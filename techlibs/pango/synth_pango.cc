
/*
 *  yosys -- Yosys Open SYnthesis Suite
 *
 *  Copyright (C) 2025  Shenzhen Pango Microsystems Co., Ltd. <marketing@pangomicro.com>
 *
 *  Permission to use, copy, modify, and/or distribute this software for any
 *  purpose with or without fee is hereby granted, provided that the above
 *  copyright notice and this permission notice appear in all copies.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 *  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 *  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 *  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 *  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 *  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

/*
Area flow based lut-mapper.
Implement of "Heuristics for Area Minimization in LUT-Based FPGA Technology Mapping".
*/

#include "kernel/celltypes.h"
#include "kernel/consteval.h"
#include "kernel/modtools.h"
#include "kernel/sigtools.h"
#include "kernel/yosys.h"
#include <queue>
#include <ranges>
#include <string.h>
#include <fstream>
#include <map>
#include <cmath>
#include <limits>
#include <bits/stdc++.h>
USING_YOSYS_NAMESPACE
using namespace std;
PRIVATE_NAMESPACE_BEGIN

// -----------------------

constexpr int LUT_SIZE = 6;
dict<SigBit, Cell*> bit2driver;
struct CutInfo {
    std::vector<SigBit> inputs;
    SigBit output;
    RTLIL::Const init;  // 存储 INIT 参数 (真值表)
};

// using sigmap to get a unique name of each signal
SigMap sigmap;
dict<std::pair<SigBit,SigBit>,vector<SigBit>> twoOutputCuts;
dict<SigBit, CutInfo> bit2cutinfo;
bool using_internel_lut_type = false;
//-------------------------------
// functions declare here
bool MapperMain(Module *module);

void SetPangoCellTypes(CellTypes *);
void BuildBit2CutInfo(Module *module);
bool Bit2oCut_FromBit2CutInfo();
bool Cone2ToLUTs(Module *module,
                dict<SigBit, pool<SigBit>> &bit2cut,
                dict<std::pair<SigBit,SigBit>, pool<SigBit>> &twoOutputCuts);
bool Bit2oCut(dict<SigBit, pool<SigBit>> &bit2cut) ;//核心



SigBit GetCellOutput(Cell *cell)
{
    // 优先使用标准输出端口 Y
    if (cell->hasPort("\\Y")) {
        RTLIL::SigSpec y = cell->getPort("\\Y");
        if (y.size() != 1)
            log_warning("GetCellOutput: cell %s output has size != 1\n", log_id(cell));
        return y.as_bit();
    }

    // 尝试常见其它单输出端口名（可按需扩展）
    static const char *fallback_ports[] = {"\\Q", "\\O", "\\OUT"};
    for (auto name : fallback_ports) {
        if (cell->hasPort(name)) {
            RTLIL::SigSpec p = cell->getPort(name);
            if (p.size() != 1)
                log_warning("GetCellOutput: cell %s output %s size != 1\n", log_id(cell), name);
            return p.as_bit();
        }
    }

    log_error("GetCellOutput: cannot find single-bit output for cell %s\n", log_id(cell));
    return SigBit(); // 永远不会到这里，log_error 会抛异常
}

//blossom_max_matching算法
pair<vector<int>, int> max_matching(const vector<vector<int>>& g) {
    int n = (int)g.size();
    vector<int> mate(n, -1), p(n), base(n);
    vector<int> q;
    vector<char> used(n), blossom(n);

    function<int(int,int)> lca = [&](int a, int b) {
        vector<char> used_path(n, false);
        while (true) {
            a = base[a];
            used_path[a] = true;
            if (mate[a] == -1) break;
            a = p[mate[a]];
        }
        while (true) {
            b = base[b];
            if (used_path[b]) return b;
            b = p[mate[b]];
        }
    };

    function<void(int,int,int)> mark_path = [&](int v, int b, int children) {
        while (base[v] != b) {
            blossom[base[v]] = blossom[base[mate[v]]] = true;
            p[v] = children;
            children = mate[v];
            v = p[mate[v]];
        }
    };

    function<int(int)> find_augmenting = [&](int root) -> int {
        used.assign(n, false);
        p.assign(n, -1);
        for (int i = 0; i < n; ++i) base[i] = i;
        q.clear();
        q.push_back(root);
        used[root] = true;
        int qh = 0;
        while (qh < (int)q.size()) {
            int v = q[qh++];
            for (int to : g[v]) {
                if (base[v] == base[to] || mate[v] == to) continue;
                if (to == root || (mate[to] != -1 && p[mate[to]] != -1)) {
                    int curbase = lca(v, to);
                    blossom.assign(n, false);
                    mark_path(v, curbase, to);
                    mark_path(to, curbase, v);
                    for (int i = 0; i < n; ++i) {
                        if (blossom[base[i]]) {
                            base[i] = curbase;
                            if (!used[i]) {
                                used[i] = true;
                                q.push_back(i);
                            }
                        }
                    }
                } else if (p[to] == -1) {
                    p[to] = v;
                    if (mate[to] == -1)
                        return to; // found augmenting path ending at 'to'
                    to = mate[to];
                    used[to] = true;
                    q.push_back(to);
                }
            }
        }
        return -1;
    };

    for (int v = 0; v < n; ++v) {
        if (mate[v] == -1) {
            int to = find_augmenting(v);
            if (to == -1) continue;
            // augment path ending at 'to'
            int cur = to;
            while (cur != -1) {
                int pv = p[cur];
                int ppv = (pv != -1 ? mate[pv] : -1);
                mate[cur] = pv;
                mate[pv] = cur;
                cur = ppv;
            }
        }
    }

    int match_pairs = 0;
    for (int i = 0; i < n; ++i)
        if (mate[i] != -1 && i < mate[i]) ++match_pairs;

    return {mate, match_pairs};
}


int IsGTP_LUT(Cell *cell)
{
	const char *type_str = cell->type.c_str();
	if (0 == strncmp(type_str, "\\GTP_LUT", 8)) {
		if (strlen(cell->type.c_str()) == 8 + 1) {
			int size = type_str[8] - '1';
			if ((size < 0 || size > 6)) {
				return 0;
			}
			return size;
		}
	}
	return 0;
}
bool IsGTP_LUT6D(Cell *cell)
{
	if (cell->type != ID(GTP_LUT6D))
		return false;
	return true;
}
bool IsGTP(Cell *cell) { return cell->type.begins_with("\\GTP_"); }


bool MapperMain(Module *module)
{		  SetPangoCellTypes(&yosys_celltypes);
        twoOutputCuts.clear();
		bit2cutinfo.clear();
        sigmap.set(module);
		BuildBit2CutInfo(module);
		Bit2oCut_FromBit2CutInfo();

        return true;
}


// ============ 全局结构与变量定义 ============
void BuildBit2CutInfo(Module *module)
{
    bit2cutinfo.clear();
    log("Building bit2cutinfo from LUTs in module: %s\n", log_id(module));

    for (auto cell : module->cells())
    {
        std::string type = cell->type.str();

        // 大小写无关匹配 LUT
        std::string type_upper;
        type_upper.reserve(type.size());
        for (char c : type) type_upper.push_back(toupper(c));

        bool is_lut = (type_upper.find("LUT") != std::string::npos)
                      || cell->type.in(ID($lut), ID($_LUT_));

        if (!is_lut)
            continue;

        // 获取输出端口
        SigBit out;
        if (cell->hasPort(ID(Y)))
            out = cell->getPort(ID(Y)).as_bit();
        else if (cell->hasPort(ID(O)))
            out = cell->getPort(ID(O)).as_bit();
        else
            continue;

        // 获取输入端口
        SigSpec inputs_spec;
        if (cell->hasPort(ID(A)))
            inputs_spec = cell->getPort(ID(A));
        else if (cell->hasPort(ID(I)))
            inputs_spec = cell->getPort(ID(I));
        else
            continue;

        int k = inputs_spec.size();
        if (k == 0 || k > 6)
            continue;

        std::vector<SigBit> inputs;
        inputs.reserve(k);
        for (auto b : inputs_spec)
            inputs.push_back(b);

        // 获取 INIT 参数
        RTLIL::Const init;
        if (cell->hasParam(ID(INIT)))
            init = cell->getParam(ID(INIT));
        else if (cell->hasParam(ID(MASK)))
            init = cell->getParam(ID(MASK));
        else if (cell->hasParam(ID(LUT)))
            init = cell->getParam(ID(LUT));
        else
            continue;

        if ((int)init.bits().size() != (1 << k)) {
            log_warning("LUT %s: INIT size mismatch (%d vs %d)\n",
                        log_id(cell), int(init.bits().size()), 1 << k);
            continue;
        }

        // 构建 CutInfo 并插入全局 bit2cutinfo
        CutInfo info;
        info.output = out;
        info.inputs = inputs;
        info.init = init;
        bit2cutinfo[out] = info;

        // 调试输出
        // log("  [%s] LUT%d output=%s inputs(%d)=",
        //     log_id(cell), k, log_signal(out), k);
        // for (auto &b : inputs)
        //     log(" %s", log_signal(b));
        // log("\n    INIT = ");
        // for (int i = (1 << k) - 1; i >= 0; --i)
        //     log("%c", init.bits()[i] == State::S1 ? '1' : '0');
        // log("\n");
    }

    log("BuildBit2CutInfo: total LUTs = %d\n", GetSize(bit2cutinfo));
}
bool Bit2oCut(dict<SigBit, pool<SigBit>> &bit2cut)
{
    bit2cut.clear();
    pool<SigBit> primary_inputs;
    pool<Cell*> topo_cells;

    // 遍历电路，找出 primary input（无驱动）
    for (auto &p : bit2cutinfo) {
        for (auto in : p.second.inputs)
            if (!bit2cutinfo.count(in))
                primary_inputs.insert(in);
    }

    // 初始化 PI cut
    for (auto &pi : primary_inputs)
        bit2cut[pi].insert(pi);

    // 拓扑遍历（假设 bit2cutinfo 已是 DAG）
    std::queue<SigBit> q;
    for (auto &p : bit2cutinfo)
        q.push(p.first);

    while (!q.empty()) {
        SigBit out = q.front();
        q.pop();

        const CutInfo &ci = bit2cutinfo[out];
        pool<SigBit> merged_cut;
        for (auto in : ci.inputs) {
            merged_cut.insert(bit2cut[in].begin(), bit2cut[in].end());
        }

        // 限制 cut 大小不超过 6
        while (merged_cut.size() > 6) {
            // 剪枝策略：删除扇出最小或随机节点
            merged_cut.erase(merged_cut.begin());
        }

        bit2cut[out] = merged_cut;
    }

    log("Bit2oCut: total cuts = %d\n", GetSize(bit2cut));
    return true;
}

State StateEval(const dict<SigBit, State> &assign, SigBit bit)
{
    // 已赋值直接返回
    if (assign.count(bit))
        return assign.at(bit);

    // 如果 bit 是常量
    if (bit == State::S0) return State::S0;
    if (bit == State::S1) return State::S1;

    // 从 bit2cutinfo 查找驱动 LUT
    if (!bit2cutinfo.count(bit))
        return State::Sx; // 无定义

    const CutInfo &ci = bit2cutinfo.at(bit);
    int k = ci.inputs.size();

    int idx = 0;
    for (int i = 0; i < k; i++) {
        State s = StateEval(assign, ci.inputs[i]);
        if (s == State::Sx) return State::Sx;
        if (s == State::S1)
            idx |= (1 << i);
    }

   return ci.init[idx];



}
std::vector<bool> GetCutInit(const std::vector<SigBit> &cut, SigBit output)
{
    int n = cut.size();
    int mask_num = 1 << n;
    std::vector<bool> init(mask_num);

    for (int m = 0; m < mask_num; ++m) {
        dict<SigBit, State> assign;
        for (int i = 0; i < n; i++)
            assign[cut[i]] = ((m >> i) & 1) ? State::S1 : State::S0;

        State val = StateEval(assign, output);
        init[m] = (val == State::S1);
    }
    return init;
}
bool IsCombinationalGate(Cell *cell)
{
    static const pool<IdString> seq_types = {
        ID(GTP_DFF), ID(GTP_DFF_C), ID(GTP_DFF_R),
        ID(GTP_DFF_S), ID(GTP_DLATCH), ID(GTP_RAM32X1SP),
        ID(GTP_RAM64X1SP)
    };

    if (seq_types.count(cell->type))
        return false;

    if (cell->type.str().find("LUT") != std::string::npos)
        return true;

    if (cell->type.str().find("BUF") != std::string::npos)
        return false;

    // 默认认为有输入输出的纯组合门
    return cell->hasPort(ID(A)) || cell->hasPort(ID(I));
}


bool HasCommonLeaf(const pool<SigBit> &cut1, const pool<SigBit> &cut2)
{
    // 选择较小和较大的集合
    const pool<SigBit> &smaller = (cut1.size() <= cut2.size()) ? cut1 : cut2;
    const pool<SigBit> &larger  = (cut1.size() <= cut2.size()) ? cut2 : cut1;

    for (auto &sig : smaller) {
        if (larger.count(sig)) {
            return true; // 有一个包含
        }
    }
    return false; 
}


int SortCutLevel(const vector<Cell*> &cells_in_level,
                  dict<SigBit, pool<SigBit>> &bit2cut,
                  vector<pair<Cell*, int>> &cell_pi_list)
{
    cell_pi_list.clear();
   int  num_six_input = 0;

    vector<pair<Cell*, int>> six_list;
    vector<pair<Cell*, int>> other_list;

    for (auto *cell : cells_in_level) {
        SigBit out = GetCellOutput(cell);
        int pi_size = bit2cut[out].size();
        if (pi_size == 6)
            six_list.push_back({cell, pi_size});
        else
            other_list.push_back({cell, pi_size});
    }

    // 6-input 的放前面
    cell_pi_list.insert(cell_pi_list.end(), six_list.begin(), six_list.end());
    cell_pi_list.insert(cell_pi_list.end(), other_list.begin(), other_list.end());

    num_six_input = six_list.size();
    return num_six_input;
}
//补充：等价验证函数，对于6picut进行验证
//枚举并按边存储可行性对
//图算法求得最优组队并存入twoOutputCuts
// 全局结构

bool Bit2oCut_FromBit2CutInfo()
{
    log("Bit2oCut_FromBit2CutInfo: entry, total LUTs = %d\n", GetSize(bit2cutinfo));
    twoOutputCuts.clear();

    // 将所有项转为 vector 方便遍历
    std::vector<std::pair<SigBit, CutInfo>> entries(bit2cutinfo.begin(), bit2cutinfo.end());
    int n = entries.size();

    int merge_cnt = 0;

    for (int i = 0; i < n; i++) {
        const auto &A = entries[i].second;
        const std::vector<SigBit> &inA = A.inputs;
		if(inA.size()==6) continue;
        for (int j = i + 1; j < n; j++) {
            const auto &B = entries[j].second;
            const std::vector<SigBit> &inB = B.inputs;
						if(inB.size()==6) continue;
            // 计算交集
            bool has_common = false;
            for (auto &a : inA)
                for (auto &b : inB)
                    if (a == b) { has_common = true; break; }
            if (!has_common) continue;

            // 构造输入并集
            std::vector<SigBit> merged = inA;
            for (auto &b : inB)
                if (std::find(merged.begin(), merged.end(), b) == merged.end())
                    merged.push_back(b);

            // 限制：最多6输入
            if ((int)merged.size() > 5) continue;

            // 记录融合
            twoOutputCuts[{A.output, B.output}] = merged;
            merge_cnt++;

            // 调试输出
          //  log("  [FUSE] %s + %s -> merged_inputs(%d):",
                // log_signal(A.output), log_signal(B.output), (int)merged.size());
            // for (auto &m : merged) log(" %s", log_signal(m));
            // log("\n");

            // 一个输出只参与一次融合（可改）
            break;
        }
    }

    log("Bit2oCut_FromBit2CutInfo: total merged pairs = %d\n", merge_cnt);
    return true;
}



void SetPangoCellTypes(CellTypes *ct)
{
	ct->setup_type(ID(GTP_ADC_E2),
		       {ID(VA), ID(VAUX), ID(DCLK), ID(DADDR), ID(DEN), ID(SECEN), ID(DWE), ID(DI), ID(CONVST), ID(RST_N), ID(LOADSC_N)},
		       {ID(DO), ID(DRDY), ID(OVER_TEMP), ID(LOGIC_DONE_A), ID(LOGIC_DONE_B), ID(ADC_CLK_OUT), ID(DMODIFIED), ID(ALARM)}, false);
	ct->setup_type(ID(GTP_APM_E2), {ID(CIN),   ID(CPI),   ID(CXI),	ID(CXBI),  ID(X),	 ID(XB),	ID(Y),	      ID(Z),	   ID(MODEIN),
					ID(MODEY), ID(MODEZ), ID(CLK),	ID(CEX1),  ID(CEX2),	 ID(CEX3),	ID(CEXB),     ID(CEY1),	   ID(CEY2),
					ID(CEZ),   ID(CEM),   ID(CEP),	ID(CEPRE), ID(CEMODEIN), ID(CEMODEY),	ID(CEMODEZ),  ID(RSTX),	   ID(RSTXB),
					ID(RSTY),  ID(RSTZ),  ID(RSTM), ID(RSTP),  ID(RSTPRE),	 ID(RSTMODEIN), ID(RSTMODEY), ID(RSTMODEZ)},
		       {ID(COUT), ID(CPO), ID(CXO), ID(CXBO), ID(P)}, false);
	ct->setup_type(ID(GTP_BUF), {ID(I)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_CFGCLK), {ID(CLKIN), ID(CE_N)}, {}, false);
	ct->setup_type(ID(GTP_CLKBUFCE), {ID(CE), ID(CLKIN)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFG), {ID(CLKIN)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFGCE), {ID(CLKIN), ID(CE)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFGMUX), {ID(CLKIN0), ID(CLKIN1), ID(SEL)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFGMUX_E1), {ID(CLKIN0), ID(CLKIN1), ID(SEL), ID(EN)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFGMUX_E2), {ID(CLKIN0), ID(CLKIN1), ID(DETECT_CLK0), ID(DETECT_CLK1), ID(SEL)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFM), {ID(CLKIN)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFMCE), {ID(CLKIN), ID(CE)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFR), {ID(CLKIN)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFX), {ID(CLKIN)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKBUFXCE), {ID(CLKIN), ID(CE)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_CLKPD), {ID(RST), ID(CLK_SAMPLE), ID(CLK_CTRL), ID(CLK_PHY), ID(DONE)}, {ID(FLAG_PD), ID(LOCK)}, false);
	ct->setup_type(ID(GTP_DDC_E2),
		       {ID(RST), ID(RST_TRAINING_N), ID(CLKA), ID(CLKB), ID(DQSI), ID(DQSIB), ID(DELAY_STEP0), ID(DELAY_STEP1), ID(DELAY_STEP2),
			ID(DELAY_STEP3), ID(DELAY_STEP4), ID(DQS_GATE_CTRL), ID(GATE_SEL), ID(CLK_GATE_CTRL), ID(CLKA_GATE)},
		       {ID(WCLK), ID(WCLK_DELAY), ID(DQSI_DELAY), ID(DQSIB_DELAY), ID(DGTS), ID(IFIFO_WADDR), ID(IFIFO_RADDR), ID(READ_VALID),
			ID(DQS_DRIFT), ID(DRIFT_DETECT_ERR), ID(DQS_DRIFT_STATUS), ID(DQS_SAMPLE)},
		       false);
	ct->setup_type(ID(GTP_DDC_E2_DFT),
		       {ID(RST), ID(RST_TRAINING_N), ID(CLKA), ID(CLKB), ID(DQSI), ID(DQSIB), ID(DELAY_STEP0), ID(DELAY_STEP1), ID(DELAY_STEP2),
			ID(DELAY_STEP3), ID(DELAY_STEP4), ID(DQS_GATE_CTRL), ID(GATE_SEL), ID(CLK_GATE_CTRL), ID(CLKA_GATE)},
		       {ID(WCLK), ID(WCLK_DELAY), ID(DQSI_DELAY), ID(DQSIB_DELAY), ID(DGTS), ID(IFIFO_WADDR), ID(IFIFO_RADDR), ID(READ_VALID),
			ID(DQS_DRIFT), ID(DRIFT_DETECT_ERR), ID(DQS_DRIFT_STATUS), ID(DQS_SAMPLE), ID(GATE_HIGHB), ID(GATE_HIGH_LATCHB)},
		       false);
	ct->setup_type(ID(GTP_DDC_E3),
		       {ID(RST), ID(RST_TRAINING_N), ID(CLKA), ID(CLKB), ID(DQSI), ID(DQSIB), ID(DELAY_STEP0), ID(DELAY_STEP1), ID(DELAY_STEP2),
			ID(DELAY_STEP3), ID(DELAY_STEP4), ID(DQS_GATE_CTRL), ID(GATE_SEL), ID(CLK_GATE_CTRL), ID(CLKA_GATE)},
		       {ID(WCLK), ID(WCLK_DELAY), ID(DQSI_DELAY), ID(DQSIB_DELAY), ID(DGTS), ID(IFIFO_WADDR), ID(IFIFO_RADDR), ID(READ_VALID),
			ID(DQS_DRIFT), ID(DRIFT_DETECT_ERR), ID(DQS_DRIFT_STATUS), ID(DQS_SAMPLE)},
		       false);
	ct->setup_type(ID(GTP_DFF), {ID(D), ID(CLK)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DFF_C), {ID(D), ID(CLK), ID(C)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DFF_CE), {ID(D), ID(CLK), ID(C), ID(CE)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DFF_E), {ID(D), ID(CLK), ID(CE)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DFF_P), {ID(D), ID(CLK), ID(P)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DFF_PE), {ID(D), ID(CLK), ID(P), ID(CE)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DFF_R), {ID(D), ID(CLK), ID(R)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DFF_RE), {ID(D), ID(CLK), ID(R), ID(CE)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DFF_S), {ID(D), ID(CLK), ID(S)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DFF_SE), {ID(D), ID(CLK), ID(S), ID(CE)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DLATCH), {ID(D), ID(G)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DLATCH_C), {ID(D), ID(G), ID(C)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DLATCH_CE), {ID(D), ID(G), ID(C), ID(GE)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DLATCH_E), {ID(D), ID(G), ID(GE)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DLATCH_P), {ID(D), ID(G), ID(P)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DLATCH_PE), {ID(D), ID(G), ID(P), ID(GE)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_DLL_E2), {ID(CLKIN), ID(SYS_CLK), ID(PWD), ID(RST), ID(UPDATE_N)}, {ID(DELAY_STEP), ID(DELAY_STEP1), ID(LOCK)}, false);
	ct->setup_type(ID(GTP_DRM18K_E1),
		       {ID(DIA), ID(DIB), ID(ADDRA), ID(ADDRA_HOLD), ID(ADDRB), ID(ADDRB_HOLD), ID(BWEA), ID(BWEB), ID(CLKA), ID(CLKB), ID(CEA),
			ID(CEB), ID(WEA), ID(WEB), ID(ORCEA), ID(ORCEB), ID(RSTA), ID(RSTB)},
		       {ID(DOA), ID(DOB)}, false);
	ct->setup_type(ID(GTP_DRM36K_E1),
		       {ID(CINA),  ID(CINB),  ID(DIA),	ID(DIB),  ID(ADDRA),	      ID(ADDRA_HOLD),	 ID(ADDRB), ID(ADDRB_HOLD), ID(CSA),
			ID(CSB),   ID(BWEA),  ID(BWEB), ID(CLKA), ID(CLKB),	      ID(CEA),		 ID(CEB),   ID(WEA),	    ID(WEB),
			ID(ORCEA), ID(ORCEB), ID(RSTA), ID(RSTB), ID(INJECT_SBITERR), ID(INJECT_DBITERR)},
		       {ID(COUTA), ID(COUTB), ID(DOA), ID(DOB), ID(ECC_SBITERR), ID(ECC_DBITERR), ID(ECC_RDADDR), ID(ECC_PARITY)}, false);
	ct->setup_type(ID(GTP_EFUSECODE), {}, {ID(EFUSE_CODE)}, false);
	ct->setup_type(ID(GTP_FIFO18K_E1), {ID(DI), ID(WCLK), ID(RCLK), ID(WCE), ID(RCE), ID(ORCE), ID(RST)},
		       {ID(ALMOST_EMPTY), ID(ALMOST_FULL), ID(EMPTY), ID(FULL), ID(DO)}, false);
	ct->setup_type(ID(GTP_FIFO36K_E1), {ID(DI), ID(WCLK), ID(RCLK), ID(WCE), ID(RCE), ID(ORCE), ID(RST), ID(INJECT_SBITERR), ID(INJECT_DBITERR)},
		       {ID(ALMOST_EMPTY), ID(ALMOST_FULL), ID(EMPTY), ID(FULL), ID(DO), ID(ECC_SBITERR), ID(ECC_DBITERR)}, false);
	ct->setup_type(ID(GTP_GPLL), {ID(CLKIN1),      ID(CLKIN2),	ID(CLKFB),	 ID(CLKIN_SEL),	  ID(DPS_CLK),	   ID(DPS_EN),
				      ID(DPS_DIR),     ID(CLKOUT0_SYN), ID(CLKOUT1_SYN), ID(CLKOUT2_SYN), ID(CLKOUT3_SYN), ID(CLKOUT4_SYN),
				      ID(CLKOUT5_SYN), ID(CLKOUT6_SYN), ID(CLKOUTF_SYN), ID(PLL_PWD),	  ID(RST),	   ID(APB_CLK),
				      ID(APB_RST_N),   ID(APB_ADDR),	ID(APB_SEL),	 ID(APB_EN),	  ID(APB_WRITE),   ID(APB_WDATA)},
		       {ID(CLKOUT0), ID(CLKOUT0N), ID(CLKOUT1), ID(CLKOUT1N), ID(CLKOUT2), ID(CLKOUT2N), ID(CLKOUT3), ID(CLKOUT3N), ID(CLKOUT4),
			ID(CLKOUT5), ID(CLKOUT6), ID(CLKOUTF), ID(CLKOUTFN), ID(LOCK), ID(DPS_DONE), ID(APB_RDATA), ID(APB_READY)},
		       false);
	ct->setup_type(ID(GTP_GRS), {ID(GRS_N)}, {}, false);
	ct->setup_type(ID(GTP_HPIO_VREF), {ID(CODE_VREF0), ID(CODE_VREF1), ID(CODE_VREF2), ID(CODE_VREF3)}, {}, false);
	ct->setup_type(ID(GTP_HSSTHP_BUFDS), {ID(COM_POWERDOWN), ID(PAD_REFCLKP), ID(PAD_REFCLKN)}, {ID(REFCLK_OUTP), ID(PMA_REFCLK_TO_FABRIC)},
		       false);
	ct->setup_type(ID(GTP_HSSTHP_HPLL),
		       {ID(P_CFG_RST_HPLL),
			ID(P_CFG_CLK_HPLL),
			ID(P_CFG_PSEL_HPLL),
			ID(P_CFG_ENABLE_HPLL),
			ID(P_CFG_WRITE_HPLL),
			ID(P_CFG_ADDR_HPLL),
			ID(P_CFG_WDATA_HPLL),
			ID(P_HPLL_POWERDOWN),
			ID(P_HPLL_RST),
			ID(P_HPLL_LOCKDET_RST),
			ID(P_RES_CAL_RST),
			ID(P_TX_SYNC),
			ID(P_HPLL_DIV_SYNC),
			ID(P_REFCLK_DIV_SYNC),
			ID(P_HPLL_VCO_CALIB_EN),
			ID(P_RESCAL_I_CODE_I),
			ID(P_HPLL_DIV_CHANGE),
			ID(P_HPLL_REFCLK_I),
			ID(REFCLK_TO_TX_SYNC_I),
			ID(REFCLK_TO_REFCLK_SYNC_I),
			ID(REFCLK_TO_DIV_SYNC_I),
			ID(ANA_TX_SYNC_I),
			ID(ANA_HPLL_REFCLK_SYNC_I),
			ID(ANA_HPLL_DIV_SYNC_I),
			ID(P_PLL_REFCLK6_I),
			ID(REFCLK0),
			ID(REFCLK1),
			ID(REFCLK0_FROM_UPPER_HSST),
			ID(REFCLK1_FROM_UPPER_HSST),
			ID(REFCLK0_FROM_LOWER_HSST),
			ID(REFCLK1_FROM_LOWER_HSST),
			ID(REFCLK_FROM_FABRIC),
			ID(TX_SYNC_FROM_UPPER_HSST),
			ID(TX_SYNC_FROM_LOWER_HSST),
			ID(REFCLK_SYNC_FROM_UPPER_HSST),
			ID(REFCLK_SYNC_FROM_LOWER_HSST),
			ID(DIV_SYNC_FROM_UPPER_HSST),
			ID(DIV_SYNC_FROM_LOWER_HSST)},
		       {ID(P_CFG_READY_HPLL),
			ID(P_CFG_RDATA_HPLL),
			ID(P_CFG_INT_HPLL),
			ID(P_RES_CAL_CODE_FABRIC),
			ID(P_HPLL_READY),
			ID(PMA_HPLL_READY_O),
			ID(PMA_HPLL_REFCLK_O),
			ID(PMA_TX_SYNC_HPLL_O),
			ID(PMA_HPLL_CK0),
			ID(PMA_HPLL_CK90),
			ID(PMA_HPLL_CK180),
			ID(PMA_HPLL_CK270),
			ID(TX_SYNC_REFSYNC_O),
			ID(REFCLK_SYNC_REFSYNC_O),
			ID(DIV_SYNC_REFSYNC_O),
			ID(LPLL_REFCKOUT_CH0),
			ID(LPLL_REFCKOUT_CH1),
			ID(LPLL_REFCKOUT_CH2),
			ID(LPLL_REFCKOUT_CH3),
			ID(REFCLK0_FOR_UPPER_HSST),
			ID(REFCLK1_FOR_UPPER_HSST),
			ID(REFCLK0_FOR_LOWER_HSST),
			ID(REFCLK1_FOR_LOWER_HSST),
			ID(TX_SYNC_FOR_UPPER_HSST),
			ID(TX_SYNC_FOR_LOWER_HSST),
			ID(ANA_TX_SYNC_O),
			ID(REFCLK_SYNC_FOR_UPPER_HSST),
			ID(REFCLK_SYNC_FOR_LOWER_HSST),
			ID(DIV_SYNC_FOR_UPPER_HSST),
			ID(DIV_SYNC_FOR_LOWER_HSST)},
		       false);
	ct->setup_type(ID(GTP_HSSTHP_HPLL_DFT),
		       {ID(P_CFG_RST_HPLL),
			ID(P_CFG_CLK_HPLL),
			ID(P_CFG_PSEL_HPLL),
			ID(P_CFG_ENABLE_HPLL),
			ID(P_CFG_WRITE_HPLL),
			ID(P_CFG_ADDR_HPLL),
			ID(P_CFG_WDATA_HPLL),
			ID(P_HPLL_POWERDOWN),
			ID(P_HPLL_RST),
			ID(P_HPLL_LOCKDET_RST),
			ID(P_RES_CAL_RST),
			ID(P_TX_SYNC),
			ID(P_HPLL_DIV_SYNC),
			ID(P_REFCLK_DIV_SYNC),
			ID(P_HPLL_VCO_CALIB_EN),
			ID(P_RESCAL_I_CODE_I),
			ID(P_HPLL_DIV_CHANGE),
			ID(P_HPLL_REFCLK_I),
			ID(REFCLK_TO_TX_SYNC_I),
			ID(REFCLK_TO_REFCLK_SYNC_I),
			ID(REFCLK_TO_DIV_SYNC_I),
			ID(ANA_TX_SYNC_I),
			ID(ANA_HPLL_REFCLK_SYNC_I),
			ID(ANA_HPLL_DIV_SYNC_I),
			ID(P_PLL_REFCLK6_I),
			ID(REFCLK0),
			ID(REFCLK1),
			ID(REFCLK0_FROM_UPPER_HSST),
			ID(REFCLK1_FROM_UPPER_HSST),
			ID(REFCLK0_FROM_LOWER_HSST),
			ID(REFCLK1_FROM_LOWER_HSST),
			ID(REFCLK_FROM_FABRIC),
			ID(TX_SYNC_FROM_UPPER_HSST),
			ID(TX_SYNC_FROM_LOWER_HSST),
			ID(REFCLK_SYNC_FROM_UPPER_HSST),
			ID(REFCLK_SYNC_FROM_LOWER_HSST),
			ID(DIV_SYNC_FROM_UPPER_HSST),
			ID(DIV_SYNC_FROM_LOWER_HSST),
			ID(P_TEST_SE_N),
			ID(P_TEST_MODE_N),
			ID(P_TEST_RSTN),
			ID(P_TEST_SI0),
			ID(P_TEST_SI1),
			ID(P_FOR_PMA_TEST_MODE_N),
			ID(P_FOR_PMA_TEST_SE_N),
			ID(P_FOR_PMA_TEST_CLK),
			ID(P_FOR_PMA_TEST_RSTN),
			ID(P_FOR_PMA_TEST_SI)},
		       {ID(P_CFG_READY_HPLL),
			ID(P_CFG_RDATA_HPLL),
			ID(P_CFG_INT_HPLL),
			ID(P_RES_CAL_CODE_FABRIC),
			ID(P_HPLL_READY),
			ID(PMA_HPLL_READY_O),
			ID(PMA_HPLL_REFCLK_O),
			ID(PMA_TX_SYNC_HPLL_O),
			ID(PMA_HPLL_CK0),
			ID(PMA_HPLL_CK90),
			ID(PMA_HPLL_CK180),
			ID(PMA_HPLL_CK270),
			ID(TX_SYNC_REFSYNC_O),
			ID(REFCLK_SYNC_REFSYNC_O),
			ID(DIV_SYNC_REFSYNC_O),
			ID(LPLL_REFCKOUT_CH0),
			ID(LPLL_REFCKOUT_CH1),
			ID(LPLL_REFCKOUT_CH2),
			ID(LPLL_REFCKOUT_CH3),
			ID(REFCLK0_FOR_UPPER_HSST),
			ID(REFCLK1_FOR_UPPER_HSST),
			ID(REFCLK0_FOR_LOWER_HSST),
			ID(REFCLK1_FOR_LOWER_HSST),
			ID(TX_SYNC_FOR_UPPER_HSST),
			ID(TX_SYNC_FOR_LOWER_HSST),
			ID(ANA_TX_SYNC_O),
			ID(REFCLK_SYNC_FOR_UPPER_HSST),
			ID(REFCLK_SYNC_FOR_LOWER_HSST),
			ID(DIV_SYNC_FOR_UPPER_HSST),
			ID(DIV_SYNC_FOR_LOWER_HSST),
			ID(P_TEST_SO0),
			ID(P_TEST_SO1),
			ID(P_FOR_PMA_TEST_SO)},
		       false);
	ct->setup_type(ID(GTP_HSSTHP_HPLL_t),
		       {ID(P_CFG_RST_HPLL),
			ID(P_CFG_CLK_HPLL),
			ID(P_CFG_PSEL_HPLL),
			ID(P_CFG_ENABLE_HPLL),
			ID(P_CFG_WRITE_HPLL),
			ID(P_CFG_ADDR_HPLL),
			ID(P_CFG_WDATA_HPLL),
			ID(P_COM_POWERDOWN),
			ID(P_HPLL_POWERDOWN),
			ID(P_HPLL_RST),
			ID(P_HPLL_LOCKDET_RST),
			ID(P_RES_CAL_RST),
			ID(P_TX_SYNC),
			ID(P_TX_RATE_CHANGE_ON_0),
			ID(P_TX_RATE_CHANGE_ON_1),
			ID(P_HPLL_DIV_SYNC),
			ID(P_REFCLK_DIV_SYNC),
			ID(P_HPLL_VCO_CALIB_EN),
			ID(P_RESCAL_I_CODE_I),
			ID(P_HPLL_REF_CLK),
			ID(P_HPLL_DIV_CHANGE),
			ID(P_FROM_LOWER_HSST_BUS),
			ID(P_FROM_UPPER_HSST_BUS),
			ID(PAD_REFCLKN_0),
			ID(PAD_REFCLKP_0),
			ID(PAD_REFCLKN_1),
			ID(PAD_REFCLKP_1)},
		       {ID(P_CFG_READY_HPLL),
			ID(P_CFG_RDATA_HPLL),
			ID(P_CFG_INT_HPLL),
			ID(P_RES_CAL_CODE_FABRIC),
			ID(P_REFCK2CORE_0),
			ID(P_REFCK2CORE_1),
			ID(P_HPLL_READY),
			ID(PMA_HPLL_READY_LEFT),
			ID(PMA_HPLL_READY_RIGHT),
			ID(PMA_HPLL_REFCLK_LEFT),
			ID(PMA_HPLL_REFCLK_RIGHT),
			ID(PMA_LPLL_REFCKOUT_CH0),
			ID(PMA_LPLL_REFCKOUT_CH1),
			ID(PMA_LPLL_REFCKOUT_CH2),
			ID(PMA_LPLL_REFCKOUT_CH3),
			ID(PMA_RES_CAL_LEFT),
			ID(PMA_RES_CAL_RIGHT),
			ID(PMA_TX_RATE_CHANGE_ON0_LEFT),
			ID(PMA_TX_RATE_CHANGE_ON0_RIGHT),
			ID(PMA_TX_RATE_CHANGE_ON1_LEFT),
			ID(PMA_TX_RATE_CHANGE_ON1_RIGHT),
			ID(PMA_TX_SYNC_HPLL_LEFT),
			ID(PMA_TX_SYNC_HPLL_RIGHT),
			ID(PMA_TX_SYNC_LEFT),
			ID(PMA_TX_SYNC_RIGHT),
			ID(PMA_HPLL_CK0_CH0),
			ID(PMA_HPLL_CK0_CH1),
			ID(PMA_HPLL_CK0_CH2),
			ID(PMA_HPLL_CK0_CH3),
			ID(PMA_HPLL_CK90_CH0),
			ID(PMA_HPLL_CK90_CH1),
			ID(PMA_HPLL_CK90_CH2),
			ID(PMA_HPLL_CK90_CH3),
			ID(PMA_HPLL_CK180_CH0),
			ID(PMA_HPLL_CK180_CH1),
			ID(PMA_HPLL_CK180_CH2),
			ID(PMA_HPLL_CK180_CH3),
			ID(PMA_HPLL_CK270_CH0),
			ID(PMA_HPLL_CK270_CH1),
			ID(PMA_HPLL_CK270_CH2),
			ID(PMA_HPLL_CK270_CH3),
			ID(P_FOR_LOWER_HSST_BUS),
			ID(P_FOR_UPPER_HSST_BUS)},
		       false);
	ct->setup_type(ID(GTP_HSSTHP_LANE),
		       {ID(P_RX_CLK_FR_CORE),
			ID(P_RCLK2_FR_CORE),
			ID(P_TX_CLK_FR_CORE),
			ID(P_TCLK2_FR_CORE),
			ID(P_PCS_RX_RST),
			ID(P_PCS_TX_RST),
			ID(P_EXT_BRIDGE_PCS_RST),
			ID(P_CFG_RST),
			ID(P_CFG_CLK),
			ID(P_CFG_PSEL),
			ID(P_CFG_ENABLE),
			ID(P_CFG_WRITE),
			ID(P_CFG_ADDR),
			ID(P_CFG_WDATA),
			ID(LANE_CIN_BUS_FORWARD),
			ID(LANE_CIN_BUS_BACKWARD),
			ID(P_TDATA),
			ID(P_PCIE_EI_H),
			ID(P_PCIE_EI_L),
			ID(P_TX_DEEMP),
			ID(P_TX_DEEMP_POST_SEL),
			ID(P_BLK_ALIGN_CTRL),
			ID(P_TX_ENC_TYPE),
			ID(P_RX_DEC_TYPE),
			ID(P_PCS_BIT_SLIP),
			ID(P_PCS_WORD_ALIGN_EN),
			ID(P_RX_POLARITY_INVERT),
			ID(P_PCS_MCB_EXT_EN),
			ID(P_PCS_NEAREND_LOOP),
			ID(P_PCS_FAREND_LOOP),
			ID(P_PMA_NEAREND_PLOOP),
			ID(P_PMA_NEAREND_SLOOP),
			ID(P_PMA_FAREND_PLOOP),
			ID(P_PCS_PRBS_EN),
			ID(P_LANE_POWERDOWN),
			ID(P_LANE_RST),
			ID(P_RX_LANE_POWERDOWN),
			ID(P_RX_PMA_RST),
			ID(P_RX_CDR_RST),
			ID(P_RX_CLKPATH_RST),
			ID(P_RX_DFE_RST),
			ID(P_RX_LEQ_RST),
			ID(P_RX_SLIDING_RST),
			ID(P_RX_DFE_EN),
			ID(P_RX_T1_EN),
			ID(P_RX_CDRX_EN),
			ID(P_RX_T1_DFE_EN),
			ID(P_RX_T2_DFE_EN),
			ID(P_RX_T3_DFE_EN),
			ID(P_RX_T4_DFE_EN),
			ID(P_RX_T5_DFE_EN),
			ID(P_RX_T6_DFE_EN),
			ID(P_RX_SLIDING_EN),
			ID(P_RX_EYE_RST),
			ID(P_RX_EYE_EN),
			ID(P_RX_EYE_TAP),
			ID(P_RX_PIC_EYE),
			ID(P_RX_PIC_FASTLOCK),
			ID(P_RX_PIC_FASTLOCK_STROBE),
			ID(P_EM_RD_TRIGGER),
			ID(P_EM_MODE_CTRL),
			ID(P_RX_CTLE_DCCAL_RST),
			ID(P_RX_SLICER_DCCAL_RST),
			ID(P_RX_SLICER_DCCAL_EN),
			ID(P_RX_CTLE_DCCAL_EN),
			ID(P_RX_SLIP_RST),
			ID(P_RX_SLIP_EN),
			ID(P_LPLL_POWERDOWN),
			ID(P_LPLL_RST),
			ID(P_LPLL_LOCKDET_RST),
			ID(P_TX_LS_DATA),
			ID(P_TX_BEACON_EN),
			ID(P_TX_SWING),
			ID(P_TX_RXDET_REQ),
			ID(P_TX_RATE),
			ID(P_TX_BUSWIDTH),
			ID(P_TX_FREERUN_BUSWIDTH),
			ID(P_TX_MARGIN),
			ID(P_TX_PMA_RST),
			ID(P_TX_LANE_POWERDOWN),
			ID(P_TX_PIC_EN),
			ID(P_RX_RATE),
			ID(P_RX_BUSWIDTH),
			ID(P_RX_HIGHZ),
			ID(P_CIM_CLK_ALIGNER_RX),
			ID(P_CIM_CLK_ALIGNER_TX),
			ID(P_ALIGN_MODE_VALID_RX),
			ID(P_ALIGN_MODE_RX),
			ID(P_ALIGN_MODE_VALID_TX),
			ID(P_ALIGN_MODE_TX),
			ID(PMA_HPLL_CK0),
			ID(PMA_HPLL_CK90),
			ID(PMA_HPLL_CK180),
			ID(PMA_HPLL_CK270),
			ID(PMA_HPLL_READY_IN),
			ID(PMA_HPLL_REFCLK_IN),
			ID(PMA_TX_SYNC_HPLL_IN),
			ID(P_LPLL_REFCLK_IN),
			ID(P_TX_RATE_CHANGE_ON_0),
			ID(P_TX_RATE_CHANGE_ON_1),
			ID(P_TX_SYNC),
			ID(P_RX_SDN),
			ID(P_RX_SDP)},
		       {ID(P_CFG_READY),
			ID(P_CFG_RDATA),
			ID(P_CFG_INT),
			ID(LANE_COUT_BUS_FORWARD),
			ID(LANE_COUT_BUS_BACKWARD),
			ID(P_RX_PRBS_ERROR),
			ID(P_PCS_RX_MCB_STATUS),
			ID(P_PCS_LSM_SYNCED),
			ID(P_RDATA),
			ID(P_RXDVLD),
			ID(P_RXDVLD_H),
			ID(P_RXSTATUS),
			ID(P_EM_ERROR_CNT),
			ID(P_LPLL_READY),
			ID(P_RX_SIGDET_STATUS),
			ID(P_RX_SATA_COMINIT),
			ID(P_RX_SATA_COMWAKE),
			ID(P_RX_LS_DATA),
			ID(P_RX_READY),
			ID(P_TEST_STATUS),
			ID(P_TX_RXDET_STATUS),
			ID(P_RCLK2FABRIC),
			ID(P_TCLK2FABRIC),
			ID(P_CA_ALIGN_RX),
			ID(P_CA_ALIGN_TX),
			ID(P_TX_SDN),
			ID(P_TX_SDP)},
		       false);
	ct->setup_type(ID(GTP_HSSTHP_LANE_DFT),
		       {ID(P_RX_CLK_FR_CORE),
			ID(P_RCLK2_FR_CORE),
			ID(P_TX_CLK_FR_CORE),
			ID(P_TCLK2_FR_CORE),
			ID(P_PCS_RX_RST),
			ID(P_PCS_TX_RST),
			ID(P_EXT_BRIDGE_PCS_RST),
			ID(P_CFG_RST),
			ID(P_CFG_CLK),
			ID(P_CFG_PSEL),
			ID(P_CFG_ENABLE),
			ID(P_CFG_WRITE),
			ID(P_CFG_ADDR),
			ID(P_CFG_WDATA),
			ID(LANE_CIN_BUS_FORWARD),
			ID(LANE_CIN_BUS_BACKWARD),
			ID(P_TDATA),
			ID(P_PCIE_EI_H),
			ID(P_PCIE_EI_L),
			ID(P_TX_DEEMP),
			ID(P_TX_DEEMP_POST_SEL),
			ID(P_BLK_ALIGN_CTRL),
			ID(P_TX_ENC_TYPE),
			ID(P_RX_DEC_TYPE),
			ID(P_PCS_BIT_SLIP),
			ID(P_PCS_WORD_ALIGN_EN),
			ID(P_RX_POLARITY_INVERT),
			ID(P_PCS_MCB_EXT_EN),
			ID(P_PCS_NEAREND_LOOP),
			ID(P_PCS_FAREND_LOOP),
			ID(P_PMA_NEAREND_PLOOP),
			ID(P_PMA_NEAREND_SLOOP),
			ID(P_PMA_FAREND_PLOOP),
			ID(P_PCS_PRBS_EN),
			ID(P_LANE_POWERDOWN),
			ID(P_LANE_RST),
			ID(P_RX_LANE_POWERDOWN),
			ID(P_RX_PMA_RST),
			ID(P_RX_CDR_RST),
			ID(P_RX_CLKPATH_RST),
			ID(P_RX_DFE_RST),
			ID(P_RX_LEQ_RST),
			ID(P_RX_SLIDING_RST),
			ID(P_RX_DFE_EN),
			ID(P_RX_T1_EN),
			ID(P_RX_CDRX_EN),
			ID(P_RX_T1_DFE_EN),
			ID(P_RX_T2_DFE_EN),
			ID(P_RX_T3_DFE_EN),
			ID(P_RX_T4_DFE_EN),
			ID(P_RX_T5_DFE_EN),
			ID(P_RX_T6_DFE_EN),
			ID(P_RX_SLIDING_EN),
			ID(P_RX_EYE_RST),
			ID(P_RX_EYE_EN),
			ID(P_RX_EYE_TAP),
			ID(P_RX_PIC_EYE),
			ID(P_RX_PIC_FASTLOCK),
			ID(P_RX_PIC_FASTLOCK_STROBE),
			ID(P_EM_RD_TRIGGER),
			ID(P_EM_MODE_CTRL),
			ID(P_RX_CTLE_DCCAL_RST),
			ID(P_RX_SLICER_DCCAL_RST),
			ID(P_RX_SLICER_DCCAL_EN),
			ID(P_RX_CTLE_DCCAL_EN),
			ID(P_RX_SLIP_RST),
			ID(P_RX_SLIP_EN),
			ID(P_LPLL_POWERDOWN),
			ID(P_LPLL_RST),
			ID(P_LPLL_LOCKDET_RST),
			ID(P_TX_LS_DATA),
			ID(P_TX_BEACON_EN),
			ID(P_TX_SWING),
			ID(P_TX_RXDET_REQ),
			ID(P_TX_RATE),
			ID(P_TX_BUSWIDTH),
			ID(P_TX_FREERUN_BUSWIDTH),
			ID(P_TX_MARGIN),
			ID(P_TX_PMA_RST),
			ID(P_TX_LANE_POWERDOWN),
			ID(P_TX_PIC_EN),
			ID(P_RX_RATE),
			ID(P_RX_BUSWIDTH),
			ID(P_RX_HIGHZ),
			ID(P_CIM_CLK_ALIGNER_RX),
			ID(P_CIM_CLK_ALIGNER_TX),
			ID(P_ALIGN_MODE_VALID_RX),
			ID(P_ALIGN_MODE_RX),
			ID(P_ALIGN_MODE_VALID_TX),
			ID(P_ALIGN_MODE_TX),
			ID(PMA_HPLL_CK0),
			ID(PMA_HPLL_CK90),
			ID(PMA_HPLL_CK180),
			ID(PMA_HPLL_CK270),
			ID(PMA_HPLL_READY_IN),
			ID(PMA_HPLL_REFCLK_IN),
			ID(PMA_TX_SYNC_HPLL_IN),
			ID(P_LPLL_REFCLK_IN),
			ID(P_TX_RATE_CHANGE_ON_0),
			ID(P_TX_RATE_CHANGE_ON_1),
			ID(P_TX_SYNC),
			ID(P_RX_SDN),
			ID(P_RX_SDP),
			ID(P_TEST_SE_N),
			ID(P_TEST_MODE_N),
			ID(P_TEST_RSTN),
			ID(P_TEST_SI0),
			ID(P_TEST_SI1),
			ID(P_TEST_SI2),
			ID(P_TEST_SI3),
			ID(P_TEST_SI4),
			ID(P_FOR_PMA_TEST_MODE_N),
			ID(P_FOR_PMA_TEST_SE_N),
			ID(P_FOR_PMA_TEST_CLK),
			ID(P_FOR_PMA_TEST_RSTN),
			ID(P_FOR_PMA_TEST_SI)},
		       {ID(P_CFG_READY),
			ID(P_CFG_RDATA),
			ID(P_CFG_INT),
			ID(LANE_COUT_BUS_FORWARD),
			ID(LANE_COUT_BUS_BACKWARD),
			ID(P_RX_PRBS_ERROR),
			ID(P_PCS_RX_MCB_STATUS),
			ID(P_PCS_LSM_SYNCED),
			ID(P_RDATA),
			ID(P_RXDVLD),
			ID(P_RXDVLD_H),
			ID(P_RXSTATUS),
			ID(P_EM_ERROR_CNT),
			ID(P_LPLL_READY),
			ID(P_RX_SIGDET_STATUS),
			ID(P_RX_SATA_COMINIT),
			ID(P_RX_SATA_COMWAKE),
			ID(P_RX_LS_DATA),
			ID(P_RX_READY),
			ID(P_TEST_STATUS),
			ID(P_TX_RXDET_STATUS),
			ID(P_RCLK2FABRIC),
			ID(P_TCLK2FABRIC),
			ID(P_CA_ALIGN_RX),
			ID(P_CA_ALIGN_TX),
			ID(P_TX_SDN),
			ID(P_TX_SDP),
			ID(P_TEST_SO0),
			ID(P_TEST_SO1),
			ID(P_TEST_SO2),
			ID(P_TEST_SO3),
			ID(P_TEST_SO4),
			ID(P_FOR_PMA_TEST_SO)},
		       false);
	ct->setup_type(ID(GTP_HSSTHP_LANE_E1),
		       {ID(P_RX_CLK_FR_CORE),
			ID(P_RCLK2_FR_CORE),
			ID(P_TX_CLK_FR_CORE),
			ID(P_TCLK2_FR_CORE),
			ID(P_PCS_RX_RST),
			ID(P_PCS_TX_RST),
			ID(P_EXT_BRIDGE_PCS_RST),
			ID(P_CFG_RST),
			ID(P_CFG_CLK),
			ID(P_CFG_PSEL),
			ID(P_CFG_ENABLE),
			ID(P_CFG_WRITE),
			ID(P_CFG_ADDR),
			ID(P_CFG_WDATA),
			ID(LANE_CIN_BUS_FORWARD),
			ID(LANE_CIN_BUS_BACKWARD),
			ID(P_TDATA),
			ID(P_PCIE_EI_H),
			ID(P_PCIE_EI_L),
			ID(P_TX_DEEMP),
			ID(P_TX_DEEMP_POST_SEL),
			ID(P_BLK_ALIGN_CTRL),
			ID(P_TX_ENC_TYPE),
			ID(P_RX_DEC_TYPE),
			ID(P_PCS_BIT_SLIP),
			ID(P_PCS_WORD_ALIGN_EN),
			ID(P_RX_POLARITY_INVERT),
			ID(P_PCS_MCB_EXT_EN),
			ID(P_PCS_NEAREND_LOOP),
			ID(P_PCS_FAREND_LOOP),
			ID(P_PMA_NEAREND_PLOOP),
			ID(P_PMA_NEAREND_SLOOP),
			ID(P_PMA_FAREND_PLOOP),
			ID(P_PCS_PRBS_EN),
			ID(P_LANE_POWERDOWN),
			ID(P_LANE_RST),
			ID(P_RX_LANE_POWERDOWN),
			ID(P_RX_PMA_RST),
			ID(P_RX_CDR_RST),
			ID(P_RX_CLKPATH_RST),
			ID(P_RX_DFE_RST),
			ID(P_RX_LEQ_RST),
			ID(P_RX_SLIDING_RST),
			ID(P_RX_DFE_EN),
			ID(P_RX_T1_EN),
			ID(P_RX_CDRX_EN),
			ID(P_RX_T1_DFE_EN),
			ID(P_RX_T2_DFE_EN),
			ID(P_RX_T3_DFE_EN),
			ID(P_RX_T4_DFE_EN),
			ID(P_RX_T5_DFE_EN),
			ID(P_RX_T6_DFE_EN),
			ID(P_RX_SLIDING_EN),
			ID(P_RX_EYE_RST),
			ID(P_RX_EYE_EN),
			ID(P_RX_EYE_TAP),
			ID(P_RX_PIC_EYE),
			ID(P_RX_PIC_FASTLOCK),
			ID(P_RX_PIC_FASTLOCK_STROBE),
			ID(P_EM_RD_TRIGGER),
			ID(P_EM_MODE_CTRL),
			ID(P_RX_CTLE_DCCAL_RST),
			ID(P_RX_SLICER_DCCAL_RST),
			ID(P_RX_SLICER_DCCAL_EN),
			ID(P_RX_CTLE_DCCAL_EN),
			ID(P_RX_SLIP_RST),
			ID(P_RX_SLIP_EN),
			ID(P_LPLL_POWERDOWN),
			ID(P_LPLL_RST),
			ID(P_LPLL_LOCKDET_RST),
			ID(P_TX_LS_DATA),
			ID(P_TX_BEACON_EN),
			ID(P_TX_SWING),
			ID(P_TX_RXDET_REQ),
			ID(P_TX_RATE),
			ID(P_TX_BUSWIDTH),
			ID(P_TX_FREERUN_BUSWIDTH),
			ID(P_TX_MARGIN),
			ID(P_TX_PMA_RST),
			ID(P_TX_LANE_POWERDOWN),
			ID(P_TX_PIC_EN),
			ID(P_RX_RATE),
			ID(P_RX_BUSWIDTH),
			ID(P_RX_HIGHZ),
			ID(P_CIM_CLK_ALIGNER_RX),
			ID(P_CIM_CLK_ALIGNER_TX),
			ID(P_ALIGN_MODE_VALID_RX),
			ID(P_ALIGN_MODE_RX),
			ID(P_ALIGN_MODE_VALID_TX),
			ID(P_ALIGN_MODE_TX),
			ID(PMA_HPLL_CK0),
			ID(PMA_HPLL_CK90),
			ID(PMA_HPLL_CK180),
			ID(PMA_HPLL_CK270),
			ID(PMA_HPLL_READY_IN),
			ID(PMA_HPLL_REFCLK_IN),
			ID(PMA_TX_SYNC_HPLL_IN),
			ID(P_LPLL_REFCLK_IN),
			ID(P_TX_RATE_CHANGE_ON_0),
			ID(P_TX_RATE_CHANGE_ON_1),
			ID(P_TX_SYNC),
			ID(P_RX_SDN),
			ID(P_RX_SDP),
			ID(P_CDR_PICTRL_OW),
			ID(P_CDR_PICTRL_OW_VAL)},
		       {ID(P_CFG_READY),
			ID(P_CFG_RDATA),
			ID(P_CFG_INT),
			ID(LANE_COUT_BUS_FORWARD),
			ID(LANE_COUT_BUS_BACKWARD),
			ID(P_RX_PRBS_ERROR),
			ID(P_PCS_RX_MCB_STATUS),
			ID(P_PCS_LSM_SYNCED),
			ID(P_RDATA),
			ID(P_RXDVLD),
			ID(P_RXDVLD_H),
			ID(P_RXSTATUS),
			ID(P_EM_ERROR_CNT),
			ID(P_LPLL_READY),
			ID(P_RX_SIGDET_STATUS),
			ID(P_RX_SATA_COMINIT),
			ID(P_RX_SATA_COMWAKE),
			ID(P_RX_LS_DATA),
			ID(P_RX_READY),
			ID(P_TEST_STATUS),
			ID(P_TX_RXDET_STATUS),
			ID(P_RCLK2FABRIC),
			ID(P_TCLK2FABRIC),
			ID(P_CA_ALIGN_RX),
			ID(P_CA_ALIGN_TX),
			ID(P_TX_SDN),
			ID(P_TX_SDP)},
		       false);
	ct->setup_type(ID(GTP_HSSTHP_LANE_t),
		       {ID(P_RX_CLK_FR_CORE),
			ID(P_RCLK2_FR_CORE),
			ID(P_TX_CLK_FR_CORE),
			ID(P_TCLK2_FR_CORE),
			ID(P_PCS_RX_RST),
			ID(P_PCS_TX_RST),
			ID(P_EXT_BRIDGE_PCS_RST),
			ID(P_CFG_RST),
			ID(P_CFG_CLK),
			ID(P_CFG_PSEL),
			ID(P_CFG_ENABLE),
			ID(P_CFG_WRITE),
			ID(P_CFG_ADDR),
			ID(P_CFG_WDATA),
			ID(LANE_CIN_BUS_FORWARD),
			ID(LANE_CIN_BUS_BACKWARD),
			ID(P_TDATA),
			ID(P_PCIE_EI_H),
			ID(P_PCIE_EI_L),
			ID(P_TX_DEEMP),
			ID(P_TX_DEEMP_POST_SEL),
			ID(P_BLK_ALIGN_CTRL),
			ID(P_TX_ENC_TYPE),
			ID(P_RX_DEC_TYPE),
			ID(P_PCS_BIT_SLIP),
			ID(P_PCS_WORD_ALIGN_EN),
			ID(P_RX_POLARITY_INVERT),
			ID(P_PCS_MCB_EXT_EN),
			ID(P_PCS_NEAREND_LOOP),
			ID(P_PCS_FAREND_LOOP),
			ID(P_PMA_NEAREND_PLOOP),
			ID(P_PMA_NEAREND_SLOOP),
			ID(P_PMA_FAREND_PLOOP),
			ID(P_PCS_PRBS_EN),
			ID(P_LANE_POWERDOWN),
			ID(P_LANE_RST),
			ID(P_RX_LANE_POWERDOWN),
			ID(P_RX_PMA_RST),
			ID(P_RX_CDR_RST),
			ID(P_RX_CLKPATH_RST),
			ID(P_RX_DFE_RST),
			ID(P_RX_LEQ_RST),
			ID(P_RX_SLIDING_RST),
			ID(P_RX_DFE_EN),
			ID(P_RX_T1_EN),
			ID(P_RX_CDRX_EN),
			ID(P_RX_T1_DFE_EN),
			ID(P_RX_T2_DFE_EN),
			ID(P_RX_T3_DFE_EN),
			ID(P_RX_T4_DFE_EN),
			ID(P_RX_T5_DFE_EN),
			ID(P_RX_T6_DFE_EN),
			ID(P_RX_SLIDING_EN),
			ID(P_RX_EYE_RST),
			ID(P_RX_EYE_EN),
			ID(P_RX_EYE_TAP),
			ID(P_RX_PIC_EYE),
			ID(P_RX_PIC_FASTLOCK),
			ID(P_RX_PIC_FASTLOCK_STROBE),
			ID(P_EM_RD_TRIGGER),
			ID(P_EM_MODE_CTRL),
			ID(P_RX_CTLE_DCCAL_RST),
			ID(P_RX_SLICER_DCCAL_RST),
			ID(P_RX_SLICER_DCCAL_EN),
			ID(P_RX_CTLE_DCCAL_EN),
			ID(P_RX_SLIP_RST),
			ID(P_RX_SLIP_EN),
			ID(P_LPLL_POWERDOWN),
			ID(P_LPLL_RST),
			ID(P_LPLL_LOCKDET_RST),
			ID(P_TX_LS_DATA),
			ID(P_TX_BEACON_EN),
			ID(P_TX_SWING),
			ID(P_TX_RXDET_REQ),
			ID(P_TX_RATE),
			ID(P_TX_BUSWIDTH),
			ID(P_TX_FREERUN_BUSWIDTH),
			ID(P_TX_MARGIN),
			ID(P_TX_PMA_RST),
			ID(P_TX_LANE_POWERDOWN),
			ID(P_TX_PIC_EN),
			ID(P_RX_RATE),
			ID(P_RX_BUSWIDTH),
			ID(P_RX_HIGHZ),
			ID(P_CIM_CLK_ALIGNER_RX),
			ID(P_CIM_CLK_ALIGNER_TX),
			ID(P_ALIGN_MODE_VALID_RX),
			ID(P_ALIGN_MODE_RX),
			ID(P_ALIGN_MODE_VALID_TX),
			ID(P_ALIGN_MODE_TX),
			ID(PMA_HPLL_CK0),
			ID(PMA_HPLL_CK90),
			ID(PMA_HPLL_CK180),
			ID(PMA_HPLL_CK270),
			ID(PMA_HPLL_READY_IN),
			ID(PMA_HPLL_REFCLK_IN),
			ID(PMA_IPN50U_IN),
			ID(PMA_LPLL_REFCLK),
			ID(PMA_RES_CAL_I),
			ID(PMA_TX_SYNC_HPLL),
			ID(PMA_TX_RATE_CHANGE_ON_0),
			ID(PMA_TX_RATE_CHANGE_ON_1),
			ID(PMA_TX_SYNC),
			ID(P_RX_SDN),
			ID(P_RX_SDP)},
		       {ID(P_CFG_READY),
			ID(P_CFG_RDATA),
			ID(P_CFG_INT),
			ID(LANE_COUT_BUS_FORWARD),
			ID(LANE_COUT_BUS_BACKWARD),
			ID(P_RX_PRBS_ERROR),
			ID(P_PCS_RX_MCB_STATUS),
			ID(P_PCS_LSM_SYNCED),
			ID(P_RDATA),
			ID(P_RXDVLD),
			ID(P_RXDVLD_H),
			ID(P_RXSTATUS),
			ID(P_EM_ERROR_CNT),
			ID(P_LPLL_READY),
			ID(P_RX_SIGDET_STATUS),
			ID(P_RX_SATA_COMINIT),
			ID(P_RX_SATA_COMWAKE),
			ID(P_RX_LS_DATA),
			ID(P_RX_READY),
			ID(P_TEST_STATUS),
			ID(P_TX_RXDET_STATUS),
			ID(P_RCLK2FABRIC),
			ID(P_TCLK2FABRIC),
			ID(P_CA_ALIGN_RX),
			ID(P_CA_ALIGN_TX),
			ID(P_TX_SDN),
			ID(P_TX_SDP)},
		       false);
	ct->setup_type(ID(GTP_IDDR_E1), {ID(D), ID(CLK), ID(CE), ID(RS)}, {ID(Q0), ID(Q1)}, false);
	ct->setup_type(ID(GTP_INBUF), {ID(I)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_INBUFDS), {ID(I), ID(IB)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_INBUFDS_E1), {ID(I), ID(IB)}, {ID(O), ID(OB)}, false);
	ct->setup_type(ID(GTP_INBUFE), {ID(EN), ID(I)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_INBUFEDS), {ID(EN), ID(I), ID(IB)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_INBUFEDS_E1), {ID(EN), ID(I), ID(IB)}, {ID(O), ID(OB)}, false);
	ct->setup_type(ID(GTP_INBUFG), {ID(I)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_INBUFGDS), {ID(I), ID(IB)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_INV), {ID(I)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_IOBUF), {ID(I), ID(T)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_IOBUF_RX_MIPI), {ID(I_LP), ID(IB_LP), ID(T), ID(TB), ID(M)}, {ID(O_LP), ID(OB_LP), ID(O_HS)}, false);
	ct->setup_type(ID(GTP_IOBUF_TX_MIPI), {ID(I_LP), ID(IB_LP), ID(I_HS), ID(T), ID(TB), ID(M)}, {ID(O_LP), ID(OB_LP)}, false);
	ct->setup_type(ID(GTP_IOBUFCO), {ID(I), ID(T)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_IOBUFCO_E1), {ID(I), ID(IB), ID(T), ID(TB)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_IOBUFDS), {ID(I), ID(T)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_IOBUFE), {ID(I), ID(EN), ID(T)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_IOBUFECO), {ID(I), ID(EN), ID(T)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_IOBUFEDS), {ID(I), ID(EN), ID(T)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_IOCLKBUF), {ID(CLKIN), ID(DI)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_IOCLKDIV_E2), {ID(CLKIN), ID(RST_N), ID(CE)}, {ID(CLKDIVOUT)}, false);
	ct->setup_type(ID(GTP_IOCLKDIV_E3), {ID(RST), ID(CLKIN)}, {ID(CLKDIVOUT)}, false);
	ct->setup_type(ID(GTP_IODELAY_E2), {ID(DI), ID(DELAY_SEL), ID(DELAY_STEP), ID(EN_N)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_IPAL_E2), {ID(RST_N), ID(CLK), ID(CS_N), ID(RW_SEL), ID(DI)},
		       {ID(DO), ID(RBCRC_ERR), ID(RBCRC_VALID), ID(ECC_VALID), ID(ECC_INDEX), ID(SERROR), ID(DERROR), ID(SEU_FRAME_ADDR),
			ID(SEU_COLUMN_ADDR), ID(SEU_REGION_ADDR), ID(SEU_FRAME_NADDR), ID(SEU_COLUMN_NADDR), ID(SEU_REGION_NADDR), ID(PRCFG_OVER),
			ID(PRCFG_ERR), ID(DRCFG_OVER), ID(DRCFG_ERR)},
		       false);
	ct->setup_type(ID(GTP_ISERDES_E2),
		       {ID(RST), ID(ICE0), ID(ICE1), ID(DESCLK), ID(ICLK), ID(ICLKB), ID(OCLK), ID(ICLKDIV), ID(DI), ID(BITSLIP), ID(ISHIFTIN0),
			ID(ISHIFTIN1), ID(IFIFO_WADDR), ID(IFIFO_RADDR)},
		       {ID(DO), ID(ISHIFTOUT0), ID(ISHIFTOUT1)}, false);
	ct->setup_type(ID(GTP_ISERDES_FIFO), {ID(DIN), ID(VALID_I), ID(RST), ID(EN), ID(WCLK), ID(RCLK)},
		       {ID(DOUT), ID(VALID_O), ID(EMPTY), ID(FULL)}, false);
	ct->setup_type(ID(GTP_JTAGIF), {ID(TCK), ID(TMS), ID(TDI)}, {ID(TDO)}, false);
	ct->setup_type(ID(GTP_KEYRAM), {ID(ERASE_KEY_N)}, {}, false);
	ct->setup_type(ID(GTP_LUT1), {ID(I0)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_LUT2), {ID(I0), ID(I1)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_LUT3), {ID(I0), ID(I1), ID(I2)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_LUT4), {ID(I0), ID(I1), ID(I2), ID(I3)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_LUT5), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_LUT6), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_LUT6CARRY), {ID(CIN), ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5)}, {ID(COUT), ID(Z)}, false);
	ct->setup_type(ID(GTP_LUT6D), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5)}, {ID(Z), ID(Z5)}, false);
	ct->setup_type(ID(GTP_LUT7), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5), ID(I6)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_LUT8), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5), ID(I6), ID(I7)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_MONITOR_E1), {ID(RST_N), ID(CLK), ID(EN), ID(SAMPLE), ID(SEL)}, {ID(DATA), ID(DATA_VALID), ID(READY)}, false);
	ct->setup_type(ID(GTP_MONITOR_E1_DFT),
		       {ID(RST_N), ID(CLK), ID(EN), ID(SAMPLE), ID(SEL), ID(TEST_BIST_MODE), ID(TEST_BIST_START), ID(TEST_TDI), ID(TEST_TCK),
			ID(TEST_FLG_JTAG), ID(TEST_SHIFTDR), ID(TEST_CLOCKDR), ID(TEST_CAPTUREDR), ID(TEST_UPDATEDR), ID(TEST_BIST_REF_MAX),
			ID(TEST_BIST_REF_MIN)},
		       {ID(DATA), ID(DATA_VALID), ID(READY), ID(TEST_TDO), ID(TEST_BIST_VALID), ID(TEST_BIST_ERROR)}, false);
	ct->setup_type(ID(GTP_MUX2LUT7), {ID(I0), ID(I1), ID(S)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_MUX2LUT8), {ID(I0), ID(I1), ID(S)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_ODDR_E1), {ID(D0), ID(D1), ID(CLK), ID(CE), ID(RS)}, {ID(Q)}, false);
	ct->setup_type(ID(GTP_ONE), {}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_OSC_E4), {ID(EN_N)}, {ID(CLKOUT)}, false);
	ct->setup_type(ID(GTP_OSERDES_E2),
		       {ID(RST), ID(OCE), ID(TCE), ID(OCLKDIV), ID(SERCLK), ID(OCLK), ID(MIPI_CTRL), ID(UPD0_SHIFT), ID(UPD1_SHIFT), ID(OSHIFTIN0),
			ID(OSHIFTIN1), ID(DI), ID(TI), ID(TBYTE_IN)},
		       {ID(OSHIFTOUT0), ID(OSHIFTOUT1), ID(TQ), ID(DO), ID(TFB), ID(TERM_FB)}, false);
	ct->setup_type(ID(GTP_OSERDES_FIFO), {ID(DIN), ID(TIN), ID(RST), ID(EN), ID(WCLK), ID(RCLK)}, {ID(DOUT), ID(TOUT), ID(EMPTY), ID(FULL)},
		       false);
	ct->setup_type(ID(GTP_OUTBUF), {ID(I)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_OUTBUFCO), {ID(I)}, {ID(O), ID(OB)}, false);
	ct->setup_type(ID(GTP_OUTBUFCO_E1), {ID(I), ID(IB)}, {ID(O), ID(OB)}, false);
	ct->setup_type(ID(GTP_OUTBUFDS), {ID(I)}, {ID(O), ID(OB)}, false);
	ct->setup_type(ID(GTP_OUTBUFT), {ID(I), ID(T)}, {ID(O)}, false);
	ct->setup_type(ID(GTP_OUTBUFTCO), {ID(I), ID(T)}, {ID(O), ID(OB)}, false);
	ct->setup_type(ID(GTP_OUTBUFTCO_E1), {ID(I), ID(IB), ID(T), ID(TB)}, {ID(O), ID(OB)}, false);
	ct->setup_type(ID(GTP_OUTBUFTDS), {ID(I), ID(T)}, {ID(O), ID(OB)}, false);
	ct->setup_type(ID(GTP_PCIEGEN2),
		       {ID(MEM_CLK),
			ID(PCLK),
			ID(PCLK_DIV2),
			ID(BUTTON_RST),
			ID(POWER_UP_RST),
			ID(PERST),
			ID(APP_INIT_RST),
			ID(DEVICE_TYPE),
			ID(RX_LANE_FLIP_EN),
			ID(TX_LANE_FLIP_EN),
			ID(APP_LTSSM_EN),
			ID(APP_REQ_RETRY_EN),
			ID(AXIS_MASTER_TREADY),
			ID(TRGT1_RADM_PKT_HALT),
			ID(AXIS_SLAVE0_TVALID),
			ID(AXIS_SLAVE0_TDATA),
			ID(AXIS_SLAVE0_TLAST),
			ID(AXIS_SLAVE0_TUSER),
			ID(AXIS_SLAVE1_TVALID),
			ID(AXIS_SLAVE1_TDATA),
			ID(AXIS_SLAVE1_TLAST),
			ID(AXIS_SLAVE1_TUSER),
			ID(AXIS_SLAVE2_TVALID),
			ID(AXIS_SLAVE2_TDATA),
			ID(AXIS_SLAVE2_TLAST),
			ID(AXIS_SLAVE2_TUSER),
			ID(DBI_ADDR),
			ID(DBI_DIN),
			ID(DBI_CS),
			ID(DBI_CS2),
			ID(DBI_WR),
			ID(APP_DBI_RO_WR_DISABLE),
			ID(SEDI),
			ID(SEDI_ACK),
			ID(SYS_INT),
			ID(VEN_MSI_REQ),
			ID(VEN_MSI_TC),
			ID(VEN_MSI_VECTOR),
			ID(CFG_MSI_PENDING),
			ID(MSIX_ADDR),
			ID(MSIX_DATA),
			ID(OUTBAND_PWRUP_CMD),
			ID(APP_REQ_ENTR_L1),
			ID(APP_READY_ENTR_L23),
			ID(APP_REQ_EXIT_L1),
			ID(APP_XFER_PENDING),
			ID(APPS_PM_XMT_TURNOFF),
			ID(APP_UNLOCK_MSG),
			ID(APPS_PM_XMT_PME),
			ID(APP_CLK_PM_EN),
			ID(SYS_AUX_PWR_DET),
			ID(APP_HDR_VALID),
			ID(APP_HDR_LOG),
			ID(APP_ERR_BUS),
			ID(APP_ERR_ADVISORY),
			ID(DIAG_CTRL_BUS),
			ID(DYN_DEBUG_INFO_SEL),
			ID(APP_RAS_DES_SD_HOLD_LTSSM),
			ID(APP_RAS_DES_TBA_CTRL),
			ID(PHY_MAC_RXELECIDLE),
			ID(PHY_MAC_PHYSTATUS),
			ID(PHY_MAC_RXDATA),
			ID(PHY_MAC_RXDATAK),
			ID(PHY_MAC_RXVALID),
			ID(PHY_MAC_RXSTATUS),
			ID(P_DATAQ_DATAOUT),
			ID(RETRYRAM_XDLH_DATA),
			ID(P_HDRQ_DATAOUT),
			ID(RAM_TEST_EN),
			ID(RAM_TEST_ADDRH),
			ID(RETRY_TEST_DATA_EN),
			ID(RAM_TEST_MODE_N)},
		       {ID(CORE_RST_N),
			ID(TRAINING_RST_N),
			ID(PHY_RST_N),
			ID(SMLH_LINK_UP),
			ID(RDLH_LINK_UP),
			ID(SMLH_LTSSM_STATE),
			ID(AXIS_MASTER_TVALID),
			ID(AXIS_MASTER_TDATA),
			ID(AXIS_MASTER_TKEEP),
			ID(AXIS_MASTER_TLAST),
			ID(AXIS_MASTER_TUSER),
			ID(RADM_GRANT_TLP_TYPE),
			ID(AXIS_SLAVE0_TREADY),
			ID(AXIS_SLAVE1_TREADY),
			ID(AXIS_SLAVE2_TREADY),
			ID(PM_XTLH_BLOCK_TLP),
			ID(LBC_DBI_ACK),
			ID(LBC_DBI_DOUT),
			ID(SEDO),
			ID(SEDO_EN),
			ID(CFG_INT_DISABLE),
			ID(INTA_GRT_MUX),
			ID(INTB_GRT_MUX),
			ID(INTC_GRT_MUX),
			ID(INTD_GRT_MUX),
			ID(VEN_MSI_GRANT),
			ID(CFG_MSI_EN),
			ID(CFG_MSIX_EN),
			ID(CFG_MSIX_FUNC_MASK),
			ID(RADM_PM_TURNOFF),
			ID(RADM_MSG_UNLOCK),
			ID(PM_STATUS),
			ID(PM_DSTATE),
			ID(AUX_PM_EN),
			ID(PM_PME_EN),
			ID(PM_LINKST_IN_L0S),
			ID(PM_LINKST_IN_L1),
			ID(PM_LINKST_IN_L2),
			ID(PM_LINKST_L2_EXIT),
			ID(WAKE),
			ID(RADM_PM_PME),
			ID(RADM_PM_TO_ACK),
			ID(PM_MASTER_STATE),
			ID(PM_SLAVE_STATE),
			ID(CFG_SEND_COR_ERR_MUX),
			ID(CFG_SEND_NF_ERR_MUX),
			ID(CFG_SEND_F_ERR_MUX),
			ID(CFG_SYS_ERR_RC),
			ID(CFG_AER_RC_ERR_MUX),
			ID(RADM_CPL_TIMEOUT),
			ID(RADM_TIMEOUT_CPL_TC),
			ID(RADM_TIMEOUT_CPL_TAG),
			ID(RADM_TIMEOUT_CPL_ATTR),
			ID(RADM_TIMEOUT_CPL_LEN),
			ID(CFG_MAX_RD_REQ_SIZE),
			ID(CFG_BUS_MASTER_EN),
			ID(CFG_MAX_PAYLOAD_SIZE),
			ID(CFG_RCB),
			ID(CFG_MEM_SPACE_EN),
			ID(CFG_PM_NO_SOFT_RST),
			ID(CFG_CRS_SW_VIS_EN),
			ID(CFG_NO_SNOOP_EN),
			ID(CFG_RELAX_ORDER_EN),
			ID(CFG_TPH_REQ_EN),
			ID(CFG_PF_TPH_ST_MODE),
			ID(CFG_PBUS_NUM),
			ID(CFG_PBUS_DEV_NUM),
			ID(RBAR_CTRL_UPDATE),
			ID(CFG_ATOMIC_REQ_EN),
			ID(CFG_ATOMIC_EGRESS_BLOCK),
			ID(CFG_EXT_TAG_EN),
			ID(RADM_IDLE),
			ID(RADM_Q_NOT_EMPTY),
			ID(RADM_QOVERFLOW),
			ID(CFG_LINK_AUTO_BW_MUX),
			ID(CFG_BW_MGT_MUX),
			ID(CFG_PME_MUX),
			ID(DEBUG_INFO_MUX),
			ID(CFG_IDO_REQ_EN),
			ID(CFG_IDO_CPL_EN),
			ID(XADM_PH_CDTS),
			ID(XADM_PD_CDTS),
			ID(XADM_NPH_CDTS),
			ID(XADM_NPD_CDTS),
			ID(XADM_CPLH_CDTS),
			ID(XADM_CPLD_CDTS),
			ID(MAC_PHY_POWERDOWN),
			ID(MAC_PHY_TXDATA),
			ID(MAC_PHY_TXDATAK),
			ID(MAC_PHY_TXDETECTRX_LOOPBACK),
			ID(MAC_PHY_TXELECIDLE_L),
			ID(MAC_PHY_TXELECIDLE_H),
			ID(MAC_PHY_TXCOMPLIANCE),
			ID(MAC_PHY_RXPOLARITY),
			ID(MAC_PHY_RATE),
			ID(MAC_PHY_TXDEEMPH),
			ID(MAC_PHY_TXMARGIN),
			ID(MAC_PHY_TXSWING),
			ID(CFG_HW_AUTO_SP_DIS),
			ID(P_DATAQ_ADDRA),
			ID(P_DATAQ_ADDRB),
			ID(P_DATAQ_DATAIN),
			ID(P_DATAQ_ENA),
			ID(P_DATAQ_ENB),
			ID(P_DATAQ_WEA),
			ID(XDLH_RETRYRAM_ADDR),
			ID(XDLH_RETRYRAM_DATA),
			ID(XDLH_RETRYRAM_WE),
			ID(XDLH_RETRYRAM_EN),
			ID(P_HDRQ_ADDRA),
			ID(P_HDRQ_ADDRB),
			ID(P_HDRQ_DATAIN),
			ID(P_HDRQ_ENA),
			ID(P_HDRQ_ENB),
			ID(P_HDRQ_WEA)},
		       false);
	ct->setup_type(ID(GTP_PCIEGEN2_DFT),
		       {ID(MEM_CLK),
			ID(PCLK),
			ID(PCLK_DIV2),
			ID(BUTTON_RST),
			ID(POWER_UP_RST),
			ID(PERST),
			ID(APP_INIT_RST),
			ID(DEVICE_TYPE),
			ID(RX_LANE_FLIP_EN),
			ID(TX_LANE_FLIP_EN),
			ID(APP_LTSSM_EN),
			ID(APP_REQ_RETRY_EN),
			ID(AXIS_MASTER_TREADY),
			ID(TRGT1_RADM_PKT_HALT),
			ID(AXIS_SLAVE0_TVALID),
			ID(AXIS_SLAVE0_TDATA),
			ID(AXIS_SLAVE0_TLAST),
			ID(AXIS_SLAVE0_TUSER),
			ID(AXIS_SLAVE1_TVALID),
			ID(AXIS_SLAVE1_TDATA),
			ID(AXIS_SLAVE1_TLAST),
			ID(AXIS_SLAVE1_TUSER),
			ID(AXIS_SLAVE2_TVALID),
			ID(AXIS_SLAVE2_TDATA),
			ID(AXIS_SLAVE2_TLAST),
			ID(AXIS_SLAVE2_TUSER),
			ID(DBI_ADDR),
			ID(DBI_DIN),
			ID(DBI_CS),
			ID(DBI_CS2),
			ID(DBI_WR),
			ID(APP_DBI_RO_WR_DISABLE),
			ID(SEDI),
			ID(SEDI_ACK),
			ID(SYS_INT),
			ID(VEN_MSI_REQ),
			ID(VEN_MSI_TC),
			ID(VEN_MSI_VECTOR),
			ID(CFG_MSI_PENDING),
			ID(MSIX_ADDR),
			ID(MSIX_DATA),
			ID(OUTBAND_PWRUP_CMD),
			ID(APP_REQ_ENTR_L1),
			ID(APP_READY_ENTR_L23),
			ID(APP_REQ_EXIT_L1),
			ID(APP_XFER_PENDING),
			ID(APPS_PM_XMT_TURNOFF),
			ID(APP_UNLOCK_MSG),
			ID(APPS_PM_XMT_PME),
			ID(APP_CLK_PM_EN),
			ID(SYS_AUX_PWR_DET),
			ID(APP_HDR_VALID),
			ID(APP_HDR_LOG),
			ID(APP_ERR_BUS),
			ID(APP_ERR_ADVISORY),
			ID(DIAG_CTRL_BUS),
			ID(DYN_DEBUG_INFO_SEL),
			ID(APP_RAS_DES_SD_HOLD_LTSSM),
			ID(APP_RAS_DES_TBA_CTRL),
			ID(PHY_MAC_RXELECIDLE),
			ID(PHY_MAC_PHYSTATUS),
			ID(PHY_MAC_RXDATA),
			ID(PHY_MAC_RXDATAK),
			ID(PHY_MAC_RXVALID),
			ID(PHY_MAC_RXSTATUS),
			ID(P_DATAQ_DATAOUT),
			ID(RETRYRAM_XDLH_DATA),
			ID(P_HDRQ_DATAOUT),
			ID(RAM_TEST_EN),
			ID(RAM_TEST_ADDRH),
			ID(RETRY_TEST_DATA_EN),
			ID(RAM_TEST_MODE_N),
			ID(TEST_SE_N),
			ID(TEST_RST_N),
			ID(TEST_MODE_N)},
		       {ID(CORE_RST_N),
			ID(TRAINING_RST_N),
			ID(PHY_RST_N),
			ID(SMLH_LINK_UP),
			ID(RDLH_LINK_UP),
			ID(SMLH_LTSSM_STATE),
			ID(AXIS_MASTER_TVALID),
			ID(AXIS_MASTER_TDATA),
			ID(AXIS_MASTER_TKEEP),
			ID(AXIS_MASTER_TLAST),
			ID(AXIS_MASTER_TUSER),
			ID(RADM_GRANT_TLP_TYPE),
			ID(AXIS_SLAVE0_TREADY),
			ID(AXIS_SLAVE1_TREADY),
			ID(AXIS_SLAVE2_TREADY),
			ID(PM_XTLH_BLOCK_TLP),
			ID(LBC_DBI_ACK),
			ID(LBC_DBI_DOUT),
			ID(SEDO),
			ID(SEDO_EN),
			ID(CFG_INT_DISABLE),
			ID(INTA_GRT_MUX),
			ID(INTB_GRT_MUX),
			ID(INTC_GRT_MUX),
			ID(INTD_GRT_MUX),
			ID(VEN_MSI_GRANT),
			ID(CFG_MSI_EN),
			ID(CFG_MSIX_EN),
			ID(CFG_MSIX_FUNC_MASK),
			ID(RADM_PM_TURNOFF),
			ID(RADM_MSG_UNLOCK),
			ID(PM_STATUS),
			ID(PM_DSTATE),
			ID(AUX_PM_EN),
			ID(PM_PME_EN),
			ID(PM_LINKST_IN_L0S),
			ID(PM_LINKST_IN_L1),
			ID(PM_LINKST_IN_L2),
			ID(PM_LINKST_L2_EXIT),
			ID(WAKE),
			ID(RADM_PM_PME),
			ID(RADM_PM_TO_ACK),
			ID(PM_MASTER_STATE),
			ID(PM_SLAVE_STATE),
			ID(CFG_SEND_COR_ERR_MUX),
			ID(CFG_SEND_NF_ERR_MUX),
			ID(CFG_SEND_F_ERR_MUX),
			ID(CFG_SYS_ERR_RC),
			ID(CFG_AER_RC_ERR_MUX),
			ID(RADM_CPL_TIMEOUT),
			ID(RADM_TIMEOUT_CPL_TC),
			ID(RADM_TIMEOUT_CPL_TAG),
			ID(RADM_TIMEOUT_CPL_ATTR),
			ID(RADM_TIMEOUT_CPL_LEN),
			ID(CFG_MAX_RD_REQ_SIZE),
			ID(CFG_BUS_MASTER_EN),
			ID(CFG_MAX_PAYLOAD_SIZE),
			ID(CFG_RCB),
			ID(CFG_MEM_SPACE_EN),
			ID(CFG_PM_NO_SOFT_RST),
			ID(CFG_CRS_SW_VIS_EN),
			ID(CFG_NO_SNOOP_EN),
			ID(CFG_RELAX_ORDER_EN),
			ID(CFG_TPH_REQ_EN),
			ID(CFG_PF_TPH_ST_MODE),
			ID(CFG_PBUS_NUM),
			ID(CFG_PBUS_DEV_NUM),
			ID(RBAR_CTRL_UPDATE),
			ID(CFG_ATOMIC_REQ_EN),
			ID(CFG_ATOMIC_EGRESS_BLOCK),
			ID(CFG_EXT_TAG_EN),
			ID(RADM_IDLE),
			ID(RADM_Q_NOT_EMPTY),
			ID(RADM_QOVERFLOW),
			ID(CFG_LINK_AUTO_BW_MUX),
			ID(CFG_BW_MGT_MUX),
			ID(CFG_PME_MUX),
			ID(DEBUG_INFO_MUX),
			ID(CFG_IDO_REQ_EN),
			ID(CFG_IDO_CPL_EN),
			ID(XADM_PH_CDTS),
			ID(XADM_PD_CDTS),
			ID(XADM_NPH_CDTS),
			ID(XADM_NPD_CDTS),
			ID(XADM_CPLH_CDTS),
			ID(XADM_CPLD_CDTS),
			ID(MAC_PHY_POWERDOWN),
			ID(MAC_PHY_TXDATA),
			ID(MAC_PHY_TXDATAK),
			ID(MAC_PHY_TXDETECTRX_LOOPBACK),
			ID(MAC_PHY_TXELECIDLE_L),
			ID(MAC_PHY_TXELECIDLE_H),
			ID(MAC_PHY_TXCOMPLIANCE),
			ID(MAC_PHY_RXPOLARITY),
			ID(MAC_PHY_RATE),
			ID(MAC_PHY_TXDEEMPH),
			ID(MAC_PHY_TXMARGIN),
			ID(MAC_PHY_TXSWING),
			ID(CFG_HW_AUTO_SP_DIS),
			ID(P_DATAQ_ADDRA),
			ID(P_DATAQ_ADDRB),
			ID(P_DATAQ_DATAIN),
			ID(P_DATAQ_ENA),
			ID(P_DATAQ_ENB),
			ID(P_DATAQ_WEA),
			ID(XDLH_RETRYRAM_ADDR),
			ID(XDLH_RETRYRAM_DATA),
			ID(XDLH_RETRYRAM_WE),
			ID(XDLH_RETRYRAM_EN),
			ID(P_HDRQ_ADDRA),
			ID(P_HDRQ_ADDRB),
			ID(P_HDRQ_DATAIN),
			ID(P_HDRQ_ENA),
			ID(P_HDRQ_ENB),
			ID(P_HDRQ_WEA)},
		       false);
	ct->setup_type(ID(GTP_PCIEGEN3),
		       {ID(PCLK),
			ID(PCLK_DIV2),
			ID(MEM_CLK),
			ID(USER_CLK),
			ID(BUTTON_RST),
			ID(POWER_UP_RST),
			ID(PERST),
			ID(APP_INIT_RST),
			ID(DEVICE_TYPE),
			ID(RX_LANE_FLIP_EN),
			ID(TX_LANE_FLIP_EN),
			ID(APP_LTSSM_ENABLE),
			ID(APP_REQ_RETRY_EN),
			ID(APP_PF_REQ_RETRY_EN),
			ID(APP_VF_REQ_RETRY_EN),
			ID(AXIS_MASTER0_TREADY),
			ID(USER_RCVD_NP_READY),
			ID(USER_RCVD_P_READY),
			ID(AXIS_SLAVE0_TDATA),
			ID(AXIS_SLAVE0_TLAST),
			ID(AXIS_SLAVE0_TUSER),
			ID(AXIS_SLAVE0_TVALID),
			ID(AXIS_SLAVE1_TDATA),
			ID(AXIS_SLAVE1_TLAST),
			ID(AXIS_SLAVE1_TUSER),
			ID(AXIS_SLAVE1_TVALID),
			ID(AXIS_SLAVE2_TDATA),
			ID(AXIS_SLAVE2_TLAST),
			ID(AXIS_SLAVE2_TUSER),
			ID(AXIS_SLAVE2_TVALID),
			ID(APP_DBI_RO_WR_DISABLE),
			ID(DBI_ADDR),
			ID(DBI_CS),
			ID(DBI_CS2),
			ID(DBI_DIN),
			ID(DBI_FUNC_NUM),
			ID(DBI_VFUNC_ACTIVE),
			ID(DBI_VFUNC_NUM),
			ID(DBI_WR),
			ID(SEDI),
			ID(SEDI_ACK),
			ID(SYS_INT),
			ID(CFG_MSI_PENDING),
			ID(CFG_VF_MSI_PENDING),
			ID(VEN_MSI_FUNC_NUM),
			ID(VEN_MSI_REQ),
			ID(VEN_MSI_TC),
			ID(VEN_MSI_VECTOR),
			ID(VEN_MSI_VFUNC_ACTIVE),
			ID(VEN_MSI_VFUNC_NUM),
			ID(MSIX_ADDR),
			ID(MSIX_DATA),
			ID(APP_CLK_PM_EN),
			ID(APPS_PM_VF_XMT_PME),
			ID(APPS_PM_XMT_PME),
			ID(APPS_PM_XMT_TURNOFF),
			ID(APP_UNLOCK_MSG),
			ID(APP_READY_ENTR_L23),
			ID(APP_REQ_ENTR_L1),
			ID(APP_REQ_EXIT_L1),
			ID(CFG_PWR_BUDGET_DATA_REG),
			ID(CFG_PWR_BUDGET_FUNC_NUM),
			ID(CFG_PWR_BUDGET_VALID),
			ID(APP_XFER_PENDING),
			ID(APP_HDR_LOG),
			ID(APP_HDR_VALID),
			ID(APP_ERR_ADVISORY),
			ID(APP_ERR_BUS),
			ID(APP_ERR_FUNC_NUM),
			ID(APP_ERR_VFUNC_ACTIVE),
			ID(APP_ERR_VFUNC_NUM),
			ID(APP_LTR_MSG_FUNC_NUM),
			ID(APP_LTR_MSG_LATENCY),
			ID(APP_LTR_MSG_REQ),
			ID(APP_FLR_PF_DONE),
			ID(APP_FLR_VF_DONE),
			ID(APP_OBFF_CPU_ACTIVE_MSG_REQ),
			ID(APP_OBFF_IDLE_MSG_REQ),
			ID(APP_OBFF_OBFF_MSG_REQ),
			ID(APP_RAS_DES_SD_HOLD_LTSSM),
			ID(DIAG_CTRL_BUS),
			ID(DYN_DEBUG_INFO_SEL),
			ID(PHY_MAC_DIRFEEDBACK),
			ID(PHY_MAC_FOMFEEDBACK),
			ID(PHY_MAC_LOCALFS),
			ID(PHY_MAC_LOCALLF),
			ID(PHY_MAC_LOCAL_TX_COEF_VALID),
			ID(PHY_MAC_LOCAL_TX_PSET_COEF),
			ID(PHY_MAC_PHYSTATUS),
			ID(PHY_MAC_RXDATA),
			ID(PHY_MAC_RXDATAK),
			ID(PHY_MAC_RXDATAVALID),
			ID(PHY_MAC_RXELECIDLE),
			ID(PHY_MAC_RXSTARTBLOCK),
			ID(PHY_MAC_RXSTATUS),
			ID(PHY_MAC_RXSYNCHEADER),
			ID(PHY_MAC_RXVALID),
			ID(PNP_RAM_RD_DATA),
			ID(RETRYRAM_XDLH_DATA),
			ID(TPH_RD_DATA_VALID),
			ID(TPH_RAM_RD_DATA),
			ID(RAM_TEST_EN),
			ID(RAM_TEST_ADDRH),
			ID(RETRY_TEST_DATA_EN),
			ID(RAM_TEST_MODE_N)},
		       {ID(USER_RST_N),
			ID(TRAINING_RST_N),
			ID(PHY_RST_N),
			ID(SMLH_LINK_UP),
			ID(RDLH_LINK_UP),
			ID(SMLH_LTSSM_STATE),
			ID(CFG_2ND_RESET),
			ID(LINK_REQ_RST),
			ID(CFG_VF_BME),
			ID(SMLH_REQ_RST),
			ID(AXIS_MASTER0_TDATA),
			ID(AXIS_MASTER0_TKEEP),
			ID(AXIS_MASTER0_TLAST),
			ID(AXIS_MASTER0_TUSER),
			ID(AXIS_MASTER0_TVALID),
			ID(CORE_AVL_NP_CNT),
			ID(CORE_AVL_P_CNT),
			ID(AXIS_MASTER1_TDATA),
			ID(AXIS_MASTER1_TKEEP),
			ID(AXIS_MASTER1_TUSER),
			ID(AXIS_MASTER1_TVALID),
			ID(AXIS_SLAVE0_TREADY),
			ID(AXIS_SLAVE1_TREADY),
			ID(AXIS_SLAVE2_TREADY),
			ID(RADM_CPL_TIMEOUT),
			ID(RADM_TIMEOUT_CPL_ATTR),
			ID(RADM_TIMEOUT_CPL_LEN),
			ID(RADM_TIMEOUT_CPL_TAG),
			ID(RADM_TIMEOUT_CPL_TC),
			ID(RADM_TIMEOUT_FUNC_NUM),
			ID(RADM_TIMEOUT_VFUNC_ACTIVE),
			ID(RADM_TIMEOUT_VFUNC_NUM),
			ID(LBC_DBI_ACK),
			ID(LBC_DBI_DOUT),
			ID(SEDO),
			ID(SEDO_EN),
			ID(CFG_INT_DISABLE),
			ID(INT_GRT),
			ID(CFG_MSI_EN),
			ID(CFG_MSI_MASK_UPDATE),
			ID(CFG_MULTI_MSI_EN),
			ID(CFG_VF_MSI_EN),
			ID(CFG_VF_MULTI_MSI_EN),
			ID(VEN_MSI_GRANT),
			ID(CFG_MSIX_EN),
			ID(CFG_MSIX_FUNC_MASK),
			ID(CFG_VF_MSIX_EN),
			ID(CFG_VF_MSIX_FUNC_MASK),
			ID(CFG_BW_MGT_MSI),
			ID(CFG_LINK_AUTO_BW_MSI),
			ID(CFG_PME_MSI),
			ID(CFG_BW_MGT_INT),
			ID(CFG_LINK_AUTO_BW_INT),
			ID(CFG_LINK_EQ_REQ_INT),
			ID(CFG_PME_INT),
			ID(CFG_NF_ERR_RPT_EN),
			ID(CFG_NO_SNOOP_EN),
			ID(CFG_OBFF_EN),
			ID(CFG_PBUS_DEV_NUM),
			ID(CFG_PBUS_NUM),
			ID(CFG_MEM_SPACE_EN),
			ID(CFG_EXT_TAG_EN),
			ID(CFG_F_ERR_RPT_EN),
			ID(CFG_ARI_FWD_EN),
			ID(CFG_ATOMIC_EGRESS_BLOCK),
			ID(CFG_ATOMIC_REQ_EN),
			ID(CFG_BUS_MASTER_EN),
			ID(CFG_COR_ERR_RPT_EN),
			ID(CFG_CRS_SW_VIS_EN),
			ID(CFG_MAX_PAYLOAD_SIZE),
			ID(CFG_MAX_RD_REQ_SIZE),
			ID(CFG_VF_EN),
			ID(CFG_TC_ENABLE),
			ID(CFG_RCB),
			ID(CFG_REG_SERREN),
			ID(CFG_RELAX_ORDER_EN),
			ID(RBAR_CTRL_UPDATE),
			ID(AUX_PM_EN),
			ID(PM_DSTATE),
			ID(PM_MASTER_STATE),
			ID(PM_PME_EN),
			ID(PM_SLAVE_STATE),
			ID(PM_STATUS),
			ID(PM_VF_DSTATE),
			ID(PM_VF_PME_EN),
			ID(PM_VF_STATUS),
			ID(PM_XTLH_BLOCK_TLP),
			ID(DPA_SUBSTATE_UPDATE),
			ID(WAKE),
			ID(CFG_PWR_BUDGET_DATA_SEL_REG),
			ID(CFG_PWR_BUDGET_SEL),
			ID(CFG_SEND_COR_ERR),
			ID(CFG_SEND_F_ERR),
			ID(CFG_SEND_NF_ERR),
			ID(CFG_AER_RC_ERR_INT),
			ID(CFG_AER_RC_ERR_MSI),
			ID(CFG_SYS_ERR_RC),
			ID(APP_LTR_MSG_GRANT),
			ID(CFG_LTR_M_EN),
			ID(CFG_DISABLE_LTR_CLR_MSG),
			ID(CFG_FLR_PF_ACTIVE),
			ID(CFG_FLR_VF_ACTIVE),
			ID(CFG_START_VFI),
			ID(CFG_NUM_VF),
			ID(APP_OBFF_MSG_GRANT),
			ID(MSG_RCVD),
			ID(MSG_RCVD_DATA),
			ID(MSG_RCVD_TYPE),
			ID(DEBUG_INFO_MUX),
			ID(CFG_IDO_CPL_EN),
			ID(CFG_IDO_REQ_EN),
			ID(XADM_CPLD_CDTS),
			ID(XADM_CPLH_CDTS),
			ID(XADM_NPD_CDTS),
			ID(XADM_NPH_CDTS),
			ID(XADM_PD_CDTS),
			ID(XADM_PH_CDTS),
			ID(RADM_Q_NOT_EMPTY),
			ID(RADM_QOVERFLOW),
			ID(MAC_PHY_BLOCKALIGNCONTROL),
			ID(MAC_PHY_DIRCHANGE),
			ID(MAC_PHY_FS),
			ID(MAC_PHY_GETLOCAL_PSET_COEF),
			ID(MAC_PHY_INVALID_REQ),
			ID(MAC_PHY_LF),
			ID(MAC_PHY_LOCAL_PSET_INDEX),
			ID(MAC_PHY_POWERDOWN),
			ID(MAC_PHY_RATE),
			ID(MAC_PHY_RXEQEVAL),
			ID(MAC_PHY_RXEQINPROGRESS),
			ID(MAC_PHY_RXPOLARITY),
			ID(MAC_PHY_RXPRESETHINT),
			ID(MAC_PHY_TXCOMPLIANCE),
			ID(MAC_PHY_TXDATA),
			ID(MAC_PHY_TXDATAK),
			ID(MAC_PHY_TXDATAVALID),
			ID(MAC_PHY_TXDEEMPH),
			ID(MAC_PHY_TXDETECTRX_LOOPBACK),
			ID(MAC_PHY_TXELECIDLE_H),
			ID(MAC_PHY_TXELECIDLE_L),
			ID(MAC_PHY_TXMARGIN),
			ID(MAC_PHY_TXSTARTBLOCK),
			ID(MAC_PHY_TXSWING),
			ID(MAC_PHY_TXSYNCHEADER),
			ID(PNP_RAM_RD_ADDR),
			ID(PNP_RAM_RD_EN),
			ID(PNP_RAM_WR_ADDR),
			ID(PNP_RAM_WR_DATA),
			ID(PNP_RAM_WR_EN),
			ID(XDLH_RETRYRAM_ADDR),
			ID(XDLH_RETRYRAM_DATA),
			ID(XDLH_RETRYRAM_EN),
			ID(XDLH_RETRYRAM_WE),
			ID(CFG_PF_TPH_ST_MODE),
			ID(CFG_TPH_REQ_EN),
			ID(CFG_VF_TPH_REQ_EN),
			ID(CFG_VF_TPH_ST_MODE),
			ID(TPH_RAM_ADDR),
			ID(TPH_RAM_FUNC_NUM),
			ID(TPH_RAM_FUNC_ACTIVE),
			ID(TPH_RAM_WR_BYTE_EN),
			ID(TPH_RAM_WR_DATA),
			ID(TPH_RAM_WR_EN)},
		       false);
	ct->setup_type(ID(GTP_PCIEGEN3_DFT),
		       {ID(PCLK),
			ID(PCLK_DIV2),
			ID(MEM_CLK),
			ID(USER_CLK),
			ID(BUTTON_RST),
			ID(POWER_UP_RST),
			ID(PERST),
			ID(APP_INIT_RST),
			ID(DEVICE_TYPE),
			ID(RX_LANE_FLIP_EN),
			ID(TX_LANE_FLIP_EN),
			ID(APP_LTSSM_ENABLE),
			ID(APP_REQ_RETRY_EN),
			ID(APP_PF_REQ_RETRY_EN),
			ID(APP_VF_REQ_RETRY_EN),
			ID(AXIS_MASTER0_TREADY),
			ID(USER_RCVD_NP_READY),
			ID(USER_RCVD_P_READY),
			ID(AXIS_SLAVE0_TDATA),
			ID(AXIS_SLAVE0_TLAST),
			ID(AXIS_SLAVE0_TUSER),
			ID(AXIS_SLAVE0_TVALID),
			ID(AXIS_SLAVE1_TDATA),
			ID(AXIS_SLAVE1_TLAST),
			ID(AXIS_SLAVE1_TUSER),
			ID(AXIS_SLAVE1_TVALID),
			ID(AXIS_SLAVE2_TDATA),
			ID(AXIS_SLAVE2_TLAST),
			ID(AXIS_SLAVE2_TUSER),
			ID(AXIS_SLAVE2_TVALID),
			ID(APP_DBI_RO_WR_DISABLE),
			ID(DBI_ADDR),
			ID(DBI_CS),
			ID(DBI_CS2),
			ID(DBI_DIN),
			ID(DBI_FUNC_NUM),
			ID(DBI_VFUNC_ACTIVE),
			ID(DBI_VFUNC_NUM),
			ID(DBI_WR),
			ID(SEDI),
			ID(SEDI_ACK),
			ID(SYS_INT),
			ID(CFG_MSI_PENDING),
			ID(CFG_VF_MSI_PENDING),
			ID(VEN_MSI_FUNC_NUM),
			ID(VEN_MSI_REQ),
			ID(VEN_MSI_TC),
			ID(VEN_MSI_VECTOR),
			ID(VEN_MSI_VFUNC_ACTIVE),
			ID(VEN_MSI_VFUNC_NUM),
			ID(MSIX_ADDR),
			ID(MSIX_DATA),
			ID(APP_CLK_PM_EN),
			ID(APPS_PM_VF_XMT_PME),
			ID(APPS_PM_XMT_PME),
			ID(APPS_PM_XMT_TURNOFF),
			ID(APP_UNLOCK_MSG),
			ID(APP_READY_ENTR_L23),
			ID(APP_REQ_ENTR_L1),
			ID(APP_REQ_EXIT_L1),
			ID(CFG_PWR_BUDGET_DATA_REG),
			ID(CFG_PWR_BUDGET_FUNC_NUM),
			ID(CFG_PWR_BUDGET_VALID),
			ID(APP_XFER_PENDING),
			ID(APP_HDR_LOG),
			ID(APP_HDR_VALID),
			ID(APP_ERR_ADVISORY),
			ID(APP_ERR_BUS),
			ID(APP_ERR_FUNC_NUM),
			ID(APP_ERR_VFUNC_ACTIVE),
			ID(APP_ERR_VFUNC_NUM),
			ID(APP_LTR_MSG_FUNC_NUM),
			ID(APP_LTR_MSG_LATENCY),
			ID(APP_LTR_MSG_REQ),
			ID(APP_FLR_PF_DONE),
			ID(APP_FLR_VF_DONE),
			ID(APP_OBFF_CPU_ACTIVE_MSG_REQ),
			ID(APP_OBFF_IDLE_MSG_REQ),
			ID(APP_OBFF_OBFF_MSG_REQ),
			ID(APP_RAS_DES_SD_HOLD_LTSSM),
			ID(DIAG_CTRL_BUS),
			ID(DYN_DEBUG_INFO_SEL),
			ID(PHY_MAC_DIRFEEDBACK),
			ID(PHY_MAC_FOMFEEDBACK),
			ID(PHY_MAC_LOCALFS),
			ID(PHY_MAC_LOCALLF),
			ID(PHY_MAC_LOCAL_TX_COEF_VALID),
			ID(PHY_MAC_LOCAL_TX_PSET_COEF),
			ID(PHY_MAC_PHYSTATUS),
			ID(PHY_MAC_RXDATA),
			ID(PHY_MAC_RXDATAK),
			ID(PHY_MAC_RXDATAVALID),
			ID(PHY_MAC_RXELECIDLE),
			ID(PHY_MAC_RXSTARTBLOCK),
			ID(PHY_MAC_RXSTATUS),
			ID(PHY_MAC_RXSYNCHEADER),
			ID(PHY_MAC_RXVALID),
			ID(PNP_RAM_RD_DATA),
			ID(RETRYRAM_XDLH_DATA),
			ID(TPH_RD_DATA_VALID),
			ID(TPH_RAM_RD_DATA),
			ID(RAM_TEST_EN),
			ID(RAM_TEST_ADDRH),
			ID(RETRY_TEST_DATA_EN),
			ID(RAM_TEST_MODE_N),
			ID(TEST_MODE_N),
			ID(TEST_RST_N),
			ID(TEST_SE_N)},
		       {ID(USER_RST_N),
			ID(TRAINING_RST_N),
			ID(PHY_RST_N),
			ID(SMLH_LINK_UP),
			ID(RDLH_LINK_UP),
			ID(SMLH_LTSSM_STATE),
			ID(CFG_2ND_RESET),
			ID(LINK_REQ_RST),
			ID(CFG_VF_BME),
			ID(SMLH_REQ_RST),
			ID(AXIS_MASTER0_TDATA),
			ID(AXIS_MASTER0_TKEEP),
			ID(AXIS_MASTER0_TLAST),
			ID(AXIS_MASTER0_TUSER),
			ID(AXIS_MASTER0_TVALID),
			ID(CORE_AVL_NP_CNT),
			ID(CORE_AVL_P_CNT),
			ID(AXIS_MASTER1_TDATA),
			ID(AXIS_MASTER1_TKEEP),
			ID(AXIS_MASTER1_TUSER),
			ID(AXIS_MASTER1_TVALID),
			ID(AXIS_SLAVE0_TREADY),
			ID(AXIS_SLAVE1_TREADY),
			ID(AXIS_SLAVE2_TREADY),
			ID(RADM_CPL_TIMEOUT),
			ID(RADM_TIMEOUT_CPL_ATTR),
			ID(RADM_TIMEOUT_CPL_LEN),
			ID(RADM_TIMEOUT_CPL_TAG),
			ID(RADM_TIMEOUT_CPL_TC),
			ID(RADM_TIMEOUT_FUNC_NUM),
			ID(RADM_TIMEOUT_VFUNC_ACTIVE),
			ID(RADM_TIMEOUT_VFUNC_NUM),
			ID(LBC_DBI_ACK),
			ID(LBC_DBI_DOUT),
			ID(SEDO),
			ID(SEDO_EN),
			ID(CFG_INT_DISABLE),
			ID(INT_GRT),
			ID(CFG_MSI_EN),
			ID(CFG_MSI_MASK_UPDATE),
			ID(CFG_MULTI_MSI_EN),
			ID(CFG_VF_MSI_EN),
			ID(CFG_VF_MULTI_MSI_EN),
			ID(VEN_MSI_GRANT),
			ID(CFG_MSIX_EN),
			ID(CFG_MSIX_FUNC_MASK),
			ID(CFG_VF_MSIX_EN),
			ID(CFG_VF_MSIX_FUNC_MASK),
			ID(CFG_BW_MGT_MSI),
			ID(CFG_LINK_AUTO_BW_MSI),
			ID(CFG_PME_MSI),
			ID(CFG_BW_MGT_INT),
			ID(CFG_LINK_AUTO_BW_INT),
			ID(CFG_LINK_EQ_REQ_INT),
			ID(CFG_PME_INT),
			ID(CFG_NF_ERR_RPT_EN),
			ID(CFG_NO_SNOOP_EN),
			ID(CFG_OBFF_EN),
			ID(CFG_PBUS_DEV_NUM),
			ID(CFG_PBUS_NUM),
			ID(CFG_MEM_SPACE_EN),
			ID(CFG_EXT_TAG_EN),
			ID(CFG_F_ERR_RPT_EN),
			ID(CFG_ARI_FWD_EN),
			ID(CFG_ATOMIC_EGRESS_BLOCK),
			ID(CFG_ATOMIC_REQ_EN),
			ID(CFG_BUS_MASTER_EN),
			ID(CFG_COR_ERR_RPT_EN),
			ID(CFG_CRS_SW_VIS_EN),
			ID(CFG_MAX_PAYLOAD_SIZE),
			ID(CFG_MAX_RD_REQ_SIZE),
			ID(CFG_VF_EN),
			ID(CFG_TC_ENABLE),
			ID(CFG_RCB),
			ID(CFG_REG_SERREN),
			ID(CFG_RELAX_ORDER_EN),
			ID(RBAR_CTRL_UPDATE),
			ID(AUX_PM_EN),
			ID(PM_DSTATE),
			ID(PM_MASTER_STATE),
			ID(PM_PME_EN),
			ID(PM_SLAVE_STATE),
			ID(PM_STATUS),
			ID(PM_VF_DSTATE),
			ID(PM_VF_PME_EN),
			ID(PM_VF_STATUS),
			ID(PM_XTLH_BLOCK_TLP),
			ID(DPA_SUBSTATE_UPDATE),
			ID(WAKE),
			ID(CFG_PWR_BUDGET_DATA_SEL_REG),
			ID(CFG_PWR_BUDGET_SEL),
			ID(CFG_SEND_COR_ERR),
			ID(CFG_SEND_F_ERR),
			ID(CFG_SEND_NF_ERR),
			ID(CFG_AER_RC_ERR_INT),
			ID(CFG_AER_RC_ERR_MSI),
			ID(CFG_SYS_ERR_RC),
			ID(APP_LTR_MSG_GRANT),
			ID(CFG_LTR_M_EN),
			ID(CFG_DISABLE_LTR_CLR_MSG),
			ID(CFG_FLR_PF_ACTIVE),
			ID(CFG_FLR_VF_ACTIVE),
			ID(CFG_START_VFI),
			ID(CFG_NUM_VF),
			ID(APP_OBFF_MSG_GRANT),
			ID(MSG_RCVD),
			ID(MSG_RCVD_DATA),
			ID(MSG_RCVD_TYPE),
			ID(DEBUG_INFO_MUX),
			ID(CFG_IDO_CPL_EN),
			ID(CFG_IDO_REQ_EN),
			ID(XADM_CPLD_CDTS),
			ID(XADM_CPLH_CDTS),
			ID(XADM_NPD_CDTS),
			ID(XADM_NPH_CDTS),
			ID(XADM_PD_CDTS),
			ID(XADM_PH_CDTS),
			ID(RADM_Q_NOT_EMPTY),
			ID(RADM_QOVERFLOW),
			ID(MAC_PHY_BLOCKALIGNCONTROL),
			ID(MAC_PHY_DIRCHANGE),
			ID(MAC_PHY_FS),
			ID(MAC_PHY_GETLOCAL_PSET_COEF),
			ID(MAC_PHY_INVALID_REQ),
			ID(MAC_PHY_LF),
			ID(MAC_PHY_LOCAL_PSET_INDEX),
			ID(MAC_PHY_POWERDOWN),
			ID(MAC_PHY_RATE),
			ID(MAC_PHY_RXEQEVAL),
			ID(MAC_PHY_RXEQINPROGRESS),
			ID(MAC_PHY_RXPOLARITY),
			ID(MAC_PHY_RXPRESETHINT),
			ID(MAC_PHY_TXCOMPLIANCE),
			ID(MAC_PHY_TXDATA),
			ID(MAC_PHY_TXDATAK),
			ID(MAC_PHY_TXDATAVALID),
			ID(MAC_PHY_TXDEEMPH),
			ID(MAC_PHY_TXDETECTRX_LOOPBACK),
			ID(MAC_PHY_TXELECIDLE_H),
			ID(MAC_PHY_TXELECIDLE_L),
			ID(MAC_PHY_TXMARGIN),
			ID(MAC_PHY_TXSTARTBLOCK),
			ID(MAC_PHY_TXSWING),
			ID(MAC_PHY_TXSYNCHEADER),
			ID(PNP_RAM_RD_ADDR),
			ID(PNP_RAM_RD_EN),
			ID(PNP_RAM_WR_ADDR),
			ID(PNP_RAM_WR_DATA),
			ID(PNP_RAM_WR_EN),
			ID(XDLH_RETRYRAM_ADDR),
			ID(XDLH_RETRYRAM_DATA),
			ID(XDLH_RETRYRAM_EN),
			ID(XDLH_RETRYRAM_WE),
			ID(CFG_PF_TPH_ST_MODE),
			ID(CFG_TPH_REQ_EN),
			ID(CFG_VF_TPH_REQ_EN),
			ID(CFG_VF_TPH_ST_MODE),
			ID(TPH_RAM_ADDR),
			ID(TPH_RAM_FUNC_NUM),
			ID(TPH_RAM_FUNC_ACTIVE),
			ID(TPH_RAM_WR_BYTE_EN),
			ID(TPH_RAM_WR_DATA),
			ID(TPH_RAM_WR_EN)},
		       false);
	ct->setup_type(ID(GTP_PPLL),
		       {ID(CLKIN1),	 ID(CLKIN2),	  ID(CLKFB),	     ID(CLKIN_SEL),   ID(CLKOUT0_SYN), ID(CLKOUT1_SYN), ID(CLKOUT2_SYN),
			ID(CLKOUT3_SYN), ID(CLKOUT4_SYN), ID(CLKOUTPHY_SYN), ID(CLKOUTF_SYN), ID(PLL_PWD),     ID(RST),		ID(APB_CLK),
			ID(APB_RST_N),	 ID(APB_ADDR),	  ID(APB_SEL),	     ID(APB_EN),      ID(APB_WRITE),   ID(APB_WDATA)},
		       {ID(CLKOUT0), ID(CLKOUT0N), ID(CLKOUT1), ID(CLKOUT1N), ID(CLKOUT2), ID(CLKOUT2N), ID(CLKOUT3), ID(CLKOUT3N), ID(CLKOUT4),
			ID(CLKOUTPHY), ID(CLKOUTPHYN), ID(CLKOUTF), ID(CLKOUTFN), ID(LOCK), ID(APB_RDATA), ID(APB_READY)},
		       false);
	ct->setup_type(ID(GTP_PPLL_DFT),
		       {ID(CLKIN1),	 ID(CLKIN2),	  ID(CLKFB),	     ID(CLKIN_SEL),   ID(CLKOUT0_SYN), ID(CLKOUT1_SYN), ID(CLKOUT2_SYN),
			ID(CLKOUT3_SYN), ID(CLKOUT4_SYN), ID(CLKOUTPHY_SYN), ID(CLKOUTF_SYN), ID(PLL_PWD),     ID(RST),		ID(APB_CLK),
			ID(APB_RST_N),	 ID(APB_ADDR),	  ID(APB_SEL),	     ID(APB_EN),      ID(APB_WRITE),   ID(APB_WDATA)},
		       {ID(PFDTOP_CLK_TEST), ID(CLKOUT0), ID(CLKOUT0N), ID(CLKOUT1), ID(CLKOUT1N), ID(CLKOUT2), ID(CLKOUT2N), ID(CLKOUT3),
			ID(CLKOUT3N), ID(CLKOUT4), ID(CLKOUTPHY), ID(CLKOUTPHYN), ID(CLKOUTF), ID(CLKOUTFN), ID(LOCK), ID(APB_RDATA), ID(APB_READY)},
		       false);
	ct->setup_type(ID(GTP_RAM32X1DP), {ID(DI), ID(RADDR), ID(WADDR), ID(WCLK), ID(WE)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_RAM32X1SP), {ID(DI), ID(ADDR), ID(WCLK), ID(WE)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_RAM32X2DP), {ID(DI), ID(RADDR), ID(WADDR), ID(WCLK), ID(WE)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_RAM32X2SP), {ID(DI), ID(ADDR), ID(WCLK), ID(WE)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_RAM32X2X4), {ID(DI0), ID(DI1), ID(DI2), ID(DI3), ID(ADDR0), ID(ADDR1), ID(ADDR2), ID(ADDR3), ID(WCLK), ID(WE)},
		       {ID(DO0), ID(DO1), ID(DO2), ID(DO3)}, false);
	ct->setup_type(ID(GTP_RAM64X1DP), {ID(DI), ID(RADDR), ID(WADDR), ID(WCLK), ID(WE)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_RAM64X1SP), {ID(DI), ID(ADDR), ID(WCLK), ID(WE)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_RAM64X1X4), {ID(DI0), ID(DI1), ID(DI2), ID(DI3), ID(ADDR0), ID(ADDR1), ID(ADDR2), ID(ADDR3), ID(WCLK), ID(WE)},
		       {ID(DO0), ID(DO1), ID(DO2), ID(DO3)}, false);
	ct->setup_type(ID(GTP_RAM128X1DP), {ID(DI), ID(RADDR), ID(WADDR), ID(WCLK), ID(WE)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_RAM128X1SP), {ID(DI), ID(ADDR), ID(WCLK), ID(WE)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_RAM256X1SP), {ID(DI), ID(ADDR), ID(WCLK), ID(WE)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_RES_CAL_E1), {ID(PCODE_IN), ID(NCODE_IN), ID(SAMPLE_IN), ID(EN), ID(CODE_SEL), ID(RST_N)},
		       {ID(PCODE_OUT), ID(NCODE_OUT), ID(CAL_DONE)}, false);
	ct->setup_type(ID(GTP_RES_CAL_E1_DFT), {ID(PCODE_IN), ID(NCODE_IN), ID(SAMPLE_IN), ID(EN), ID(CODE_SEL), ID(RST_N)},
		       {ID(PCODE_OUT), ID(NCODE_OUT), ID(CAL_DONE)}, false);
	ct->setup_type(ID(GTP_RES_CAL_E2), {ID(PCODE_IN), ID(NCODE_IN), ID(SAMPLE_IN), ID(EN), ID(CODE_SEL), ID(RST_N)},
		       {ID(PCODE_OUT), ID(NCODE_OUT), ID(CAL_DONE)}, false);
	ct->setup_type(ID(GTP_ROM32X1), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_ROM32X2), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_ROM64X1), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_ROM128X1), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5), ID(I6)}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_ROM256X1), {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5), ID(I6), ID(I7)}, {ID(Z)}, false);
	ct->setup_type(
	  ID(GTP_SCANCHAIN_E1), {ID(TCK), ID(TDI), ID(TMS), ID(TDO_USER)},
	  {ID(TDO), ID(CAPDR), ID(SHFTDR), ID(UPDR), ID(JCLK), ID(RST), ID(JRTI), ID(FLG_USER), ID(TCK_USER), ID(TDI_USER), ID(TMS_USER)}, false);
	ct->setup_type(ID(GTP_START_E1), {ID(CLK), ID(GOE), ID(GRS_N), ID(GWE)}, {ID(WAKEUP_OVER)}, false);
	ct->setup_type(ID(GTP_UDID), {ID(DI), ID(LOAD), ID(SE), ID(CLK)}, {ID(DO)}, false);
	ct->setup_type(ID(GTP_ZERO), {}, {ID(Z)}, false);
	ct->setup_type(ID(GTP_ZEROHOLDDELAY), {ID(DI)}, {ID(DO)}, false);
}



RTLIL::Cell *addLut(Module *module, const pool<SigBit> &cut, const RTLIL::SigBit &sig_z)
{
	log_assert(cut.size() <= LUT_SIZE && cut.size() >= 1);
	vector<SigBit> vcut;
	for (auto bit : cut) {
		vcut.push_back(bit);
	}
	vector<bool> cut_init_bools = GetCutInit(vcut, sig_z);
	Cell *drv = bit2driver[sig_z];
	log_assert(drv);
	IdString name = drv->name;
	string new_name = string(name.c_str()) + "_lut";
	Cell *cell = nullptr;
	if (using_internel_lut_type) {
		// instantiate $lut, need call techmap pass
		cell = module->addCell(IdString(new_name), ID($lut));
		cell->parameters[ID::WIDTH] = RTLIL::Const(vcut.size());
		cell->parameters[ID::LUT] = RTLIL::Const(cut_init_bools);

		cell->setPort(ID(A), vcut);
		cell->setPort(ID(Y), sig_z);
		cell->set_src_attribute(drv->get_src_attribute());
	} else {
		IdString types[] = {ID(GTP_LUT1), ID(GTP_LUT2), ID(GTP_LUT3), ID(GTP_LUT4), ID(GTP_LUT5), ID(GTP_LUT6)};
		cell = module->addCell(RTLIL::IdString(new_name), types[vcut.size() - 1]);
		cell->parameters[ID::INIT] = RTLIL::Const(cut_init_bools);
		for (size_t i = 0; i < vcut.size(); ++i) {
			string pin_name = "\\I" + to_string(i);
			cell->setPort(RTLIL::IdString(pin_name), vcut[i]);
		}
		cell->setPort(ID(Z), sig_z);
		cell->set_src_attribute(drv->get_src_attribute());
	}
	return cell;
}

RTLIL::Cell* addDualOutputLut(Module *module,
                              const pool<SigBit> &merged_cut,
                              const std::vector<SigBit> &outputs,
                              dict<SigBit, pool<SigBit>> &bit2cut)
{
    if (outputs.size() != 2 || merged_cut.empty() || merged_cut.size() > 6)
        return nullptr;

    SigBit out_a = outputs[0];
    SigBit out_b = outputs[1];
    pool<SigBit> cut_a = bit2cut[out_a];
    pool<SigBit> cut_b = bit2cut[out_b];

    auto sanitize_cut = [&](const pool<SigBit> &cut) {
        pool<SigBit> clean;
        for (auto bit : cut) {
            if (bit == out_a || bit == out_b)
                continue;
            clean.insert(bit);
        }
        return clean;
    };

    cut_a = sanitize_cut(cut_a);
    cut_b = sanitize_cut(cut_b);

    bool a_is6 = cut_a.size() == 6;
    bool b_is6 = cut_b.size() == 6;

    SigBit out_z = out_b;
    SigBit out_z5 = out_a;
    pool<SigBit> cut_z = cut_b;
    pool<SigBit> cut_z5 = cut_a;

    if (a_is6 && !b_is6) {
        out_z = out_a;
        out_z5 = out_b;
        cut_z = sanitize_cut(cut_a);
        cut_z5 = sanitize_cut(cut_b);
    } else if (b_is6 && !a_is6) {
        out_z = out_b;
        out_z5 = out_a;
         cut_z = sanitize_cut(cut_b);
        cut_z5 = sanitize_cut(cut_a);
    } else if (a_is6 && b_is6) {
        // 当前不支持 6PI + 6PI 融合
        return nullptr;
    }

    dict<SigBit, int> pin_index;
    std::vector<SigBit> pins;

    bool share_inputs = false;
    for (auto bit : cut_z5)
        if (cut_z.count(bit)) {
            share_inputs = true;
            break;
        }

    if (!share_inputs)
        return nullptr;

     auto add_cut_bits_unique = [&](const pool<SigBit> &cut) {
        for (auto bit : cut) {
            if (!pin_index.count(bit)) {
                pin_index[bit] = pins.size();
                pins.push_back(bit);
            }
        }
    };

    add_cut_bits_unique(cut_z5);
    add_cut_bits_unique(cut_z);

    if (pins.size() > 6)
        return nullptr;

    size_t input_count = pins.size();
    bool z_uses_6_inputs = cut_z.size() == 6;

    std::vector<bool> init64(64, false);
    size_t comb_limit = 1ull << input_count;
    for (size_t mask = 0; mask < comb_limit; ++mask) {
        dict<SigBit, State> assignment;
        for (size_t i = 0; i < input_count; ++i)
            assignment[pins[i]] = ((mask >> i) & 1) ? State::S1 : State::S0;

        State val_z = StateEval(assignment, out_z);
        State val_z5 = StateEval(assignment, out_z5);
        if (val_z == State::Sx || val_z5 == State::Sx)
            log_error("Cannot evaluate dual-output LUT for %s/%s.\n",
                      log_signal(out_z), log_signal(out_z5));

        int idx_base = 0;
        int o5_inputs = std::min<size_t>(5, input_count);
        for (int i = 0; i < o5_inputs; ++i)
            if ((mask >> i) & 1)
                idx_base |= (1 << i);

        if (!z_uses_6_inputs) {
            // 两个输出均为 ≤5 输入：I5 接 1，使主输出使用上半区
            init64[idx_base | (1 << 5)] = (val_z == State::S1);
        } else {
            // 6 输入 + 5 输入：第 6 个输入位作为 I5
            bool i5_val = (mask >> 5) & 1;
            init64[idx_base | (i5_val << 5)] = (val_z == State::S1);
        }

        // 仅在 I5=0 平面写入 5 输入函数
        if (!z_uses_6_inputs || ((mask >> 5) & 1) == 0)
            init64[idx_base] = (val_z5 == State::S1);
    }

    auto eval_from_init = [&](bool is_z, size_t mask) {
        size_t idx_base = 0;
        size_t o5_inputs = std::min<size_t>(5, input_count);
        for (size_t i = 0; i < o5_inputs; ++i)
            if ((mask >> i) & 1)
                idx_base |= (1 << i);

        if (!is_z)
            return init64[idx_base];

        if (!z_uses_6_inputs)
            return init64[idx_base | (1 << 5)];

        bool i5_val = (mask >> 5) & 1;
        return init64[idx_base | (i5_val << 5)];
    };

    auto depends_on_pin = [&](bool is_z, size_t pin_idx) {
        if (!is_z && pin_idx >= 5)
            return false;

        for (size_t mask = 0; mask < comb_limit; ++mask) {
            if ((mask >> pin_idx) & 1)
                continue;

            size_t mask1 = mask | (1ull << pin_idx);
            bool v0 = eval_from_init(is_z, mask);
            bool v1 = eval_from_init(is_z, mask1);
            if (v0 != v1)
                return true;
        }
        return false;
    };

    size_t shared_dependencies = 0;
    for (size_t idx = 0; idx < input_count; ++idx) {
        bool z_dep = depends_on_pin(true, idx);
        bool z5_dep = depends_on_pin(false, idx);
        if (z_dep && z5_dep)
            ++shared_dependencies;
    }

    if (shared_dependencies == 0)
        return nullptr;

    Cell *drv = bit2driver[out_z];
    log_assert(drv);
    IdString drv_name = drv->name;
    std::string new_name = std::string(drv_name.c_str()) + "_lut6d";
    Cell *cell = module->addCell(RTLIL::IdString(new_name), ID(GTP_LUT6D));

    for (size_t i = 0; i < input_count && i < 6; ++i) {
        std::string pin_name = "\\I" + std::to_string(i);
        cell->setPort(RTLIL::IdString(pin_name), pins[i]);
    }

    for (size_t i = input_count; i < 5; ++i) {
        std::string pin_name = "\\I" + std::to_string(i);
        cell->setPort(RTLIL::IdString(pin_name), RTLIL::SigSpec(RTLIL::Const(0, 1)));
    }

    if (z_uses_6_inputs) {
        log_assert(input_count == 6);
        cell->setPort(RTLIL::IdString("\\I5"), pins[5]);
    } else {
        cell->setPort(RTLIL::IdString("\\I5"), RTLIL::SigSpec(RTLIL::Const(1, 1)));
    }

    cell->setPort(ID(Z), out_z);
    cell->setPort(ID(Z5), out_z5);
    cell->parameters[ID::INIT] = RTLIL::Const(init64);
    cell->set_src_attribute(drv->get_src_attribute());

    return cell;
}

bool Cone2ToLUTs(Module *module,
                dict<SigBit, pool<SigBit>> &bit2cut,
                dict<std::pair<SigBit,SigBit>, pool<SigBit>> &twoOutputCuts)
{
    static size_t vf_count = 0;

 pool<SigBit> fused_outputs;

    // 1️⃣ 先处理 dual-output cut
    for (auto &p : twoOutputCuts) {
        vf_count++;
        SigBit out1 = p.first.first;
        SigBit out2 = p.first.second;
        pool<SigBit> &inputs = p.second;

        std::vector<SigBit> outputs = {out1, out2};
        Cell *lut = addDualOutputLut(module, inputs, outputs,bit2cut); // 需要实现 dual-output addLut
        if (!lut) {
           // log_warning("Bit2oCut: dual-output LUT build failed for %s/%s, falling back to single-output mapping.\n",
            //            log_signal(out1), log_signal(out2));
            continue;
        }

        fused_outputs.insert(out1);
        fused_outputs.insert(out2);
    }

    // 2️⃣ 再处理剩余的单输出 cut
    for (auto &p : bit2cut) {
        SigBit out = p.first;

        // 已经在 dual-output 中处理过的输出就跳过
        if (fused_outputs.count(out))
            continue;

        pool<SigBit> cut = p.second;
        Cell *lut = addLut(module, cut, out); // 单输出 LUT
        log_assert(lut);
    }

    // 3️⃣ 删除原始组合逻辑门
    for (auto c : module->cells_) {
        if (IsCombinationalGate(c.second)) {
            module->remove(c.second);
        }
    }

    return true;
}




struct testpass:public ScriptPass{
	testpass():ScriptPass("test_pango","just test"){}
	bool Score =false;
	string top_module_name;
	void execute(std::vector<std::string> args, RTLIL::Design *design) override{
		for (auto &arg : args) {
			 if (arg == "-h") {
                log("\nUsage: test_pango [options]\n");
                log("  -t     Enable TEST mode\n");
                log("  -h     Show this help message\n\n");
                return; // 打印帮助后退出，不执行脚本
            }
			else if (arg == "-c") {
				Score =true;
            }
        }
		log_header(design, "Start test_pango\n");
		log_push();
		run_script(design);
		log_pop();
	}
	void script() override{
	Module *module = active_design->top_module();
	if (!module) log_cmd_error("No top module found\n");
	top_module_name = module->name.str();
	std::string input_file ="/home/yujingning/ziguang/test/"+ top_module_name.substr(1)+ ".v";
			MapperMain(module);
			if(Score){
				run("check -mapped");
			run(stringf("write_verilog -noexpr -noattr %s_syn.v", top_module_name.c_str() + 1));
			run(stringf("score -before %s -after %s_syn.v", input_file.c_str(), top_module_name.c_str() + 1));
			}
			
	}
}testpass;



PRIVATE_NAMESPACE_END
