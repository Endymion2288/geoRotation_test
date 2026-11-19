# 专项审计报告：2025 年 `p0013` Preshower `GFIT_FAILED` Run 的根因定位、`ps_dig1_vetonu` 影响与全量 `FASER_DQ_Analyis` 检查

- **日期**：2026-11-19
- **主机与环境**：`lxplus901.cern.ch`（28 cores），LCG_107（ROOT `RDataFrame`）+ Athena `24.0.41`（`xAOD` 波形读取）
- **分析与绘图脚本**：[dq_audit_2025_preshower_gfit.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/dq_audit_2025_preshower_gfit.py)、[plot_dq_audit_2025_p0013.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/plot_dq_audit_2025_p0013.py)
- **审计汇总数据**：[dq_audit_2025_summary.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/dq_audit_2025_p0013/dq_audit_2025_summary.json)、[dq_audit_2025_full.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/dq_audit_2025_p0013/dq_audit_2025_full.json)
- **完整诊断图集（3 页，PNG + 单图 PDF + 汇总 PDF）**：[dq_audit_2025_p0013_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/dq_audit_2025_p0013/figures/dq_audit_2025_p0013_report.pdf)（索引：[figure_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/dq_audit_2025_p0013/figures/figure_index.json)）

---

## 1. 阶段目标

回答关于 2025 年 Preshower `GFIT_FAILED` 异常 run 的四个具体问题：
1. 这些 run 在 `t_Preshower(Dig1) - t_νveto`（`fig13_04`）中有无体现，定量占比与形态如何？
2. 这批 run 到底出了什么问题（物理硬件还是生产重建配置/数据库缺陷），影响范围是仅限抽样到的 6 个 run 还是全部 `p0013`？
3. 这批 run 的事例数、完整性、亮度归一化产额以及其他所有探测器的物理量读数是否正常？
4. 按照 [FASER_DQ_Analyis](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/README.md) 的完整 DQ 规范（`FASER_DQ_RDF.py`、7 套 `histograms/*.yaml`、`internal_note/analysis.py`、`beam_scan.py`、`timing_scan.py`、`yield_checks/thresholds.json`）逐项核查，排查除 Preshower 外的全部隐藏差异。

---

## 2. 工作内容

1. **全量 2025 普查（Pass 1，20 workers，耗时 `80.8 s`）**：扫描 `FASER_DQ_Analyis/faser_filelists/2025_filelist.txt` 中全部 **199 个** 2025 GRL run（每 run 取首个有效文件），统计每个 run 的生产标签（`p0013` vs `p0014`）、`GoodTracks`、`Preshower0/1/2` `GFIT_FAILED` 占比、由 `localtime - triggertime - 200` 反推的 COOL `/WAVE/DAQ/TimingOffset`、以及 `VetoNu0 / Veto20 / Timing0 / CaloLo0` 正常率与 `distanceToCollidingBCID`。
2. **全 run 完整 `FASER_DQ_Analyis` 审计（Pass 2，8 workers × 3 ROOT threads，耗时 `633.0 s`）**：使用全部有效 PHYS 文件（共处理约 `1.07 × 10⁸` 事例），对以下 **16 组 `(run, tag)`** 执行完整 `FASER_DQ_RDF.py` cutflow、产额 `m ± 5s` 阈值检验、全部 7 套 YAML 分布、16 个闪烁体/量能器通道完整 `status` 位分解、`Track` 运动学与外推位置、以及 `crossingAngle`：
   - 6 个焦点 `p0013` run：`19415, 19524, 19866, 19923, 19970, 20081`
   - 3 个 `p0013 / p0014` 同数据 A/B 共有 run：`19336, 19878, 20012`（逐事件对标 `20012` 35,021 个共有事例）
   - 4 个 `p0014` 基准 run：`20170, 20279, 20481, 20591`
3. **生产日志、COOL 数据库与真实 xAOD 波形溯源**：对比 `/eos/experiment/faser/rec/2025/{p0013_log,p0014_log}/020012/`、查询 CVMFS `DBRelease/1.6.0`、`1.6.1`、`current` 的 `ALLP200.db` `/WAVE/DAQ/TimingOffset`，并核对 [WaveRecAlgsConfig.py](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Waveform/WaveRecAlgs/python/WaveRecAlgsConfig.py) 与 [WaveformReconstructionTool.cxx](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Waveform/WaveRecTools/src/WaveformReconstructionTool.cxx) 源码及 `20012 p0014` 真实 xAOD `ch 12` 波形。

---

## 3. 详细数据

### 3.1 `ps_dig1_vetonu`（`fig13_04`）中的定量体现

在 `timing_histograms.json` 的 2025 `ps_dig1_vetonu` 分布（总计 `6,837,434` 事例）中：
- 主峰 `[5, 22]` ns：`6,741,650` 事例（**`98.60%`**，来自 `p0014` 正常脉冲，峰值 `13.27 ns`）。
- **次峰 `[22, 42]` ns：`88,895` 事例（`1.300%`），中心在 `30–34 ns`**，峰值概率密度约 `9.6 × 10⁻⁴ ns⁻¹`（仅为主峰高度 `~0.45 ns⁻¹` 的 `1/470`，因此在 `fig13_04` 线性纵轴下几乎贴近零线不可见，但在对数纵轴图 `figDQ1` 中与 6 个 `p0013` run 的加权贡献逐 bin 严格重合）。
- 6 个 `p0013` run 的 `3.13M` 个 `GoodTracks` 中仅有 `~3%`（`95,161` 事例，全文件为 `699,594` 事例）侥幸通过 `(status & ~512) == 0 && raw_charge > 2.5 pC`，但它们的 `localtime` 均被卡在 `p0013` 拟合窗口左边界 `178 ns` 右侧（中位数 `180.1 ns`），给出假的 `t_PS(Dig1) - t_VN ≈ 30–34 ns`（以及假的 `t_PS(Dig2) - t_PS(Dig1) ≈ -14.5 ns`）。

![Fig DQ1: 2025 ps_dig1_vetonu linear vs log showing p0013 bump](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/dq_audit_2025_p0013/figures/figDQ1_ps_dig1_vetonu_p0013_contribution.png)

### 3.2 根因确认：COOL `/WAVE/DAQ/TimingOffset = +18 ns` 导致的波形窗口截断

- **`p0013` 生产配置**：`Calypso-22.0.49`（`reco/p0013`，2025-05-04 构建）、`DBRelease/1.6.0`、全局标签 `OFLCOND-FASER-04`（`TRACKER-ALIGN-04`）。其 COOL 数据库 `/WAVE/DAQ/TimingOffset` 对 `ch 12, 13`（`Preshower0/1`）仍沿用 2022–2024 旧值 **`+18.0 ns`**（对新增 `ch 24, 25` `Preshower2/3` 已设为 **`-35.0 ns`**）。
- **`p0014` 生产配置**：`Calypso-22.0.49`（`reco/p0014`，2025-08-27 构建）、`DBRelease/1.6.1` + 在线 `COOLOFL_DIGITIZER/CONDBR3`、全局标签 `OFLCOND-FASER-05`（`TRACKER-ALIGN-05`）。`/WAVE/DAQ/TimingOffset` 更新为 `ch 12, 13 = -40.0 ns`、`ch 24, 25 = -42.0 ns`。
- **精确窗口算术**（[WaveRecAlgsConfig.py:47](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Waveform/WaveRecAlgs/python/WaveRecAlgsConfig.py#L47) 与 [WaveformReconstructionTool.cxx:98-99](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Waveform/WaveRecTools/src/WaveformReconstructionTool.cxx#L98-L99)：`NominalTriggerTime = 200 ns`，`FitWindowStart = -20` samples，`FitWindowWidth = 100` samples，`BaselineSampleHi = 60` 即 `0–120 ns`）：
  - `p0013` `Preshower0/1` 拟合窗：`lo_edge = (200+18)/2 - 20 = 89` -> **`[178.0, 378.0]` ns**。
  - `p0014` `Preshower0/1` 拟合窗：`lo_edge = (200-40)/2 - 20 = 60` -> **`[120.0, 320.0]` ns**（已在 `Faser-Physics-020012-00000-p0014-xAOD.root` 中实测 `len(time_vector) = 101`, `[120.0 .. 320.0]` ns 验证）。
  - 2025 升级后的 `Preshower0/1` 真实信号到达时刻 `localtime` 位于 `158–169 ns`（峰值在 `< 178 ns`），导致 `p0013` 窗口 **`[178, 378]` ns 完全切掉上升沿和主峰**，仅剩下降拖尾：
    - `GFIT_FAILED (16)` 飙升至 **`84.4%–92.0%`**（`p0014` 为 `0.5%–0.7%`）；
    - `Preshower0_charge` 中位数由 **`5.56 pC`** 跌至 **`1.36–2.26 pC`**，`raw_charge` 由 `6.16 pC` 跌至 `2.56–3.46 pC`（直接解释了 `FASER_DQ_Analyis` 既有图 [Fig5_j_preshower0.png](file:///eos/home-x/xcheng/FASER/DQ_internal_note_2022_2026_L10/figures/Fig5_j_preshower0.png) 中 2025 绿点分裂为 `0–3.5 pC` 与 `5.2 pC` 两支）；
    - 基线窗口为 `0–120 ns`，因此 `Preshower0/1_baseline`（`1891.38 / 1909.22 mV`）与 `baseline_rms`（`0.584 mV`）在 `20012 p0013 vs p0014` 35,021 个同事例上 **`max |Δ| = 0.000`**；
    - `Preshower2/3`（`ch 24, 25`）因 `-35 ns` 对应窗口 `[124, 324]` ns，在 `p0013` 中保持 **`100%` 正常**。
- **全量 2025 普查结果（199 个 GRL run）**：
  - 全部 **`67 / 67`** 个非零 `p0013` run（`19413–20142`）的 `Preshower0 TriggerOffset` 均为 **`+18.0 ns`**、`GFIT_FAILED` 均为 **`84.5%–92.0%`**；
  - 全部 **`94 / 94`** 个非零 `p0014` run（`19336, 19878, 20012, 20170–20677`）的 `TriggerOffset` 均为 **`-40.0 ns`**、`GFIT_FAILED` 最高仅 **`0.67%`**；
  - 另有 `19` 个早期 `p0013` run（`19116–19411` + `19336_p0013`）及 `2` 个 `p0014` run（`20281–20282`）的 `distanceToCollidingBCID = -3564`（`Colliding` 归零）；以及年末 4 个 `p0014` run（`20667, 20668, 20675, 20677`）`VetoNu0` 有效率为 `0.0%`。

![Fig DQ2: Root cause mechanism and all-2025 survey](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/dq_audit_2025_p0013/figures/figDQ2_p0013_preshower_window_bug_and_all2025.png)

### 3.3 全量 `FASER_DQ_Analyis` 逐项审计结果（16 组全 run 数据集）

#### 表 1：事例数、文件完整性、Cutflow 与 `m ± 5s` 产额检查（2025 阈值：Events/pb⁻¹ ∈ `[5502, 12455]`，GoodTracks/pb⁻¹ ∈ `[4641, 10330]`）

| Run & Tag | 有效/总文件 | Lumi [pb⁻¹] | 输入总事例 | `GoodTimes` | `Colliding` | 选中事例 (`Yield`) | `GoodTrkYield` | `Evt/pb⁻¹` | `Trk/pb⁻¹` | `Pos/Neg` | `CaloHi/pb⁻¹` | `CaloLo/pb⁻¹` | DQ 产额判定 |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `19415_p0013` | `11/11` | `769.2` | `6,762,511` | `100.0%` | `99.5%` | `6,728,179` | `5,614,107` | `8,747.5` | `7,299.0` | `0.7923` | `5.59` | `64.7` | `normal` |
| `19524_p0013` | `6/6` | `403.1` | `3,766,024` | `100.0%` | `99.5%` | `3,746,764` | `3,110,267` | `9,295.4` | `7,716.3` | `0.7745` | `3.69` | `64.0` | `normal` |
| `19866_p0013` | `1/1` | `12.2` | `122,094` | `100.0%` | `98.7%` | `120,451` | `100,238` | `9,864.1` | `8,208.8` | `0.7831` | `5.08` | `70.2` | `normal` |
| `19923_p0013` | `14/14` | `989.2` | `8,597,147` | `100.0%` | `99.6%` | `8,561,787` | `7,130,862` | `8,655.6` | `7,209.0` | `0.7752` | `5.59` | `62.6` | `normal` |
| **`19970_p0013`** | **`1/3`** | `213.3` | `202,672` | `100.0%` | `99.5%` | `201,692` | `167,618` | **`945.6`** | **`785.9`** | `0.7700` | **`0.66`** | **`6.9`** | **`incomplete_input`** |
| `20081_p0013` | `14/14` | `949.5` | `7,523,653` | `100.0%` | `99.7%` | `7,498,735` | `6,246,639` | `7,897.6` | `6,578.9` | `0.7741` | `5.01` | `57.3` | `normal` |
| **`19336_p0013`** | `12/12` | `816.5` | `12,852,713` | `100.0%` | **`0.0%`** | **`0`** | **`0`** | **`0.0`** | **`0.0`** | — | `0.00` | `0.0` | **`low` (`BCID=-3564`)** |
| `19878_p0013` | `15/15` | `1051.6` | `9,322,960` | `100.0%` | `99.5%` | `9,278,813` | `7,729,757` | `8,823.7` | `7,350.6` | `0.7777` | `5.62` | `64.4` | `normal` (缺末片 `00750-00772`) |
| `20012_p0013` | `7/7` | `478.8` | `3,396,570` | `100.0%` | `99.5%` | `3,380,631` | `2,808,496` | `7,060.7` | `5,865.8` | `0.7717` | `4.59` | `50.3` | `normal` (缺 `00300-00349` 等片) |
| `19336_p0014` | `12/12` | `816.5` | `7,318,163` | `100.0%` | `99.3%` | `7,263,020` | `6,050,709` | `8,895.0` | `7,410.3` | `0.8375` | `5.76` | `66.3` | `normal` |
| `19878_p0014` | `16/16` | `1051.6` | `9,634,295` | `100.0%` | `99.6%` | `9,589,284` | `7,991,309` | `9,118.9` | `7,599.3` | `0.8045` | `5.85` | `66.7` | `normal` |
| `20012_p0014` | `8/8` | `478.8` | `4,426,415` | `100.0%` | `99.5%` | `4,405,745` | `3,661,767` | `9,201.7` | `7,647.9` | `0.7981` | `5.96` | `65.6` | `normal` |
| `20170_p0014` | `16/16` | `1033.4` | `9,559,689` | `100.0%` | `99.7%` | `9,526,434` | `7,942,562` | `9,218.5` | `7,685.8` | `0.7960` | `5.94` | `66.7` | `normal` |
| `20279_p0014` | `14/14` | `940.5` | `8,695,329` | `100.0%` | `99.7%` | `8,663,646` | `7,219,550` | `9,211.8` | `7,676.3` | `0.7976` | `7.18` | `66.6` | `normal` |
| `20481_p0014` | `9/9` | `734.4` | `4,492,651` | `100.0%` | `99.7%` | `4,476,854` | `3,726,851` | `6,096.4` | `5,075.0` | `0.7853` | `3.87` | `43.3` | `normal` |
| `20591_p0014` | `10/10` | `624.5` | `5,697,626` | `100.0%` | `99.7%` | `5,682,276` | `4,734,213` | `9,098.6` | `7,580.5` | `0.8031` | `5.82` | `66.1` | `normal` |

#### 表 2：全子系统物理量对照（`p0013` 焦点 6 run + 同数据 `19878 / 20012` A/B vs `p0014` 基准 4 run）

| 子系统 / DQ 物理量 | 6 个焦点 `p0013` run | `19878 / 20012` (`p0013`) | `19878 / 20012` (`p0014`) | 4 个 `p0014` 基准 run | 结论 |
| --- | --- | --- | --- | --- | --- |
| **`Preshower0/1` OK% (`GFIT%`)** | **`3.7–6.5%`** (`84.4–92.0%`) | **`3.8–4.0%`** (`91.5–91.8%`) | **`97.0%`** (`0.6%`) | **`97.0–97.1%`** (`0.5–0.7%`) | **严重损坏（COOL Offset `+18 ns`）** |
| **`Preshower0` charge / raw_q [pC]** | **`1.36–2.26`** / `2.56–3.46` | **`2.26`** / `3.16–3.46` | **`5.56`** / `6.16` | **`5.56`** / `6.16` | **脉冲截断（丢前半峰）** |
| `Preshower0` baseline [mV] | `1891.31–1891.51` | `1891.38–1891.40` | `1891.38–1891.40` | `1891.33–1891.40` | 正常（`0–120 ns` 未受影响） |
| `Preshower2/3` (Dig2) OK% / `dt_VN` | `65.2–66.6%` / `17.5 ns` | `66.0–66.5%` / `17.5 ns` | `66.9–67.3%` / `17.5 ns` | `64.3–66.8%` / `17.5 ns` | **正常**（`-35 ns` 窗口覆盖脉冲） |
| `VetoNu0/1` OK% / charge med | `97.5–97.9%` / `99.0 pC` | `97.6–97.8%` / `99.0 pC` | `97.6–97.8%` / `99.0 pC` | `97.6–97.8%` / `99.0 pC` | **正常且完全一致** |
| `Veto10/11, Veto20/21` OK% / charge | `96.1–97.3%` / `87.0 pC` | `96.1–97.2%` / `87.0 pC` | `96.1–97.2%` / `87.0 pC` | `96.1–97.2%` / `87.0 pC` | **正常且完全一致** |
| `Timing0–3` OK% / charge med | `61.6–62.4%` (`50.7–51.4%`) / `21 pC` | `61.9–62.2%` / `21 pC` | `61.9–62.2%` / `21 pC` | `61.8–62.3%` / `21 pC` | **正常且完全一致** |
| `CaloLo` / `CaloHi` OK% / `E_EM` | `29.7–30.1%` / `0.94–0.97` / `223–265 GeV` | `29.8–29.9%` / `0.96–0.97` / `243–248 GeV` | `29.9–30.0%` / `0.96–0.97` / `240–247 GeV` | `29.7–30.2%` / `0.96–0.97` / `246–252 GeV` | **正常**（`19524` `CaHi/pb` 略低 `3.69`） |
| **`Track_Chi2` 均值** | **`17.28–17.36`** | **`17.30–17.33`** | **`16.56–16.60`** | **`16.51–16.57`** | **轻微系统差（`TRACKER-ALIGN-04 -> 05`）** |
| **`Pos/Neg` / `θ_y0` / `Y_atTrig`** | `0.770–0.792` / `-0.06 mrad` / `-1.0 mm` | `0.772–0.778` / `-0.06 mrad` / `-1.1 mm` | `0.798–0.805` / `-0.03 mrad` / `-0.4 mm` | `0.785–0.803` / `-0.03 mrad` / `-0.3 mm` | **同上对齐更新（`+0.026`、`+0.03 mrad`、`+0.7 mm`）** |
| `Track_nLayers` / `nHits` / `p0` med | `8.597` / `16.62` / `110.0 GeV` | `8.597` / `16.61` / `110.0 GeV` | `8.599` / `16.62` / `110.0 GeV` | `8.599` / `16.62` / `110.0 GeV` | **完全一致** |

![Fig DQ3: Full FASER_DQ_Analyis multi-subsystem audit](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/dq_audit_2025_p0013/figures/figDQ3_faser_dq_all_subsystems_check.png)

---

## 4. 阶段性结论

### 已证实的事实
1. **在 `ps_dig1_vetonu`（`fig13_04`）中确有体现**：6 个抽样 `p0013` run 在 `[22, 42]` ns 贡献了 **`88,895` 事例（占 2025 曲线的 `1.30%`）** 的假次峰（中心 `30–34 ns`），因线性纵轴下仅为 `1/470` 主峰高度而被压缩在零线附近，对数纵轴（`figDQ1`）与 `FASER_DQ_Analyis` `Fig5_j/k` 中清晰可见。
2. **不是硬件故障，而是 100% 覆盖全部 `67/67` 个非零 `p0013` run 的生产 COOL 数据库配置缺陷**：`p0013` 使用 `DBRelease/1.6.0` 中过期的 `Preshower0/1 TriggerOffset = +18 ns`，使波形拟合窗口 `[178, 378]` ns 晚于真实脉冲主峰，切掉了前半峰并导致 `84.4%–92.0%` 的 `GFIT_FAILED`；`p0014` 将其更正为 `-40 ns`（窗口 `[120, 320]` ns）后，全部 `94/94` 个 `p0014` run 恢复至 `97.0%` 正常率与 `5.56 pC` MIP 峰。
3. **完整 `FASER_DQ_Analyis` 检查定位的其余 3 类差异**：
   - **`19970_p0013` 文件不完整**：4 个 PHYS slice 中 `0000` 未生成、`0001/0002` 因 `019970-00077-p0013-xAOD.root` 损坏崩溃生成 `555` 字节空壳，仅剩 `0003`（`167,618` GoodTracks，`785.9 trk/pb⁻¹`）；其余 5 个焦点 run（`19415, 19524, 19866, 19923, 20081`）文件完整、产额均在 `m ± 5s` 正常带内。
   - **全部其他闪烁体与量能器通道 100% 正常**：`VetoNu0/1`、`Veto10/11`、`Veto20/21`、`Timing0–3`、`Preshower2/3`、`CaloLo/CaloHi` 的状态正常率、电荷、时间与基线在 `p0013` 与 `p0014` 间高度一致（仅 `19524_p0013` 的 `PS0` 振幅略低导致 `THR_FAILED` 为 `8.3%` 而 `GFIT_FAILED` 略低至 `84.4%`）。
   - **Tracker 对齐更新差异**：`p0013`（`TRACKER-ALIGN-04`）相比 `p0014`（`TRACKER-ALIGN-05`）存在可观测的系统偏移：`Track_Chi2` 均值偏高 `+0.74`（`17.3 vs 16.6`），`Pos/Neg` 偏低 `-0.026`，`θ_y0` 偏 `-0.03 mrad`，`Y_atTrig` 偏 `-0.7 mm`。此外早期 19 个 `p0013` run（`19116–19411` + `19336_p0013`）缺失 `LHC/BCIDData`（`distanceToCollidingBCID = -3564`）。

---

## 5. 下一步计划

1. 任何涉及 `Preshower0/1` 电荷或时间、以及精确 tracker 对齐（`TRACKER-ALIGN-05`）的 2025 物理分析，均应优先使用 `p0014`（或待重处理的 `r0023`）生产版本；若必须读取存量 `p0013` 数据，需依赖未受影响的 `Preshower2/3`（Dig2 `ch 24, 25`，配合 `PreshowerDig2ExtraOffset = 3.97 ns`）作为 Preshower fallback，并排除已确认文件残缺的 `19970` 等 `incomplete_input` run。
2. `2025` 年末 4 个 `p0014` run（`20667, 20668, 20675, 20677`）实测 `VetoNu0_ok = 0.0%`，后续做 `VetoNu` 相关分析时应单列核查其高压/运行状态。
