# Stage 04 — 2022–2026 Data 与 MC24 多时段 Fig. 12 风格波形对照与 Reduced Charge 积分时段优化研究

- **日期**：2026-11-20
- **关联 Calypso 构建**：`/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1`（branch `fix/latest-seed-and-momentum-v1`，commit `1c0f3f45`；两级 Preshower 回退 `Dig1 -> Dig2`）
- **参考规范**：`/afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/docs/Electronic_2D_measurement_internal_note-6.pdf`（p. 21 Fig. 12、p. 22 Fig. 13）
- **分析与绘图脚本**：
  - 数据扫描与窗口优化：[fig12_waveform_and_window_study.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/fig12_waveform_and_window_study.py)
  - 绘图脚本：[plot_fig12_waveforms_and_window.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/plot_fig12_waveforms_and_window.py)
- **产物目录**：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/`
  - 完整合并报告 PDF：[fig12_all_periods_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_all_periods_report.pdf)
  - 冻结统计 JSON：[fig12_window_summary.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/fig12_window_summary.json)

---

## 阶段目标

1. 严格仿照 `Electronic_2D_measurement_internal_note-6.pdf` 第 21 页 **Fig. 12** 的版式（共享时间轴的三站堆叠波形、2 ns 采样点与连线、Veto / FASERν Veto / 预期时间竖线、$\Delta\bar{t}$ 箭头、灰色 `+35 ns` 最大积分区间），基于真实 `xAOD` 波形数据绘制 2022、2023、2024（pre / early / during CaloNu `status=0` 与 `status=512` / post-CaloNu）、2025、2026 及 MC24（Muon + GENIE FaserNu + GENIE CaloNu）全部时段的单事例波形图。
2. 明确列出每个时段应使用的 `VetoVetoNuTimeDifference`、`PreshowerVetoNuTimeDifference`、`PreshowerDig2ExtraOffset` 及对应 `NtupleDumper` 配置预设，并在波形图上逐一给出物理与读出时序对应关系。
3. 在全部时段 51,122 个真实通过 μ 子与 1,503 个 GENIE 中微子事例上定量检验 `ReducedChargeIntegrationTime`（标称 `35 ns`）以及积分下界、基线、双边截断、Dig1/Dig2 差异在不同年份和 MC 中是否需要调整或附加额外处理。

---

## 工作内容

1. **真实 xAOD 多进程波形扫描**（11 worker 并行，Athena `24.0.41` PyROOT 直接读取 `xAOD::WaveformHitContainer` 与 `xAOD::WaveformClock`）：
   - 扫描 11 份覆盖 2022–2026 与 MC24 的真实 `xAOD` 文件（总计 231,000 事例，筛选出 **51,122** 个全站干净通过 μ 子 + **1,503** 个 GENIE 中微子事例，其中 571 个具明显 VetoNu 延迟背散射 $Q_{\mathrm{raw}} \ge 40\ \mathrm{pC}$），全程耗时 `391.1 s`。
   - 对每个时段提取：
     - 具有完整 3/4 站波形（`VetoNu ch 4/5` + `Veto2 ch 6/7` 或 `22/23` + `Preshower Dig1 ch 12/13` [+ 2025–2026 `Dig2 ch 24/25`]）的典型通过 μ 子事例；2024 CaloNu 额外提取 `Veto2 status == 512` 真实事例；MC24 GENIE FaserNu 与 CaloNu 提取延迟背散射被 `+35 ns` 窗口拯救（$Q_{\mathrm{raw}} \ge 40\ \mathrm{pC} \to Q_{\mathrm{red}} < 30\ \mathrm{pC}$）的典型事例。
     - 系综平均的对齐 `VetoNu` 脉冲形状（以 $t - t_{\mathrm{exp}}$ 对齐，计算 `FWHM`、Crystal-Ball 拟合宽度 $\sigma_{\mathrm{CB}}$、峰值相对时间 $t_{\mathrm{peak}} - t_{40\%}$、脉冲前读出窗口长度 $t_{\mathrm{exp}} - t_{\mathrm{start}}$、脉冲前噪声电荷 $Q(t < t_{\mathrm{exp}} - 15\ \mathrm{ns})$、脉冲后最小电压下冲 $\min V([+20, +120]\ \mathrm{ns})$）。
     - 积分上限 $\Delta T \in \{15, 20, 25, 30, 35, 40, 45, 50, 60, 80\}\ \mathrm{ns}$ 的全网格扫描，同步对比：单边窗口 $t < t_{\mathrm{exp}} + \Delta T$（现行算法）与双边窗口 $t \in [t_{\mathrm{exp}} - 20, t_{\mathrm{exp}} + \Delta T]\ \mathrm{ns}$，以及 `Veto2 (fallback=1)`、`Preshower Dig1 (fallback=2)`、`Preshower Dig2 (fallback=3)` 三级参考源的 μ 子漏放率 $P(Q_{\mathrm{red}} < 30\ \mathrm{pC} \mid Q_{\mathrm{raw}} > 40\ \mathrm{pC})$ 与 GENIE 背散射拯救率 $P(Q_{\mathrm{red}} < 30\ \mathrm{pC} \mid Q_{\mathrm{raw}} \ge 40\ \mathrm{pC})$。
2. **Fig. 12 风格成套出图与目视核验**：
   - 产出 7 张多时段波形对照图（`fig12_01`–`fig12_07`）与 2 张定量脉冲形状及窗口优化图（`fig12_08`–`fig12_09`），同时导出 300 dpi PNG、逐图 PDF、9 页汇报 PDF 及 `figure_index.json`，并逐图目视检查确认无标题裁切、无刻度重叠。

执行命令：

```bash
source /eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/setup_environment.sh
python3 /eos/home-x/xcheng/FASER/geoRotation_test/analysis/fig12_waveform_and_window_study.py

source /afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/setup.sh
python3 /eos/home-x/xcheng/FASER/geoRotation_test/analysis/plot_fig12_waveforms_and_window.py
```

---

## 详细数据

### 1. 全时期 Reduced Charge 时序参数与波形对应总表

下表列出 2022–2026 Data 与 MC24 每个时期的通道配置、读出窗口、应应用的 `NtupleDumper` 预设与精确参数值、本批真实 xAOD 实测时间差以及对应波形图编号：

| 时段（代表 Run / DSID） | 配置预设 (`--reduced-charge-preset`) | VetoNu 读出窗口 `[t_start, t_end]` [ns] | Veto2 通道 | Preshower 通道 (`Dig1` / `Dig2`) | `VetoVetoNuTimeDifference` [ns] | `PreshowerVetoNuTimeDifference` [ns] | `PreshowerDig2ExtraOffset` [ns] | 实测中位数 $\Delta t_{\mathrm{V2}}$ / $\Delta t_{\mathrm{PS1}}$ [ns] | 建议积分窗 `ReducedChargeIntegrationTime` [ns] | 对应波形图 |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | --- | ---: | --- |
| **2022** (`008090`, `r0021`) | `data_default` | `[744, 944]`（前余量 `42.5 ns`） | `6, 7` (Dig1) | `12, 13` / — | **33.165** | **53.625** | 0.00 | `32.83 ± 4.11` / `53.52 ± 1.87` | **35.0** | `fig12_01` (左) |
| **2023** (`011703`, `r0021`) | `data_default` | `[144, 344]`（前余量 `42.4 ns`） | `6, 7` (Dig1) | `12, 13` / — | **33.165** | **53.625** | 0.00 | `32.88 ± 1.85` / `53.32 ± 0.59` | **35.0** | `fig12_01` (右) |
| **2024 pre-CaloNu** (`015387`, `r0022`) | `data_default` | `[124, 324]`（前余量 `23.4 ns`） | `6, 7` (Dig1) | `12, 13` / — | **33.165** | **53.625** | 0.00 | `33.10 ± 0.71` / `53.53 ± 0.63` | **35.0** | `fig12_02` (左) |
| **2024 early-CaloNu** (`015838`, `r0022`, runs `15821–15844`) | `data_default` | `[124, 324]`（前余量 `24.2 ns`） | `6, 7` (Dig1) | `12, 13` / — | **33.165** | **53.625** | 0.00 | `33.03 ± 0.83` / `53.48 ± 1.64` | **35.0** | `fig12_02` (右) |
| **2024 during-CaloNu** (`015849`, `r0022`, runs `15845–16540`) | `data_caloNu` | `[124, 324]`（前余量 `24.2 ns`） | **`22, 23` (Dig2)**（22.3% 需 `+16 ns`） | `12, 13` / — | **35.783** | **53.512** | 0.00 | `35.76 ± 0.73` / `53.48 ± 0.63` | **35.0** | `fig12_03` (左/右) |
| **2024 post-CaloNu** (`016929`, `r0022`) | `data_default` | `[124, 324]`（前余量 `24.4 ns`） | `6, 7` (Dig1) | `12, 13` / — | **33.165** | **53.625** | 0.00 | `33.05 ± 1.62` / `53.50 ± 2.53` | **35.0** | `fig12_04` (左) |
| **2025** (`020481`, `p0014`；`p0013` 需 `p0014` 重建) | `data_2025` | `[124, 324]`（前余量 `26.2 ns`） | `6, 7` (Dig1) | `12, 13` / **`24, 25`** | **33.165** | **13.27** | **3.97** | `33.08 ± 2.77` / `13.30 ± 2.99`（PS2 对齐后 `13.27 ± 6.26`） | **35.0** | `fig12_04` (右), `fig12_05` (左) |
| **2026** (`022404`, `p0015`) | `data_2026` | `[124, 324]`（前余量 `25.6 ns`） | `6, 7` (Dig1) | `12, 13` / **`24, 25`** | **33.165** | **13.27** | **3.64** | `33.11 ± 1.02` / `13.34 ± 1.86`（PS2 对齐后 `13.69 ± 4.65`） | **35.0** | `fig12_05` (右) |
| **MC24 Muon** (`100044`) | `mc` | `[780, 980]`（前余量 `12.2 ns`） | `6, 7` | `12, 13` / — | **27.21** | **51.41** | 0.00 | `27.215 ± 0.073` / `51.392 ± 0.593` | **35.0** | `fig12_06` (左) |
| **MC24 GENIE FaserNu** (`200093`) | `mc` | `[780, 980]` | `6, 7` | `12, 13` / — | **27.21** | **51.41** | 0.00 | 共享 MC24 模拟时序模板 | **35.0** | `fig12_06` (右), `fig12_07` (左) |
| **MC24 GENIE CaloNu** (`200139`) | `mc` | `[780, 980]` | `6, 7` | `12, 13` / — | **27.21** | **51.41** | 0.00 | 共享 MC24 模拟时序模板 | **35.0** | `fig12_07` (右) |

> **注**：2022–2026 全部 51,122 个扫描事例中，`VetoNu (ch 4, 5)` 出现 `status & 512` 的计数均为 **0**，说明 VetoNu 始终位于主数字化仪（Dig1），不会触发跨卡 16 ns 跳变。

---

### 2. 各时段 Fig. 12 风格单事例波形图与参数逐一对应

#### (1) 2022 vs 2023：`NominalTriggerTime` 整体平移 `-600 ns`，但站间相对时间差不变

- 2022（run `8090`, entry `26287`）：`VetoNu = 793.7 ns`，`Veto2 = 826.8 ns`（回推 `exp_V2 = 793.6 ns`，误差 `0.1 ns`），`Preshower = 847.2 ns`（回推 `exp_PS = 793.6 ns`）；`Raw Q = 87.0 pC -> Red Q(35 ns) = 75.5 pC (87%)`。
- 2023（run `11703`, entry `2389`）：触发窗口由 `[744, 944]` 平移至 `[144, 344] ns`；`VetoNu = 187.3 ns`，`Veto2 = 220.5 ns`，`Preshower = 241.0 ns`，用完全相同的 `33.165 / 53.625 ns` 回推均精准落在 `187.3–187.4 ns`；`Raw Q = 94.1 pC -> Red Q(35 ns) = 84.3 pC (90%)`。

![2022 vs 2023 Fig.12 waveforms](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_01_2022_vs_2023_waveforms.png)

#### (2) 2024 pre-CaloNu vs 2024 early-CaloNu（runs `15821–15844`）：Veto2 仍在 Dig1 `ch 6,7`

- 2024 读出窗口前移至 `[124, 324] ns`（脉冲前余量由 `42 ns` 缩至 `23–24 ns`），但 early-CaloNu（run `15838`）在 Veto2 迁至 Dig2 之前仍保持 `Veto2 ch 6, 7` 与 `33.165 / 53.625 ns`，`exp_V2` 和 `exp_PS` 与实测 `VetoNu` 误差 `< 0.1 ns`。

![2024 pre vs early CaloNu Fig.12 waveforms](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_02_2024_pre_vs_early_caloNu_waveforms.png)

#### (3) 2024 during-CaloNu（runs `15845–16540`）：Veto2 迁至 Dig2 `ch 22,23`（`+2.62 ns` 硬件偏移 + `status==512` `+16 ns` 还原）

- 左图（run `15849`, entry `2565`, `status=0`）：`Veto2 ch 22` 到达时间相对 VetoNu 从 `33.165 ns` 延迟至 `35.783 ns`（`+2.618 ns`），而 `Preshower` 仍停在 `Dig1 ch 12,13`（`53.512 ns`）。
- 右图（run `15849`, entry `14522`, **`status=512`**）：原始 `Veto2` 波形（中图灰色虚线）整体提前 `16 ns` 位于 `174.4 ns`；若不加 `+16 ns` 修正，积分截止线会错误前移到 `173.6 ns`（顶图灰点线），直接截掉 μ 子脉冲主体后半段！加上 `if (status & 512) t += 16.0` 后，`Veto2` 回到 `190.4 ns`，减去 `35.783 ns` 得 `exp_V2 = 154.6 ns`，与实测 `VetoNu = 154.58 ns` 精确重合（误差 `0.01 ns`），`Red Q(35 ns) = 74.9 pC (85%)`。

![2024 CaloNu status=0 vs status=512 Fig.12 waveforms](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_03_2024_caloNu_status0_vs_status512.png)

#### (4) 2024 post-CaloNu vs 2025 & 2026：升级 4-PMT Preshower 导致 `PreshowerVetoNu` 跳变 `-40.3 ns` 并出现 Dig1/Dig2 双卡时间差

- **fig12_04**（2024 post-CaloNu run `16929` vs 2025 run `20481 p0014`）：2024 post-CaloNu Veto2 迁回 Dig1，恢复 `33.165 / 53.625 ns`；2025 升级 4-PMT Preshower 后，`Preshower Dig1 ch 12,13` 的信号时间骤降 `-40.3 ns`，直接跳到 **`Veto2` 之前约 20 ns**（`t_PS1 = 165.4 ns < t_V2 = 185.2 ns`），对应 `PreshowerVetoNuTimeDifference` 必须由 `53.625 ns` 更新为 **`13.27 ns`**，否则 Preshower 回退会将积分窗整体左移 40 ns、彻底漏掉 VetoNu 脉冲。
- **fig12_05**（2025 run `20481` vs 2026 run `22404`）：同屏叠加紫色 `Preshower Dig2 ch 24,25` 波形，清晰可见两处硬件特征：(a) `Dig2` 相对 `Dig1` 晚 `+3.97 ns`（2025）/ `+3.64 ns`（2026），需要 `PreshowerDig2ExtraOffset` 对齐；(b) `Dig2` 通道幅度显著低于 `Dig1`（约 3–4 mV vs 20 mV），信噪比低，这正是我们在 `1c0f3f45` 中改为“优先 `Dig1 (fb=2)`、失效才回退对齐后的 `Dig2 (fb=3)`”的直接波形证据。

![2024 post-CaloNu vs 2025 Fig.12 waveforms](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_04_2024_post_vs_2025_preshower_shift.png)

![2025 vs 2026 Fig.12 waveforms](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_05_2025_vs_2026_4pmt_waveforms.png)

#### (5) MC24 Muon 与 GENIE 中微子（FaserNu & CaloNu）：标准 Fig. 12 物理机制对照

- **fig12_06**（左：MC24 Muon `100044`, entry `4235`；右：MC24 GENIE FaserNu `200093`, entry `2334`）：MC24 使用 `27.21 / 51.41 ns`。左图通过 μ 子脉冲完全落在 `[780, 827.2] ns` 灰色窗口内（`Raw Q = 95.9 pC -> Red Q = 80.9 pC, 84%`）；右图中微子事例在 `836 ns` 出现主背散射峰（位于 `828.8 ns` 截止线之后），被 `+35 ns` 窗口成功截除（`Raw Q = 56.3 pC -> Red Q = 14.8 pC, 26%, rescued <30 pC`）。
- **fig12_07**（左：GENIE FaserNu `200093`；右：GENIE CaloNu `200139`, entry `6983`）：CaloNu 几何下延迟背散射同样在 `818.1 ns` 才过 40% 上升沿（相对 `exp_V2 = 792.2 ns` 延迟 `+25.9 ns`），`+35 ns` 窗口将其在峰顶处截断（`Raw Q = 61.8 pC -> Red Q = 20.5 pC, 33%, rescued <30 pC`）。

![MC24 Muon vs GENIE FaserNu Fig.12 waveforms](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_06_mc24_muon_vs_genie_fasernu.png)

![MC24 GENIE FaserNu vs CaloNu Fig.12 waveforms](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_07_mc24_genie_fasernu_vs_calonu.png)

---

### 3. VetoNu 脉冲形状、脉冲前基线与下冲跨时期检验

对 9 组 μ 子样本（共 51,122 事例）的 VetoNu 波形几何特征定量汇总如下：

| 时段 | 扫描 / 选中 μ 子数 | 前窗口长度中位数 (`q05`) [ns] | 脉冲 `FWHM` [ns] | CB 宽度 $\sigma$ [ns] | $t_{\mathrm{peak}} - t_{40\%}$ [ns] | 脉冲前电荷 $Q(t < t_{\mathrm{exp}}-15)$ 中位数 (`q99`) [pC] | $Q_{\mathrm{pre15}} > 5\ \mathrm{pC}$ 占比 | 脉冲后最小电压 `[+20, +120]` [mV] | 原始总电荷中位数 [pC] | `@35 ns` 电荷俘获中位数 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2022 (`8090`) | 25,000 / 166 | 42.46 (33.71) | **20.0** | 6.00 | 7.84 | −0.003 (0.13) | 0.000% | +0.46 | 92.8 | **83.66%** |
| 2023 (`11703`) | 25,000 / 401 | 42.35 (34.72) | **20.0** | 6.03 | 8.14 | +0.008 (0.14) | 0.000% | +0.53 | 100.3 | **84.92%** |
| 2024 pre (`15387`) | 25,000 / 6,118 | 23.36 (16.16) | **20.0** | 6.08 | 8.00 | +0.000 (0.08) | 0.020% | +0.49 | 100.1 | **85.06%** |
| 2024 early (`15838`) | 25,000 / 6,949 | 24.20 (17.09) | **20.0** | 6.04 | 8.02 | +0.000 (0.08) | 0.030% | +0.48 | 100.7 | **85.13%** |
| 2024 CaloNu (`15849`) | 30,000 / 8,901 | 24.19 (17.09) | **20.0** | 6.03 | 7.97 | +0.001 (0.08) | 0.010% | +0.48 | 100.6 | **85.26%** |
| 2024 post (`16929`) | 25,000 / 6,893 | 24.42 (17.13) | **20.0** | 6.03 | 7.99 | +0.001 (0.09) | 0.000% | +0.47 | 100.6 | **85.23%** |
| 2025 (`20481`) | 25,000 / 6,419 | 26.17 (19.03) | **20.0** | 6.06 | 8.05 | +0.001 (0.09) | 0.050% | +0.47 | 100.2 | **85.53%** |
| 2026 (`22404`) | 25,000 / 9,770 | 25.55 (18.20) | **20.0** | 6.02 | 7.98 | +0.000 (0.09) | 0.010% | +0.29 | 100.2 | **85.63%** |
| **MC24 Muon (`100044`)** | 6,000 / 5,505 | **12.16 (12.15)** | **18.0** | **5.06** | **5.85** | +0.000 (0.00) | 0.000% | +0.10 | 101.8 | **84.28%** |

关键物理观测（见 **fig12_08**）：

1. **2022–2026 真实数据 VetoNu 脉冲形状高度不变**：七个真实数据段的归一化平均波形完全重合（`FWHM = 20.0 ns`，$\sigma_{\mathrm{CB}} = 6.00\text{–}6.08\ \mathrm{ns}$，峰值滞后 $t_{40\%}$ 恒为 `7.8–8.1 ns`）。
2. **MC24 脉冲模板略窄但在 `35 ns` 处与数据高度吻合**：MC24 模拟脉冲的 `FWHM = 18.0 ns`（$\sigma_{\mathrm{CB}} = 5.06\ \mathrm{ns}$，上升更快），但尾部稍缓，使得在 $\Delta T = 35\ \mathrm{ns}$ 处中位电荷俘获率为 `84.28%`，恰好位于真实数据 `83.66%–85.63%` 区间之内（差异 $\le 1.3\%$）。
3. **脉冲前窗口长度变化 3.5 倍但对单边积分无基线漂移污染**：尽管脉冲前长度从 `42.5 ns`（2022–2023）变为 `23–26 ns`（2024–2026）再到 `12.2 ns`（MC24），所有时期脉冲前 15 ns 之外的累积电荷中位数均为 `-0.003..+0.008 pC`（99% 分位 $\le 0.14\ \mathrm{pC}$，$> 5\ \mathrm{pC}$ 的早到堆积占比 $\le 0.05\%$），脉冲后无任何负电压下冲（最小电压 `+0.10..+0.53 mV`）。

![VetoNu pulse shape and charge capture across all periods](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_08_vetonu_pulse_shape_and_capture_all_periods.png)

---

### 4. 积分窗口 $\Delta T$ 扫描：μ 子漏放率 vs GENIE 中微子背散射拯救率

下表汇总不同积分上限 $\Delta T$ 下各时段通过 μ 子的误漏放率 $P(Q_{\mathrm{red}} < 30\ \mathrm{pC} \mid Q_{\mathrm{raw}} > 40\ \mathrm{pC})$，以及 MC24 GENIE 中微子背散射事例的拯救率 $P(Q_{\mathrm{red}} < 30\ \mathrm{pC} \mid Q_{\mathrm{raw}} \ge 40\ \mathrm{pC})$（见 **fig12_09**）：

| 样本 / 参考路径 | `15 ns` | `20 ns` | `25 ns` | `30 ns` | **`35 ns`（标称）** | `40 ns` | `45 ns` | `50 ns` | `60 ns` | `80 ns` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **2022 Veto2 (`fb=1`) 漏放** | 15.060% | 7.229% | 1.205% | 0.000% | **0.000%** | 0.000% | 0.000% | 0.000% | 0.000% | 0.000% |
| **2023 Veto2 (`fb=1`) 漏放** | 2.743% | 0.998% | 0.000% | 0.000% | **0.000%** | 0.000% | 0.000% | 0.000% | 0.000% | 0.000% |
| **2024 pre / early / CaloNu Veto2 漏放** | 0.17–0.25% | 0.03–0.08% | 0.000% | 0.000% | **0.000%** | 0.000% | 0.000% | 0.000% | 0.000% | 0.000% |
| **2024 post Veto2 (`fb=1`) 漏放** | 0.276% | 0.058% | 0.029% | 0.015% | **0.015%** | 0.015% | 0.015% | 0.015% | 0.015% | 0.015% |
| **2025 Veto2 (`fb=1`) / PS Dig1 (`fb=2`) 漏放** | 0.203% / 0.214% | 0.078% / 0.092% | 0.016% / 0.061% | 0.016% / 0.031% | **0.000% / 0.000%** | 0.000% / 0.000% | 0.000% / 0.000% | 0.000% / 0.000% | 0.000% / 0.000% | 0.000% / 0.000% |
| **2025 Veto2 双边 `[-20, ΔT]` 漏放** | 0.234% | 0.109% | 0.047% | 0.047% | **0.031%** | 0.031% | 0.031% | 0.031% | 0.031% | 0.031% |
| **2025 PS Dig2 (`fb=3`) 漏放** | 1.809% | 1.166% | 1.125% | 0.884% | **0.764%** | 0.683% | 0.683% | 0.442% | 0.161% | 0.000% |
| **2026 Veto2 (`fb=1`) / PS Dig1 (`fb=2`) 漏放** | 0.133% / 0.296% | 0.041% / 0.020% | 0.020% / 0.000% | 0.010% / 0.000% | **0.010% / 0.000%** | 0.010% / 0.000% | 0.000% / 0.000% | 0.000% / 0.000% | 0.000% / 0.000% | 0.000% / 0.000% |
| **2026 PS Dig2 (`fb=3`) 漏放** | 1.282% | 1.015% | 0.828% | 0.641% | **0.427%** | 0.321% | 0.321% | 0.187% | 0.080% | 0.000% |
| **MC24 Muon Veto2 / PS1 漏放** | 0.000% | 0.000% | 0.000% | 0.000% | **0.000%** | 0.000% | 0.000% | 0.000% | 0.000% | 0.000% |
| **GENIE FaserNu (`200093`, N=273 背散射) 拯救率（单边 = 双边）** | 69.60% | 60.81% | 55.31% | 47.62% | **42.12%** | 36.26% | 32.97% | 29.67% | 20.51% | 9.16% |
| **GENIE CaloNu (`200139`, N=298 背散射) 拯救率（单边 = 双边）** | 43.96% | 40.60% | 37.25% | 34.23% | **31.54%** | 27.52% | 25.17% | 22.15% | 18.46% | 12.42% |

> **GENIE 总体通过率与背散射延迟**：FaserNu（785 事例，原始 $Q_{\mathrm{raw}} < 40\ \mathrm{pC}$ 占 `65.22%`，延迟背散射中位数 `+17.82 ns`）在 `35 ns` 处将总通过率提高至 **`79.87%`**（拯救 `115/273`）；CaloNu（718 事例，原始占 `58.50%`，延迟背散射中位数 `+11.21 ns`）在 `35 ns` 处将总通过率提高至 **`71.59%`**（拯救 `94/298`）。

![Integration window optimization: muon leakage vs GENIE backsplash rescue](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_09_window_optimization_leakage_vs_rescue.png)

---

### 5. 完整图表产物索引

- 多页合并报告 PDF：[fig12_all_periods_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_all_periods_report.pdf)
- 产物元数据索引：[figure_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/figure_index.json)

| 编号 | 图表内容 | PNG | PDF |
| --- | --- | --- | --- |
| Fig.12-01 | 2022 vs 2023 真实 μ 子波形（触发窗 820 → 220 ns，相同 $\Delta t$） | [fig12_01_2022_vs_2023_waveforms.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_01_2022_vs_2023_waveforms.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_01_2022_vs_2023_waveforms.pdf) |
| Fig.12-02 | 2024 pre-CaloNu vs early-CaloNu 真实 μ 子波形（均在 Dig1 `ch 6,7`） | [fig12_02_2024_pre_vs_early_caloNu_waveforms.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_02_2024_pre_vs_early_caloNu_waveforms.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_02_2024_pre_vs_early_caloNu_waveforms.pdf) |
| Fig.12-03 | 2024 during-CaloNu `Veto2 ch 22,23`：`status=0` vs `status=512`（`+16 ns` 修正） | [fig12_03_2024_caloNu_status0_vs_status512.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_03_2024_caloNu_status0_vs_status512.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_03_2024_caloNu_status0_vs_status512.pdf) |
| Fig.12-04 | 2024 post-CaloNu vs 2025（Preshower 升级跳变 `−40.3 ns` 领先 Veto2） | [fig12_04_2024_post_vs_2025_preshower_shift.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_04_2024_post_vs_2025_preshower_shift.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_04_2024_post_vs_2025_preshower_shift.pdf) |
| Fig.12-05 | 2025 vs 2026 升级 4-PMT Preshower（`Dig1` 与低增益 `Dig2 +3.97/+3.64 ns` 同屏对照） | [fig12_05_2025_vs_2026_4pmt_waveforms.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_05_2025_vs_2026_4pmt_waveforms.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_05_2025_vs_2026_4pmt_waveforms.pdf) |
| Fig.12-06 | MC24 模拟 μ 子 vs GENIE FaserNu 延迟背散射（标准 Fig. 12 对照） | [fig12_06_mc24_muon_vs_genie_fasernu.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_06_mc24_muon_vs_genie_fasernu.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_06_mc24_muon_vs_genie_fasernu.pdf) |
| Fig.12-07 | MC24 GENIE FaserNu vs CaloNu 中微子背散射对比 | [fig12_07_mc24_genie_fasernu_vs_calonu.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_07_mc24_genie_fasernu_vs_calonu.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_07_mc24_genie_fasernu_vs_calonu.pdf) |
| Fig.12-08 | 全时段 VetoNu 对齐脉冲形状与累积电荷俘获曲线 | [fig12_08_vetonu_pulse_shape_and_capture_all_periods.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_08_vetonu_pulse_shape_and_capture_all_periods.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_08_vetonu_pulse_shape_and_capture_all_periods.pdf) |
| Fig.12-09 | 积分窗 $\Delta T$ 优化：μ 子漏放率 vs GENIE 背散射拯救率 | [fig12_09_window_optimization_leakage_vs_rescue.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_09_window_optimization_leakage_vs_rescue.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/fig12_waveform_window_study/figures/fig12_09_window_optimization_leakage_vs_rescue.pdf) |

---

## 阶段性结论

### 事实结论（数据直接支持）

1. **所有时期（2022–2026 Data 及 MC24）均应保持统一的 `ReducedChargeIntegrationTime = 35.0 ns`，无需分时期更改积分时长**：
   - 2022–2026 全部真实数据的 VetoNu 脉冲形状完全一致（`FWHM = 20.0 ns`，$\sigma_{\mathrm{CB}} = 6.00\text{–}6.08\ \mathrm{ns}$），`@35 ns` 中位电荷俘获率为 `83.66%–85.63%`，与略窄的 MC24 模拟波形（`FWHM = 18.0 ns`，俘获率 `84.28%`）在 `35 ns` 处高度一致。
   - **不可缩窄至 `< 30 ns`**：2022 与 2023 数据受较低时序分辨与少量早/晚散点影响，在 `25 ns` / `20 ns` / `15 ns` 处 μ 子漏放率分别升至 `1.21% / 7.23% / 15.06%`（2022）和 `0.00% / 1.00% / 2.74%`（2023）；在 `35 ns` 处，全时期 `Veto2 (fb=1)` 与 `Preshower Dig1 (fb=2)` 的 μ 子漏放率均降至 **`0.000%–0.015%`**，且相对 `30 ns` 留有 `+5 ns` 安全余量。
   - **不应放宽至 `> 35–40 ns`**：由于中微子延迟背散射中位延迟仅 `+17.8 ns`（FaserNu）和 `+11.2 ns`（CaloNu），将 $\Delta T$ 从 `35 ns` 放宽至 `40 / 50 / 60 ns` 会使 FaserNu 背散射拯救率从 `42.12%` 跌至 `36.26% / 29.67% / 20.51%`，CaloNu 从 `31.54%` 跌至 `27.52% / 22.15% / 18.46%`，却无法带来任何额外 μ 子拒斥收益。
2. **现行单边积分条件 `t < exp_time + 35.0 ns` 最优，无需增加下界截断或基线重扣**：
   - 尽管不同年份脉冲前读出余量相差 3.5 倍（`42.5 ns` vs `23–26 ns` vs `12.2 ns`），各时期脉冲前基线漂移极小（$Q_{\mathrm{pre15}}$ 中位数 $\le 0.008\ \mathrm{pC}$，$q_{99} \le 0.14\ \mathrm{pC}$），且不存在负电压下冲。
   - 改为双边窗口 `[-20, +35] ns` 对 GENIE 拯救率无任何增益（`42.12%` vs `42.12%`、`31.54%` vs `31.54%` 逐项一致），反而在 2025/2026 将少量参考时间偏晚事例的 μ 子漏放率从 `0.000%` / `0.010%` 微升至 `0.031%` / `0.020%`。
3. **唯一需要特别处理的跨时期效应不在积分窗宽，而在三类硬件/条件变化（已全部由两级回退配置与状态位逻辑覆盖）**：
   - **2024 CaloNu `Veto2 ch 22,23`**：使用 `data_caloNu`（`35.783 / 53.512 ns`）并由 `status & 512` 实时加 `+16.0 ns`；无须硬编码 run 号。
   - **2025–2026 升级 4-PMT Preshower**：使用 `data_2025` / `data_2026`（`33.165 / 13.27 ns`，`Dig2Extra = 3.97 / 3.64 ns`）；低增益 `Dig2 (ch 24,25)` 即使放宽到 `45–50 ns` 仍有 `0.68% / 0.44%`（2025）与 `0.32% / 0.19%`（2026）远尾时间离群点漏放，无法靠调宽积分窗解决，必须依赖已实现的“优先 `Dig1 (fb=2)`、失效才回退对齐 `Dig2 (fb=3)`”两级逻辑，并在物理分析中可按 `VetoNu_reduced_charge_fallback == 3` 单独监控该极少数次级回退事例。
   - **2025 `p0013` 重建产物**：因 COOL `/WAVE/DAQ/TimingOffset` 滞后 `+58 ns` 导致 Preshower 波形拟合窗错位，必须使用 `p0014` 及以上重建（见前序 DQ 审计报告 [2026-11-19_stage04_2025_p0013_preshower_dq_audit.md](file:///eos/home-x/xcheng/FASER/geoRotation_test/workbook/2026-11-19_stage04_2025_p0013_preshower_dq_audit.md)）。

### 待验证假设

- 当前 GENIE 背散射拯救率基于 MC24 模拟波形模板（`FWHM = 18.0 ns`）测得；真实数据中由 2025–2026 束流背景及中微子候选事例产生的延迟背散射时间分布可能与 MC24 存在 1–2 ns 级的脉冲宽度差异，但因 `35 ns` 对 18 ns 与 20 ns 模板均处于 `84–85.6%` 平台区，预期该差异不会改变 `35 ns` 的最优点。

---

## 下一步计划

1. 若后续对 2025–2026 全量数据执行批量 `faser_ntuple_maker.py` 生产，按上表对应 `--reduced-charge-preset`（`data_default` / `data_caloNu` / `data_2025` / `data_2026` / `mc`）保持 `ReducedChargeIntegrationTime = 35.0 ns` 统一输出。
2. 在真实 2025–2026 数据中微子分析选例中，单独抽查 `VetoNu_reduced_charge_fallback == 3`（即 Veto2 与 Preshower Dig1 同时饱和/失效、仅依赖低增益 Dig2 的极少数事例）的波形与电荷分布，核实无需为 `fb == 3` 引入额外电荷系统误差。
