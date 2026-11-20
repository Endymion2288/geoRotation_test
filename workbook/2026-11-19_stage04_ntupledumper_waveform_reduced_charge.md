# 阶段报告：`geoRotation_calypso_latest_fix_v1` 中 2022–2026 数据及 MC 的 Waveform Reduced VetoNu Charge 移植、去硬编码重构与 Fig. 13 真实数据时间差研究

- **日期**：2026-11-19
- **主机与环境**：`lxplus901.cern.ch`（28 cores），Athena `24.0.41`（`x86_64-el9-gcc13-opt`）与 LCG_107（ROOT RDataFrame）
- **基线源码（未改动）**：`/eos/home-x/xcheng/FASER/geoRotation_calypso_latest`，HEAD `7fd40df8cd3b60e81b4095a4f84310c4a54700a0`
- **工作树源码与提交**：`/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1`（分支 `fix/latest-seed-and-momentum-v1`，提交 `7729ddf0` -> 两级 Dig1/Dig2 升级 `1c0f3f45`）
- **可复用补丁**：[0051_ntupledumper_waveform_reduced_charge_2022_2026_mc.patch](file:///eos/home-x/xcheng/FASER/geoRotation_test/patches/0051_ntupledumper_waveform_reduced_charge_2022_2026_mc.patch)
- **编译产物（Overlay）**：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/software_ntuple_fix_v1_clean/`
- **严格仿 Fig. 13 全年份时间差图集（11 页，PNG + 单图 PDF + 汇总 PDF）**：[fig13_all_timing_diffs_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_all_timing_diffs_report.pdf)（索引：[figure_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/figure_index.json)）
- **多面板诊断辅助图集（6 页）**：[timing_study_fig13_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/figures/timing_study_fig13_report.pdf)（索引：[figure_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/figures/figure_index.json)）
- **端到端验证完整图集**：[validation_v2_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v2/figures/validation_v2_report.pdf)（索引：[figure_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v2/figures/figure_index.json)）

---

## 1. 阶段目标

在 `geoRotation_calypso_latest_fix_v1` 的 `PhysicsAnalysis/NtupleDumper` 中正确移植并验证适用于 2022–2026 数据与 MC 的 waveform 及 reduced VetoNu charge 计算，同时满足以下架构与物理验收目标：

1. **彻底消除 C++ 硬编码 run 编号**：严格遵循 Calypso（`NtupleDumperAlg`）与 `/eos/home-x/xcheng/FASER/El9CalypsoForWaveForm`（`MyWaveformTimeRangeAlg`、`MyWaveformTimeRangeConfig.py`）的代码风格。C++ 算法仅依据 `xAOD::WaveformStatus` 状态位（`TIMING_SHIFTED = 512`、`SECONDARY = 2`）和 `IWaveformCableMappingTool`（`m_wave_map`）工作，**不包含任何按 run 编号分支的逻辑**；所有随时期或 MC 变化的站点时间差参数均放在易于维护的 Python 配置与命令行接口（`NtupleDumperConfig.py`、`faser_ntuple_maker.py`）中传入。
2. **阐明各站点时间差偏移的物理与硬件根源**：仿照 [Electronic_2D_measurement_internal_note-6.pdf](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/docs/Electronic_2D_measurement_internal_note-6.pdf) 第 22 页 Fig. 13，在 [FASER_DQ_RDF.py](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/FASER_DQ_RDF.py) 定义的精确 `GoodTracks` 与 Data Quality 选择下，扫描 2022–2026 全时期官方 PHYS 数据及 MC24，绘制 `Veto2 - VetoNu`、`Preshower (Dig1) - VetoNu`、`Preshower (Dig2) - VetoNu`、`Preshower (Dig2 - Dig1)` 及 `Timing / Veto1 - VetoNu` 分布，用真实数据定量说明为什么需要 `VetoVetoNuTimeDifference`、`PreshowerVetoNuTimeDifference` 和 `PreshowerDig2ExtraOffset`。
3. **在官方默认通量 GENIE 中微子 MC 与 2022–2026 真实数据上端到端闭合验证**：
   - 使用 [FASER_FaserMC.pdf](file:///eos/home-x/xcheng/FASER/nu_muAnalysis/docs/FASER_FaserMC.pdf) 中 **Genie Neutrino Samples** 三个章节的全部 9 个 **default flux prediction** run（`2022/23 Crossing angle distribution`：`200082/200083/200089`；`2024 Crossing angle distribution, FaserNu geometry`：`200093/200094/200100`；`2024 Crossing angle distribution, CaloNu geometry`：`200139/200140/200146`）以及穿行 muon MC（`100044`）。
   - 在 2022、2023、2024（pre-CaloNu / early-CaloNu / during-CaloNu / post-CaloNu）、2025、2026 真实 xAOD 上运行编译后的 `faser_ntuple_maker.py`，并与 Ben Wilson 的 5 套参考目录（`/eos/user/b/bewilson/El9CalypsoForWaveForm/`）逐事件对标。

---

## 2. 工作内容

### 2.1 为什么存在 `VetoVetoNuTimeDifference`、`PreshowerVetoNuTimeDifference` 与 `PreshowerDig2ExtraOffset`（物理、硬件与波形重建机制）

逐行审计 [WaveformReconstructionTool.cxx](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Waveform/WaveRecTools/src/WaveformReconstructionTool.cxx)、[WaveformDigitisationTool.cxx](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Waveform/WaveDigiTools/src/WaveformDigitisationTool.cxx)、[RawWaveformDecoderTool.cxx](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Waveform/WaveEventCnv/WaveByteStream/src/RawWaveformDecoderTool.cxx) 及 COOL SQLite 数据库 `CABP200.db` / `ALLP200.db` 后，完整机理如下：

1. **为什么 `reduced_charge` 必须在 `localtime` 坐标系下使用站点时间差（而非 `triggertime`）**：
   - 原始采样波形 `hit->time_vector()` 存储的是 500 MHz 数字化仪原始采样时刻 $t_i = 2i\ \mathrm{ns}$（`WaveformReconstructionTool.cxx:416`）；Crystal Ball 拟合得到的 40% 恒比定时 `hit->localtime()` 同样处于原始采样时刻坐标系，**并未扣除** COOL 数据库 `/WAVE/DAQ/TimingOffset` 中的 `TriggerOffset` 或第 2 数字化仪的 `-16 ns` 偏移（`WaveformReconstructionTool.cxx:331-345` 仅在计算 `trigger_time` 和 `bcid_time` 时扣除 `TriggerOffset + digi_offset`）。
   - 同时，MC 中不填 `triggertime`/`bcidtime`（`NtupleDumperAlg.cxx:255` 跳过 `if (!isMC)`），且真实数据中 `/WAVE/DAQ/TimingOffset` 仅用于粗略对准重建窗口（例如 2025 年 `Preshower` 设 `-40/-42 ns` 导致 `Preshower - VetoNu` 的 `triggertime` 差值反而跳变至 `+18.3 / +24.3 ns`）。
   - 因此，对正穿光速 muon（$\beta \simeq 1$），任意下游参考站 $R$ 与 `VetoNu` 的 `localtime` 差值由四项叠加决定（其中全局 `NominalTriggerTime` 在 2022 的 `820 ns`、2023 的 `220 ns`、2024–2026 的 `200 ns` 间变化，但在同事件两站相减中**精确抵消**）：
     $$\Delta t_{\mathrm{local}}(R - \mathrm{VetoNu}) = \frac{z_R - z_{\mathrm{VetoNu}}}{c} + \left(\tau_{\mathrm{PMT+cable}, R} - \tau_{\mathrm{PMT+cable}, \mathrm{VetoNu}}\right) + \left(\delta_{\mathrm{board}, R} - \delta_{\mathrm{board}, \mathrm{VetoNu}}\right) + \delta_{16\,\mathrm{ns}}(\mathrm{bit\,512})$$
2. **为什么需要 `VetoVetoNuTimeDifference` 区分标准时期（`33.165 ns`）、2024 CaloNu（`35.783 ns`）与 MC24（`27.21 ns`）**：
   - **标准真实数据**（2022、2023、2024 pre-CaloNu、2024 early-CaloNu `15821–15844`、2024 post-CaloNu、2025、2026）：`VetoNu`（`z = -3112 mm`，`ch 4,5`）与 `Veto St2`（`z = -1610 mm`，`ch 6,7`）均位于第 1 块数字化仪（Board 0，通道 `0–15`），光速飞行时间 `+5.0 ns` 加 PMT/电缆延迟差（`ALLP200.db` 标称 `TriggerOffset`: VetoNu `-35 ns`, Veto2 `0 ns`），实测稳定在 **`33.08–33.23 ns`**。
   - **2024 CaloNu 时期**（`15845–16927`）：为腾出 Board 0 通道给 CaloNu 触发，`Veto St2` 迁至第 2 块数字化仪（Board 1，`ch 22, 23`），而 `VetoNu` 留在 Board 0（`ch 4, 5`）。双数字化仪引入**两项独立效应**：
     - **离散 `-16.0 ns` 触发量化跳变**：约 `22.3%` 的事例中两块数字化仪触发时间戳相差 1 个 8-sample 时钟节拍（`wave.trigger_offset() >= 1`），`WaveformReconstructionTool` 将采集窗口前移 8 个点（`-16 ns`）并置位 `xAOD::WaveformStatus::TIMING_SHIFTED`（bit 9 = `512`）。未校正时在 `19.87 ns` 处形成次峰；施加 `if (status & 512) t += 16.0;` 后与 `status == 0` 峰精确重合至 **`0.00 ns`** 以内（见图 `fig02`）。
     - **固定 `+2.7 ns` 跨板/重接电缆延迟**：即使完成 `+16.0 ns` 校正，Board 1 (`ch 22, 23`) 的 `Veto2 - VetoNu` 峰值仍比 Board 0 (`ch 6, 7`) 晚 **`+2.69 ns`**（`33.18 -> 35.87 ns`，Ben Wilson 参考值 `35.783 ns`）。该固定板间硬件延迟不带状态位标记，必须通过 Python 参数配置。
   - **MC24（FaserNu 与 CaloNu 几何一致）**：根据 [makeDigiDB.py](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Waveform/WaveformConditions/WaveCondUtils/scripts/makeDigiDB.py) 与 `WaveformDigitisationTool.cxx:274-295`，MC 数字化在卷积前**已扣除 `tof = z/c`**，并直接按 Crystal Ball 模板中心 `cb_mean`（VetoNu `800 ns, σ=5.0 ns`；Veto `825 ns, σ=3.35 ns`；Timing `830 ns, σ=3.21 ns`；Preshower `850 ns, σ=4.0 ns`）生成波形。40% 上升沿恒比定时近似为 $t_{40\%} \approx \mu_{\mathrm{CB}} - \sqrt{2\ln 2.5}\,\sigma_{\mathrm{CB}} = \mu_{\mathrm{CB}} - 1.354\,\sigma_{\mathrm{CB}}$，解析给出 $\Delta t_{\mathrm{Veto2-VN}} \approx (825-800) - 1.354(3.35-5.0) = 27.23\ \mathrm{ns}$、$\Delta t_{\mathrm{PS-VN}} \approx (850-800) - 1.354(4.0-5.0) = 51.35\ \mathrm{ns}$，与 MC24 实测（`27.28 / 51.36 ns`）及官方常数（`27.21 / 51.41 ns`）严格自洽！
3. **为什么需要 `PreshowerVetoNuTimeDifference` 在 2025–2026 降至 `13.27 ns`**：
   - 在中微子 CC 相互作用中，前向簇射使 **`36%–72%`** 的事例发生 `Veto2` 波形溢出饱和（`status & 4`），从而触发 `Preshower` fallback（`fallback_mode == 2`，占比高达 `29%–66%`）。
   - 2022–2024 年旧 Preshower（`ch 12, 13`，`TriggerOffset = +18 ns`）实测 `t_PS(Dig1) - t_VN = 53.49–53.62 ns`。
   - 2025–2026 年 Preshower 升级（更换新 PMT 与更短信号电缆，并增设底部两路 PMT `ch 24, 25`；`ALLP200.db` 中 `TriggerOffset` 由 `+18 ns` 改为 `-40 / -42 ns`），实测 `Preshower0/1`（Dig1 `ch 12, 13`）绝对 `localtime` 由 `201.7 ns` 前移约 `40 ns` 至 `162.8 ns`，导致 `t_PS(Dig1) - t_VN` 从 `53.62 ns` 骤降 **`-40.3 ns`** 至 **`13.27 ns`**（2025 与 2026 一致至 `0.01 ns`）。若沿用旧值 `53.625 ns`，fallback 积分截止点会提前 `40.3 ns`，直接截断 `VetoNu` 正常 muon 脉冲，将 `88%–97%` 的穿行 muon 误判为 `< 40 pC` 假中微子候选！
4. **为什么 2025–2026 需要 `PreshowerDig2ExtraOffset`（2025: `+3.97 ns`；2026: `+3.64 ns`）**：
   - 2025–2026 的 4 路 Preshower 跨接在两块数字化仪上：顶部 `Preshower0/1`（`ch 12, 13`）在 Board 0，底部低增益 `Preshower2/3`（`ch 24, 25`）在 Board 1。
   - 经 `+16.0 ns`（`status & 512`）校正后，`Preshower2/3` 仍比 `Preshower0/1` 晚 **`+4.04 ns`**（2025 核心均值 `+3.97 ± 1.75 ns`）或 **`+3.53 ns`**（2026 核心均值 `+3.64 ± 1.76 ns`，因 2026 runs `22408–23172` 中 `Preshower3` 有效率由 `29.5%` 降至 `1.2%`，以 `ch 24` 为主）。
   - 若不对齐 Dig2 相对 Dig1 的偏移，4 路平均时间会在 `13.3 ns`（仅 Dig1 命中）、`15.0 ns`（4 路同中）与 `17.3 ns`（仅 Dig2 命中，如 2025 `p0013` 中 `PS0/1` 发生 `84%–92%` `GFIT_FAILED` 时）之间系统漂移并展宽（`core_std` 从 `1.57 ns` 改善至 `1.23 / 1.09 ns`，模式峰值重新严格锚定 `13.21–13.31 ns`，见图 `fig04`）。

### 2.2 去硬编码重构（严格沿用 Calypso / `El9CalypsoForWaveForm` 风格）

修改集中在 `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/PhysicsAnalysis/NtupleDumper/` 的 4 个文件（提交 `7729ddf0` + 两级回退升级 `1c0f3f45`，相对 `geoRotation_latestMaster` 仅 `+267/-68` 行）：

- **[NtupleDumperAlg.h](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/PhysicsAnalysis/NtupleDumper/src/NtupleDumperAlg.h)**（基于 `official/master` 并保留 `m_rotatedPhysicalAngles`）：
  - **完全删除 `resolveReducedChargeTiming()`** 及任何 `-999` 哨兵值。
  - 声明清晰的 Gaudi `DoubleProperty`：`VetoVetoNuTimeDifference`（默认 `33.165`）、`PreshowerVetoNuTimeDifference`（默认 `53.625`）、`PreshowerDig2ExtraOffset`（默认 `0.0`）、`ReducedChargeIntegrationTime`（默认 `35.0`）。
  - 保留与 `MyWaveformTimeRangeAlg` 对应的诊断分支：`VetoNu{0,1}_reduced_charge_valid`、`VetoNu_reduced_charge_fallback`（`0` = 无有效定时/回退 `raw_charge`，`1` = `Veto2`，`2` = `Preshower Dig1 (ch < 16)`，`3` = 对齐后的 `Preshower Dig2 (ch >= 16)`）、`VetoNu_expected_localtime`。
- **[NtupleDumperAlg.cxx](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/PhysicsAnalysis/NtupleDumper/src/NtupleDumperAlg.cxx)**（**零 run 编号硬编码**）：
  - `get_mean_scintillator_time()` 严格仿照 `MyWaveformTimeRangeAlg::get_local_time()`：使用命名位 `kTimingShiftedBit = 1U << xAOD::WaveformStatus::TIMING_SHIFTED`，以 `(status & ~kTimingShiftedBit) != 0U` 筛选好击中；若 `status & kTimingShiftedBit` 则 `t += 16.0`；若配置了 `dig2ExtraOffset != 0.0` 且 `hit->channel() >= kDig2ChannelStart (16)`（CAEN Board 1），则 `t -= dig2ExtraOffset`。
  - **两级 Preshower 回退逻辑（`1c0f3f45`）**：当 `vetoLocalTime` 为 `NaN` 时，按板卡边界 `kDig2ChannelStart = 16` 将 `m_wave_map["preshower"]` 分为 `Dig1 (ch < 16)` 与 `Dig2 (ch >= 16)`；**优先仅使用同板、定时更锐利的 `Dig1`（`fallback = 2`）**，避免正常事例混入 `Dig2` 较宽抖动或额外偏移误差；**仅当 `Dig1` 全部失效（如大簇射饱和或 2025 `p0013` `GFIT_FAILED`）时，才回退到施加 `+16 ns`（bit 512）与 `-PreshowerDig2ExtraOffset` 对齐后的 `Dig2`（`fallback = 3`）**。同时正常 `Veto2` 有效事例不再多余遍历 `preshowerContainer`。
  - 直接从 `IWaveformCableMappingTool` 填充的 `m_wave_map["veto2"]`、`m_wave_map["preshower"]` 与 `m_wave_map["vetonu"]` 获取各时期正确通道（COOL 数据库 `CABP200.db` 自动给出 CaloNu `ch 22, 23` 与 2025–2026 `ch 12, 13, 24, 25`，无需在 C++ 硬编码）。
  - 保留全部程序缺陷修复：`VetoNu` 按 `hit->channel()` 映射至层下标 `0/1` 彻底消除 `size() >= 3` 时数组越界；空 `wave_vector()` 保持 `m_wave_raw_charge[ch]` 回退不归零；容器 `isValid()` 保护；blinding 改用 `m_wave_map` 动态读取 veto 通道；在 `initialize()` 中打印一次当前生效的 Python 配置参数。
- **[NtupleDumperConfig.py](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/PhysicsAnalysis/NtupleDumper/python/NtupleDumperConfig.py)** 与 **[faser_ntuple_maker.py](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/PhysicsAnalysis/NtupleDumper/scripts/faser_ntuple_maker.py)**：
  - 将全部时期/MC 标定值集中管理在 Python 字典 `REDUCED_CHARGE_TIMING_PRESETS`（`data_default`、`data_caloNu`、`data_2025`、`data_2026`、`mc`）和辅助函数 `getReducedChargeTimingKwargs(preset, **overrides)` 中，后续微调无需修改或重编任何 C++。
  - `faser_ntuple_maker.py` 暴露 `--timingPreset`（默认 `--isMC` 取 `mc`，否则 `data_default`）及 `--vetoVetoNuTimeDiff`、`--preshowerVetoNuTimeDiff`、`--preshowerDig2ExtraOffset`、`--reducedChargeWindow` 四个浮点覆写参数。

---

## 3. 详细数据

### 3.1 2022–2026 真实数据与 MC24 站点时间差测量（4,040 万 `GoodTracks` 穿行 Muon，Fig. 13 协议）

使用多核 ROOT `RDataFrame` 脚本 [waveform_timing_study.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/waveform_timing_study.py)（22 并发进程，原生 `[-60, 60]` ns 分箱、重跑耗时 `174.0 s`），在 [datasets.json](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/faser_filelists/datasets.json) 对应的官方 2022–2026 PHYS 文件及 MC24 上严格执行 [FASER_DQ_RDF.py](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/FASER_DQ_RDF.py) 选择：
- 年度 GRL `GoodTimes && !ExcludedTimes`（2022–2024 v8，2025–2026 v10）
- `(Timing0..3_status & 4) == 0`、`distanceToCollidingBCID == 0`、`(TAP & 4) != 0`
- `GoodTracks`：`Track_nLayers >= 7 && (Track_Chi2/Track_nDoF) < 25 && (Track_nDoF+5) >= 12 && Track_pz0 > 20000 && RemoveDuplicates(Track_p0)`，`NGoodTracks >= 1`
- 闪烁体击中：`VetoNu/Veto2 raw_charge > 40 pC`、`Preshower0/1 raw_charge > 2.5 pC`、`(status & ~512) == 0`（统计文件：[timing_summary.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/timing_summary.json)）。

| 时期 / 样本 | 扫描 Run 数 / 文件数 | `GoodTracks` 事例数 | `Veto2 - VetoNu` (+16 ns) 峰值 / 核心均值±σ [ns] | `PS Dig1 (12,13) - VetoNu` 峰值 / 核心均值±σ [ns] | `PS Dig2 (24,25) - VetoNu` (+16 ns) [ns] | `PS Dig2 - Dig1` (+16 ns) [ns] | `Timing - VetoNu` [ns] | `Veto20 st==512` 占比 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **2022 Data** | 11 / 17 | 540,303 | `33.08` (`33.02 ± 0.55`) | `53.49` (`53.46 ± 0.57`) | — | — | `34.70` | `0.0%` |
| **2023 Data** | 11 / 20 | 1,819,191 | `33.09` (`33.02 ± 0.56`) | `53.51` (`53.48 ± 0.57`) | — | — | `34.70` | `0.0%` |
| **2024 pre-CaloNu** (`<15821`) | 13 / 24 | 7,126,068 | `33.18` (`33.11 ± 0.56`) | `53.61` (`53.58 ± 0.54`) | — | — | `34.86` | `0.0%` |
| **2024 early-CaloNu** (`15821–15844`) | 4 / 8 | 1,978,799 | `33.15` (`33.08 ± 0.56`) | `53.59` (`53.55 ± 0.54`) | — | — | `34.84` | `0.0%` |
| **2024 during-CaloNu** (`15845–16924`) | 17 / 29 | 7,199,520 | **`35.87`** (`35.80 ± 0.56`) | `53.59` (`53.55 ± 0.54`) | — | — | `34.84` | **`22.3%`** |
| **2024 post-CaloNu** (`>=16925`) | 12 / 23 | 6,578,042 | `33.18` (`33.10 ± 0.56`) | `53.62` (`53.58 ± 0.54`) | — | — | `34.87` | `0.0%` |
| **2025 Data** (`p0013`+`p0014`) | 21 / 37 | 10,695,828 | `33.21` (`33.13 ± 0.56`) | **`13.27`** (`13.36 ± 1.15`) | **`17.29`** (`17.29 ± 1.60`) | **`+4.04`** (`+3.97 ± 1.75`) | `34.88` | `0.0%` |
| **2026 Data** (`p0015`) | 17 / 31 | 4,462,527 | `33.23` (`33.16 ± 0.56`) | **`13.27`** (`13.38 ± 1.15`) | **`16.79`** (`16.98 ± 1.61`) | **`+3.53`** (`+3.64 ± 1.76`) | `34.99` | `0.0%` |
| **MC24 Muon FaserNu** (`100124`) | 1 / 4 | 301,366 | **`27.28`** (`27.25 ± 0.10`) | **`51.36`** (`51.35 ± 0.26`) | — | — | `32.47` | `0.0%` |
| **MC24 Muon CaloNu** (`100050`) | 1 / 4 | 180,981 | **`27.27`** (`27.23 ± 0.14`) | **`51.36`** (`51.36 ± 0.26`) | — | — | `32.46` | `0.0%` |
| **MC24 GENIE 22/23** (`200082/83/89`) | 3 / 18 | 16,419 | `22.67` (反向背散射迟到) | `45.55` (反向背散射迟到) | — | — | `28.89` | `0.0%` |
| **MC24 GENIE 24 FaserNu** (`200093/94/100`) | 3 / 18 | 16,275 | `22.13` (反向背散射迟到) | `45.69` (反向背散射迟到) | — | — | `26.50` | `0.0%` |
| **MC24 GENIE 24 CaloNu** (`200139/40/46`) | 3 / 18 | 11,292 | `25.49` (近靶背散射迟到) | `46.71` (近靶背散射迟到) | — | — | `27.74` | `0.0%` |

> **逐 run 巡检定位的两项 2025/2026 生产与硬件现象（见图 `fig06`）**：
> 1. **2025 `p0013` vs `p0014` 差异**：所有 `p0013` 生产 run（`19415, 19524, 19866, 19923, 19970, 20081`）均出现 `Preshower0 GFIT_FAILED` 达 **`84.4%–92.0%`** 且残余 `PS Dig1 - VN` 跳至 `30–33 ns`；而在同一 run 区间内经 **`p0014`** 重建的 `19336, 19878, 20012` 及全部后续 `p0014` run（`20170–20610`），`PS0 GFIT_FAILED` 降至 **`0.4%–0.7%`** 且 `PS Dig1 - VN` 稳定回归 **`13.24–13.31 ns`**，`PS Dig2` 在全部 `p0013/p0014` run 中始终稳定于 `17.1–17.4 ns`。
> 2. **2026 `Preshower3` 有效率变化**：runs `22132–22404` 及 `23176–23241` 中 `PS3_ok ≈ 29.4%`、`PS Dig2 - Dig1 = +4.00..+4.42 ns`（与 2025 `p0014` 完全一致）；runs `22408–23172` 中 `PS3_ok` 降至 `1.2%`（仅 `ch 24` 活跃），`PS Dig2 - Dig1` 同步平移 `~0.6 ns` 至 `+3.20..+3.50 ns`，均值对应 Python 预设 `data_2026: PreshowerDig2ExtraOffset = 3.64 ns`。

#### 施加 Python 预设后的全时期对齐闭合残差（对应 `fig13_10` 与 `fig13_11`）

将最终 Python 预设（`data_default: 33.165 / 53.625 / 0.0`、`data_caloNu: 35.783 / 53.512 / 0.0`、`data_2025: 33.165 / 13.27 / 3.97`、`data_2026: 33.165 / 13.27 / 3.64`、`mc: 27.21 / 51.41 / 0.0`）与 `status & 512` 的 `+16.0 ns` 校正逐事件作用于上述 4,040 万真实 muon 与 48.2 万 MC24 muon 后，全部真实数据时期（2022–2026）及两套 MC24 几何的 `Veto2` 与 `Preshower`（含 4-PMT Dig1/Dig2 对齐）残差均严格收敛至 `0.0 ns`（最大偏差 `≤ 0.15 ns`）；GENIE 中微子 MC 因 VetoNu 仅由反向背散射迟到粒子触发，残差按物理预期呈 `-3.2..-6.5 ns` 负偏移：

| 时期 / 样本 | 生效 Python Preset | `(t_Veto2 - t_VetoNu) - Δt_preset` 峰值 / 核心均值±σ [ns] | `(t_Preshower - t_VetoNu) - Δt_preset` 峰值 / 核心均值±σ [ns] |
| --- | --- | ---: | ---: |
| **2022 Data** | `data_default` | `-0.084` (`-0.150 ± 0.554`) | `-0.140` (`-0.163 ± 0.567`) |
| **2023 Data** | `data_default` | `-0.081` (`-0.146 ± 0.556`) | `-0.114` (`-0.145 ± 0.566`) |
| **2024 pre-CaloNu** | `data_default` | `+0.015` (`-0.059 ± 0.559`) | `-0.013` (`-0.048 ± 0.543`) |
| **2024 early-CaloNu** | `data_default` | `-0.002` (`-0.084 ± 0.559`) | `-0.040` (`-0.076 ± 0.541`) |
| **2024 during-CaloNu** | `data_caloNu` | `+0.087` (`+0.015 ± 0.559`) | `+0.081` (`+0.042 ± 0.542`) |
| **2024 post-CaloNu** | `data_default` | `+0.015` (`-0.061 ± 0.558`) | `-0.006` (`-0.043 ± 0.542`) |
| **2025 Data** | `data_2025` | `+0.040` (`-0.030 ± 0.555`) | `+0.059` (`+0.109 ± 1.227`) |
| **2026 Data** | `data_2026` | `+0.063` (`-0.004 ± 0.556`) | `+0.042` (`+0.146 ± 1.092`) |
| **MC24 Muon FaserNu** | `mc` | `+0.044` (`+0.017 ± 0.115`) | `-0.052` (`-0.055 ± 0.263`) |
| **MC24 Muon CaloNu** | `mc` | `+0.036` (`+0.010 ± 0.140`) | `-0.049` (`-0.052 ± 0.263`) |

#### 严格仿 Fig. 13 风格的全年份时间差分布图集（共 11 图，完整汇总见 [fig13_all_timing_diffs_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_all_timing_diffs_report.pdf)，脚本 [plot_fig13_all_timing_diffs.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/plot_fig13_all_timing_diffs.py)）

全套图严格复现 [Electronic_2D_measurement_internal_note-6.pdf](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/docs/Electronic_2D_measurement_internal_note-6.pdf) 第 22 页 Fig. 13 的单图布局、官方 FASER logo、`FASER Internal` 标识、四边内向主/次刻度、浅灰实线网格、`[-60, 60]` ns `1.0 ns` 分箱面积归一化（`Normalized`）及年份配色命名（并附 `fig13_03` 局部放大与 6 页多面板诊断汇总 [timing_study_fig13_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/figures/timing_study_fig13_report.pdf)）：

![Fig. 13-01: All years Veto St2 - VetoNu (+16 ns corrected)](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_01_veto2_vetonu_all_years.png)

![Fig. 13-02: All years Veto St2 - VetoNu raw (before +16 ns fix)](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_02_veto2_vetonu_raw_uncorrected.png)

![Fig. 13-04: All years Preshower Dig1 - VetoNu](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_04_preshower_dig1_vetonu_all_years.png)

![Fig. 13-05: 2025-2026 Preshower Dig1 vs Dig2 - VetoNu](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_05_preshower_dig2_vs_dig1_vetonu_2025_2026.png)

![Fig. 13-06: 2025-2026 Preshower Dig2 - Dig1 direct offset](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_06_preshower_dig2_minus_dig1_2025_2026.png)

![Fig. 13-07: 2025-2026 Preshower 4-PMT alignment effect](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_07_preshower_4pmt_dig2_alignment_effect.png)

![Fig. 13-08: All years Timing - VetoNu control](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_08_timing_vetonu_all_years.png)

![Fig. 13-09: All years Veto St1 - VetoNu (+16 ns corrected)](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_09_veto1_vetonu_all_years.png)

![Fig. 13-10: All years Veto St2 post-preset alignment closure](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_10_veto2_after_preset_alignment.png)

![Fig. 13-11: All years Preshower post-preset alignment closure](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_11_preshower_after_preset_alignment.png)

> 注：内部笔记 Fig. 13 图注文字称 CaloNu 期间“FASERν veto 更换了数字化仪”，但由 `fig13_08`（`t_Timing - t_νveto` 全数据时期平直）与 `fig13_04`（`t_Preshower(Dig1) - t_νveto` 在 2024 CaloNu 期间平直）对比 `fig13_01`、`fig13_09`（`Veto2` 与 `Veto1` 同步跳变 `+2.7..+2.8 ns` 且出现 `22.3%` 的 `status == 512` 事例）可直接证明：**实际迁移至 Digitizer 2（`ch 20–23`）的是 Veto St1/St2，而 VetoNu 始终留在 Digitizer 1（`ch 4, 5`）**。局部放大图见 [fig13_03_veto2_vetonu_zoom.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_study_fig13/fig13_exact/fig13_03_veto2_vetonu_zoom.png)。

---

### 3.2 端到端编译验证：2022–2026 数据与 Ben Wilson 5 套参考目录对标

使用 [validate_reduced_charge_all_periods.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/validate_reduced_charge_all_periods.py)（16 并发运行 41 个真实 `faser_ntuple_maker.py` 作业，全部 `rc = 0`，耗时 `663.0 s`；汇总：[validation_summary_v2.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v2/validation_summary_v2.json)）：

| 真实数据样本 | Python `--timingPreset` | 事例数 | vs Ben 参考目录 `\|Δq0\| < 10⁻³ pC` / `fallback` 一致率 | Preshower fallback muon (`N`)：基线 vs 重构中位数 [pC] | Preshower fallback `< 40 pC` 误判率（基线 -> 重构） | `Veto2 st==512` muon (`N`)：基线 vs 重构中位数 [pC] |
| --- | --- | ---: | --- | --- | --- | --- |
| `2022_008090` (`r0021`) | `data_default` | 5,000 | **`100.00%`** (`max 7.6e-5`) / **`100.00%`** | `N=0` | — | `N=0` |
| `2023_011703` (`r0021`) | `data_default` | 5,000 | **`100.00%`** (`max 1.2e-4`) / N/A | `N=9`: `78.01` vs `78.01` | `0.00%` -> `0.00%` | `N=0` |
| `2024_pre_015387` (`r0022`) | `data_default` | 5,000 | **`100.00%`** (`max 3.7e-4`) / N/A | `N=128`: `97.70` vs `97.70` | `0.00%` -> `0.00%` | `N=0` |
| `2024_early_015838` (`r0022`) | `data_default` (覆写 `35.783/53.512` 验证 Ben) | 5,000 | 覆写 **`100.00%`** (`max 3.8e-4`)；`data_default` 匹配真实 Dig1 `33.15 ns` | `N=146`: `91.57` vs `91.57` | `0.00%` -> `0.00%` | `N=0` |
| `2024_CaloNu_015849` (`r0022`) | `data_caloNu` | 8,000 | **`99.98%`** (`100% < 1.1e-3`) / N/A | `N=247`: `91.49` vs `91.49` | `0.00%` -> `0.00%` | **`N=620`**: `89.28` -> **`83.86`**（恢复 `512` +16 ns） |
| `2024_post_016929` (`r0022`) | `data_default` | 5,000 | **`100.00%`** (`max 4.9e-4`) / **`100.00%`** | `N=101`: `86.88` vs `86.88` | `0.00%` -> `0.00%` | `N=0` |
| `2025_20012` (`r0023`) | `data_2025` (`13.27 / 3.97`) | 8,000 | 无 Ben 参考 | **`N=67`**: `11.30` -> **`87.78`** | **`58.21%`** -> **`1.49%`** | `N=0` |
| `2025_20481` (`p0014`) | `data_2025` (`13.27 / 3.97`) | 8,000 | 无 Ben 参考 | **`N=67`**: `0.18` -> **`87.57`** | **`95.52%`** -> **`0.00%`** | `N=0` |
| `2026_22404` (`p0015`) | `data_2026` (`13.27 / 3.64`) | 8,000 | 无 Ben 参考 | **`N=135`**: `0.33` -> **`92.26`** | **`97.04%`** -> **`0.00%`** | `N=0` |
| `2026_23176` (`p0015`) | `data_2026` (`13.27 / 3.64`) | 8,000 | 无 Ben 参考 | **`N=97`**: `0.23` -> **`87.40`** | **`91.75%`** -> **`0.00%`** | `N=0` |
| `MC24_Muon_100044` (`s0013`) | `mc` (`27.21 / 51.41 / 0.0`) | 3,000 | — | `N=165`: `89.52` vs `89.52` | `0.00%` -> `0.00%` | `N=0` |

#### 两级 Preshower 回退（`1c0f3f45`：优先 `Dig1 fb=2`，失效再回退对齐 `Dig2 fb=3`）复测结果

重新编译 `software_ntuple_fix_v1_clean` 后对 11 组代表样本（22 个真实 `faser_ntuple_maker.py` 作业，16 并发耗时 `330.9 s`，全部 `rc = 0`；汇总：[validation_summary_v2.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v3_twotier/validation_summary_v2.json)）执行端到端复测：
- **2022–2024 数据与 MC24（Muon + 3 几何代表 GENIE `200082/200093/200139`）**：`fb == 3` 恒为 `0`，Ben Wilson 对标（`2022: 100.00%`, `2024_post: 100.00%`, `2024_CaloNu: 99.98%`）与 GENIE 挽救率（`35.94% / 40.74% / 23.33%`）与之前逐事例完全一致。
- **2025–2026 4-PMT Preshower 实测 `fb=2`（Dig1）vs `fb=3`（Dig2）分工**：
  - `2025_20481` (`p0014`)：`fb = {0: 3529, 1: 3956, 2: 509, 3: 6}`；穿行 muon 中 `66/67`（`98.5%`）走纯 Dig1（`fb=2`，中位数 `87.68 pC`，`<40 pC = 0.00%`），`1` 个 Dig1 失效事例由对齐 Dig2 兜底（`fb=3`，`92.68 pC`）。
  - `2026_22404` (`p0015`)：`fb = {0: 2178, 1: 5130, 2: 689, 3: 3}`；穿行 muon 中 `134/135`（`99.3%`）走纯 Dig1（`fb=2`，中位数 `91.79 pC`，`<40 pC = 0.00%`），`1` 个 Dig1 失效事例由对齐 Dig2 兜底（`fb=3`，`104.01 pC`）。
  - `2026_23176` (`p0015`)：`fb = {0: 2131, 1: 5244, 2: 620, 3: 5}`；穿行 muon 中 `92/97`（`94.8%`）走纯 Dig1（`fb=2`，中位数 `87.61 pC`，`<40 pC = 0.00%`），`5` 个 Dig1 失效事例全部由对齐 Dig2 兜底（`fb=3`，中位数 `84.42 pC`，`<40 pC = 0.00%`）。
  - `2025_20012` (`r0023`，`PS0/1` 高失效率样本)：`fb = {0: 2270, 1: 5314, 2: 69, 3: 347}`；`7` 个 Dig1 正常 muon 走 `fb=2`（`<40 pC = 0.00%`），`60` 个 Dig1 失效 muon 自动落入对齐 Dig2 `fb=3`（中位数 `85.65 pC`，`<40 pC` 由基线 `58.21%` 降至 `1.67%`）。

---

### 3.3 端到端编译验证：`FASER_FaserMC.pdf` 全部 9 个 Official Default-Flux GENIE 中微子 MC Run

在 `/eos/experiment/fasernu-data0/faser/sim/mc24/genie/<run>/rec/{s0014,s0015}/` 中对三类几何/交叉角下的全部 9 个 default flux prediction run 各处理 `3,000` 个完整 xAOD 事例（共 `27,000` 事例，Python `--isMC` 自动选择 `--timingPreset mc`：`Veto2-VN = 27.21 ns`, `PS-VN = 51.41 ns`, `Dig2Extra = 0.0 ns`）。下游中微子事例选择：`longTracks > 0 && Veto20/21_raw_charge > 40 pC && Preshower0/1_raw_charge > 2.5 pC`：

| MC24 GENIE 类别（`FASER_FaserMC.pdf`） | Run 编号与产生模式 | 选中事例 `N_sel` | `Veto2` 饱和占比 (`status & 4`) | Timing 来源 `fb1(Veto2)` / `fb2(PS)` / `fb0` | `raw >= 40 pC` 背散射事例数 | `reduced < 30 pC` 挽救事例数（挽救率） | `pass raw<40` -> `pass red<30` | `max \|q_fix - q_base\|` [pC] |
| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| **2022/23 Crossing angle** (`s0014`, p. 11) | `200082` (`epos_pi`) | 326 | `35.9%` | `230 / 95 / 1` (`fb2 = 29.1%`) | 64 | `23` (**`35.9%`**) | `262 -> 285` | `0.0001` |
| **2022/23 Crossing angle** (`s0014`, p. 11) | `200083` (`epos_k`) | 311 | `50.5%` | `167 / 142 / 2` (`fb2 = 45.7%`) | 78 | `25` (**`32.1%`**) | `233 -> 258` | `0.0001` |
| **2022/23 Crossing angle** (`s0014`, p. 11) | `200089` (`charm_central`) | 186 | `59.7%` | `85 / 101 / 0` (`fb2 = 54.3%`) | 53 | `25` (**`47.2%`**) | `133 -> 158` | `0.0001` |
| **2024 FaserNu geometry** (`s0014`, p. 13–14) | `200093` (`epos_pi`) | 334 | `38.3%` | `224 / 106 / 4` (`fb2 = 31.7%`) | 54 | `22` (**`40.7%`**) | `280 -> 302` | `0.0002` |
| **2024 FaserNu geometry** (`s0014`, p. 13–14) | `200094` (`epos_k`) | 277 | `52.7%` | `153 / 124 / 0` (`fb2 = 44.8%`) | 72 | `27` (**`37.5%`**) | `205 -> 232` | `0.0001` |
| **2024 FaserNu geometry** (`s0014`, p. 13–14) | `200100` (`charm_central`) | 187 | `54.6%` | `96 / 91 / 0` (`fb2 = 48.7%`) | 57 | `26` (**`45.6%`**) | `130 -> 156` | `0.0002` |
| **2024 CaloNu geometry** (`s0015`, p. 15–16) | `200139` (`epos_pi`) | 238 | `55.0%` | `124 / 114 / 0` (`fb2 = 47.9%`) | 60 | `14` (**`23.3%`**) | `178 -> 192` | `0.0002` |
| **2024 CaloNu geometry** (`s0015`, p. 15–16) | `200140` (`epos_k`) | 201 | `65.2%` | `83 / 116 / 2` (`fb2 = 57.7%`) | 58 | `9` (**`15.5%`**) | `143 -> 152` | `0.0002` |
| **2024 CaloNu geometry** (`s0015`, p. 15–16) | `200146` (`charm_central`) | 127 | `72.4%` | `42 / 84 / 1` (`fb2 = 66.1%`) | 37 | `4` (**`10.8%`**) | `90 -> 94` | `0.0002` |

#### 端到端验证主图（完整 3 图见 [validation_v2_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v2/figures/validation_v2_report.pdf)）

![End-to-End 2024 CaloNu and 2025-2026 Preshower Fallback Recovery](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v2/figures/figV1_data_caloNu_and_2025_2026_recovery.png)

![Official Default-Flux MC24 GENIE 9-Run Validation](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v2/figures/figV2_genie_default_flux_9runs_validation.png)

![Per-Run GENIE 9-Run Fallback and Backsplash Summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v2/figures/figV3_genie_9runs_fallback_and_oob_summary.png)

### 3.4 执行命令、回归测试与产物清单

- **构建命令**：
  ```bash
  source /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/cached_athena_env.sh
  python3 /eos/home-x/xcheng/FASER/geoRotation_test/scripts/build_fix_v1_ntuple.py \
    --output /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/software_ntuple_fix_v1_clean \
    --jobs 4
  ```
- **Fig. 13 真实数据时间差扫描与绘图**：
  ```bash
  source /afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/setup.sh
  python3 /eos/home-x/xcheng/FASER/geoRotation_test/analysis/waveform_timing_study.py --workers 22 --files-per-run 2 --runs-per-period 12
  python3 /eos/home-x/xcheng/FASER/geoRotation_test/analysis/plot_fig13_all_timing_diffs.py
  python3 /eos/home-x/xcheng/FASER/geoRotation_test/analysis/plot_waveform_timing_fig13.py
  ```
- **端到端验证与绘图**（自动清理临时 `.root` 输出）：
  ```bash
  python3 /eos/home-x/xcheng/FASER/geoRotation_test/scripts/validate_reduced_charge_all_periods.py \
    --overlay /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/software_ntuple_fix_v1_clean \
    --env-script /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/cached_athena_env.sh \
    --outdir /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v2 \
    --workers 16
  python3 /eos/home-x/xcheng/FASER/geoRotation_test/analysis/plot_validation_v2.py
  ```
- **回归单元测试**：`python3 -m unittest discover -s tests -v`（全部 `121/121` 通过，其中 [test_reduced_charge_logic.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/tests/test_reduced_charge_logic.py) 显式断言 C++ 中零硬编码 run 编号并校验全部 Python 预设与边界条件）。

---

## 4. 阶段性结论

### 已由实测数据证实的事实
1. **零硬编码 C++ 架构**：`NtupleDumperAlg.{h,cxx}`（`7729ddf0`）不包含任何 run 编号分支，完全依靠 `xAOD::WaveformStatus::TIMING_SHIFTED (512)`、`SECONDARY (2)`、`IWaveformCableMappingTool`（`m_wave_map`）与标准 Gaudi `DoubleProperty` 工作；全部时期和 MC 时间差配置清晰集中于 Python `NtupleDumperConfig.REDUCED_CHARGE_TIMING_PRESETS` 与 `faser_ntuple_maker.py` CLI 参数，代码简洁易维护。
2. **4,040 万真实 `GoodTracks` 穿行 muon 的 Fig. 13 实测完全解释了三项参数的必要性**：
   - `VetoVetoNuTimeDifference`：2022–2024 pre/early/post 及 2025–2026（Dig1 `ch 6,7`）稳定在 `33.08–33.23 ns`；2024 CaloNu（Dig2 `ch 22,23`）除 `22.3%` 事例出现 `-16.0 ns`（`status == 512`）跳变外，校正后仍有 **`+2.69 ns`** 跨板固定偏移（`35.87 ns`，对应 `data_caloNu: 35.783 ns`）；MC24 因 `WaveformDigitisationTool` 已扣 TOF 并按 CB 模板中心（`825 vs 800 ns`）数字化而位于 `27.27–27.28 ns`（对应 `mc: 27.21 ns`）。
   - `PreshowerVetoNuTimeDifference`：2022–2024 稳定于 `53.49–53.62 ns`，MC24 为 `51.36 ns`；2025–2026 升级新 PMT/电缆后整体前移 **`-40.3 ns`** 至 **`13.27 ns`**（Dig1 `ch 12,13`）。
   - `PreshowerDig2ExtraOffset`：2025–2026 新增底部 Dig2 `ch 24,25` 在 `+16.0 ns` 校正后仍晚于 Dig1 `ch 12,13` 约 **`+3.97 ns`**（2025）与 **`+3.64 ns`**（2026，受 `22408–23172` 中 `PS3_ok` 降至 `1.2%` 影响）；扣除该偏移后 4 路 Preshower 均值峰值重新稳定收敛至 `13.21–13.31 ns`。
3. **2022–2026 真实数据端到端验收**：与 Ben Wilson 5 套参考目录逐事件一致率达 **`100.00%`**（`15849` 为 `99.98% < 10⁻³ pC`, `100.00% < 1.1×10⁻³ pC`）；2025–2026 的 366 个 Preshower-fallback 穿行 muon 的 `< 40 pC` 误判率从基线的 **`88.25%`** 降至 **`0.27%`**（中位数从 `0.49 pC` 恢复至 **`87.4–92.3 pC`**）。
4. **9 个官方默认通量 GENIE MC 端到端验收**：`2022/23 Crossing angle`（`200082/83/89`）、`2024 FaserNu`（`200093/94/100`）与 `2024 CaloNu`（`200139/40/46`）共 `27,000` 事例全部零错误运行，与基线在浮点精度内严格闭合（`max |Δq| ≤ 0.0002 pC`）；同时实测揭示 `35.9%–72.4%` 的选中中微子事例因前向簇射使 `Veto2` 饱和而必须依赖 Preshower fallback（`29.1%–66.1%`），并证实 `reduced_charge < 30 pC` 可在 FaserNu 几何挽救 **`32%–47%`**、在 CaloNu 几何挽救 **`11%–23%`** 的 `raw >= 40 pC` 背散射事例。

### 适用边界说明
- `2025 p0013` 早期生产文件中 `Preshower0/1` 存在 `84%–92%` 的 `GFIT_FAILED`（重建配置缺陷，已在 `p0014` 中完全消除）；若在旧 `p0013` xAOD 上运行，因 `PS2/3`（Dig2 `ch 24,25`）保持正常且通过 `PreshowerDig2ExtraOffset = 3.97 ns` 对齐，fallback 仍保持无偏。
- 历史 `MC22`（`WAVE-Digitization-01`，`cb_mean` 均为 `815 ns`）不适用 `MC24` 预设，可通过 `--vetoVetoNuTimeDiff 2.34 --preshowerVetoNuTimeDiff 32.44` 从 Python 显式传入。

---

## 5. 下一步计划

1. 在批量生产（HTCondor）调用 `faser_ntuple_maker.py` 时，按数据时期在提交脚本中显式传入 `--timingPreset`（`data_default` / `data_caloNu` / `data_2025` / `data_2026`，MC 自动使用 `mc`），与 `El9CalypsoForWaveForm/submit_waveform_time_range_cfg.py` 的提交方式保持完全一致。
2. 若将 `fix/latest-seed-and-momentum-v1` 的跟踪修复（`e47b85fd`）与本次 `NtupleDumper` 改进（`7729ddf0`）合并进入主干构建，可在完整 worktree 中执行一次 `build_local.sh` 与 `verify_local.sh` 集成验收。
