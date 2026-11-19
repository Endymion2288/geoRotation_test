# 阶段报告：在 `geoRotation_calypso_latest_fix_v1` 中合并 2022–2026 数据及 MC 的 Waveform Reduced VetoNu Charge 处理程序

- **日期**：2026-11-19
- **主机与环境**：`lxplus901.cern.ch`，Athena `24.0.41`（`x86_64-el9-gcc13-opt`）
- **基线源码（未改动）**：`/eos/home-x/xcheng/FASER/geoRotation_calypso_latest`，HEAD `7fd40df8cd3b60e81b4095a4f84310c4a54700a0`
- **工作树源码与提交**：`/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1`（分支 `fix/latest-seed-and-momentum-v1`，提交 `a4ecefeaa151040f57aa20be041133bb3c1f0e16`）
- **可复用补丁**：[0051_ntupledumper_waveform_reduced_charge_2022_2026_mc.patch](file:///eos/home-x/xcheng/FASER/geoRotation_test/patches/0051_ntupledumper_waveform_reduced_charge_2022_2026_mc.patch)
- **编译产物（Overlay）**：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/software_ntuple_fix_v1/`
- **完整图集与汇报 PDF**：[waveform_reduced_charge_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/figures/waveform_reduced_charge_report.pdf)（索引：[figures_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/figures/figures_index.json)）

---

## 1. 阶段目标

参考 `/eos/home-x/xcheng/FASER/El9CalypsoForWaveForm`（`MyWaveformTimeRangeAlg.cxx`、`MyWaveformTimeRangeConfig.py`、`submit_waveform_time_range_cfg.py`）、通道映射文档 `/eos/home-x/xcheng/FASER/nu_muAnalysis/docs/FASER_CaloScintillator.pdf` 与 COOL 数据库 `CABP200.db`/`ALLP200.db`，在 `geoRotation_calypso_latest_fix_v1` 的 `PhysicsAnalysis/NtupleDumper` 中完成正确、通用的 waveform 与 reduced VetoNu charge 处理程序移植与升级，实现以下验收目标：

1. `NtupleDumperAlg.h` 以 `official/master` 版本为基准（保留 `geoRotation_calypso_latest` 已有且 `NtupleDumperAlg.cxx` 依赖的 `m_rotatedPhysicalAngles` 开关）。
2. 正确处理第 2 数字化仪（Digitizer 2）的 `TIMING_SHIFTED` 状态位（bit 9 = `512`）及 `+16.0 ns` localtime 校正，并覆盖 2022–2026 全部 Veto St2（`ch 6, 7, 22, 23`）与 Preshower（`ch 12, 13, 24, 25`）通道。
3. 在 2022、2023、2024（pre-CaloNu / 早期 CaloNu `15821–15844` / CaloNu `15845–16927` / post-CaloNu）、2025、2026 真实数据及 MC 上直接测量穿行 muon 峰值时间差，实现按 `(isMC, run_number, cable_map)` 自动解析且支持命令行/属性显式覆盖。
4. 修复原 `NtupleDumperAlg.cxx` 中除 `status==512` 之外的全部已知程序缺陷（`VetoNu` 容器索引越界与通道错位、空 `wave_vector()` 归零、CaloNu blinding 硬编码旧通道等）。
5. 使用真实编译后的 Calypso `faser_ntuple_maker.py` 在 2022–2026 全部关键时期及 MC 上完成端到端运行验证，并与 Ben Wilson 的 5 套参考目录（`/eos/user/b/bewilson/El9CalypsoForWaveForm/`）进行逐事件数值交叉核对。

---

## 2. 工作内容

### 2.1 原 `NtupleDumper` 缺陷审计（基线 `7fd40df8`）

1. **`get_mean_scintillator_time()` 丢弃 Digitizer 2 `TIMING_SHIFTED (512)` 击中且缺少 `+16.0 ns` 修正**（原 `NtupleDumperAlg.cxx:2638`）：原条件 `hit->hit_status() == 0` 拒绝全部 `status == 512` 击中。根据 `WaveformReconstructionTool.cxx:111-114,331-345,446-448`，Digitizer 2 在 `trigger_offset >= 1` 时将采集窗口前移 8 个采样点（`16 ns`）并置位 `xAOD::WaveformStatus::TIMING_SHIFTED`（bit 9 = `512`），`trigger_time`/`bcid_time` 已扣减 `-16 ns`，而 `localtime()` 保留原始窗口时间。因此必须以 `(status & ~512U) == 0U` 筛选好击中，并在 `status & 512U` 时执行 `t += 16.0`。
2. **`VetoNu` 循环用容器下标 `idx` 直接索引固定数组 `m_vetonu_reduced_charge[2]` 与 `good_vetoNu_hit[2]`**（原 `NtupleDumperAlg.cxx:1618-1642`）：xAOD 扫描证实 2022–2026 与 MC 中每 40,000 事例存在 `140–956` 个事例满足 `VetoNuWaveformHits->size() >= 3`（含 secondary hit）。原代码以 `idx = 0..size()-1` 越界访问 `bool[2]` 与 `float[2]`（紧邻成员 `m_clock_phase`），且在通道乱序时将 `VetoNu0`/`VetoNu1` 对调。
3. **空 `wave_vector()` 保护缺失**（原 `NtupleDumperAlg.cxx:1630`）：若某 VetoNu hit 未保存原始波形（`WaveformCleanupTool` prescale 或 `WAVEFORM_MISSING`），原代码先置零再跳过空循环，将缺失波形误判为 `0 pC` 中微子信号候选。
4. **单一硬编码时间差在 2024 CaloNu 与 2025–2026 严重失配**（原 `faser_ntuple_maker.py:384-386`）：对所有数据统一硬编码 `VetoVetoNuTimeDifference=33.165 ns`、`PreshowerVetoNuTimeDifference=53.625 ns`。
   - 2024 CaloNu（`15845–16927`）Veto St2 迁至 Dig2 `ch 22, 23` 后实测 `Veto2 - VetoNu = 35.79–35.87 ns`（偏置 `+2.62 ns`）。
   - 2025–2026（`run >= 18194`）Preshower 升级为上下 4 PMT（Dig1 `ch 12, 13` + Dig2 `ch 24, 25`）且硬件触发偏移由 `+18 ns` 调整至 `-40/-42 ns`，实测 `Preshower - VetoNu` 骤降约 `40 ns` 至 `13.28 ns`（Dig1）与 `17.25 ns`（Dig2）！沿用 `53.625 ns` 会使 Preshower fallback 积分截止时间提前约 `40 ns`，将 **88–97% 的穿行 muon 电荷截断为 `< 40 pC` 假中微子事例**。
5. **缺少诊断分支与容器有效性检查**：原代码未检查 `vetoContainer.isValid()`/`preshowerContainer.isValid()`，且未输出 `fallback_mode`、`reduced_charge_valid` 与 `expected_localtime`。
6. **Blinding 逻辑硬编码旧 Veto 通道**（原 `NtupleDumperAlg.cxx:1550`）：硬编码 `ch 4,5,6,7,14`，在 2024 CaloNu 期间 `ch 6,7,14` 实为 CaloNu 而非 Veto（`ch 20,21,22,23`）。

### 2.2 硬件通道与时间差全景核查（2022–2026 数据 + MC）

1. 直接查询 `/cvmfs/faser.cern.ch/repo/sw/database/DBRelease/current/sqlite200/CABP200.db`（2026-10-28）与 `ALLP200.db`（2026-11-10）的 `/WAVE/DAQ/CableMapping`、`/WAVE/DAQ/TimingOffset` 与 `/WAVE/Digitization`，完整冻结 2022–2026 数据（`CONDBR3`）与 MC（`OFLP200`）全部 IOV。
2. 编写多核 PHYS 扫描脚本 [measure_scint_timing_diffs.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/measure_scint_timing_diffs.py)（16 并发扫描 82 个 PHYS 文件）与原生 C++ (`g++ -O2`) xAOD 扫描器 [scan_xaod_waveform.cxx](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/scan_xaod_waveform.cxx) + [measure_xaod_timing_diffs.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/measure_xaod_timing_diffs.py)（12 并发直读 55 个 2022–2026 及 MC xAOD 文件），对所有时期测量穿行 muon 峰值时间差，并分别对比 `status == 0` 与 `status == 512`（加/不加 `+16.0 ns`）。

### 2.3 代码修改与编译验证

修改均位于 `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/PhysicsAnalysis/NtupleDumper/`（仅 3 个文件，`+300/-63` 行）：

- **`src/NtupleDumperAlg.h`**（基于 `official/master` + 保留 `m_rotatedPhysicalAngles`）：
  - 将 `VetoVetoNuTimeDifference`、`PreshowerVetoNuTimeDifference` 与新增 `PreshowerDig2ExtraOffset` 默认值设为哨兵 `-999.0`（`<= -900.0` 启用自动周期查找，`> -900.0` 立即采用用户显式覆写）。
  - 声明 `resolveReducedChargeTiming()` 与支持 `dig2ExtraOffset` 的 `get_mean_scintillator_time()`；新增树分支变量 `m_vetonu_reduced_charge_valid[2]`、`m_vetonu_reduced_charge_fallback`、`m_vetonu_expected_localtime`。
- **`src/NtupleDumperAlg.cxx`**：
  - 引入 `xAODFaserWaveform/WaveformDefs.h`，使用命名位 `xAOD::WaveformStatus::TIMING_SHIFTED`（`512`）与 `SECONDARY`（`2`）。
  - `get_mean_scintillator_time()`：接受 `(status & ~512U) == 0U`；对 `status & 512U` 加 `+16.0 ns`；对 `channel >= 16` 扣除 `dig2ExtraOffset`（用于将 2025–2026 Preshower Dig2 `ch 24, 25` 对齐至 Dig1 `ch 12, 13` 参考系），双精度求平均。
  - `resolveReducedChargeTiming()`：按 `(isMC, run_number, mappedVeto2Channels)` 自动给出 2022–2026 及 MC24 实测标定时间差，首事例输出一条 `ATH_MSG_INFO` 日志记录生效参数。
  - Reduced VetoNu charge 计算：扩展 Veto2 通道 `{6, 7, 22, 23}` 与 Preshower 通道 `{12, 13, 24, 25}`；通过 `vetoNuChannels`（`4 -> 0, 5 -> 1`）映射 VetoNu 层下标彻底消除越界与对调；保持 primary-first 规则且空 `wave_vector()` 回退至 `m_wave_raw_charge[ch]`；移除提前 `break` 并在 `double` 中积分 `< m_reducedChargeIntegrationTime`；填充 `VetoNu{0,1}_reduced_charge`、`VetoNu{0,1}_reduced_charge_valid`、`VetoNu_reduced_charge_fallback`（`0/1/2`）与 `VetoNu_expected_localtime`。
  - Blinding 改为按 `m_wave_map` 动态读取 `vetonu`/`veto1`/`veto2` 通道，修复 CaloNu 期间通道错位。
- **`scripts/faser_ntuple_maker.py`**：移除硬编码 `33.165 / 53.625` 与 `27.21 / 51.41`，新增 `--vetoVetoNuTimeDiff`、`--preshowerVetoNuTimeDiff`、`--preshowerDig2ExtraOffset`、`--reducedChargeWindow` 可选 CLI 覆写参数。
- 使用 [build_fix_v1_ntuple.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/build_fix_v1_ntuple.py) 复用 `geoRotation_calypso_latest/run` 依赖编译 `libNtupleDumper.so` 并运行 `listcomponents`、`genconf` 与 `mergeConfdb2.py`，随后用 [validate_reduced_charge_all_periods.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/validate_reduced_charge_all_periods.py) 并发执行 35 个真实 `faser_ntuple_maker.py` 作业（覆盖 12 个 2022–2026 与 MC 样本的固定自动模式、未改基线对比及 Ben 精确覆写对齐）。

---

## 3. 详细数据

### 3.1 COOL 数据库 IOV 与 2022–2026 / MC 穿行 Muon 时间差测量

所有时间差均定义为 $\Delta t = t_{\mathrm{station,corr}} - t_{\mathrm{VetoNu,corr}}$（单位：ns，局部时间 `localtime` 坐标系；`status & 512` 击中已加 `+16.0 ns`；筛选条件为单长径迹或 VetoNu 两层 `raw_charge > 40 pC` 且 `(status & ~512) == 0`）：

| 时期 / 样本 | Run 范围 (`CABP200.db` IOV) | Veto St2 通道 | Preshower 通道 | 实测 `Veto2 - VetoNu` [ns] | 自动查找 `Veto2 - VetoNu` [ns] | 实测 `PS Dig1 (12,13) - VetoNu` [ns] | 实测 `PS Dig2 (24,25) - VetoNu` [ns] | 自动查找 `PS(Dig1)` / `Dig2Extra` [ns] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2022–2023 Data | `6525..13846` | `6, 7` (Dig1) | `12, 13` (Dig1) | `33.03–33.09` (`σ≈0.55`) | `33.055` | `53.46–53.51` | N/A | `53.490` / `0.000` |
| 2024 pre-CaloNu | `13847..15820` | `6, 7` (Dig1) | `12, 13` (Dig1) | `33.06–33.21` | `33.165` | `53.47–53.61` | N/A | `53.550` / `0.000` |
| 2024 早期 CaloNu（Veto 仍在 Dig1） | `15821..15844` | `6, 7` (Dig1) | `12, 13` (Dig1) | `33.05–33.07`（4 runs: `15838–15844`） | `33.165` | `53.47–53.49` | N/A | `53.550` / `0.000` |
| 2024 CaloNu（Veto 移至 Dig2） | `15845..16927` | **`22, 23` (Dig2)** | `12, 13` (Dig1) | `35.78–35.87` (`st0` 与 `st512+16` 差 `<0.02`) | `35.783` | `53.49–53.57` | N/A | `53.512` / `0.000` |
| 2024 post-CaloNu | `16928..18192` | `6, 7` (Dig1) | `12, 13` (Dig1) | `33.07–33.18` | `33.165` | `53.51–53.60` | N/A | `53.550` / `0.000` |
| 2025 启动调试（单 run） | `18193` | `22, 23` (Dig2) | `12, 13, 24, 25` | 未测（无公开束流文件；按 `TimingOffset` `-28 vs -40` 推算） | `35.783` | 推算 `25.28` | 推算 `29.25` | `25.280` / `3.970` |
| 2025 Data | `18194..20701` | `6, 7` (Dig1) | **`12, 13` + `24, 25` (Dig2)** | `33.07–33.21` | `33.165` | `13.22–13.38`（注：`19480–20080` `PS0,1` 多为 `GFIT_FAILED`） | `17.09–17.39` (`st0` 与 `st512+16` 差 `<0.07`) | `13.280` / `3.970` |
| 2026 Data（质子物理段） | `22126..23176` | `6, 7` (Dig1) | **`12, 13` + `24, 25` (Dig2)** | `33.12–33.19` | `33.165` | `13.23–13.36` | `16.47–17.41`（`ch24≈16.6, ch25≈18.6`） | `13.280` / `3.970` |
| 2026 末期（`>= 23359`） | `23359..24068` | `6, 7` (Dig1) | `12, 13, 24, 25` | 实测 `VetoNu ch 4,5` 无信号（`q_p95 < 0.8 pC`，HV 关或拆除） | `33.165` | N/A | N/A | `13.280` / `3.970` |
| MC24（`WAVE-Digitization-02+`） | `isMC = true` (`100012/43/44/49/50`, `210010`) | `6, 7` | `12, 13` | `27.115`（标准）/ `27.314`（CaloNu `s0015`） | `27.210` | `51.33–51.39` | N/A | `51.410` / `0.000` |
| MC22（历史 `WAVE-Digitization-01`） | `210001` (`r0013`) | `6, 7` | `12, 13` | `2.335`（`cb_mean` 均为 `815 ns`） | 需 CLI 显式覆写 | `32.442` | N/A | 需 CLI 显式覆写 |

> **Digitizer 2 `TIMING_SHIFTED (512)` `+16.0 ns` 实测闭合验证**：
> - 2024 CaloNu Run `15849`（Veto2 `ch 22`）：`status == 0` 中位数 **`35.168 ns`**（`N=2057`）；`status == 512` 修正前 **`19.160 ns`**（`N=643`），加 `+16.0 ns` 后 **`35.160 ns`**（差仅 **`0.008 ns`**）。
> - 2026 Run `22404`（Preshower `ch 24`）：`status == 0` 中位数 **`16.747 ns`**（`N=828`）；`status == 512` 修正前 **`0.739 ns`**（`N=354`），加 `+16.0 ns` 后 **`16.739 ns`**（差仅 **`0.008 ns`**）。
> - 所有 55 个 xAOD 扫描中，`VetoNuWaveformHits` 的 `wave_vector()` 非空覆盖率均为 **`100.0%`**。

### 3.2 与 Ben Wilson 参考文件逐事件数值交叉验证

在 5 个参考目录（及 `15838` 早期 CaloNu 过渡 run）的完全相同 xAOD 事例上进行 `(run, eventID)` 逐事件匹配（输出：`validation_v1/validation_summary.json`）：

| 验证样本 | Ben 参考目录 | 匹配事例数 / 有效时间窗事例数 | `fallback_mode` 一致率 | 精确覆写模式 `\|Δq0\| < 10⁻³ pC` 比例 | 精确覆写模式 `max \|Δq0\|` [pC] | `valid` 分支一致率 (`< 10⁻³ pC`) | 自动模式 `med \|Δq0\|` [pC] / `< 1 pC` 比例 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `2022_008090` (`r0021`) | `waveforms-2022` | `5000 / 548` | **`100.00%`** (`5000/5000`) | **`100.00%`** | `7.6 × 10⁻⁵` | **`100.00%`** | `0.000` / `99.82%` |
| `2023_011703` (`r0021`) | `Waveforms-2023` | `5000 / 626` | N/A（旧版无该分支） | **`100.00%`** | `1.22 × 10⁻⁴` | N/A | `0.000` / `99.36%` |
| `2024_pre_015387` (`r0022`) | `waveforms-5FilesPerChunk` | `5000 / 2714` | N/A | **`100.00%`** | `3.66 × 10⁻⁴` | N/A | `0.000` / `99.82%` |
| `2024_cn_early_015838` (`r0022`) | `WaveformsCaloNuPeriod` | `5000 / 3036` | N/A | **`100.00%`**（覆写 `35.783`） | `3.81 × 10⁻⁴` | N/A | 注：`15838` 实为 Dig1 Veto (`33.06 ns`)，自动模式 `33.165` 更准确 |
| `2024_CaloNu_015849` (`r0022`) | `WaveformsCaloNuPeriod` | `8000 / 4837` | N/A | **`99.98%`**（`100% < 1.1×10⁻³`） | `1.10 × 10⁻³` | N/A | **`0.000`** / **`100.00%`** |
| `2024_post_016929` (`r0022`) | `waveforms_postCaloNu` | `5000 / 3008` | **`100.00%`** (`5000/5000`) | **`100.00%`** | `4.88 × 10⁻⁴` | **`100.00%`** | `0.000` / `99.83%` |

### 3.3 修复版（`fix_v1` 自动模式）vs 未改基线（`latest 7fd40df8`）关键改进数据

| 测试样本与子集 | 事例数 `N` | `Raw` 中位数 [pC] | 基线 `latest` 中位数 [pC] | 修复 `fix_v1` 中位数 [pC] | 基线 `< 40 pC` 误判率 | 修复 `< 40 pC` 误判率 |
| --- | --- | --- | --- | --- | --- | --- |
| `2024_CaloNu_015849`：`Veto2 st==512` 且无好 Preshower | `299` | `77.65` | `77.65`（丢弃 `512`，退化为 raw） | **`66.02`**（恢复 Veto2 `+16 ns`） | — | — |
| `2024_CaloNu_015849`：`Veto2 st==0` muon | `2056` | — | 较修复偏高 `+1.87 pC`（`88.5%` 偏差 `>1 pC`） | 与 Ben 参考 `0.000 pC` 对齐 | — | — |
| `2025_20012`：Preshower fallback muon（`PS0,1` 失效、仅 `PS2,3`） | `67` | — | `11.30` | **`87.10`** | **`58.21%`** | **`1.49%`** |
| `2025_20481`：Preshower fallback muon（4-PMT 正常） | `67` | — | `0.18` | **`87.57`** | **`95.52%`** | **`0.00%`** |
| `2026_22404`：Preshower fallback muon（4-PMT 正常） | `135` | — | `0.33` | **`92.26`** | **`97.04%`** | **`0.00%`** |
| `2026_23176`：Preshower fallback muon（2026 后期质子束） | `97` | — | `0.23` | **`87.40`** | **`91.75%`** | **`0.00%`** |
| **2025–2026 Preshower fallback muon 合并** | **`366`** | **`102.82`** | **`0.49`** | **`88.81`** | **`88.25%`**（`323/366`） | **`0.27%`**（`1/366`） |
| `MC24_PG_100044`（穿行 muon，`fallback > 0`） | `2976` | `103.32` | `87.02` | `87.02` | `0.00%` | `0.00%` |
| `MC24_Genie_200139`（中微子相互作用，`fallback > 0`） | `718` | `<40 pC: 94.57%` | `<40 pC: 95.82%` | `<40 pC: 95.82%` | — | — |

### 3.4 关键验证图表

![2022-2026 & MC Scintillator Timing Timeline](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/figures/fig01_timing_diff_timeline_2022_2026_mc.png)

![2nd-Digitizer status==512 +16 ns Alignment](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/figures/fig02_status512_16ns_alignment.png)

![Event-by-Event Parity vs Ben Wilson Reference Files](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/figures/fig03_parity_vs_ben_wilson.png)

![Fix vs Baseline Smoking Guns across 2024 CaloNu, 2025-2026, and MC24](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/figures/fig04_fix_vs_baseline_smoking_guns.png)

### 3.5 执行命令、单元测试与产物清单

- **构建命令**：
  ```bash
  python3 /eos/home-x/xcheng/FASER/geoRotation_test/scripts/build_fix_v1_ntuple.py \
    --source /eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1 \
    --latest /eos/home-x/xcheng/FASER/geoRotation_calypso_latest \
    --output /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/software_ntuple_fix_v1 \
    --jobs 4
  ```
- **端到端验证命令**：
  ```bash
  python3 /eos/home-x/xcheng/FASER/geoRotation_test/scripts/validate_reduced_charge_all_periods.py \
    --overlay /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/software_ntuple_fix_v1 \
    --env-script /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/cached_athena_env.sh \
    --outdir /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v1 \
    --workers 10
  ```
- **回归单测**：`python3 -m unittest tests/test_reduced_charge_logic.py -v`（6/6 通过）。
- **紧凑产物路径**（总计 `20 MB`，已清理临时验证 ROOT）：
  - [timing_scan_phys.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_scan_phys.json)
  - [timing_scan_xaod.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/timing_scan_xaod.json)
  - [validation_summary.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/validation_v1/validation_summary.json)
  - [software.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/software_ntuple_fix_v1/software.json)
  - [waveform_reduced_charge_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/figures/waveform_reduced_charge_report.pdf)

---

## 4. 阶段性结论

### 已由实测数据证实的事实
1. 修复后的 `NtupleDumper` 在 `geoRotation_calypso_latest_fix_v1`（`a4ecefea`）中编译通过且在 2022–2026 与 MC24 全部 12 组真实 xAOD 样本上零错误运行。
2. 在精确时间差覆写下，与 Ben Wilson 5 套参考目录的逐事件数值一致性达到 **`100.00%`**（`< 10⁻³ pC`，2024 CaloNu 15849 为 `99.98% < 10⁻³ pC`、`100.00% < 1.1×10⁻³ pC`），`fallback_mode` 与 `valid` 分支一致率均为 **`100.00%`**。
3. 在 2024 CaloNu（`15845–16927`），修复版正确恢复全部 `status == 512`（`+16.0 ns`）的 Veto St2（`ch 22, 23`）击中，消除 `status == 0` 事例 `+1.87 pC` 系统偏置及无好 Preshower 时的 raw-charge 退化；同时自动区分 `15821–15844`（Veto 仍在 Dig1 `ch 6, 7`，实测 `33.06 ns`）。
4. 在 2025–2026 数据中，通过引入 Preshower 新时间差（Dig1 `13.28 ns`）与 Dig2 通道相对偏移（`+3.97 ns`，即 `17.25 ns`），即使在 `19480–20080` 期间 `PS0,1` 发生 `GFIT_FAILED` 仅剩 `PS2,3` 时也能保持无偏；366 个 2025–2026 Preshower-fallback 穿行 muon 的 `< 40 pC` 误判率由基线的 **`88.25%`** 降至 **`0.27%`**，中位数电荷由 `0.49 pC` 恢复至 **`88.81 pC`**。
5. 通过通道映射（`4 -> 0, 5 -> 1`）与空波形保护，彻底消除了原代码在 `VetoNuWaveformHits->size() >= 3`（实测占比 `0.5%–2.4%`）时的数组越界未定义行为。

### 说明与边界条件（非推测，明确标注适用范围）
1. `run == 18193`（2025 启动单 run 调试期）在 EOS 公开 `phys`/`rec` 中无含穿行 muon 的文件，其自动查找值（`Veto2 = 35.783 ns`, `PS = 25.28 ns`）系依据 `CABP200.db` 与 `ALLP200.db` 中 `TriggerOffset` 跳变（`-28 ns` vs `-40 ns`）推导，标注为未直接测量项。
2. `MC22` 旧数字化（`WAVE-Digitization-01`，`cb_mean = 815/815/846 ns`）的实测时间差为 `2.34 / 32.44 ns`，与 `MC24`（`27.21 / 51.41 ns`）不同；当前自动模式默认对应 `MC24`，若处理历史 `MC22` 样本可通过 `--vetoVetoNuTimeDiff 2.34 --preshowerVetoNuTimeDiff 32.44` 显式指定。
3. 2026 末期 `run >= 23359` 实测 `VetoNu ch 4, 5` 无束流信号（`q_p95 < 0.8 pC`，疑似 VetoNu 断高压或移除），不属于软件缺陷。

---

## 5. 下一步计划

1. 若将 `fix/latest-seed-and-momentum-v1` 的 `FaserActsKalmanFilter` 与 `NtupleDumper` 变更合并至主干生产构建，可在完整 `geoRotation_calypso_latest_fix_v1` worktree 中执行一次全量 `build_local.sh` 与 `verify_local.sh` 集成验收。
2. 如需对 2025–2026 全年大规模样本批量生成带新 `VetoNu{0,1}_reduced_charge` 与 `VetoNu_reduced_charge_fallback` 分支的 PHYS ntuple，可直接通过 HTCondor 调用 `software_ntuple_fix_v1`（无需任何手工时间差参数）。
