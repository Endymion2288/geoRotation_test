# 阶段 08 报告：最新 GeoRotation Calypso 双径迹 Overlay 追踪性能验证与旋转种子 `covPhi` 缺陷定位

- **日期**：2026-11-25
- **参考材料**：[Overlay_TrackingMeeting_9June_Pawan.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test/docs/Overlay_TrackingMeeting_9June_Pawan.pdf)（Pawan，2026-06-09 FASER Tracking Meeting，Slides 1–16）
- **核心软件与 Commit**：
  - **Latest GeoRotation (`latest_georot` / `fix_v1`)**：`/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1`，分支 `fix/latest-seed-momentum-for-sazhang`，HEAD `218813c2bcfc656fb8c395a167affcc033250867`（clean），预编译库 `/eos/home-x/xcheng/FASER/geoRotation_test_results/calypso_sazhang_mr/run`（Athena `24.0.41`，`x86_64-el9-gcc13-opt`）
  - **Pre-fix GeoRotation (`t2b_7fd40df8`)**：`/eos/home-x/xcheng/FASER/geoRotation_calypso_latest`，分支 `geoRotation_latestMaster`，HEAD `7fd40df8cd3b60e81b4095a4f84310c4a54700a0`，预编译库 `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/run`
  - **旋转材料图**：[material-maps-alma9-rotated-rowmajor.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/calypso_master_sync_v1/material/material-maps-alma9-rotated-rowmajor.json)（与 `218813c2` 和 `7fd40df8` 配套的 row-major 解压版）
  - **未旋转基线（EOS 冻结 ntuples）**：
    - `Baseline-May30`（`/eos/experiment/faser/data0/overlay/<dataset>/OverlaidNTPs_bug/`）：`Calypso-22.0.49` `two-track-overlay-2025/54a9026b`（2026-05-30，June 9 报告原始样本，旧版 `longTracks` 未过滤 1/2 站径迹）
    - `Baseline-Jul8`（`/eos/experiment/faser/data0/overlay/<dataset>/OverlaidNTPs/`）：同一批 May 30 xAOD（`54a9026b`），于 2026-07-08 使用 `811b81b8` 重导并启用 3-station `checkIsLongTrack`（与 `218813c2` 严格同口径，作为**主基线**）
  - **真实信号参考 (`110038`, `A' 100 MeV, ε=1e-5`)**：未旋转参考 `/eos/experiment/faser/data0/sim/mc24/foresee/110038/phy/master_dev/`（10 文件，20,000 事件，2025-07-01 `master_dev`），旋转版由 `rdo/master_dev/` 用 `218813c2` 原位重建
- **仓库脚本、测试与修复补丁**：
  - 重建与坐标闭合：[overlay_reco_ntuple.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/configs/overlay_reco_ntuple.py)、[overlay_acts_to_faser.cpp](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_acts_to_faser.cpp)、[overlay_worker.sh](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/overlay_worker.sh)、[overlay_signal_worker.sh](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/overlay_signal_worker.sh)
  - 任务冻结、提交与编译：[submit_overlay_2track.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/submit_overlay_2track.py)、[build_overlay_analysis.sh](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/build_overlay_analysis.sh)
  - C++17 多线程分析与 22 图绘制：[overlay_analyze.cpp](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_analyze.cpp)、[overlay_signal_110038.cpp](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_signal_110038.cpp)、[overlay_tier2_analyze.cpp](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_tier2_analyze.cpp)、[overlay_figures.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_figures.py)、[overlay_tier2_figures.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_tier2_figures.py)、[overlay_plot_style.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_plot_style.py)
  - 单元测试：[test_overlay_frame.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/tests/test_overlay_frame.py)、[test_overlay_metrics.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/tests/test_overlay_metrics.py)、[test_overlay_figures.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/tests/test_overlay_figures.py)（13/13 PASS）
  - 最小修复补丁：[stage08_ckf2_rotated_seed_covphi.patch](file:///eos/home-x/xcheng/FASER/geoRotation_test/patches/stage08_ckf2_rotated_seed_covphi.patch)（在 `CKF2Config.py` 非 `legacySwappedCovLoc` 旋转分支将 `sigma_phi_rotated = sigma_theta = 2.9e-4`，使 `covPhi = 8.41e-6 rad²`）
- **EOS 结果目录**：
  - 主运行（323 overlay + 10 `110038` + 22 页图集）：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/`
  - Tier-2 归因运行（各 36 shards）：`run_20261125_v1_t2a_covphi/`（`218813c2 + covPhi=8.41e-6`）、`run_20261125_v1_t2b_7fd40df8/`（`7fd40df8`）、`run_20261125_v1/tier2_attribution/`（4 路合并结果）、`run_20261125_v1/tier2_110038_check/`（`110038` 200 事件对照）
- **HTCondor 提交目录**：`~/condor/geoRotation_test/stage08_overlay_2track_run_20261125_v1{,_t2a_covphi,_t2b_7fd40df8}/`
- **独立审计报告**：Gate 1 闭合、Gate 2 完备性及 Gate 3 全部 22 图与数值复核均通过（[completeness_audit.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/gate2_audit/completeness_audit.json)、[signal_audit.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/gate1_closure/signal_110038_20evt/signal_audit.json)）。

---

## 阶段目标

1. **完整复现与直接对比**：在最新旋转版 Calypso（`fix/latest-seed-momentum-for-sazhang`，`218813c2`）上全量重建 EOS `OverlaidRDOs`（`008749`、`010698`、`Pos_008749-Neg_010698` 与 `10004X` 共 323 文件、6,354,483 事件）及真实信号样本 `110038`（`A' 100 MeV, ε=1e-5`，20,000 事件），按 Pawan（2026-06-09）全部核心物理量（$\Delta R_1, \Delta R_0, |\Delta Y_0|, |\Delta Y_1|, |\Delta X_0|, \Delta X_0, |\Delta X_1|, \Delta X_1, \Delta Y_1, \Delta Y_0, \text{Track } P_{z0}$）与未旋转基线（`Baseline-Jul8` 及 `Baseline-May30`）执行逐事件精确配对比较。
2. **解答报告遗留的追踪机制问题**：定量分解（a）高原区（$\sim 92$–$94\%$）与单径迹效率（$\sim 99\%$）的差距来源；（b）大 $|\Delta X|$ 效率下降的几何与重建成因；（c）小 $\Delta R_1 < 5\text{ mm}$ 效率下降与径迹多重数分裂；（d）MC/Data `Pz0` 谱差异重加权对高原比值的影响。
3. **异常归因与最小修复验证**：定位全量运行中新发现的 $|\Delta Y_0| \approx 16\text{ mm}$ 局部效率深坑（$-20.9\text{ pp}$ Data / $-27.2\text{ pp}$ MC）的精确代码根因，并通过 4 路受控 HTCondor 归因实验（Tier-2，36 文件对、908,731 事件）与 `110038` 回归检查验证最小修复的因果闭合效果。

---

## 工作内容

### 1. 坐标系闭合与单步零 xAOD 重建管线

- **坐标系问题与闭合方案**：检查 `TwoTrackOverlayNTuplerAlg.cxx` 确认，在旋转几何（`geometryRotation = True`）下该算法将重建径迹 `Track_{x,y,z,px,py,pz}{0,1}` 直接以旋转 ACTS 坐标系 $(x_A, y_A, z_A) = (z_F, y_F, -x_F)$ 输出，而原始单径迹真值 `og_Track_*` 与 `TrackSegment_*` 位于标准 FASER 坐标系（$z$ 沿束流，$x$ 水平，$y$ 竖直，磁场沿 $+x$）。为在不重编 Calypso 的前提下保持 100% 真实组件运行，设计了单步零 xAOD 管线：
  1. [overlay_reco_ntuple.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/configs/overlay_reco_ntuple.py) 在同一 Athena 进程内执行 `FaserActsKalmanFilter`（`--noBackward`，前向 CKF + refit）并直接将 `TwoTrackOverlayNTuplerAlg` 接入内存 `StoreGate`，不落盘任何中间 xAOD；
  2. [overlay_acts_to_faser.cpp](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_acts_to_faser.cpp) 在 worker 本地 scratch 内对 `Track_{x,y,z,px,py,pz}{0,1}` 施加精确旋转逆变换 $(x_F, y_F, z_F) = (-z_A, y_A, x_A)$，同时对全量径迹施加物理边界硬断言（$\text{med}(z_0) \in [-20, 80]\text{ mm}$、$\text{med}(z_1) \in [2350, 2500]\text{ mm}$、`0` 条径迹满足 $|z_1| < 500\text{ mm}$、$\ge 99.9\%$ 径迹满足 $z_0 \in [-100, 200]\text{ mm}$ 且 $z_1 \in [2200, 2650]\text{ mm}$ 与 $p_{z0}, p_{z1} > 0$，外加 `og_longTracks == 2`、零非有限值与 `CollectionTree` 逐条对齐），并写入 `meta_frame="FASER"` 及软件 provenance 元数据。
- **Conditions 与上游代码一致性**：核对发现 `218813c2` 将 `TI12MC04` 默认 conditions 标签改为 `OFLCOND-FASER-06`，而 EOS `OverlaidRDOs` 生产日志使用 `OFLCOND-FASER-05`。本阶段全部 overlay 任务均显式传参 `--cond OFLCOND-FASER-05`，确保与基线对齐至完全相同的 `TRACKER-ALIGN-05`（Data，`CONDBR3`）与 `TRACKER-ALIGN-02`（MC，`OFLP200`）。另对 `git diff 54a9026b..b28f02da` 审计确认，2026-05 至 2026-11 之间 `CKF2`、`CircleFitTrackSeedTool`、`KalmanFitterTool`、`SegmentFitAlg`、`GhostBusters` 与 `FaserActsGeometry` 零算法改动（实测 `107,421` 个 segment 零差异）。
- **Gate 1 闭合验证**：在 `008749` 与 `10004X` 各 100 事件（共 200 事件）上对比单步零 xAOD（`--noBackward`）与标准两步（写盘 xAOD 再 dump ntuple），48/48 分支在 `float32` 精度下逐比特一致（`float64` 最大相对差 $< 5.9 \times 10^{-8}$，源自 POOL xAOD T/P `float32` 序列化），坐标转换零越界；另对 `110038` 20 事件 dry-run 核对确认 `20/20` `eventID` 与 `phy/master_dev` 精确对齐、`16/20` `longTracks == 2` 一致、全部 `44` 条重建径迹均由 `NtupleDumper` 自动输出于标准 FASER 坐标系（$z_0 \in [9.50, 21.42]\text{ mm}$、$z_1 \in [2398.33, 2465.30]\text{ mm}$），且零残留 xAOD。

### 2. 数据冻结、效率定义与符号约定

- **冻结掩码**：EOS `OverlaidRDOs/` 共 323 个文件（`008749` 81、 `010698` 81、`Pos_008749-Neg_010698` 81、`10004X` 80，共 `6,354,483` 事件）。其中 `010698/PosBin2_NegBin3`（`14,884` 事件）在 `OverlaidNTPs` / `OverlaidNTPs_bug` 中缺失，因此**主 3 版本严格配对掩码**冻结为 322 个共同文件对（`6,339,599` 事件：`All3Data` = 242 文件、`3,487,558` 事件；`10004X` MC = 80 文件、`2,852,041` 事件），`010698_PosBin2_NegBin3` 单独作为辅助项报告。
- **物理量、坐标系与效率定义**（全部在标准 **FASER 坐标系**，位置单位 **mm**，ntuple 动量分支单位 **MeV**，动量重加权与动量谱图转换为 **GeV**）：
  - 相对分离量按原始单径迹定义（`og_Track[0] - og_Track[1]` = 正径迹减负径迹，Station 1 记下标 `0`，Station 3 记下标 `1`）：$\Delta X_{0,1}, \Delta Y_{0,1}, \Delta R_{0,1} = \sqrt{\Delta X_{0,1}^2 + \Delta Y_{0,1}^2}$。
  - **主效率（Pawan 口径 `eff_pawan_c2R`）**：分母为全部合法 overlay 事件（原始单径迹满足 $|R| < 95\text{ mm}$）；分子要求 `longTracks == 2`（3 站长径迹）且两条重建径迹均满足 $\chi^2/\text{ndof} \le 25$ 且 $R_0, R_1 < 95\text{ mm}$。同时并行计算纯多重数指标 `eff_lt2`（`longTracks == 2`，无 $\chi^2/R$ 截断）、`eff_ge2`（`longTracks >= 2`，无 $\chi^2/R$ 截断）以及几何匹配指标 `eff_matched`（`longTracks >= 2` 且存在两条不同重建径迹同时通过 `c2R` 质控并与两条 `og_Track` 在 St.1 和 St.3 均满足 1-to-1 $\Delta r < 5\text{ mm}$ 匹配）。
  - 逐事件配对差值 $\Delta \varepsilon = (g - l)/N$ 与 McNemar 标准误差 $\sigma = \sqrt{(g + l) - (g - l)^2 / N} / N$（$g$ 为新版通过/旧版失败事件数，$l$ 为旧版通过/新版失败事件数）。

### 3. HTCondor 批量执行、C++17 多线程分析与 Tier-2 归因实验

- **批量提交策略**：每次提交前查询普通两位 EL9 `bigbirdNN` schedd 队列负载，按 `TotalRunningJobs + TotalIdleJobs` 最小原则选中 `bigbird20.cern.ch`（提交文件全部 $< 1\text{ KB}$ 置于 `~/condor/`，xAOD 零落盘，紧凑结果直接回传 EOS）。4 个 HTCondor cluster 共 **405/405 作业一次性全部成功（0 失败、0 重试）**。
- **高性能分析与可视化**：
  - [overlay_analyze.cpp](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_analyze.cpp)（C++17 + OpenMP 16 线程 + ROOT `TTreeReader`）：一次扫描 `322 × 3` 个 ROOT 文件（19,018,797 条目），在 `44.55 s` 内完成 `6,339,599/6,339,599` 逐事件 `og_Track` 身份硬校验（0 失败）、全曲线与 McNemar 差值、逐阶段损失（`stage_of_loss`）、大 $|\Delta X|$ 分层、`Pz0` 重加权及动量残差归约。
  - [overlay_tier2_analyze.cpp](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_tier2_analyze.cpp)（C++17 + OpenMP 16 线程）：在 `8.60 s` 内完成 4 路（`Baseline-Jul8`、`t2b_7fd40df8`、`Latest-218813c2`、`t2a_covphi`）36 文件对（`908,731` 匹配四元组，0 身份失败）的逐事件配对归因。
  - [overlay_figures.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_figures.py) 与 [overlay_tier2_figures.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/overlay_tier2_figures.py)：产出 22 张英文发布级图（PNG + 单图 PDF + 22 页统一报告 PDF + `bin_counts/` CSV + `figure_index.json`）。

关键复现与验证命令：

```bash
cd /eos/home-x/xcheng/FASER/geoRotation_test
# 1. 提交主运行（323 overlay + 10 signal 110038）与 Tier-2 归因运行（36 shards × 2）
python3 scripts/submit_overlay_2track.py --run-id run_20261125_v1 --include-signal-110038 --submit
python3 scripts/submit_overlay_2track.py --run-id run_20261125_v1_t2a_covphi --subset tier2_36 --seed-cov-phi 8.41e-6 --submit
python3 scripts/submit_overlay_2track.py --run-id run_20261125_v1_t2b_7fd40df8 --subset tier2_36 \
    --calypso-source /eos/home-x/xcheng/FASER/geoRotation_calypso_latest \
    --calypso-run /eos/home-x/xcheng/FASER/geoRotation_calypso_latest/run \
    --calypso-commit 7fd40df8cd3b60e81b4095a4f84310c4a54700a0 --submit

# 2. 编译 C++17 多线程分析器并生成全量 22 页图集与单元测试验收
bash analysis/build_overlay_analysis.sh
python3 analysis/overlay_figures.py \
    --summary /eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/analysis/overlay_summary.json \
    --signal /eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/analysis/signal_110038.json \
    --out-dir /eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures --mode full
python3 analysis/overlay_tier2_figures.py \
    --tier2-summary /eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/tier2_attribution/tier2_attribution_summary.json \
    --out-dir /eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures
python3 -m unittest tests/test_overlay_frame.py tests/test_overlay_metrics.py tests/test_overlay_figures.py -v
```

---

## 详细数据

### 1. HTCondor Schedd 选择与执行资源汇总

全部提交均在 `lxplus-gpu` 发起，提交前实时查询 EL9 `bigbirdNN` schedd 候选：主运行提交时刻（`2026-11-25 16:01 UTC`）最空闲前三位为 `bigbird20`（`275 run / 557 idle / 832 total`）、`bigbird21`（`2015/266/2281`）、`bigbird17`（`306/2268/2574`）；Tier-2 提交时刻（`18:52 UTC`）最空闲前三位为 `bigbird20`（`169/655/824`）、`bigbird21`（`1291/237/1528`）、`bigbird17`（`235/2076/2311`），均选择 `bigbird20.cern.ch`（容量 `10,055`）。

| 任务批次 | Cluster ID | Shards (成功/总) | 事件总数 | Wall 中位数 / 最大值 | RSS 中位数 / 最大值 | EOS 紧凑产物 | 提交文件大小 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 主运行 Overlay (`218813c2`, 默认 `covPhi=0.1089`) | `1587793` | `323 / 323` | `6,354,483` | `122.0 s` / `668.0 s`（全批 `13m49s`） | `1,059.8 MB` / `1,733.0 MB` | `7.654 GB` | `909 B` |
| 主运行 Signal `110038` (`218813c2`, `reco+ntp`) | `1587794` | `10 / 10` | `20,000` | `5,120.0 s` / `7,607.0 s` | `1,343.9 MB` / `1,431.4 MB` | `119.25 MB` | `945 B` |
| Tier-2 Arm A (`218813c2 + covPhi=8.41e-6`) | `1587796` | `36 / 36` | `908,731` | `193.0 s` / `538.0 s` | `1,269.0 MB` / `1,710.0 MB` | `1.067 GB` | `963 B` |
| Tier-2 Arm B (`7fd40df8` pre-fix GeoRot) | `1587797` | `36 / 36` | `908,731` | `238.0 s` / `851.0 s` | `1,274.7 MB` / `1,733.6 MB` | `1.064 GB` | `973 B` |

Gate 2 核对结果：主运行 `323/323` + `10/10` 及 Tier-2 `72/72` ROOT 文件无 zombie、零残留 xAOD、`actsToFaser` 坐标越界 `0 / 12,175,413` 条径迹；独立 PyROOT `RDataFrame` 复核与 C++17 汇总 JSON 逐 bin 差值全部为 **`0`**。

### 2. 全量 322 文件对（6,339,599 事件）总览与配对对比

数据来源：[overlay_summary.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/analysis/overlay_summary.json)。主指标为 `eff_pawan_c2R`（`lt==2`, $\chi^2/\text{ndof}\le 25$, $R_{0,1}<95\text{ mm}$）；$\Delta$ 与误差为逐事件配对 McNemar 统计。

| 数据集 | 文件对 | 事件数 ($N$) | Baseline-May30 (%) | Baseline-Jul8 (%) | Latest (`218813c2`) (%) | Latest − Jul8 (pp) | Gain ($g$) / Loss ($l$) | Latest − Jul8 (`ge2`, pp) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `008749` (Data) | 81 | 1,171,836 | 92.2211 | 92.0538 | 92.5052 | **`+0.4513 ± 0.0189`** | 27,060 / 21,771 | +0.2815 |
| `010698` (Data) | 80 | 1,143,938 | 92.8532 | 92.7008 | 93.2193 | **`+0.5185 ± 0.0190`** | 26,593 / 20,662 | +0.3422 |
| `Pos_008749-Neg_010698` | 81 | 1,171,784 | 92.5750 | 92.4158 | 92.8291 | **`+0.4134 ± 0.0188`** | 26,733 / 21,889 | +0.2449 |
| `008749 + 010698` | 161 | 2,315,774 | 92.5333 | 92.3734 | 92.8579 | **`+0.4845 ± 0.0134`** | 53,653 / 42,433 | +0.3115 |
| **`All3Data` (合并 Data)** | **242** | **3,487,558** | **92.5473** | **92.3877** | **92.8483** | **`+0.4606 ± 0.0109`** | **80,386 / 64,322** | **+0.2891** |
| **`10004X` (MC)** | **80** | **2,852,041** | **93.3298** | **93.3028** | **91.9505** | **`-1.3523 ± 0.0094`** | **17,118 / 55,685** | **-1.4104** |

- 辅助单文件 `010698/PosBin2_NegBin3`（仅 `Latest` 有，无基线对应）：`13,900 / 14,884 = 93.3889%`。
- `Baseline-May30` 与 `Baseline-Jul8` 差值极小且集中在小 $\Delta R$（Data 全局 $+0.1596\text{ pp}$，MC 全局 $+0.0270\text{ pp}$），确认 `811b81b8` 的 3-station 修正仅去除少量短径迹污染（详见图 20）。

### 3. 关键相空间区域与 `$|\Delta Y_0| \approx 16\text{ mm}$` 效率深坑

全量曲线检查立即暴露出一个未见于未旋转基线的尖锐特征：在 `Latest-218813c2` 中，双径迹效率在 $|\Delta Y_0| \in [13, 18]\text{ mm}$（及对应的 $\Delta R_0 \approx 16$–$17\text{ mm}$）处出现深度超过 **20 pp** 的局部效率深坑；而在坑外区域（如大 $|\Delta X_0|$ 与小 $|\Delta Y_0| < 2\text{ mm}$），`Latest-218813c2` 反而明显优于基线：

| 相空间区域 | 样本 | 事件数 ($N$) | Baseline-Jul8 (%) | Latest (`218813c2`) (%) | Latest − Jul8 (pp) | 特征说明 |
| --- | --- | --- | --- | --- | --- | --- |
| 小间距 $\Delta R_1 < 5\text{ mm}$ | `All3Data` / `10004X` | 11,887 / 10,504 | 90.250 / 89.071 | 88.458 / 85.558 | `-1.7919 ± 0.2935` / `-3.5129 ± 0.3248` | `lt0` 与 `ckf_lost` 上升，`ge3` 反而减少 |
| 高原区 $\Delta R_1 \in [30, 70]\text{ mm}$ | `All3Data` / `10004X` | 1,301,156 / 1,084,164 | 92.965 / 94.074 | 93.148 / 92.541 | `+0.1838 ± 0.0192` / `-1.5329 ± 0.0160` | MC 受深坑拖累下降，Data 因大 $|\Delta X|$ 收益微升 |
| 大水平间距 $|\Delta X_0| \in [70, 100]\text{ mm}$ | `All3Data` / `10004X` | 490,258 / 393,193 | 90.524 / 90.861 | 91.662 / 90.262 | `+1.1374 ± 0.0253` / `-0.5992 ± 0.0198` | Data `ckf_lost` 减半（`1.90% → 0.91%`） |
| **深坑极小 bin ($|\Delta Y_0| = 15.25$ / $15.75\text{ mm}$)** | `All3Data` / `10004X` | 22,046 / 18,195 | **91.87 / 93.46** | **71.01 / 66.24** | **`-20.86 ± 0.30` / `-27.22 ± 0.33`** | 下文 Tier-2 证实由 `CKF2Config.py` `covPhi` 缺陷唯一导致 |
| $\Delta R_0$ 深坑极小 bin ($17.25$ / $16.25\text{ mm}$) | `All3Data` / `10004X` | — | — | — | `-11.98 ± 0.43` / `-20.18 ± 0.52` | $|\Delta Y_0| \approx 16\text{ mm}$ 坑在 $\Delta R_0$ 上的投影 |

### 4. 逐阶段损失分解与 Pawan 报告四个遗留问题的定量回答

下表给出全量 322 文件对在四个代表区域的逐阶段损失占分母百分比（`pass_c2R + lt2_fail_c2R + lt_ge3 + lt0 + miss_seg + gb + ckf_lost + unm = 100%`）：

| 版本 / 样本 | 区域 | `pass_c2R` (%) | `c2R` 失败 (%) | `lt>=3` (%) | `lt==0` (%) | `miss_seg` 总 (St1 / St2 / St3 / 多站) (%) | `GhostBusters` (%) | `all3gb → ckf_lost` (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `Jul8` / `All3Data` | `plateau [30,70]` | 92.965 | 0.274 | 0.070 | 0.106 | **3.985** (0.10 / 1.48 / 2.07 / 0.33) | 0.042 | **2.469** |
| `Jul8` / `10004X` | `plateau [30,70]` | 94.074 | 0.079 | 0.112 | **1.619** | **3.504** (0.11 / 1.53 / 1.65 / 0.21) | 0.042 | **0.489** |
| `Latest` / `All3Data` | `plateau [30,70]` | 93.148 | 0.094 | 0.071 | 0.490 | **3.934** (0.10 / 1.44 / 2.10 / 0.29) | 0.041 | 2.194 |
| `Latest` / `10004X` | `plateau [30,70]` | 92.541 | 0.033 | 0.113 | 2.251 | **3.434** (0.12 / 1.47 / 1.66 / 0.19) | 0.042 | 1.561 |
| `Jul8` / `All3Data` | `small ΔR1<5` | 90.250 | 0.454 | **4.408** | 0.126 | 1.312 (0.21 / 1.03 / 0.05 / 0.02) | 0.042 | 3.373 |
| `Jul8` / `10004X` | `small ΔR1<5` | 89.071 | 0.181 | **6.093** | 1.447 | 1.285 (0.23 / 0.97 / 0.07 / 0.02) | 0.029 | 1.828 |
| `Latest` / `All3Data` | `small ΔR1<5` | 88.458 | 0.067 | **2.751** | 1.775 | 1.296 (0.19 / 1.05 / 0.03 / 0.02) | 0.025 | 5.578 |
| `Latest` / `10004X` | `small ΔR1<5` | 85.558 | 0.086 | **4.227** | 3.646 | 1.314 (0.30 / 0.90 / 0.10 / 0.02) | 0.019 | 5.093 |
| `Jul8` / `All3Data` | `large \|ΔX0\| [70,100]` | 90.524 | 0.385 | 0.017 | 0.125 | **6.958** (0.07 / 2.68 / 3.58 / 0.63) | 0.031 | 1.900 |
| `Jul8` / `10004X` | `large \|ΔX0\| [70,100]` | 90.861 | 0.084 | 0.023 | 1.595 | **7.018** (0.07 / 3.17 / 3.33 / 0.45) | 0.041 | 0.319 |
| `Latest` / `All3Data` | `large \|ΔX0\| [70,100]` | 91.662 | 0.282 | 0.020 | 0.092 | **6.948** (0.06 / 2.62 / 3.69 / 0.57) | 0.046 | **0.915** |
| `Latest` / `10004X` | `large \|ΔX0\| [70,100]` | 90.262 | 0.088 | 0.027 | 1.690 | **6.969** (0.07 / 3.08 / 3.39 / 0.43) | 0.060 | 0.872 |

由此直接回答 Pawan 报告中的四个关键问题：
1. **为何双径迹高原效率（$\sim 93$–$94\%$）低于单径迹效率（$\sim 99\%$，Slide 15）**：在未旋转基线高原区，总损失（`7.04%` Data / `5.93%` MC）由三大项完全闭合：（i）**上游 `SegmentFit` 缺段（`miss_seg`）**占 `3.99%`（Data）/ `3.50%`（MC），集中在 Station 2（`1.48–1.53%`）与 Station 3（`1.65–2.07%`）；（ii）**全 3 站过 `GhostBusters` 后的 CKF 丢失（`ckf_lost`）**占 `2.47%`（Data）/ `0.49%`（MC），这正是未旋转基线中 Data 低于 MC 约 `1.1 pp` 的真正原因；（iii）`10004X` 存在约 `1.62%` 的 `lt==0` 基底（Data 仅 `0.11%`）。`GhostBusters` 仅剔除 `0.042%`，绝非高原损失主因。
2. **为何大 $|\Delta X_0|$ 处效率下降（Slides 11–12）**：在大 $|\Delta X_0| \in [70, 100]\text{ mm}$ 区间，`miss_seg` 跃升至 **`6.96%`（Data）/ `7.02%`（MC）**（且在 $|\Delta X_0| \in [95, 100]\text{ mm}$ 达 `~8.1%`），几乎完全由 Station 2/3 缺段驱动。按最大单径迹横向位置 $\max(|x_0|)$ 分层检验（图 18）证实：在核心区 $\max(|x_0|) \in [0, 33]\text{ mm}$ 内，效率保持平坦且高达 **`95.37%`（Jul8 Data）/ `96.55%`（Jul8 MC）**（`Latest` Data 更达 **`96.13%`**）；但进入外缘区 $\max(|x_0|) \in [33, 66]\text{ mm}$ 与 $[66, 95]\text{ mm}$ 后，效率整体跌至 `91.06%` 与 `90.81%`（在 $|\Delta X_0|<20\text{ mm}$ 且两径迹均在外缘时仅 `83.88%` Data / `83.25%` MC）。因为大 $|\Delta X_0|$ 在几何上必然强迫至少一条径迹落入外层模块/芯片边界，从而表现为随 $|\Delta X_0|$ 的下降。
3. **小间距 $\Delta R_1 < 5\text{ mm}$ 下降机制**：在基线中主要由 `lt>=3` 径迹分裂/假径迹（`4.41%` Data / `6.09%` MC）造成（使宽松口径 `ge2` 在基线保持约 `95%`）。在 `Latest-218813c2` 中，`lt>=3` 显著压低至 `2.75%` / `4.23%`，但默认 `covPhi` 缺陷导致 `ckf_lost` 与 `lt0` 激增，掩盖了该项改善（Tier-2 修复后即完全反转）。
4. **`Pz0` 谱差异与重加权闭合（Slide 16）**：按 **40 个**对数等距 `og_p0` 分箱（`[50, 5000] GeV`）将 `10004X` MC 重加权至 `All3Data` 谱形（有效样本量 $N_{\text{eff}} = 750,516.34$，占 $2,852,041$），高原区 $\Delta R_1 \in [30, 50]\text{ mm}$ 的 MC/Data 效率比值变化**小于 $0.0013$**：
   - `Baseline-May30`：`1.01149`（未加权）→ `1.01210`（重加权）
   - `Baseline-Jul8`：`1.01333`（未加权）→ `1.01387`（重加权）
   - `Latest-218813c2`：`0.99429`（未加权）→ `0.99557`（重加权）
   因此，**基线中 ~1.3% 的 MC/Data 高原差值并非动量谱差异所致**，而是 Data 侧更高的 `ckf_lost`（`2.47%` vs `0.49%`）。在动量分辨率方面，匹配相对残差 $(P_{z0,\text{reco}}/P_{z0,\text{og}} - 1)$ 的 `iqr68`：`All3Data` 为 `0.0685 (Jul8) → 0.0753 (Latest)`，且在 `008749` 高动量区（`0.5–5 TeV`）旋转版显著收窄（`0.330→0.253`、`0.455→0.350`、`0.479→0.377`）；`10004X` 因基线与自身 `og_Track` 同引擎自比（`iqr68 = 0.0071`）而旋转版跨引擎对比为 `0.1103`（且在 Tier-2 `t2a_covphi` 中稳定为 `0.1104`）。

### 5. 真实信号样本 `110038`（`A' 100 MeV, ε=1e-5`）验证

数据来源：[signal_110038.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/analysis/signal_110038.json)（全量 10 文件、20,000 事件，其中 `17,724` 事件满足 truth St.1 & St.3 均在 $R < 95\text{ mm}$ 接受度内）及 [summary_110038_200evt.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/tier2_110038_check/summary_110038_200evt.json)（file0 前 200 事件 `covPhi` 单因素控制）：

| 样本与版本 | 总事件 | Truth 接受度 ($R<95$) | `2trk + OS + R<95` | `2trk + OS + c2R` | 相对总事件效率 (%) | 相对 Truth 接受度效率 (%) | $\chi^2/\text{ndof}>25$ 失败数 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 全量 20k：`master_dev` (2025-07 未旋转) | 20,000 | 17,724 | 9,362 | 9,244（其中 9,224 在 truth 接受度） | 46.81 ± 0.35 | **52.04 ± 0.38** | 118 |
| 全量 20k：`Latest-218813c2` (默认 `covPhi`) | 20,000 | 17,724 | 11,141 | 11,141（其中 11,122 在 truth 接受度） | 55.71 ± 0.35 | **62.75 ± 0.36**（**`+10.71 pp`**） | **0** |
| 200 事件对照：`ref_master_dev` | 200 | — | — | 97 (`lt==2`: 138, `c2R`: 126, `SS`: 29) | 48.5 | — | median $\chi^2/\text{ndof}=1.292$ |
| 200 事件对照：`Latest-218813c2` (默认) | 200 | — | — | 115 (`lt==2`: 142, `c2R`: 129, `SS`: 14) | 57.5 (`+9.0 ± 3.1 pp`) | — | median $\chi^2/\text{ndof}=1.047$ |
| **200 事件对照：`Latest + covPhi=8.41e-6`** | **200** | — | — | **119** (`lt==2`: 146, `c2R`: 134, `SS`: 15) | **59.5 (`+2.0 ± 1.2 pp` vs 默认, `+11.0 ± 3.1 pp` vs ref)** | — | median $\chi^2/\text{ndof}=1.049$ |

关键发现：`110038` 中 **98% 以上事件**具有 $|\Delta Y_0| < 2\text{ mm}$（位于 `16 mm` 效率深坑之外），因此在默认 `218813c2` 下就已展现旋转版在双径迹重建与异号电荷识别上的显著增益（`+10.71 pp`，其中 2-track 产额增加 `+570` 事件、异号识别率从 `73.6%` 升至 `83.8%`、`chi2>25` 失败清零）；叠加 `covPhi=8.41e-6` 后，在 200 事件样本上进一步净增 `+4` 个 `OS+c2R` 通过事件（`gain=5, loss=1`，无任何退化）。

### 6. Tier-2 根因定位与 4 路 908,731 事件归因对照

#### 6.1 `$|\Delta Y_0| \approx 16\text{ mm}$` 深坑的代码根因与数学机理

针对全量结果中的 $|\Delta Y_0| \approx 16\text{ mm}$ 深坑，我们开展了 9×9 `PosBin × NegBin` 失效率矩阵、Event 8（`10004X_PosBin0_NegBin0`，$|\Delta Y_0| = 16.0\text{ mm}$）逐层 `trackstates` 追踪、6 路单因素消融与 24 点 `(covPhi, covTheta)` 参数扫描，将缺陷唯一定位于 [CKF2Config.py](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/python/CKF2Config.py#L50-L68)：

```python
sigma_loc0 = 1.9e-2
sigma_loc1 = 9e-1
sigma_phi = 3.3e-2          # <--- 未旋转 1/sin(theta_F) 极点奇异遗留值！
sigma_theta = 2.9e-4        # <--- 真实物理角度分辨率
...
if legacy_swapped_cov_loc:
    default_cov_loc0 = initial_variance_inflation[0] * sigma_loc1 * sigma_loc1  # 81.0 mm^2
    default_cov_loc1 = initial_variance_inflation[1] * sigma_loc0 * sigma_loc0  # 0.0361 mm^2
else:
    default_cov_loc0 = initial_variance_inflation[0] * sigma_loc0 * sigma_loc0  # 0.0361 mm^2 (y_F)
    default_cov_loc1 = initial_variance_inflation[1] * sigma_loc1 * sigma_loc1  # 81.0 mm^2 (x_F)

seed_tool_kwargs.setdefault("covLoc0", default_cov_loc0)
seed_tool_kwargs.setdefault("covLoc1", default_cov_loc1)
seed_tool_kwargs.setdefault("covPhi", initial_variance_inflation[2] * sigma_phi * sigma_phi)      # = 100*(3.3e-2)^2 = 0.1089 rad^2
seed_tool_kwargs.setdefault("covTheta", initial_variance_inflation[3] * sigma_theta * sigma_theta) # = 100*(2.9e-4)^2 = 8.41e-6 rad^2
```

- **物理机制**：在未旋转 FASER 坐标系中，束流沿 $+z_F$（$\theta_F \approx 0.001$–$0.01\text{ rad}$ 极点奇异区），因此 $\sigma_{\phi,F} = \sigma_{\text{angle}} / \sin\theta_F \approx 3.3\times 10^{-2}\text{ rad}$ 仅因极点奇异而偏大。在旋转 ACTS 坐标系 $(x_A, y_A, z_A) = (z_F, y_F, -x_F)$ 中，束流沿 $+x_A$（$\sin\theta_A \approx 1$，极点奇异彻底消除），**$\theta_A \approx \pi/2 - t_x$ 为水平斜率，而 $\phi_A = \text{atan2}(p_{y,F}, p_{z,F}) \approx t_y$ 为精密竖直弯曲斜率**。沿用未旋转 `sigma_phi = 3.3e-2 rad` 导致竖直弯曲角方差被放大 **12,950 倍**（`covPhi = 0.1089 rad²`，$\sigma_\phi = 330\text{ mrad} = 18.9^\circ$，比真实 $\sigma_\theta = 2.9\text{ mrad}$ 宽 **114 倍**）。
- **Event 8（`10004X_PosBin0_NegBin0`，$|\Delta Y_0| = 16.0\text{ mm}$）因果闭合与深坑形成机理**：
  1. `218813c2`（`fix_v1`）修正了旋转 seed 平面上的 `covLoc0/covLoc1` 互换（使竖直方向 `covLoc0` 从 `81.0 mm²` 收紧至物理值 `0.0361 mm²`），但未同步修正 `covPhi = 0.1089 rad²`，导致 `covLoc0 = 0.0361 mm² ≪ L²·covPhi = 10.89 mm²` 暴露出来。
  2. 种子从起点（$z_F = -0.5\text{ mm}$）外推 $L \approx 10.0$–$10.89\text{ mm}$ 到达 Station 1 Plane 0（$z_F = 9.50 / 10.39\text{ mm}$）时，竖直预测不确定度为 $\sigma_y \approx L\cdot\sigma_\phi \approx 3.3$–$3.6\text{ mm}$，在 `chi2Max = 15` 下张开高达 $\pm\sqrt{15}\,\sigma_y \approx \pm 13$–$16\text{ mm}$ 的**竖直搜索窗**——精确解释为何效率深坑严格出现在 $|\Delta Y_0| \approx 15$–$17\text{ mm}$（而非 $|\Delta X_0|$）。
  3. 当两条径迹穿过同一 SCT module 且竖直相距 $|\Delta Y_0| \approx 16\text{ mm}$ 时，CKF 在首层拉入邻近径迹的 cluster（$\Delta y \approx \pm 16\text{ mm}$）；由于 `covLoc0 ≪ L²·covPhi`，Kalman 角度更新增益高达 **$K_\phi = L\cdot\text{covPhi} / (\text{covLoc0} + L^2\cdot\text{covPhi}) \approx 1/L \approx 0.10\text{ rad/mm}$**，一举将滤波角踢偏 **$\Delta\phi \approx \pm 16 \times 0.10 = \pm 1.60\text{ rad} \ge \pi/2$**，使 $p_{x,A} \le 0$ 反向翻转并在 **`nStates = 1, nMeasurements = 1`** 处中止传播（`Input track has only 1 measurements and cannot be fitted`），两条种子对称失败从而出现 **`lt = 2 → lt = 0`** 崩溃。
  4. 而在 `legacySwappedCovLoc=True`（如 `7fd40df8`，`covLoc0 = 81.0 mm²`）下，$K_\phi$ 下降约 10 倍至 `0.012 rad/mm`，虽避免 $|\phi_{\text{flt}}| \ge \pi/2$ 反向翻转（因此无 `16 mm` 深坑），但 `phi_flt` 仍存在 `-38` 至 `-46 mrad` 的非物理跳变，且错误的 `covLoc1 = 0.0361 mm²` 在大 $|x|$ 处损失效率。
- **最小修复**（[stage08_ckf2_rotated_seed_covphi.patch](file:///eos/home-x/xcheng/FASER/geoRotation_test/patches/stage08_ckf2_rotated_seed_covphi.patch)）：在 `CKF2Config.py` 非 `legacySwappedCovLoc`（`else`）旋转分支将 `sigma_phi_rotated` 设为与 `sigma_theta` 对称的 `2.9e-4 rad`（即 `covPhi = 100 × (2.9e-4)² = 8.41e-6 rad²`），并在 `legacySwappedCovLoc=True` 分支保留 `sigma_phi_legacy = 3.3e-2`；修复后 Event 8 两条径迹均以 `18/18` measurements 完美重建（`chi2Sum = 18.23 / 23.10`），`git apply --check` 在洁净 `218813c2` 源码树通过。

#### 6.2 Tier-2 四路归因全量结果（36 文件对，908,731 匹配事件）

样本构成：`008749` Data 18 文件（`265,606` 事件 = 9 对角 `130,204` + 9 非对角 `135,402`）与 `10004X` MC 18 文件（`643,125` 事件 = 9 对角 `340,463` + 9 非对角 `302,662`），覆盖全部 9 个对角深坑 bin 与 9 个非对角高统计对照 bin。数据来源：[tier2_attribution_summary.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/tier2_attribution/tier2_attribution_summary.json)。

| 样本 / 区域 | `Baseline-Jul8` (%) | `T2-B: 7fd40df8` (%) | `Latest: 218813c2` (%) | **`T2-A: 218813c2 + covPhi=8.41e-6` (%)** | `T2-A − Base` (pp) | `T2-A − Latest` (pp) |
| --- | --- | --- | --- | --- | --- | --- |
| **`008749` 全局 (`265,606`)** | 91.327 | 88.196 (`-3.13`) | 90.913 (`-0.41`) | **93.677** | **`+2.35 ± 0.04`** | **`+2.76`** |
| `008749` 深坑 `[12.5, 18.5]` (`29,085`) | 91.374 | 88.891 (`-2.48`) | **76.658 (`-14.72`)** | **93.846** | **`+2.47 ± 0.11`** | **`+17.19`**（对角 `+24.65`） |
| `008749` 小间距 `\|ΔY0\|<2` | 86.699 | 80.017 (`-6.68`) | 89.939 (`+3.24`) | **91.138** | **`+4.44 ± 0.27`** | `+1.20` |
| `008749` 小间距 `ΔR1<5` (`1,759`) | 88.857 | 85.787 (`-3.07`) | 86.754 (`-2.10`) | **93.121** | **`+4.26 ± 0.65`** | **`+6.37`** |
| `008749` 高原 `ΔR1 [30,70]` (`93,964`) | 91.777 | 87.515 (`-4.26`) | 90.904 (`-0.87`) | **94.586** | **`+2.81 ± 0.06`** | **`+3.68`** |
| `008749` 大水平 `\|ΔX0\|>70` | 90.661 | 90.306 (`-0.35`) | 92.020 (`+1.36`) | **92.651** | **`+1.99 ± 0.08`** | `+0.63` |
| **`10004X` 全局 (`643,125`)** | 92.584 | 88.191 (`-4.39`) | 90.378 (`-2.21`) | **92.945** | **`+0.36 ± 0.01`** | **`+2.57`** |
| `10004X` 深坑 `[12.5, 18.5]` (`70,934`) | 92.669 | 89.063 (`-3.61`) | **77.400 (`-15.27`)** | **92.771** | **`+0.10 ± 0.04`** | **`+15.37`**（对角 `+20.97`） |
| `10004X` 小间距 `\|ΔY0\|<2` | 87.931 | 79.227 (`-8.70`) | 89.359 (`+1.43`) | **90.184** | **`+2.25 ± 0.16`** | `+0.83` |
| `10004X` 小间距 `ΔR1<5` (`5,120`) | 88.105 | 85.078 (`-3.03`) | 84.238 (`-3.87`) | **89.648** | **`+1.54 ± 0.37`** | **`+5.41`** |
| `10004X` 高原 `ΔR1 [30,70]` (`234,372`) | 93.220 | 87.860 (`-5.36`) | 90.818 (`-2.40`) | **93.648** | **`+0.43 ± 0.02`** | **`+2.83`** |
| `10004X` 大水平 `\|ΔX0\|>70` | 91.262 | 89.929 (`-1.33`) | 90.701 (`-0.56`) | **91.494** | **`+0.23 ± 0.03`** | `+0.79` |

**逐阶段损失与物理闭合验证（Tier-2）**：
- **深坑彻底消除**：在深坑 `[12.5, 18.5] mm` 内，`008749` 的 `lt0` 由 `5.821%` 骤降至 **`0.045%`**，`ckf_lost` 由 `12.261%` 骤降至 **`0.289%`**（下降 **42 倍**）；`10004X` 的 `lt0` 由 `9.672%` 回归基线基底 **`1.678%`**，`ckf_lost` 由 `8.165%` 骤降至 **`0.187%`**（下降 **44 倍**）。
- **CKF 损失全面低于未旋转基线**：高原区 `ckf_lost` 在 `008749` 由基线 `2.780%` 降至 **`0.325%`**（降低 **8.6 倍**），在 `10004X` 由基线 `0.484%` 降至 **`0.180%`**（降低 **2.7 倍**）；同时上游 `miss_seg`（`4.71%` Data / `4.25%` MC）在四路间基本不变，确认增益 100% 源自 CKF 寻迹质量提升。
- **MC/Data 高原比值演变**：`ΔR1 ∈ [30, 50] mm` 的 MC/Data 比值从 `Baseline-Jul8` 的 `1.0190` → `T2-B (7fd40df8)` `1.0056` → `Latest (218813c2)` `1.0009` → **`T2-A (covPhi=8.41e-6)` `0.9902`**（即 Data `94.59%` 反超 MC `93.65%`）。这是因为 `T2-A` 几乎消除了 Data 的额外 `ckf_lost`（两者均 `< 0.33%`），余下的 `-0.94 pp` 差值完全对应 `10004X` 跨四路恒定的 `lt0 ≈ 1.58%` 基底（Data 仅 `0.04%`）；若扣除 `lt0`，MC 的 `miss_seg + ckf_lost`（`4.43%`）仍略低于 Data（`5.04%`）。

---

### 7. 主图展示与 22 页统一图集索引

完整 22 页英文发布级汇报 PDF：[stage08_overlay_2track_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/stage08_overlay_2track_report.pdf)（`24.68 MB`，22 页；另含独立 Tier-2 2 页摘编 [stage08_tier2_attribution_figures.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/stage08_tier2_attribution_figures.pdf)）；完整图表元数据索引：[figure_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/figure_index.json) 与 [tier2_figure_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/tier2_figure_index.json)；全部原始分箱计数 CSV 目录：[analysis/bin_counts/](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/analysis/bin_counts/) 与 [tier2_attribution/bin_counts/](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/tier2_attribution/bin_counts/)。

以下嵌入 8 张核心主图（覆盖 Pawan 复现、Baseline vs Latest 对比、动量与残差、高原/大 $|\Delta X|$ 机制分解以及 Tier-2 `covPhi` 因果闭合）：

#### 图 1（Pawan Slide 7 对应）：Latest GeoRotation 效率 vs $\Delta R_1$ 与 `110038` `A'` 信号谱
![fig01_dR1_efficiency](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/png/fig01_dR1_efficiency.png)

#### 图 9：Baseline-Jul8 vs Latest GeoRotation 全量逐事件配对对比（$\Delta R_1$）
![fig09_cmp_dR1_baseline_vs_latest](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/png/fig09_cmp_dR1_baseline_vs_latest.png)

#### 图 11：Baseline-Jul8 vs Latest GeoRotation 水平间距对比（$|\Delta X_0|, |\Delta X_1|$）
![fig11_cmp_abs_dX_baseline_vs_latest](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/png/fig11_cmp_abs_dX_baseline_vs_latest.png)

#### 图 14：动量谱与匹配相对动量残差对比（Baseline-Jul8 vs Latest）
![fig14_cmp_momentum_and_residual](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/png/fig14_cmp_momentum_and_residual.png)

#### 图 17：高原效率逐阶段损失漏斗与依赖关系（回答 Slide 15）
![fig17_plateau_stage_of_loss](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/png/fig17_plateau_stage_of_loss.png)

#### 图 18：大 $|\Delta X|$ 诊断分层（按 $\max(|x_0|)$）与 2D 探测器损失分布（回答 Slides 11–12）
![fig18_large_dX_and_2d_loss_maps](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/png/fig18_large_dX_and_2d_loss_maps.png)

#### 图 21（Tier-2 归因主图 1）：$|\Delta Y_0| \approx 16\text{ mm}$ 效率深坑消除与小 $\Delta R_1$ 效率反转（908,731 事件四路对照）
![fig21_tier2_covphi_crater_attribution](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/png/fig21_tier2_covphi_crater_attribution.png)

#### 图 22（Tier-2 归因主图 2）：四路全相空间总结、水平间距与逐阶段损失归因
![fig22_tier2_covphi_regimes_and_loss](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/png/fig22_tier2_covphi_regimes_and_loss.png)

其余 14 张图（`fig02`–`fig08` 覆盖 $\Delta R_0, |\Delta Y_0|, |\Delta Y_1|, \Delta X_0, \Delta X_1, \Delta Y_1, \Delta Y_0, P_{z0}$ 完整复现；`fig10, fig12, fig13, fig15, fig16, fig19, fig20` 覆盖完整配对对比、5×5 多重数迁移、重加权闭合与 May30 vs Jul8 校验）全部收录于 [stage08_overlay_2track_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/stage08_overlay_2track_report.pdf) 与 [figure_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/figures/figure_index.json)。

---

## 阶段性结论

### 事实（实测与因果闭合）

1. **全量管线与坐标系闭合**：单步零 xAOD `overlay_reco_ntuple.py` + `overlay_acts_to_faser` 完整重构了 `323` overlay 文件（`6,354,483` 事件）与 `10` `110038` 文件（`20,000` 事件），`405/405` HTCondor 作业一次性完成；`6,339,599` 共同事件 `og_Track` 逐条校验零失败，坐标转换 `0 / 12,175,413` 越界，独立 PyROOT 复核零计数偏差。
2. **默认 `Latest-218813c2` 表现（修复 `covPhi` 前）**：
   - 全量 `All3Data` 净升 **`+0.4606 ± 0.0109 pp`**（`92.39% → 92.85%`），`10004X` MC 净降 **`-1.3523 ± 0.0094 pp`**（`93.30% → 91.95%`）。
   - 性能呈显著两极分化：在真实 `110038` `A'` 信号（98% 位于 $|\Delta Y_0| < 2\text{ mm}$）上效率飙升 **`+10.71 pp`**（`52.04% → 62.75%`），在大 $|\Delta X_0|$ 处 Data 提升 `+1.14 pp`；但在 $|\Delta Y_0| \approx 16\text{ mm}$ 处存在深达 `-20.9 pp`（Data）/ `-27.2 pp`（MC）的剧烈效率深坑，且在小 $\Delta R_1 < 5\text{ mm}$ 下滑 `-1.79` / `-3.51 pp`。
3. **Pawan 报告四个遗留问题的实测答案**：
   - **高原 6–7% 损失**：上游 St.2/3 `SegmentFit` 缺段占 `3.5–4.0%`、`all3gb → ckf_lost` 在基线占 `2.47%`（Data）/ `0.49%`（MC），`GhostBusters` 仅占 `0.042%`。
   - **大 $|\Delta X|$ 下降**：核心区 $\max(|x_0|)<33\text{ mm}$ 效率平坦（`95.4–96.6%`）；大 $|\Delta X|$ 强制至少一条径迹进入外缘区域（$|x|>33\text{ mm}$），使 St.2/3 `miss_seg` 增至 `7–8%`。
   - **小 $\Delta R_1 < 5\text{ mm}$ 下降**：基线以 `lt>=3` 径迹分裂（`4.4–6.1%`）为主；旋转版有效压制 `lt>=3`。
   - **`Pz0` 重加权**：重加权改变高原 MC/Data 比值不足 `0.0013`，证实基线 ~1.3% 偏移并非动量谱差异所致，而是 Data 侧更高的 `ckf_lost`。
4. **精确代码缺陷与最小修复因果闭合（Tier-2 908,731 事件）**：
   - `CKF2Config.py:52,67` 在 `fix_v1` 旋转分支中仍沿用未旋转 $1/\sin\theta_F$ 极点奇异遗留的 `sigma_phi = 3.3e-2 rad`（`covPhi = 0.1089 rad²`，$\sigma_\phi = 330\text{ mrad}$，比旋转系中竖直弯曲斜率 $\phi_A \approx t_y$ 的真实误差大 114 倍、方差放大 12,950 倍）；在 `covLoc0` 收紧至 `0.0361 mm²` 后，其导致 St.1 首层竖直搜索窗扩至 $\pm 13$–$16\text{ mm}$ 且 Kalman 增益 $K_\phi \approx 0.10\text{ rad/mm}$ 将 $\phi_{\text{flt}}$ 踢翻过 $\pm\pi/2$，是 $|\Delta Y_0| \approx 16\text{ mm}$ 深坑、`lt=2→lt=0` 崩溃、小 $\Delta R_1$ 下滑及大部分 `ckf_lost` 的唯一根因。
   - 在非 `legacySwappedCovLoc` 旋转分支将 `sigma_phi_rotated` 对齐至 `sigma_theta = 2.9e-4 rad`（`covPhi = 8.41e-6 rad²`，补丁 [stage08_ckf2_rotated_seed_covphi.patch](file:///eos/home-x/xcheng/FASER/geoRotation_test/patches/stage08_ckf2_rotated_seed_covphi.patch)），在 `908,731` 事件上彻底消除深坑（`008749` 深坑 `76.66% → 93.85%`、`10004X` `77.40% → 92.77%`），并使 **Data 与 MC 在全 7 个相空间区域同步显著超越未旋转基线**：`008749` 全局 **`+2.35 ± 0.04 pp`**（`91.33% → 93.68%`）、小 $\Delta R_1 < 5\text{ mm}$ **`+4.26 ± 0.65 pp`**、高原 **`+2.81 ± 0.06 pp`**；`10004X` 全局 **`+0.36 ± 0.01 pp`**、小 $\Delta R_1 < 5\text{ mm}$ **`+1.54 ± 0.37 pp`**、高原 **`+0.43 ± 0.02 pp`**；且 `110038` 200 事件进一步提升至 `59.5%`（`+2.0 ± 1.2 pp` vs 默认，`+11.0 ± 3.1 pp` vs `master_dev`）。

### 假设、范围说明与未测量项

- **2022 历史基线对比（Slide 14）**：EOS 上不存在 Pawan Slide 14 所用的 2022 原始 overlay ntuples，因此本阶段 `Then vs Now` 使用真实冻结的 2026 未旋转基线（`Baseline-Jul8` / `Baseline-May30`）与 `Latest GeoRotation` 的逐事件配对对比，不凭空重构 2022 曲线。
- **`covPhi = 8.41e-6` 的全量与单径迹回归（未测量）**：该最小修复已在 36 shards（`908,731` 双径迹事件）及 `110038` 200 事件完成验证；但**全量 323 shards 重跑、全量 20,000 事件 `110038` 以及单径迹 MC24（`100044 / 100116 / 100117 / 110033 / 110051`）上的动量尾率、角度结构、电荷误判率与大角度低动量径迹效率回归均未测量**，必须作为合入分支前的验收门槛。
- **`10004X` 恒定 `lt0 ≈ 1.60%` 基底成因（待查）**：`10004X` MC 在全部 4 路版本（含未旋转基线与 `T2-A`）均存在约 `1.60%` 的 `lt==0` 基底（显著高于 Data 的 `0.04–0.12%`），该比例不受 `covPhi` 或旋转影响；其具体来源（RDO 底层输入、单径迹过滤或外缘 hit 丢失）尚未逐事件拆解。

---

## 下一步计划

1. **在 `fix/latest-seed-momentum-for-sazhang` 落地 `covPhi = 8.41e-6` 并完成单径迹 + 全量信号回归**：
   - **输入**：洁净 `218813c2` + [stage08_ckf2_rotated_seed_covphi.patch](file:///eos/home-x/xcheng/FASER/geoRotation_test/patches/stage08_ckf2_rotated_seed_covphi.patch)；样本为单径迹 MC24（`100044`, `100116`, `100117`, `110033`, `110051`）及全量 `110038`（20k）与全量 323 overlay shards（均走 HTCondor）。
   - **验收标准**：单径迹动量残差尾率 `f(p_reco/p_truth<0.5)`、角度核心/肩峰、μ± 电荷误判率与单径迹重建效率（尤其低动量/大斜率尾部）无退化；全量 323 overlay 确认全曲线深坑清零、Data/MC 全区域正增益；`110038` 20k 全量效率 $\ge 62.75\%$。
2. **修复 `TwoTrackOverlayNTuplerAlg` 旋转坐标系输出**：
   - 将本阶段已验证的 `actsToFaser` 坐标逆变换合入 `TwoTrackOverlayStudies` 源码（或通过 `geometryRotation` flag 控制），使旋转模式下直接输出标准 FASER 坐标系的 `Track_*` 分支，并与本阶段 `overlay_acts_to_faser` 产物逐值闭合。
3. **上游 `SegmentFit` 外缘缺段与 `10004X` `lt0 ≈ 1.60%` 基底专项排查**：
   - 针对高原区剩余最大的 `miss_seg`（`4.2–4.7%`，集中在 St.2/3 与 $|x|>33\text{ mm}$）及 `10004X` `lt0` 基底，抽取代表事件审计原始 clusters、module 死区/几何间隙与 `SegmentFitAlg` 阈值，定量评估上限可恢复空间。
