# S04：native truth、isolated smoothing 与 final purity 的几何前后对比

结果目录：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_native_purity/20261112_full110051`，下文记为 `R`。

**全量结果：两组均已开启 isolated smoothing 并应用严格 final purity > 0.7，只使用 native truth。afterGeo 的 S1 L0 smoothed 负尾从仅 isolation 时的 22.2013% 降至 10.5236%；本次 beforeGeo 为 36.5994%。**

## 阶段目标

在 110051 全部 20 个 RDO 文件、20,000 个独立源事件上，比较 beforeGeo 与 afterGeo。两组均开启 candidate-isolated smoothing，在 GitHub 原宏 cuts 上增加严格的 `final purity > 0.7`；truth 只使用原 writer 的 native `t_eQOP`。

参考为 [geometryRotation.pdf](../docs/geometryRotation.pdf) 第17页及 [GitHub test_pResolution 原宏](https://github.com/zsm2024/FASER_tracking/tree/main/test_pResolution)，冻结源 commit 为 `33ca2281b961935456dfd4dc9fc06235e24b551e`。

## 工作内容与冻结设置

两组沿用[上一轮全量 isolation 对照](2026-11-12_stage04_smoothing_covariance.md)的原生 Calypso 软件、材料、事件分组和 covariance：

| 设置 | beforeGeo | afterGeo |
| --- | --- | --- |
| Calypso commit | `1eb885ae22299fb39bb352bead6005d41b0f953e` | `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` |
| 冻结库 | 上轮 `software_before/lib` | 上轮 `software_after/lib` |
| initial variance inflation | `[100,100,100,100,1000]` | `[100,100,100,10000,100000]` |
| `covTheta` / rad² | `8.41e-6` | `8.41e-4` |
| `covQOverP` / GeV⁻² | `1e-5` | `1e-3` |
| material JSON | `stage02_closure/material-current-source.json` | `/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json` |
| candidate-isolated smoothing | 开 | 开 |
| audit / reset-event-maps | 关 / 关 | 关 / 关 |

三项共同 covariance 为 `covLoc0=81 mm²`、`covLoc1=0.0361 mm²`、`covPhi=0.1089 rad²`。这是两套原生配置的版本比较；两组的 covariance 尚未通过坐标 Jacobian 匹配为相同物理先验。

共同设置为 3station forward、mask `[0,1,2,3,4,5]`、FASERNU-04、OFLCOND-FASER-05。每个文件完整处理 entries `[0,1000)`，每文件一个作业，保留 `reset-event-maps=False` 下的原事件分组。输入目录清单与文件大小重新核对，worker 再读取各文件实际 `CollectionTree` entries。

GitHub 原宏的 selection 保留为：

```text
chi2Sum/NDF < 20 && rMax_smt < 90 && NDF > 9
```

本次只增加 `&&final_purity>0.7`。保留原 residual 公式、100 bins、`[-1.2,1.2]`、raw Entries、P=green/F=red/S=blue 和九层布局。原 residual 为 `p_reco/p_native - 1`；writer 的 reconstructed `q/p` 单位为 GeV⁻¹，native truth `q/p` 为 MeV⁻¹，原宏进行相应的 1,000 单位换算。NaN、flow 和除零保护沿用原宏。

isolation 使用 [patch0010](../patches/0010_candidate_isolated_smoothing.patch)：在候选深拷贝后调用实际 ACTS `GainMatrixSmoother`。两套实际 worktree 分别为 `geoRotation_calypso_smoothing_before` 和 `geoRotation_calypso_smoothing_after`，复用上轮构建，没有重新编译 tracking。精确源码 diff、构建与共享依赖记录位于上轮两套 `software.json`。

### final purity 的实现

观察器读取实际 `CKFTrackCollectionWithoutIFT` 的每条最终输出轨迹：refit 成功时为 refit 轨迹，失败时为原轨迹。收集 `measurementsOnTrack()` 中的 `FaserSCT_ClusterOnTrack::prepRawData()`，调用**实际加载的 Calypso 库**所导出的 `identifyContributingParticles(SCT_SDO_Map, clusters, counts)`。源码为 `FaserActsKalmanFilter/src/TrackClassification.cxx:77`。

```text
final_purity = max_barcode(number of associated SCT clusters) / all final-track SCT clusters
```

同一 barcode 在同一 cluster 内只计一次，一个 cluster 可以关联多个 barcode；不做 deposit 能量加权。没有 SDO 关联的 cluster 仍在分母。严格 `>0.7` 的整数交叉检查为 `10*majority_hits > 7*clusters`，恰好 `7/10` 被排除。

实现见 [reference_truth_observer.py](../scripts/reference_truth_observer.py) 与 [native_purity.py](../analysis/native_purity.py)。观察器只读 final tracks 和 SDO，不更改重建。final track 与 selected CKF writer 按“输入文件 + source entry + track ordinal”对应，并检查 writer event_nr 和 track_nr；重复的生产 run/event 不用于去重。purity 仅作 track gate，图中 native `t_eQOP` 原值保留，不替换为 final majority 的参考动量。计数同时记录 native writer / final majority 不同的通过轨迹。

原宏差异保存在 `R/reference/purity_gate.diff`；原 GitHub 来源与既有适配记录保留于 `R/reference/`。最终 histogram 由冻结宏直接生成，Python 另行统计有限 residual 与左尾，用 raw Entries/flow closure 检查两者分母关系。worker 在 scratch 内完成归约，最终仅保留 compact ROOT histograms、逐事件尾部计数、purity 分子/分母频数与日志。单独导图宏为每层创建独立 legend 并显示三条曲线的 raw Entries，差异保存在 `R/reference/plot_labels.diff`；worker 的计数宏保持冻结原值。

## 详细数据：小样本验证

两组均重跑第一个冻结文件的 entries `[0,20)`。额外生成未加 purity 的 native 控制图，与上轮同配置 smoke 比较：**54 个直方图，5,508 个 bins（含 flow）和 raw Entries 全部逐项相等**。见 `R/smoke_validation.json`。这验证了新增观察器不改变原生 writer 分布。

| 配置 | writer tracks | GitHub cuts | GitHub cuts + purity > 0.7 |
| --- | ---: | ---: | ---: |
| beforeGeo，20 events | 44 | 30 | 25 |
| afterGeo，20 events | 42 | 30 | 23 |

阈值回归检查覆盖 `7/10` 拒绝、`13/18` 通过、`12/18` 拒绝和缺失关联留在分母。`python3 -m unittest discover -s tests -v`：38 tests passed，日志 `R/unittest.log`。五页小样本图集已导出并检查主图与左尾图布局；其标题明确标为 1 文件、20 events。

首轮 smoke clusters `1573970/1573971` 因 ROOT interpreter 未配置 HepMC3 而失败。按 Athena24.0.41 实际 compile recipe 设置 `HEPMC3=1` 与对应 include 后，实际冻结库的 helper 编译成功；重跑 `1573974/1573975` 两组通过。失败日志保存在 `R/smoke_failed_v1/`，最终验证基于重跑结果。

## 详细数据：全量覆盖与轨迹选择

40/40 个全量重建任务全部成功。每组覆盖20文件、20,000个独立输入 entries；这是同一批输入的两次处理。逐文件 `CollectionTree=1000`、范围 `[0,1000)` 和文件/entry 唯一键均已检查，无遗漏或重叠。每个 job 的实际加载三个 ACTS 库、五项 seed covariance、isolation、audit/reset 设置均与冻结任务相符。

| 配置 | writer tracks | GitHub cuts | 加 purity > 0.7 | smoothed raw states | 有限 smoothed residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| beforeGeo | 44,238 | 36,178 | 26,530 | 462,514 | 461,404 |
| afterGeo | 41,940 | 33,597 | 25,306 | 441,209 | 439,921 |

两组 writer tracks 与未加 purity 的上一轮全量结果分别完全一致，GitHub cuts 通过数也分别完全一致。新增 purity gate 后，保留的是各自 GitHub cohort 的子集。

两组 smoothed invalid reconstruction 均为0，缺失 native truth 分别为1,110和1,288个 states；它们仍按原宏保留在 raw Entries 和 flow 中。最终 gate 通过轨迹中，native writer majority 与 final majority 不同的分别为63和10条，native truth 均保持原值。全部 writer tracks 中 final purity 恰为0.7的分别为9和14条，严格 gate 不接受它们。

所有54个 layer/stage直方图的 raw Entries、flow 与独立有限 residual 计数闭合，见 `R/reductions/summary.json` 的 `histogram_flow_checks`。`bin_counts.json` 和 `histogram_schema.json` 保存原始 counts、sumw、sumw2、共同边界及 flow 索引；本次均为单位权重。

### smoothed 左尾

负尾定义为 residual < −0.5，近−1区域定义为 residual < −0.8，分母均为该层有限 reconstructed/native truth residual。Δ为 afterGeo − beforeGeo，pp为百分点。95%区间以同一20,000个源事件配对，用 source-event cluster ratio delta method 计算；事件内多个 tracks/states 不当作独立事件。两组的通过候选可以不同。

| 层 | beforeGeo：<−0.5 / % | afterGeo：<−0.5 / % | Δ / pp | Δ的95%区间 / pp | beforeGeo：<−0.8 / % | afterGeo：<−0.8 / % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| S1 L0 | 36.5994 | 10.5236 | -26.0758 | [-26.6688, -25.4828] | 12.1204 | 4.0172 |
| S1 L1 | 36.5342 | 10.5031 | -26.0311 | [-26.6239, -25.4384] | 12.0263 | 4.0093 |
| S1 L2 | 36.1458 | 10.4108 | -25.7350 | [-26.3242, -25.1459] | 11.8875 | 3.9247 |
| S2 L0 | 36.0680 | 10.3535 | -25.7146 | [-26.2994, -25.1297] | 11.8887 | 3.9093 |
| S2 L1 | 36.2620 | 10.3575 | -25.9045 | [-26.4939, -25.3151] | 11.9062 | 3.9127 |
| S2 L2 | 35.9678 | 10.3185 | -25.6492 | [-26.2373, -25.0611] | 11.8267 | 3.8458 |
| S3 L0 | 36.0087 | 10.1526 | -25.8561 | [-26.4411, -25.2711] | 11.7447 | 3.7256 |
| S3 L1 | 36.0515 | 10.1942 | -25.8573 | [-26.4464, -25.2682] | 11.7333 | 3.7299 |
| S3 L2 | 35.8108 | 10.1199 | -25.6909 | [-26.2778, -25.1040] | 11.6290 | 3.6351 |

同一组中，相对上一轮“GitHub cuts + isolation”的变化直接量化了新增 purity gate 的总体选择效应：

| 配置/层 | 原 <−0.5 / % | 加 purity 后 / % | Δ / pp | Δ的95%区间 / pp |
| --- | ---: | ---: | ---: | ---: |
| beforeGeo S1 L0 | 45.0805 | 36.5994 | -8.4811 | [-8.8240, -8.1382] |
| beforeGeo S3 L2 | 43.2185 | 35.8108 | -7.4077 | [-7.7476, -7.0677] |
| afterGeo S1 L0 | 22.2013 | 10.5236 | -11.6777 | [-12.0931, -11.2624] |
| afterGeo S3 L2 | 19.8993 | 10.1199 | -9.7793 | [-10.1707, -9.3880] |

afterGeo S1 L0 的 residual < −0.8 同时由12.6646%降至4.0172%，减少8.6474 pp，95%区间为[−9.0078, −8.2870] pp；S3 L2由10.6772%降至3.6351%。这说明仅修正 smoothing 后仍包含的低 purity 轨迹，对靠近−1的尾部贡献很大。purity 是对轨迹样本的 truth 选择，不是新的动量更新；本结果不测量固定候选/hits 下的动量修复量。

### 正式图集

主图保留原宏 raw Entries、共同分箱和可比的轴范围；原生 beforeGeo 与旋转 afterGeo 的 P/F/S 叠图中，实线为 beforeGeo、虚线为 afterGeo：

![Full 110051 native geometry comparison](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_native_purity/20261112_full110051/figures/110051_geometry_overlay.png)

仅 smoothed 的左尾形状另用各自 raw Entries 归一化，方便观察样本数不同下的分布：

![Full 110051 smoothed native left tail](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_native_purity/20261112_full110051/figures/110051_smoothed_left_tail.png)

[完整五页PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_native_purity/20261112_full110051/figures/native_purity_comparison.pdf) · [beforeGeo原宏九层图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_native_purity/20261112_full110051/figures/110051_beforeGeo_native_nine_layers.png) · [afterGeo原宏九层图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_native_purity/20261112_full110051/figures/110051_afterGeo_native_nine_layers.png) · [左尾log图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_native_purity/20261112_full110051/figures/110051_smoothed_left_tail_log.png) · [图集索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_native_purity/20261112_full110051/figures/figure_index.md)。每图同时保存PNG与单页PDF，全部只使用native truth。

最终导图 exit 0。PDF中两组共54个逐层 N 标注逐项等于冻结 histogram 的 raw Entries；五个单页与合并PDF的content、media/crop boxes及rotation一致，验证见 `R/figures/figure_validation.json`。五张PNG及合并PDF的五页渲染已目视检查，PDF预览保存在 `R/pdf_review_contact.png`。本仓库源码与实际冻结执行副本的6组比较均一致，见 `R/code_validation.json`。

## 执行与复现

入口为 [prepare_native_purity.py](../scripts/prepare_native_purity.py)、[run_reference_task.py](../scripts/run_reference_task.py) 和 [native_purity_figures.py](../analysis/native_purity_figures.py)。准备命令：

```bash
python3 scripts/prepare_native_purity.py \
  --reference /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051 \
  --output /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_native_purity/20261112_full110051
```

所有重建通过 HTCondor。提交前重新查询 schedd 广告与队列，选择兼容的最低负载 `bigbird20.cern.ch`，排除既有 EOS 兼容问题的 bigbird21；UTC 时间、候选负载、选择依据、任务、库和资源记录保存在各 `*_worker/submission.json` 与相邻 load/queue 文件。beforeGeo/afterGeo 全量 clusters 分别为 `1573978/1573979`，各20个作业，申请每作业1 CPU、3,500 MB、12 GB scratch、workday。线程限制为1。交互节点为 `lxplus901.cern.ch`，检测到 Tesla T4；本次 CPU 重建没有使用 GPU。

冻结配置、逐文件实际 entries、宏 diff 与 smoke 验证均保存在 `R`；`configuration_validation.json` 核对40个任务相对上轮只更换观察器/宏入口、输出与新增 purity 开关。

全量重建实际资源，CPU hours包括 user+system，wall为 `acc.run`：

| 配置 | CPU hours | wall中位数 / s | wall范围 / s | 最大RSS / MB |
| --- | ---: | ---: | ---: | ---: |
| beforeGeo | 7.824 | 1507.2 | 1186.6–1615.7 | 2329.2 |
| afterGeo | 8.028 | 1564.7 | 1054.1–1802.8 | 2435.1 |

逐文件 hostname、CPU、wall与内存保存在 `R/resource_summary.json`。全量合并/导图 cluster 为 `1573987`，补齐逐层legend的raw Entries后最终导图为 `1573994`，命令分别位于 `R/merge_tasks.json` 与 `R/plot_final_tasks.json`。导图调试中 `1573990` 缺少宏头文件而失败，补齐后 `1573991` 完成；目视核对发现原宏共用legend会重复第一层N，最终版本改为逐层legend并增加PDF文本与冻结counts逐项核对。

```bash
python3 R/code/analysis/native_purity_figures.py --root R --phase all
```

这里的 `R` 需替换为本文开头的绝对路径，使用相同 Athena24.0.41 环境。

## 阶段性结论

本次全量 beforeGeo/afterGeo 图已按指定的 GitHub cuts、candidate-isolated smoothing、final purity严格>0.7和native truth生成。新增purity选择明显减少afterGeo左尾；S1 L0的有限负尾比例为10.5236%，beforeGeo为36.5994%。这些数值包括各自原生covariance与通过候选集合的差别，属于指定配置的整体对照。

## 下一步计划

如继续研究几何差异的物理来源，需固定共同候选/hits，并在不同frame间匹配seed covariance；验收应同时报告固定cohort的分布变化与候选迁移。
