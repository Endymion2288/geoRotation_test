# S04：独立 smoothing 与初始 covariance 的全量动量对照

最新复核：[GitHub cuts 下 smoothing OFF/ON 与 writer 逐参数检查](2026-11-12_stage04_github_smoothing_recheck.md)。正式主图采用 GitHub cuts + candidate-isolated smoothing；全量 afterGeo S1 L0 左尾 42.63%→22.20%，S3 基本不变。本报告全量 ON 数据已与旧实现交叉验证。

结果目录：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051`，下文记为 `R`。

**全量结果：恢复原 inflation 明显减少 S2 predicted/filtered 的负尾，但未改善 smoothed 左尾。smoothed residual<−0.8 的比例在九层均增加 0.47–1.13 个百分点；S3 L2 从 10.6772% 升至 11.8046%。beforeGeo/afterGeo 全量比较与 covariance 单因素图均已交付。**

## 阶段目标

在 110051 全部 20 个 RDO 文件、20,000 个源事件上，开启 `IsolateCandidateSmoothing` 重画 beforeGeo/afterGeo。afterGeo 内部只将 seed variance inflation 从 `[100,100,100,10000,100000]` 恢复到 `[100,100,100,100,1000]`，检验动量分布、尤其左侧尾部的变化。

上一轮[鼓包定位](2026-11-12_stage04_smoothing_bump_cause.md)测得 CKF 共享 smoothed slots 覆写会改变 writer 中的逐层状态。此次补丁在候选深拷贝后调用实际 ACTS `GainMatrixSmoother`，修正候选的逐层状态。它位于 CKF2，而原 ROOT residual 公式、cuts 和分箱保持原样。

## 工作内容

| 设置 | beforeGeo | afterGeo | afterGeo_originalVariance |
| --- | --- | --- | --- |
| Calypso commit | `1eb885ae22299fb39bb352bead6005d41b0f953e` | `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` | 同 afterGeo |
| inflation | `[100,100,100,100,1000]` | `[100,100,100,10000,100000]` | `[100,100,100,100,1000]` |
| `covTheta` / rad² | `8.41e-6` | `8.41e-4` | `8.41e-6` |
| `covQOverP` / GeV⁻² | `1e-5` | `1e-3` | `1e-5` |
| 材料 | 原 F frame JSON | AFS rotated JSON | 同 afterGeo |
| isolate smoothing | 开 | 开 | 开 |
| audit / reset-event-maps | 关 / 关 | 关 / 关 | 关 / 关 |

三组均使用 3station forward、mask `[0,1,2,3,4,5]`、FASERNU-04、OFLCOND-FASER-05，原 GitHub 宏的 `chi2Sum/NDF<20 && rMax_smt<90 && NDF>9`、100 bins、`[-1.2,1.2]`、raw Entries 和 P=green/F=red/S=blue 九层布局。主图 truth 仍为原 `t_eQOP`；同 majority barcode 的 production truth 作为原有控制保留。它们均表示 HepMC production momentum。

native 分支的两套源码位于 `geoRotation_calypso_smoothing_before` / `geoRotation_calypso_smoothing_after` 隔离 worktree，每套仅增加相同 smoothing 补丁。`FaserActsGeometry` 与 `FaserActsKalmanFilter` 按 Athena 24.0.41 的既有 compile/link recipes 构建，各 53 个 translation units；未改依赖复用同 release。源码 patch、compile/link recipes、Gaudi metadata、实际加载库和运行 hostname 随 `R` 保存。

前 3 个 covariance 均为 `covLoc0=81 mm²`、`covLoc1=0.0361 mm²`、`covPhi=0.1089 rad²`。恢复原 inflation 只将最后两项 variance 除以 100，相应标准差除以 10；不改 seed 动量、测量误差或 refit 设置。beforeGeo 与 afterGeo 原生角参数属于不同 frame；相同数组不表示已用 Jacobian 构成相同物理 covariance。

每个文件独立处理 1,000 个事件，共 60 个重建任务，保持上一轮的事件分组和 `reset-event-maps=False`。每个 worker 检查实际 `CollectionTree` entries，并在 scratch 内完成 writer truth 关联和原宏 histogram 归约；仅保留 compact counts、统计、日志及逐源事件的三项尾部计数，不保存大规模 states ROOT。

## 详细数据：覆盖、选择与资源

60/60 个全量 Calypso 重建任务成功，三组各覆盖全部 20 文件、20,000 源事件；这是同一批 20,000 个独立输入 entries 的三次处理。每个文件的实际 `CollectionTree` 均为 1,000 entries，处理范围为 `[0,1000)`；各组无重复、遗漏或重叠。Athena event/source entry 保留在紧凑事件计数中，重复的生产 run/event 未被删除。

| 配置 | writer tracks | 通过原宏 cuts | cuts / writer | smoothed raw states | 有限 smoothed residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| beforeGeo | 44,238 | 36,178 | 81.7804% | 623,584 | 617,104 |
| afterGeo 增大 variance | 41,940 | 33,597 | 80.1073% | 579,318 | 573,170 |
| afterGeo 原 variance | 42,434 | 34,829 | 82.0781% | 603,095 | 597,150 |

`cuts / writer` 是 writer 集合内的通过比例。恢复原 variance 后，writer 增加 494 条，cuts 通过数增加 1,232 条（3.6670%）。本次报告总体分布及选择变化；固定相同 hits/candidate 的 covariance 更新效应与候选迁移尚未拆分。逐 seed 失败、refit fallback、SDO purity matching 去重：未另行测量，audit 按原设置关闭。

三组 smoothed 的 invalid reconstruction 均为 **0**；缺失 truth 分别为 **6,480 / 6,148 / 5,945** 个 states。它们仍保留在原宏 raw Entries 与 flow 统计中。下面的尾率使用有限 residual 分母；本次 smoothed 左尾变化不来自原宏将 invalid reconstructed q/p 映射到 −1。

运行软件为 Athena24.0.41、ACTS32.0.2、EL9/GCC13。60 个运行逐一核对实际加载的三个 ACTS 库和五项 seed covariance；全部符合冻结设置。提交前均重新选择 schedd，实际选择为 `bigbird20.cern.ch`；负载广告、队列、UTC 时间、资源与选择依据保存在每个 `*_worker/submission.json`、`schedd_load.txt`、`schedd_queue.txt`。submit 文件位于 `/afs/cern.ch/user/x/xcheng/condor/geoRotation_test/` 对应任务目录。

| 阶段 | cluster ID | 任务数与申请资源 |
| --- | --- | --- |
| before/after 原生构建 | 1573894 / 1573897 | 各1；4 CPUs、10,000 MB、12 GB disk |
| 三组20事件检查 | 1573900 / 1573901 / 1573902 | 各1；1 CPU、3,500 MB |
| 三组全量重建 | 1573904 / 1573905 / 1573906 | 各20；1 CPU、3,500 MB、12 GB disk |
| 小样本图表检查 | 1573907 / 1573914 | 各1；1 CPU、6,000 MB |
| 全量合并与图集 | 1573931 | 1；1 CPU、6,000 MB、12 GB disk |

CPU 重建设置 `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`，未申请或使用 GPU。交互节点实测 `lxplus901.cern.ch`；GPU 查询记录在 `R/interactive_environment.json`。实际 worker hostname、逐库路径和逐文件资源记录见 `R/resource_summary.json`。

| 配置 | 重建 CPU hours | acc.run wall 中位数 / s | acc.run wall 范围 / s | 最大 RSS / MB |
| --- | ---: | ---: | ---: | ---: |
| beforeGeo | 8.697 | 1645.8 | 1159.7–2219.9 | 2338.4 |
| afterGeo 增大 variance | 9.369 | 1826.8 | 1415.0–2108.6 | 2341.2 |
| afterGeo 原 variance | 7.671 | 1507.8 | 1296.3–1953.7 | 2259.8 |

这些任务使用不同 worker，表中时间记录实际资源消耗。源码来源及主 Calypso 树已有本地改动清单见 `R/source_provenance.json`；本次两套 tracking 软件的精确 diff、53 个 compile recipes、link recipes、生成的 Gaudi metadata 与 `software.json` 位于 `R/software_before/`、`R/software_after/`。通用补丁为 [0010_candidate_isolated_smoothing.patch](../patches/0010_candidate_isolated_smoothing.patch)，已在两个实际 worktree 上分别通过反向 `git apply --check`。

## 详细数据：smoothed 左尾

负尾定义为 `p_reco/p_truth<0.5`（residual<−0.5）；靠近 −1 的区域定义为 `p_reco/p_truth<0.2`（residual<−0.8）。`p_truth` 是原 writer 的 SiHit-linked HepMC production momentum。所有表格按 measurement state 统计、保留两 strip sides；95% 区间按源事件聚类并配对，用 ratio delta method 计算，不把同事件内多个 states 当作独立事件。Δ 始终为 **原 variance − 增大 variance**，单位 pp 为百分点。

| 层 | beforeGeo 负尾 / % | 增大 variance / % | 原 variance / % | Δ / pp | Δ 的95%区间 / pp |
| --- | ---: | ---: | ---: | ---: | ---: |
| S1 L0 | 45.0805 | 22.2013 | 22.3193 | +0.1179 | [-0.2509, +0.4868] |
| S1 L1 | 44.7047 | 22.0247 | 22.0794 | +0.0547 | [-0.3158, +0.4253] |
| S1 L2 | 44.1926 | 21.5799 | 21.6978 | +0.1179 | [-0.2465, +0.4823] |
| S2 L0 | 44.2701 | 21.7064 | 21.8238 | +0.1174 | [-0.2409, +0.4757] |
| S2 L1 | 44.2713 | 21.6488 | 21.7918 | +0.1429 | [-0.2224, +0.5083] |
| S2 L2 | 43.7875 | 21.1840 | 21.4457 | +0.2617 | [-0.0994, +0.6229] |
| S3 L0 | 43.6490 | 20.6879 | 21.3768 | +0.6889 | [+0.3444, +1.0334] |
| S3 L1 | 43.5598 | 20.4905 | 21.1816 | +0.6911 | [+0.3382, +1.0441] |
| S3 L2 | 43.2185 | 19.8993 | 20.6371 | +0.7378 | [+0.3891, +1.0865] |

S1/S2 的 smoothed 负尾变化为 +0.05–0.26 pp，逐层区间均包含 0；S3 为 +0.69–0.74 pp。恢复原 variance 没有使这些 smoothed 负尾减少。

| 层 | 有限N：增大 variance | 有限N：原 variance | residual<−0.8：增大 / % | 原 / % | Δ / pp | Δ 的95%区间 / pp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| S1 L0 | 63,879 | 66,279 | 12.6646 | 13.1686 | +0.5040 | [+0.1621, +0.8459] |
| S1 L1 | 63,347 | 65,835 | 12.5184 | 12.9901 | +0.4717 | [+0.1326, +0.8108] |
| S1 L2 | 63,851 | 66,555 | 12.0969 | 12.6256 | +0.5287 | [+0.1974, +0.8601] |
| S2 L0 | 64,474 | 66,936 | 12.2856 | 12.7614 | +0.4759 | [+0.1469, +0.8048] |
| S2 L1 | 64,045 | 66,548 | 12.2320 | 12.7547 | +0.5227 | [+0.1891, +0.8562] |
| S2 L2 | 64,020 | 66,582 | 11.8494 | 12.4433 | +0.5939 | [+0.2649, +0.9228] |
| S3 L0 | 63,530 | 66,591 | 11.3663 | 12.4041 | +1.0378 | [+0.7263, +1.3493] |
| S3 L1 | 63,283 | 66,180 | 11.1973 | 12.2968 | +1.0994 | [+0.7835, +1.4154] |
| S3 L2 | 62,741 | 65,644 | 10.6772 | 11.8046 | +1.1274 | [+0.8170, +1.4377] |

靠近 −1 的区域在九层方向一致地增加；S1 L0 增加 0.5040 pp，S3 L2 增加 1.1274 pp。相同数值 inflation 数组下，beforeGeo 与 afterGeo 的 smoothed 分布差别仍很明显，例如 S1 L0 负尾为 45.0805% 与 22.3193%。这仍是各自原生 frame/parameterization 下的版本比较；相同数组不代表经过 covariance Jacobian 匹配后的相同物理先验。

## 详细数据：predicted/filtered 的变化

下表同样使用 residual<−0.5 与各自有限分母。全部27个 layer/stage的计数和区间见 `R/reductions/tail_summary.json`、`tail_comparisons.csv`。

| 层与 stage | 增大 variance / % | 原 variance / % | Δ / pp | Δ 的95%区间 / pp |
| --- | ---: | ---: | ---: | ---: |
| S1 L0 predicted | 27.8295 | 22.6739 | -5.1556 | [-5.5770, -4.7341] |
| S1 L0 filtered | 27.8245 | 22.6633 | -5.1612 | [-5.5826, -4.7398] |
| S2 L0 predicted | 49.4610 | 32.4225 | -17.0385 | [-17.4303, -16.6467] |
| S2 L0 filtered | 66.3182 | 42.8095 | -23.5087 | [-24.0270, -22.9903] |
| S2 L1 predicted | 65.0043 | 43.0306 | -21.9737 | [-22.4522, -21.4952] |
| S2 L1 filtered | 63.0338 | 43.1794 | -19.8544 | [-20.3280, -19.3809] |
| S3 L2 predicted | 19.9359 | 20.6477 | +0.7118 | [+0.3648, +1.0588] |
| S3 L2 filtered | 19.9168 | 20.6355 | +0.7187 | [+0.3717, +1.0658] |

S2 P/F 对这两项 seed variance 的联合改变很敏感，例如 S2 L0 filtered 减少 23.5087 pp，而相同层 smoothed 仅增加 0.1174 pp。前向阶段的分布变化没有转化为九层 smoothed 低动量尾的改善。本次未分别改变 covTheta 与 covQOverP，二者的独立贡献未测量；final refit 动量分布亦未在本次重新测量。

## 全量图集与数值产物

主图保留原宏 raw Entries、共同分箱、P/F/S颜色及九层布局。左尾补图使用相同原宏 TH1，纵轴为各 bin entries / 全部 raw Entries（包含 NaN/flow），不按显示窗口重新归一化。图中 Raw N 按 legend 的 beforeGeo / afterGeo 增大 variance / afterGeo 原 variance 顺序排列。

**beforeGeo / afterGeo，均开启候选独立 smoothing**

![110051_geometry_overlay](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/figures/110051_geometry_overlay.png)

**afterGeo 两组 covariance 对照**

![110051_variance_overlay](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/figures/110051_variance_overlay.png)

**smoothed 左尾归一化 log 叠图**

![110051_smoothed_left_tail](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/figures/110051_smoothed_left_tail.png)

[完整16页PDF图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/figures/smoothing_covariance_comparison.pdf)；[逐图索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/figures/figure_index.md)。逐图 PNG/PDF、原宏 standalone 九层图、同数值 inflation 的 before/after 叠图及原有 truth 控制均保存在 `R/figures/`。

原始 bin counts、sumw/sumw2、flow 存于 `R/reductions/*_3station.root`；smoothed 三组 JSON bin 表为 `smoothed_bin_counts.json`。`tail_summary.json` 保存有限分母、原始尾部计数、两项阈值与全部配对区间；逐文件 `event_counts.npz` 保留源 entry 和三阶段逐事件计数，可重算区间。原 ROOT writer states 和临时 truth trees 在 worker scratch 完成归约并随作业清理。

## 验证与复现

- 全覆盖：`R/reductions/coverage.json`，三组均20文件、20,000事件、20 shards；每文件输入 entries 及事件序列通过核对。
- 原宏核对：三组×4 truth references×27 TH1，共324个 entries/flow 闭合检查通过；逐事件 counts 的全样本和与冻结 summary 一致。
- 运行核对：60个作业全部成功，实际库路径、五项 seed covariance、isolate=True、audit=False、reset=False 均通过。
- 图表核对：已目视检查全量 geometry/variance 叠图与左尾图；16页汇编与逐图PDF的内容、MediaBox、CropBox一致，记录见 `R/figures/combined_figure_validation.json`。
- 回归：`python3 -m unittest discover -s tests -v`，27/27通过；`git diff --check`通过。日志见 `R/unit_tests.log`。

实际入口：

```bash
python3 scripts/prepare_smoothing_covariance.py \
  --output /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051 \
  --reference /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison
python3 scripts/drive_smoothing_covariance.py --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051
```

构建命令、source/patch/输出库目录及4-CPU申请保存在 `R/before_build_tasks.json`、`after_build_tasks.json`；具体 compile/link recipes 与库副本已冻结。驱动程序使用这些已完成构建，先验收20事件检查，再提交60个全量任务和实际分析任务。全部提交记录为 `R/all_submissions.json`；全量分析实际命令为：

```bash
python3 /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/code/analysis/smoothing_covariance.py --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051
```

该分析命令由 HTCondor cluster 1573931 执行。复现绘图时可使用已合并 ROOT counts；统计重算入口为 `analysis/smoothing_covariance.py` 中 `statistics(root)`，读取20,000事件的紧凑计数。

## 阶段性结论

1. 用户要求的全110051、isolate-smoothing开启、其他名义配置保持的 beforeGeo/afterGeo 对比已完成。
2. 两项 inflation 恢复原值会明显改变动量分布，最明显的是 S2 P/F 负尾及候选/选样数。smoothed 的 residual<−0.5 在 S1/S2 变化小、S3略增；靠近−1的区域在九层均增。因此本次数据不支持“调大的最后两项造成了 smoothed 左尾鼓包、恢复原值即可消除”的解释。
3. 目前保留增大的两项，相比恢复原值具有略小的 smoothed 低动量极端尾；这是一项相对 HepMC production momentum 的分布测量。sensor-local G4动量、固定hits下的单独增益分解及最终refit分辨率未在本次测量。

## 下一步计划

保留本次三个全量 baseline。若继续拆解初始误差的作用，分别只恢复 covTheta 或 covQOverP，并增加固定共同 hits/candidate 的对照；验收同时报告 P/F/S 左尾、候选迁移和相同源事件配对区间。剩余 smoothed 左尾的物理与程序来源另设受控检查，不从一次初始误差调整直接归因。
