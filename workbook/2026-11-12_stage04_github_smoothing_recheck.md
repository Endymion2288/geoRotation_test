# S04：GitHub cuts 下 candidate-isolated smoothing 的直接复核

旧筛选逐项实测已完成，见[isolated smoothing 后的 selection 消融](2026-11-12_stage04_selection_ablation.md)：同一GitHub+isolation样本，仅加**final purity≥0.7**就使S1 L0左尾21.91%→10.40%，左窗/肩部计数比2.50→0.82；每truth排序的作用与purity重叠，Good95在GitHub已筛集合中仅小幅补充。

**全量 afterGeo 的 S1 L0 smoothed 左尾率从 42.6341% 降到 22.2013%；S3 L2 从 19.7706% 到 19.8993%。此前图中的 smoothing 已生效，但缺少 OFF/ON 直接叠图，未清楚展示减少发生在哪些层。旧图更小的剩余鼓包还包含 truth purity、final Good95 等额外选择的作用。**

## 阶段目标

回答用户追问：其他筛选完全按 GitHub，candidate-isolated smoothing 使用旧图的实际逐候选重算，为什么 afterGeo 仍有左侧鼓包？检查实际 writer 参数归属，并补充全量 smoothing OFF/ON 分布，区分 smoothing、cut membership 和旧图额外选样的作用。

结果目录 `Q=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_github_recheck`。

## 工作内容

全量叠图复用两次已经完成的实际 Calypso 重建，读取原 GitHub 宏产生的 compact TH1，不重做 tracking 或更改 residual。OFF 来自 `stage05_geometry_comparison`；ON 来自 `stage04_smoothing_covariance/20261112_full110051`。每组均为 110051 全部 20 文件，每文件独立处理 entries `[0,1000)`，合计 20,000 个不同源 entries。文件与 entry 定义事件身份，保留重复生产 run/event。

逐文件核对 80 份运行 manifest：OFF/ON 的输入、event range、mask、材料、conditions 和五项 seed covariance 相同，audit/reset-event-maps 均关闭。各文件实际 CollectionTree entries=1,000，依据上轮 worker 保存的输入检查记录核对。完整定位见 `Q/source_validation.json` 和 `Q/full_summary.json`。

beforeGeo 始终采用原 inflation `[100,100,100,100,1000]` 与原 F frame 材料；afterGeo 始终采用增大后的 `[100,100,100,10000,100000]` 与 AFS rotated 材料。此次 OFF/ON 不同时修改 covariance；恢复原 inflation 的单因素结果见[上一轮全量报告](2026-11-12_stage04_smoothing_covariance.md)。

所有全量图保持 GitHub 的 `chi2Sum/NDF<20 && rMax_smt<90 && NDF>9`、实际启用的 residual、100 bins、`[-1.2,1.2]`、raw Entries 和九层布局。P/F/S 仍用 green/red/blue；OFF 实线、ON 虚线。另附 smoothed-only 线性与 log 归一化图，归一化分母为 raw Entries，finite tail 统计另列。truth 为 writer 原生 `t_eQOP`，即 SiHit 关联 HepMC 粒子的 production momentum，不称为 sensor 局部 G4 动量。

## 详细数据：全量 OFF/ON

| 配置 | writer tracks OFF / ON | GitHub cuts 通过 OFF / ON |
| --- | ---: | ---: |
| beforeGeo | 44,238 / 44,238 | 36,158 / 36,178 |
| afterGeo | 41,940 / 41,940 | 33,178 / 33,597 |

cut 表达式完全一致，但重算改变 smoothed 位置及少量 `rMax_smt`，因此通过候选分别增加 20 / 419 条。下表是各自相同 cuts 下的总体分布，不能把其差值解释为全量固定候选的单因素效应。固定候选检查另列。

负尾定义 `(p_reco-p_truth)/p_truth < -0.5`，分母为该层有限 residual，单位为百分比；与图中 raw Entries 分母分别保存。OFF 历史统计未记录精确 `<−0.8` 数，故不据分箱反推该阈值计数。

| 配置 | 层 | finite OFF / ON | 负尾率 OFF → ON / % | 变化 / 百分点 |
| --- | --- | ---: | ---: | ---: |
| beforeGeo | S1 L0 | 68,406 / 68,442 | 38.5551 → 45.0805 | +6.5254 |
| beforeGeo | S1 L1 | 67,822 / 67,852 | 38.2280 → 44.7047 | +6.4766 |
| beforeGeo | S1 L2 | 68,167 / 68,197 | 37.8350 → 44.1926 | +6.3575 |
| beforeGeo | S2 L0 | 69,309 / 69,347 | 42.0479 → 44.2701 | +2.2222 |
| beforeGeo | S2 L1 | 68,953 / 68,986 | 42.1766 → 44.2713 | +2.0947 |
| beforeGeo | S2 L2 | 69,048 / 69,086 | 41.6898 → 43.7875 | +2.0976 |
| beforeGeo | S3 L0 | 68,682 / 68,714 | 43.6199 → 43.6490 | +0.0292 |
| beforeGeo | S3 L1 | 68,368 / 68,398 | 43.5394 → 43.5598 | +0.0204 |
| beforeGeo | S3 L2 | 68,049 / 68,082 | 43.1953 → 43.2185 | +0.0231 |
| afterGeo | S1 L0 | 63,149 / 63,879 | 42.6341 → 22.2013 | -20.4327 |
| afterGeo | S1 L1 | 62,621 / 63,347 | 32.8883 → 22.0247 | -10.8636 |
| afterGeo | S1 L2 | 63,109 / 63,851 | 31.6215 → 21.5799 | -10.0416 |
| afterGeo | S2 L0 | 63,713 / 64,474 | 21.9704 → 21.7064 | -0.2640 |
| afterGeo | S2 L1 | 63,283 / 64,045 | 21.4323 → 21.6488 | +0.2165 |
| afterGeo | S2 L2 | 63,248 / 64,020 | 21.0378 → 21.1840 | +0.1462 |
| afterGeo | S3 L0 | 62,762 / 63,530 | 20.5507 → 20.6879 | +0.1372 |
| afterGeo | S3 L1 | 62,503 / 63,283 | 20.3350 → 20.4905 | +0.1555 |
| afterGeo | S3 L2 | 61,981 / 62,741 | 19.7706 → 19.8993 | +0.1287 |

![afterGeo GitHub cuts: source-style P/F/S OFF/ON](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_github_recheck/figures/afterGeo_github_smoothing_off_on.png)

![afterGeo smoothed OFF/ON, raw-Entries normalization](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_github_recheck/figures/afterGeo_smoothed_off_on_linear.png)

## 详细数据：旧实现与 writer 逐参数验证

重新运行第一文件 entries `[0,20)` 的实际 Calypso。设置与上一轮 20-event ON smoke 相同，仅换用旧图冻结的 audit-enabled `software_controls_v3/lib/libFaserActsKalmanFilter.so`，开启 selected-state 观察；geometry 两个库及 Python 配置继续使用上一轮冻结的 afterGeo 软件。实际加载路径、commit/patch 来源、五项 covariance、hostname 均保存于 `Q/software_audit/software.json` 与 `Q/audit_check/run.json`。reset-event-maps 和接口修正仍关闭。

旧实现是在深拷贝 selected candidate 后调用 ACTS `GainMatrixSmoother`；audit 分别记录重算前 `selected_shared_state` 与重算后 `selected_state`。与上一轮只加 smoothing 最小补丁的 20-event smoke 对照，四种 truth 控制共 **108 个直方图**的全部 102 个 bins（含 flow）及 raw Entries 逐项完全相等，见 `Q/audit_check/audit_vs_minimal.json`。

writer 归属直接检查得到：

- 20 源事件、42 条 writer tracks、718 个 measurement states；30 条通过原 GitHub cuts。
- 718×3×6=**12,924 个 writer P/F/S 六维参数**，全部等于旧实现重算后快照转换为 float32 的数值。最大相对存储舍入误差 `5.91173e-8`；writer 没有读回重算前的 shared smoothed 值。
- 重算前后 hit、geometry ID、predicted 和 filtered 参数逐项相同；718 个 states 中 156 个 smoothed q/p 改变。
- 从 audit 保存的真实 sensor surface transform 和 smoothed 局部位置重建 writer 的 `rMax_smt`，与 writer 最大差为 `1.42109e-14 mm`。没有用 audit 中基于 hit 位置的 `rmax_mm` 替代 GitHub 的 cut 字段。

固定这 30 条 ON writer 通过候选、同一 hits 和同一 writer truth 后，S1 L0 的60个有限 states 左尾 **20/60→10/60**，S3 L2 的59个有限 states为 **10/59→10/59**。这是精确 GitHub cohort 的小样本单因素验证；全量尾率仍以上节20文件结果为准。检查数据见 `Q/audit_check/writer_smoothing_check.json`。

![Same candidates: old shared and isolated snapshots](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_github_recheck/figures/afterGeo_same_candidates_audit_check.png)


## 旧图为何更干净

[鼓包定位报告](2026-11-12_stage04_smoothing_bump_cause.md)中较干净的旧图右列使用前 10 文件、10,000 源事件、strict material、reset-event-maps 开启，并额外要求 signal e±、truth fiducial+H18、SDO majority fraction≥0.7、final unique matching 和 final Good95。本次按用户要求保留 GitHub cuts，不加入这些额外 truth 选择。

旧报告已在固定同一批候选/hits/truth 后测得：S1 L0 左尾由 37.4410% 降至 10.7539%，S3 L2 的 smoothed q/p 则逐 state 完全相同。固定 shared states 的选择控制中，仅加入 purity≥0.7 就将 S3 L2 左尾从 23.0126% 降至 11.5911%，靠近−1的 `<−0.8` 比例从 13.0142% 降至 4.2873%。因此，不能期待只开启 isolation 就在 GitHub 选样中重现旧图全部净化程度。

## 验证与复现

实现见 [smoothing_recheck.py](../analysis/smoothing_recheck.py)。`reference_tracking.py --audit-states` 是显式 opt-in 观察入口，普通全量入口的 audit 默认仍关闭。audit 原始 states、writer ROOT 和临时宏输出全部在 worker scratch，最终只保留 compact histogram、逐层计数与一致性检查。

全部工作通过 HTCondor 执行，提交前重新查询 schedd 广告和队列，选择兼容的最低负载 `bigbird20.cern.ch`，排除已有 EOS 兼容问题记录的 bigbird21。audit 检查 cluster **1573942**、全量导图 **1573943**、图集合并 **1573945**、诊断图分母标注修订后最终导图 **1573949** 均 exit 0。请求均为1 CPU，audit 3,500 MB、3 GB scratch，绘图2,500 MB、2 GB scratch，BLAS/OpenMP线程=1；实际分配资源见各 Condor event log。负载、队列、UTC时间、submit路径和软件目录保存于 `Q/*_worker/submission.json` 及相邻记录。

audit 重建 hostname=`b9g03p4168.cern.ch`，`acc.run` wall=207.38 s，含初始化的实测wall=219.07 s，CPU user+system=121.91 s，最大RSS=2,036,560 KiB。全量导图wall=35.85 s、CPU=25.59 s、最大RSS=661,764 KiB。运行 Athena24.0.41，实际加载三个ACTS库均已检查；交互节点与GPU查询另存 `Q/source_validation.json`，CPU任务未使用GPU。

7张PNG完整性、逐图PDF与汇报PDF全部page contents、page boxes、rotation一致性检查通过；主图和诊断图已目视检查。72组全量 smoothed 计数数组的100个共同bins、2个flow与raw Entries总和一致，保存sumw2；诊断图明确采用finite residual分母。`python3 -m unittest discover -s tests -v` 当次27项通过，日志见 `Q/unittest.log`；三个修改的Python文件 `py_compile`、34个本地文档链接与 `git diff --check` 通过。

复现全量导图命令保存在 `Q/plot_tasks.json`，逐参数检查完整配置在 `Q/audit_tasks.json`，最终图集合并在 `Q/finish_v2_tasks.json`。例如通过仓库批量入口提交最终任务：

```bash
python3 scripts/submit_momentum.py \
  --tasks /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_github_recheck/finish_v2_tasks.json \
  --library /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/software_after/lib \
  --name <unique_task_name> --stage S04-smoothing-recheck \
  --memory-mb 2500 --disk-gb 2 --exclude-schedd bigbird21.cern.ch
```

全量 counts 和归一化图的逐 bin `sumw/sumw2`、underflow/overflow、raw Entries 见 `Q/bin_counts.json`；P/F/S 原始 counts 见叠图 ROOT 和源 reductions；逐候选配对图 counts 见 `Q/paired_bin_counts.json`。完整 [7页图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_github_recheck/figures/smoothing_recheck.pdf)、[图表索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_github_recheck/figure_index.json)、[验证记录](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_github_recheck/validation.json) 与源 provenance 随结果保存。beforeGeo/afterGeo 开启 smoothing 的正式主图继续使用[全量 covariance 报告](2026-11-12_stage04_smoothing_covariance.md)中的同一份 ON 数据。

## 阶段性结论

全量事实：candidate isolation 显著减少 afterGeo 的 S1 负尾；S2/S3 的剩余左尾基本保留。旧报告的层依赖结论与本次 GitHub cuts 下结果一致。共享 slots 的程序问题解释早期层的一部分额外鼓包；剩余鼓包的主要物理/候选混合来源不能仅凭此次 OFF/ON 整体图定量分解。

## 下一步计划

继续保留 GitHub cuts + isolated smoothing 作为用户指定的主图定义。若进一步追查剩余左尾，固定通过集合、majority barcode 和 hits，单列 purity 分组及 sensor G4 truth 控制，验收各层尾率变化；保留完整 GitHub 样本分布作为基准。
