# S04：isolated smoothing 后旧筛选对左侧鼓包的逐项作用

**同一份10,000事件数据、GitHub cuts及candidate-isolated smoothing下，仅加入final SDO purity≥0.7，S1 L0左尾率从21.9076%降至10.3974%，左端抬升比从2.4966降至0.8183；S3 L2从20.4697%降至9.9984%，抬升比从2.3370降至0.7535。每truth选一条也有作用，但与purity重叠。Good95在GitHub cuts之后的额外作用较小。上传图的54个P/F/S直方图已精确复现。**

## 阶段目标

用户提供的 `110051_smoothing_only.png` 显示，在完整旧选择内重算 candidate-isolated smoothing 后，左端几乎不再抬升；使用 GitHub cuts 的全量重画图仍有左峰。本次在上传图对应的同一批实际 Calypso states 上，测量每个旧条件的增量作用、独立添加作用及从完整选择中移除后的作用，并直接测量“左端是否翘起”。

结果目录 `Q=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts`。

## 工作内容与固定条件

本次使用上传图的原始 10 文件、10,000 个源事件，每文件 entries `[0,1000)`。不是全量20文件的新重建；上一轮全量20文件、20,000源事件结果仍见[GitHub smoothing 复核](2026-11-12_stage04_github_smoothing_recheck.md)。事件身份使用冻结输入路径和 entry，未用重复的生产 run/event 去重。

数据是 `stage04_momentum/operations/strict/110051_f00...f09_e000000_01000_current/stages.jsonl.gz` 的实际 Calypso selected shared/isolated snapshots，来源为 commit `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` 加保存于 `software_controls_v3/tracking.patch` 的诊断改动。全部源配置见 `Q/sources/f00.json...f09.json`，其中保留 original task、输入路径、软件目录与 runner。该批固定 strict rotated material、`reset-event-maps=True`、候选 isolation 开启；此次只改选样，不更改这些重建设置。没有将 strict/reset 与 AFS/reset-off 的配置差值混入 cut 增量。

所有普通对照固定每条 selected candidate 的同一 majority barcode、production momentum、hits、P/F/S 数值与100个共同bins `[-1.2,1.2]`。同时检查 shared/isolated 的 hit、geometry ID、predicted 和 filtered 参数逐项一致。唯一单列的 `original_figure` 采用原图的 final-majority truth，以验证上传图的全部54个 shared/isolated P/F/S直方图；与 fixed-selected-truth 的切换影响另行核对。

图中 residual 为 `(p_reco-p_production)/p_production`，逐 sensor measurement state，旋转 ACTS frame。诊断图采用有限 residual 分母，finite、invalid 与flows分别保存；不将其冒称原 GitHub writer 的 raw Entries。单位为无量纲，truth 是 HepMC production momentum。

## 精确实现 GitHub cuts 与旧筛选

先用 ON states 冻结 `chi2Sum/NDF<20 && NDF>9 && rMax_smt<90 mm` 通过集合，再对 shared/isolated 比较固定同一集合；其他筛选对照全部保持 isolation ON。原图并未包含这组 GitHub cuts，因此另外保留原图的完整旧集合，并明确标注与 GitHub 集合的交集。

`rMax_smt` 由每个 smoothed bound 局部位置和实际 sensor transform 计算。transform 来自真实 Calypso GeoModel 导出 `stage05_slopes/diagnostics_v1/geometry/geometry.json`，使用相同 FASERNU-04 / OFLCOND-FASER-05。另运行同一 reference 配置的20事件 audit 检查，将 GeoModel→旋转 ACTS 的完整4×4矩阵与实际 audit surface 比较，全部元素误差为0；所有718个measurement states的writer参数检查也继续通过。旧10文件所有 measurement states 的 sensor z 与归档surface z一致，并检查 wafer→geometry ID 的映射稳定。未使用 audit 中基于 cluster hit 位置的 `rmax_mm` 代替 GitHub smoothed 半径。

旧选择按实际 `pilot.process_event` 的顺序复现：

1. final track 必须有有效 momentum、可关联的 majority truth。
2. signal e±：truth `signal=True`、`status=1`、`abs(pdg)=11`。
3. truth fiducial：生产vertex及ST1–3平均hit位置的 F-frame横向半径均 `<100 mm`。
4. H18：truth在ST1–3的SiHit总数 `≥18`。
5. **final** SDO majority fraction `≥0.7`，分母为该final track的实际hits数。selected purity仅作为另列控制。
6. 每个final-majority truth按 majority hit数、总hit数、chi2和candidate ID选一条；改变purity时重新执行排序。
7. 对已选中的final track应用Good95：`ndof>0 && chi2/ndof<25`、`layers≥7`、`ndof+5≥12`、`rmax_final<95 mm`。**先排序，再应用Good95**，不能事后挑一条通过Good95的备选替代失败的最佳候选。

前四个truth/覆盖条件与candidate质量分别处理。逐项累加说明指定顺序下的贡献；独立添加说明只加这一项时的效果；去掉一项说明其余条件保持时的作用。由于候选重排与条件关联，三类差值不要求可加。

## 鼓包与误差定义

- 普通负尾：residual `<−0.5`。
- 靠近−1的极端负尾：residual `<−0.8`。
- 左端抬升比 `R_left`：原100-bin分箱中，左窗 `[-0.984,-0.792)` 的计数除以肩部窗 `[-0.600,-0.408)` 的计数。两窗均为8个bins，宽度均0.192；`R_left>1` 表示左窗计数高于肩部，`<1` 表示该左窗不再高于肩部。这是分布形状诊断，不是拟合的峰显著性。

所有比例和差值使用源事件聚类的delta method给出95%区间，保持同一event内多个tracks/states的相关性。全量逐层区间见 `Q/summary.json`、`Q/paired_cut_effects.csv`。

## 详细数据：smoothing 与选择的交互

| 选择 | 层 | finite states | shared左尾 / % | isolated左尾 / % | shared R_left | isolated R_left |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| GitHub cuts | S1L0 | 32,733 | 41.8385 | 21.9076 | 4.1386 | 2.4966 |
| GitHub cuts | S3L2 | 31,549 | 20.4697 | 20.4697 | 2.3370 | 2.3370 |
| GitHub ∩ 完整旧选择 | S1L0 | 24,067 | 37.1089 | 10.3544 | 3.4917 | 0.8032 |
| GitHub ∩ 完整旧选择 | S3L2 | 23,735 | 9.9937 | 9.9937 | 0.7457 | 0.7457 |
| 上传图：完整旧选择，原truth | S1L0 | 26,065 | 37.4410 | 10.7539 | 3.5137 | 0.7758 |
| 上传图：完整旧选择，原truth | S3L2 | 25,691 | 10.4122 | 10.4122 | 0.7234 | 0.7234 |

同一GitHub通过集合中，S1 L0重算后仍有 `R_left=2.50>1`；同样的重算在完整旧选择内将其降至 `<1`。S3 L2在每个固定cohort内的smoothed值均未因isolation改变，但加入旧选择后左峰明显减少。这两项控制直接区分了共享states污染与剩余选样尾部。

## 详细数据：逐项累加旧条件

下表始终固定isolation ON、selected-majority truth及GitHub cuts，按所列顺序增加旧条件。候选数仅计具有可用selected-majority truth的候选，finite states分母与原宏raw Entries分开。

| 累加条件 | 可用truth候选 | S1 L0左尾 / % | S1 L0 <−0.8 / % | S1 R_left | S3 L2左尾 / % | S3 L2 <−0.8 / % | S3 R_left |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GitHub cuts基线 | 16,760 | 21.9076 | 12.5653 | 2.4966 | 20.4697 | 11.4393 | 2.3370 |
| +final有效且truth可关联 | 16,760 | 21.9076 | 12.5653 | 2.4966 | 20.4697 | 11.4393 | 2.3370 |
| +signal e± | 16,301 | 22.3426 | 12.8869 | 2.5477 | 20.8946 | 11.7432 | 2.3846 |
| +truth fiducial | 16,292 | 22.3412 | 12.8871 | 2.5465 | 20.8922 | 11.7455 | 2.3846 |
| +H18 | 16,158 | 22.2887 | 12.8290 | 2.5323 | 20.8432 | 11.6923 | 2.3701 |
| +final purity≥0.7 | 12,323 | 10.5289 | 4.1074 | 0.8332 | 10.1440 | 3.7562 | 0.7663 |
| +每truth一条 | 12,291 | 10.4020 | 3.9785 | 0.8059 | 10.0807 | 3.6978 | 0.7547 |
| +final Good95 | 12,238 | 10.3544 | 3.9307 | 0.8032 | 9.9937 | 3.6444 | 0.7457 |

加入final purity这一步，S1 L0左尾变化为 **-11.7597 [-12.3576, -11.1618] 个百分点**，S3 L2为 **-10.6992 [-11.2763, -10.1222] 个百分点**。purity后再选一条，S1仅再减少0.1269个百分点；Good95仅再减少0.0476个百分点。

![Isolated smoothing with cumulative legacy selections](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/figures/110051_cumulative_isolated.png)

## 详细数据：独立添加与去掉一项

每次仅在GitHub+isolation基线上添加所列条件，其余旧条件关闭。Δ为“添加后−基线”，单位为百分点。

| 独立添加条件 | 可用truth候选 | S1 L0左尾 / % | S1 Δ / 百分点 | S3 L2左尾 / % | S3 Δ / 百分点 | S1 / S3 R_left |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| signal e± | 16,301 | 22.3426 | +0.4350 | 20.8946 | +0.4249 | 2.548 / 2.385 |
| truth fiducial | 16,624 | 22.0413 | +0.1337 | 20.6032 | +0.1334 | 2.497 / 2.339 |
| H18 | 16,326 | 22.1481 | +0.2406 | 20.6993 | +0.2295 | 2.508 / 2.346 |
| final purity≥0.7 | 12,705 | 10.3974 | -11.5101 | 9.9984 | -10.4714 | 0.818 / 0.754 |
| 每truth选一条 | 14,799 | 15.6249 | -6.2827 | 15.5579 | -4.9119 | 1.546 / 1.566 |
| 完整Good95 | 16,630 | 21.6691 | -0.2384 | 20.2840 | -0.1857 | 2.465 / 2.317 |
| 仅final chi2/ndof<25 | 16,751 | 21.8731 | -0.0345 | 20.4401 | -0.0296 | 2.493 / 2.336 |
| 仅final layers≥7 | 16,647 | 21.7277 | -0.1799 | 20.3220 | -0.1477 | 2.474 / 2.322 |
| 仅final ndof+5≥12 | 16,712 | 21.8260 | -0.0815 | 20.3673 | -0.1025 | 2.493 / 2.328 |
| 仅final rmax<95 mm | 16,758 | 21.8980 | -0.0095 | 20.4597 | -0.0101 | 2.495 / 2.334 |

从GitHub∩完整旧选择出发，每次移除一项，重新执行其余匹配/排序。下表为左尾率的变化和95% source-event配对区间，单位为百分点；参考完整选择为S1 L0=10.3544%、S3 L2=9.9937%。

| 移除条件 | 可用truth候选 | S1 L0 Δ [95%区间] | S3 L2 Δ [95%区间] | S1 / S3 R_left |
| --- | ---: | --- | --- | --- |
| signal e± | 12,378 | -0.0595 [-0.1026, -0.0163] | -0.0674 [-0.1071, -0.0278] | 0.799 / 0.742 |
| truth fiducial | 12,241 | +0.0057 [-0.0091, +0.0205] | -0.0004 [-0.0012, +0.0004] | 0.803 / 0.746 |
| H18 | 12,320 | -0.0086 [-0.0467, +0.0295] | +0.0000 [-0.0390, +0.0390] | 0.801 / 0.745 |
| final purity≥0.7 | 14,163 | +5.4695 [+5.0326, +5.9065] | +5.7019 [+5.2658, +6.1379] | 1.566 / 1.576 |
| 每truth选一条 | 12,263 | +0.0850 [+0.0311, +0.1390] | +0.0522 [+0.0135, +0.0909] | 0.821 / 0.755 |
| 完整Good95 | 12,291 | +0.0476 [+0.0019, +0.0932] | +0.0870 [+0.0325, +0.1415] | 0.806 / 0.755 |
| 仅final chi2/ndof<25 | 12,239 | -0.0009 [-0.0025, +0.0008] | -0.0004 [-0.0012, +0.0004] | 0.803 / 0.746 |
| 仅final layers≥7 | 12,267 | +0.0197 [-0.0144, +0.0539] | +0.0210 [-0.0092, +0.0513] | 0.801 / 0.746 |
| 仅final ndof+5≥12 | 12,241 | +0.0103 [-0.0061, +0.0267] | +0.0029 [-0.0047, +0.0106] | 0.805 / 0.747 |
| 仅final rmax<95 mm | 12,238 | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] | 0.803 / 0.746 |

完整其他条件保持时，移除purity将S1/S3的R_left重新推到 **1.57 / 1.58**，左端再度翘起。移除每truth排序的变化仅0.0850 / 0.0522个百分点，说明其独立添加时的较大作用已大部分被purity覆盖。条件相关且改变排序集合，因此不能把独立添加的purity与unique降幅直接相加。

![Individual selections with isolated smoothing](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/figures/110051_individual_isolated.png)

![Remove one condition with isolated smoothing](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/figures/110051_remove_one_isolated.png)

在本次GitHub已筛集合中，Good95半径条件是冗余的，移除后所有九层尾率变化均为0；其余Good95分量也很小。这里比较的是不同stage的cut：GitHub用selected smoothed rmax和chi2，Good95用最终refit/fallback track；不预设它们逻辑等价。

在**不加GitHub cuts的完整旧流程**中，Good95仍有可见作用：每truth一条后的S1 L0左尾12.3648%→10.7539%，S3 L2为11.5482%→10.4122%。所以“Good95额外作用小”特指GitHub已经筛过的样本。

## 上传图与重画图为何不同

上传图的两列曲线始终已经采用完整旧选择，其中purity和每truth排序已大量排除混合候选。此时原shared slots还会给S1带来额外负尾；isolation移除这一部分后，剩余分布的左窗不再高于肩部。图中“只改smoothing”描述该固定cohort内的比较，并不表示这个cohort只有GitHub cuts。

重画图采用GitHub cuts，保留了更多低purity及重复匹配候选。此次同一批真实states上，**GitHub+isolation本身就重现左端翘起，加入purity后左窗立即降到肩部以下**。因此不需要假设画图漏用了isolation来解释两张图的形状。S1确实有smoothing污染的作用；S3末端的额外区别主要体现在选样中。

本次精确复现上传图的54个shared/isolated P/F/S直方图全部102个bins、finite/invalid及两项tail计数；另外selected-all基线的54个直方图也逐项复现旧summary。完整旧集合13,242条候选中只有1条final-majority与selected-majority不同：切换为固定selected truth后尾率四位小数不变，S1 R_left仅0.7758→0.7765、S3为0.7234→0.7249，不能解释2.5到0.8量级的差别。

本次旧10文件的GitHub+isolation S1左尾21.91%，与此前AFS材料、reset-off、全20文件重画的22.20%处于相近范围；该跨批次差值未在此次实验中拆分，不能把其0.29个百分点全部分配给选样。此处的逐项筛选贡献来自固定旧数据内的受控比较。

![Exact uploaded-figure histogram reproduction](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/figures/110051_uploaded_figure_reproduced.png)

## 验证、执行与产物

实现见 [selection_ablation.py](../analysis/selection_ablation.py)；回归检查见 [test_selection_ablation.py](../tests/test_selection_ablation.py)。五个有针对性的检查覆盖final purity的分母、排序先于Good95、移除purity时重新排序、Good95严格边界与signal/H18分离，全部通过。

几何/writer交叉检查cluster **1573955**、10个归约任务 **1573957**、合并与绘图 **1573959** 均exit0。提交前比较schedd广告和队列，选择兼容的最低负载bigbird20，排除已有EOS兼容问题记录的bigbird21；UTC时间、load、queue、submit路径和资源申请保存在Q各`*_worker/submission.json`及相邻文件。请求1 CPU，geometry检查3,500 MB/3 GB scratch，归约2,500 MB/2 GB，合并4,000 MB/3 GB，线程数=1；实际分配资源见Condor event logs。

归约合计CPU=520.88 s，最大RSS=447,324 KiB；合并绘图wall=40.11 s、CPU=29.53 s、最大RSS=697,640 KiB。实际小样本重建hostname=`b9p02p8097.cern.ch`、acc.run wall=207.79 s，Athena24.0.41；加载的ACTS库与源软件清单均保存于Q/geometry_check，sensor矩阵最大误差=0，旧源数据覆盖144个wafer的z验证最大误差=0。资源汇总见`Q/resource_summary.json`。

所有6张PNG完整性与6页汇报PDF逐页content、page boxes、rotation一致性检查通过，图已目视检查。54个原图直方图和54个selected-all基线精确复现；58个cohorts、两种smoothing模式和九层的smoothed计数与事件充分统计逐项核对。3,132个P/F/S直方图的finite+flow总和逐项一致，unit-weight sumw/sumw2按同一bin count保存。5个新增回归检查和3个既有smoothing统计检查全部通过；Python编译、41个本地文档链接及`git diff --check`通过。

冻结的归约任务为`Q/tasks.json`，合并任务为`Q/merge_tasks.json`，可按仓库`submit_momentum.py`入口以唯一任务名重新提交；分析脚本支持`--phase reduce/merge/plot/tables`。CSV及schema由已冻结的紧凑summary和histogram格式导出，不再读逐state数据。

```bash
python3 scripts/submit_momentum.py \
  --tasks /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/merge_tasks.json \
  --library /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/software_after/lib \
  --name <unique_task_name> --stage S04-selection \
  --memory-mb 4000 --disk-gb 3 --exclude-schedd bigbird21.cern.ch
```

所有worker只读实际Calypso归档，使用批量数组累计histogram及source-event充分统计，不长期保存逐state副本。`histograms.npz` 存共同bins、flow和原始计数；unit weights下同一数组同时为sumw与sumw2。轴/字段顺序见 `histogram_schema.json`，`event_counts.npz` 保存计算配对区间需要的紧凑事件计数。详细表见 [selection_metrics.csv](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/selection_metrics.csv)、[paired_cut_effects.csv](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/paired_cut_effects.csv)。

[6页完整图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/figures/selection_ablation.pdf)、[图表索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/figure_index.json)、[冻结summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_selection_ablation/20261112_legacy_cuts/summary.json) 和验证记录均在Q。英文图包含原图P/F/S九层复现、累加/独立/移除筛选的九层形状、各层leave-one-out热图和带source-event区间的尾率/左端抬升比。

## 阶段性结论

事实：在固定真实Calypso isolated states及GitHub cuts后，final purity是将左窗从高于肩部压到低于肩部的主要旧条件；每truth排序有明显的独立作用，但与purity高度重叠。signal、fiducial、H18并不单独压低这一鼓包，Good95在GitHub已经筛过的集合中只作小幅补充。上传图与当前重画图的剩余左峰差别，已在同一份states中通过添加/移除purity直接重现。

算法层面：isolation修正的是候选共享smoothed slots的污染；它不会自动移除在自己hit路径上仍产生低动量结果的候选。purity改变保留的候选集合，未改变任何候选的P/F/S拟合值。剩余低purity候选中具体的误关联、物理能损和refit来源比例，本次未逐事件分解。

## 下一步计划

正式主图继续保留用户指定的 GitHub cuts + isolated smoothing；完整旧选择、仅purity控制分别作为附图。若后续要降低完整GitHub样本的剩余尾部，需要对被purity控制移除的候选做逐hit归属和最终动量来源诊断，并用相同输入、候选与truth验证算法变化；本次选样结果不代替该重建研究。
