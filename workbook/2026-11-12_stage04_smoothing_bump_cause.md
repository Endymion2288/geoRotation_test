# S04复核：smoothed动量在−1附近鼓包的来源

最新复核：[GitHub cuts 下 smoothing OFF/ON 与 writer 逐参数检查](2026-11-12_stage04_github_smoothing_recheck.md)。正式主图采用 GitHub cuts + candidate-isolated smoothing；全量 afterGeo S1 L0 左尾 42.63%→22.20%，S3 基本不变。本文历史关闭开关的描述仅指对应旧批次。

日期：2026-11-12。结果目录为`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/bump_diagnosis_20261112`，下文记为`B`。

**结论：差别由候选smoothing处理和选样共同造成，不是单一cut。S1的突出鼓包主要来自原生shared smoothed slots；旧图开启候选隔离后的smoothing重算，使其显著减少。S3末端的差别主要来自SDO purity及final选样，该层的状态完全不受这次重算影响。**

## 阶段目标

解释用户指定的两套图为何在smoothed residual≈−1处差别很大：

- 新图：`stage05_geometry_comparison/figures/*nine_layers.png`，实际Calypso原生CKF writer、GitHub cuts。
- 旧图：`stage04_momentum/reference_figures/figures/110051_source_binning_nine_layers.png`，其右列使用strict材料、S04 final Good95祖先及candidate-isolated smoothing。

验收为固定相同输入、几何和材料，分别测量smoothing单因素与选择条件的影响；给出与原图对应的九层分布，并逐bin复现旧图验证归约。

## 工作内容与受控比较

本次使用旧图原有的前10个110051 RDO文件、10,000源事件，以及原任务保存的真实Calypso `selected_shared_state_trace`和`selected_state_trace`。二者分别在同一候选复制后、隔离smoothing重算前后记录。未运行独立tracking，也没有重新生成拟合状态。

几何固定geoRotation；材料固定`stage02_closure/material-current-rotated.json`。新图选`afterGeo_strict`相同前10文件，避免混入AFS材料、beforeGeo主分支或额外10文件。旧图左列是材料坐标修正前，不与beforeGeo几何对齐。

| 项目 | 新图原生口径 | 旧图口径 |
| --- | --- | --- |
| 画的对象 | native CKF selected tracks的measurement states | matched final Good95的selected CKF祖先states |
| 主要选择 | `chi2Sum/NDF<20 && rMax_smt<90 && NDF>9` | signal truth、fiducial+H18、SDO匹配≥0.7、每truth选一条final、final Good95 |
| smoothed slots | 保留原生finder输出 | 候选复制后独立重算 |
| truth | 本次统一为同一匹配粒子的production p | production p |
| 分箱 | 100 bins、[-1.2,1.2]、raw Entries | 相同 |

旧图的truth fiducial要求vertex及station1/2/3平均truth位置的横向半径均小于100 mm；H18要求下游三站SiHit计数和≥18。Good95作用于final track，要求χ²/ndof<25、至少7层、ndof+5≥12及rmax<95 mm。它们不同于原宏直接作用于CKF selected track的三项cuts。旧运行还开启`--reset-event-maps`；本次smoothing单因素只比较同一次运行保存的前后快照，不把两套运行间全部差值归为cut效应。

两份原运行的全部selected candidates均为20,928条；新图通过GitHub cuts为16,557条，旧图最终选择为13,242条。每个state的前后配对要求同一seed/candidate、相同hit及measurement flag；全部predicted、filtered六维向量逐值一致，变化只发生于smoothed。

## 详细数据：固定旧图tracks和cuts，只替换smoothed states

定义负尾为`p_reco/p_prod<0.5`，即residual<−0.5；“靠近−1区域”固定为`p_reco/p_prod<0.2`，即residual<−0.8。分母为各层有限smoothed measurement states，两strip sides均保留。

| 层 | 相同states数 | 负尾：原生shared→isolated / % | 靠近−1：原生shared→isolated / % |
| --- | ---: | ---: | ---: |
| S1 L0 | 26,065 | 37.4410→10.7539 | 24.3814→3.9517 |
| S1 L1 | 25,783 | 28.4451→10.7629 | 15.4986→4.0026 |
| S1 L2 | 25,885 | 27.0156→10.7398 | 14.2438→3.9251 |
| S2 L0 | 25,793 | 17.0783→10.6657 | 8.3976→3.9197 |
| S2 L1 | 25,736 | 16.1758→10.6777 | 7.8917→3.9478 |
| S2 L2 | 25,812 | 15.9732→10.6617 | 7.6399→3.8238 |
| S3 L0 | 25,706 | 10.5306→10.5306 | 3.7890→3.7929 |
| S3 L1 | 25,749 | 10.5713→10.5713 | 3.7827→3.7904 |
| S3 L2 | 25,691 | 10.4122→10.4122 | 3.6667→3.6667 |

S1 L0有14,085/26,065个states的smoothed q/p改变；其中7,441个从负尾恢复、485个新进入负尾，净减少6,956个。靠近−1区域从6,355降到1,030，净减少5,325个。S3 L2的25,691个smoothed q/p则逐个完全相同，排除该层鼓包差别由候选隔离重算造成。

九层合计44,681个smoothed q/p改变；原shared快照有1,933个measurement states的smoother递推恒等式误差超过`10⁻⁶×|q/p_prod|`，隔离重算后为0。误差检查只在其对应递推边界上定义；不将1,933解释为全部受影响states数。

同一批tracks、同一cuts，只改变smoothing的叠图：

![Same candidates, only smoothing changes](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/bump_diagnosis_20261112/figures/110051_smoothing_only.png)

## 详细数据：固定原生shared states，只改变选样

S3 L2完全不受本次smoothing重算影响，适合隔离选样效应。下表始终读取同一份shared states，按所列条件选track；阈值变化没有改动任何拟合值。

| 选择 | states数 | 负尾数 / % | 靠近−1数 / % |
| --- | ---: | ---: | ---: |
| 有可关联majority truth的全部selected candidates | 38,513 | 8,781 / 22.8001 | 4,977 / 12.9229 |
| 限signal e± | 37,222 | 8,695 / 23.3598 | 4,957 / 13.3174 |
| 再要求truth fiducial+H18 | 36,706 | 8,447 / 23.0126 | 4,777 / 13.0142 |
| 再要求selected SDO majority fraction≥0.7 | 27,780 | 3,220 / 11.5911 | 1,191 / 4.2873 |
| 按final matching为每truth选一条，未加Good95 | 27,762 | 3,206 / 11.5482 | 1,184 / 4.2648 |
| 完整旧图：final Good95祖先 | 25,691 | 2,675 / 10.4122 | 942 / 3.6667 |

前四行是逐级增加条件；后两行按旧图的final匹配/排序定义重新选祖先，不预设与前一行严格嵌套。SDO purity这一项本身就将负尾从23.01%降至11.59%，将靠近−1区域从13.01%降至4.29%。它主要移除了包含较多其他粒子hits的混合候选；这是一项truth辅助选样，不能把更小的尾部当作对全部原生CKF tracks的重建改善。

本批13,242条旧图轨迹中，只有1条final匹配barcode与selected-majority不同；排除此条后主要比例不变。因而大幅鼓包差别也不能归因于旧图普遍改用了另一个truth粒子的动量。

相同10文件下，新图与旧图两步变化的概览如下，数值均为smoothed负尾率：

| 层 | 新图：GitHub cuts、原生slots / % | 旧图选样、原生slots / % | 旧图选样、隔离重算 / % |
| --- | ---: | ---: | ---: |
| S1 L0 | 41.4126 | 37.4410 | 10.7539 |
| S1 L1 | 32.2672 | 28.4451 | 10.7629 |
| S2 L0 | 21.6708 | 17.0783 | 10.6657 |
| S3 L2 | 20.0603 | 10.4122 | 10.4122 |

第一、二列是各自选择的样本，不能按同state解释其差值；第二、三列为严格配对的smoothing单因素。旧图主体形状改变的位置依赖性由这两种作用共同解释。

两种选择均读取原生smoothed slots的叠图如下。此图比较各自raw Entries；上面的固定shared状态cutflow给出purity本身的比例变化。

![Native GitHub and S04 selection](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/bump_diagnosis_20261112/figures/110051_selection_comparison.png)

## 源码原因与结果的含义

控制开关定义在[CKF2.h:85](../../geoRotation_calypso/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.h:85)，默认`false`。实际重算位于[CKF2.cxx:383](../../geoRotation_calypso/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx:383)：

```cpp
destProxy.copyFrom(tracks.getTrack(selected.index), true);
if (m_isolateSmoothing) {
  // selected_shared records the copied original states here
  Acts::GainMatrixSmoother()(gctx, selectedTracks.trackStateContainer(),
                            destProxy.tipIndex(), *m_logger);
}
// selected records the independently smoothed candidate here
```

CKF不同候选可能共享早期祖先states；对不同后续分支进行smoothing时，共享的smoothed slots可被其他分支的后向结果覆写。候选的target fitted parameters已经单独保存，但事后读取早期state的smoothed值时，可能拿到与该候选后续hit路径不一致的值。越靠近共同祖先的S1影响越明显；S3末端分支的states大多独立。复制候选后，用该候选自己的predicted/filtered/Jacobian重算smoothing，恢复该路径的后向解。

因此`--isolate-smoothing`是会改变输出states的实际操作，不是单纯的绘图设置。旧图标注的candidate-isolated smoothing确实使曲线不同；新图为忠实复现GitHub及原生Calypso而关闭该开关。此前S04的100,000次事件配对已确认final p/charge和Good95判定不变、少量rmax元数据变化，见[原S04报告](2026-11-12_stage04_momentum.md)。本次解释逐层图的变化，不能据此宣称最终track动量问题已修复。

## 验证、运行与产物

实现为[smoothing_bump.py](../analysis/smoothing_bump.py)。全部10个归约jobs使用已保存的Calypso数据，经HTCondor执行；申请1 CPU、2300 MB、2 GB scratch，BLAS/OpenMP线程=1。每次提交比较实际schedd广告及队列，选兼容的最低负载bigbird20，排除已记录EOS兼容问题的bigbird21。归约cluster1573840、合并1573841、最终仅导图1573842，均exit0；任务、source manifests、软件快照和资源日志在`B`中。

实际提交入口如下；三次提交分别使用`tasks.json`、`merge_tasks.json`、`plot_tasks.json`及相应唯一`--name`：

```bash
python3 scripts/submit_momentum.py \
  --tasks /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/bump_diagnosis_20261112/plot_tasks.json \
  --library /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/software_after/lib \
  --name s04_smoothing_bump_plot_20261112 --stage S04-bump \
  --cpus 1 --memory-mb 2300 --disk-gb 2 --flavour microcentury \
  --exclude-schedd bigbird21.cern.ch
```

完整旧图右列27个P/F/S直方图的100个可见bins、finite/invalid数和tail数，均逐项精确复现。配对全过程断言hit、measurement flag、predicted和filtered向量一致。新图同文件的smoothed invalid reco为0；此前TTreeFormula零除保护产生的可见invalid entries仅在predicted，不解释此处smoothed鼓包。truth数值及AFS/strict敏感性结论沿用[truth单因素报告](2026-11-12_stage04_native_geometry_truth.md)。

首轮合并的统计和导出成功，但ROOT即时读取刚写入EOS的counts文件时报告`TNetXNGFile::ReadBuffer`错误。导出改为worker scratch中写入、读取和关闭全部ROOT/PDF，再复制最终文件至EOS；仅重做绘图，不改归约数据。

最终导图日志无ROOT读写错误；实测wall 84.25 s，user+system CPU 13.58 s，最大RSS 653,040 KiB。三组九层图共81个TH1的全部flow bins、raw Entries、分箱及颜色通过独立读取核对；shared/isolated对照summary，native对照相同10个worker的原宏计数。五张PNG通过完整性检查；5页汇报PDF与逐图PDF的content、各page box和rotation完全一致。五张PNG及PDF第5页渲染已目视检查，验证记录见[validation.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/bump_diagnosis_20261112/validation.json)。`python3 -m py_compile analysis/smoothing_bump.py`、116个本地文档链接和`git diff --check`通过。

主产物为[5页对照PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/bump_diagnosis_20261112/figures/smoothing_bump_comparison.pdf)、[汇总与逐bin counts](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/bump_diagnosis_20261112/summary.json)。其中包含原生选择、旧图选择+原生slots、旧图选择+隔离重算三组GitHub布局九层图，以及两步叠图。原图和控制图均使用英文标注。

## 阶段性结论

事实：S1的主要额外鼓包在固定旧图全部选样后，仍随候选隔离smoothing显著变化；S3 L2对该操作逐state不变，其主要变化由purity/final选样产生。文件覆盖、truth值和材料图差别均不是本次受控结果的原因。旧图更干净的形状不能直接代表全部原生CKF tracks具有同样表现。

## 下一步计划

保留两种用途的图：历史复现继续使用原生slots和GitHub cuts；逐运算物理诊断另列候选隔离smoothing及truth purity分组。若将隔离smoothing作为正式Calypso修复，独立验收其对rMax_smt/GitHub选择、refit、效率、最终动量、电荷和CPU的影响，并保留原生基线。局部物理偏差继续以同sensor G4 truth判别。

用户随后指定正式比较图开启 isolate-smoothing，其全110051复核已完成，见[全量 smoothing 与 seed covariance 对照](2026-11-12_stage04_smoothing_covariance.md)。新报告在原宏 cuts 下比较 beforeGeo/afterGeo，并将 afterGeo inflation 从 `[100,100,100,10000,100000]` 恢复到 `[100,100,100,100,1000]`：S2 P/F 负尾明显减少，smoothed 靠近−1的尾率在九层均略增。完整16页图集、选择迁移和事件配对区间随报告保存。
