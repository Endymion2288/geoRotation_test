# S04补充：参考图对照与材料旋转前后动量分布

> 本报告保留S04材料诊断结果；其Good95/SDO选择与GitHub原宏不同。用户最新指定的整体beforeGeo/afterGeo和仅替换truth的对照，按[更新后的规范](../roadmap/05_reference_figures.md)执行。

日期：2026-11-12。结果根目录为 `/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum/reference_figures`，下文记为 `$F`。

110051逐层图在本文使用前10个文件、10,000源事件，占该DSID目录50%；independent final图使用其后10文件。用户指定全文件后，已另行交付[20文件、20,000事件的完整图集与P/F/S研究](2026-11-12_stage04_full110051_pfs.md)，完整110051分布请使用该版本。

## 阶段目标

补齐`docs/problem.pdf`第3页风格的动量分布，展示材料图坐标严格旋转前后的seed、CKF逐层和final形状；核对用户提供的绘图源码；为S05–S07固定参考图、PNG/PDF汇编和数值核对要求。

比较的两端均为旋转后的探测器几何。`current`使用未旋转current22材料图，`strict`使用保持geometry IDs/bin data不变的严格旋转材料图。本次图中的before/after指材料图坐标，不能解读为历史beforeGeo/afterGeo的几何旋转净效应。

所有图的seed、CKF和final参数均来自geoRotation Calypso作业；本仓库脚本负责matching、归约和绘图。独立验证与逐层诊断使用的具体Calypso构建、patch、配置差别及材料问题源码路径，见[Calypso来源报告](2026-11-12_stage04_calypso_provenance.md)。独立ACTS toy没有进入本图集。

## 原绘图源码核对

已读取用户提供的[GitHub目录](https://github.com/zsm2024/FASER_tracking/tree/main/test_pResolution)，冻结commit为`33ca2281b961935456dfd4dc9fc06235e24b551e`。5个原文件与来源表保存在`$F/reference_scripts/`；下列链接固定到该commit。

机器可读核对记录见[source_compatibility.json](../../geoRotation_test_results/stage04_momentum/reference_figures/source_compatibility.json)。

| 项目 | GitHub脚本实际行为 | 本次补图的处理 |
| --- | --- | --- |
| Seed横轴 | [seeding_initialParams_residual.C L50–63](https://github.com/zsm2024/FASER_tracking/blob/33ca2281b961935456dfd4dc9fc06235e24b551e/test_pResolution/seeding_initialParams_residual.C#L50)：启用`abs(track_QoverP)-abs(t_QoverP)`，100 bins、[-0.03,0.03]；relative-p公式被注释 | 同时给PDF中的relative-p residual和当前源码启用的绝对curvature差；正确标出GeV⁻¹单位 |
| Seed取样与truth | 原宏遍历seed ROOT vectors，只跳过NaN；本地RootSeedWriter中`t_QoverP/t_theta`来自production particle传播到seed plane后的truth mean | S04主图保留fiducial+H18、SDO majority≥0.7的matched seed candidates，以production p/方向为参考；这是明确的truth/selection差异，不宣称与原seed ROOT逐bin相同 |
| CKF逐层横轴与单位 | [ckf_residual_perLayer.C L165–167](https://github.com/zsm2024/FASER_tracking/blob/33ca2281b961935456dfd4dc9fc06235e24b551e/test_pResolution/ckf_residual_perLayer.C#L165)：`(1/abs(eQOP)-1/abs(t_eQOP)/1000)/(1/abs(t_eQOP)/1000)`，reco为GeV、truth为MeV | 使用compact trace的GeV⁻¹ q/p和production truth GeV，计算等价的无量纲relative-p residual；不再次除1000。原writer的`t_eQOP`从该track majority粒子的SiHit link momentum构造，缺对应truth hit时为NaN；本图保留matched粒子的全部selected measurement states，以该粒子production p为参考。该参考不属于实测sensor局部p |
| 九层对应 | ACTS volumes=25/27/29，layers=2/4/6；每个物理layer可有两个strip sides的measurement states | 从hit Identifier解码station/layer，并逐state核对geometry ID的volume/layer与原宏映射一致 |
| CKF图色 | smoothed=blue、filtered=red、predicted=green | 保持相同颜色 |
| CKF分箱/纵轴 | [ckf_residual_perLayer.C L233–244](https://github.com/zsm2024/FASER_tracking/blob/33ca2281b961935456dfd4dc9fc06235e24b551e/test_pResolution/ckf_residual_perLayer.C#L233)：100 bins、[-1.2,1.2]、原始Entries，无Scale归一化 | 补齐原宏binning和raw counts的视图；另保留S04共同[-2,2)分箱、包含overflow的百分比视图，便于跨配置比较形状 |
| CKF质量选择 | 原宏为`chi2Sum/NDF<20 && rMax_smt<90 && NDF>9`，筛的是CKF track，未要求S04 final Good95 | 本次保留S04 fiducial+H18、matching与final Good95：χ²/ndof<25、layers≥7、ndof+5≥12、rmax<95 mm。原`rMax_smt`是smoothed positions的最大物理半径，不能拿hit radius或final rmax直接替代 |
| Good-seed开关 | [beforeVsAfter_trackWhichSeed_residual.C L525–530](https://github.com/zsm2024/FASER_tracking/blob/33ca2281b961935456dfd4dc9fc06235e24b551e/test_pResolution/beforeVsAfter_trackWhichSeed_residual.C#L525)：三项cuts当前均为注释；charge行即使启用也比较seed与自身 | 本次显式执行`abs(pseed/pprod−1)<0.1`、seed/truth同charge、同A-frame`abs(theta_seed−theta_prod)<1e−4 rad`。这是PDF所述条件的一套明确生产点定义，不等同原writer的propagated-truth theta |
| Good-seed逐层取值 | [同脚本L752–773](https://github.com/zsm2024/FASER_tracking/blob/33ca2281b961935456dfd4dc9fc06235e24b551e/test_pResolution/beforeVsAfter_trackWhichSeed_residual.C#L752)：`info.residual[idx]`在同layer每次写入，最后一个state覆盖前一个；当前h4填充也被注释 | 本次给出每个物理layer内全部measurement states的P/F/S分布，保留两侧。原宏每track/layer最终只保留一个state的覆盖行为未作为物理定义复制 |
| 事件与seed关联 | 原脚本使用ROOT entry的`i+1`与seed number在before/after间关联；原始路径在`/home/pub/...` | 本次使用DSID/file_index/source_entry/barcode定位，good-seed对照还要求真实seed mean逐值相同、selected hit ID序列相同；未把复用event number当作唯一身份 |
| 当前硬编码样本与标签 | `ckf_residual_perLayer.C`当前启用的文件是afterGeo的100043 μ−、4-station forward，但canvas/legend仍有110033、110044 μ+、beforeGeo文本 | 本次从冻结DSID和task配置生成标签，未复制这些硬编码文本；GitHub当前main也不能直接当作PDF当时完整运行配置 |

因此当前GitHub main、PDF第3页文字和S04冻结分析不是同一个完整运行定义。已对齐可验证的公式、颜色、九层映射，并提供原脚本的raw-count/binning视图；图中保留“S04 selection / production truth”标注。历史私有ROOT输入、本地保存的旧脚本版本及原writer实际配置未恢复，本次不报告历史图数值的精确复现。

两个seed横轴不可互换。例如p_truth=100 GeV、p_seed=50 GeV时，relative-p residual为−0.5，`abs(q/p_seed)-abs(q/p_truth)`为+0.01 GeV⁻¹；低估在两种图中位于相反一侧。

## 工作内容与样本

- 新增`analysis/momentum_reference_figures.py`，流式读取已有compact traces；没有重新模拟或重建。
- 五类DSID各10000个发现样本源事件，current/strict共100000次已有重建；独立预留样本另各50000源事件/配置。与S04相同，A′每个源事件可有两个signal truth。
- 初始seed图保留fiducial+H18 signal truth的全部多数匹配seed候选；CKF/final图使用各配置自己的Good95 matched final粒子，CKF取其真实祖先candidate。
- Good-seed九层图进一步要求两配置均为Good95、相同seed及selected measurement ID序列，且满足上表三项生产点参考cuts；这是条件分布，不用于估计T1母样本效率。含failure/pruning的T1效率仍使用S04既有定义与完整分母。
- 每层曲线的entries是measurement states；两侧pooled，不把states数量当作独立事件数。直方图保存原始counts、sumw/sumw2、invalid和under/overflow；百分比以全部有效值为分母。

## 详细数据与主图

独立验证的final Good95样本、tail=`p_reco/p_production<0.5`计数与S04冻结summary逐项核对：

| DSID | current Good95 / tail | strict Good95 / tail | tail比例current→strict / % |
| --- | --- | --- | --- |
| 100044 | 6959 / 774 | 6960 / 40 | 11.122→0.575 |
| 100116 μ+ | 4902 / 72 | 4918 / 23 | 1.469→0.468 |
| 100117 μ− | 5038 / 120 | 5041 / 20 | 2.382→0.397 |
| 110033 | 12254 / 1510 | 12702 / 903 | 12.323→7.109 |
| 110051 | 12509 / 4272 | 13288 / 1401 | 34.151→10.543 |

发现样本的matched初始seed候选在current/strict下逐bin相同。以下为每配置相同的计数；overflow定义为relative-p residual≥2，留在百分比分母中。

| DSID | finite seed candidates | seed tail<0.5 | invalid momentum | residual≥2 |
| --- | --- | --- | --- | --- |
| 100044 | 8168 | 110 | 16 | 162 |
| 100116 | 5764 | 17 | 3 | 61 |
| 100117 | 5950 | 42 | 2 | 51 |
| 110033 | 15045 | 1141 | 17 | 314 |
| 110051 | 16190 | 2224 | 46 | 808 |

对共同Good95 truth，先固定相同seed和selected hit序列，再施加上述good-seed条件。表中tail均相对production p，分母在两配置间相同；这组条件样本不承担总体效率的估计。

| DSID | 相同seed/hits粒子 | 其中good seed粒子 | CKF target tail current→strict | final tail current→strict |
| --- | --- | --- | --- | --- |
| 100044 | 6856 | 185 | 36→0 | 35→0 |
| 100116 | 4721 | 69 | 2→0 | 2→0 |
| 100117 | 4876 | 86 | 4→0 | 4→0 |
| 110033 | 10361 | 323 | 15→1 | 13→1 |
| 110051 | 9260 | 276 | 52→3 | 52→3 |

![Independent final residuals](../../geoRotation_test_results/stage04_momentum/reference_figures/figures/independent_final_residuals.png)

110051的seed、selected CKF target和final分布：

![110051 seed CKF final](../../geoRotation_test_results/stage04_momentum/reference_figures/figures/110051_seed_ckf_final_residuals.png)

110051九层P/F/S，左侧current，右侧strict：

![110051 nine layers](../../geoRotation_test_results/stage04_momentum/reference_figures/figures/110051_source_binning_nine_layers.png)

相同seed/hits且满足good-seed条件的九层图：

![110051 good seed nine layers](../../geoRotation_test_results/stage04_momentum/reference_figures/figures/110051_source_binning_good_seed_nine_layers.png)

与GitHub seed宏的两个横轴对照：

![110051 source seed views](../../geoRotation_test_results/stage04_momentum/reference_figures/figures/110051_source_seed_views.png)

完整产物：[汇报PDF](../../geoRotation_test_results/stage04_momentum/reference_figures/figures/momentum_material_comparison.pdf)、[逐图索引](../../geoRotation_test_results/stage04_momentum/reference_figures/figures/figure_index.md)、[分箱与统计](../../geoRotation_test_results/stage04_momentum/reference_figures/figures/summary.json)。

## 运行与验证

所有批量统计与导图均通过HTCondor执行，申请1 CPU、3500 MB memory、8 GB scratch，OpenMP/BLAS线程数=1。每次提交重新比较schedd负载，排除已证实不兼容当前EOS executable/log布局的bigbird21。

两次归约分别保存S04共同分箱与用户新提供源码的binning/raw-count视图；task、冻结worker、submission与resource位于`$F`。各组内部只流式读取一次，不保存额外states副本。源码文件仅作核对来源，不执行下载的宏。

```bash
python3 analysis/momentum_reference_figures.py --phase reduce --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum --dsid 110051 --output <new-reduction-directory>
python3 analysis/momentum_reference_figures.py --phase reference-reduce --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum --dsid 110051 --output <new-source-view-reduction-directory>
python3 analysis/momentum_reference_figures.py --phase plot --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum --reductions <all-five-reductions> --source-reductions <all-five-source-view-reductions> --output <new-figure-directory>
```

以上为worker执行的分析入口，完整批量提交使用`reduction_tasks.json`、`source_reduction_tasks.json`和`plot_tasks.json`。三次提交均选择`bigbird20.cern.ch`；实际资源如下，wall列为各组作业运行时间之和，不含排队时间，不能作为并行任务的端到端耗时。

| 作业组 | Cluster / 作业数 | CPU合计 / s | 作业wall合计 / s | 最大单作业RSS / MiB |
| --- | --- | --- | --- | --- |
| S04共同分箱归约 | 1573793 / 5 | 1700.76 | 1762.05 | 477.06 |
| 原脚本分箱归约 | 1573794 / 5 | 1549.39 | 1613.95 | 469.56 |
| 最终导图 | 1573795 / 1 | 182.96 | 198.07 | 206.51 |

11个作业全部exit 0。两组归约各覆盖100000次已有事件读取；逐组、逐DSID和配置的事件数、seed分布、good-seed人数及逐层P/F/S总体计数已交叉核对。共1170个直方图满足`finite entries = bin counts + underflow + overflow`，单位权重下`sumw=sumw2=finite entries`。独立样本的final entries和tail与S04冻结统计逐项一致，geometry ID与hit Identifier的九层映射检查通过。

最终实际导出31张PNG、31份逐图PDF及31页汇编，图索引与summary齐全。已目视检查五类独立样本总览，以及110051的seed/CKF/final、原脚本分箱九层图与good-seed九层图，标题、图例、英文标注、before/after顺序和页边裁切正常。布局验证另检查了原脚本seed双横轴图。

新脚本的Python语法与`--help`检查通过，已有11项回归测试通过；Markdown本地链接与`git diff --check`通过。验证记录：[统计守恒](../../geoRotation_test_results/stage04_momentum/reference_figures/statistics_validation.json)、[资源与导图结果](../../geoRotation_test_results/stage04_momentum/reference_figures/resources_summary.json)、[Markdown链接](../../geoRotation_test_results/stage04_momentum/reference_figures/markdown_links_validation.json)、[最终验收](../../geoRotation_test_results/stage04_momentum/reference_figures/validation.json)。

## 阶段性结论与下一步

事实：动量图对应的独立样本tail计数已在S04得到验证；材料坐标修正显著减少后续拟合的负尾，而相同输入seed不随材料查表坐标改变。good-seed条件分布与总体发生率分别展示，production参考仍包含真实沿程能损。

绘图定义方面，原GitHub脚本有启用公式、cuts和逐层覆盖与PDF文字不同的地方；本次在图注中明确实际定义。保留的S04选择让补图与已报告物理数值可互相核对。

S05、S06、S07的最终产物按[参考图清单](../roadmap/05_reference_figures.md)验收，分别覆盖角度、signed curvature/charge和完整修复前后图组。阶段开始先检查原脚本实际启用的公式/cuts及writer语义，冻结图清单；阶段完成时输出主图、PDF汇编、bin counts与中文报告。
