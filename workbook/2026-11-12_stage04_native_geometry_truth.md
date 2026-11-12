# S04复核：原生Calypso几何对照与truth单因素检查

日期：2026-11-12。结果目录为`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison`，下文记为`G`。该目录保存S04动量复核，不表示已开始路线图S05的100044角度研究。

**仅替换同一majority粒子的truth动量来源，不能解释旧图与参考图的大幅差别。** afterGeo全20文件的566,369个共同有效states上，两种q/p的最大相对差为`5.924649615×10⁻⁸`，原宏27个直方图的所有bin完全一致。production reference补足5,822个原来没有SiHit truth的states，改变约1.02%的覆盖；负尾率最大变化为0.679个百分点。两种reference均为关联粒子的HepMC产生点动量，并非sensor处的局部p。

## 阶段目标

在110051样本全部20个RDO文件、20,000源事件上，使用GitHub `test_pResolution`原宏复核动量分布。beforeGeo/afterGeo指整体坐标旋转前后；AFS/strict材料图单列对照。用户最新提出的truth检查固定同一条writer track、同一majority barcode及所有其他绘图设置，仅替换truth动量，并将共同有效states与缺少SiHit truth的覆盖变化分别统计。

第17页是本次动量图的规范；第19页的100044 μ+500k角度图留到后续阶段。本次不追加100044重建。

## 软件、配置与绘图来源

所有tracking均由Calypso原生组件执行，配置入口为[reference_tracking.py](../configs/reference_tracking.py)。没有独立tracking或简化detector/field实现。

| 项目 | beforeGeo | afterGeo |
| --- | --- | --- |
| Calypso源码 | `/eos/home-x/xcheng/FASER/geoRotation_calypso_main`，本地master的隔离worktree | `/eos/home-x/xcheng/FASER/geoRotation_calypso` |
| Commit | `1eb885ae22299fb39bb352bead6005d41b0f953e` | `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` |
| 软件 | Athena24.0.41、ACTS32.0.2、EL9/GCC13 | 同左 |
| geometry/conditions | FASERNU-04 / OFLCOND-FASER-05 | 同左 |
| tracking模式 | 3station forward，mask0–5 | 同左 |
| 材料 | `stage02_closure/material-current-source.json`，F frame | AFS `material-maps-alma9-rotated.json`，A frame；strict JSON为独立控制 |
| seed variance inflation | `[100,100,100,100,1000]` | `[100,100,100,10000,100000]` |
| audit/reset-event-maps/isolate-smoothing | 全部关闭 | 全部关闭 |

beforeGeo重编真实Calypso的`FaserActsGeometry`、`FaserActsKalmanFilter`两个packages，包括geometry interfaces头文件；共53个translation units，产出三个原生库。复用相同Athena release中的未变依赖，源码tracked clean。脚本为[build_native_calypso.py](../scripts/build_native_calypso.py)，完整compile/link recipes与`software.json`位于`G/software_before/`。移除debug symbols以减少产物大小，保留原优化及物理代码。afterGeo使用S02保留的未patch原版KalmanFilter库，Build ID为`15dc7462b839ae6cce6254e5e309b6a7a2005ad6`；Python configs从上述commit导出。原生两个分支的covariance不同，因此这是一项整体版本比较，不能将其全部差别解释为保持先验一致的纯数学旋转效应。

原宏来源为[GitHub/test_pResolution](https://github.com/zsm2024/FASER_tracking/tree/33ca2281b961935456dfd4dc9fc06235e24b551e/test_pResolution)，冻结commit为`33ca2281b961935456dfd4dc9fc06235e24b551e`。主图直接调用`ckf_residual_perLayer.C`，保留实际启用的设置：

- `chi2Sum/NDF<20 && rMax_smt<90 && NDF>9`，逐层`volume_id=25/27/29`、`layer_id=2/4/6`。
- 原始relative-p residual公式，100 bins、[-1.2,1.2]、raw Entries；P=green、F=red、S=blue，3×3九层布局。
- `eQOP_prt/flt/smt`单位为GeV⁻¹，truth branch单位为MeV⁻¹；保留原公式中的`/1000`。

[prepare_reference_macro.py](../scripts/prepare_reference_macro.py)仅适配输入、输出、标签、truth alias、共同纵轴和合并TH1的读取。原宏完成九层canvas之后有`hists[3]`等越界的废旧绘图块，适配时删除该部分；不改active cuts或residual公式。原版、适配版和完整diff保存在`G/reference/`。worker用原宏的`TTree::Draw`生成27个TH1，再合并counts并调用相同绘图代码。全量ROOT states和xAOD只留worker scratch。

## 两种truth实际表示什么

原writer的详细路径为[RootTrajectoryStatesWriterTool.cxx](../../geoRotation_calypso/Tracking/Acts/FaserActsKalmanFilter/src/RootTrajectoryStatesWriterTool.cxx:425)。它先从SDO找到track的majority barcode，再按`(barcode, waferId)`找SiHit；找到后累加`siHit->particleLink()->momentum().rho()`和charge并构造`t_eQOP`，缺少对应SiHit时写NaN。

`HepMcParticleLink::operator->()`返回HepMC `GenParticle`，不是G4 step。接口完整路径为`/cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/x86_64-el9-gcc13-opt/include/GeneratorObjects/HepMcParticleLink.h:476`；Calypso的[FaserSiHit.h](../../geoRotation_calypso/Tracker/TrackerSimEvent/TrackerSimEvent/FaserSiHit.h:105)保存start/end position、energy deposit、time及particle link，没有保存sensor处的动量。因此，现有`t_eQOP`不能称为“实测sensor局部p”。多SiHit时writer的`nHits=count(key)`为0/1；同一粒子的charge与p同时累加通常相消，该代码也没有恢复实际局部能损后的p。

本次[reference_truth_observer.py](../scripts/reference_truth_observer.py)只观察原生Calypso的xAOD `TruthParticles`和最终collection顺序；不改变tracking。production q/p取相同majority barcode的`charge/sqrt(px²+py²+pz²)`，仍以MeV⁻¹写ROOT friend。关联使用输入文件、source entry和selected-track ordinal，逐条断言writer event/track编号；不使用会重复的event number作唯一key。主、次级truth分别计数。[reference_truth.py](../analysis/reference_truth.py)另生成`link_common`与`prod_common`，只在两种参考都有限且非零时保留值。

四个绘图输入为原`t_eQOP`、全production、共同有效states的link和共同有效states的production。原宏未增加finite cut；NaN等invalid另行统计，ROOT的flow bins按原行为保留。有限state分母与图中的raw Entries分别报告。

## 小样本核对

原生afterGeo重建file0前20个事件成功，wall time397.112 s；42条writer tracks，其中30条通过GitHub cuts。9层共525个selected measurement states，两种truth均有效；majority均为generator电子或正电子，14条e−、16条e+。

| 检查 | 实测结果 |
| --- | --- |
| 原`t_eQOP`相对production q/p最大差 | `4.587465674×10⁻⁸` |
| 共同有效states中相对q/p差>10⁻⁶ | 0/525 |
| 原truth与production的27个TH1 | 全部bin counts一致，含under/overflow |
| 四种truth×27个TH1的独立核对 | 108/108全部flow bins一致 |
| xAOD production与旧compact trace的HepMC production | 前20事件61个可关联粒子的p逐值相同，最大相对差0 |

最后一项使用同一输入的S04 `TrackingAudit.h`输出，仅比较truth，不混用该旧重建的reco状态。它确认此次production reference与用户提到的compact trace引用同一类产生点动量。小样本结果支持writer的float `t_eQOP`与production q/p只有数值精度差异；全样本结果另列。

验证产物为`G/smoke_after/truth_comparison.json`、`macro_histogram_comparison.json`、`independent_histogram_validation.json`和`production_encoding_check.json`。独立核对命令：

```bash
python3 scripts/validate_reference_output.py /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/smoke_after
```

## 全量truth单因素结果

主设置afterGeo使用用户指定的AFS JSON。strict材料图仅为控制，truth比较在每个设置内部使用同一份writer ROOT，完全固定reco、majority barcode、cuts与分箱。

| 指标 | afterGeo / AFS JSON | afterGeo / strict JSON |
| --- | ---: | ---: |
| 输入文件 / 源事件 | 20 / 20,000 | 20 / 20,000 |
| 原生CKF writer tracks | 41,940 | 42,160 |
| 通过GitHub cuts的tracks | 33,178 | 33,324 |
| majority generator / secondary或缺失 | 32,025 / 1,153 | 32,169 / 1,155 |
| majority e− / e+ / π− / π+ / 缺失 | 16,786 / 16,387 / 1 / 1 / 3 | 16,865 / 16,455 / 1 / 1 / 2 |
| 九层selected measurement states | 572,234 | 574,685 |
| 共同有效link / production states | 566,369 | 568,746 |
| production有效states | 572,191 | 574,659 |
| production-only / link-only states | 5,822 / 0 | 5,913 / 0 |
| production-only占全部states | 1.017416% | 1.028911% |
| 共同states相对q/p差>10⁻⁶ | 0 | 0 |
| 最大相对q/p差 | 5.924649615×10⁻⁸ | 5.924649615×10⁻⁸ |
| common link / production逐bin差 | 27/27图完全相同，含flow bins | 27/27图完全相同，含flow bins |

原始`t_eQOP`为float，production friend为double；观察到的差异与float精度相符。共同states上不只是主峰接近，而是100个可见bins及under/overflow逐bin一致。全production图的变化来自新增truth覆盖。AFS设置中，27图合计1,603个bins变化、bin差绝对值总和22,002；该数包含NaN从overflow移入有限分布以及同一state的P/F/S三次计数，不是22,002条独立轨迹。

下表的负尾为`p_reco/p_truth<0.5`，分母为该层、该stage的有限residual数；它与原宏画布中的raw Entries不同。表中每格为原truth→production truth，单位为%。

| AFS设置的物理层 | predicted | filtered | smoothed |
| --- | ---: | ---: | ---: |
| S1 L0 | 27.6190→27.1336 | 27.6141→27.1289 | 42.6341→41.9555 |
| S2 L0 | 49.3973→49.1730 | 66.2769→65.9751 | 21.9704→21.8769 |
| S2 L1 | 64.9416→64.7521 | 62.9822→62.7910 | 21.4323→21.3859 |
| S2 L2 | 59.8264→59.7279 | 55.6018→55.5060 | 21.0378→21.0244 |
| S3 L2 | 19.8141→20.2354 | 19.7915→20.2194 | 19.7706→20.1985 |

最大变化发生在S1 L0 smoothed，为−0.678543个百分点；strict控制的最大变化为−0.689147个百分点。S2中P/F的大负尾仍然存在，切换truth没有消除它。各层state数与候选之间相关，以上比例是描述性统计，不将数十万states当作独立事件计算显著性。

原truth九层主图：

![110051 afterGeo native truth](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/figures/110051_afterGeo_native_nine_layers.png)

仅切换truth的叠图，实线为SiHit link、虚线为production：

![110051 truth-only comparison](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/figures/110051_afterGeo_truth_overlay.png)

[共同有效states叠图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/figures/110051_afterGeo_common_truth_overlay.png)中两组逐bin重合。详细数据为`G/reductions/afterGeo_truth_summary.json`、`afterGeo_layer_metrics.csv`；strict控制对应同目录`afterGeo_strict_*`。原宏下的S1双峰、S2大负尾和S3的P/F/S逐渐重合，与第17页右图的主要结构相符；历史私有ROOT及完整运行manifest未取得，未作历史逐bin相等声明。

## 整体坐标旋转前后的全样本结果

beforeGeo与afterGeo分别使用两个原生版本自己的配置及GitHub cuts，输入完全相同，均为20文件、20,000事件。这里比较每个配置实际选中的分布，未强制使用共同最终tracks。

| 指标 | beforeGeo / master | afterGeo / geoRotation、AFS |
| --- | ---: | ---: |
| writer tracks | 44,238 | 41,940 |
| GitHub selected tracks | 36,158 | 33,178 |
| 九层measurement states | 623,285 | 572,234 |
| 共同有效link / production states | 616,804 | 566,369 |
| production-only states | 6,448（1.034519%） | 5,822（1.017416%） |

beforeGeo上616,804个共同有效states的最大相对q/p差同样为`5.924649615×10⁻⁸`，27个common link/production直方图也逐bin完全相同。truth定义等价的结论在两个原生版本中均成立。

下表固定原`t_eQOP`；每格为beforeGeo→afterGeo的有限residual负尾率，单位为%。

| 物理层 | predicted | filtered | smoothed |
| --- | ---: | ---: | ---: |
| S1 L0 | 25.7876→27.6190 | 25.7755→27.6141 | 38.5551→42.6341 |
| S1 L1 | 25.8456→27.8533 | 25.8367→27.8884 | 38.2280→32.8883 |
| S1 L2 | 25.4933→28.9515 | 25.4786→31.7958 | 37.8350→31.6215 |
| S2 L0 | 31.6471→49.3973 | 38.9531→66.2769 | 42.0479→21.9704 |
| S2 L1 | 40.0969→64.9416 | 40.1607→62.9822 | 42.1766→21.4323 |
| S2 L2 | 39.8896→59.8264 | 40.2256→55.6018 | 41.6898→21.0378 |
| S3 L0 | 40.8477→37.2550 | 42.5395→20.7562 | 43.6199→20.5507 |
| S3 L1 | 43.1108→20.7366 | 43.3273→20.3526 | 43.5394→20.3350 |
| S3 L2 | 43.0366→19.8141 | 43.1204→19.7915 | 43.1953→19.7706 |

整体版本变化确实显著改变了P/F/S之间的结构。afterGeo的S2 P/F尾率更大，smoothed尾率更小；末端S3 L2三个stage接近。该事实不能单独区分几何、seed covariance及候选选择的贡献，需要固定物理先验和共同seed/hits的后续消融。整体旋转图不等同于过去只改变材料坐标、同时启用诊断修正的S04对照。

beforeGeo原truth主图：

![110051 beforeGeo native truth](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/figures/110051_beforeGeo_native_nine_layers.png)

整体版本叠图，实线为beforeGeo、虚线为afterGeo：

![110051 geometry comparison](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/figures/110051_geometry_overlay.png)

predicted使用当前测量更新前的状态，filtered加入当前测量，smoothed再利用后续测量进行反向修正。S2的P/F相当于只使用较短的前向信息，动量受seed及早期曲率约束影响；smoothed具有下游测量的信息。但“后向信息更多”不能证明所有尾部都正常：S1 L0 smoothed尾率反而更高，必须检查同一state的P→F→S迁移、covariance/gain及同面局部truth。此前相应逐运算研究见[全样本P/F/S报告](2026-11-12_stage04_full110051_pfs.md)，其诊断选择与这里的原宏样本不同，不能直接合并比例。

## 历史复现中的已知差异

旧`110051_nine_layers.png`来自matched final Good95祖先，并使用fiducial+H18、SDO≥0.7、每truth选一条、production truth以及reset-event-maps/isolate-smoothing。GitHub原宏直接画native CKF states writer中的selected tracks，只加上述三个cuts。两者所选tracks、states和诊断控制不同；旧S04图不能直接当作第17页的原宏复现。本报告的truth实验只固定native writer的barcode；后续[同10文件鼓包复核](2026-11-12_stage04_smoothing_bump_cause.md)已测量smoothing及purity/final选择的作用，并检查旧图final/selected的barcode归属：13,242条轨迹中仅1条不同。S1的主要额外鼓包随候选独立smoothing消失，S3末端主要随选样变化；其余配置因素仍按对应控制的覆盖范围解释。

报告附带的公开脚本位于`/afs/cern.ch/work/s/sazhang/public/FASER/xai_alma9/check_geometryRotation/faser_reco.py`，已冻结到`G/reference/historical_faser_reco.py`。其中所谓3station的mask为`[0,1,2]`，而seed mask按`6*station+2*layer+side`编码，完整IFT排除应为0–5。当前GitHub宏实际启用的hardcoded输入又是100043 μ−4station；注释才包含110051。这些是可核对的配置差异，PDF截图不能确认当时私有ROOT具体使用哪一版脚本。控制样本使用file0的1,000事件，分别测试mask012与4station，不能将其标为20文件全量。

该1,000事件控制实测三个mask设置均为2,046条writer tracks、1,664条GitHub selected tracks，27个原truth直方图含flow bins逐bin相同。因此，这个mask差异在本控制样本上没有改变本次图，不能将其列为已证明的复现差异来源。

AFS与strict材料设置的20文件结果有小差异：selected tracks为33,178与33,324，相差146条（相对AFS为0.4401%）。27个原truth直方图共有2,166个bins不同，bin差绝对值总和11,169，最大单bin差98 entries。共同12个材料surfaces的86,400个材料bins相同，但strict还包含10个空surfaces及法向平移±2 mm，详见原材料图核对报告。变化确实存在；这些额外surfaces具体如何通过navigation或candidate路径影响结果尚未定位。更换JSON不能解释旧S04图与原宏图之间的大幅结构差异。

## 执行与验证记录

完整输入清单和配置矩阵为`G/inputs.json`、`beforeGeo_tasks.json`、`afterGeo_tasks.json`、`afterGeo_strict_tasks.json`。复用S04已逐文件测量的CollectionTree entries，并重新核对目录完整文件清单和文件大小。每个全量配置20 jobs，每job一文件1,000 events；reconstruction单event slot，BLAS/OpenMP线程=1。申请1 CPU、3500 MB、8 GB、workday；CERN实际资源调整另见job广告和resource。

提交前分别查询schedd负载与队列，选择兼容的普通bigbird中running+idle最少的`bigbird20.cern.ch`；已知EOS兼容性问题的bigbird21排除。每次查询、选择依据和资源保存在相应`*_worker/submission.json`、`schedd_load.txt`与`schedd_queue.txt`。native build最终cluster1573823，afterGeo/afterGeo_strict为1573820/1573821，mask/4station控制为1573822，beforeGeo为1573827。

构建适配最初在package筛选及relative library path解析时失败，第三次把`-fplugin=...so`与`-Wl,-soname,...so`误当作文件路径。修正为区分linker flags与文件，并在编译前验证link recipes；v4完成53 units及三个库的链接。日志为`G/build_job_v*/reconstruction.log.gz`。绘图worker没有`pdfunite`；最终使用pypdf 5.9.0合并原始单页ROOT vector PDFs，保留每页内容与crop boxes，软件冻结于`G/software_tools/`。PNG、逐图PDF和counts独立核对。

第一次worker导出没有给`TH1::Write`显式key，ROOT保留原对象名`hist_station_layer`及S/F/P对应的cycles 1/2/3。合并器读取该实际格式时同时校验line colors 600/632/416，最终counts使用明确`h_component_station_layer` keys；后续worker源脚本已改为显式key。AFS file0的native ROOT另因EOS `TFile::Flush: No such device`只写出63 bytes，导致analysis cluster1573832失败。该作业的`common_native` ROOT完好，且selected missing production=0、每层link-only=0，原始link值和common link值对所有states完全相同。逐层P/F/S的finite、sum、sum2、tail、under/overflow和invalid-reco统计逐值相等，按这一等价关系恢复native导出；27个histogram的binning、颜色和entries再次通过读取检查。原截断文件、错误与全部相等断言保存在`G/export_recovery/afterGeo_f00/`。没有重建或修改任何物理值。源worker及合并器改为先在worker scratch关闭ROOT，再复制到EOS，写出失败由文件系统异常报告；AFS复算cluster1573836成功。

15项现有回归测试通过，记录在`G/unit_tests.log`；新增配置和脚本通过`py_compile`，改动通过`git diff --check`。

最终合并/绘图cluster1573838成功，60个全文件重建jobs及2个1,000事件控制jobs均exit0。冻结代码为`G/code/analysis/reference_figures_final.py`，任务入口为`G/final_analysis_tasks.json`。本次实际提交命令：

```bash
python3 scripts/submit_momentum.py \
  --tasks /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/final_analysis_tasks.json \
  --library /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/software_after/lib \
  --name s05_complete_reference_analysis --stage S04-reference \
  --cpus 1 --memory-mb 3500 --disk-gb 2 --flavour microcentury \
  --exclude-schedd bigbird21.cern.ch
```

六个合并case（含file0的nominal alias）共648个TH1通过raw Entries、有限residual及flow bins的统计闭合。原宏的`TTreeFormula`除零返回0；因此reco q/p=0或±∞可被原式填到residual=−1。原truth下，beforeGeo/afterGeo/strict分别有33/16/17个这样的invalid reconstructed states进入可见bins，均在predicted；主图原样保留，而有限residual指标排除。此次统计闭合将这部分单列，没有新增cut。独立NumPy校验也保留该原公式语义，smoke的108个TH1全部flow bins再次一致。

| 资源实测 | beforeGeo | afterGeo / AFS | afterGeo / strict |
| --- | ---: | ---: | ---: |
| 20 jobs重建wall之和 / s | 36,618.819 | 31,667.967 | 29,092.252 |
| 单job重建wall范围 / s | 1,311.434–2,142.603 | 1,289.109–2,219.132 | 1,153.959–1,666.263 |
| 20 jobs含truth归约的worker wall之和 / s | 39,237.047 | 34,130.517 | 31,497.976 |
| 单job最大RSS的最大值 / KiB | 2,360,744 | 2,330,792 | 2,342,040 |

以上wall之和不是并行任务的实际历时，也不直接用作版本CPU性能比较；运行节点与I/O不同。Condor history确认三个主配置的file0均实际分配2 CPU、6000 MB、8 GiB disk；reconstruction仍为单event slot、单计算线程。每配置约57 GB的临时production JSON以及完整ROOT states和xAOD在worker scratch完成归约，不写入EOS长期保存。

所有14张PNG已经完整性检查和contact sheet目视检查；14页汇报PDF逐页验证与单页ROOT vector PDF的content、MediaBox、CropBox和rotation一致，并渲染核对afterGeo主图。最终产物：

- [14页完整汇报PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/figures/reference_comparison.pdf)，包含三个主case各4图、整体geometry叠图及材料叠图。
- [逐图PNG/PDF/ROOT索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/figures/figure_index.md)；bin counts为`G/reductions/*_3station.root`，逐层表为`*_layer_metrics.csv`。
- [全量覆盖与验证记录](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/validation.json)及`G/figures/figure_validation.json`；[输入与配置manifest](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/inputs.json)。

## 阶段性结论

事实：相同majority barcode下，原`t_eQOP`与显式production truth在两个原生版本中等价，剩余差异为float精度和约1%的truth覆盖。使用GitHub原宏及native CKF writer后，可重现第17页左右两图的主要结构。旧S04图存在已核对的选样、stage与诊断设置差异；AFS/strict还有小的产额变化，file0的三种mask控制则逐bin一致。

待检验的解释：旧图中的candidate/final-truth归属及各项cuts各自贡献尚未量化；AFS/strict的小差异如何经navigation或候选选择产生尚未定位；两个原生版本的几何与covariance贡献尚未分离。它们不能由本次truth值等价实验替代。

## truth选择建议与后续物理检验

| 分析目的 | 建议reference | 原因与本次实测 |
| --- | --- | --- |
| 复现geometryRotation第17页和GitHub原图 | 原`t_eQOP` | 完整保留历史宏口径；本次主图采用这一选择 |
| 评价相对产生点的整体动量响应 | 显式production truth，保持相同majority barcode | 参考含义更清楚，并补足约1%的states；同时报告缺失、混合hits、majority purity和共同覆盖对照 |
| 判断某层prediction/update/smoothing本身是否有bias | 同sensor、同次穿越的G4 pre/post-step momentum | 分离实际沿程能损与reco误差；上述两种reference都不满足这个局部定义 |

对于电子样本，相对production p的负尾包含真实沿程能损和拟合偏差，不能直接命名为“该层KF动量低估率”。production补足的states也可能对应majority粒子并未在该wafer留下hit的混合轨迹；覆盖更完整本身不证明这些states是纯净匹配。用于定位算法时保留共同有效states及truth purity分组。110051现存RDO没有sensor局部p；局部偏差的定量检验使用Calypso内带G4动量记录的诊断simulation，按既有sensor audit方案关联。

## 下一步计划

1. 使用本次冻结的20文件和Calypso native writer，逐项控制旧图的candidate/final-truth归属及cuts，报告配对轨迹数、有限分母、tail迁移和原宏图，区分选样与数值变化。
2. 在Calypso内固定共同seed/hits及物理covariance，检查P→F→S首次tail迁移；以同sensor G4 pre/post-step p分解真实能损与拟合偏差，给出每类事件的比例和可复查输入定位。AFS/strict导航差异作为独立小效应控制。
3. 第19页角度和μ±误判研究按[路线图](../roadmap/README.md)推进，分别验收同面角度结构、谱匹配的误判率及区间。本次完成范围为110051；材料图bug作为已解决的历史控制。
