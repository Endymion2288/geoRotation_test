# S06：保留strip量化的SDO坐标控制与无truth-link的真实沉积

**固定完整seed、原cluster宽度/error和CKF measurement σ=0.0231 mm，仅将测量X换为majority SDO正沉积所支持的实际strip中心均值，六条原pure wrong中五条恢复selected正确；60条normal全部保持selected正确。439421在同一ordered hits下恢复。51839仍错，但final q/p从−1.53412×10⁻³变为−2.26365×10⁻⁵ GeV⁻¹；SiHit centroid才恢复其符号。无可用SDO贡献不等同于电子学噪声：29条pure无SDO strips全部与同wafer barcode=0的真实G4 SiHit端点X范围重叠。**

## 阶段目标

承接[实际measurement坐标与误差控制](2026-11-14_stage06_measurement_coordinates_and_errors.md)，区分宽cluster中缺少可用SDO贡献的strips对坐标的影响和剩余strip量化。验收要求保留原RDO标识、cluster宽度/covariance、完整seed均值及36项covariance，实际CKF消费坐标与固定measurement V逐值闭合，并同时报告路径迁移、质量失败和正常控制。

结果根目录`C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`，本控制记为`S=C/sdo_strip_coordinate_v1`，上一轮五模式控制记为`V=C/ckf_measurement_v1`。

## 工作内容与实际Calypso配置

输入沿用100116 μ+、100117 μ−的冻结RDO file/entry：六条pure wrong（两条μ+、四条μ−）和每电荷30条normal，共66源事件，分布于10个RDO文件。该cohort按既有重建结果选取，统计回答机制控制，不用于估计总体misID改善率。

实际源码基于geoRotation commit `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`，新增隔离worktree：

```text
/eos/home-x/xcheng/FASER/geoRotation_calypso_charge_sdo_position/
  Tracker/TrackerRecTools/FaserSiClusterizationTool/src/TrackerClusterMakerTool.cxx
  Tracker/TrackerRecTools/FaserSiClusterizationTool/src/TrackerClusterMakerTool.h
```

在旧cluster audit源码上新增默认关闭的`UseSDOSupportedStripPosition`。先按cluster的SDO电荷求majority barcode，再对该barcode有正deposit的实际strip，调用原生sensor `rawLocalPositionOfCell()`取中心，使用非加权X均值并加原Lorentz shift；本样本shift为0。缺少支持strip时保留native坐标。该控制与SiHit centroid互斥。增量见[0032 patch](../patches/0032_sdo_supported_strip_position.patch)，应用于0009+0030 cluster audit之后。

真实cluster package编译产物为`S/software`；实际tracking/geometry复用`V/software`、observer复用`C/observer_hitcontrol_v2`，完整组合见`S/overlay/software.json`。材料固定已修正decoder编码的`C/material_encoding_control_v1/corrected_afs_map.json`，isolated smoothing始终ON，CKF采用原生固定V，不切换cluster V。所有原始cluster IDs、width、error、localY、SiHit列表和channel deposits保持不变。

四事件native/SDO各一模式，共八项smoke重建；四条native在新cluster库下均与corrected参考的58项数组逐值相同。smoke validator通过后只提交剩余62项SDO控制。新运行合计70项tracking作业＝4 native基线+66 SDO；native和centroid的完整66事件结果复用上一轮已验收数据，没有重复运行该330项矩阵。

全量归约对66事件均检查完整seed的六维均值、36项C、reference surface及hit multiset一致，并核对实际加载库。1159次实际measurement updates的X等于新cluster X、V等于`0.0231² mm²`，更新恒等式通过。固定native selected measurement路径用于坐标统计，SDO新selected路径迁移另行计数。

## 详细数据：电荷与选样

| 模式 | 六条pure：correct / wrong / quality failed | pure保留ordered hits | 60条normal：correct / wrong / quality failed | normal保留ordered hits |
| --- | ---: | ---: | ---: | ---: |
| native | 0 / 6 / 0 | 6 | 60 / 0 / 0 | 60 |
| SDO-supported strip centres，固定seed | 5 / 1 / 0 | 2 | 60 / 0 / 0 | 59 |
| SiHit centroid，固定seed | 5 / 0 / 1 | 2 | 59 / 0 / 1 | 53 |

下表q/p单位为`10⁻³ GeV⁻¹`；same hits比较本模式与native的完整ordered measurement IDs。

| DSID / sample entry | native q/p | SDO q/p | SDO结果 / same hits | centroid q/p | centroid结果 |
| --- | ---: | ---: | --- | ---: | --- |
| 100116 / 51839 | −1.53412433 | −0.02263651 | wrong / 是 | +1.88565176 | correct |
| 100116 / 224394 | −1.46960213 | +0.99159060 | correct / 否 | +1.19396473 | correct |
| 100117 / 227965 | +0.26443571 | −1.93121400 | correct / 否 | −2.08683550 | correct |
| 100117 / 390966 | +0.40725288 | −1.91101843 | correct / 否 | −1.87325222 | correct |
| 100117 / 439421 | +2.59627373 | −2.10158410 | correct / 是 | −2.17096846 | rmax质量失败 |
| 100117 / 458033 | +0.32730181 | −1.69645689 | correct / 否 | −1.66515169 | correct |

439421的native/SDO/centroid最大半径分别为94.7058204、94.7469119、95.0914453 mm。SDO仍通过原`rmax<95 mm`，centroid符号虽正确却越界。四条其他SDO恢复伴随selected hit变化，不能将它们解释为固定同一路径上的单次局部修正。正常控制只有100116/309070发生ordered hit迁移，仍通过全部选择且符号正确。

51839保留原12个ordered measurements。最大偏移cluster的X由−28.48改为−28.44 mm，相对primary SiHit centroid的误差从−0.0453821004降至−0.0053821004 mm。其剩余错号最后发生于S3 L2的actual update：predicted q/p=`+0.0003603503 GeV⁻¹`，filtered=`−0.0003078224 GeV⁻¹`，residual=`+0.01985155 mm`，gain=`−0.03365846 GeV⁻¹/mm`；`Δq/p=K r`误差为`3.53×10⁻¹⁷ GeV⁻¹`。predicted q/p仅有0.1681σ的正号显著性，filtered σκ=`0.0018231741 GeV⁻¹`。这说明剩余测量/散射不确定性下的低显著曲率能够穿零；尚不能将这一例命名为新的程序缺陷，也没有将其从待解释类别移除。

## 详细数据：固定native路径的实际坐标

| 固定cohort | measurement states | bit级X变化 / 超过10⁻¹² mm | 无支持strip的fallback | native绝对误差median / max，mm | SDO绝对误差median / max，mm |
| --- | ---: | ---: | ---: | --- | --- |
| pure wrong | 85 | 28 / 7 | 0 | 0.0176344490 / 0.2793710187 | 0.0163345337 / 0.0374594116 |
| normal | 1048 | 326 / 30 | 1 | 0.0174216461 / 0.1250849517 | 0.0169790649 / 0.0926288308 |

误差参考为同primary、同sensor的energy-deposit加权SiHit端点中点，不是HepMC production momentum。纯strip均值的重新求和会带来约10⁻¹⁵ mm舍入差，不能把28/326个bit变化都解释为物理位移。实际pitch为0.08 mm，SDO控制仍保留strip量化；normal的多strip和无primary支持fallback不保证其误差小于半pitch。

全部37个超过10⁻¹² mm的坐标变化均发生于含无可用SDO贡献strip的cluster，没有其他非零barcode支持。pure的7个此类cluster全部改变；normal共31个，其中30个改变，另1个没有majority支持而fallback。

六条pure各自native偏移最大的cluster如下。表内SDO后的误差仍使用同一个原cluster和同primary centroid；barcode0的能量为同wafer全部此类SiHits之和，不擅自分配给某条strip。

| DSID / entry | strips / primary支持 / 无SDO贡献 | native / SDO误差，mm | 同wafer barcode0 SiHits | barcode0总Edep，MeV | 无SDO strip与barcode0端点X重叠数 |
| --- | --- | --- | ---: | ---: | ---: |
| 100116 / 51839 | 2 / 1 / 1 | −0.0453821004 / −0.0053821004 | 3 | 0.1609099992 | 1 |
| 100116 / 224394 | 8 / 1 / 7 | −0.2793710187 / +0.0006289813 | 20 | 0.3162399990 | 7 |
| 100117 / 227965 | 6 / 1 / 5 | +0.2270492363 / +0.0270492363 | 21 | 0.4533599969 | 5 |
| 100117 / 390966 | 7 / 1 / 6 | +0.2319251368 / −0.0080748632 | 16 | 0.3226300022 | 6 |
| 100117 / 439421 | 4 / 1 / 3 | −0.1102360312 / +0.0097639688 | 14 | 0.2648800001 | 3 |
| 100117 / 458033 | 7 / 2 / 5 | +0.2098450955 / +0.0098450955 | 10 | 0.1954099999 | 5 |

全部85个pure native states中，11个所在wafer有barcode0 SiHits，29条无可用SDO贡献的strips全部与这些step端点的X范围重叠。1048个normal states中108个有barcode0 SiHits，40条无SDO strips也全部有此X重叠。重叠使用实际strip中心±半pitch的区间；这是几何兼容性证据，没有逐strip charge的独占来源分解。

## 无SDO贡献与真实G4沉积的源码证据

这里“无SDO贡献”严格指audit未得到可用deposit列表；旧audit没有保存`has_sdo_entry`，不能据此断言map条目不存在。Calypso的实际机制为：

1. [FaserTrackHelper.cxx:32](../../geoRotation_calypso_charge_sdo_position/Simulation/G4Sim/FaserMCTruth/src/FaserTrackHelper.cxx:32)：没有trackInfo或HepMC particle时`GetBarcode()`返回0。
2. [FaserSctSensorSD.cxx:53](../../geoRotation_calypso_charge_sdo_position/Tracker/TrackerG4/FaserSCT_G4_SD/src/FaserSctSensorSD.cxx:53)：SiHit保存真实G4能量沉积、sensor局部start/end及particle link。
3. [FaserSCT_SurfaceChargesGenerator.cxx:410](../../geoRotation_calypso_charge_sdo_position/Tracker/TrackerDigitization/FaserSCT_Digitization/src/FaserSCT_SurfaceChargesGenerator.cxx:410)：barcode0 SiHit仍使用`SiCharge::track`，不自动丢弃其沉积。
4. [FaserSCT_DigitizationTool.cxx:731](../../geoRotation_calypso_charge_sdo_position/Tracker/TrackerDigitization/FaserSCT_Digitization/src/FaserSCT_DigitizationTool.cxx:731)：SDO创建的`real_particle_hit`要求非零barcode且type为track；[CreateNoiseSDO默认false](../../geoRotation_calypso_charge_sdo_position/Tracker/TrackerDigitization/FaserSCT_Digitization/src/FaserSCT_DigitizationTool.h:125)。只有无truth-link charge的strip可以影响实际读出，却没有进入已识别truth贡献。

例如51839最大偏移cluster的无SDO strip27中心为−28.52 mm；同wafer三个barcode0 step的X范围为[−28.5090580,−28.4345989] mm，Edep为0.16091 MeV。primary的X范围仅[−28.4346924,−28.4345589] mm、Edep=0.07767 MeV。227965的primary仅支持strip217，另外五条strips218–222的X范围均有barcode0 step重叠。strip是一维读出，局部Y不同的沉积仍可能进入相同X读出cluster。

RDO没有保存这些barcode0 steps的PDG、G4 trackID和parentID；它们与δ电子等无truth-link次级沉积相容，但尚未逐条确认粒子种类及parent。电子学噪声和crosstalk的逐strip贡献也未测量。实测足以否定“没有可用SDO贡献就一定是电子学噪声”，不支持把全部无SDO charge唯一归为某种secondary。

## 执行、资源与产物

全部批量build/reconstruction、归约、绘图及PDF汇编通过HTCondor。每次提交重新查询schedd并按兼容普通bigbird的running+idle最低选择bigbird20，排除已有兼容问题的bigbird15/21。build为1578821，八项smoke为1578824，validator为1578825，62项扩展为1578826，smoke归约为1578827，full归约为1578829，资源归约为1578830，均exit0。1578826的62条history均为Completed、ExitCode=0、ExitBySignal=false。

66项SDO tracking合计CPU=792.10 s，单项CPU median=11.525 s、P90=14.245 s、max=16.41 s；最大RSS=1,020,092 KiB。申请1 CPU、3500 MB、4 GB scratch，线程数为1。这些资源不包含四项native smoke或build。实际任务、选schedd依据、库和范围保存在各`*_worker/submission.json`及run manifests。

重建和分析入口为[prepare_sdo_strip_coordinates.py](../scripts/prepare_sdo_strip_coordinates.py)，归约为[sdo_strip_coordinates.py](../analysis/sdo_strip_coordinates.py)，英文图表为[sdo_strip_figures.py](../analysis/sdo_strip_figures.py)。已执行的phase不重复生成；复现实验使用新根目录并冻结相同上游输入与库。完整归约的提交命令为：

```bash
python3 scripts/submit_momentum.py \
  --tasks /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/sdo_strip_coordinate_v1/full_analysis_tasks.json \
  --library /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/sdo_strip_coordinate_v1/overlay/lib \
  --name s06_sdo_strip_full_analysis_20261114 --stage S06-measurement \
  --cpus 1 --memory-mb 3500 --disk-gb 4 \
  --exclude-schedd bigbird15.cern.ch --exclude-schedd bigbird21.cern.ch
```

主数据为[summary.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/sdo_strip_coordinate_v1/analysis_v1/full/summary.json)、同目录的`events.json`、`coordinates.json`及原始同wafer证据`cluster_evidence.json`；[resources.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/sdo_strip_coordinate_v1/resources_v1/results/resources.json)保存逐任务资源。

![SDO strip control outcomes](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/sdo_strip_coordinate_v1/figures_v1/results/sdo_strip_outcomes.png)

![Same-wafer actual G4 segments](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/sdo_strip_coordinate_v1/figures_v1/results/sdo_strip_g4_evidence.png)

[四页英文图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/sdo_strip_coordinate_v1/figures_v1/results/sdo_strip_controls_report.pdf)包含outcomes、六条signed q/p、固定native paths的坐标误差及六个最大偏移cluster的原始G4 step图；每图保存PNG、单页PDF和逐值plot数据。绘图1578831、独立validator1578833均exit0；plot_values与冻结归约精确一致，四页normalized vector content及page boxes与逐图PDF一致，四张PNG均原尺寸目视通过。G4图逐段绘制，不将同barcode0的不同SiHit擅自连接成同一G4 track。

[S06完整105页汇报PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v9/charge_stage06_report.pdf)已汇入原101页及本次4页，组装1578836 exit0。全部vector contents/page boxes和来源单页一致，前101页来源索引精确继承，后4页与本次图索引精确一致；[汇编验收](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v9/validation.json)保存原101页与新增4页的目视记录。

[代码与来源验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/sdo_strip_coordinate_v1/verification_v1/validation.json)确认四份Python代码编译通过，0032在确切旧cluster基底上应用后与实际编译的两个源文件逐字节相同，完整source.patch与当前冻结worktree差异精确，reconstruction/analysis snapshots与仓库代码精确。新增文档和入口链接、`git diff --check`通过；本次以实际smoke/full gates验证新物理控制，没有无变更地重复上一轮既有回归矩阵。

## 阶段性结论与下一步计划

事实：保留strip量化、完整seed及原固定CKF误差后，移除无majority SDO支持的strips对X均值的影响即可恢复五条pure selected正确，其中439421保留同一条hit路径。四条包含路径迁移。剩余51839的极小曲率仍对量化后的measurement updates敏感，SiHit centroid在同hits下恢复。仅用SDO majority purity=1不能排除无truth-link沉积引起的cluster中心偏移。

下一步在真实Calypso中继续固定hit路径、measurement残差与过程噪声的可判别控制，分开测量51839的strip量化和散射对穿零的贡献；避免由单条低显著曲率推断全样本偏差。Q5继续完成共同有限物理先验下的非线性CKF坐标控制，并与已完成的7,680点transport闭合明确区分。生产cluster/error修复与全量效率、pull及留出μ±/A′验收在S07执行；S06保持未完成。
