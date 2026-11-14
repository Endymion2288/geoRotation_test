# S06：实际strip坐标与CKF测量误差的独立控制

**完整66事件、330次真实Calypso重建已通过归约。CKF原生直接使用固定σ=0.0231 mm，没有消费cluster自身随宽度变化的local-X方差。六条原始pure wrong中，CKF V×4与cluster V均恢复4条通过质量选择的正确轨迹，余2条仍错；SiHit centroid恢复5条，另1条符号正确但rmax失败。全部完整seed保持一致。四条误差控制恢复都伴随selected hit路径迁移；51839仅对坐标替换恢复，说明一个全局误差缩放不能统一解释六条。**

## 阶段目标

推进Q4中剩余pure轨迹的measurement机制：区分digitized cluster坐标、strip量化、同barcode的实际SiHit局部轨迹片段，以及CKF采用的measurement covariance。固定此前六条pure wrong与60条正常事件、corrected material encoding和isolated smoothing ON。结论须同时给出符号、原quality判据、完整seed及selected ordered hits的变化。

结果根目录`C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`；只读cluster控制记为`A=C/cluster_measurement_v1`，CKF误差控制记为`V=C/ckf_measurement_v1`。输入均为已有100116/100117 RDO，沿用原input file与entry身份；不重新生成模拟truth。

## 工作内容与实际源码

原生[CircleFitTrackSeedTool.cxx:99](/eos/home-x/xcheng/FASER/geoRotation_calypso/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:99)虽然读取`cluster->localCovariance()`，构造ACTS measurement时仍硬编码`0.0231*0.0231 mm²`。因此仅修改cluster error会影响其他组件，不能作为CKF V的单因素控制。本轮直接在实际Calypso该工具中提供默认关闭的`UseClusterMeasurementVariance`和默认1的`MeasurementVarianceScale`，不改seed C。

| 控制 | 实际Calypso位置与构建 | 默认行为与单因素 |
| --- | --- | --- |
| cluster/RDO/SiHit坐标记录 | `/eos/home-x/xcheng/FASER/geoRotation_calypso_charge_measurements/Tracker/TrackerRecTools/FaserSiClusterizationTool/src/TrackerClusterMakerTool.{h,cxx}`；`A/software` | 只读记录实际strip IDs/centres、SDO deposits、sensor-local SiHit起止点、width与error；默认关闭 |
| CKF测量方差 | `/eos/home-x/xcheng/FASER/geoRotation_calypso_charge_measurement_ckf/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.{h,cxx}`；`V/software` | 保留原固定V；独立切换V×0.25、V×4或cluster local-X V |
| tracking与observer | `V/overlay`组合上述两库与原`observer_hitcontrol_v2` | 完整native clustering、segment/seed、CKF、conversion、refit及观察器 |

两worktree均基于geoRotation commit`88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。CKF继承冻结`software_hitcontrol/software.json`的完整patch；所有新增修改保存在新worktree和新输出路径，没有覆盖既有库、任务或日常run。增量patch为[0030](../patches/0030_cluster_coordinate_error_audit.patch)，应用于已有[0009](../patches/0009_truth_cluster_position.patch)之后；CKF误差增量为[0031](../patches/0031_ckf_measurement_variance_control.patch)。两项patch分别在对应旧源码上`git apply --check`通过。

## 详细数据：只读开关与小样本门槛

cluster package重编4个translation units，Condor1578797 exit0。选择每电荷一个pure wrong及一个原normal，每事件audit OFF/ON各一组，共8个实际Calypso jobs，cluster1578802全部exit0。每组与corrected-material reference的58个数组字段精确相同；四个ON/OFF配对的58字段也逐值相同。记录的clusters分别为22、30、31、51个，全部cluster coordinates和error未因只读audit改变。

首轮validator1578803在解析JSON时失败：C++输出`Gaudi::Property<double>`时附带了字面属性名`'MeasurementVarianceScale':`。原数值token完整保留。[专用decoder](../scripts/cluster_coordinate_io.py)只去除该固定标签，严格JSON解析其余内容；没有改写原audit。新validator1578807 exit0，8组全部通过；不是重跑tracking修正数据。

CKF+geometry两个真实package重编53个translation units，Condor1578804 exit0。五模式为native、SiHit centroid+完整native seed、CKF V×0.25、CKF V×4及CKF cluster V。两电荷各wrong/normal一例，共20项重建1578810全部exit0；validator1578811与归约1578812均exit0。

20组均检查完整seed mean、36项C、surface及hit多重集合；四组native与corrected reference的58字段精确。逐selected measurement的实际V与请求factor/cluster V精确一致，actual measurement坐标与实际cluster local-X精确一致。局部G4信息是RDO携带的`FaserSiHit`起止位置，不含sensor处momentum；没有从HepMC production p冒称局部G4 p。

小样本两条wrong的分类如下，质量判据保持原选择：

| DSID/entry | native | centroid固定seed | CKF V×0.25 | CKF V×4 | CKF cluster V |
| --- | --- | --- | --- | --- | --- |
| 100116/51839 | wrong | correct | wrong | wrong | wrong |
| 100117/227965 | wrong | correct | wrong | correct | correct |

本表先报告分类；逐值q/p和selected hit路径见`V/analysis_v3/smoke/events.json`。两个normal在五模式下均correct。不能用20项配置作业当成20个独立物理事件。

## 完整66事件的五模式结果

全量归约为`V/analysis_v3/full`，cluster1578814 exit0；此前20项smoke与310项扩展重建1578813均exit0。沿用六条原始pure wrong（μ+2、μ−4）及60条normal（每电荷30），不以330次配置作业作为独立事件分母。原Good95、truth fiducial及charge-independent matching全部保持不变。

| 模式 | 六条pure：通过且correct / wrong / quality failed | pure中最终符号正确数 | pure保留原ordered hits | 60 normal：通过且correct / quality failed | normal保留原ordered hits |
| --- | --- | ---: | ---: | --- | ---: |
| native | 0 / 6 / 0 | 0 | 6 | 60 / 0 | 60 |
| centroid + fixed seed | 5 / 0 / 1 | 6 | 2 | 59 / 1 | 53 |
| CKF V×0.25 | 0 / 6 / 0 | 0 | 5 | 59 / 1 | 54 |
| CKF V×4 | 4 / 2 / 0 | 4 | 2 | 60 / 0 | 57 |
| CKF cluster V | 4 / 2 / 0 | 4 | 2 | 60 / 0 | 59 |

每个native与corrected-material reference的58个重建数组字段精确相同。全部330组的seed数量、六维mean、完整36项covariance、surface及hit多重集合精确；actual CKF measurement坐标与cluster used local-X逐值相同，实际V与请求的fixed/factor/cluster值逐值相同。三种V控制的所有cluster audit records与native逐值相同，因此没有将cluster误差或位置变化混入CKF V单因素。

下表为实际signed final q/p，单位GeV⁻¹，末位按显示精度取整。μ+的正确符号为正，μ−为负；星号表示符号正确但原质量选择失败。

| DSID/entry | native | centroid固定seed | CKF V×0.25 | CKF V×4 | CKF cluster V |
| --- | ---: | ---: | ---: | ---: | ---: |
| 100116/51839 | −0.00153412433 | +0.00188565176 | −0.00152459760 | −0.00153394311 | −0.00155926058 |
| 100116/224394 | −0.00146960213 | +0.00119396473 | −0.00146656435 | +0.000631693405 | +0.000643833690 |
| 100117/227965 | +0.000264435713 | −0.00208683550 | +0.000263736695 | −0.00160175301 | −0.00160981828 |
| 100117/390966 | +0.000407252880 | −0.00187325222 | +0.000404619750 | −0.00151662578 | −0.00152348415 |
| 100117/439421 | +0.00259627373 | −0.00217096846* | +0.00236711014 | +0.00278148955 | +0.00246252986 |
| 100117/458033 | +0.000327301812 | −0.00166515169 | +0.000322434333 | −0.00139499275 | −0.00138874140 |

V×4和cluster V恢复的是同四条：224394、227965、390966、458033；四条全部更换selected ordered hits，不能解释为在同一路径上单纯缩小曲率残差。51839在所有V控制中保留同一路径并仍错，centroid控制保留同一路径后恢复，是坐标偏移影响该轨迹的直接证据。439421的centroid控制保留同一路径、符号恢复，但`rmax<95 mm`失败；不能记为selected正确。

normal中的100116/109670仅在centroid控制质量失败，100116/405779仅在V×0.25质量失败；均缺少完整下游三站、至少七层及ndof+5≥12，符号仍正确。四种控制没有在这60条normal中新增错号。该小型、按结果选取的mechanism cohort不提供总体改善率或生产error最优值。

## 实际strip与SiHit证据

下表选各wrong轨迹native selected states中cluster-centroid偏移最大的measurement，所有坐标由实际Calypso sensor transform计算。

| DSID/entry；S/L/side | strips | primary SDO strips / 无SDO strips | cluster−primary centroid / mm | cluster σ / mm | primary SiHit数 / 连续片段数 | primary全部端点X范围 / mm |
| --- | ---: | --- | ---: | ---: | --- | ---: |
| 100116/51839；1/2/0 | 2 | 1 / 1 | −0.0453821004 | 0.0124707658 | 2 / 1 | 0.0001335144 |
| 100117/227965；2/0/0 | 6 | 1 / 5 | +0.2270492363 | 0.1385640646 | 1 / 1 | 0.0000715256 |

两处primary centroid对应的native cell都包含在该cluster中。将偏移按实际cell centre拆开：51839为`−0.040000−0.0053821004 mm`；227965为`+0.200000+0.0270492363 mm`。前项是cluster中心相对truth所在strip中心的偏移，后项是该strip中心相对truth centroid的量化位置差。实际pitch为0.08 mm；这不是以手写探测器模型代替native geometry。

全量六条pure的最大偏移measurement如下。所有native selected states共85个，其中7个cluster含29条无SDO strips；没有其他barcode的SDO支持。全部85个measurement均有同primary SiHit，所保存端点在10 nm连接标准下均为单一连续片段，X跨度最大0.0005502701 mm。端点连续性不包含RDO未保存的G4 track ID信息。

| DSID/entry；S/L/side | strips | primary支持 / 无SDO | 最大偏移 / mm | cluster σ / mm | primary端点X跨度 / mm |
| --- | ---: | --- | ---: | ---: | ---: |
| 100116/51839；1/2/0 | 2 | 1 / 1 | −0.0453821004 | 0.0124707658 | 0.0001335144 |
| 100116/224394；2/0/0 | 8 | 1 / 7 | −0.2793710187 | 0.1847520861 | 0.0000629425 |
| 100117/227965；2/0/0 | 6 | 1 / 5 | +0.2270492363 | 0.1385640646 | 0.0000715256 |
| 100117/390966；2/0/1 | 7 | 1 / 6 | +0.2319251368 | 0.1616580754 | 0.0002584457 |
| 100117/439421；1/0/1 | 4 | 1 / 3 | −0.1102360312 | 0.0923760431 | 0.0000419617 |
| 100117/458033；2/0/0 | 7 | 2 / 5 | +0.2098450955 | 0.1616580754 | 0.0004692078 |

60条normal的1048个native selected measurement中，31个cluster含40条无SDO strips；有1个measurement没有primary SDO支持，但同wafer存在同primary SiHit。两组绝对cluster-centroid偏移的中位数分别为pure 0.0176344490 mm、normal 0.0174216461 mm；最大值分别0.2793710187与0.1250849517 mm。常见单strip量化尺度相近，而pure组出现更宽、SDO支持稀疏的cluster。此为固定结果cohort的实际组成，不作为总体关联强度估计。

连续片段仅按保存SiHit的end→start距离`<1e−5 mm`连接，不将断开的同barcode片段自动认作secondary。当前两处端点跨度远小于cluster偏移，不支持由同primary在sensor内横向扫过0.2 mm来解释该偏移。track/cluster的majority purity=1也没有排除无SDO strips。

实际[FrontEnd串扰阈值检查:833](/eos/home-x/xcheng/FASER/geoRotation_calypso/Tracker/TrackerDigitization/FaserSCT_Digitization/src/FaserSCT_FrontEnd.cxx:833)可生成noise/crosstalk diodes，[addNoiseDiode:1011](/eos/home-x/xcheng/FASER/geoRotation_calypso/Tracker/TrackerDigitization/FaserSCT_Digitization/src/FaserSCT_FrontEnd.cxx:1011)使用`SiCharge::noise`；[addSDO:731](/eos/home-x/xcheng/FASER/geoRotation_calypso/Tracker/TrackerDigitization/FaserSCT_Digitization/src/FaserSCT_DigitizationTool.cxx:731)仅在real-particle hit或createNoiseSDO开启时保存对应SDO。该源码支持无SDO读出可以参与cluster的机制；当前RDO没有每个读出strip的噪声/串扰历史，不能把这五条strip唯一分配给某一来源。

## 执行、验证与产物

准备入口为[cluster小样本](../scripts/prepare_cluster_measurements.py)、[CKF五模式](../scripts/prepare_ckf_measurements.py)，构建入口为[cluster实际package](../scripts/build_cluster_measurements.py)与已有native-package builder。归约为[actual坐标/更新检查](../analysis/charge_measurement_coordinates.py)，图表为[测量误差图](../analysis/charge_measurement_error_figures.py)。全部批量运行、归约与绘图走HTCondor；每次重新比较schedd广告及队列，排除bigbird15/21，选择兼容且running+idle最少的bigbird20。

```bash
python3 scripts/prepare_cluster_measurements.py --root "$C" --source <new-cluster-worktree> --phase build
python3 scripts/prepare_cluster_measurements.py --root "$C" --phase smoke
python3 scripts/prepare_ckf_measurements.py --root "$C" --source <new-CKF-worktree> --phase build
python3 scripts/prepare_ckf_measurements.py --root "$C" --phase smoke
python3 scripts/prepare_ckf_measurements.py --root "$C" --phase replay
python3 scripts/submit_momentum.py --tasks <frozen-tasks.json> --library <matching-overlay/lib> \
  --name <unique-task-name> --stage S06-measurement --cpus 1 --memory-mb 3500 --disk-gb 4 \
  --exclude-schedd bigbird15.cern.ch --exclude-schedd bigbird21.cern.ch
```

完整五模式的66事件共330次重建已完成；全部完整seed、native58字段与5722个actual CKF measurement updates的坐标/V及运算门槛通过。独立统计1578816 exit0，数值与资源在[details.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/ckf_measurement_v1/details_v1/results/details.json)。330次tracking合计CPU 4054.21 s；每job CPU中位数13.005 s、P90 14.632 s、最大17.11 s；最大RSS 1023752 KiB。均为1CPU/3500MB申请，OpenMP/BLAS线程为1。

三页机制图1578815 exit0，已逐张全尺寸目视核对；独立validator1578820 exit0，plot values与冻结summary/events/coordinates精确相同，PNG完整解码、逐页vector内容/page boxes/rotation一致。[三页机制图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/ckf_measurement_v1/figures_v2/results/measurement_controls_report.pdf)与[逐图索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/ckf_measurement_v1/figures_v2/results/figure_index.json)均已交付。

汇编1578822 exit0，[101页S06汇报PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v8/charge_stage06_report.pdf)保留前版98页来源索引，追加本轮三页；合并vector内容与boxes逐页相同，视觉来源记录在[report validation](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v8/validation.json)。12项charge matching、seed lineage与paired regression检查通过；相关Python编译、局部文档链接、增量patch应用及`git diff --check`通过，记录在`C/verification_v13`。

![Five measurement modes](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/ckf_measurement_v1/figures_v2/results/measurement_outcomes.png)

![Actual cluster coordinates and errors](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/ckf_measurement_v1/figures_v2/results/measurement_coordinates.png)

![Six pure wrong curvatures](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/ckf_measurement_v1/figures_v2/results/measurement_pure_curvature.png)

## 阶段性结论与下一步计划

事实：CKF固定measurement error与实际cluster error不同；pure truth matching仍会保留含无SDO strips的cluster。四条pure错号对V放宽响应并更换selected路径，两条未被误差控制恢复；其中51839在固定同一路径后仅改变坐标即可恢复。六条的响应不支持一个全局error缩放统一解释所有wrong。假设：读出cluster中心偏移及其测量权重共同影响部分branch选择；无SDO strip相对strip量化的独立作用还需下一项控制。

下一项只用majority SDO正沉积支持的实际strip中心均值替换坐标，保留strip量化、cluster宽度/error/IDs及完整seed，以区分无SDO读出对cluster中心的作用和剩余量化误差；该truth辅助控制不作为生产修复。Q5完整共同物理先验仍未完成，S06保持未完成。
