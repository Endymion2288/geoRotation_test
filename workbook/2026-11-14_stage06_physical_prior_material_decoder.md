# S06：共同有限物理先验、原生材料解码与transport闭合

**材料开启时的坐标非闭合已定位到JSON旋转数组的存储顺序。旧脚本按列展开，实际ACTS按行读取，因而加载了预期旋转的转置。只修正编码后，176次原生查询全部恢复相同bin及原始材料；全部7,680个共同物理点的transport位置、方向和动量闭合至浮点精度，终点电荷差异为0。完整CKF的角度及误判率影响另作全量单因素测量。**

## 阶段目标

推进[Q5纯坐标旋转检查](../roadmap/04_charge_misidentification.md)：使用真实Calypso检查完整有限物理先验在两个坐标表达下的transport；定位先前材料开启时的非闭合，记录实际材料矩阵、bin和slab，更正历史“strict物理严格旋转”的表述。

结果根目录`C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`。旧配置和产物保留原样；本文`legacy`指原先加载的材料编码，`corrected`指本轮按实际decoder修正后的新文件。

## 工作内容与输入定义

输入来自100044、100116、100117各一份RDO的entries `[0,20)`，每种样本按实际记录seed方向的极角排序选择5个代表先验，每个先验抽样512次，共每种2,560点、三种7,680点。RNG为`2026111400+DSID`。抽样保留记录的6×6 covariance全部36项；这些记录的矩阵本身为对角矩阵，没有据此验证任意非零相关矩阵。

共同输入面明确为`z_F=0 mm`：保留seed横向mean及方向、q/p、time和记录的C，不将它声称为原seed surface上的同一Gaussian。每个完整Gaussian只抽样一次，再将每个点用三角函数转成实际方向；两版Calypso接收同一物理position、momentum和q/p的精确旋转表达，不用一阶Jacobian近似有限方向分布。time在静态transport中边缘化；每个点不再携带Gaussian covariance。

| DSID | 先验数×抽样数 | θ抽样穿越0的点数 | q/p与该seed mean反号的点数 |
| --- | ---: | ---: | ---: |
| 100044 | 5×512 | 573 | 840 |
| 100116 | 5×512 | 348 | 67 |
| 100117 | 5×512 | 183 | 455 |

输入电荷取每个draw的q/p符号。因此下文检查输入→终点和坐标表达之间的sign一致性，不以名义PDG统计draw的charge misID。物理点、完整mean/C、选中seed索引和随机种子均冻结在[prior manifest](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/physical_prior_transport_v1/manifest.json)。

实际transport分别到共同输入面、`VetoNu z_F=−3112 mm`及对应primary的首个ST1 SiHit plane。使用Calypso原生geometry、field、extrapolator和材料相互作用；energy loss/scattering均开启。ACTS deterministic material transport没有随机G4散射抽样，本实验不测量CKF measurement update或完整posterior的坐标等价性。

## 实际软件与诊断闭合

原生beforeGeo/master与afterGeo的geometry/field库继续使用已冻结版本；新增observer位于真实Calypso worktree：

```text
/eos/home-x/xcheng/FASER/geoRotation_calypso_material_bins_v4/
  Tracking/Acts/FaserActsKalmanFilter/src/CalypsoSlopeDiagnosticsAlg.cxx
```

实际编译副本为`C/software_material_bins_v4`，各版组合库见`C/physical_prior_material_bins_v4/software_beforeGeo`、`software_afterGeo`。原生`propagate()`仍提供endpoint；另外调用原生`propagationSteps()`记录材料。在实际`BinnedSurfaceMaterial`上读取已解码的4×4 matrix、native local coordinates、`binTriple(hit.position)`和`materialSlab(hit.position)`原始厚度/五个材料参数。不是用Python重算结果代替runtime查询。

增量[0029–0031 patches](../patches/README.md)已在各自正确基底上顺序应用，产出源码与实际编译源码逐字节相同，见[patch验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/verification_v9/patches.json)。六组before/legacy-after的原有18个endpoint字段全部逐值复现此前v3；新增读取没有改变这些transport结果。

## 根因：JSON数组顺序与实际读取规则不一致

F→A的预期旋转为：

```text
R = [[ 0, 0, 1],
     [ 0, 1, 0],
     [-1, 0, 0]]
```

旧`rotate_material.py`及`compare_material_maps.py`按`order="F"`解读/写出数组。ACTS 32.0.2的`Transform3::from_json`按逻辑行执行`rot << rdata[0], ..., rdata[8]`。旧文件实际加载的是`Rᵀ`，将局部横向u反号；原材料数据未随之镜像，查到了另一个bin。CVMFS头文件给出32.0.2及其安装hash，官方[AlgebraJsonConverter.cpp](https://github.com/acts-project/acts/blob/v32.0.2/Plugins/Json/src/AlgebraJsonConverter.cpp)说明读取语义；实际加载对象中的矩阵和查询值才是runtime证明。[源码摘录与获取记录](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/physical_prior_material_bins_v4/source_evidence/sources.json)保存来源。

此前4096点自检在生成和验证两端重复使用同一错误解码假设，只证明脚本内部自洽，没有验证Calypso实际消费后的材料。现已修正两份脚本，并加入3个独立按ACTS逻辑行规则读取的回归检查；历史JSON不覆盖。

用户提供的AFS材料图共同12项也保存了与旧strict相同的旋转数组，因此“AFS与旧strict共同bin数据和旋转编码一致”仍成立，**“已经正确物理旋转”“旧strict是纯旋转控制”撤回**。共同面只有u/v分箱，两文件的±2 mm平移仍只改变法向w，不改变它们彼此的横向查表；但二者共同使用了错误的u符号。旧S04性能数字是实际两配置的测量，保留为历史配置比较，不能解释为严格物理旋转的净效应。

## 详细数据：实际bin与原始材料

从每个DSID的2,560点中，取VetoNu和ST1终点`|Δp/p|`最大的各5点，再加每个seed的第0个draw，去重后共44点。查询覆盖176次材料交互：

| DSID | 交互数 | legacy改变bin | legacy改变原始材料 | corrected改变bin / 原始材料 |
| --- | ---: | ---: | ---: | ---: |
| 100044 | 60 | 60 | 59 | 0 / 0 |
| 100116 | 60 | 60 | 60 | 0 / 0 |
| 100117 | 56 | 56 | 56 | 0 / 0 |

全部paired geometry IDs相同。修正后最大实际物理交点及local-coordinate差均为`2.0605739337×10⁻¹³ mm`；原始slab参数与JSON对应cell逐值相同。

例子：geometry ID=`144115600661151744`，F坐标`(−60.8920685443,−78.8851533336,−1820.1075) mm`。before的local坐标为`(−60.8921,−78.8852,+2) mm`，legacy-after为`(+60.8921,−78.8852,−2) mm`；bin从`[15,22,0]`变为`[44,22,0]`。原始thickness从`811.1181640625`变为`841.239501953125 mm`，X0从`16.2975864410`变为`16.4145965576 mm`。corrected-after恢复原bin和全部材料参数。此为选中先验点在对应材料面查询的raw slab，不当作tracker各层的典型厚度。

## 详细数据：全7,680点transport

三个目标、两种表达下全部点成功。旧编码导致的最大绝对relative-p差如下；修正后两个非平凡目标的relative q/p差均精确为0。

| DSID | legacy VetoNu max\|Δp/p\| | legacy首ST1 max\|Δp/p\| | corrected VetoNu / ST1 max\|Δp/p\| |
| --- | ---: | ---: | ---: |
| 100044 | 0.0020012941 | 0.0005489024 | 4.4409×10⁻¹⁶ / 4.4409×10⁻¹⁶ |
| 100116 | 0.0022792312 | 0.0000213756 | 4.4409×10⁻¹⁶ / 4.4409×10⁻¹⁶ |
| 100117 | 0.0089433635 | 0.0025710953 | 4.4409×10⁻¹⁶ / 4.4409×10⁻¹⁶ |

修正后全目标最大位置差`4.6895820560×10⁻¹³ mm`，最大方向夹角差`2.2074381395×10⁻¹⁶ rad`。包含输入面Cartesian norm舍入后的最大relative-p差为`6.6613381478×10⁻¹⁶`。实际endpoint `charge()`与q/p sign全部一致，before/after终点sign差异为0，输入→终点sign改变为0。

![Native material decoder transport comparison](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v13/material_bins/100117_material_decoder_transport.png)

[3页完整材料解码图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v13/material_bins/physical_prior_transport_report.pdf)同时包含100044、100116、100117。旧map的material-on/off图另保留于[transport图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v11/transport/physical_prior_transport_report.pdf)，用于复核非闭合来源。

[84页S06汇报图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v5/charge_stage06_report.pdf)汇入此前58页、固定seed2页、联合refit18页、legacy transport3页与decoder3页；组装cluster1578725 exit0，全部页的vector content、page boxes与来源单页精确一致，PNG完整性通过。各来源主图/全页contact及本次3张decoder图已有目视验收记录。

## 执行、验证与产物

全部批量transport、归约、绘图和PDF组装经HTCondor执行。每次提交重新查询schedd，按兼容普通bigbird的running+idle最低选择bigbird20，排除已有兼容问题的bigbird15/21。12个v4查询/transport jobs为clusters1578710–1578713，全部exit0，合计CPU163.11 s，最大RSS960,068 KiB。全量归约1578714、最终绘图1578717、图表核对1578721均exit0；最终绘图wall22.54 s、CPU21.27 s、最大RSS122,136 KiB。

复现入口如下，prepare需使用新的run根目录及其中匹配的上游manifest：

```bash
python3 scripts/prepare_charge_prior_transport.py --root <new-root>
python3 scripts/prepare_charge_transport_sign.py --root <new-root>
python3 scripts/prepare_charge_material_trace.py --root <new-root>
python3 scripts/prepare_charge_material_bins.py --root <new-root>
python3 analysis/charge_material_bins_analysis.py --root <root> --output <new-summary>
python3 analysis/charge_prior_transport_analysis.py --phase plot --input <summary.json> --output <new-figures>
```

实际任务、命令、library和resources由各worker的`submission.json`、`source_task.json`、`run.json`保存。主要入口为[查询/transport manifest](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/physical_prior_material_bins_v4/manifest.json)、[完整统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_bins_analysis_v1/summary.json)、[图表验收](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v13/material_bins/validation.json)。54个绘制直方图含flow和sumw2与summary精确相同，3页PDF vector streams/page boxes与单页一致，PNG完整解码；全部3张PNG已全尺寸目视检查。首版legend遮挡underflow的布局已在v13修正，统计未变。

57个unit tests通过，记录为`C/verification_v9/unit_tests.log`，其中材料3项检查独立验证decoder顺序、非identity原transform及null translation。实现为[任务准备](../scripts/prepare_charge_material_bins.py)、[归约](../analysis/charge_material_bins_analysis.py)及[transport分析](../analysis/charge_prior_transport_analysis.py)。

## 阶段性结论

事实：旧map实际加载了转置旋转，改变176次查询中的全部bin；按实际decoder修正后，raw材料及有限物理点transport均闭合。旧4096点自检未覆盖runtime decoder，这一验证缺口已补齐。尚未完成：完整非线性CKF先验及剩余错号measurement机制；transport的零sign变化不代替重建误判率。

## 重建影响控制与下一步计划

真实CKF单因素控制为`C/material_encoding_control_v1`：转置AFS原12个rotation数组，保持IDs、translations、bin axes、data、metadata和surface presence不变，原路径不修改。复用原始runner、worker、实际四个library及所有输入/ranges，isolated smoothing始终ON。

三DSID的legacy/corrected各20事件smoke均成功。legacy replay分别有58/58/56个数组字段与旧完整shard前20事件逐值一致；corrected的原始truth和segments不变，允许truth-model propagation随材料变化。[smoke验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/smoke_validation.json)已通过。μ±各5文件500k、100044单文件500k，共300个5k-event任务已完成，clusters1578722–1578724；[全量覆盖验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/full_validation.json)核对实际库、flags和原始ranges一致。后续[全量统计与12页图验收](2026-11-14_stage06_material_encoding_reconstruction.md)现已完成：μ±仍是同一27/32个wrong，100044为292→291，角峰窗口及核心变化较小；三样本truth-model momentum/slope/validity精确不变。

300任务逐range/库/开关闭合、同母样本的角度/momentum/charge/efficiency及迁移、共同truth与truth-model变化均已验收，详见上述后续报告。继续修正map下的详细audit、Q5完整CKF非线性先验和剩余measurement机制。本文证明材料解码及确定性transport的坐标闭合，S06保持未完成。
