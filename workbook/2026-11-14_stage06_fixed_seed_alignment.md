# S06：固定完整 seed 的测量控制与 nominal alignment

**六条原 purity-one 错号在恢复完整 native seed 后，仍有五条通过选择且符号恢复正确；第六条符号正确，但 rmax=95.0914 mm 超过95 mm cut。恢复不需要重估 seed，测量坐标及其引起的 CKF 响应已足以改变这些结果。当前 MC conditions 下两版全部192个sensor的实际transform与nominal完全相同；匹配G4、digitization与重建链通过，default/nominal六组58字段也逐值相同。**

## 阶段目标

推进[Q4/Q5](../roadmap/04_charge_misidentification.md)：区分SiHit-centroid控制引起的seed重估和测量/branch响应；解释quality失败；在实际Calypso中验证nominal alignment的conditions链与sensor transforms。Q5完整共同物理先验及留出μ±/A′联合回归继续执行，不将本轮控制等同于生产修复。

结果根目录`C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`。所有重建、simulation和digitization均为真实Calypso组件，isolated smoothing ON。

## 工作内容与控制定义

固定集合为[measurement/seed报告](2026-11-13_stage06_measurement_seed_controls.md)中六条原native final purity=1、primary-SDO过滤后仍错的轨迹，加最初冻结的60个正常事件；没有根据本控制的outcome重新挑选。输入以RDO路径和entry定位，μ+为100116、μ−为100117。原500k/电荷的全部59个selected错号不等于本轮六条pure子集。

真实cluster maker继续开启`UseTruthClusterPosition`：locX改为同wafer的majority-SiHit能量加权midpoint，保留measurement error、width、noise、IDs及缺truth fallback。实际CKF的`SeedOverrideMode=mean`将六维mean恢复为原native，完整C和surface逐项精确，seed hit多重集合相同。66个seeds中16个hit vector顺序改变，单独保存；没有声称其顺序逐元素不变。真实CKF仍可改变branch及所选measurement路径。

完整mean的identical override已用100044/100116/100117各20事件验证：58字段精确，20/17/18个seed的mean/C/surface/hits相同，记录在`C/seed_mean_smoke_v1/validation.json`。固定seed重放66/66 jobs exit0，实际两cluster库与冻结版本全匹配。

## 详细数据：固定完整 seed 后的结果

correct和wrong要求通过原charge-independent matching、purity>0.7和Good95条件；符号正确而quality失败仍单列。表中是冻结目标集合的事件计数，不能当作总体misID率。

| 集合 | charge | native correct / wrong / quality_failed | centroid、自由seed | centroid、完整native seed |
| --- | --- | --- | --- | --- |
| 原pure wrong | μ+，2事件 | 0 / 2 / 0 | 2 / 0 / 0 | 2 / 0 / 0 |
| 原pure wrong | μ−，4事件 | 0 / 4 / 0 | 3 / 0 / 1 | 3 / 0 / 1 |
| 原normal | μ+，30事件 | 30 / 0 / 0 | 30 / 0 / 0 | 29 / 0 / 1 |
| 原normal | μ−，30事件 | 30 / 0 / 0 | 30 / 0 / 0 | 30 / 0 / 0 |

![Fixed native seed outcomes](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v8/fixed_seed/fixed_seed_outcomes.png)

六条pure的signed q/p如下，单位GeV⁻¹。charge定义为μ+正、μ−负。

| DSID / sample entry | native q/p | centroid自由seed q/p | centroid固定seed q/p | 固定seed选择与路径 |
| --- | ---: | ---: | ---: | --- |
| 100116 / 51839 | −0.001534125092 | +0.001992525791 | +0.001885651571 | correct；ordered selected hits相同 |
| 100116 / 224394 | −0.001469602076 | +0.001194047528 | +0.001193964721 | correct；selected路径改变 |
| 100117 / 227965 | +0.000264475596 | −0.002086818750 | −0.002086832738 | correct；selected路径改变 |
| 100117 / 390966 | +0.000407198602 | −0.001873240217 | −0.001873238149 | correct；selected路径改变 |
| 100117 / 439421 | +0.002789309002 | −0.002172017002 | −0.002170998023 | quality_failed；ordered selected hits相同 |
| 100117 / 458033 | +0.000327301885 | −0.001665168309 | −0.001665151196 | correct；selected路径改变 |

![Fixed seed signed curvature](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v8/fixed_seed/fixed_seed_pure_curvature.png)

439421只有半径条件失败：rmax=95.0914009080 mm、chi2=0.1101484373、ndof=13、18 hits、9层、station_bits=14、purity=1。没有把略超cut的轨迹记为恢复成功，也没有改cut获得成功。

正常控制中100116/109670只剩station1：station_bits=2、6 hits、3层、ndof=1、chi2=0.7100260，q/p仍为正。失败条件是三下游站覆盖、至少7层和ndof+5≥12。它在自由centroid seed下通过，说明固定旧seed与新测量之间也会造成质量损失。

共同已审计measurement的最大位移为0.04538–0.27937 mm；这是本控制的实际locX差，不能直接把它分配给单一strip量化或alignment原因。六条pure的每个measurement update均验证`Δ(q/p)=(K r)_{q/p}`，最大绝对恒等式误差3.925×10⁻¹⁶ GeV⁻¹。

## 逐运算证据与解释范围

51839的完整seed、seed hit顺序及最终selected ordered hits均相同，因此它提供了measurement响应的直接对照。native最后进入错号的更新位于S3 L0、z_F=2389.5 mm：predicted q/p=+5.66006×10⁻⁵、residual=+1.18814 mm、K_q/p=−0.00224278 GeV⁻¹/mm，filtered q/p=−0.00260815。固定seed centroid同hit的predicted q/p=+0.00406765、residual=+1.59703 mm、K_q/p=−0.00224787，filtered仍为+0.000477725。不是简单地要求每次residual都减小；前序测量改变了进入该面的曲率估计。

439421也保持相同ordered hits。native在S2 L0、z_F=1211.415 mm从−0.00115469更新为+0.00944911；centroid固定seed在同hit由−0.00170806更新为−0.000138221，符号保持正确。其最终淘汰来自半径cut。

其余四条的selected路径改变，不能声称固定了全部CKF hits；这些结果证明seed重估不是恢复的必要条件，同时包含measurement引起的branch响应。只替换production seed q/p的旧控制则六条全仍错，见上述measurement/seed报告。measurement模型、局部G4轨迹及散射的进一步拆分仍待完成。

## nominal alignment：实际源码、构建与初步实测

Calypso现有`flags.GeoModel.Align.Disable=True`会跳过SCT alignment和detector-element conditions algorithms；只关该flag不能建立可用的nominal reconstruction链。因此在实际`FaserSCT_AlignCondAlg`中加入默认关闭的`UseNominalAlignment`：不读取alignment folder、不调用conditions displacement，仍生成真实`SCTAlignmentStore`并沿原detector-element/ACTS conditions链消费。另有只读audit保存每个sensor的实际/default 4×4矩阵。

- 实际源码：`/eos/home-x/xcheng/FASER/geoRotation_calypso_nominal_alignment/Tracker/TrackerConditions/FaserSCT_ConditionsAlgorithms/src/FaserSCT_AlignCondAlg.{h,cxx}`，基底`88c457c1`。
- 补丁：[0028](../patches/0028-nominal-alignment-store-control.patch)。该conditions package在master`1eb885ae`和旋转基底的源码diff为空，复用同一package控制库；其余before/after几何、tracking库分别保持原生来源。
- 编译入口：[build_calypso_alignment.py](../scripts/build_calypso_alignment.py)，复用既有CMake compile/link recipes，在Condor scratch重编11个translation units，用Athena genconf/listcomponents生成匹配metadata；冻结库在`C/software_nominal_alignment_v1`。
- 配置和验证：[calypso_alignment.py](../scripts/calypso_alignment.py)、[prepare_charge_alignment.py](../scripts/prepare_charge_alignment.py)、[validate_charge_alignment.py](../scripts/validate_charge_alignment.py)。

6组兼容验证为before/after×100044/100116/100117，各20事件。全部58字段与原native数据精确相同，实际条件库加载匹配。每组192个sensor的actual/default矩阵差精确为0，六组sensor矩阵也逐元素相同。这是在本次OFLCOND-FASER-05 MC条件及所查IOV的实测；不能从flag开启推断有非零misalignment。

匹配full-G4/digitization/reco的4个20-event控制全部完成：两电荷各default与nominal，三个重建模式为beforeGeo、beforeGeo_regular_both与afterGeo，HITS/RDO仅在scratch。首轮两个nominal作业完成G4事件后，审计在application teardown后读取已释放GeoModel而失败；将只读审计移到首事件后，仅恢复这两个terminal jobs。default两项复用原成功结果。

最终逐元素验证通过：两电荷的production位置与动量相同，三种重建的truth表相同；nominal simulation、digitization与两版重建的192个sensor实际/default矩阵差均为0，跨步骤实际矩阵在1e−9绝对容差内相同。default重建六组58字段与原symmetric-gun smoke精确一致；同charge、同模式default/nominal的58字段也全部逐值相同。在该OFLCOND-FASER-05与所查IOV下，nominal控制就是同一实际几何，alignment位移不能解释本次剩余错号；这项等价性控制不需扩大至80k重新测同一分布。

## 执行、验证与产物

fixed-seed真实重放cluster1578090共66作业全部exit0；归约v3为1578174，wall24.95 s、CPU7.90 s、最大RSS77,712 KiB。首轮归约1578170遗漏冻结依赖，未进入数据处理；v2的1578173因把hit vector顺序当成hit集合而assert失败。核对六条原audit后改为hit多重集合精确匹配，数值mean/C/surface仍保持严格逐元素检查；保留两次失败及冻结代码，未重跑tracking。

两页新图cluster1578179 exit0，10项源表/图计数核对、2页PDF vector/boxes与PNG完整性通过，两张分别目视检查。源summary与逐运算表在`C/fixed_seed_analysis_v3`，图counts与validation在`C/figures_v8/fixed_seed`。

nominal库编译cluster1578172 exit0，wall91.30 s、CPU239.18 s、最大RSS1,018,932 KiB，请求4CPU/14GB；兼容控制1578175/1578176共6作业exit0。gun首轮1578177中两个default成功；恢复1578180两个nominal均exit0，原失败目录保留。四个成功full-chain jobs的wall合计1164.12 s、CPU合计746.50 s、最大RSS1,857,812 KiB，另两次失败的消耗未计入该成功合计。每次提交重新查询schedd广告及queue，选择兼容且running+idle最少的bigbird20；排除实际拒绝EOS提交路径的bigbird15及已知不兼容的bigbird21。其余请求1CPU，tracking4–4.5GB，归约3GB，BLAS/OpenMP线程与资源一致。

本轮51项单测全部通过（2.027s），新增检查区分“符号正确但半径失败”和“失去下游站覆盖”；日志在`C/verification_v6/unit_tests.log`。Python编译与文档链接在本阶段提交前核对，结果保存在同verification目录。

```bash
python3 scripts/prepare_charge_seed_qop.py --root "$C" --phase mean-smoke
python3 scripts/prepare_charge_seed_qop.py --root "$C" --phase centroid-fixed
python3 scripts/validate_charge_alignment.py --root "$C" --phase compatibility
python3 scripts/validate_charge_alignment.py --root "$C" --phase gun-smoke
python3 scripts/prepare_charge_alignment.py --root "$C" --phase gun-recovery
```

以上prepare命令创建唯一任务目录，复现需使用新的run目录并按冻结manifest恢复上游输入，不覆盖当前已完成结果。实际提交使用`submit_momentum.py --tasks <frozen tasks.json> --library <matching overlay/lib> --stage S06`，完整命令、库来源、资源及schedd选择随各worker记录保存。

[fixed-seed summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/fixed_seed_analysis_v3/summary.json) · [逐事件/逐update结果](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/fixed_seed_analysis_v3/events.json) · [2页PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v8/fixed_seed/fixed_seed_report.pdf) · [alignment兼容验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/alignment_compatibility_v1/validation.json)。此前[58页完整S06图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v4/charge_stage06_report.pdf)仍有效。

## 阶段性结论

事实：六条pure的measurement控制在完整seed固定后仍消除最终错误符号，其中五条满足原选择。仅改seed q/p不能达到这一结果。两条同ordered hits的逐update证据直接显示measurement改变前序曲率及穿零响应；四条包含branch/path变化。正常控制有一例质量损失，不能把控制当作无代价修复。

事实：当前MC conditions的192个sensor实际alignment为nominal；条件库默认行为和只读audit不改变本次58字段。匹配G4/digitization/reconstruction链通过，六组default/nominal重建表完全相同。[完整链验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/alignment_gun_smoke_v1/validation.json)保存了逐模式几何及字段检查。假设：measurement error、数字化坐标和沿程散射共同影响剩余错误；本控制尚未给出三者的独立比例。

## 下一步计划

1. nominal完整链等价性已完成；后续更换conditions或IOV时重新核对实际sensor transforms。
2. 对measurement响应保留同完整seed/相同ordered hits控制；区分SiHit-centroid位置定义、量化与误差模型，不将truth辅助控制部署为生产重建。
3. 完整共同物理先验和留出μ±/A′联合回归继续按Q5执行，覆盖动量、角度、电荷、效率、fallback与CPU；S06保持未完成。
