# S06：测量位置、seed 曲率与联合角度回归

**afterGeo 原59条错号中，SiHit-centroid测量控制使23条通过选择且恢复正确；只把seed q/p替换为production q/p则没有恢复任何通过选择的正确轨迹。原hit过滤后仍错的6条purity-one轨迹，centroid后5条恢复正确，另1条符号正确但质量失败。** beforeGeo refit的updater-only虽然减少电荷错误，却增加角度大尾；同时改updater和smoother才使本次电荷与角度同时改善。

## 阶段目标

推进[电荷路线图Q4/Q5](../roadmap/04_charge_misidentification.md)：在真实Calypso、isolated smoothing ON下区分hit可用性、测量位置和seed曲率的作用，并给出电荷改善对应的角度、效率和资源变化。完整纯旋转、nominal几何与留出μ±/A′联合验收继续执行；本报告记录本轮已测量的子实验。

结果根目录 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`。

## 工作内容与输入

重放母集合为两电荷各500k原生afterGeo测量中全部59条selected错号（μ+27、μ−32），加原来冻结的60个正常控制（各30）。事件以RDO文件和entry定位；使用原`charge_replay_v1/targets.json`，不重新挑选易恢复事件。matching、候选排序、strict purity>0.7与r95质量条件沿用[S06基线](2026-11-13_stage06_charge_baseline.md)。这是一组针对已知outcome的控制，所列比例不能当作新配置的全样本误判率。

所有tracking为真实Calypso CKF/refit。hit控制库来自`geoRotation_calypso_charge_hitcontrol`，冻结在`C/overlay_hitcontrol_v2`；centroid overlay额外使用S05已验证的实际`FaserSiClusterizationTool`两库，位于`stage05_slopes/software_cluster_control_v2`。未另写tracking或简化detector/field。

| 控制 | 实际变化 | 保持条件与解释 |
| --- | --- | --- |
| primary hits | CKF只接收primary-SDO支持的source links | 原seed六维均值与完整C不变；分支及winner可以变化 |
| SiHit centroid | 实际cluster maker将measured locX替换为同wafer majority-SiHit的能量加权midpoint | 保留width/error、noise、IDs及缺truth fallback；会影响segment/seed，不称digitization off |
| centroid + primary hits | 同时施加上述两项 | 用于检查测量位置与hit可用性的交互 |
| production seed q/p | 实际CKF `SeedOverrideMode=qop`仅替换seed第4分量为primary生产点charge/p | 其他五个均值、完整C、surface及seed hits逐项不变；production p不是sensor局部G4动量 |

centroid三类各20事件smoke的default58字段与基线完全相同，实际cluster两库和替换计数通过。q/p控制三类各20事件的identical override同样逐字段复现基线；production override只改变q/p。全部控制使用实际原生finding与refit，保留质量失败及purity失败。

## 详细数据：59条原错号与60条正常控制

表中“正确/错误”均要求通过matching、quality、purity及有效曲率选择；quality或purity失败单列，不计为恢复。全部为源事件原始计数。

| 原cohort | 控制 | correct | wrong | quality_failed | low_purity |
| --- | --- | ---: | ---: | ---: | ---: |
| 59 wrong | native | 0 | 59 | 0 | 0 |
| 59 wrong | primary hits | 17 | 6 | 36 | 0 |
| 59 wrong | centroid | 23 | 30 | 6 | 0 |
| 59 wrong | centroid + primary hits | 28 | 0 | 31 | 0 |
| 59 wrong | production seed q/p | 0 | 57 | 1 | 1 |
| 60 normal | native | 60 | 0 | 0 | 0 |
| 60 normal | primary hits | 59 | 0 | 1 | 0 |
| 60 normal | centroid | 60 | 0 | 0 | 0 |
| 60 normal | centroid + primary hits | 59 | 0 | 1 | 0 |
| 60 normal | production seed q/p | 60 | 0 | 0 | 0 |

按电荷拆分，centroid为μ+9 correct/15 wrong/3 quality_failed、μ−14/15/3；centroid+primary为μ+12 correct/15 quality_failed、μ−16/16。production q/p为μ+26 wrong/1 quality_failed、μ−31 wrong/1 low_purity。119个q/p作业的126个seeds全部核对其他mean/C/hits/surface不变。

centroid在共同seed ID/hits/surface上改变所有可配对seed均值，完整C保持相同。59个wrong events中，μ+有1个seed hit set被替换，μ−移除5个、增加4个；正常60个events的seed hit sets全部保留。因而“centroid恢复23条”同时包含测量位置及其引起的seed/branch响应，尚不能全部分配给某一次measurement update。

三页新图：outcome counts、六条pure cases的final signed q/p、共同measurement最大位移。首轮绘图的PNG/PDF已经产生，但JSON归约遇到NumPy int64不能序列化；仅将flow计数转为Python int，冻结新revision重新导出，未改重建数据或统计定义。最终v3仅改pure页坐标留白与legend位置，完整数值与v2精确相同，68项summary计数、3页vector及PNG核对通过，三图已目视检查。

![Measurement and seed outcomes](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v7/measurements/measurement_control_outcomes.png)

[六条pure cases](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v7/measurements/measurement_pure_wrong_cases.png) · [测量位移](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v7/measurements/measurement_centroid_shifts.png) · [3页PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v7/measurements/charge_measurement_report.pdf)。

## 详细数据：六条primary过滤后仍错的purity-one轨迹

单位均为GeV⁻¹。列中的q/p保留实际符号，μ+正确为正、μ−正确为负。centroid+primary与centroid的这六个结果相同。

| DSID / sample entry | native final q/p | centroid final q/p | centroid outcome | production-seed final q/p |
| --- | ---: | ---: | --- | ---: |
| 100116 / 51839 | −0.001534125092 | +0.001992525791 | correct | −0.001533917591 |
| 100116 / 224394 | −0.001469602076 | +0.001194047528 | correct | −0.001469617746 |
| 100117 / 227965 | +0.000264475596 | −0.002086818750 | correct | +0.000264474376 |
| 100117 / 390966 | +0.000407198602 | −0.001873240217 | correct | +0.000407177579 |
| 100117 / 439421 | +0.002789309002 | −0.002172017002 | quality_failed | +0.002806015411 |
| 100117 / 458033 | +0.000327301885 | −0.001665168309 | correct | +0.000327302527 |

这六个events的centroid seed数量、ID、hits、surface和完整C均与native相同，但seed mean改变。仅改seed q/p全部仍错；改变测量位置能使全部六条最终符号正确，但只有五条通过质量条件。该证据排除了“只要给对seed p就可恢复这六条”的解释，支持测量位置及其传播响应的作用。

进一步的66事件控制已完成：上述6条加原60条正常控制，centroid测量仍开启，同时通过实际 `SeedOverrideMode=mean`恢复原生六维seed mean。完整数值mean、C和surface逐项精确，hit多重集合相同；16个seed的hit vector顺序改变，未称逐元素顺序相同。六条pure仍为5条selected正确、1条符号正确但rmax=95.0914 mm未通过95 mm cut。原60条正常控制中59条正确、1条失去下游两站而质量失败。完整结果与逐hit路径见[固定完整seed报告](2026-11-14_stage06_fixed_seed_alignment.md)。

## 详细数据：beforeGeo refit的联合角度回归

复用真实100k/电荷/模式的native、regular smoother、regular updater、both四组输出。所有CKF记录逐项相同。实际Calypso将track投影到`z_F=0 mm`；truth为首st1 SiHit chord附上production p后传播的模型，和独立gun的同面实际launch truth分别解释。

| charge | 模式 | |Δsy|的未分箱68%分位数 / μrad | |Δsy|>1 mrad / finite selected |
| --- | --- | ---: | ---: |
| μ+ | native | 74.6225 | 403 / 50335 |
| μ+ | smoother | 67.3225 | 410 / 50342 |
| μ+ | updater | 71.1521 | 913 / 50289 |
| μ+ | both | 47.8394 | 48 / 50318 |
| μ− | native | 76.7799 | 499 / 52030 |
| μ− | smoother | 69.6590 | 615 / 52031 |
| μ− | updater | 73.0600 | 1116 / 51928 |
| μ− | both | 43.6194 | 50 / 51999 |

**updater-only的charge改善有角度大尾代价**，不能只看wrong-sign数验收。both在相同candidate ordinal和ordered hits的共同selected集合中也改善：μ+50,096条的core74.3422→47.8388 μrad、1mrad尾385→46；μ−51,734条为76.3211→43.6447 μrad、478→50。共同集合排除了仅因换候选或selection迁移导致的表观改善。

![Joint angle regression](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v5/regression/refit_angle_core_and_tail.png)

[共同candidate角度图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v5/regression/refit_common_candidate_angles.png) · [CPU/RSS图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v5/regression/refit_cpu_and_memory.png)。

GNU time的user+system覆盖完整tracking进程；按source events归约CPU，单位ms/event：μ+ native/smoother/updater/both为11.5590/11.2271/12.7739/10.9547，μ−为11.2561/12.5677/14.2231/13.9982。最大RSS分别约989–991 MiB与1377–1391 MiB，原始KiB保存在`resource_analysis_v1/summary.json`。不同worker/batch，这些是实测资源描述，不给严格加速比。`acc.run` wall包含Athena初始化与结束，排除之前的Python配置阶段。

## seed-hit分类覆盖修正

旧归约把“没有出现在任何candidate measurement audit中的seed hit”归入`no_SDO`。现在单列`not_in_candidate_measurement_audit`，只有确实审计且SDO为空才记`no_SDO`。119事件中34个seed分类改变；summary、pruning及除该字段外的全部event记录与v1精确相同，56次最终CKF翻转的31 no-SDO/25 primary-only结论不变。

六条pure中的224394、227965、390966各有5个seed hits缺measurement audit，458033有13个；不能把它们解释为已证实的seed污染。修正与逐项相等检查保存在`C/replay_analysis_v2/classification_validation.json`。

## 执行与验证

centroid实际重放cluster1577639的238/238 jobs exit0；归约1577949 exit0，实际cluster两库全部匹配。归约wall143.24 s、CPU24.97 s、最大RSS75,744 KiB。production seed q/p实际重放cluster1578082的119/119 exit0；归约1578083 exit0，wall44.19 s、CPU11.77 s、最大RSS76,288 KiB。refit regression、资源归约1577643/1577645以及3页绘图1578084均exit0。

新提交时bigbird15的running+idle短暂低于bigbird20，但该schedd实际拒绝submit文件中的EOS路径；两次提交均未进入queue。保存原提交文件与日志后排除bigbird15及已知不兼容的bigbird21，重新比较负载，选bigbird20，恢复提交成功。失败来自提交环境，不影响已完成的tracking数据。

分析与准备入口：[charge_centroid_analysis.py](../analysis/charge_centroid_analysis.py)、[charge_seed_qop_analysis.py](../analysis/charge_seed_qop_analysis.py)、[prepare_charge_seed_qop.py](../scripts/prepare_charge_seed_qop.py)、[charge_refit_regression.py](../analysis/charge_refit_regression.py)、[charge_resource_analysis.py](../analysis/charge_resource_analysis.py)、[charge_measurement_figures.py](../analysis/charge_measurement_figures.py)。冻结任务、软件与source manifests分别在对应`*_tasks.json`、`*_code`和运行目录。

```bash
python3 scripts/prepare_charge_seed_qop.py --root "$C" --phase replay
python3 scripts/submit_momentum.py --tasks "$C/seed_qop_analysis_v1_tasks.json" \
  --library "$C/overlay_hitcontrol_v2/lib" --name <unique_name> --stage S06 \
  --cpus 1 --memory-mb 3000 --disk-gb 2 --flavour microcentury \
  --exclude-schedd bigbird21.cern.ch --exclude-schedd bigbird15.cern.ch
```

[centroid summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/centroid_analysis_v1/summary.json) · [seed q/p summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/seed_qop_analysis_v1/summary.json) · [联合角度summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/refit_regression_v1/summary.json)。

[58页完整S06图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v4/charge_stage06_report.pdf)由原36页、机制5页、gun11页、角度/资源3页及measurement3页组成。cluster1578091 exit0，58页vector内容及page boxes逐页相同，三张contact sheets已目视检查；完整索引与validation在同目录。fixed-seed新增图另列，不改变这份58页历史图集。

## 阶段性结论

事实：原afterGeo错号对cluster坐标控制的响应明显，对production seed q/p-only的响应不足以恢复任何selected正确轨迹。hit过滤与centroid的作用有交互，全部原59条在两项共同控制后无selected wrong，但31条质量失败，不能当作“59条全部修复”。正常控制也保留1例质量损失。before refit updater与smoother需要联合看角度和电荷，updater-only不是可接受的联合改善。

事实补充：恢复原seed完整均值后，五条selected恢复仍存在，六条pure的符号均正确；其中两条selected ordered hits与native相同。这个结果说明seed重估不是这些恢复的必要条件；其余四条仍包含measurement引起的branch/path变化。measurement error模型、strip量化及沿程散射的各自贡献与完整共同物理先验的坐标效应继续判别。生产修复的全效率、pull及留出A′回归按Q5验收。

## 下一步计划

1. 固定完整seed已完成；继续对固定seed及相同ordered hits的事件检查measurement innovation与curvature响应，保持质量损失单列。
2. 真实nominal conditions库与6组兼容验证已完成，默认58字段精确；全部192个sensor的conditions transform等于nominal。匹配G4/digitization/reco的小样本链已提交，surface与全表闭合后再判断是否需要统计扩展。
3. 完整共同物理先验、空间适用性和留出μ±/A′的动量/角度/电荷/效率联合回归仍须完成，不将一阶J C Jᵀ称为全链纯旋转。
