# S05：A1–A5 证据与验收审计

本表逐项核对 [角度路线图](../roadmap/03_slope_structure.md)，测量、输入和范围见[主报告](2026-11-12_stage05_slope_structure.md)及[最终结论](2026-11-13_stage05_angle_structure_conclusions.md)。结果根目录为 `/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes`，下文记为 `S`。beforeGeo主要近零增强已通过联合regular update/smoothing及独立μ± full-G4定位到near-pole方向参数化；afterGeo主要结构具有独立strip控制与预留样本response证据。下表核对的是成因定位与所列窗口的剩余余量；补丁的完整生产修复验收、效率优化和pull校准尚未测量。

验收结论：S05成因定位与A1–A5交付通过。290页统一图集、52组报告、75项证据和逐图来源页码已核对，记录为 `S/final_report_v1/validation.json`。该结论适用于下表明确的输入、cohort与窗口范围。

## A1：物理定义与 writer

| 路线图要求 | 已检查的实际证据与适用范围 | 判定 |
| --- | --- | --- |
| A1.1，F/A 的 slope、符号、逆变换和已知方向 | `main_summary_v1/{beforeGeo,afterGeo}_summary.json` 的 frame/momentum 误差均为0；`controls_analysis_v1/direction/summary.json` 使用294个等价free states/版、实际Calypso传播 | 已测；自由传播闭合不代替近极点拟合验证 |
| A1.2，seed→CKF→conversion→refit→Ntuple 路径 | `native_stage_analysis_v3`、`stage_analysis_v2`、`writer_checks_v2`；真实native writers、afterGeo audit及Ntuple开关对照 | 实际路径已覆盖；seed仅population，不声明同ancestor seed |
| A1.3，st1/VetoNu reference plane | observer实际传播至 `z_F=0/-3112 mm`；geometry probe确认VetoNu两板及间隙；全部原upstream surfaces不同于固定st1 | 已测；VetoNu真值模型与局部G4真值分开 |
| A1.4，精确表达式、单位、binning、分母、collection | `slope_reference.json`、`main_summary_v1`、`figures_v2/*/index.json`、`angle_macro_provenance.json`；原宏实际公式另做literal ROOT核对 | 已测；原作者私有ROOT及完整历史manifest未取得，不声明历史逐bin相等 |
| A1.5，无效值、占位符与失败 | 全量两版propagation success/failed/missing/invalid及writer占位值对照；before/after VetoNu失败5/12、sensor缺truth15/20 | 已测，保留原分母；flow验证见下文 |

## A2：结构参数和选样

| 路线图要求 | 实际证据 | 判定 |
| --- | --- | --- |
| 分箱宽度、半bin移动、ECDF/KDE bandwidth | `main_summary_v1` 的 `histograms/structure`、`figures_v2/100044`；共3280个冻结histogram经 `figure_counts_validation_v1` 核对 | 已测 |
| 精确零、unique spacing、float/default/clamp | 全量正确物理slope均unique、exact zero=0；`zero_analysis_v1`固定25/50/100/200µrad窗口；native raw angles保留float writer值 | 已排除exact-zero堆积；near-zero另有实际chart因果控制 |
| 全谱、reco-p、truth-p、legacy-pz、A′窗口 | 五个DSID的 `figures_v2/*/index.json` 和冻结summary；主100044每版500k，μ±辅助各10k，110033全10k，110051全20k | 已测；不混用各自不同母分母 |
| 同面truth/reco/residual/pull和二维图 | 主summary三面数据及参考图；fixed-plane truth是production-p大小加SiHit chord的Calypso传播模型，sensor方向来自实际SiHit chord | 已测；旧RDO没有实际sensor局部动量 |
| 多峰检验、峰位/宽度/占比/区间和预留事件 | afterGeo discovery冻结窗口；validation窗口占比/条件矩及区间在 `main_summary_v1/*_summary.json`；before另用物理beam-centred窗口 | 已测；窗口条件矩不当作去背景的Gaussian峰参数，检验不单独判物理来源 |

## A3：条件分布与响应重组

| 路线图分解项 | 实际证据与剩余范围 | 判定 |
| --- | --- | --- |
| kinematics、位置、module、hits、cluster width、purity | 主summary的 `groups/maps`；`origin_analysis_v2` 的三站segments、truth效率和response map；`topology_analysis_v2` 的两类A′全量decay z/opening/ST1-centroid separation；110051稀疏opening合并在v3 | 主100044及A′条件分布已测；拓扑是描述性关联，不作为独立消融或局部G4动量真值 |
| 材料/场与动量标度 | 实际G4 world直线材料积分约250.31X0；`controls_analysis_v2/gun` 的μ±×100/250GeV full-G4响应；S02实际field/gradient检查 | 支持VetoNu宽度主要来自靶散射；未将直线积分冒充逐track累计材料 |
| 测量离散性及最终两肩 | `centroid_windows_v1/independent40k` 固定同粒子、同完整hits及native-p窗口11511 pairs；两肩窗口和几何谐波均有配对区间；另有独立G4 segment控制 | 已有独立因果证据；不归因全部肩峰面积 |
| algorithm 与问题交集 | native逐stage、scale1/10/100、周期φ、canonical、spatial basis/plane；`zero_analysis_v1`及实际refit outcomes；`chart_update_analysis_v1`、`chart_gun_analysis_v2`和质量原数组重放 | 联合chart控制RDO近零109→61/1758、truth58；truth-p56→31/1595、truth30。独立μ+1299→765/3906、truth783，μ−1321→753/3912、truth749；所列剩余差区间含0，same-hit结果一致。RDO角误差RMS0.467→0.329mrad、>1mrad38→8。近零全部refit成功，不外推全500k修复效果 |
| 独立response重组、覆盖和系统项 | `phase_response_v2`；`response_systematics_v1` 五种固定分组、discovery学习/validation应用、共同covered支持域与完整uncovered计数 | 已测；afterGeo两肩的固定窗口保守正过量余量为全validation的0.4504%/0.6048%，仅适用于声明过的分组与窗口 |
| response counts、sumw2、flow完整性 | `response_flow_validation_v1/{beforeGeo,afterGeo}/validation.json` 从实际完整输入表重放24组冻结曲线，prediction及observed/covered/uncovered均闭合 | 已测；没有以完整性检查代替物理验收 |

## A4：首次出现位置与最小控制

| 路线图试验 | 实际证据 | 判定 |
| --- | --- | --- |
| truth→segment→seed→CKF P/F/S→conversion→refit→writer/propagation | `stage_analysis_v2`、`native_stage_analysis_v3`、`native_transport_v2`；固定final候选并固定p窗口，两个版本各10k | beforeGeo near-zero在refit F→S增强；同sensor和truth-p控制保留该增强 |
| A4.1，只改变输出解释 | `writer_checks_v2` 实际Ntuple flag false/true，20-event及500-event检查 | 保存的reconstruction字段不变，物理表达式闭合 |
| A4.2，方向稳定性、角度周期与normalization | free-direction scan、refit covariance、ACTS32 periodic/canonical及0015–0017；`chart_math_validation_v3`、`chart_update_refit_v1/smoke_validation.json`；RDO和独立G4四模式 | periodic/canonical及canonical update不足以解释。0017真实53 units编译、关闭逐字段不变及smoother-only复现通过；联合update/smoothing同时改善近零、角误差核心和大尾；pull RMS RDO2.03、独立1.373/1.331，尚未校准 |
| A4.3，strip pitch/stereo/lever arm与digitization | 实际geometry pitch0.08mm、lever63mm、stereo约±0.02；保留完整digitization、只替换locX的真实Calypso诊断 | 实际控制已验证；该操作保留误差/噪声且有fallback，不称为digitization关闭或正式修复 |
| A4.4，material/field | 实际Calypso/G4材料量与100/250GeV响应，旧RDO与新gun严格区别 | VetoNu宽度证据成立；未测的逐track材料或任意关场/关材料结果不作结论 |
| A4.5，选择/refit/covariance/smoothing | 同候选、same-hit、native-p固定；shared-slot、refit covariance、periodic/canonical、spatial basis、input plane及regular chart四模式 | 空间basis109→109、plane109→101；联合chart候选/majority全9707不变、38998 seed/CKF rows不变、76条final hit sets变化；same-hit1744条近零108→61、truth58，RMS0.462→0.282mrad。actual成功/fallback9110/597不变。独立gun same-hit也重复改善 |
| A4.6，geometry/conditions/local basis | production 22.0.49/tag04与运行24.0.41/tag05分开记录；S02 sensor/field闭合与真实geometry probe；module响应保持实际配置 | conditions差异独立记录；结论基于冻结条件中的实际控制 |

## A5：最终交付与证据范围

| 必交付项 | 当前证据 | 验收范围 |
| --- | --- | --- |
| 每个主要结构的来源表、位置/宽度/占比/区间 | 最终结论来源表；主报告五个afterGeo冻结窗口条件矩/区间；联合chart的RDO/gun剩余差及残差核心/大尾 | near-zero取100µrad窗口与声明cohort，RDO+0.1706 [−0.4077,+0.7928]pp、独立μ+−0.4608 [−1.1779,+0.2753]pp/μ−+0.1022 [−0.6480,+0.8544]pp；after两肩保守正过量余量0.4504%/0.6048%。不是任意窄峰的全局排除 |
| 程序问题的最小控制及独立输入 | 实际Calypso原生update/smoother及完整covariance/Jacobian变换；10kRDO与独立μ±各4k；fixed-native/same-hit对照 | 来源定位成立。loader20-event0→20tracks的程序实证及严格成套加载gun已完成，native/smoother-only before字段精确复现；旧after空参考不作物理证据 |
| 同步效率、p-tail、charge、失败/fallback | 全量原生指标、10k控制、独立gun及三面失败；实际refit返回值 | 独立pure-long效率下降2.000/1.925pp，全部候选仍保留且固定native/same-hit角度改善；pull未校准。完整生产修复尚未验收；全500k actual outcome字段未保存，不用10k计数代替 |
| 参考PDF对应主图、legacy/physical、truth/residual、分组、逐stage、放大图 | `figures_v2`五DSID105图；`refit_display_v2`69图/1104hist原843452数值字段精确复现；RDO update/both108hist、gun三模式324hist；最终质量18图/498hist原数组重放通过 | 已按参考页与各自分母交付；新gun包含有效afterGeo own-cohort参考。旧regular及更早nested gun只使用经验证的before控制，不以失效after参考作结论 |
| PNG、逐图PDF、统一汇报PDF、bin counts及index | 52组报告/75项证据；290页统一PDF及290张PNG/单页PDF；冻结counts与独立重放 | `final_report_v1/results/index.json`及`validation.json`：汇编290页decoded streams/text/boxes/rotation精确保留源报告；单页resource/指令及全部bookmarks/index通过；主要原图、质量图及汇编页已目视 |
| 中文报告、可复现配置/源码、阶段commit | 本审计、主报告、最终结论、patch0013–0017和实际worktree、冻结任务/库maps/schedd/资源记录 | 文档与代码同阶段提交，最终检查结果见结论报告 |

数学验证直接使用实际Calypso header、1887组真实writer角度，方向往返6.94e-18、逆Jacobian2.18e-11、差分3.60e-11。0013+0014共同对应basis worktree，0015/0016/0017各自为master-based完整替代patch；apply及源码一致性由 `S/patch_validation_v2/validation.json`核对。后续生产优化使用新的保留输入检验效率与pull；本表的成因验收不将该工作标为已测量。
