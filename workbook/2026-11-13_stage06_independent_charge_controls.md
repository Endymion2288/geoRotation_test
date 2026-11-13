# S06：独立 full-G4 电荷控制与模式回退计数

**相同生成条件的1 TeV μ± 中，beforeGeo 原生误判率为12.5585%/13.1033%；仅启用 regular refit updater+smoother 后为0.1392%/0.1783%，afterGeo为0.0271%/0.0162%。角度核心和大尾同步改善。100 GeV 的 before regular 控制损失1.755/1.900个百分点正确效率，必须与误判率一起判断。各组μ−−μ+的95%区间均包含0。**

## 阶段目标

用独立 full-Geant4 μ±、相同实际 production 条件检验 near-pole refit 的电荷归因，完成 Q3/Q5 的新样本与真实空间分组；补齐24个实际 direction/station/merge 配置的全输出与 selected refit fallback。完整纯坐标先验、剩余纯轨迹机制及μ±/A′联合回归仍待完成，本报告不将这些子项记为通过。

## 工作内容

结果根目录 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`。真实 Calypso 的 full-G4→digitization→native reconstruction 共32个作业、80,000个G4源事件及240,000次重建事件，全部步骤exit0。HITS/RDO只在worker scratch；保留紧凑表、实际日志、配置、库加载和统计。

每个charge每个100/1000 GeV点20,000事件：launch在`z_F=0 mm`，x/y均匀±30 mm，两个方向斜率范围±0.5 mrad。μ±同shard使用相同gun seed；归约逐元素验证全部实际production位置/动量相同，三个重建版本的truth表相同。使用真实FASERNU-04 apparatus、磁场和材料，不假定正负粒子经历相同随机相互作用，也不假定真实探测器精确镜像对称。请求nominal p与保存的实际p最大相对差约3.9×10⁻⁸；residual使用实际保存的production p。

三个真实重建为beforeGeo native、beforeGeo仅refit updater+smoother使用regular chart、afterGeo native，全部isolated smoothing ON。before基底master 1eb885ae+0022+0027，库`C/software_refit_before`；after基底88c457c1+0021，库`C/software_afterGeo`，observer为`C/observer`。实际加载四个相关库全量核对一致。seed、finding、geometry、material和前向CKF在before native与regular控制间不同时修改。主FLUKA控制在[逐运算报告](2026-11-13_stage06_charge_controls.md)。

## 详细数据：独立 gun 主率

每行源事件与fiducial truth分母均为20,000。conditional分母为相应模式选中的matched track；95%使用exact binomial。

| p / GeV | charge | 配置 | wrong / selected | misID / % [95%] | correct / fiducial / % |
| ---: | --- | --- | ---: | --- | ---: |
| 100 | μ+ | beforeGeo | 0 / 19595 | 0.000000 [0.000000,0.018824] | 97.975 |
| 100 | μ+ | beforeGeo_regular_both | 1 / 19245 | 0.005196 [0.000132,0.028948] | 96.220 |
| 100 | μ+ | afterGeo | 0 / 18564 | 0.000000 [0.000000,0.019869] | 92.820 |
| 1000 | μ+ | beforeGeo | 2255 / 17956 | 12.558476 [12.077086,13.052028] | 78.505 |
| 1000 | μ+ | beforeGeo_regular_both | 25 / 17960 | 0.139198 [0.090101,0.205416] | 89.675 |
| 1000 | μ+ | afterGeo | 5 / 18431 | 0.027128 [0.008809,0.063297] | 92.130 |
| 100 | μ− | beforeGeo | 1 / 19633 | 0.005093 [0.000129,0.028376] | 98.160 |
| 100 | μ− | beforeGeo_regular_both | 1 / 19253 | 0.005194 [0.000132,0.028936] | 96.260 |
| 100 | μ− | afterGeo | 1 / 18566 | 0.005386 [0.000136,0.030006] | 92.825 |
| 1000 | μ− | beforeGeo | 2351 / 17942 | 13.103333 [12.612807,13.605851] | 77.955 |
| 1000 | μ− | beforeGeo_regular_both | 32 / 17947 | 0.178303 [0.121990,0.251618] | 89.575 |
| 1000 | μ− | afterGeo | 3 / 18480 | 0.016234 [0.003348,0.047435] | 92.385 |

100 GeV的原生before几乎没有错号，zero count保留非零上限；regular控制没有可测的电荷改善，却降低正确效率。因此不把regular chart称为已通过生产验收。1 TeV两电荷的错误显著下降且正确效率从78.505%/77.955%提高到89.675%/89.575%；after为92.130%/92.385%。

## 详细数据：相同生成条件后的正负差

下表区间为两份同时97.5% Clopper–Pearson区间构成的保守95%差值界，保留两个zero-count样本的不确定度；不将独立μ±的重建随机过程按event号强行配对。

| p / GeV | 配置 | μ−−μ+ / 百分点 [95%] |
| ---: | --- | --- |
| 100 | beforeGeo | +0.005093 [-0.022296,+0.032497] |
| 100 | beforeGeo_regular_both | -0.000002 [-0.033086,+0.033073] |
| 100 | afterGeo | +0.005386 [-0.023534,+0.034364] |
| 1000 | beforeGeo | +0.544857 [-0.579731,+1.669229] |
| 1000 | beforeGeo_regular_both | +0.039105 [-0.100371,+0.178448] |
| 1000 | afterGeo | -0.010894 [-0.066679,+0.045347] |

事实：这些生成条件下未检出剩余charge差；1 TeV after差值95%为[−0.066679,+0.045347]个百分点。区间描述本gun及统计精度，不等同于主FLUKA谱下的界，也不证明真实apparatus严格对称。

## 详细数据：同面角度、pull与资源

实际Calypso将track传播到launch平面`z_F=0 mm`，以该面实际G4生成方向为truth。此gun没有production→st1的上游靶材传播；与主样本的SiHit chord truth model严格区分。角度细图为5 μrad/bin；下表是未分箱的绝对残差68%分位数及|Δsy|>1 mrad原始计数。

| p / GeV | charge | 配置 | |Δsy| 68% / μrad | >1 mrad / finite selected | 实测reco wall / ms·event⁻¹ |
| ---: | --- | --- | ---: | ---: | ---: |
| 100 | μ+ | beforeGeo | 191.316 | 242 / 19595 | 9.939 |
| 100 | μ+ | beforeGeo_regular_both | 29.682 | 2 / 19245 | 9.442 |
| 100 | μ+ | afterGeo | 29.476 | 0 / 18564 | 14.024 |
| 1000 | μ+ | beforeGeo | 177.121 | 536 / 17956 | 14.014 |
| 1000 | μ+ | beforeGeo_regular_both | 20.918 | 5 / 17960 | 12.274 |
| 1000 | μ+ | afterGeo | 20.472 | 1 / 18431 | 12.326 |
| 100 | μ− | beforeGeo | 187.153 | 197 / 19633 | 13.393 |
| 100 | μ− | beforeGeo_regular_both | 30.306 | 1 / 19253 | 13.122 |
| 100 | μ− | afterGeo | 30.026 | 3 / 18566 | 14.919 |
| 1000 | μ− | beforeGeo | 182.978 | 511 / 17942 | 14.556 |
| 1000 | μ− | beforeGeo_regular_both | 20.957 | 10 / 17947 | 14.800 |
| 1000 | μ− | afterGeo | 20.595 | 0 / 18480 | 14.894 |

before原生角误差核心约177–191 μrad，regular后约21–30 μrad，after约20–30 μrad。大角尾同样减少，支持方向chart的联合机制；formal pull仍需误差模型校准。wall包含各自worker与启动开销，不能据此宣称严格CPU加速比。

每个launch x/y象限的misID、效率及其exact区间保存在`C/gun_analysis_v2/summary.json`的`production_regions`。分组检验实际apparatus，不用假定的镜像关系消去实际场、材料及接受度效应。

## 详细数据：24组 native refit fallback

100k/电荷/模式的400份只读counter重放中，391份首次成功，9份因EOS库/NPZ读写失败。恢复仅重跑4份实际未完成重建，另外5份复用已完成表在scratch重新验证；原失败日志和退出码保留。全部400份58字段与原生baseline精确相同。每个native counter的source entry、ordinal和q/p与独立方向track逐值闭合，merge内部两方向status与独立方向一致；输出用完整hits/qop/upstream参数溯源。

| charge | 版本 | 模式 | fallback / 全输出 | fallback / selected |
| --- | --- | --- | ---: | ---: |
| μ+ | beforeGeo | 3forward | 14224 / 90890 | 0 / 50335 |
| μ+ | beforeGeo | 3backward | 15280 / 89139 | 0 / 51052 |
| μ+ | beforeGeo | 3merge | 25278 / 103839 | 0 / 51835 |
| μ+ | beforeGeo | 4forward | 19775 / 97716 | 9 / 46251 |
| μ+ | beforeGeo | 4backward | 16466 / 92014 | 16 / 50027 |
| μ+ | beforeGeo | 4merge | 33012 / 114052 | 21 / 50787 |
| μ+ | afterGeo | 3forward | 14225 / 90893 | 0 / 49468 |
| μ+ | afterGeo | 3backward | 14575 / 89135 | 0 / 49906 |
| μ+ | afterGeo | 3merge | 24455 / 104532 | 0 / 51760 |
| μ+ | afterGeo | 4forward | 27260 / 97689 | 3 / 43066 |
| μ+ | afterGeo | 4backward | 15415 / 91962 | 13 / 48947 |
| μ+ | afterGeo | 4merge | 39891 / 121530 | 14 / 50573 |
| μ− | beforeGeo | 3forward | 13492 / 90771 | 0 / 52030 |
| μ− | beforeGeo | 3backward | 14510 / 89201 | 0 / 52790 |
| μ− | beforeGeo | 3merge | 24153 / 103256 | 0 / 53582 |
| μ− | beforeGeo | 4forward | 19045 / 97464 | 13 / 47841 |
| μ− | beforeGeo | 4backward | 15572 / 91777 | 13 / 51767 |
| μ− | beforeGeo | 4merge | 31575 / 112892 | 22 / 52594 |
| μ− | afterGeo | 3forward | 13551 / 90775 | 0 / 51197 |
| μ− | afterGeo | 3backward | 13636 / 89199 | 0 / 51784 |
| μ− | afterGeo | 3merge | 23340 / 103930 | 0 / 53522 |
| μ− | afterGeo | 4forward | 26799 / 97427 | 3 / 44725 |
| μ− | afterGeo | 4backward | 14577 / 91738 | 11 / 50850 |
| μ− | afterGeo | 4merge | 38769 / 120588 | 12 / 52302 |

三站selected的fallback全部为0，但全输出回退大量存在。四站selected也有3–22例，不能将三站的0推广到全部模式。merge表中全输出分母是merge最终collection；其内部forward/backward的全部调用分母另外保存在summary。未出现相同输出但来源refit status冲突。

## 执行与验证

真实full-G4为cluster1576723，gun归约v1为1577334，v2也已exit0。v2仅增加fine角分箱、未分箱quantile和1mrad尾数；本轮逐项确认v1所有已有group字段、主率、迁移与charge差值完全不变。原smoke归约先后暴露inner wrapper退出码路径和过严nominal-p精度要求，修正后smoke v3通过，均保留原失败记录。

模式计数原集群1576726/1576727，恢复1577030/1577032，归约1577336，恢复及归约全部exit0；读取有效来源表`C/mode_counter_recovery_v2/effective_sources.json`，不修改冻结原运行。所有提交均实际比较schedd负载，选兼容最低负载bigbird20，排除bigbird21，依据随submission.json保存。

分析代码：[charge_gun_analysis.py](../analysis/charge_gun_analysis.py)、[charge_mode_counter_analysis.py](../analysis/charge_mode_counter_analysis.py)。输入准备与恢复：[prepare_charge_gun.py](../scripts/prepare_charge_gun.py)、[prepare_charge_mode_recovery.py](../scripts/prepare_charge_mode_recovery.py)。

[gun全量summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/gun_analysis_v2/summary.json) · [模式精确计数](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/mode_counter_analysis_v1/summary.json)。

11页gun图最终版为`C/figures_v5/gun`，绘图1577947 exit0。92项plot-summary核对、11页vector PDF与PNG完整性检查通过；完整contact sheet及主率、角度core/tail和wall图已目视。v1/v2全部数值payload相同，v2仅改版式与fallback纵轴。新增5页机制图已经生成并目视；独立核对全部主率、refit计数、动量hist及primary-filter统计与冻结summary相同。PDF五页的graphics操作在解析ExtGState资源后均与单页相同，page boxes一致；动量图仅有等价graphics资源编号重命名，另外四页vector字节也相同。验证在各图目录的`validation.json`。

seed-hit分类发现一个覆盖问题：未出现在任何candidate measurement audit中的seed hit，旧归约将其误记为`no_SDO`。现在区分`not_in_candidate_measurement_audit`与已经审计且SDO为空的hit，并加回归测试。119事件归约v2任务1577640已exit0：34个seed分类改变；summary、pruning及除此字段外的全部event记录精确相同，31 no-SDO/25 primary-only的actual crossing结论不变。主样本同面/共同候选角度回归也已完成，见[测量与seed控制报告](2026-11-13_stage06_measurement_seed_controls.md)：updater-only的大角尾增加，both联合改善。这两项不重复运行tracking。

前次实际运行49项单测全部通过，2.137 s；`analysis/scripts/configs/tests` compileall成功，84个相关本地文档链接有效，`git diff --check`通过，日志在`C/verification_v5`。centroid smoke的三类各20事件，default58字段全部精确复现，实际cluster两库及替换计数通过；validator首轮因Gaudi属性在日志中带空格而未匹配，修正日志解析后通过，不修改重建。238个centroid/centroid+primary-filter真实重放及归约全部exit0；原59错号中centroid恢复23条selected正确、30仍错、6质量失败。只替换production seed q/p的119重放也完成：57仍错、1质量失败、1低purity，没有selected恢复。完整数值、六条pure cases及后续fixed-seed边界见上述新报告。

## 阶段性结论

事实：独立相同生成条件复现before near-pole refit导致的电荷与角度问题，1 TeV最明显；只改actual refit chart已足以大幅降低错误。主FLUKA未加权正负差有production谱贡献，独立gun也未见显著charge剩余差。100 GeV正确效率损失说明修复验收必须同时看角度、动量、效率及fallback。

假设与待判别：after剩余wrong、100 GeV效率损失、formal covariance校准及完整共同先验的非线性坐标效应尚未全部解释。这些结果支持特定操作的因果归因，不替代全链G0/G1纯旋转或所有条件下的生产修复验收。

## 下一步计划

1. 在6个pure wrong及60个正常控制中，固定原生完整seed mean并替换centroid测量，分离seed响应与CKF更新响应。
2. 已完成的100k同面/共同候选角度和资源回归与独立gun分别解释；合并已验证的新图与原参考图。
3. 完成Q5的完整共同物理先验、适用的空间/nominal控制及μ±/A′留出回归。目标保持S05/S06完整研究，不以本报告的已完成部分替代剩余验收。
