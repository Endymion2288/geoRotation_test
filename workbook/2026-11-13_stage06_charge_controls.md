# S06：production谱、方向与逐运算控制

**production-only谱平衡后，beforeGeo的μ−−μ+误判率差从原始+0.150833个百分点变为−0.018691个百分点，95%区间[−0.044406,+0.007024]；afterGeo为−0.001480个百分点，95%区间[−0.007554,+0.004594]。原始正负差异包含明显的输入谱贡献。全部59个afterGeo错误及60个正常控制的实际Calypso重放逐字段闭合，56次最终CKF翻转中31次发生在无SDO关联hit、25次在仅primary支持hit。**

## 阶段目标

在[S06全量基线](2026-11-13_stage06_charge_baseline.md)上推进roadmap Q3–Q5：区分production输入谱与入口条件，测量真实forward/backward/IFT/merge，逐hit闭合符号穿零，并检验物理seed covariance。hit/qB和100k refit控制已完成；[独立full-G4与模式fallback报告](2026-11-13_stage06_independent_charge_controls.md)补充相同生成条件的μ±及精确回退率。完整纯旋转与剩余纯轨迹机制继续按文末验收。

结果根目录 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`。全部重建、拟合、传播均为真实Calypso；两版均开启candidate-isolated smoothing。主样本、charge无关SDO matching和r95质量定义与基线相同。

## 工作内容：production-only与入口条件分别控制

production-only使用全部500,000 primary/charge；不先加truth fiducial，也不把st1位置或sensor方向当作生成条件。特征为production log p、x/y和sx/sy。discovery文件0/1/2拟合，held-out文件3/4检验，保留共同支持域、sumw/sumw2、有效N及冻结模型。

预留样本五项最大|SMD|=0.00596986，最大weighted marginal KS=0.00687522，均通过预先冻结的0.05/0.05标准。此结果与另一个包含真实上游场弯曲及接受度的入口条件overlap分开：后者held-out最大|SMD|=0.0192431、KS=0.0426882。

| 谱控制 / 样本 | before μ−−μ+ / 百分点 [95%] | after μ−−μ+ / 百分点 [95%] |
| --- | --- | --- |
| production-only，全部文件 | −0.018691 [−0.044406,+0.007024] | −0.001480 [−0.007554,+0.004594] |
| production-only，held-out | −0.026085 [−0.067004,+0.014833] | +0.001575 [−0.007485,+0.010635] |
| 入口条件，held-out | +0.021769 [−0.053414,+0.096951] | +0.002284 [−0.009223,+0.013791] |

production-only全部文件加权conditional率为before μ+0.219378%、μ−0.200687%，after μ+0.012315%、μ−0.010835%。held-out分别before0.227295%/0.201209%、after0.009971%/0.011546%。全部文件selected有效N为before248100/252858，after243638/248759（μ+/μ−）。加权truth fiducial接受度59.3611%/59.2791%；after correct/all-primary为49.9591%/49.8603%。这些outcome区间条件于冻结的discovery模型，不包含训练模型不确定度。

事实：只平衡production谱已经使before原始正负差异区间包含0，说明不能直接从未加权总率推断算法对电荷不对称。入口条件校正回答另一问题，不能将其消去的接受度效应也归为generator谱。加权结果尚不能单独判定每一个错误事件的程序或物理原因。

[production平衡与率summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/production_overlap_v1/summary.json) · [入口条件summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/overlap_v1/summary.json)。production特征与权重分别在`C/production_features_v1`、`C/production_overlap_v1`，分析入口为[charge_production_features.py](../analysis/charge_production_features.py)、[charge_overlap.py](../analysis/charge_overlap.py)。

![Production charge differences](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v3/controls_v2/production_charge_difference.png)

## 详细数据：实际方向、IFT与merge

每电荷、每版本、每模式100,000源事件，取5文件各前20,000entries；三个forward复用主样本的相同20个shards。其余实际重建400jobs全部exit0，before cluster1574281、after1574282；归约1574591成功。与主测量相同matching/r95条件，并要求下游ST1–3覆盖。4station仍是原生mask012，不冒称所有IFT六层都启用。

下表误判率与效率均用源事件/primary分母；95%为exact binomial区间。三/四站merge是实际MergeTracksAlg输出。

### 100116，beforeGeo

| 模式 | wrong / selected | misID / % [95%] | correct / fiducial / % |
| --- | ---: | --- | ---: |
| 3forward | 82 / 50335 | 0.162909 [0.129586,0.202173] | 85.577807 |
| 3backward | 120 / 51052 | 0.235054 [0.194921,0.281003] | 86.734103 |
| 3merge | 55 / 51835 | 0.106106 [0.079943,0.138089] | 88.178196 |
| 4forward | 32 / 46251 | 0.069188 [0.047329,0.097658] | 78.708150 |
| 4backward | 55 / 50027 | 0.109941 [0.082833,0.143079] | 85.099281 |
| 4merge | 28 / 50787 | 0.055132 [0.036638,0.079672] | 86.439495 |

### 100116，afterGeo

| 模式 | wrong / selected | misID / % [95%] | correct / fiducial / % |
| --- | ---: | --- | ---: |
| 3forward | 6 / 49468 | 0.012129 [0.004451,0.026398] | 84.230782 |
| 3backward | 10 / 49906 | 0.020038 [0.009609,0.036847] | 84.969858 |
| 3merge | 9 / 51760 | 0.017388 [0.007951,0.033005] | 88.128810 |
| 4forward | 0 / 43066 | 0.000000 [0.000000,0.008565] | 73.338783 |
| 4backward | 10 / 48947 | 0.020430 [0.009798,0.037569] | 83.336739 |
| 4merge | 5 / 50573 | 0.009887 [0.003210,0.023071] | 86.114233 |

### 100117，beforeGeo

| 模式 | wrong / selected | misID / % [95%] | correct / fiducial / % |
| --- | ---: | --- | ---: |
| 3forward | 176 / 52030 | 0.338266 [0.290203,0.391994] | 85.370431 |
| 3backward | 280 / 52790 | 0.530403 [0.470228,0.596113] | 86.450445 |
| 3merge | 148 / 53582 | 0.276212 [0.233553,0.324391] | 87.971683 |
| 4forward | 84 / 47841 | 0.175582 [0.140074,0.217337] | 78.625288 |
| 4backward | 116 / 51767 | 0.224081 [0.185197,0.268704] | 85.036220 |
| 4merge | 61 / 52594 | 0.115983 [0.088729,0.148960] | 86.488311 |

### 100117，afterGeo

| 模式 | wrong / selected | misID / % [95%] | correct / fiducial / % |
| --- | ---: | --- | ---: |
| 3forward | 5 / 51197 | 0.009766 [0.003171,0.022790] | 84.280540 |
| 3backward | 13 / 51784 | 0.025104 [0.013368,0.042925] | 85.233783 |
| 3merge | 12 / 53522 | 0.022421 [0.011586,0.039161] | 88.096806 |
| 4forward | 1 / 44725 | 0.002236 [0.000057,0.012457] | 73.631874 |
| 4backward | 12 / 50850 | 0.023599 [0.012194,0.041219] | 83.697728 |
| 4merge | 10 / 52302 | 0.019120 [0.009169,0.035159] | 86.091538 |

4forward的低misID伴随明显效率损失：after正确效率约73.3–73.6%，3forward约84.2–84.3%；零错误保留上限，不能单独以错误数最少选最佳模式。3merge提高效率到约88.1%，同时保留各模式的配对7×7迁移和truth-p分箱。

每条merge输出的完整hits、q/p及upstream参数都精确对应独立forward/backward输出，全部native输出及selected/wrong来源分别统计；没有近似hash匹配。selected主要来自backward，例如after μ+三站merge为44367 backward、7393 forward。初版100k生产未开启refit计数，MessageSvc抑制重复warning；后续400份只读计数重放的58字段均精确复现。归约1577336成功：所有三站selected fallback为0；四站selected为3–22条，详细全输出/selected分母及恢复记录见[模式计数报告](2026-11-13_stage06_independent_charge_controls.md)。

[模式summary与迁移](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/mode_analysis_v1/summary.json) · [charge_modes.py](../analysis/charge_modes.py)。

![Native mode misidentification](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v3/controls_v2/mode_misid.png)

![Native mode efficiency](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v3/controls_v2/mode_efficiency.png)

## 详细数据：119事件逐hit重放与pruning

全部59个afterGeo最终wrong及60个正常控制已重放（每charge每文件最早6个满足sample_entry mod997=0的正常事件）。119/119作业exit0；每event全部58重建数组与原5000-event源shard精确相同。归约1574592成功，全部measurement hit的filtered q/p及首次/最后stage标签与原lineage闭合。

直接从actual predicted mean/C、measurement H/V计算：`r=m-Hx`，`S=HCHᵀ+V`，`K=CHᵀS⁻¹`，并验证`Δκ=(Kr)κ`。两电荷wrong更新的最大绝对恒等式误差为1.164×10⁻¹⁵/1.658×10⁻¹⁵ GeV⁻¹，normal最大1.961×10⁻¹⁴ GeV⁻¹。

| 最后一次进入最终错误的CKF update | μ+ | μ− | 合计 |
| --- | ---: | ---: | ---: |
| hit无SDO关联 | 15 | 16 | 31 |
| hit仅primary支持 | 11 | 14 | 25 |
| ST3 | 23 | 26 | 49 |
| ST2 | 2 | 4 | 6 |
| ST1 | 1 | 0 | 1 |

各翻转innovation χ²均在0–25内，符合当前measurement selector；每次r/S/K、完整Cκ相关、前后κ与significance保存在`events.json`。无SDO关联不等于已经证明电子噪声；primary支持也不证明seed/路径正确。59条中53条final purity<1，6条purity=1（μ+2、μ−4），不能将全部剩余错误归为其他粒子混入。

有10个wrong事件存在primary-supported正确符号candidate被pruned，两电荷各5。全部这些备选都没有majority tie。淘汰reason bits：1为nMeasurements≤12，2为shared hits>6，3为两者。完整候选和winner见`pruning.json`。存在正确备选不证明可以无损恢复，因为改变候选规则也可能改变正常轨迹和duplicates。

[逐运算summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/replay_analysis_v2/summary.json) · [119事件详细表](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/replay_analysis_v2/events.json) · [pruning表](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/replay_analysis_v2/pruning.json)。v2单列未进入measurement audit的seed hits，34个events的该字段改变，其余event内容、summary及pruning与v1精确相同。源matrix audit为各重放目录的`charge_audit.jsonl.gz`，来源与58字段验证随输出保存。[charge_replay_analysis.py](../analysis/charge_replay_analysis.py)为实际归约入口。

![Actual charge-crossing innovations](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v3/controls_v2/charge_crossing_innovations.png)

![Six purity-one wrong paths](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v3/controls_v2/pure_wrong_sign_paths.png)

## 详细数据：一阶物理seed covariance匹配

真实Calypso worktree `geoRotation_calypso_charge_prior`以88c457c1+0021+0024构建。默认关闭的MatchFaserSeedCovariance保留native mean，用实际F/A seed surface basis和角度Jacobian变换完整6×6 covariance：`C_A=J C_F Jᵀ`。before配置的local swap在该转换中显式处理；mean、field、材料、finding与refit不同时修改。

100044/100116/100117的20-event默认行为与baseline58字段精确相同，audit ON/OFF精确相同，seed mean逐项不变；实际seed行数20/17/18。中心差分Jacobian最大scaled误差9.37×10⁻¹⁰，C映射误差0，同一阶physical covariance相对误差≤1.054×10⁻¹⁵。记录先验的纯代数抽样显示，有限宽度球角Gaussian经过非线性映射，其物理covariance可差约2.8–15.7%；因此**一阶匹配不等于完整非线性先验分布相同**，也不作为全链纯旋转闭合。

全5文件500k/charge控制重建1574895，200/200成功；归约1575201成功。全部919232个seed mean与native逐项相同，q/p variance均1×10⁻⁵ GeV⁻²。主选择结果：

| charge / prior | wrong / selected | misID / % [95%] | correct / fiducial / % |
| --- | ---: | --- | ---: |
| 100116 after matched | 19 / 254297 | 0.007472 [0.004498,0.011668] | 86.606653 |
| 100117 after matched | 26 / 262987 | 0.009886 [0.006458,0.014486] | 86.546164 |

相对native after，μ+ conditional变化−0.003452个百分点，95%[−0.007139,+0.000235]；μ−−0.002633个百分点，95%[−0.006922,+0.001656]。两者区间均含0。correct效率从84.1806%/84.1169%提高至86.6067%/86.5462%，说明原生prior选择显著影响效率，不能把原生版本全部变化归于坐标旋转。

配对迁移中，native-after wrong→matched correct为7/13，wrong→wrong为12/13，wrong→quality_failed为8/6；另有native correct→matched wrong为3/1。这些迁移分别保存，不把失去轨迹算成符号修复。全native-final的refit fallback为69544/65980，分母453936/453896；这是所有输出轨迹的计数，不能与主selected分母混用。补充归约`prior_analysis_v2`确认主selected的254297/262987条全部refit成功、fallback均0；新增计数与v1的全部率、histograms、迁移逐项相同。

[prior验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/prior_validation_v1/summary.json) · [全量prior结果与迁移](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/prior_analysis_v1/summary.json) · [patch0024](../patches/0024_matched_seed_prior_after.patch)。

## 详细数据：primary hit 与 q/B 单因素控制

238个实际重放全部成功，修正后的归约1576722成功。source-link过滤保持所有119个事件的seed mean与完整covariance不变，继续使用native组合CKF；它只移除没有指定primary SDO支持的source links。所有保留的selected measurement均有primary支持。该控制允许branches和winners改变，不能称为固定完整hit路径。

| 原cohort | 源事件 | 恢复/保持correct | 仍wrong | quality failed |
| --- | ---: | ---: | ---: | ---: |
| μ+，原wrong | 27 | 8 | 2 | 17 |
| μ−，原wrong | 32 | 9 | 4 | 19 |
| μ+，原normal | 30 | 30 | 0 | 0 |
| μ−，原normal | 30 | 29 | 0 | 1 |

59个wrong中17个恢复correct，6个仍wrong，36个quality failed；本次没有unreconstructed/unmatched、low-purity或invalid-curvature outcome。后36个不能计作符号恢复。这是刻意选择的全部错误与冻结正常控制，不能将17/59外推成对全500k样本的无损改善率。

剩余6个正好是原final purity=1的事件，最终q/p逐值不变。它们的sample_entry为μ+51839、224394，μ−227965、390966、439421、458033。`q_truth κ/σ_κ`分别为−0.8122、−6.0746、−1.0876、−1.6768、−1.4242、−1.3515。一个事件有较大的formal wrong-sign significance，不能把全部六例直接解释为统计上不可分辨的曲率；当前covariance是否校准与逐事件物理原因仍需检验。

q与B同时反号的119事件中，全部6902项数组比较通过；126个seed及6122个路径states的P/F/S mean与完整6×6 covariance按q/p行列反号闭合。整数身份/命中计数精确相同；浮点表使用rtol/atol=10⁻¹⁰，state covariance使用10⁻⁹，见实际[validator](../scripts/validate_charge_hitcontrol.py)。该控制支持seed/field/refit/conversion的qB符号协变性，不代表真实μ+、μ−随机G4过程或apparatus具有精确对称性。

[238事件summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/hitcontrol_analysis_v2/summary.json) · [逐事件表](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/hitcontrol_analysis_v2/events.json)。分析入口[charge_hitcontrol_analysis.py](../analysis/charge_hitcontrol_analysis.py)。

## 详细数据：beforeGeo refit chart 的实际电荷控制

每charge100k、native/smoother/updater/both四组；native复用原同输入shards。120个新增重建全部成功，归约1576420成功。两电荷的6173405/6190328条非final CKF记录在三种控制下均与native逐项相同；变化仅发生于refit及其后输出。

| charge / 控制 | wrong / selected | misID / % [95%] | correct / fiducial / % | selected refit correct→wrong |
| --- | ---: | --- | ---: | ---: |
| μ+ native | 82 / 50335 | 0.162909 [0.129586,0.202173] | 85.577807 | 49 |
| μ+ smoother | 81 / 50342 | 0.160899 [0.127797,0.199944] | 85.591431 | 48 |
| μ+ updater | 5 / 50289 | 0.009943 [0.003228,0.023201] | 85.630598 | 1 |
| μ+ both | 5 / 50318 | 0.009937 [0.003227,0.023188] | 85.679984 | 1 |
| μ− native | 176 / 52030 | 0.338266 [0.290203,0.391994] | 85.370431 | 123 |
| μ− smoother | 175 / 52031 | 0.336338 [0.288417,0.389924] | 85.373724 | 122 |
| μ− updater | 14 / 51928 | 0.026960 [0.014740,0.045231] | 85.469213 | 4 |
| μ− both | 14 / 51999 | 0.026924 [0.014720,0.045169] | 85.586105 | 4 |

both−native conditional误判率的配对差：μ+−0.152972个百分点，95%[−0.187115,−0.118829]；μ−−0.311343个百分点，95%[−0.360386,−0.262299]。变化显著，且correct/fiducial效率略升，没有通过大幅丢弃正常轨迹制造此改善。只改smoother的误判率变化区间包含0。

预留文件3/4同样成立：μ+ native35/20084→both2/20076，μ−68/20756→4/20747。这里的每文件范围仍为原20k，held-out合计40k/charge；不称为完整500k修复验证。

所有主selected refit均成功、fallback=0。全native-final输出的fallback则为μ+14224/90890、μ−13492/90771，四控制相同。两个分母分开保存；未用质量选择后的0替代全输出回退率。完整7×7迁移、动量响应、held-out、成功/回退与新生/恢复wrong见[refit summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/refit_analysis_v1/summary.json)。

事实：beforeGeo额外错号的主要来源被单因素定位到near-pole方向参数化下的refit measurement update。该结论由固定CKF、实际native updater与regular-chart updater比较支持；smoother-only不足以移除主要错号。它与S05角度结论使用同一真实refit组件，但不把方向分布变好单独当作电荷归因，也不将该100k结果扩称为所有版本差的纯旋转效应。

## 执行、图表与验证

所有批量任务先检查schedd负载，选择兼容且running+idle最少的bigbird20，排除bigbird21；各submission.json保留UTC查询、负载、请求资源及cluster。主重建和归约冻结源码与软件后运行。production-only1574283、mode归约1574591、replay归约1574592、prior验证1574590、prior全量1574895与归约1575201均exit0。

新增9页控制图cluster1575504成功，输出[C/figures_v3/controls_v2](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v3/controls_v2/figure_index.md)；修正seed小柱标注的7页诊断在[C/figures_v3/diagnostics](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v3/diagnostics/figure_index.md)。新增3页prior图cluster1576113成功。全部36张图的三张contact sheets及主要新图已目视检查；PNG完整性、合并PDF逐页vector内容及page boxes一致性通过。完整[36页更新图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v3/charge_stage06_report.pdf)及[逐页索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v3/figure_index.md)包含12主图、7诊断、5历史、9控制及3prior。原主图及历史Alma图保持原定义；新控制图全部英文，并保存PNG、逐图PDF、counts/sumw2/flows和index。

本轮修正的q/B observer只影响新增控制：首版TruthChargeScale同时乘了已反号的reco charge，导致projection不同；原生seed/P/F/S及完整C均闭合。默认58字段与基线精确一致。失败输出和日志保留在`C/hitcontrol_smoke_v1`，v2只在原始production-truth传播处翻号，不对reco charge重复操作。修正的observer build1575808成功；新9组smoke1575809的默认58字段及全部seed/P/F/S/完整C和投影均闭合，验证在`C/hitcontrol_smoke_v2/validation.json`。238个全119事件hit过滤/qB重放1576115全部exit0，机制归约1576722成功，实际数据见上述hit/qB表。

before refit机制控制用master+0022实际库，并复制S05已验证的5个refit chart文件；不替换CKF finding。增量[patch0027](../patches/0027_charge_refit_chart_before.patch)，source/provenance在`C/refit_control_provenance.json`；build1575807成功，8个native/smoother/updater/both小样本1575810全部exit0，默认58字段精确复现，CKF所有非final行不变。100k/charge三控制共120jobs1576116全部exit0，归约1576420成功，实际结果见上表。

hit/qB首轮归约1576276因raw worker表缺少`sample_entry`而exit1，失败输出保存在`C/hitcontrol_analysis_v1_job`。已用冻结task的skip、sample_entry、file_index补齐分析索引，不更改重建数组或238份实际Calypso重放。修正版冻结于`C/hitcontrol_analysis_v2_code`，归约1576722成功。

独立full-G4 μ± ON的2×20-event smoke1576117全部步骤exit0；[实际验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/symmetric_gun_smoke_v1/validation.json)确认两电荷实际production position/momentum精确相同，三种重建truth相同，四份实际加载库与各自配置一致。100/1000 GeV各20k/charge、32jobs的full-G4控制1576723全部exit0，归约v1/v2也均成功，原有率、迁移、直方图在细化角分箱后完全不变；[独立控制报告](2026-11-13_stage06_independent_charge_controls.md)记录全部12组实际数据。使用真实apparatus，不假定空间反射对称。

实际执行`python3 -m unittest discover -s tests -v`，48项全部通过（2.152s），日志`C/verification_v4/unit_tests.log`。`compileall`覆盖analysis/scripts/configs/tests成功，91个本地文档链接全部有效，记录在`C/verification_v4/compile_and_links.json`；本轮`git diff --check`通过。图集实际合并入口为`analysis/charge_report.py --figures ... --output C/report_v3 --pdf-tools <stage05_geometry_comparison/software_tools>`，逐页内容与box校验保存在`C/report_v3/validation.json`。

## 阶段性结论

事实：production谱平衡后，before/after两版的正负差异都与0相容；真实模式变化带来可量化效率/误判权衡。after剩余符号翻转的actual update恒等式已闭合，31无SDO与25仅primary支持说明需要hit和seed的独立控制。matched prior显示原生covariance是版本比较中的实际混杂因素。

事实：100k/charge的实际refit控制把before额外错号主要定位到near-pole chart的measurement update；primary source-link控制恢复17个wrong，但另有36个quality failed，6个纯轨迹错误完全不变。q/B全119事件路径与完整C闭合。

后续[测量、seed与角度回归](2026-11-13_stage06_measurement_seed_controls.md)已测量centroid及production seed q/p-only：6个pure wrong在centroid后5条selected正确、1条符号正确但质量失败，q/p-only全部仍错。centroid也改变seed完整均值，下一步固定这些均值以继续归因。无SDO hits可能含噪声或缺关联，不能从audit缺失反推污染。全500k的834个refit新生错误不按100k比例直接外推。

## 下一步计划

1. 在已完成的centroid与q/p-only控制基础上，固定原生完整seed均值重放6个pure wrong及60正常控制，进一步分离测量与seed响应；保留完整质量迁移。
2. 补before refit控制的角度/CPU联合回归及新机制图；在新样本和共同候选中报告效果，保留全输出与selected的两套fallback分母。
3. 完成独立full-Geant4 μ± gun图组验证和主报告汇编；保留100 GeV效率损失及1 TeV残余，量化空间分组与区间。
4. 完成纯坐标控制的非线性prior及剩余parameterization路径、留出μ±与A′回归。模式fallback精确计数已完成。Q5尚未完成；本阶段不宣称纯旋转或完整生产修复通过。
