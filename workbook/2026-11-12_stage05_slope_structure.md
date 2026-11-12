# S05：角度结构定位与实测进展

本阶段尚未通过 A1–A5 完整验收。已完成 100044 的 500,000-event beforeGeo/afterGeo 重建、物理角度与原 writer 闭合、预留样本响应检验及真实 Calypso/G4 控制。**历史表达式存在 frame/选样问题；strip coordinate 造成 seed 之前的 segment 周期结构，并对最终 st1 两个肩峰有独立40k事件的配对证据；VetoNu 残差的主要宽度与实际靶材料散射相符。旋转前近零增强已定位到 refit smoothing，具体数值原因正在用原生 Calypso 最小控制验证。**

结果根目录：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes`，下文记为 `S`。本报告保存已经测得的结果；后续补充不修改原 discovery 窗口和运行快照。

## 阶段目标

按 [A1–A5](../roadmap/03_slope_structure.md) 回答剩余角度结构首先出现在哪个环节、是否由 truth 谱/接受度/测量离散性和响应混合解释，以及最小修改是否在独立输入上消除对应异常。按[参考图规范](../roadmap/05_reference_figures.md)输出 `problem.pdf` p4、`geometryRotation.pdf` p19–20、p23 对应图和逐阶段诊断。

## 输入、软件与定义

主输入为第 19 页指定的 100044 μ+ 单个 500k production unit，完整文件：

```text
/eos/experiment/faser/data0/sim/mc24/particle_gun/100044/rdo/s0013/FaserMC-MC24_PG_mupl_fasernu_5mrad_flukaE-100044-00000-00099-s0013-RDO.root
```

输入按 entry `[0,250000)` discovery、`[250000,500000)` validation 分开。每版本 100 个 5000-event shards 全部成功，覆盖统计保存在 `S/main_summary_v1/coverage.json`。生产 run/event 在不同 unit 复用，身份使用输入文件和 source entry。

| 项目 | beforeGeo | afterGeo |
| --- | --- | --- |
| Calypso | `geoRotation_calypso_main`，master `1eb885ae22299fb39bb352bead6005d41b0f953e` | `geoRotation_calypso`，base `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` |
| 实际软件 | `S/software_before` | `S/software_after` |
| 材料 | `stage02_closure/material-current-source.json` | `/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json` |
| 重建 | 原生 CKF2、conversion 和 KalmanFitterTool；3station forward、mask 0–5 | 同样 collection、方向和 mask，保留该分支原生 covariance |

两版使用 Athena 24.0.41、ACTS 32.0.2、FASERNU-04、OFLCOND-FASER-05。历史 RDO metadata 则是 Athena 22.0.49、OFLCOND-FASER-04。本次为实际版本整体比较；两分支 seed covariance 不同，不能把净变化称为匹配物理先验的纯数学旋转。所有重建和传播都调用真实 Calypso，模拟调用其完整 G4/digitization；没有独立 tracking 或简化 detector。

物理 frame 的 `s_y=py_F/pz_F=py_A/px_A`，旋转为 `p_A=(pz_F,py_F,-px_F)`。fixed st1 为 `z_F=0 mm`；VetoNu reference 为 `z_F=-3112 mm`。实际 geometry 的两块 VetoNu 板中心是 −3122.5、−3101.5 mm，厚 20 mm，reference 在间隙内，不能称为板中心。

SiHit truth 是同 sensor 中心处的 G4 chord 方向。fixed-plane model 将 **production momentum 大小附在首个 st1 SiHit midpoint/chord 上，再经 Calypso 外推**；不是实测 sensor 局部动量。旧 RDO 没有 VetoNu truth hits，因此用新的完整 G4 gun 检查该面的真实响应。`segments.parameters_F` 的 momentum 单位为 MeV，大小固定 100000 MeV，只用于方向，未把这个约定大小当成真实粒子动量。

## A1：frame、writer 和参考宏

全量 long tracks 的 raw→F momentum 和 F/A 两种 slope 算法最大差均为 0。全部 upstream stored surfaces 均不同于 fixed st1 z=0，主图使用实际 Calypso 外推后的同一物理面。

| 计数 | beforeGeo | afterGeo |
| --- | ---: | ---: |
| 源事件 | 500000 | 500000 |
| 原生 final tracks | 484845 | 484931 |
| 三站 long tracks | 436859 | 432124 |
| 物理 upstream p=50–250 GeV | 176702 | 157132 |
| 历史 event-first raw pz cut | 176714 | 76 |
| 原始 st1 legacy angle 超出 ±0.01 | 16156 | 429790 |
| fixed st1 外推失败 | 0 | 0 |
| VetoNu 外推失败 | 5 | 12 |
| sensor 缺少 truth reference | 15 | 20 |

afterGeo 的 legacy `atan(raw_py/raw_pz)` 分母不再是束流方向；沿用 raw pz cut 会把 157132 条物理动量窗口 tracks 换成仅 76 条历史 cut tracks。该变化属于输出解释和选样，不能拿 legacy 曲线直接判断旋转后的物理角度。

`S/writer_checks_v2/results/summary.json` 比较 20-event smoke 和另外 500 events 的实际 NtupleDumper flag false/true：events/truth/tracks/segments 的全部保存字段逐值相同，保存的非角 Ntuple branches 逐值相同，正确角度与物理公式闭合至约 1e−17。未保存全部 Ntuple branches，验证范围为 summary 列明的 branches。

原 GitHub `test_tanTheta.C` 冻结 commit `33ca2281b961935456dfd4dc9fc06235e24b551e`，路径为 `root_txt/geoRotation/after_vs_before/test_geoMerge/afterGeo_afterModifyErr_test/test_tanTheta.C`。实际启用表达式 before 是 theta_x 加逐 track `100<pz<250 GeV`，after 是位置比 `Y/sqrt(X²+Y²)`、没有 p cut；源代码只 draw h1。用真实 writer 保存的 branches 组成小 trees，逐字未改的原宏产生 h1=7、h2=14 entries，全部 flow bins 与 Python literal 计算一致；无用的 truth/p0 branches 用空 vector，此项验证算术、cuts、ROOT Fill，没有恢复原作者全部 ROOT 输入。

实际 rotated writer 的 X 为 `z_A=-x_F`，Y 为 `y_A=y_F`，所以该位置比分母确为 `sqrt(x_F²+y_F²)`。曾在 `main_analysis_v4` 错改 literal 坐标，原宏核对发现后已纠正。`S/main_summary_v1/provenance.json` 记录最终汇总采用 v4 物理统计、恢复 v3 正确 after literal histogram 的来源；不声称所有版本 JSON 逐字相同。

## 全量效率、动量和电荷伴随指标

匹配不依赖电荷或 reconstructed p。效率母分母是 414417 个 fiducial H18 truth μ+；numerator 要求 long、Good95、purity≥0.7，并按 truth 去重。misID 和 p-tail 则以 matched primary long tracks 为分母，未额外施加上述效率 cuts。

| 指标 | beforeGeo | afterGeo |
| --- | ---: | ---: |
| truth 效率 | 351504/414417 = 84.8189% [84.7078,84.9343]% | 347615/414417 = 83.8805% [83.7668,83.9976]% |
| charge misID | 18765/436634 = 4.29765% [4.23509,4.35897]% | 3524/431840 = 0.816043% [0.786854,0.842346]% |
| p_reco/p_production<0.5 | 41373/436634 = 9.47544% [9.38709,9.56301]% | 7887/431840 = 1.82637% [1.78604,1.86713]% |

区间为 1000 次 source-event bootstrap 的 95% 区间。效率降低约 0.94 个百分点，需与曲线改善同时报告。此表比较各版自己的最终样本，未作为固定同候选的旋转因果估计。

## A2–A3：结构、谱与响应

主统计已保存共同分箱及半 bin offset、ECDF、KDE bandwidth scan、未舍入值的 unique spacing、truth/reco/residual/pull、p_truth/位置/module/hit 数/cluster width/purity 分组。两版本在 st1、VetoNu、sensor 的有限正确 slopes 全部 unique，exact zero=0。**没有 exact-zero 堆积不等于没有近零异常结构。**

旧峰窗口仅在 afterGeo discovery 找出，然后冻结用于两个版本。afterGeo st1 p_reco=50–250 GeV 的窗口中心分别是 −1.89、−0.71、+0.59 mrad，半宽分别 0.5、0.107219、0.1 mrad。不能用它们替代 beforeGeo 明显近零结构的独立窗口；已另完成 beam-centred 25/50/100/200 µrad 固定半宽的 validation 检查，见 `analysis/slope_zero.py`、`S/zero_analysis_v1/`。

近零统计与 response 都要求 primary 匹配且 purity≥0.7 的 long tracks，另加标明的 p 窗口；分母因此小于未加 truth 匹配的主图。±100 µrad 的 st1 窗口中，beforeGeo reco-p 选样为 5574/88182=6.3210%，相同选中粒子的 truth 为 3222/88182=3.6538%；truth-p 选样仍为 reco 2413/78940=3.0568%、truth 1727/78940=2.1877%。因此 p 选样迁移不能单独解释近零尖峰。phase response 对 beforeGeo 预测 4287.96 条，实际多出 1286.04 [1106.37,1441.49] 条；这个差值是该条件模型的验证偏差，仍需处理未覆盖部分和模型系统项。afterGeo 相同 reco-p 窗口为 1852/78360=2.3635%、truth 为 1854/78360=2.3660%；phase response 预测 1809.80 条，差 42.20 [−36.49,109.14]，与零相容。sensor 的 truth-p 对照同样保留 beforeGeo 的近零过量，排除仅 fixed-plane 外推制造结构的解释。

经验 response 使用 discovery 的 residual 和 track multiplicity，按 production p、同面 truth slope、sensor x/y 分组，应用于 held-out truth。原 coarse response 的 afterGeo st1 reco-p shape TV 为 2.4413%，VetoNu 为 1.3566%；低响应分组的覆盖不足单列，不能把全部结构称为已解释。

实际 geometry 为 strip pitch=0.08 mm、三层 lever arm=63 mm、stereo≈±0.01999867 rad，几何周期 `pitch/lever_arm=1.269841 mrad`。`S/phase_response_v1` 使用此固定周期的 8 phase bins，没有自由扫描周期或拟合 Gaussian 分量。afterGeo held-out st1 p_reco 窗口如下：

| 窗口中心 | observed | phase prediction | 差值和 source-event bootstrap 95% 区间 |
| --- | ---: | ---: | --- |
| −1.89 mrad | 11534 | 11333.32 | 200.68 [64.57,302.98] |
| −0.71 mrad | 2552 | 2573.00 | −21.00 [−88.55,51.61] |
| +0.59 mrad | 1728 | 1655.44 | 72.56 [8.70,134.02] |

st1 shape TV=2.3462%；VetoNu=1.6119%。该 st1 响应未覆盖 96926/244741 个 validation truth，对应 2486/78360 条实际选中 reco；VetoNu 未覆盖 92522/244636 truth、1549/78360 reco。未覆盖 reco counts 已逐 bin 保存。+0.59 mrad 的差值不能直接换算成全样本额外峰的严格上限，也未做多窗口 familywise 显著性声明。−0.71 mrad 结构与几何 phase 条件响应相容，支持测量/响应相关来源；经验模型仍从同一模拟重建中学习，不能据此断言全部为物理效应。

上述历史差值将全部 validation reco 与仅有 covered truth cells 的预测比较。`S/phase_response_v2` 保留全部原字段，另补相同覆盖范围的比较；两版两个面的旧字段均与 v1 逐项相同，covered+uncovered 的 bins、总数及 unbinned frozen-window counts 全部闭合。afterGeo 的结果如下，区间仍为 source-event bootstrap 95%，未包含响应分组选择等模型系统项：

| st1 窗口中心 | covered observed | prediction | covered 差值及区间 | 该窗口 uncovered reco |
| --- | ---: | ---: | --- | ---: |
| −1.89 mrad | 11455 | 11333.32 | 121.68 [−9.05,223.28] | 79 |
| −0.71 mrad | 2508 | 2573.00 | −65.00 [−135.57,8.04] | 44 |
| +0.59 mrad | 1688 | 1655.44 | 32.56 [−29.25,98.16] | 40 |

三窗口的 covered 差值、均值差和宽度差区间均包含零，st1 covered shape TV=1.4293%，VetoNu=1.5282%。这支持已覆盖条件中的 shoulders 与测得响应相容；未覆盖部分仍单列，不能把其数量当成已解释或删去后宣称全样本没有额外结构。完整新图和逐窗口数据见 `S/phase_response_v2/`，核对记录为其 `validation.json`。

## A4：真实 Calypso/G4 控制

### 自由方向与逐阶段审计

每版本 294 个等价 free states，扫描正负近零 slopes、50/100/1000 GeV、q±，传播到 z_F=±100 mm：各 588 次成功、0 失败。before/after physical mean 差最大约 2e−13 GeV、1e−14 mm，physical Var(sy) 比在 1±3e−15 内。before 的 6 个 exact-pole states 无定义的 phi covariance 单列，只传播 mean。结果证明此次自由传播数值闭合，不证明整个近极点 CKF 更新正常或两版重建 covariance 相同。

afterGeo 已审计 4000 events（`[0,2000)` 与 `[250000,252000)`），3452 条 matched pure long tracks。从真实 seed、CKF P/F/S、converted、refit、final 经实际 Calypso 投到 st1/VetoNu/sensor，final projection 与 observer 最大误差为 0。结果位于 `S/stage_analysis_v2/results/`。部分窄窗峰与匹配 sideband 统计少，不能据少数事件终止该分支。

另做 2000-event 仅切换 candidate-isolated smoothing 的配对：P/F、seed、selected target、converted、refit、final 的 means 均逐值不变。共享 slots 的角度 smoother 恒等式超阈值项在两 shards 为 1/7，隔离后均为 0；final mean、charge、hits、projection 全部不变，Good95 不变。故它可污染逐 state 诊断，但不解释这批最终角度结构。详见 `S/checks_analysis_v1/smoothing/summary.json`。

为定位 beforeGeo 的近零尖峰，开启其本来就有的 `RootTrajectoryStatesWriterTool`（CKF/refit）和 `RootSeedWriterTool`，无需另一套 tracking 或重编 master。20-event/版本 smoke 全部成功；各自的 5/15/30/6 个保存 event/truth/track/segment branches 与原基线逐值相同。refit 分别 19/15 条，均以同一事件、完整 module/local measurement signature 唯一匹配，无缺失或歧义。源事件身份还以 observer counts 和 CKF 顺序交叉核对。`S/native_states_v1` 每版 10000 events 均已 exit0：beforeGeo refit 9110 条中 149 条未匹配，afterGeo 8928 条中 65 条未匹配，均无歧义；未匹配者不进入 same-hit stage 对照。seed 暂按 population 保存，未伪称已关联同一 ancestor seed。`S/native_transport_v2` 四个 jobs 已全部 exit0，实际 Calypso 传播后由 `S/native_stage_analysis_v3` 归约。

首轮 transport_v1 的四个 jobs 在配置阶段退出：主重建的旧 observer 库没有注册后来新增的 `CalypsoSlopeTransportProbeAlg`。未处理事件，也未损坏已有 states。v2 保持 inputs 不变，改用此前 direction/stage 控制已成功的 `S/diagnostics_v1/software_{beforeGeo,afterGeo}` transport 组合，不重建 tracking；失败日志及恢复来源均保留。

逐阶段归约 v2 在身份检查退出：native 输入按 stage 分块，transport 按 event 输出。修订 v3 按唯一 `(source_entry, track ordinal, stage)` 对齐，检查完整身份集合及 barcode/layer 后读取状态，不重新传播。两个新的回归检查覆盖顺序交换及等行数但身份不同的情况。10k/版本的全部 event/truth/track/segment 保存字段均与原全量基线逐值相同，final 同面投影最大绝对差为 0。

固定同一 final candidate、按其 native final p=50–250 GeV 选样，再限制 validation entries≥250000、primary μ+、purity≥0.7 和 long，完整 CKF/refit 的同面比较如下。数字为 `|sy|<100 µrad` 的条数/共同分母，未按各中间 stage 的 p 重选：

| stage | beforeGeo /1735 | afterGeo /1569 |
| --- | ---: | ---: |
| CKF predicted | 45 | 34 |
| CKF filtered | 45 | 34 |
| CKF smoothed | 50 | 31 |
| refit predicted | 51 | 32 |
| refit filtered | 56 | 30 |
| refit smoothed | 108 | 36 |
| final | 108 | 36 |

beforeGeo 的 refit filtered→smoothed 由 3.2277% 增至 6.2248%，afterGeo 为 1.9120%→2.2945%。beforeGeo 完整 cohort 覆盖 1735/1758 条该选样 final、108/109 条近零 final；afterGeo 为 1569/1590、36/36。改用 truth-p 窗口时 beforeGeo 同一分母1579中的近零条数为35→55，afterGeo 同一分母1559为26→32；sensor 面也有对应增强。这定位到 refit smoothing，不能归因于 final writer 或 st1 外推。尚不能仅凭 stage 计数将全部近零增强归为某个实现缺陷。完整分母、1000次 event-bootstrap 区间及 seed population 在 `S/native_stage_analysis_v3/{beforeGeo,afterGeo}/summary.json`；每版有10张PNG/逐图PDF和10页汇编。

后续单因素为原生 `KalmanFitterTool1.SeedCovarianceScale` 的1/10/100对照；CKF配置、输入、geometry和material固定。两版 CKF2Config 实际将此值设为10，C++类默认100。最初 smoke 错将100当成名义配置，逐值验证发现 final 差异后已在 `S/refit_covariance_controls_v1/baseline_correction.json` 更正：100是敏感性变化，原配置为10。补做的 scale10 20-event/版本 smoke 的 events/truth/tracks/segments/native_states 全部保存字段与原基线逐值一致，见 `native_smoke_validation.json`。生产对照每个配置为10000 events，12个jobs已全部exit0；分析入口为 [slope_refit_covariance.py](../analysis/slope_refit_covariance.py)。

固定原scale10候选、原majority、原native-p=50–250 GeV窗口、primary/purity/long和validation范围，不按改变后的动量重选。这个cohort包含前述未获得完整refit state链的少量候选，因此分母为1758/1590，与逐stage表1735/1569不同：

| refit scale | beforeGeo近零数/1758 | 相对原scale10变化及95%区间 / 百分点 | afterGeo近零数/1590 |
| --- | ---: | ---: | ---: |
| 1 | 67 | −2.3891 [−3.5502,−1.1434] | 35 |
| 10（原配置） | 109 | 0 | 36 |
| 100 | 103 | −0.3413 [−1.5093,+0.8766] | 35 |

scale1的beforeGeo完整same-hit子组1739条中近零64条，配对变化−2.5302 [−3.7321,−1.3019]个百分点；truth-p窗口同分母1595条56→40，变化−1.0031 [−1.8150,−0.1952]个百分点。故scale1差别在固定原动量选择及相同hits时仍存在。三组全部seed/CKF保存向量、事件、truth和segments逐值相同；scale10全部final字段与原10k逐值闭合。afterGeo两个非原scale的变化均为−0.0629个百分点，区间[−0.1946,0]。

原scale10的beforeGeo cohort中，186条有refit smoothed φ超出[−π,π]、154条有θ超出[0,π]、54条有|φ_smt−φ_prt|>π；109条近零track中对应为46/46/15条。afterGeo三个scale的这些flags均为0。越界/相关性用于定位数值问题，不把所有近零tracks归给单一flag，也不将缩小covariance的敏感性结果当作正式修复。

在同一原native-p cohort内，beforeGeo的misID counts（scale1/10/100）为39/37/87，动量负尾为88/191/268；afterGeo为3/5/11及35/32/32。缩小scale1时近零峰及动量负尾减少，但没有证据说明电荷误判同时改善。

结果在 `S/refit_covariance_analysis_v2/{beforeGeo,afterGeo}`，各有1张PNG/逐图PDF及报告。首轮归约1573944用目录名把诊断s0001关联到不存在的原生s0001，实际同一source range在原任务名为s0050，因而在读数据阶段失败；v2按输入路径、skip和events关联源shard，保留全部state/candidate断言，只重做统计（1573947），没有重复重建。

新增 raw angle 诊断保留真实 native writer 的逐measurement φ、θ、q/p及对角误差，未在分析中重新拟合。scale100 smoke 中 beforeGeo refit 出现 φ超出[−π,π]、θ<0，afterGeo未出现相同现象；这是扩大 covariance 的受控配置中的线索，不当作原配置的定量结论。ACTS32.0.2 对应发布源码 [GainMatrixSmoother.cpp](https://github.com/acts-project/acts/blob/v32.0.2/Core/src/TrackFitting/GainMatrixSmoother.cpp) 使用普通参数差分；当前 upstream 已使用 periodic bound-parameter subtraction 和 normalization。这提供检查角度 chart crossing 的具体路径，仍需实际 Calypso 单因素和独立输入验证。

为验证该具体假设，在同一 master commit 的实际 Calypso worktree `geoRotation_calypso_slope_refit` 实现 [0011](../patches/0011_refit_periodic_phi.patch)。新属性 `KalmanFitterTool1.WrapRefitSmoothingPhi` 默认 false；true 保留 ACTS32 原 `GainMatrixSmoother::calculate` 的 gain/covariance 运算，只将 smoother 中的 φ 差分改为2π周期差分；不改 θ normalization、CKF、refit forward 或 covariance scale。编译 cluster1573937 成功，53个实际 translation units，构建和完整 patch 冻结在 `S/periodic_refit_v1/build`。软件由原 `S/software_before` 叠加新编译的实际组件和 Gaudi metadata 组成，逐文件来源见 `software/software.json`。

20-event off/on smoke（1573939）均exit0：关闭属性时，events/truth/tracks/segments/native_states 全部保存字段与原20-event基线逐值相同；开启时 events/truth/segments、80条 seed/CKF states 和38条 refit P/F states仍逐值相同。19条 refit S 中1条变化；这是受控开关工作的验证，不能由20事件判断峰因果。对应验证为 `native_smoke_validation.json`、`periodic_smoke_validation.json`。相同两个5k ranges的10k/模式原生生产控制1573940四jobs及归约1573950均exit0，结果在 `S/periodic_refit_analysis_v1/rdo`。原native-p同一1758条中近零109→104，变化−0.2844个百分点 [−0.6488,0]；truth-p同一1595条56→56，变化0 [−0.2456,+0.2500]个百分点；全p同一4348条320→318，变化−0.0460 [−0.2930,+0.1661]个百分点。**周期φ差分不能解释大部分近零增强。** 所有off final字段及原生states与原10k基线精确闭合，on的seed/CKF及refit P/F原生state也逐值相同。

独立 full G4/digi 的2-event pipeline smoke1573941全部步骤exit0，随后μ±各2000events三模式生产1573948及归约1573952全部完成。在各自固定native p窗口中，μ+同一1974条近零54→54，变化0 [−0.1514,+0.1511]个百分点；μ−同一1973条51→50，变化−0.0507 [−0.1563,0]个百分点。两模式的效率均为μ+1982/2000、μ−1980/2000；misID counts为1/0、动量负尾counts为4/3，前后分别相同。独立输入同样不支持周期φ差分是主要成因。结果在 `S/periodic_gun_analysis_v1/results`。统计入口为 [slope_periodic_refit.py](../analysis/slope_periodic_refit.py)，已用同一真实20-event数据完成3页分析smoke。

下一最小控制在实际 master Calypso 的独立worktree `geoRotation_calypso_slope_canonical` 中增加 `CanonicalRefitSmoothingAngles=false`，补丁[0012](../patches/0012_refit_canonical_angles.patch)为可直接应用于master的完整替代patch。开启时每次实际smoother更新后调用ACTS已有 `normalizePhiTheta`；θ跨pole反射时同步翻转smoothed covariance的θ行/列，保持该等价chart的一阶误差变换。保留真实gain/covariance计算、scale10、前向拟合和CKF，不截断物理方向。默认off、periodic-only及canonical三组先做20-event闭合，再做固定原输入及独立G4。编译1573951进行中；尚无该控制的物理结果，不称为修复。

### VetoNu 实际材料与完整 G4

真实 Calypso Geant4 world 内，以 G4Navigator 对六个 x/y 位置作只读直线积分，Veto reference→st1 为 250.31–250.33 X0；中央包括 Wolfram 226.436、Lead 17.820、Emulsion 1.378、SSteel 3.957 X0。这是实际 world 的直线材料量，不是每条弯曲散射轨迹的精确累计值。

full G4 μ±×100/250 GeV，每 case 2000 events、总 8000；afterGeo VetoNu residual central68 half-width 如下：

| case | 全 truth 重建效率 | VetoNu half-width |
| --- | ---: | --- |
| μ+ 100 GeV | 53.70% | 2.443 mrad [2.314,2.610] |
| μ+ 250 GeV | 53.25% | 0.9522 mrad [0.9049,1.0163] |
| μ− 100 GeV | 54.15% | 2.489 mrad [2.328,2.621] |
| μ− 250 GeV | 52.70% | 0.9513 mrad [0.9068,1.0118] |

100/250 宽度比 μ+ 为 2.5655 [2.3683,2.8028]，μ− 为 2.6168 [2.3811,2.7846]，与多重散射 1/p 的 2.5 相容。约 250 X0 的 Highland 量级也约为 2.6 mrad@100 GeV。sensor width 为约 22–34 µrad；VetoNu residual 包含真实靶散射，不能称为拟合本征分辨率。这批 gun 的发射位置/方向和生产样本不同，约 54% 效率不代替主样本效率。数据、95% event bootstrap 宽度区间和图在 `S/controls_analysis_v2/gun/`。

### strip coordinate 单因素

实际 Calypso `TrackerClusterMakerTool::UseTruthClusterPosition` 默认 false。开启时根据 cluster RDO 的 SDO majority 与同 wafer SiHits，以 energy-weighted midpoint 和真实 transformHit 转换 local，**仅替换 measured locX**；locY、宽度、误差、噪声和身份继续使用原处理。无相应 truth 的 cluster 保留原坐标并计 fallback。此为诊断 ideal centroid 控制，没有关闭完整 digitization，也不是正式修复。

false 的 20-event smoke 与原 afterGeo observer 数据逐值一致。生产 RDO 对照为 `[0,5000)`、`[250000,255000)`，native 复用已有主重建，true 重新调用相同 Calypso。主结果：

| 指标 | native | truth coordinate |
| --- | ---: | ---: |
| matched pure long | 8617 | 8901 |
| 全 10000 truth 的效率 | 86.16% [85.50,86.87]% | 89.01% [88.41,89.63]% |
| charge misID | 62/8617 = 0.7195% | 31/8901 = 0.3483% |
| p_reco/p_prod<0.5 | 141/8617 = 1.6363% | 75/8901 = 0.8426% |
| station-1 segment 几何谐波幅度 | 0.18847 | 0.02010 |

8604 个共同 particles 中，7763 对完整 cluster identifiers 相同；固定 native p=50–250 GeV 剩 2869 对。st1 harmonic amplitude 0.03879→0.00598，配对变化 −0.03281，95% bootstrap 区间 [−0.04009,+0.00664]，仍含零。VetoNu 的区间也含零。总样本的效率/候选迁移与 same-hit 控制分别报告，不能用 shape 变化掩盖丢失的困难事件。

独立 full G4 在 st1 z=0、x/y±30 mm、slopes±5 mrad 均匀发射 100 GeV μ±，各 2000 events。同一次 G4/digi 的 RDO 各做两种坐标重建。station-1 segment 谐波 μ+ 0.27034→0.01685、μ− 0.26915→0.02887；效率 μ+ 93.85%→99.90%、μ− 94.40%→99.95%。完整同 hits 的 1613/1617 对中，final st1 harmonic 变化 95% 区间分别 [−0.00538,+0.00227]、[−0.00688,+0.00369]，均含零。

`S/cluster_rdo_v2` 将生产RDO配对范围扩展到 `[0,25000)` 和 `[250000,275000)` 共50000 events；原10000 truth-coordinate结果和全量native基线全部复用，只新增40000 events的同一Calypso重建。8个jobs（1573932）均exit0。50k中 native / truth-coordinate 的 matched pure long 为43186/44619，truth效率86.370% [86.076,86.682]% / 89.238% [88.974,89.524]%；misID 为305/43186 / 132/44619，动量负尾为747/43186 / 401/44619。此为各配置自己的总体分母；下面的方向控制另固定同一粒子和完整hit集合。

50k有14380对 same-hit、native p=50–250 GeV tracks，st1 几何谐波幅度0.0297546→0.00637219，配对变化−0.0233824，95% event-bootstrap [−0.0282167,−0.00810037]。为避免最初10k发现样本进入验证，只用新增40k（与原两个5k ranges不重叠）另做独立归约：34504个共同particles，31176对完整hits相同，native-p窗口11511对；谐波0.0275040→0.00663039，变化−0.0208736 [−0.0267572,−0.00461191]。VetoNu对应变化−0.00086197 [−0.00741790,+0.00708119]，区间含零。

冻结主分析 discovery 窗口后，只增加窗口统计而保持上述全部既有字段逐项相同，验证见 `S/centroid_windows_v1/validation.json`。下表为独立40k中的同一11511对tracks，所有窗口在看该对照之前已冻结；变化以**百分点**表示，CI为500次source-event bootstrap：

| st1窗口中心 | native entries | truth-coordinate entries | 占比变化及95%区间 / 百分点 |
| --- | ---: | ---: | ---: |
| −1.89 mrad，半宽0.5 mrad | 1611 | 1619 | +0.0695 [−0.1079,+0.2617] |
| −0.71 mrad，半宽0.107219 mrad | 358 | 329 | −0.2519 [−0.4256,−0.0653] |
| +0.59 mrad，半宽0.1 mrad | 260 | 234 | −0.2259 [−0.3929,−0.0602] |

50k的对应两肩窗口453→417、323→293，与独立40k方向一致。−1.89 mrad宽核窗口及VetoNu两个窗口的变化区间均含零。**证据支持 measured strip coordinate 既产生 segment 周期，又对最终 st1 两肩有可重复的贡献；不将该小幅窗口变化解释为肩峰全部面积或全部物理成因。** 生产50k与独立40k结果位于 `S/centroid_windows_v1/{full50k,independent40k}`，各有两张PNG/逐图PDF及2页图集；原独立G4结果保留于 `S/centroid_analysis_v1/gun/`。

![Independent 40k same-hit centroid control](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/centroid_windows_v1/independent40k/production_muplus_same_hits.png)

## 图表、执行与验证

已交付主样本 57 张 PNG、逐图 PDF、57 页[主图汇编](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/figures_v2/100044/slope_reference_report.pdf)和 index。各图含样本、cuts、frame、reference、entries/invalid/flows；bin counts/sumw2 在 summary 保存。当前主图：

![Physical st1 before/after](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/figures_v2/100044/100044_st1_physical_reference.png)

![Physical VetoNu before/after](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/figures_v2/100044/100044_VetoNu_physical_reference.png)

辅助输入：100116/100117 各 100k 文件的前 10k；110033 全 5 文件 10k；110051 全 20 文件 20k。native 重建与分析全部成功；110051 首次 after JSON 在 EOS 截断，已用原输入和原窗口仅重做归约，完整解析两版 JSON（14.98/14.61 MB），没有重跑重建。`write_json` 改为 worker local 关闭后分块复制并检查 size。`S/figures_v2` 已成功重导五组共 105 张 PNG/逐图 PDF，统一非主面板纵轴、maps 色标，并从恢复后的 110051 counts 导图。105 张 PNG 完整性、105 份单页 PDF 和五份 57/12/12/12/12 页汇编已核对，记录在 `S/figures_v2/file_validation.json`；已目视主样本 truth/residual、肩峰 zoom 和 110051 全量主图。尚待阶段来源表及完整物理验收。

入口为 [slope_tracking.py](../configs/slope_tracking.py)、[run_slope_task.py](../scripts/run_slope_task.py)、[slope_analysis.py](../analysis/slope_analysis.py)、[slope_figures.py](../analysis/slope_figures.py)。原生 writer 归约为 [slope_native_states.py](../scripts/slope_native_states.py)，受控位置替换 patch 为 [0009](../patches/0009_truth_cluster_position.patch)。每次运行的全部命令在对应 task JSON 和日志中；例：

```bash
python3 scripts/submit_momentum.py \
  --tasks /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/native_states_v1/beforeGeo_tasks.json \
  --library /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/software_before/lib \
  --name s05_native_states_before_v1 --stage S05 \
  --cpus 1 --memory-mb 3500 --disk-gb 5 --flavour microcentury \
  --exclude-schedd bigbird21.cern.ch
```

批量全部经 Condor；每次提交重新比较 schedd，选择兼容的最低 running+idle bigbird20，排除已有 EOS 兼容失败的 bigbird21。源任务、schedd 广告/队列、资源和 cluster ID 保存于 `S` 下对应 `*_worker/submission.json`；home 只保存小型 submit。已完成归约包括 production centroid 1573899、独立 gun 分析1573903、图1573910、近零统计1573911、native stage1573922和phase response1573925；最后两组全部exit0。其24张PNG完整性、24份逐图PDF和10/10/2/2页汇编已验证，记录在 `S/native_stage_analysis_v3/file_validation.json`。已目视beforeGeo st1逐阶段图和afterGeo新phase response图。50k及独立40k归约1573933/1573935、冻结窗口1573938、refit covariance归约1573947、周期φ差分生产1573940及归约1573950、独立G4生产1573948、响应分组系统项1573946和独立周期φ分析1573952均exit0；canonical编译1573951正在执行。

基础及 S05 共25个单测通过，另有3个S04 covariance单测同样通过；新增common-validation-support回归确认统一验证覆盖不裁掉训练样本，全部28个测试通过。native writer归约已用两版真实20-event原配置smoke验证，保存字段逐值闭合。新增Python分析的py_compile和git diff --check通过。`S/patch_validation_v1/validation.json` 记录0006–0009可应用且与实际源码对应内容相同。五类共3280个frozen histograms已通过bin counts、sumw2、flows/invalid、共同物理分箱、完整覆盖检查，记录于 `S/figure_counts_validation_v1/results/count_validation.json`。50k/40k新增4张PNG、4份逐图PDF及两份2页图集验证于 `S/centroid_windows_v1/figure_validation.json`，独立40k同hit图已目视。A5最终来源表尚未完成。

本轮refit covariance、周期φ RDO/独立gun、响应系统项共15张PNG、15份逐图PDF及1/1/3/6/2/2页图集均通过完整性/页数检查；114个新control histograms的counts、sumw2、flows/invalid闭合。各响应模型的共同观测逐bin相同，窗口common+uncovered与full counts精确相加，验证位于 `S/refit_response_validation_v1/validation.json`。已目视beforeGeo covariance、production st1周期φ和afterGeo st1响应系统项。

响应分组系统项使用 [slope_response_systematics.py](../analysis/slope_response_systematics.py)：phase bins=4/8/16、phase8细p分组、phase8细x/y分组，五种选择均只用discovery学习response/效率；在其共同covered validation truth上比较相同观测，再把未覆盖窗口entries全部算入未解释余量上界。两个版本任务1573946均exit0，结果位于 `S/response_systematics_v1/{beforeGeo,afterGeo}`。每版本仍以完整500k source events构造discovery/validation，只有条件response的验证支持域取交集；训练样本不裁掉。

afterGeo st1共同覆盖66781/78360条，VetoNu64236/78360条；fine-position分组的稀疏cells是交集缩小的主要来源，不能隐藏其未覆盖部分。下表为afterGeo冻结窗口；CI采用各模型逐窗口event-bootstrap区间的包络，分母均为该plane全部78360条validation reco。"保守余量"再把未覆盖窗口全部entries分配给额外正过量，不将其解释为已经发现的异常：

| plane与窗口中心 | common covered observed | 五模型prediction范围 | covered差值95%包络 / entries | 未覆盖窗口entries | 全validation保守正过量余量 |
| --- | ---: | ---: | ---: | ---: | ---: |
| st1 −1.89 mrad，半宽0.5 | 10732 | 10642.25–10691.29 | [−60.90,+182.86] | 802 | ≤1.2568% |
| st1 −0.71 mrad，半宽0.107219 | 2306 | 2278.38–2353.50 | [−100.74,+106.90] | 246 | ≤0.4504% |
| st1 +0.59 mrad，半宽0.1 | 1378 | 1297.97–1368.39 | [−39.62,+123.92] | 350 | ≤0.6048% |
| VetoNu −0.37 mrad，半宽0.164045 | 3596 | 3562.77–3587.90 | [−102.49,+141.26] | 319 | ≤0.5874% |
| VetoNu +0.43 mrad，半宽0.5 | 10945 | 10946.43–11001.33 | [−223.18,+175.33] | 936 | ≤1.4182% |

这些是声明过的五种分组系统变化和固定窗口内的实测范围，不是全局/同时置信上限，也不覆盖任意更窄的未知峰模型。结合独立40k coordinate控制，它们支持"truth/效率/测量响应的混合可重现主要肩部，strip coordinate提供其中一项周期贡献"；不声称所有未覆盖条件已被解释。beforeGeo对应两肩的保守余量为0.3221%/0.3754%，但它们不是其近零异常的窗口，不用该响应自洽性排除已定位的refit问题。

最终汇编入口为 [slope_final_artifacts.py](../analysis/slope_final_artifacts.py) 和 [图组catalog](../configs/slope_report_catalog.json)，要求25组参考/控制报告及相应证据全部存在后才生成统一PDF和逐图索引；canonical控制尚未完成，当前没有声称最终图集或A5验收通过。

## 阶段性结论与剩余工作

| 结构 | 已测证据 | 当前判断与未完成验收 |
| --- | --- | --- |
| legacy angle 大范围失真、历史 pz 窗口骤减 | 实际 Ntuple 与物理 frame/原宏闭合 | writer/表达式问题已定位；主物理图使用同面正确 slope |
| segment 周期梳状结构、final st1两肩 | full RDO 与独立 G4 的 locX 单因素，独立40k same-hit final谐波及冻结窗口变化 | measured coordinate 来源有独立验证；已量化对两肩的贡献，不归因其全部面积 |
| VetoNu 的 mrad 宽度 | 真实 geometry≈250 X0，full G4 的 100/250 宽度比与 1/p 相符 | 靶散射是宽度主要来源；不是本征 tracker 分辨率 |
| 逐 state shared smoothing 异常 | 恒等式在隔离后恢复，但 final 逐值不变 | 解释逐 state 污染，不解释该批 final 结构 |
| beforeGeo 近零尖峰 | 同候选10k/版本定位到refit smoothing；scale1使固定cohort近零109→67；周期φ仅109→104，独立gun也无主要效果 | 原始covariance/chart对其有作用；仅φ差分解释不足；θ规范化及其C变换控制进行中，程序原因未完成 |
| afterGeo −0.71/+0.59 mrad shoulders | 相同covered分组的phase response；独立40k strip-coordinate配对；五种分组与未覆盖余量 | 支持谱/效率/测量response混合及strip周期贡献；冻结窗保守额外正过量余量为全validation的0.4504%/0.6048%，不声称全条件或任意窄峰均已排除 |

下一步验证canonical角度chart控制的真实Calypso编译、off/on闭合、固定10k和独立full-G4结果；若仍不能解释近零峰，继续检查ACTS→Trk→refit covariance basis，不把本次φ负结果当作终点。补齐最终来源表、统一总PDF/index和A1–A5逐项审计。本报告不将尚未解释的结构转移到后续阶段以宣称S05完成。
