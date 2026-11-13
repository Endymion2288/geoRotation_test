# S05：角度结构定位与实测进展

已完成 100044 的 500,000-event beforeGeo/afterGeo 重建、物理角度与原 writer 闭合、预留样本响应检验及真实 Calypso/G4 成因控制。**历史表达式存在 frame/选样问题；strip coordinate 造成 seed 之前的 segment 周期结构，并对最终 st1 两个肩峰有独立40k事件的配对证据；VetoNu 残差的主要宽度与实际靶材料散射相符。旋转前额外近零结构的来源是近极点方向参数化下的 refit update/smoothing，联合 regular chart 控制在生产RDO及独立μ±中移除主要增强并改善角误差。** 最终结果、剩余余量与效率/pull副作用见[结论报告](2026-11-13_stage05_angle_structure_conclusions.md)，逐项验收见[审计](2026-11-12_stage05_acceptance.md)。本文保留各次试验的发现顺序及冻结结果。

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

为汇报固定主要结构的参数，下面直接摘录 `S/main_summary_v1/afterGeo_summary.json` 的 validation 统计。此表采用参考主图的原生 long tracks、物理 p=50–250 GeV，分母78,428；它没有额外施加 response 的 primary/purity 条件，因此不同于下文78,360条 response cohort。均值和标准差均是在冻结窗口内的条件矩，**不是去背景后的 Gaussian 峰中心/分辨率**；区间为 source-event delta method。窗口占比的区间为1000次 source-event bootstrap。

| 物理面、窗口中心 / mrad | 窗内均值及95%区间 / mrad | 窗内标准差及95%区间 / µrad | 窗口entries / 78428 | 窗口占比及95%区间 / % |
| --- | --- | --- | ---: | --- |
| st1 −1.89，半宽0.5 | −1.87986 [−1.88501,−1.87471] | 282.36 [279.97,284.75] | 11541 | 14.7154 [14.4495,14.9726] |
| st1 −0.71，半宽0.107219 | −0.708446 [−0.710816,−0.706077] | 61.064 [59.987,62.142] | 2552 | 3.25394 [3.12995,3.38204] |
| st1 +0.59，半宽0.1 | +0.590863 [+0.588182,+0.593543] | 56.927 [55.689,58.166] | 1731 | 2.20712 [2.10275,2.31532] |
| VetoNu −0.37，半宽0.164045 | −0.371183 [−0.374134,−0.368232] | 94.213 [92.878,95.547] | 3916 | 4.99311 [4.83518,5.15140] |
| VetoNu +0.43，半宽0.5 | +0.420757 [+0.415591,+0.425922] | 287.40 [285.07,289.72] | 11892 | 15.1630 [14.8945,15.4222] |

beforeGeo的近零结构另使用事先固定的beam-centred窗口：同一primary/purity cohort中，±100µrad占比为6.3210% [6.1604,6.4816]%，afterGeo为2.3635% [2.2571,2.4698]%；不把afterGeo discovery定义的−0.71/+0.59mrad窗口用于替代近零峰。完整参数和分母保存在 `S/zero_analysis_v1/{beforeGeo,afterGeo}/summary.json`。

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

下一最小控制在实际 master Calypso 的独立worktree `geoRotation_calypso_slope_canonical` 中增加 `CanonicalRefitSmoothingAngles=false`，补丁[0012](../patches/0012_refit_canonical_angles.patch)为可直接应用于master的完整替代patch。开启时每次实际smoother更新后调用ACTS已有 `normalizePhiTheta`；θ跨pole反射时同步翻转smoothed covariance的θ行/列，保持该等价chart的一阶误差变换。保留真实gain/covariance计算、scale10、前向拟合和CKF，不截断物理方向。编译1573951已exit0，53个translation units，实际构建diff与仓库补丁完全相同；逐库来源在 `S/canonical_refit_v1/software/software.json`。

三模式20-event检查1573954全部exit0，[validate_refit_control_smoke.py](../scripts/validate_refit_control_smoke.py)核对得到：native的全部events/truth/tracks/segments/native_states与原生基线逐值相同；native与periodic模式还分别精确复现上一轮的25个native_angles字段。canonical的events/truth/segments、80条seed/CKF及38条refit P/F states仍逐值相同，19条refit S中2条改变；312个有限refit measurement states没有超出规范角度范围。检查保留float32存储容差，见 `S/canonical_refit_v1/smoke_validation.json`。

固定原两个5k ranges的canonical生产1573958和归约1573961均exit0。原native comparator复用 `periodic_refit_v1/production/native`，全部events/truth/tracks/segments及native states与最初10k基线精确闭合。固定原候选、majority、primary/purity/long及validation entries≥250000，st1近零结果如下；区间为source-event配对bootstrap：

| 固定选择 | 共同分母 | native→canonical近零数 | 变化及95%区间 / 百分点 |
| --- | ---: | ---: | --- |
| 原native p=50–250 GeV | 1758 | 109→105 | −0.2275 [−0.9090,+0.4116] |
| truth p=50–250 GeV | 1595 | 56→54 | −0.1254 [−0.4657,+0.1964] |
| 不限p | 4348 | 320→313 | −0.1610 [−0.6021,+0.2541] |

所有9707条候选的完整hit集合和majority均不变；38998条seed/CKF与17922条refit P/F保存states逐值不变。原1758 cohort的φ/θ越界由186/154条降至0/0，但大φ差分flag由54条升至213条，其中55条近零。规范化等价chart消除了范围违规，不能据此声称物理角度结构消失。全部10k的matched pure long8680、efficiency8680/10000、wrong charge317/8680和动量负尾773/8680均相同；VetoNu传播失败由0变1，其余两面为0。完整结果为 `S/canonical_refit_analysis_v1/rdo/summary.json`。

独立full-G4 μ±各2000events使用50681/50682新随机种子，区别于上一轮50581/50582。pipeline smoke1573956、生产1573960和分析1573963均exit0。实际完整G4/digitization不变，在同一次RDO上比较beforeGeo native/canonical，afterGeo单列参考。固定原native-p cohort的st1近零：μ+为46→49/1976，变化+0.1518 [0,+0.3427]个百分点；μ−为53→50/1972，变化−0.1521 [−0.3943,+0.0511]个百分点。两电荷VetoNu近零数分别为41/1976及23/1972，控制前后不变。μ+效率1980/2000、wrong charge1/1980、p-tail1/1980；μ−为1975/2000、0/1975、2/1975，均不变，三面传播失败均0。完整结果为 `S/canonical_gun_analysis_v1/results/summary.json`。因此φ周期差分及θ规范化不能解释beforeGeo主要近零增强。

下一控制使用真实master worktree `geoRotation_calypso_slope_basis`、相同commit/条件、patches [0013](../patches/0013_refit_covariance_basis.patch)和[0014](../patches/0014_native_refit_outcomes.patch)。0013分别开启ACTS surface→Trk curvilinear→refit plane的空间local covariance basis变换，以及将refit reference从原input.z−10mm移回input.z；保留native scale10及真实ACTS拟合。basis-only仍保留原mean displacement，作为独立敏感性；basis+plane才在同点使用两次basis变换。本补丁未修正raw theta反射时的angular covariance，不能将空间控制解释为完整角度chart修复。0014只记录实际fit返回的success/fallback，按source entry/run/event/ordinal与observer逐项核对。

编译1573962 exit0，53个translation units；实际source diff、冻结source.patch完全一致，18项software overlay来源在 `S/basis_refit_v1/software/software.json`。五模式20-event smoke1573964全部exit0：native与原始基线逐字段相同，observed与native的events/truth/tracks/segments/native_states/native_angles全部相同；所有控制保持80条seed/CKF states不变。每模式实际20条selected tracks中refit成功19、fallback1，观察身份全闭合。四模式各10k生产1573965及三组归约1573968均exit0，原生observer全部字段与最初10k基线闭合。

| 空间接口控制 | native-p近零数/1758 | 变化及95%区间 / 百分点 | truth-p近零数/1595 | all-p近零数/4348 |
| --- | ---: | --- | ---: | ---: |
| native（仅观察） | 109 | — | 56 | 320 |
| basis | 109 | 0，配对计数不变 | 56 | 320 |
| input plane | 101 | −0.4551 [−1.0592,+0.1098] | 51 | 303 |
| basis + input plane | 101 | −0.4551 [−1.0592,+0.1098] | 51 | 303 |

basis的9707个完整hit集合及majority不变；plane/basis+plane有9671/9707个完整hit集合不变，majority仍全部相同。固定same-hit native-p子组1753条，近零仍109→101。truth-p控制变化−0.3135 [−0.7685,+0.1258]个百分点；all-p为−0.3910 [−0.7416,−0.0235]个百分点。空间basis在所测近零计数上没有作用，起点修正有小幅敏感性，仍保留绝大部分增强。

实际refit返回值排除了一项具体解释：原9707条候选中成功9110、fallback597；原native-p 1758中成功1753、fallback5；**109条近零轨迹全部refit成功，fallback为0**。plane下101条近零同样全部成功，所有模式的成功/fallback总数一致。这里的0是该固定样本的实测计数，不将bootstrap的退化[0,0]当作总体误率上限；全500k未保存该观察字段，不作全量fallback归因。native/basis的matched pure long8680、wrong charge317、p-tail773；plane/basis+plane为8679、326、796。各自分母分别8680/8679，效率分母10000，三面传播失败均0。各组完整计数、source-event区间和9张PNG/逐图PDF见 `S/basis_refit_analysis_v1`。

独立full-G4 basis/plane控制使用新种子50781/50782、μ±各2000。2-event smoke1573967七个实际步骤全部exit0，同一digitized RDO各控制的events/truth/segments精确相同，成功/fallback身份闭合；生产1573973及归约1573980均exit0。固定native-p cohort，空间basis使st1近零数μ+46→46/1965、μ−58→58/1968；input plane与basis+plane均为μ+46→50/1965（+0.2036 [−0.0998,+0.4760]个百分点）、μ−58→58/1968（0 [−0.2484,+0.2527]个百分点）。basis完整hit集合2000/2000不变，两个plane模式各1993/2000不变，majority全部不变。native/basis效率为1967/2000、1972/2000，wrong charge为0/1967、0/1972，动量负尾1/1967、3/1972；plane两模式的效率不变，μ+两项不变，μ−wrong charge1/1972、动量负尾1/1972，三面失败均0。该独立输入同样不支持空间basis或起点为主要近零原因。18张PNG/逐图PDF、3份6页报告和216个histogram的counts/sumw2/flow/invalid闭合，见 `S/basis_gun_analysis_v2/artifact_validation.json`；μ+ st1 plane图已目视。

后续线索在真实20-event native writer中可见：CKF filtered有4/317个states的θ越界，refit filtered为4/312，因此进一步验证ACTS32 update输出mean/covariance到下一次transport之前的chart一致性。这不同于前面的smoother-only规范化。实际master worktree `geoRotation_calypso_slope_update` 已保存完整替代patch [0015](../patches/0015_refit_canonical_update.patch)，新增默认false的 `CanonicalRefitUpdateAngles`：先调用原生 `GainMatrixUpdater`，仅将其成功后的filtered mean和theta covariance一起转到等价规范chart，再交还实际Calypso传播。CKF保持原配置，smoother独立开关；不以此代替新的tracking实现。编译1573972实际exit0、53个translation units。四模式20-event smoke1573977均exit0：native及canonical分别精确复现原基线和上一轮smoother-only的所有保存字段；80条seed/CKF rows不变，新控制312个有限refit filtered states的角度范围违规为0，见 `S/update_refit_v1/smoke_validation.json`。updater/update_canonical各10k生产1573981四jobs及分析1573984均exit0。独立full-G4新种子50881/50882的2-event七步骤smoke1573982已exit0，source字段、加载库、开关及outcome身份由 [validate_slope_gun.py](../scripts/validate_slope_gun.py)核对；μ±各2000生产1573983及三模式分析1573986均exit0。

### forward-update 结果与方向 chart 控制

固定原native candidate/majority/p窗口及validation范围，两个forward控制仍保留主要近零增强：

| 模式 | st1 native-p近零 /1758 | 配对变化及95%区间 / 百分点 | truth-p近零 /1595 | all-p近零 /4348 |
| --- | ---: | --- | ---: | ---: |
| native | 109 | 参考 | 56 | 320 |
| updater | 102 | −0.3982 [−1.2095,+0.3516] | 54 | 299 |
| updater + canonical smoother | 105 | −0.2275 [−0.9811,+0.5566] | 53 | 309 |

两个控制9694/9707条候选的完整final-hit集合不变，majority全部不变，实际success/fallback仍9110/597；控制后的102/105条近零全部refit成功。seed/CKF fields保持不变。matched pure long8680→8681；wrong charge317/8680→295/8681、294/8681，动量负尾773/8680→731/8681；VetoNu传播失败0→1，另外两面仍0。近零结构没有随这些伴随指标的变化消失。

独立新种子50881/50882、μ±各2000的st1固定native-p结果如下；分母分别1974/1975：

| 模式 | μ+近零 | μ−近零 | μ+配对变化95%区间 / 百分点 | μ−配对变化95%区间 / 百分点 |
| --- | ---: | ---: | --- | --- |
| native | 48 | 46 | 参考 | 参考 |
| canonical smoother | 49 | 54 | +0.0507 [−0.1977,+0.2584] | +0.4051 [+0.1516,+0.7039] |
| updater | 58 | 52 | +0.5066 [+0.1024,+0.9472] | +0.3038 [0,+0.6569] |
| updater + canonical smoother | 50 | 48 | +0.1013 [−0.1510,+0.3527] | +0.1013 [−0.1539,+0.4017] |

独立μ+的效率1977/2000、wrong charge2/1977、动量负尾0/1977在三个控制中均不变；μ−native/canonical为1978/2000、1/1978、2/1978，两个updater模式为1979/2000、0/1979、1/1979。canonical完整hits全部不变，两个updater模式μ+2000/2000、μ−1998/2000相同，所有majority不变、三面失败均0。此处零计数不作为总体率上限。24张PNG/逐图PDF、2份3页RDO报告和3份6页gun报告，以及288个histogram counts/sumw2/flows均通过检查，分别记录在 `S/update_refit_analysis_v2/artifact_validation.json` 和 `S/update_gun_analysis_v1/artifact_validation.json`。因此forward canonicalization也不是主要近零增强的解释。

下一试验区分“角度范围规范化”和“在近极点chart中线性化”。原生ACTS32 smoother在参数向量上作差；即使phi周期、theta反射均规范，near-pole的phi仍会随微小物理方向变化产生很大变化。实际master worktree `geoRotation_calypso_slope_chart_centered` 使用完整替代补丁 [0016](../patches/0016_refit_regular_direction_chart.patch)，新增默认false的 `RegularRefitSmoothingChart`，在真实refit smoother内部将方向改写为A-frame longitude与latitude；geometry、B/material、measurements、native forward fit保持原配置。

每个实际RTS相邻states分别使用其filtered、predicted和smoothed mean处的解析Jacobian：`C_f′=J_f C_f J_fᵀ`、`C_p′=J_p C_p J_pᵀ`、`C_s′=J_s C_s J_sᵀ`、`F′=J_p F J_f⁻¹`。保留全部position/angular/q/p/time交叉项，由实际 `Acts::GainMatrixSmoother::calculate`完成gain与后向更新，再将结果及covariance变回native chart。这是在Calypso中对数值参数化的控制，不是另一套transport或tracking。

数学验证直接include实际Calypso的 `RefitDirectionChart.h`，输入1887组真实native-writer角度，见 [C++检查](../configs/refit_direction_chart_check.cxx) 和 [入口](../scripts/validate_refit_direction_chart.py)。最初使用接近π/2的polar角，在逆Jacobian检查中得到4.03×10⁻⁶误差；对应冻结v1未用于物理结论。改为以0为中心的latitude后，方向往返最大误差6.94×10⁻¹⁸、逆Jacobian误差2.18×10⁻¹¹、中心差分scaled误差3.60×10⁻¹¹，1887组全部通过，记录 `S/chart_math_validation_v3/validation.json`。第一次本地编译未加载GCC13的libstdc++，checker plugin报GLIBCXX_3.4.31缺失；v2/v3显式采用同编译器runtime，失败日志未覆盖。

ACTS的实现与参数/曲面设计参考 [A Common Tracking Software Project](https://arxiv.org/abs/2106.13593)；使用局部chart检验非线性方向变量的依据参考 [Integrating Generic Sensor Fusion Algorithms with Sound State Representations through Encapsulation of Manifolds](https://arxiv.org/abs/1107.1119)。这些文献支持数值试验的设计，不是本样本原因的证据。中心化版本实际编译1573988成功，53个translation units；20-event smoke1573993均exit0。native全部events/truth/tracks/segments/native_states/native_angles/refit_outcomes精确复现旧基线；80条seed/CKF与38条refit P/F rows逐值不变。记录为 `S/chart_refit_v2/smoke_validation.json`。

### regular direction-chart：实测作用与剩余质量问题

固定10k生产1573995及分析1573997全部exit0。全部9707条候选的完整hits和majority不变；38998条seed/CKF与17922条refit P/F rows逐值不变。下面仅改变实际refit后向smoother的方向chart；近零窗仍为预先固定的±100µrad，选样固定native的primary muon、三站、purity≥0.7及validation entry≥250000。

| 同面与固定选样 | 分母 | native近零 | regular chart近零 | 同面truth近零 | control−native及95%区间 / 百分点 |
| --- | ---: | ---: | ---: | ---: | --- |
| st1，native-p | 1758 | 109 | 43 | 58 | −3.7543 [−4.8711,−2.7087] |
| st1，truth-p | 1595 | 56 | 32 | 30 | −1.5047 [−2.2743,−0.7205] |
| sensor，native-p | 1758 | 110 | 41 | 58 | −3.9249 [−5.0141,−2.8993] |
| sensor，truth-p | 1595 | 56 | 30 | 30 | −1.6301 [−2.3780,−0.8280] |

st1 truth-p子组的control−truth为+0.1254 [−0.5563,+0.7881]个百分点；这是此固定1595条cohort的局部剩余差值，不外推为全500k或任意峰宽的上限。全部9707候选的actual refit success/fallback均为9110/597；native-p1758中1753/5，控制后的43条近零全部refit成功。效率8680/10000、wrong charge317/8680、动量左尾773/8680及三面传播失败0均不变。st1全p近零320→199/4348，VetoNu在固定native-p中39→38/1758。完整结果与54个histograms的counts/sumw2/flows、3张PNG/逐图PDF检查为 `S/chart_refit_analysis_v1/{rdo/summary,artifact_validation}.json`。

独立full-G4种子50981/50982、μ±各4000，100GeV、st1 z_F=0、x/y±30mm、theta_x/y输入均匀±0.5mrad，生产1573996及分析1573999全部exit0。st1参考直接用该面生成动量，不再用首sensor传播model。μ+固定native-p近零1299→686/3906，精确生成truth783；μ−1321→660/3912，truth749。完整hits两电荷均4000/4000不变；效率分别3920/4000、3935/4000，wrong charge1/3920、1/3935，动量左尾13/3920、20/3935，三面失败0，控制前后均相同。所有4000候选的refit成功/回退每电荷3965/35；近零候选均成功。108个histograms、6张PNG/逐图PDF验证在 `S/chart_gun_analysis_v1/artifact_validation.json`。

这一控制在同输入及独立两电荷样本上显著移除主要额外近零结构，支持将其来源定位到近beam极点下的refit smoothing方向参数化。但它尚未通过整体拟合质量验收：[slope_chart_diagnostics.py](../analysis/slope_chart_diagnostics.py)的1574000在相同cohort上测得，生产st1 residual RMS从0.46667增至1.49993mrad，absolute-error68从108.71降至100.06µrad；独立μ+ RMS从0.29319增至0.30528mrad、μ−0.28731增至0.31804mrad。独立两电荷control−generated近零仍为−2.4834及−2.2751个百分点，区间均低于0。pull亦未达到标准normal校准。因此当前补丁是成因控制，尚不能作为完整生产修复；下一步检查原生前向update在polar chart中的线性化与后向regular chart的一致性。不能以近零数下降代替残差尾部验收。

同批新gun的afterGeo参考出现0条long tracks。其segments计数与beforeGeo完全相同（μ+12216、μ−12201），seed总数4001/4000，损失发生于CKF。已查到启动路径的程序原因：outer worker使用beforeGeo overlay，旧 `run_slope_simulation.py`在afterGeo子进程只prepend其局部库，保留了beforeGeo fallback；而旧 `software_after`只提供KalmanFilter库，geometry会从父beforeGeo目录解析。已在 [driver](../scripts/run_slope_simulation.py)隔离两个版本的库/Python搜索路径，在 [slope_tracking.py](../configs/slope_tracking.py)新增实际`/proc/self/maps`记录，并冻结完整afterGeo geometry overlay `S/software_after_complete_v1`。同20-event old/clean启动方式对照1574001全部exit0：旧方式实际加载master的两份geometry库，0条track；清理搜索路径后实际加载rotated geometry，20条track。生成粒子、SiHit truth、segments及seed covariance逐字段相同，记录为 `S/loader_audit_v1/validation.json`。beforeGeo native/regular控制使用完整同版overlay，其上述配对数字不依赖这条afterGeo参考。受影响的旧nested gun afterGeo参考不作为物理结果；主500k两版分别按其overlay提交，不属于这种子进程跨版本fallback。

随后在真实master worktree `geoRotation_calypso_slope_chart_update` 增加独立开关 `RegularRefitUpdateChart`，完整替代patch为 [0017](../patches/0017_refit_regular_update_chart.patch)。开启时先将predicted mean及完整6×6 covariance转到同一regular longitude/latitude chart，调用实际 `Acts::GainMatrixUpdater`，再把filtered mean/covariance变回native chart交给原生transport；恢复predicted slots，保留原measurement和projector。与regular smoother构成native、smoother-only、update-only、both四个模式。实际编译1574002已exit0，53个translation units与Gaudi metadata产出成功；18项软件来源在 `S/chart_update_refit_v1/software/software.json`。四模式20-event smoke1574003通过：native关闭逐字段不变，smoother-only精确复现上一版本；实际加载库、seed/CKF及refit outcome身份闭合。10k生产1574005和分析1574007、独立μ±各4k的1574004/1574006及分析1574010均exit0。

联合regular update/smoothing的validation近零为61/1758，同面truth58；truth-p窗口为31/1595、truth30。固定same-hit子组仍为108→61/1744、truth58，残差RMS为0.46156→0.28220mrad。联合控制在固定原native-p cohort上的RMS为0.46667→0.32907mrad，绝对误差68%分位108.71→41.76µrad，超过1mrad的误差38→8；单独改update或smoother反而增大大尾。全部9707候选身份/majority不变，38998条seed/CKF rows逐值相同；76条完整final hit sets变化，actual成功/fallback9110/597不变。

严格成套加载的独立gun验证在 `S/chart_gun_v2`，每个电荷native与smoother-only逐字段精确复现原有效beforeGeo结果，修正后的afterGeo均4000条tracks。联合控制使固定native-p近零μ+1299→765/3906（truth783）、μ−1321→753/3912（truth749）；both−truth的source-event Poisson bootstrap95为−0.4608 [−1.1779,+0.2753]和+0.1022 [−0.6480,+0.8544]个百分点，same-hit子组的剩余差也与0相容。残差RMS降至0.10404/0.09288mrad。pure-long效率3920→3840、3935→3858/4000，下降2.000/1.925个百分点；pull RMS仍1.373/1.331，生产RDO为2.03。成因验证成立，生产效率优化和误差校准尚未测量，不能据此宣布完整生产修复通过。完整数字、同面truth定义和区间口径见[最终结论](2026-11-13_stage05_angle_structure_conclusions.md)。

最终展示重导1573998的15个jobs全部exit0；没有重跑任何tracking。`S/refit_display_v2/validation.json`核对原843452个数值字段全部精确相同，1104个histograms的counts/sumw2/flows/invalid闭合，69张PNG可完整解码，69份单页PDF与15份报告逐页文本一致。full、zoom、residual每个绘图面板均保存计数。已目视regular RDO/gun及重导updater sensor图。新增2个环境回归测试通过，确认跨版本fallback不会残留；原38个仓库测试亦通过。

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

辅助输入：100116/100117 各 100k 文件的前 10k；110033 全 5 文件 10k；110051 全 20 文件 20k。native 重建与分析全部成功；110051 首次 after JSON 在 EOS 截断，已用原输入和原窗口仅重做归约，完整解析两版 JSON（14.98/14.61 MB），没有重跑重建。`write_json` 改为 worker local 关闭后分块复制并检查 size。`S/figures_v2` 已成功重导五组共 105 张 PNG/逐图 PDF，统一非主面板纵轴、maps 色标，并从恢复后的 110051 counts 导图。105 张 PNG 完整性、105 份单页 PDF 和五份 57/12/12/12/12 页汇编已核对，记录在 `S/figures_v2/file_validation.json`；已目视主样本 truth/residual、肩峰 zoom 和 110051 全量主图。后续来源表与物理验收见最终结论报告及A1–A5审计。

入口为 [slope_tracking.py](../configs/slope_tracking.py)、[run_slope_task.py](../scripts/run_slope_task.py)、[slope_analysis.py](../analysis/slope_analysis.py)、[slope_figures.py](../analysis/slope_figures.py)。原生 writer 归约为 [slope_native_states.py](../scripts/slope_native_states.py)，受控位置替换 patch 为 [0009](../patches/0009_truth_cluster_position.patch)。每次运行的全部命令在对应 task JSON 和日志中；例：

```bash
python3 scripts/submit_momentum.py \
  --tasks /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/native_states_v1/beforeGeo_tasks.json \
  --library /eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/software_before/lib \
  --name s05_native_states_before_v1 --stage S05 \
  --cpus 1 --memory-mb 3500 --disk-gb 5 --flavour microcentury \
  --exclude-schedd bigbird21.cern.ch
```

批量全部经 Condor；每次提交重新比较 schedd，选择兼容的最低 running+idle bigbird20，排除已有 EOS 兼容失败的 bigbird21。源任务、schedd 广告/队列、资源和 cluster ID 保存于 `S` 下对应 `*_worker/submission.json`；home 只保存小型 submit。已完成归约包括 production centroid 1573899、独立 gun 分析1573903、图1573910、近零统计1573911、native stage1573922和phase response1573925；最后两组全部exit0。其24张PNG完整性、24份逐图PDF和10/10/2/2页汇编已验证，记录在 `S/native_stage_analysis_v3/file_validation.json`。已目视beforeGeo st1逐阶段图和afterGeo新phase response图。50k及独立40k归约1573933/1573935、冻结窗口1573938、refit covariance归约1573947、周期φ差分生产1573940及归约1573950、独立G4生产1573948、响应分组系统项1573946和独立周期φ分析1573952均exit0；canonical编译1573951及三模式smoke1573954亦已exit0。

基础及 S05 共25个单测通过，另有3个S04 covariance单测同样通过；新增common-validation-support回归确认统一验证覆盖不裁掉训练样本，全部28个测试通过。native writer归约已用两版真实20-event原配置smoke验证，保存字段逐值闭合。新增Python分析的py_compile和git diff --check通过。`S/patch_validation_v1/validation.json` 记录0006–0009可应用且与实际源码对应内容相同。五类共3280个frozen histograms已通过bin counts、sumw2、flows/invalid、共同物理分箱、完整覆盖检查，记录于 `S/figure_counts_validation_v1/results/count_validation.json`。50k/40k新增4张PNG、4份逐图PDF及两份2页图集验证于 `S/centroid_windows_v1/figure_validation.json`，独立40k同hit图已目视。A5最终来源表已在结论报告汇总。

本轮refit covariance、周期φ RDO/独立gun、响应系统项共15张PNG、15份逐图PDF及1/1/3/6/2/2页图集均通过完整性/页数检查；114个新control histograms的counts、sumw2、flows/invalid闭合。各响应模型的共同观测逐bin相同，窗口common+uncovered与full counts精确相加，验证位于 `S/refit_response_validation_v1/validation.json`。已目视beforeGeo covariance、production st1周期φ和afterGeo st1响应系统项。

[slope_response_flows.py](../analysis/slope_response_flows.py)从完整输入表重放冻结RNG的8次response抽样，验证nominal phase及五种系统分组，不重新重建或做额外拟合。1573953两jobs均exit0：两版×两面×六模型共24组prediction的逐bin counts、Monte Carlo sumw2、under/overflow及expected total全部与冻结值闭合；完整/covered/uncovered的observed counts、sumw2和flows也逐项相加闭合。计数遵循各原算法的右边界规则，本批没有恰落右端导致歧义的entries。数据位于 `S/response_flow_validation_v1/{beforeGeo,afterGeo}/validation.json`，补足了response图的flow sumw2证据。

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

## A′ 双track拓扑的条件分布

[slope_topology.py](../analysis/slope_topology.py)复用两版真实Calypso已保存的全量结果：110033的5文件/10000源事件，110051的20文件/20000源事件。以输入file+entry身份的每个事件配对status=1的e+/e−，要求共同production vertex；所有10000/20000事件均得到唯一的两粒子配对，错误multiplicity/charge pair/vertex均0。原生truth中分别211/28个事件缺少至少一个粒子的station-1 SiHit centroid，相关条件列为missing，不将其填0。

三个因素独立分组：共同production z；产生点两动量的 `atan2(|p1×p2|,p1·p2)` 张角；实际station-1 SiHit平均位置的横向两粒子间距。最后一项是同station的truth hit中心间距，没有将其称为共同z外推距离，也没有重算任何tracking/transport。重建选择为primary e±、三站long、SDO purity≥0.7，保留重复候选；分别保存all-p、reco-p=50–100GeV及truth-p=50–100GeV的st1/VetoNu/sensor truth、reco、residual、近零/p-tail/charge/失败计数。图采用reco-p窗口，归一化以各组有限entries为分母，单位mrad；truth仍为明确声明的SiHit-chord/production-p传播模型。

以st1 `|sy|<100µrad`为例，以下每格是近零条数/共同有限分母。区间使用source-event cluster delta method，保持同事件两粒子和重复候选的相关性；各版按自己的reco-p窗口选样，因此before/after列不作同候选的因果差值：

| ST1 truth centroid间距 / mm | 110033 before | 110033 after | 110051 before | 110051 after |
| --- | ---: | ---: | ---: | ---: |
| [0,0.1) | 195/436 | 134/376 | 179/578 | 91/337 |
| [0.1,0.5) | 93/565 | 39/515 | 129/823 | 32/477 |
| [0.5,2) | 35/1041 | 10/935 | 13/525 | 7/429 |
| [2,∞) | 15/1374 | 4/1337 | 5/135 | 1/120 |

afterGeo 110051小间距组为27.003% [22.272%,31.734%]，大间距组0.833% [0,2.460%]；110033相应为35.638% [30.760%,40.517%]和0.299% [0.006%,0.592%]。同图的truth本身已在小间距组呈现更窄的st1中心分布，reco保留这个趋势；这支持将拓扑与输入谱的混合纳入双track角图解释，不把上述近零条数全部认定为algorithm peak。

衰变位置也有明显关联：afterGeo 110051在z<−1000、[−1000,−500)、[−500,0)mm的近零分别16/421、28/478、87/464；110033为5/1091、8/1065、174/1035。该关联同时改变输入方向/间距和接受度，没有以单因素描述性分组冒充独立物理消融。opening方面，110033使用0/0.05/0.1/0.25/∞ mrad；110051首次0.05mrad分界以上仅6条reco-p tracks，细分检查的0–5µrad也仅4条，因此合并为0/10/20/∞µrad，避免解释稀疏曲线形状。探索和合并顺序在v1/v2/v3 manifest完整记录，不作为预留样本的显著性检验。

首轮分组1573966、细化1573969均exit0，110051稀疏组的最终合并1573976已exit0。110033正式结果为 `S/topology_analysis_v2/110033`，110051最终入口为 `S/topology_analysis_v3/110051`。v2四组12张PNG/逐图PDF及1296个histogram counts、sumw2、flow/invalid已核对，已目视110051间距条件图。新增三个回归检查验证file-entry身份下的配对、opening单位、缺hit不变成零间距和共同vertex要求，均通过。v3三张PNG/逐图PDF及每版297个histogram已核对，合并后的opening bins精确相加，其他两个条件组与v2逐项相同。afterGeo的opening PNG出现IDAT CRC损坏；从原封存PDF渲染恢复，不重算统计，也不修改原产物；恢复图及完整展示入口在 `S/topology_analysis_v3/export_recovery/afterGeo`，来源记录为同级 `recovery.json`，修复图已目视。后续新增图的PNG/逐图PDF改为worker scratch完整关闭后复制至EOS并verify PNG，旧冻结输出保持原样。

最终汇编入口为 [slope_final_artifacts.py](../analysis/slope_final_artifacts.py) 和 [图组catalog](../configs/slope_report_catalog.json)，列出52组参考/控制报告及75项证据；regular update/both、严格成套加载的独立gun及18张最终质量图已补入。旧图15组、69张的单位及panel counts重导已完成；最新质量展示6组、498个histograms、18张图从冻结counts重导，494991个原数值字段精确复现。原始数组重放及count/sumw2/flows见 `S/chart_diagnostics_display_v1/validation.json`。最终汇编与逐页核对结果记录在[结论报告](2026-11-13_stage05_angle_structure_conclusions.md)。

## 阶段性结论与剩余工作

| 结构 | 已测证据 | 当前判断与未完成验收 |
| --- | --- | --- |
| legacy angle 大范围失真、历史 pz 窗口骤减 | 实际 Ntuple 与物理 frame/原宏闭合 | writer/表达式问题已定位；主物理图使用同面正确 slope |
| segment 周期梳状结构、final st1两肩 | full RDO 与独立 G4 的 locX 单因素，独立40k same-hit final谐波及冻结窗口变化 | measured coordinate 来源有独立验证；已量化对两肩的贡献，不归因其全部面积 |
| VetoNu 的 mrad 宽度 | 真实 geometry≈250 X0，full G4 的 100/250 宽度比与 1/p 相符 | 靶散射是宽度主要来源；不是本征 tracker 分辨率 |
| 逐 state shared smoothing 异常 | 恒等式在隔离后恢复，但 final 逐值不变 | 解释逐 state 污染，不解释该批 final 结构 |
| beforeGeo 近零尖峰 | 同候选10k/版本定位到refit smoothing；周期φ、canonical、空间basis/input plane不足以解释；联合regular update/smoothing近零109→61、truth58；独立μ+1299→765、truth783，μ−1321→753、truth749；same-hit子组重复改善 | 近极点方向参数化的因果证据成立；所列cohort/window的剩余差与0相容，角误差核心和大尾同时改善。效率副作用及pull未校准另报，未外推全500k修复验收 |
| afterGeo −0.71/+0.59 mrad shoulders | 相同covered分组的phase response；独立40k strip-coordinate配对；五种分组与未覆盖余量 | 支持谱/效率/测量response混合及strip周期贡献；冻结窗保守额外正过量余量为全validation的0.4504%/0.6048%，不声称全条件或任意窄峰均已排除 |
| A′双track的近零结构 | 全110033/110051按decay z/opening/ST1 truth centroid间距；小间距组truth与reco均较窄 | 明确的输入拓扑/接受度关联；不把此分组结果当作μ+ refit增强的解释 |

空间接口独立full-G4、拓扑合并、canonical及regular forward-update/smoothing的固定10k和独立输入均已完成，nested gun的afterGeo参考已经更正。主要结构来源、固定窗口余量与副作用已汇总至[结论报告](2026-11-13_stage05_angle_structure_conclusions.md)。后续生产修复须在新保留输入上验收效率、pull及完整tracking指标；本阶段的正控制用于定位成因。
