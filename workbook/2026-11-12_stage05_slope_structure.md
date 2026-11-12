# S05：角度结构定位与实测进展

本阶段尚未通过 A1–A5 完整验收。已完成 100044 的 500,000-event beforeGeo/afterGeo 重建、物理角度与原 writer 闭合、预留样本响应检验及真实 Calypso/G4 控制。**历史表达式存在 frame/选样问题；segment 在 seed 之前已有 strip 周期结构；VetoNu 残差的主要宽度与实际靶材料散射相符。旋转前近零尖峰及旋转后剩余肩峰仍需逐阶段和最小控制完成定位。**

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

±100 µrad 的 st1 窗口中，beforeGeo reco-p 选样为 5574/88182=6.3210%，相同选中粒子的 truth 为 3222/88182=3.6538%；truth-p 选样仍为 reco 2413/78940=3.0568%、truth 1727/78940=2.1877%。因此 p 选样迁移不能单独解释近零尖峰。phase response 对 beforeGeo 预测 4287.96 条，实际多出 1286.04 [1106.37,1441.49] 条；这个差值是该条件模型的验证偏差，仍需处理未覆盖部分和模型系统项。afterGeo 相同 reco-p 窗口为 1852/78360=2.3635%、truth 为 1854/78360=2.3660%；phase response 预测 1809.80 条，差 42.20 [−36.49,109.14]，与零相容。sensor 的 truth-p 对照同样保留 beforeGeo 的近零过量，排除仅 fixed-plane 外推制造结构的解释。

经验 response 使用 discovery 的 residual 和 track multiplicity，按 production p、同面 truth slope、sensor x/y 分组，应用于 held-out truth。原 coarse response 的 afterGeo st1 reco-p shape TV 为 2.4413%，VetoNu 为 1.3566%；低响应分组的覆盖不足单列，不能把全部结构称为已解释。

实际 geometry 为 strip pitch=0.08 mm、三层 lever arm=63 mm、stereo≈±0.01999867 rad，几何周期 `pitch/lever_arm=1.269841 mrad`。`S/phase_response_v1` 使用此固定周期的 8 phase bins，没有自由扫描周期或拟合 Gaussian 分量。afterGeo held-out st1 p_reco 窗口如下：

| 窗口中心 | observed | phase prediction | 差值和 source-event bootstrap 95% 区间 |
| --- | ---: | ---: | --- |
| −1.89 mrad | 11534 | 11333.32 | 200.68 [64.57,302.98] |
| −0.71 mrad | 2552 | 2573.00 | −21.00 [−88.55,51.61] |
| +0.59 mrad | 1728 | 1655.44 | 72.56 [8.70,134.02] |

st1 shape TV=2.3462%；VetoNu=1.6119%。该 st1 响应未覆盖 96926/244741 个 validation truth，对应 2486/78360 条实际选中 reco；VetoNu 未覆盖 92522/244636 truth、1549/78360 reco。未覆盖 reco counts 已逐 bin 保存。+0.59 mrad 的差值不能直接换算成全样本额外峰的严格上限，也未做多窗口 familywise 显著性声明。−0.71 mrad 结构与几何 phase 条件响应相容，支持测量/响应相关来源；经验模型仍从同一模拟重建中学习，不能据此断言全部为物理效应。

## A4：真实 Calypso/G4 控制

### 自由方向与逐阶段审计

每版本 294 个等价 free states，扫描正负近零 slopes、50/100/1000 GeV、q±，传播到 z_F=±100 mm：各 588 次成功、0 失败。before/after physical mean 差最大约 2e−13 GeV、1e−14 mm，physical Var(sy) 比在 1±3e−15 内。before 的 6 个 exact-pole states 无定义的 phi covariance 单列，只传播 mean。结果证明此次自由传播数值闭合，不证明整个近极点 CKF 更新正常或两版重建 covariance 相同。

afterGeo 已审计 4000 events（`[0,2000)` 与 `[250000,252000)`），3452 条 matched pure long tracks。从真实 seed、CKF P/F/S、converted、refit、final 经实际 Calypso 投到 st1/VetoNu/sensor，final projection 与 observer 最大误差为 0。结果位于 `S/stage_analysis_v2/results/`。部分窄窗峰与匹配 sideband 统计少，不能据少数事件终止该分支。

另做 2000-event 仅切换 candidate-isolated smoothing 的配对：P/F、seed、selected target、converted、refit、final 的 means 均逐值不变。共享 slots 的角度 smoother 恒等式超阈值项在两 shards 为 1/7，隔离后均为 0；final mean、charge、hits、projection 全部不变，Good95 不变。故它可污染逐 state 诊断，但不解释这批最终角度结构。详见 `S/checks_analysis_v1/smoothing/summary.json`。

为定位 beforeGeo 的近零尖峰，现已开启其本来就有的 `RootTrajectoryStatesWriterTool`（CKF/refit）和 `RootSeedWriterTool`，无需另一套 tracking 或重编 master。20-event/版本 smoke 全部成功；各自的 5/15/30/6 个保存 event/truth/track/segment branches 与原基线逐值相同。refit 分别 19/15 条，均以同一事件、完整 module/local measurement signature 唯一匹配，无缺失或歧义。源事件身份还以 observer counts 和 CKF 顺序交叉核对。`S/native_states_v1` 每版 10000 events 均已 exit0：beforeGeo refit 9110 条中 149 条未匹配，afterGeo 8928 条中 65 条未匹配，均无歧义；未匹配者不进入 same-hit stage 对照，峰区的覆盖率需单列。seed 暂按 population 保存，未伪称已关联同一 ancestor seed。`S/native_transport_v2` 正在用实际 Calypso 同面传播，完整逐阶段物理结论尚未测量。

首轮 transport_v1 的四个 jobs 在配置阶段退出：主重建的旧 observer 库没有注册后来新增的 `CalypsoSlopeTransportProbeAlg`。未处理事件，也未损坏已有 states。v2 保持 inputs 不变，改用此前 direction/stage 控制已成功的 `S/diagnostics_v1/software_{beforeGeo,afterGeo}` transport 组合，不重建 tracking；失败日志及恢复来源均保留。

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

**证据支持 segment 的周期首先来自 measured strip coordinate；尚不能把最终肩峰的全部幅度归因于该机制。** 生产及独立 gun 的结果在 `S/centroid_analysis_v1/{rdo,gun}/`，包含同 hits 配对、逐事件统计、PNG/PDF 和 source manifests。

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

批量全部经 Condor；每次提交重新比较 schedd，选择兼容的最低 running+idle bigbird20，排除已有 EOS 兼容失败的 bigbird21。源任务、schedd 广告/队列、资源和 cluster ID 保存于 `S` 下对应 `*_worker/submission.json`；home 只保存小型 submit。最新归约为 production centroid 1573899、独立 gun 分析 1573903、图 1573910、近零统计 1573911；native writer smoke 1573908/1573909 全部 exit0。

已有 24 个本地单测通过；新增 native writer 归约已用两版真实 20-event smoke 验证，保存字段逐值闭合。`S/patch_validation_v1/validation.json` 记录 0006–0009 可应用且与实际源码对应内容相同。当前图表总验收与 A5 最终来源表尚未完成。

## 阶段性结论与剩余工作

| 结构 | 已测证据 | 当前判断与未完成验收 |
| --- | --- | --- |
| legacy angle 大范围失真、历史 pz 窗口骤减 | 实际 Ntuple 与物理 frame/原宏闭合 | writer/表达式问题已定位；主物理图使用同面正确 slope |
| segment 周期梳状结构 | full RDO 与独立 G4 的 locX 单因素均大幅抑制 | measured coordinate 来源有独立验证；对最终 shoulders 的量化映射未闭合 |
| VetoNu 的 mrad 宽度 | 真实 geometry≈250 X0，full G4 的 100/250 宽度比与 1/p 相符 | 靶散射是宽度主要来源；不是本征 tracker 分辨率 |
| 逐 state shared smoothing 异常 | 恒等式在隔离后恢复，但 final 逐值不变 | 解释逐 state 污染，不解释该批 final 结构 |
| beforeGeo 近零尖峰 | truth-p 选样仍可见；无 exact-zero 堆积 | 需 native seed/CKF/refit 同面定位、固定候选最小实验；未完成 |
| afterGeo −0.71/+0.59 mrad shoulders | phase response 分别相容/仍有 72.56 条局部差值 | 不宣称全部已解释；需计入未覆盖 reco、same-hit 控制和模型系统项 |

下一步是完成 native writer 的 10k/版本诊断与真实 Calypso 传播，先定位 beforeGeo 近零结构首个 stage；在主样本固定窗口量化其 p-tail、charge、truth-slope 组成。随后针对来源做最小控制，在独立输入上检验；扩充窄肩峰的 same-hit 有效统计，完成条件响应与未覆盖部分的余量表。最后完成统一轴图集、总 PDF/index、counts/flows/invalid 验证和 A1–A5 逐项审计。本报告不将尚未解释的结构转移到后续阶段以宣称 S05 完成。
