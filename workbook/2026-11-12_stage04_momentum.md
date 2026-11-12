# S04：动量低估的逐运算定位与材料坐标一致性验证

日期：2026-11-12（Europe/Zurich）。阶段结果目录：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum`，下文以 `$R` 表示。

## 阶段目标

按 roadmap M1–M6 区分 seed 的初始低估、真实上游能损、CKF 内部更新、转换/refit 与候选选择。验收采用固定输入/候选的逐运算证据、总体发生率、独立配对验证，以及效率、charge、angle、failure 和资源回归。

## 已确认的主要结果

当前旋转几何与未旋转 current22 材料查表坐标不一致，是本配置中后续拟合大负尾的重要程序原因。严格旋转同一份 map，保留全部 geometry IDs 和 bin data，在五类独立样本上均显著降低 production-referenced tail；A′ 110051 从 4272/12509（34.151%）降至 1401/13288（10.543%），Good95 效率提高 4.018 个百分点，bootstrap95 为 [3.544,4.529] 个百分点。

生产点 residual 包含真实上游能损。历史 RDO 只能从 continuation vertices 给局部 p 的上下界；本阶段另生成与 seed plane 同面的 G4 样本，执行确定性 q/p、完整 mean 和纯 truth hits 非组合 KF 实验。

另确认 shared ancestor 的 smoothed slots 可被其他分支覆写，污染逐层诊断。独立重算使 smoother 恒等式恢复闭合。全量100000次事件处理的136540条final tracks中，p/charge和Good95判定均不变；部分rmax元数据会变。该诊断覆写没有解释本批最终动量负尾。

## 工作内容与冻结定义

- 发现样本沿用 S03 五类各 10000 事件，共 50000 事件、70000 个 signal truth。事件身份为 DSID/file_index/source_entry；A′ 两个子粒子共同按源事件 bootstrap，不使用会复用的 run/lumi/event 去重。

- 主选择：truth fiducial + 下游station1/2/3的SiHit计数和≥18（H18），SDO majority≥0.7，最终 Good95：χ²/ndof<25、至少7层、ndof+5≥12、rmax<95 mm。匹配不使用动量或电荷。分别保留所有 truth、未重建、未匹配、重复、失败和 fallback。

- 三站 forward 主实验冻结 FASERNU-04、OFLCOND-FASER-05、maskedLayers=0..5；所有候选间比较通过 final ancestry 的 seed/candidate ID。参考动量默认是 production truth；同面真值和上下界另标。所有 p ratio 均使用有效正的 Cartesian norm。

- `current`为冻结的未旋转current22，`strict`为同一map的严格刚体旋转current22。历史AFS rotated12与strict的共同材料数据相同，缺少的10项均为空材料，部分binUtility平移沿法向相差2 mm；本阶段保留22项以固定原map全部IDs。逐项差异见[材料图比较](2026-11-12_stage04_material_maps.md)。本阶段改善不是AFS旋转图→strict的测量。

- 主材料配对和独立验证均使用 S03 的冻结重建库，未同时启用 chart/refit plane 等修正。新 observer/control 的普通开关默认保持原行为。

- 独立样本每类另取 10000 事件：单粒子文件 entries10000–19999；110033 使用预留后5个文件，每个2000事件；110051 使用预留后10个文件，每个1000事件。具体文件/range 见 `$R/independent_inputs.json`，不与 S03 重叠。

- 全部批量 reconstruction/simulation/逐层统计使用 HTCondor。每次提交前重新记录 schedd 负载；bigbird21 实际拒绝当前 EOS executable/log 布局，后续排除该不兼容节点，再在兼容普通 bigbird 中按 running+idle 最小选择。实际选中 bigbird20，依据保存在各 `*_worker/schedd_load.txt`、`schedd_queue.txt` 和 `submission.json`。

## 详细数据：材料配对

下表的 tail 为 `p_reco/p_production<0.5`。比例分母是该配置自己的 Good95 匹配粒子，效率分母为共同 fiducial+H18 truth。

### 发现样本

| DSID | current Good95 / tail / wrong | strict Good95 / tail / wrong | tail 比例 current→strict / % | Δtail / pp（event bootstrap95） | ΔGood95 效率 / pp（event bootstrap95） |
| --- | --- | --- | --- | --- | --- |
| 100044 | 6965 / 791 / 369 | 6963 / 39 / 8 | 11.357→0.560 | -10.797 [-11.479,-10.023] | -0.024 [-0.145,0.085] |
| 100116 | 4868 / 75 / 5 | 4884 / 16 / 0 | 1.541→0.328 | -1.213 [-1.550,-0.904] | 0.276 [0.069,0.501] |
| 100117 | 5047 / 110 / 34 | 5052 / 32 / 0 | 2.180→0.633 | -1.546 [-1.932,-1.190] | 0.084 [-0.151,0.318] |
| 110033 | 12325 / 1510 / 265 | 12778 / 879 / 37 | 12.252→6.879 | -5.373 [-5.849,-4.860] | 2.707 [2.336,3.048] |
| 110051 | 12460 / 4182 / 1695 | 13242 / 1405 / 305 | 33.563→10.610 | -22.953 [-23.757,-22.140] | 4.043 [3.535,4.574] |

| DSID | current p ratio 16/50/84% | strict p ratio 16/50/84% |
| --- | --- | --- |
| 100044 | 0.59278 / 0.92623 / 1.18978 | 0.89556 / 0.95727 / 1.02208 |
| 100116 | 0.76789 / 0.91782 / 1.10296 | 0.87659 / 0.93391 / 0.97328 |
| 100117 | 0.75103 / 0.91964 / 1.10943 | 0.87589 / 0.93662 / 0.97881 |
| 110033 | 0.58345 / 0.93260 / 1.14641 | 0.75807 / 0.97926 / 1.01537 |
| 110051 | 0.22493 / 0.78105 / 1.24140 | 0.66480 / 0.98201 / 1.22487 |



### 独立验证

| DSID | current Good95 / tail / wrong | strict Good95 / tail / wrong | tail 比例 current→strict / % | Δtail / pp（event bootstrap95） | ΔGood95 效率 / pp（event bootstrap95） |
| --- | --- | --- | --- | --- | --- |
| 100044 | 6959 / 774 / 371 | 6960 / 40 / 5 | 11.122→0.575 | -10.548 [-11.308,-9.868] | 0.012 [-0.096,0.121] |
| 100116 | 4902 / 72 / 21 | 4918 / 23 / 1 | 1.469→0.468 | -1.001 [-1.348,-0.711] | 0.274 [0.103,0.462] |
| 100117 | 5038 / 120 / 32 | 5041 / 20 / 1 | 2.382→0.397 | -1.985 [-2.408,-1.581] | 0.050 [-0.133,0.250] |
| 110033 | 12254 / 1510 / 297 | 12702 / 903 / 31 | 12.323→7.109 | -5.213 [-5.671,-4.698] | 2.685 [2.321,3.054] |
| 110051 | 12509 / 4272 / 1815 | 13288 / 1401 / 311 | 34.151→10.543 | -23.608 [-24.447,-22.856] | 4.018 [3.544,4.529] |

| DSID | current p ratio 16/50/84% | strict p ratio 16/50/84% |
| --- | --- | --- |
| 100044 | 0.60063 / 0.92176 / 1.19379 | 0.88875 / 0.95593 / 1.02055 |
| 100116 | 0.76850 / 0.91677 / 1.11369 | 0.87087 / 0.93449 / 0.97241 |
| 100117 | 0.75763 / 0.91937 / 1.12283 | 0.87432 / 0.93822 / 0.98058 |
| 110033 | 0.57832 / 0.92786 / 1.13888 | 0.74606 / 0.97669 / 1.01473 |
| 110051 | 0.22100 / 0.76986 / 1.24215 | 0.66110 / 0.98026 / 1.20796 |



![Paired material validation](../../geoRotation_test_results/stage04_momentum/validation/paired_material_validation.png)

五类独立样本的 tail 改善区间均在 0 以下；效率保持或提高。材料修正还减少 wrong charge，μ± 独立样本各剩1条 wrong。正负电荷残余差异须在 S06 做共同谱/位置/方向支持域比较；角度多峰的显著性与机制由 S05 检验。

材料图旋转前后的seed、九层predicted/filtered/smoothed、good-seed与final分布，以及用户提供绘图源码的逐项核对见[S04补图报告](2026-11-12_stage04_reference_figures.md)。该报告提供PNG、逐图PDF和完整汇报图集，并区分原脚本与S04的选择/truth定义。

## 逐层材料的固定轨迹证据

100044/file0/entry21，production p=1708.985 GeV，seed p=1436.63 GeV。current 下 found=258.788 GeV、final=262.020 GeV；strict 下 final=1301.198 GeV。相同18个 measurement IDs，初始 seed 不变。selected 路径的完整 slab 查询和由 55.05665 X0 降至 0.0486853 X0；path correction 均约1。

| z_F / mm | current slab X/X0 | current thickness / mm | strict slab X/X0 | strict thickness / mm |
| --- | --- | --- | --- | --- |
| 38.858 | 0.0074 | 0.830 | 0.006086 | 0.570 |
| 70.358 | 0.0584 | 5.750 | 0.006086 | 0.570 |
| 1197.358 | 1.9461 | 45.429 | 0.006086 | 0.570 |
| 1228.858 | 1.9213 | 44.759 | 0.006086 | 0.570 |
| 1260.358 | 50.9495 | 1040.008 | 0.006086 | 0.570 |
| 2387.358 | 0.0577 | 5.194 | 0.006086 | 0.570 |
| 2418.858 | 0.0600 | 5.408 | 0.006086 | 0.570 |
| 2450.358 | 0.0563 | 5.070 | 0.006086 | 0.570 |

上述量是在每个 predicted 位置查询完整 slab；ACTS 可在传播/更新前后分裂施加材料，不能将该和解释成已测得的逐次能损积分。巨大厚度差来自查表，不是本例的斜入射路径修正。关闭 scattering 与关闭 energy loss 的单因素结果用于区分机制。

## 历史 RDO 的局部真值界

使用同 PDG、`barcode % 1000000 == primary barcode` 的 continuation chain，按 vertex.z 排序，要求 forward 且 p 单调。参考面前最近 p 为上界、后最近 p 为下界，无后续点时下界为0。不做连续能损插值。分类：保证为 local tail、排除 local tail、区间尚不能区分。

| DSID | production tail | guaranteed local tail | excluded local tail | interval unresolved |
| --- | --- | --- | --- | --- |
| 100044 | 791 | 452 | 75 | 264 |
| 100116 | 75 | 17 | 28 | 30 |
| 100117 | 110 | 37 | 36 | 37 |
| 110033 | 1510 | 730 | 344 | 436 |
| 110051 | 4182 | 3152 | 391 | 639 |

这是 current、主选择的 final tail 分类。100044/entry670 在 production 为336.625 GeV，但 z_F=−2021 mm 已只有90.3835 GeV，最终121.85 GeV，不能把 production tail 全归为重建错误。110051/entry1/barcode2 从508.684 GeV 在 z≈15.3 mm 发生大幅能损，后续 continuation 为91.655 GeV，最终98.683 GeV。精确 sensor p 在这些历史 RDO 中未测量。

## 全量逐运算审计与 seed 分类

另对同一发现样本运行 current/strict 各 50000 事件，保存每个选中路径的 predicted/filtered/smoothed、signed κ、innovation、Kκ、材料 slab、candidate ancestry 与 refit；完整矩阵仅留 worker scratch，结果为 compact traces。`operations_tasks.json` 为100个任务，全部成功。`operations_analysis/{current,strict}` 保存粒子表、分类、异常位置与英文逐层图。

按最终 Good95 track 的祖先路径寻找最早已经低于0.5的 stage；stage 内再定位正常进入、低比值离开的操作。暂时偏低后恢复的 `first_new_tail_operation` 与最终负尾的 `tail_entry_operation` 分开，不能将任意一次阈值穿越直接当成最终根因。

| DSID | current tail：seed / found / refit | strict tail：seed / found / refit | strict local：guaranteed / excluded / unresolved |
| --- | --- | --- | --- |
| 100044 | 33 / 726 / 32 | 24 / 14 / 1 | 7 / 29 / 3 |
| 100116 | 9 / 59 / 7 | 10 / 6 / 0 | 1 / 15 / 0 |
| 100117 | 21 / 84 / 5 | 20 / 12 / 0 | 0 / 31 / 1 |
| 110033 | 841 / 621 / 48 | 614 / 262 / 3 | 212 / 340 / 327 |
| 110051 | 2065 / 2014 / 103 | 967 / 435 / 3 | 519 / 426 / 460 |

逐层 κ update 恒等式 `κfiltered−κpredicted = Kκ r` 在两份全量 Good95 队列中的最大归一化误差为 7.87×10⁻¹⁰（以 |κtruth| 归一）；独立重算后的 smoother 最大为 5.25×10⁻⁸，阈值10⁻⁶以上不闭合计数均为0。ACTS→Trk 与 Trk→refit-input 的动量模长最大相对改变量均为6.67×10⁻¹⁶。接口 covariance 问题仍见S02；这些数值只验证本阶段的动量 mean，不能代替角度/covariance 的 S05/S07 验收。

原 shared ancestor smoother slot 不闭合的计数，current 五类依次为144/104/97/3785/7392，strict 为137/92/86/3758/7651。独立重算用于得到同一 candidate 内一致的诊断 states。其对全部输出字段的影响单独保存于 `observer_integrity/`，最终物理统计以未启用该开关的冻结名义重建为准。

`observer_integrity` 对每个stage逐字段核对全部100000次事件处理。current/strict分别有70170/66370条final，动量、电荷、位置、bound、κ variance、hits、匹配、χ²/ndof、层数、refit状态和Good95判定均一致；唯一差异是 `rmax_mm`：converted分别8832/6625条、final分别37/36条。例100044/entry27的converted rmax从83.6141变为76.5917 mm。独立smoothing会改变用于计算rmax的state位置，refit/fallback后少量rmax仍不同，但本批没有跨越Good95半径阈值。观察开关的适用范围由这份完整差异表限定。

下表是发现样本在Good95之前的完整输出计数，每格为current→strict，每类均10000源事件。无final按事件计数，其他列按track或finder错误记录计数；finder错误不等于整事件失败，refit fallback表示保留CKF结果。`unmatched_to_signal`只表示没有匹配到本研究的signal truth，包含其他真实粒子，未把它直接解释成物理fake率。

| DSID | final tracks | 无final事件 | unmatched to signal | signal重复tracks | finder错误记录 | refit fallback tracks |
| --- | --- | --- | --- | --- | --- | --- |
| 100044 | 9687→9684 | 324→324 | 72→67 | 0→0 | 4→4 | 723→727 |
| 100116 | 9104→9103 | 897→897 | 64→64 | 0→0 | 0→1 | 1374→1389 |
| 100117 | 9068→9067 | 934→934 | 70→70 | 0→0 | 1→1 | 1265→1288 |
| 110033 | 19556→17588 | 360→333 | 4947→2574 | 56→20 | 617→416 | 487→529 |
| 110051 | 22755→20928 | 438→371 | 8866→6192 | 73→42 | 1806→1298 | 315→318 |


### Seed：数值圆拟合与实际输入分开

实际生产 seed 使用 segment center 和沿 segment direction 的±30 mm fake points，在 FASER 的 y–z 平面做 Taubin circle；诊断 spacepoint 列表并非其实际输入，该列表还有 static map 跨事件/多候选覆写。清空诊断 maps 不改变实际 fake-point 圆拟合输入。

理想圆 p=10/50/100/250/1000 GeV、两电荷的 radius→p closure 最大相对误差1.31×10⁻⁹；fake tangent points 的偏差约4.0–4.83×10⁻⁷。57个代表 seed 的同输入 Taubin 重算相对差≤1.301×10⁻¹²，独立 geometric circle 最大半径改变量0.0314%。这些测量支持继续检查 segment 来源、方向和场近似，未发现足以解释巨大负尾的圆方程数值误差。

以下分母是最终 Good95 的祖先 seed，包含后来被 CKF 恢复的 seed tail；因此 seed tail 数可大于最终 tail。互斥标签优先级为 segment majority不属于该truth、局部界排除低估、direction替换敏感、guaranteed local、局部界待定。direction 替换只保留原中心、使用 station 的平均 truth 方向，是敏感性诊断，未当成精确同点修复。

| DSID | current seed tail / ancestry seeds | strict seed tail / ancestry seeds | strict seed tail：other segment / local excluded / direction sensitive / local guaranteed / unresolved |
| --- | --- | --- | --- |
| 100044 | 45 / 6965 | 43 / 6963 | 3 / 20 / 16 / 3 / 1 |
| 100116 | 9 / 4868 | 10 / 4884 | 0 / 10 / 0 / 0 / 0 |
| 100117 | 24 / 5047 | 23 / 5052 | 0 / 20 / 3 / 0 / 0 |
| 110033 | 1889 / 12325 | 1462 / 12778 | 928 / 146 / 43 / 173 / 172 |
| 110051 | 3791 / 12460 | 3020 / 13242 | 2346 / 153 / 146 / 161 / 214 |

110051 strict 的3020条祖先 seed tail 中，2346条至少一个 segment 的 majority barcode 属于别的粒子；另146条对 segment direction 替换敏感。它们支持 cross-track segment 组合和方向估计作为后续重点，但不能仅凭标签给出2346条“已修复”结论。固定B=0.55 T与真实场积分、完整 segment covariance 和 charge 分组继续在S07检验。无有效有限 seed p 的计数 current为0/2/0/2/4，strict为0/2/0/4/7，未把无穷半径写成物理动量。

候选 pruning 也有独立贡献。定义 accurate-found 为 production p ratio∈[0.8,1.2]、同电荷且至少12 measurements；current 中“有 accurate-found 但没有 accurate-selected”的五类计数为51/47/0/1205/2133，strict 为12/5/0/459/1717，完整分母见 `candidate_counts`。这是选择损失诊断，尚未把 production 相对 accurate 等同局部 truth 正确。

## 同 seed、同 hits 的总体配对与 core

只保留两配置均为 Good95、相同 selected measurement ID 序列和 seed 的粒子，避免把改善全部归于换候选。下表的材料量仍为 predicted 位置完整 slab 查询和。

| DSID | common Good95 | same seed+selected hits | tail recovered / introduced | material X/X0 median current→strict |
| --- | --- | --- | --- | --- |
| 100044 | 6953 | 6856 | 743 / 5 | 56.15395→0.09007 |
| 100116 | 4855 | 4721 | 59 / 3 | 56.06396→0.08441 |
| 100117 | 5024 | 4876 | 82 / 9 | 56.06075→0.08639 |
| 110033 | 12081 | 10361 | 446 / 88 | 56.05510→0.08180 |
| 110051 | 11895 | 9260 | 1799 / 116 | 56.17820→0.07690 |

同一批 hits 并不保证拟合在完全相同的空间位置查表；这正是材料/导航模型敏感性的组成部分。多数配对路径的材料 slab 从约56 X0 降至约0.08–0.09 X0。不同轨迹仍可查到局部较厚材料，不能把严格旋转解释成所有轨迹固定为零材料。单独的 gain 最大值分位数没有统一单调下降，关键是 `Kκ×innovation` 在具体错误更新上的变化。

![Same selected paths](../../geoRotation_test_results/stage04_momentum/validation_complete/same_path_material.png)

独立样本另预定义 core 窗口0.8≤p ratio≤1.2，报告窗口内 median 与(84%−16%)/2；它不是完整分布分辨率。event bootstrap 512 replicas 的95%区间及0.8/0.5/0.2三档 tail 原始计数均见 `validation_complete/independent_core.json`。

| DSID | core median current→strict | core half68 width current→strict | strict median bootstrap95 | strict width bootstrap95 |
| --- | --- | --- | --- | --- |
| 100044 | 0.97032→0.95758 | 0.11820→0.05348 | [0.95614,0.95911] | [0.05212,0.05525] |
| 100116 | 0.94054→0.93757 | 0.11124→0.04410 | [0.93610,0.93900] | [0.04236,0.04563] |
| 100117 | 0.94465→0.94076 | 0.11327→0.04660 | [0.93927,0.94237] | [0.04502,0.04848] |
| 110033 | 0.98768→0.98875 | 0.11351→0.04999 | [0.98806,0.98970] | [0.04844,0.05182] |
| 110051 | 0.98737→0.99273 | 0.12074→0.09535 | [0.99026,0.99481] | [0.09332,0.09732] |

## 同面 truth seed 与纯 hits：T0–T4

PDG=±13/±11、p=100/1000 GeV，每格300个tracker+Dipole-only Geant4事件，共2400个母事件；每个母事件只simulate/digitize一次，再运行12种受控模式，共28800次模式事件。gun=(20,−20,11.415) mm、θ=0.002 rad、φ=0，seed plane与gun同面。T0/T2/T3/T4默认使用strict材料；同母样本的current变体另标。

T2只替换signed q/p；T3替换完整6维局部mean（含position、direction、q/p、time=0 ns），保留原C；T4在T3基础上只选唯一贡献barcode属于该primary的纯clusters，实际调用非组合ACTS KF。共2406个候选seed记录，其中2400个满足primary与同面注入条件，另6个不注入；同面距离最大3.815×10⁻⁸mm。12000次实际mean/C检查均一致。没有把带smearing/固定正电荷的TruthBasedInitialParameterTool当作truth注入。

| PDG / p [GeV] | T0 strict Good95 / tail / wrong | T2 strict | T3 strict | T4 strict | T4 current |
| --- | --- | --- | --- | --- | --- |
| -13 / 100 | 299 / 0 / 0 | 299 / 0 / 0 | 297 / 0 / 0 | 300 / 0 / 0 | 300 / 1 / 0 |
| -13 / 1000 | 296 / 1 / 0 | 296 / 1 / 0 | 295 / 1 / 0 | 297 / 0 / 0 | 297 / 246 / 1 |
| 13 / 100 | 298 / 0 / 0 | 298 / 0 / 0 | 298 / 0 / 0 | 300 / 0 / 0 | 300 / 0 / 1 |
| 13 / 1000 | 297 / 0 / 0 | 297 / 0 / 0 | 296 / 0 / 0 | 298 / 0 / 0 | 298 / 15 / 51 |
| -11 / 100 | 296 / 6 / 0 | 296 / 6 / 0 | 296 / 6 / 0 | 295 / 6 / 0 | 295 / 10 / 0 |
| -11 / 1000 | 299 / 8 / 0 | 299 / 8 / 0 | 299 / 8 / 0 | 300 / 9 / 0 | 300 / 255 / 3 |
| 11 / 100 | 296 / 3 / 0 | 296 / 3 / 0 | 295 / 3 / 0 | 295 / 4 / 0 | 295 / 3 / 0 |
| 11 / 1000 | 299 / 10 / 0 | 299 / 10 / 0 | 297 / 10 / 0 | 298 / 10 / 0 | 298 / 16 / 58 |

T1在finding前按 `abs(pseed/ptruth−1)<0.1`、同电荷、同A-frame `abs(θseed−θtruth)<10⁻⁴ rad` 定义合格seed；母truth须至少一个合格seed，final须来自这些seed。失败、被prune或改由其他seed选中仍留分母，另给any-seed结果；以下真实母样本均无该类失败/换seed，但回归测试覆盖它们。

| PDG / p [GeV] | 合格seed母truth | T1 current Good95 / tail / wrong | T1 strict Good95 / tail / wrong |
| --- | --- | --- | --- |
| -13 / 100 | 0 | 无合格seed | 无合格seed |
| -13 / 1000 | 59 | 59 / 58 / 0 | 59 / 0 / 0 |
| 13 / 100 | 0 | 无合格seed | 无合格seed |
| 13 / 1000 | 0 | 无合格seed | 无合格seed |
| -11 / 100 | 4 | 4 / 0 / 0 | 4 / 0 / 0 |
| -11 / 1000 | 70 | 70 / 64 / 1 | 70 / 0 / 0 |
| 11 / 100 | 1 | 1 / 0 / 0 | 1 / 0 / 0 |
| 11 / 1000 | 2 | 2 / 0 / 0 | 2 / 0 / 0 |

1 TeV μ+ 的good-seed母样本有59个，current下58/59为负尾，strict下0/59；1 TeV e+对应64/70→0/70。它们直接复现good seed后仍低估的现象，并表明该现象受材料坐标控制。其他格T1的样本量为0–4，属于强烈条件选择，空母样本没有效率点；不据此比较整体重建性能。

| 1 TeV PDG | T4 current Good95 / tail / wrong | strict | current no scattering | current no energy loss | Δtail strict−current / pp（paired bootstrap95） |
| --- | --- | --- | --- | --- | --- |
| -13 | 297 / 246 / 1 | 297 / 0 / 0 | 297 / 0 / 0 | 297 / 246 / 1 | -82.828 [-87.162,-78.783] |
| 13 | 298 / 15 / 51 | 298 / 0 / 0 | 298 / 0 / 0 | 298 / 15 / 51 | -5.034 [-7.744,-2.685] |
| -11 | 300 / 255 / 3 | 300 / 9 / 0 | 300 / 9 / 0 | 300 / 255 / 3 | -82.000 [-86.333,-78.000] |
| 11 | 298 / 16 / 58 | 298 / 10 / 0 | 298 / 10 / 0 | 298 / 16 / 58 | -2.013 [-4.013,0.000] |

这里材料修正没有减少各格T4的Good95计数。truth完整mean、纯hits仍保留current的大负尾，strict/no-scattering使其大幅消退，no-energy-loss保留大尾；结合逐层slab和update恒等式，证据指向材料散射协方差→transport covariance→Kκ×innovation 的错误放大，而非单位转换或仅seed p错误。每格约300事件的零观测仍有约1.3%的95%上界；总体发生率使用上面的独立历史样本配对。

![Truth seed and fixed hits](../../geoRotation_test_results/stage04_momentum/simulation_analysis_verified/pdg-13_p1000_injections.png)

严格材料下T4 μ±四格合计1195条Good95，无ratio<0.5、无wrong sign；e±四格仍有29条production tail，其中23条由sensor界保证为local tail、6条区间待定。更换electron hypothesis、chart/refit接口、seed C×0.1没有系统消除此残余。该结论定位到固定真hits的KF/材料模型，G4真实能损与拟合模型需进一步区分；没有将它们全部归为production-reference定义效应。

新sensor truth覆盖2400个Athena事件、81153个step，和原模拟除事件标识外逐值相同；59个measurement-state关联因surface/primary覆盖不能找到0.6mm内sensor界，单列为missing，不用它们构造精确truth。`sensor_layers.json`保存逐层predicted/filtered/smoothed对sensor p上界的分位数，final粒子另保留上下界分类。

## 定向单因素控制

从 S03 冻结的30个事件、32条无匹配歧义代表 track 出发，baseline 重放的 final p/charge 与 S03 完全一致。它们按 normal、seed/found/refit tail 与 wrong-charge 刻意抽取，以下结果仅用于机制判别，不是总体发生率。

| 受控配置 | matched / Good95 | tail <0.5 / wrong |
| --- | --- | --- |
| baseline | 32 / 32 | 21 / 14 |
| no_scatt | 31 / 30 | 6 / 1 |
| no_eloss | 32 / 32 | 21 / 14 |
| no_material | 31 / 30 | 6 / 1 |
| bypass | 32 / 32 | 16 / 10 |
| chart | 32 / 32 | 21 / 14 |
| plane | 32 / 32 | 21 / 14 |
| isolated | 32 / 32 | 21 / 14 |
| electron | 31 / 31 | 21 / 12 |
| refit_all1 | 32 / 32 | 16 / 11 |
| refit_diag1 | 32 / 32 | 19 / 11 |
| seedC01 | 31 / 30 | 18 / 11 |
| step10 | 31 / 31 | 20 / 13 |

`no_scatt` 同时关闭 CKF 与 refit 的实际 scattering；`no_eloss` 同时关闭实际 energy loss；没有用 `resolveMaterial=False` 代替实际物理开关。关闭 scattering 大幅改变选中路径的负尾/电荷，关闭 energy loss 在这些案例不改变结果，符合错误材料厚度主要通过过程噪声和 gain 影响 curvature 的机制。refit bypass、covariance inflation 和 step=10 mm 也会改变部分案例，不能据此替代材料修复。`fixed_hits` 案例使用普通 seed，并非正式 T4，正式结论使用前节的确定性 truth mean+纯 hits。

此前运行过`configs/acts_momentum_toy.cxx`的30000次独立ACTS toy fits，它使用均匀场和自建strip planes，不属于Calypso实际重建。根据用户本次明确要求，该项退出本报告的物理结论和T5验收；旧代码和运行记录仅保留历史来源。T5在Calypso实际几何/组件内的完整统计closure未测量，列入S07。主材料配对、逐层分布和上表单因素结果来自实际Calypso作业，来源详见[专项核查报告](2026-11-12_stage04_calypso_provenance.md)。

## 固定真 hits 的电子剩余尾：局部能损与运算

对T4 strict全部Good95再用同sensor的pre/post p界追踪。23条guaranteed local tail中，11条首次可由同sensor界证实的穿越发生在measurement update，12条发生在backward smoothing；它们都出现在第一次固定hits KF中。这里的“首次”限定于有sensor界的measurement update/smoothing，不把两个不同面的状态直接作同点比较。

| PDG / p [GeV] | local tail / Good95 | first measurement update / smoothing | p_last,upper / p_first,lower 的 min / median / max |
| --- | --- | --- | --- |
| pdg-11_p100 | 4 / 295 | 0 / 4 | 0.0864 / 0.2659 / 0.2957 |
| pdg-11_p1000 | 7 / 300 | 5 / 2 | 0.1038 / 0.2146 / 0.5709 |
| pdg11_p100 | 3 / 295 | 0 / 3 | 0.2209 / 0.3367 / 0.4030 |
| pdg11_p1000 | 9 / 298 | 6 / 3 | 0.0693 / 0.2107 / 0.3983 |

23条电子局部负尾的末端/首端真实p比均≤0.571，显示显著沿程能损；四格μ±则末端/首端比最低仍约0.999。量化对照支持在辐射能损后的轨迹模型中研究这些残余，单凭这份相关性尚未确定全部过程来源。

- e+、100 GeV、source_entry124：首sensor p∈[99.999904,99.999997] GeV，末sensor约24.1535 GeV；z_F=22.300 mm 的filtered p=99.999962 GeV，smoothing后24.263338 GeV，最终24.263571 GeV。后部低p的信息被带回首端，直接证实同面低估的数值运算。

- e+、1000 GeV、source_entry36：首sensor约1000 GeV，末sensor118.434 GeV；z_F=1211.415 mm局部truth约245.0086 GeV，update前491.2643 GeV、后86.6665 GeV，r=−0.946240 mm、Kκ=−0.0100428 GeV⁻¹mm⁻¹；最终首sensor p=156.6963 GeV。另在首sensor smoothing中可见约1000→156.641 GeV的变化。

相应完整事件、前后参数和truth界见 `g4_residuals/summary.json`。这既区分了“production点到首sensor已经掉能”和“在tracking区掉能后，拟合将错误曲率传播回首sensor”，也解释了为什么强制正确seed仍不足以保证最终第一面p正确。

后续工作假设是Gaussian KF的能损/过程噪声模型不能充分描述这些stochastic electron trajectories。检索了W. Adam等《Reconstruction of electrons with the Gaussian-sum filter in the CMS tracker at LHC》，[arXiv:physics/0306087](https://arxiv.org/abs/physics/0306087)、[DOI:10.1088/0954-3899/31/9/N01](https://doi.org/10.1088/0954-3899/31/9/N01)。论文用Gaussian mixture描述非Gaussian bremsstrahlung能损，适合这里的固定真hits电子残余；Calypso的具体场/材料与测量模型仍需单独验证。S07先做关闭G4 bremsstrahlung的匹配生成对照、记录radiation位置与分段p，再比较piecewise q/p或GSF；验收看首sensor local residual/pull及效率，不把electron mass开关等同于GSF。网页metadata/摘要保存于gsf_reference.json。

## Backward、4-station 与实际 merge

五类样本各取500源事件，对 current/strict 分别重建，共5000次源事件处理。每次实际运行 forward3/backward3/forward4/backward4，并分别执行真正的 MergeTracksAlg；没有用离线任选一条track代替merge。所有任务成功，5000事件、15712条merged tracks 的动量/位置/charge/hits与其输入来源逐项一致，来源歧义与不一致均为0。

下表展示110051（每配置500源事件、1000个signal truth）的 fiducial+H18、Good95且最终hits包含downstream station1/2/3 的结果。全五类和配对bootstrap区间见 `modes_analysis_complete/summary.json`，本表不承担独立验证的最终精度。

| finding / merge mode | current Good95 / tail / wrong | strict Good95 / tail / wrong |
| --- | --- | --- |
| forward3 | 625 / 220 / 104 | 668 / 71 / 18 |
| backward3 | 545 / 176 / 90 | 604 / 57 / 21 |
| merge3 | 659 / 235 / 118 | 711 / 78 / 24 |
| forward4 | 625 / 220 / 104 | 668 / 71 / 18 |
| backward4 | 543 / 178 / 89 | 603 / 57 / 21 |
| merge4 | 656 / 233 / 116 | 711 / 78 / 24 |

4-station finding 确实启用IFT，并非mask配置无效。以100116 μ+ strict forward4为例，345个seed、382个found、334个converted tracks带全四站hits；final的全四站track为0，347条带station1/2/3。原因是既有 `TrackFitterFunction.cxx:4–11` 在 `cluster_z<-10000` 的refit模式下将F frame `globalPosition().z()<-100 mm` 的IFT measurements判作outlier；`KalmanFitterTool::fitAudited` 显式设置了该模式。refit fallback仍可能保留部分IFT，但本例没有全四站final。

因此“4-station finding”的结果不能解释成四站共同进入最终refit。A′ 在IFT后衰变时本就没有IFT seed，forward3/4可相同。仅Good95的≥7层proxy还可能接纳IFT+两个downstream stations；本轮额外要求最终1/2/3站后，current forward4 μ+从244变为243，其他所列模式/样本计数不变。station位解码在真实truth_hit的Identifier与已存station之间逐项验证。

## 转换、refit 与 Ntuple 实际输出

转换/refit-input 的 p mean 单位与模长已由全量逐运算审计闭合。另用真实 NtupleDumperAlg 对100044 entries21–22与110051 entries1–2分别运行current/strict，检查 ROOT Tree 的 Track_p0/px0/py0/pz0，而非只审读源码。

| material / DSID | events | Ntuple tracks / Good95 final tracks | max norm relative error | max final p relative error |
| --- | --- | --- | --- | --- |
| current / 100044 | 2 | 2 / 2 | 0 | 1.11e-16 |
| current / 110051 | 2 | 5 / 5 | 0 | 2.22e-16 |
| strict / 100044 | 2 | 2 / 2 | 0 | 0 |
| strict / 110051 | 2 | 4 / 4 | 0 | 1.11e-16 |

共13条tracks的p分量norm完全一致；audit的GeV与Ntuple的MeV换算后最大相对差2.23×10⁻¹⁶，所有检查样本中的Good95 final tracks均在Ntuple找到。此入口创建空 waveform containers，以满足tracking-only RDO未运行waveform reco的接口；calorimeter能量和时间未测量。ROOT只在scratch中验证，不保存中间xAOD。

历史RDO regeneration界的源码依据是 `FaserTruthSvc.cxx:262–276` 对同PDG延续粒子将barcode增加 `HepMC::SIM_REGENERATION_INCREMENT`；Athena24.0.41 `TruthUtils/MagicNumbers.h:30` 将其定义为1000000，相关same-generator比较也使用此模数。局部界同时要求forward与p单调，不能仅按模数把任意secondary并入。

## 剩余类别与下一项区分实验

严格旋转材料后，不将所有 production tail 当成残余程序错误。下表只统计 guaranteed local tail（最终p低于局部p下界的一半），分母仍为该配置全部 Good95；区间按这些选中粒子的源事件聚类 bootstrap。零观测的 bootstrap 区间退化为0，不据此声称真实概率严格为0。

| DSID | guaranteed local tail / Good95 | 比例 / %（bootstrap95） | production tail 中局部界仍待定数 |
| --- | --- | --- | --- |
| 100044 | 7 / 6963 | 0.101 [0.029,0.173] | 3 |
| 100116 | 1 / 4884 | 0.020 [0.000,0.061] | 0 |
| 100117 | 0 / 5052 | 0.000 [0.000,0.000] | 1 |
| 110033 | 212 / 12778 | 1.659 [1.419,1.861] | 327 |
| 110051 | 519 / 13242 | 3.919 [3.570,4.259] | 460 |

110051 的519条 guaranteed local tail，按首个异常 stage/seed标签拆分：232条带 other-majority segment 的 seed，175条 CKF measurement update，83条 seed 自身已保证为local tail，16条 seed 对direction替换敏感，10条 seed 当时局部界待定，3条 seed 当时被界排除但最终变为local tail。这些是观测到的起点与状态变化，不将相关标签直接提升为完整物理解释。435条 CKF-first production tail 中，175/126/134条最终分别 guaranteed/excluded/unresolved；另外3条 refit-first 中1条 excluded、2条 unresolved。完整交叉表和正常进入/异常离开的事件见 `mechanism_tables/summary.json`。

| 待查类别 | 已有区分证据 | 下一项实验与验收 |
| --- | --- | --- |
| 混合 segment 的 seed | 110051 strict 最终tail中484条 seed 标记other segment，其中232条最终guaranteed local tail | 固定当前segment组合与排序，替换为同truth segment或hit-pure组合；记录κ、方向、charge及效率，若同组合恢复而圆数值解不变，归因到segment association |
| pure/同truth seed 的方向或场近似 | 理想circle闭合，部分真实seed在只替换segment方向后越过0.5 | 使用同segment center处的局部truth方向与真实field-integral拟合，分离散射、0.55T近似和±30mm fake-point权重；在独立样本给配对差值 |
| strict下CKF更新新tail | 110051仍435条，其中175条保证为local tail；update恒等式闭合 | 固定同seed/同hits，依次控制measurement C、material scattering、branch/χ²；对首个更新比较r、S、Kκ与hit purity，区分观测错误和gain模型 |
| shared-hit pruning | strict下有accurate-found而无accurate-selected计数已量化 | 保留被淘汰候选，做同truth hit purity排序或受控branch上限扫描；同时报告效率/fakes，不通过丢弃困难tracks改善尾部 |
| 历史局部truth区间待定 | 110033/110051仍327/460条production tail无法由vertex界判定 | 在匹配拓扑的新G4样本记录同sensor p；既有历史RDO不插值成“精确truth” |
| covariance chart/refit/IFT | mean转换closure；S02已证实chart和10mm起点问题，4station refit丢IFT | S07在修正材料后分别启用chart与同点plane，显式控制IFT outlier策略；检验p/方向/pull与效率，而非仅看本轮定向案例 |

S04交付的是主要材料机制的因果闭环、完整负尾定位和剩余类别的可检验划分；残余 seed/measurement 模型优化在S07按单因素推进。S05的角度多峰显著性和S06的谱匹配电荷不对称仍未测量，本阶段的总misID改善不替代这两个问题的验收。

## 角度监测与阶段验收

独立样本的 physical slope_y(reco)−slope_y(production truth) 分位数如下，单位为10⁻³。该量含生产点到最终参考面的真实偏转，作为材料变更的监测量；S05会把truth/reco统一到同物理面再研究峰和分辨率。

| DSID | current 16 / 50 / 84% [10⁻³] | strict 16 / 50 / 84% [10⁻³] |
| --- | --- | --- |
| 100044 | -2.85820 / -0.64747 / 0.01332 | -2.86085 / -0.66579 / -0.02588 |
| 100116 | -8.55470 / -2.73396 / -0.09306 | -8.61024 / -2.74213 / -0.11736 |
| 100117 | 0.05967 / 2.08791 / 7.87447 | 0.08659 / 2.13357 / 7.79352 |
| 110033 | -1.68286 / 0.00711 / 1.71744 | -1.54518 / 0.01049 / 1.57637 |
| 110051 | -0.33399 / -0.00830 / 0.32360 | -0.22173 / 0.00318 / 0.23532 |

| S04验收项 | 实际交付与结论 |
| --- | --- |
| M1 定义/图/真值 | production与local两套定义、5类cutflow、逐层英文图、历史vertex界和2400事件sensor界；负p/invalid单列，所有有效p比使用正norm |
| M2 seed | fake-point实际输入、circle/helix闭合、独立geometric解、全量seed标签与发生数；主要混合segment类和direction敏感类已量化，精确修复留S07 |
| M3 T0–T4 | Calypso内同母事件good-seed、q/p/full-mean/纯hits KF；12000次注入检查；good seed之后负尾与材料因果控制闭合。T5完整Calypso统计closure未测量，独立toy不计入验收 |
| M4 运算 | 两份50k全量selected/refit compact traces、Kκr与smoother恒等式、材料查询与first-stage/operation标签；23条剩余electron local tail补充same-sensor运算证据 |
| M5 转换/选择/output | p mean闭合、refit/bypass/C/hypothesis/step消融、pruning损失、真实merge来源、IFT outlier源码及Ntuple branch验证 |
| M6 独立验证 | 预留50k母事件的current/strict配对，五类tail均改善且效率保持/提高；core/charge/production-slope/失败及资源表；所有剩余类有发生数、具体运算与下个实验 |

S04材料配对、T0–T4和逐运算研究已交付。已证实的主因是冻结current22材料与旋转几何的坐标不一致，主要通过散射协方差及后续gain更新制造大负尾；严格旋转同一map在独立Calypso样本得到确认。剩余seed组合、辐射能损后的KF模型和历史truth界待定类已经分别定位，尚未形成全部残余的唯一解释。独立toy退出后，T5完整Calypso统计closure属于S07待完成项。

## 资源与运行记录

主批量均申请1 CPU、3500 MB memory、8 GB worker scratch，`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`；Calypso重建和G4的单个事件串行，事件分块通过Condor并发。下表汇总成功Calypso主矩阵的GNU time实测值，wall为各任务墙钟时间之和，不是用户等待时长；不包含独立toy、被替换失败任务、额外debug、分析作业和编译。

| 矩阵 | tasks | 源事件处理或fits | sum CPU / h | sum wall / h | max RSS / MiB |
| --- | --- | --- | --- | --- | --- |
| 初轮材料+代表事件 | 80 | 50000+30 | 2.7730 | 3.5477 | 1485.7 |
| 独立材料配对 | 100 | 100000 | 8.9469 | 9.4000 | 1248.0 |
| 全量逐层重放 | 100 | 100000 | 6.9409 | 7.2959 | 1493.0 |
| G4+digitization+12模式 | 8 | 2400源事件、28800模式事件 | 1.0349 | 1.3413 | 1735.9 |
| 四种finding+两种真实merge | 10 | 5000源事件 | 0.8413 | 1.0214 | 2078.4 |
| 代表事件单因素 | 390 | 390 | 1.3166 | 2.3002 | 1017.7 |

主要submission：初轮材料/代表事件1573761（initial_submission.json）、独立验证1573764、全量operations1573769、G4矩阵1573770、modes1573773、gzip定点retry1573775、完整operations统计1573776、旧metadata模拟统计1573779、最终完整审计1573780、mechanism表1573781、Ntuple验证1573783、sensor metadata重生成1573784、资源统计1573785、正确metadata/T1统计1573786、电子剩余尾审计1573787。每次选择均保存实时schedd load/queue及理由，实际选中bigbird20；bigbird21的EOS路径兼容性拒绝已记录并排除。各任务表和 `*_worker/submission.json` 为精确定位入口。

## 失败、修正与验证边界

| 失败步骤 | 实际现象与已知原因 | 处理和结果 |
| --- | --- | --- |
| 第一轮controls metadata | 新库属性不在genconf结果中；使用basename时载入run中的旧plugin | genconf绑定新build库绝对路径，冻结v2/v3；失败smoke不计入物理分母 |
| 补充G4退出134 | Pythia与FaserGeant4的静态G4ios符号抢占导致finalizer重复delete，gdb和minimal load pair复现 | Pythia链接使用 --exclude-libs,ALL；minimal pair和实际完整模拟均正常退出，没有跳过finalize |
| 诊断ParticleGun单位 | HepMC3默认GeV而填入MeV数值 | 显式MeV/mm；100GeV配置实际truth=99.9999969GeV，同seed plane距离3.815×10⁻⁸mm。历史RDO不据此追加1000倍修正 |
| EOS两个audit截断 | pdg-11_p100的T4_no_scatt_current和pdg13_p1000的T3_current在event52/110截断，初轮worker仍exit0；94/96 gzip完整 | 仅重跑两组，用scratch完成压缩、EOS .partial+flush/fsync、完整gzip读取验证再rename；simulation_verified_tasks.json关联最终96个合法模式与修正sensor metadata。底层截断原因未测量 |
| sensor event key | G4 event ID在300个Athena事件中恒为1，旧sidecar按此聚合不成立 | patch0005写Athena event/source_entry；相同随机种子重生成2400事件，81153个物理steps逐项相同。只重做metadata与统计，保留原reconstruction |
| 初版T1效率分母 | 以final ancestry筛good seed会遗漏finding失败或被prune的good seed | finding前确定母样本，final须来自合格seed；失败/换seed留在分母。回归测试覆盖该失效情形 |
| tracking-only Ntuple | 缺Veto/Preshower WaveformHits时，Ntuple内部reduced-charge计算解引用空handle | 专用验证入口提供空waveform containers；8事件13条track branch闭合。worker对非零reco返回码先保存失败，不再解析未闭合audit |
| SD增量构建的全量configure | installed run/src重复参与package扫描，CaloGeoModelUtils重复声明 | 使用既有单object build.make、两个既有link.txt与package install；两个事件probe及2400事件逐step同值验证加载的是新代码。一般全量configure问题仍记录在构建说明 |

所有控制保持实际材料选项可见。既有refit factory会将m_scattering/m_energyLoss传入fitter，成员的初始false并不代表生产refit关闭材料；新属性默认true保持原行为，全量观察配对提供实际输出验证。

## 复现入口与产物

按源码基底顺序应用patches/0001–0005，tracking与G4按变更组件编译。`patch_validation_final/validation.json`验证15个文件（12 tracking、2生成、1 sensor）与实际源码逐字节一致。S03冻结软件用于纯材料对照；controls_v3用于逐层与T矩阵，controls_v5仅加入Ntuple专用Python入口和sensor event metadata/G4库。

```bash
cd /eos/home-x/xcheng/FASER/geoRotation_test
georotation_stage04=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum
python3 scripts/rotate_material.py --source /eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/material-current-source.json --output "$georotation_stage04/reproduced-current22-rotated.json" --report "$georotation_stage04/reproduced-material-report.json"
python3 -m unittest discover -s tests -v
git diff --check
```

批量实际入口是 `scripts/submit_momentum.py --tasks <冻结task表> --library <冻结lib目录> --name <新任务名> --exclude-schedd bigbird21.cern.ch`，会在提交时重新比较schedd并记录选择。task表中的output必须为新的结果目录，已有输出用作只读输入。

| 内容 | 冻结输入/任务与产物 |
| --- | --- |
| S03母样本与材料发现 | stage03_baseline/audit_current；initial_tasks.json；material_consistent；merged_consistent；paired_material.json |
| 独立验证 | independent_inputs.json、independent_tasks.json；independent_merged；paired_independent.json；validation_complete |
| 逐运算和观察量一致性 | operations_tasks.json；operations_analysis；observer_integrity；mechanism_tables |
| 受控seed/KF | simulation_verified_tasks.json关联原96个模式与已验证sensor metadata；simulation_analysis_verified；simulation_sensor_validation.json |
| Calypso定向单因素、seed | case_controls_tasks.json；case_analysis；seed_ideal/seed_cases |
| 四种finding/实际merge | modes_tasks.json；mode_matrix；modes_analysis_complete |
| 实际Ntuple与资源 | ntuple_v2_tasks.json；ntuple_checks_v2；resources/summary.json |

所有xAOD/HITS/RDO临时文件在worker scratch中使用；EOS保留compact审计、必要的逐层诊断、输入定位、英文图和日志。已对7个字节相同的冻结库副本改用相对symlink，减少75212193 bytes；另清理旧sensor metadata统计产生的10个已替代衍生文件，减少16269088 bytes。旧simulation_analysis仅保留summary和冻结worker等来源记录，并用NOT_USED.md指向权威结果simulation_analysis_verified。整理记录见deduplication.json；未执行大型数据SHA扫描。

## 下一步计划

1. S05使用已冻结的current/strict配对，先在同一st1/VetoNu物理面定义tanθ_y与truth，分解reco-p窗口的选择迁移，再量化峰位置/宽度/占比与hit pattern、segment、材料和covariance的对应关系；验收为独立样本上峰结构的forward解释或最小修复。

2. S06先对已有μ±做共同truth p/位置/方向支持域和谱平衡，再决定追加规模；当前独立strict各1条wrong不足以回答正负不对称。保留原current、strict和完整G0/G1的区别；未测量的旋转净效应不以材料修正效果代替。

3. S07分别处理同truth segment组合、同点方向/field-integral seed、measurement C/branch策略、辐射后的KF模型、chart/refit起点和显式IFT refit策略。每次只改一项，报告p/angle/charge/效率/fakes/失败率并使用新的预留验证文件或随机种子；本阶段独立样本不再用作盲验证集。

验证：11项unittest、37个Python文件编译、两个新增shell入口语法和git diff --check通过；15份Markdown的本地链接检查通过；补丁0001–0005顺序应用后的15个源码文件与实际Calypso一致。数据覆盖与检查记录见[final_validation.json](../../geoRotation_test_results/stage04_momentum/final_validation.json)。
