# 材料证据、软件起点与源码检查点

核查日期：2026-11-11。页码均指 PDF 物理页码，从 1 开始。历史图表只作为问题证据；本次没有重新测量 tracking 的发生率或分布。

## 三份报告分别说明了什么

| 材料 | 关键页 | 可支持的事实与尚待确认事项 |
| --- | --- | --- |
| [problem.pdf](../docs/problem.pdf)（4 页，Xiaocong Ai，2026-11-11） | 3 | seed 初始动量残差已有负侧结构；对 good seed 筛选后，CKF 逐层分布仍有低估。筛选为 `abs(residualP)<0.1`、同电荷、`abs(residualTheta)<0.0001`；角度单位与确切定义需由旧代码核对。图示并非六维参数逐项替换为 truth 的证明。用户另述曾强制 truth，本研究将把该旧配置与新的精确注入实验分别复现。 |
| 同上 | 4 | 旋转后零附近尖峰有所改善，要求确认剩余结构；提出 μ− 高于 μ+ 的 charge misID。st1 图标题实际为 `atan(Track_py0/Track_pz0)`，须区别于斜率与 ACTS 的 spherical theta。 |
| 同上 | 2 | 讨论 IFT alignment 参数异常；本次三问题研究中将 alignment/conditions 作为系统效应检查，完整 alignment 自由度优化另立研究阶段。 |
| [geometryRotation.pdf](../docs/geometryRotation.pdf)（24 页，Saimeng Zhang、Xiaocong Ai，2026-10-30） | 2–9 | 绕全局 y 旋转 90°，将 beam +z 变为 ACTS +x；涉及 detector surfaces、layers/volumes、binning、磁场/gradient、seed 参数/covariance、CKF target 和 refit。covariance 也作了调优。 |
| 同上 | 10–15 | 给出旧运行脚本及样本；报告中效率使用 truth fiducial、truth hits≥18，以及事件级双 track 条件。A′ forward/backward 单独指标有变化，merged 表现可接近；不同指标分母不能混用。 |
| 同上 | 16–18 | 110033、110051 的 CKF predicted/filtered/smoothed 动量残差负尾仍在；110051 更明显；100044 μ+ 的负尾改善。逐层 truth 的确切物理定义须核对 writer。 |
| 同上 | 19–20、23 | 100044 使用报告标注 500k events；VetoNu/st1 图有全谱和标为 50–250 GeV 的切片；A′ 附图使用不同切片（50–100 GeV）。图中的切选表达式实际使用 `Track_pz0`，不等价于按动量模长选择。 |
| 同上 | 24 | A′ seed residual 的 before/after 曲线和图中统计相同；μ+ seed 宽度改善。相同统计需要原始 ROOT/绘图对象复核。 |
| [Alma 9 Validation Overview.pdf](../docs/Alma%209%20Validation%20Overview.pdf)（27 页，S. Eley、M. D'Onofrio、C. Gwilliam） | 2–3、5 | 单 μ± particle gun，2D FLUKA 谱；比较 s0012-dev 与 s0012-r0019；fiducial truth at vertex 与 3 stations、`r_atMaxRadius<100`、恰好 1 Good Track。Good Track 只引用 Twiki，PDF 未给全部数值；位置 truth 使用 station 平均，与 reco reference plane 不同。 |
| 同上 | 16–18、23–25 | 新 tracking 在 wrong-sign q/p 尾部更多；μ− 的误判率高于 μ+；有 old-correct/new-wrong 事件。PDF 第 3 页提到用 truth PDG/reco charge matching，应复查旧分析是否对误判样本选择产生影响。 |
| 同上 | 19–22、27 | truth/reco 动量相关与 high-p 分布发生变化。第 27 页有宽范围 Track_p0 比较，不能由一张总谱图定位哪个拟合环节改变了动量。 |

### 可直接转录的历史数值

以下均为 PDF 图表的转录，未获得原始分子分母与误差：

| 指标 | before/old | after/new | 来源 |
| --- | --- | --- | --- |
| A′（10 MeV，ε=10⁻⁴）seed relative-p residual mean / std. dev. | −0.199088 / 0.592149 | −0.199088 / 0.592149 | rotation p24 |
| μ+ seed relative-p residual mean / std. dev. | 0.071982 / 0.302100 | 0.069308 / 0.237745 | rotation p24 |
| 100044 st1 全谱图 entries | 436859 | 395679 | rotation p20 / problem p4 |
| 同图低动量切片 entries | 176715 | 142671 | 同上；图表达式为 pz 切选 |
| μ+ charge misID | 0.101% | 0.197% | Alma p17–18 |
| μ− charge misID | 0.217% | 0.359% | Alma p17–18 |
| old correct & new wrong，μ+ / μ− | — | 0.122% / 0.223% | Alma p17–18；需核对配对分母 |

新 tracking 的 μ−−μ+ 差为 0.162 个百分点，旧 tracking 为 0.116 个百分点，均为表值相减；显著性未测量。图中的 entries 变化也不能直接当作旋转引起的效率损失，因为事件配对、文件覆盖和完整选择尚未闭合。

## 软件和现成构建

- 源码：`/eos/home-x/xcheng/FASER/geoRotation_calypso`，分支 `geoRotation_rebasedOfficialMaster`。
- HEAD：`88c457c1e4bf4154e67b96114ee8e7020d9b1a97`；前一旋转提交为 `73c97758`，其父提交 `c95170f2` 可作为建立同版未旋转对照的起点。
- `faser-common`：`1770d15df9311f6535c09909edf9ad64a5779e32`。
- 已有 `build/`、`run/`；环境为 Athena 24.0.41、x86_64-el9-gcc13-opt。实际库加载和组件配置在本次通过，见 [阶段报告](../workbook/2026-11-11_stage00_material_review.md)。
- 本地 `README.local.md`、`setup_environment.sh`、`build_local.sh`、`verify_local.sh` 为已有 helper；本次没有改动外部 Calypso 源码或重建软件。

### 当前入口与旧报告的差别

正式入口为 `Control/CalypsoExample/Reconstruction/scripts/faser_reco.py`。`--help` 已实际验证，支持 `--isMC`、`--skip`、`-n`、`--geom`、`--noForward`、`--noBackward`、`--noIFT`，没有 `-o`。

默认 runtype 是 TI12Data04，映射到 FASERNU-04 / OFLCOND-FASER-05；TI12MC03 映射到 FASERNU-03 / OFLCOND-FASER-03。不能凭文件名直接确定使用哪一组，S01 必须读取生产 metadata。`--cond` 分支打印了不存在的 `args.config`（第 190 行），该路径需先修复或用独立配置设置条件。

材料文件默认是 `/cvmfs/faser.cern.ch/repo/sw/database/DBRelease/current/acts/material-maps-alma9.json`；运行需记录 `current` 最终解析路径与文件版本。默认多个 CKF 实例均关闭 diagnostics；refit fitter 也单独关闭 diagnostics。旧脚本位置 `/afs/cern.ch/work/s/sazhang/public/FASER/xai_alma9/check_geometryRotation` 已确认可访问，但不视为当前源码的默认配置。

## 已确认存在的输入目录

共同前缀 `/eos/experiment/faser/data0/sim/mc24/`。本次仅检查目录与文件名；每文件 event 数、branch 内容、生产 geometry/conditions 和逐层真实动量均未测量。

| 用途 | DSID / 粒子 | 已存在目录（相对共同前缀） | 首轮选择依据 |
| --- | --- | --- | --- |
| 复现旋转角度图、单 μ 控制 | 100044 / μ+（5 mrad FLUKA E） | `particle_gun/100044/phy/s0013-r0022/`、`s0013-r0023/`；RDO 在 `particle_gun/100044/rdo/s0013/` | 历史图用 r0022 PHYS；重新 reconstruction 用对应 s0013 RDO，不能以新 PHYS 当原始 RDO |
| Alma 正电荷 | 100116 / μ+（2D FLUKA） | `particle_gun/100116/phy/s0012-dev/`、`s0012-r0019/`；RDO 在 `particle_gun/100116/rdo/s0012/`、`s0013/` | 先使用 s0012 与历史样本配对；s0013 作为独立系统检查 |
| Alma 负电荷 | 100117 / μ−（2D FLUKA） | `particle_gun/100117/phy/s0012-dev/`、`s0012-r0019/`；RDO 在 `particle_gun/100117/rdo/s0012/`、`s0013/` | 与 100116 在相同生产版本、共同谱域比较 |
| A′→e+e−，10 MeV，ε=10⁻⁵ | 110033 | `foresee/110033/rdo/master_dev/`，另有 dev、s0008、s0013 | 报告使用 master_dev；避免混入 MC22 或另一批 simulation |
| A′→e+e−，10 MeV，ε=10⁻⁴ | 110051 | `foresee/110051/rdo/master_dev/`、`s0008/` | seed/CKF 异常较明显，作为优先复现样本 |

已确认的示例：

```text
particle_gun/100044/rdo/s0013/FaserMC-MC24_PG_mupl_fasernu_5mrad_flukaE-100044-00000-00099-s0013-RDO.root
particle_gun/100116/rdo/s0012/FaserMC-MC24_PG_mupl_fasernu_2D_fluka-100116-00000-00009-s0012-RDO.root
particle_gun/100117/rdo/s0012/FaserMC-MC24_PG_mumi_fasernu_2D_fluka-100117-00000-00009-s0012-RDO.root
foresee/110033/rdo/master_dev/FaserMC-MC24_FS_Aee_010MeV_1Em5-110033-00000-00001-master_dev-RDO.root
foresee/110051/rdo/master_dev/FaserMC-MC24_FS_Aee_010MeV_1Em4-110051-00000-00009-master_dev-RDO.root
```

相同 DSID、文件序号或 TTree entry 不保证同事件。用 dataset/production ID + run/event + truth barcode 验证 old/new 和 before/after 交集；另存输入文件与 entry 位置供重放。

## 从输入到最终图的源码地图

下表路径相对 Calypso 根目录。行号属于当前 HEAD；后续修改按函数名和 commit 重新定位。

| 环节 | 路径、起始行/函数 | 需要回答的问题 |
| --- | --- | --- |
| cluster、space point、segment | `Control/CalypsoExample/Reconstruction/scripts/faser_reco.py:270`；`Tracker/TrackerRecAlgs/TrackerSegmentFit/src/SegmentFitAlg.cxx:692` | stereo/局部测量、segment 方向、协方差与质量选择是否已产生偏差；ghost removal 前后是否换了输入组合 |
| seed 组合和筛选 | `Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:110` | 4/3/2/1 station fallback、排序和 shared hits 是否改变 seed 组成 |
| seed 数值 | 同文件 `Segment:235`、`Seed:276`、`fakeFit:338`、`get_params:378`；`CircleFit.cxx` | ±30 mm fake points、固定 B=0.55 T、circle radius/center 符号、线拟合方向、近共线条件数与能标 |
| seed covariance | `Tracking/Acts/FaserActsKalmanFilter/python/CKF2Config.py:47` | 当前 inflation 为 [100,100,100,10000,100000]；分量位置及旋转前后物理等价性 |
| CKF 初始化与层更新 | `Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx:127`、`TrackFinderFunction.cxx:14` | target surface、measurement selector、navigation、predicted→filtered→smoothed、material 与 Jacobian |
| 候选排序与选中 | `CKF2.cxx:254` | nMeasurements 优先、χ² tie break、nMeasurements≤12/shared>6 的移除是否造成选择偏差 |
| ACTS→Trk | `CreateTrkTrackTool.cxx:104` | frame、GeV→MeV、plane/curvilinear 参数与 covariance 是否一致 |
| refit | `CKF2.cxx:321`；`KalmanFitterTool.cxx:343`、`getParametersFromTrack:581`；`TrackFitterFunction.cxx` | 改 reference plane、测量重建、outlier、能损/散射、二次 smoothing、失败 fallback |
| CKF/seed truth 与状态输出 | `RootTrajectoryStatesWriterTool.cxx:405`；`RootSeedWriteTool.cxx`；`RootTrajectorySummaryWriterTool.cxx` | truth 来源、state 有效性、单位、每层/每 track 对齐和写出对象是否真为最终 track |
| truth matching | `TrackTruthMatchingTool.cxx:16`、`identifyContributingParticles:77` | SDO/hit 的多数 barcode、匹配纯度、重复候选、次级粒子，不依赖待检验 charge |
| 合并 | `MergeTracksAlg.cxx`、`MergeTracksConfig.py`；Ntuple 配置 | 仅实际启用时审计 forward/backward 合并；记录选中来源与被替换候选 |
| 最终 Ntuple 与外推 | `PhysicsAnalysis/NtupleDumper/src/NtupleDumperAlg.cxx:1625`、`1858`、`1924` | Track_p0/py0/pz0 的 frame 与 surface、VetoNu 外推、Good Track 和失败占位值 |
| 几何、binning 和 alignment | `FaserActsGeometryContext.h:21`；`FaserActsDetectorElement.cxx:90`；`FaserActsLayerBuilder.cxx`、`CuboidVolumeBuilder.cxx` | surface 法向/local basis、construction/alignment transforms、gluing、binning、material bins 是否一致旋转 |
| 磁场 | `Tracking/Acts/FaserActsGeometry/FaserActsGeometry/FASERMagneticFieldWrapper.h:56` | 位置反变换、B 正变换、gradient 双侧变换、kT/T 单位和查询点 |

### 源码已经显示、需要优先做实验的具体问题

这些是实现事实或明确的不一致候选，尚未建立它们与三种总体异常的因果关系。

| 编号 | 实现事实 | 首个验证与判别 |
| --- | --- | --- |
| C1 | `geometryRotation=true` 位于 context；field wrapper、seed 参数、CKF target、refit 和 Ntuple 有不依赖该成员的旋转代码 | 单改成员不能构成完整 off 对照；用未旋转共同父版或覆盖全链的配置实现对照 |
| C2 | `CKF2.cxx:226` 定义 `nSeed=0`，循环内 `trackGroup=nSeed-1`，未见递增 | 诊断 seed ID 会无符号下溢；用显式 iseed 和 event/collection key 建立可验证 lineage，检查其下游影响 |
| C3 | `CircleFitTrackSeedTool.cxx:320` 单 station seed 固定 p=9999999、charge=+1 | 按 seed station 数统计两电荷 fallback、纠正率与最终误判；不能把默认正号当随机无偏先验 |
| C4 | seed 使用固定 `0.0231² mm²` measurement variance，refit 使用 `0.08²/12 mm²`；cluster 原 covariance 未直接使用 | 两者数值接近但来源与含义不同，分别对比 cluster 条件与 refit 更新，避免把它们称为量级差异 |
| C5 | ACTS→Trk 创建 CurvilinearParameters；refit 再以 plane 参数重建，并复制 covariance 的 5×5 block | 在同一物理 state 检验完整 Jacobian、local basis 和单位转换，不能仅凭对角线合理就视为闭合 |
| C6 | refit 的 `getParametersFromTrack` 对角乘 `SeedCovarianceScale`（配置 10），`fit` 中又将整矩阵乘 10 | 当前两步使对角累计×100、非对角×10；记录实际矩阵与 correlations，分开控制两步并测影响 |
| C7 | CKF writer 最后写 `selectedTracks`；final collection 写 refit 的 `trk2`，失败时写 `trk`；refit writer 默认 `noDiagnostics=True` | 需要同时捕获 refit，不能把 CKF states 图当作最终 Ntuple 的完整轨迹 |
| C8 | states writer 用 SiHit 起终点给方向，但 `truthP += particleLink()->momentum().rho()`；源码还有 “get truthQOP for each state” FIXME | 判断该 p 是否为生产顶点动量；没有逐层真实 p 时补充 G4 记录，分别画生成点与同面 residual |
| C9 | `NtupleDumperAlg.cxx:1653` 直接写 track momentum 为 px0/py0/pz0；`1936` 仍用外推 ACTS 的 [0]/[2]、[1]/[2] 生成角度 | 旋转关系要求 FASER y-slope 为 py_ACTS/px_ACTS；用一条已知方向验证实际执行路径和输出 frame，并核对 x 的反变换符号 |
| C10 | `s_spacePointMap` 为 static map，run 中更新；末尾只清空 `s_indexMap`；Segment 的 clPositions/clMomentum 收集在外层 state 循环内重复遍历参数 | 比较单独事件/连续事件/改变顺序及不同 CKF 实例；检查陈旧指针、重复诊断点与额外内存，测量它们是否影响实际计算 |
| C11 | Ntuple 将外推失败角度预填 `999999`；`Track_py0/pz0` 不一定对应严格的 st1 同一平面 | 区分成功外推与无效值；显式匹配 reference surface，完整报告丢失和失败计数 |
| C12 | states writer 在第 426/468 行由 HepMC momentum 构造 truthQOP，第 630/642 行存储 raw ACTS q/p 并直接相减，所读代码未显式作两者单位转换 | 用输入 HepMC 的 momentum unit 与真实值检查 t_eQOP/eQOP/residual/pull；即使图中 p residual 在外部修正，writer 自带 residual 仍需单独核验 |

## 尚待取得或核对的输入

1. 报告中原始 ROOT、绘图脚本、before/after 确切 commits、geometry/conditions/material map、forward/backward/merged 选择；保留 `Track_pz0` 切选的历史复现，同时新增正确物理 p 的对照。
2. 用户提到的强制 truth seed 测试配置：替换哪些分量、在哪个 surface、是否 smearing、covariance、charge、matching 与 candidate identity。
3. Alma Good Track 的 Twiki/原脚本；当前 Ntuple 的 goodTracks 判据还使用 `r<95 mm` 和可配置质量阈值，不自动视为历史 `r_atMaxRadius<100` 的同一选择。
4. 各输入是否存有同层真实动量，以及 hit/SDO 的可用性。若缺失，在最小 G4 样本记录 pre/post-step 动量、能损过程和穿面位置，不用 generator p 填充“逐层 truth”。

上述缺项影响历史数字的精确复现；frame、单位和候选来源等独立检查可以先开展。

## 方法参考及在本研究中的用途

2026-11-11 已通过 arXiv/Crossref 核对以下题名与摘要/书目信息；方法细节在实现相应实验时继续查阅。

- X. Ai et al., [A Common Tracking Software Project](https://arxiv.org/abs/2106.13593)：用于梳理 ACTS geometry、track states、propagation 和线程模型。落实到本地版本实际 API，作为逐层观测设计依据。
- R. Frühwirth, [Application of Kalman filtering to track and vertex fitting](https://doi.org/10.1016/0168-9002(87)90887-4)（1987）：用于拆分 prediction/update/smoothing、残差与 covariance 的因果链。直接适用于固定 hit 的 KF closure。
- J. A. Hartigan、P. M. Hartigan, [The Dip Test of Unimodality](https://doi.org/10.1214/aos/1176346577)（1985）：用于检验条件角度样本的多峰性。它只能描述分布结构，不能判定非物理；需结合 truth 条件谱、detector 分组与 forward model。

优先使用可解释的 curvature/update 统计和条件模型。Gaussian mixture、聚类或异常检测可以辅助归类，因果结论仍由固定输入的受控实验给出；若 electron bremsstrahlung 被证实是主要机制，再检索并比较适合非高斯能损的拟合方法。
