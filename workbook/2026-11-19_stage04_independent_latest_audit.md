# S04 latest 动量审计的独立复核

## 阶段目标

复算 `110051` 的 seed、Station 2 和最终动量结果，检查五类样本的逐级修复证据，并审查 `geoRotation_calypso_latest_fix_v1` 中会影响找迹、refit 和输出的实现。[动量深查](2026-11-18_stage04_latest_momentum_deep_dive.md) · [修复验证](2026-11-18_stage04_latest_worktree_fixes_and_checks.md)

## 复核结论

**Station 2 的前向负尾与 good-seed 样本的 smoothed 左移都能复现。** 单独缩小初始 `q/p` 方差显著压低前向负尾。单独标定 seed 场常数让通过同一 truth gate 的轨迹群体改变，几乎消除这一群体的 smoothed 左移。固定 measurement path 时，两个干预对 smoothed 中心的影响都远小于原来约 12% 的偏移。结论适用于 `110051` 的三站正向、指定 truth matching 和当前 CKF 设置。五项修复整体上改善该样本的单轨和双轨效率，逐级累加实验给出每一步在前序组合下的边际效应。[独立统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/independent_summary.json)

**优先复核 refit 的协方差坐标与起点。** 输出转换把 ACTS 的 surface-bound 协方差直接放进 Trk curvilinear 参数，refit 则把这份矩阵复制到提前 10 mm 的新平面并再膨胀。现有动量对照约束了选中样本的中心偏移，协方差、角度和边缘轨迹的偏差仍待测量。[输出转换](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CreateTrkTrackTool.cxx:105) · [refit 起点和协方差](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/KalmanFitterTool.cxx:585)

## 工作内容与详细数据

在 `lxplus901.cern.ch` 上直接读冻结任务、运行清单、worker NPZ 和 ROOT，未重新提交批量任务。两轮分别核对 120、330 个 shard 的退出状态、输入和事件范围，覆盖 120,000、600,000 次重建。分组验证包含 `110051` 全 20 文件，标定后的决策量采用 `file01..file19` 的 19,000 事件。校验五样本轮次时核对各配置使用同一输入范围、`FASERNU-04 / OFLCOND-FASER-05` 的库配置和三站正向 mask。冻结库的 build manifest 以 `7fd40df8cd3b60e81b4095a4f84310c4a54700a0` 加当时工作区 diff 记录来源。这份 diff 与 `e47b85fd74ecb785dcc0ece6b8f61022075b2e7d` 的 12 个已跟踪文件 diff 一致。编译使用的审计头当时另加，结果顶层早期同名头缺少后来的中性粒子子链处理。仓库补丁目录的头与当前 commit 一致。[构建清单](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/software_fix_v1/software.json) · [审计头](../patches/LatestFastValidationAudit.h)

### S2 与最终动量

对每条满足 GitHub quality、purity `>0.7` 和 seed/track 同多数 barcode 的候选，取实际第一条 S2 measurement。相对 production HepMC 动量，`δp = 1/(|q/p| p_prod) − 1`。正向路径的配对键为输入文件、entry、barcode、按传播顺序的 `(geometry ID, cluster ID)` 全序列，去掉多重键，使用 baseline good-seed gate。区间为以源事件为单位的 2,000 次 Poisson bootstrap。这里的第一条 S2 measurement 视实际路径而定，层内序号随已测量状态计数。基线在 S2 L0 可见 696 个首序号、679 个第二序号以及 7 个第三序号，因此 `s0/s1` 汇总应描述为测量顺序；固定物理 wafer 需用 geometry ID 分组。[独立统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/independent_summary.json) · [历史归约](../analysis/latest_deep_reduce.py:243)

| 110051 独立 19 文件 | good-seed 有效轨迹 | First S2 predicted `δp < −0.10` | First S2 filtered 中位数 | First S2 smoothed 中位数 | 无 seed gate 的 smoothed 中位数 |
| --- | ---: | ---: | ---: | ---: | ---: |
| baseline | 702 | 35.90% | −74.01% | −11.67% | −2.44% |
| 仅降低 `covQOverP` | 833 | 0.60% | −16.84% | −11.20% | −1.71% |
| 降低 `covQOverP` 并修正 loc 协方差 | 799 | 0.50% | −5.63% | −11.08% | −1.60% |
| 仅改 seed B | 1,465 | 32.01% | −46.24% | −0.49% | −2.42% |
| 三项联合 | 1,666 | 0.60% | −1.20% | −0.50% | −1.71% |

固定路径后，仅降曲率方差的配对中位变化是 **+0.0285 pp**，95% CI `[+0.0244,+0.0341] pp`，共 675 对。与 loc 一起改时，607 对的变化是 **+0.0581 pp**，95% CI `[+0.0425,+0.0764] pp`。只改 B 时，692 对为 **−0.00080 pp**，95% CI `[−0.00131,−0.00011] pp`。这是前向宽曲率先验对原 smoothed 中心影响很小的直接证据，也给出了有限区间内仍然存在的效应和各组的路径保留范围。`1/|q/p|` 会把宽曲率分布的正负大模值折进低动量侧，S2 同层合并时混入了第一次 S2 测量更新之后的 predicted 状态。[独立统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/independent_summary.json)

![First S2 momentum residuals on held-out files](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/station2_independent.png)

全 20 文件、无 seed gate 的 17,125 条有限 smoothed 状态相对 production 中位数是 **−2.41%**。低记录损失组 6,624 条是 **+1.65%**。记录到大于 5% 损失的 7,376 条从 **−16.99%** 移至 **−0.52%**，S2 前记录损失超过 10% 的 3,566 条从 **−37.11%** 移至 **+0.33%**。只变换分母就产生这一改善，支持 production 与局部动量差解释了大量负尾。延续节点通过同 PDG 和 barcode 关联，以固定 `z_F` 上限处的最小记录动量作估计。它既无法识别全部连续损失，也未读取逐 step G4 过程标签。因而 `−9.2 pp` 是两个不同选择分布中位数的差，`−2.4 pp` 是其中一个总体中位数，不能作为两个独立因果量。PG μ 子约 `0.8–1.1 GeV` 连续损失仍是物理解释，未以逐 step 记账确认。[真值提取](../configs/latest_truth_vertices.py:24) · [延续链统计](../analysis/latest_deep_study.py:395)

原基线 725 条 good-seed 同真值轨迹最终 refit 成功率 **99.86%**，最终 `Trk::Track` 与 S2 smoothed 的逐轨迹动量相对差中位数 **+0.0334%**。该比较跨参考面，支持中心左偏在 refit 之前已经形成。第二份报告的列 `final_refit` 则读取 refit 自己的首条测量 `w0`，并在缺失时回退 CKF。同为 `C5` 的 23,556 条 ungated 同真值轨迹中有 **33** 条 refit 与 CKF 的首面不同、**27** 条用了 fallback。五参数 pull 代码对 refitted wafer ID 施加了同面要求，`q/p` 的 truth 仍取 `p_st1` 延续节点 proxy。另一个口径差异是 deep dive 的 seed gate 角度来自 ROOT seed writer，修复验证用 SiHit chord 作 seed 方向参考。相同文件上 `C0` 的 702 与 711 条分别对应这两种定义，不能直接混用。[refit 输出核对](../analysis/latest_validation_study.py:263) · [seed gate](../analysis/latest_validation_study.py:214)

### 五项修复、目录构建与探索组

- **Seed 组成。** 重算得到 52,338 条选中 seed，50,364 条有效，17,736 条 `δp < −0.5`（35.22%）。占尾部 59.76% 的是 10,599 条跨正负主粒子组合，32.74% 来自其余混合组合，7.27% 来自同粒子、记录硬损失组。低损失同粒子 seed 中位数 **+17.20%**，仍有 41 条低于 −0.5。21,681 条与前序 seed 恰好共享六个 cluster。三站中心的圆拟合确实无自由度，但实际排序使用 9 个相关的 segment 伪点以及 `x_F(z_F)` 残差。应以实测分类说明组合问题，避免把实际排序统计量描述为零自由度检验。[seed 排序与分数](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:326) · [独立统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/independent_summary.json)
- **C1、C2。** 旋转 seed 平面局部轴为 `(+y_F,+x_F)`，位置方差确实反序。修复验证中的 `110051` 串扰率 `32.54% → 13.53% → 9.37%` 复现。`C0→C1` 的 `74.90% → 71.57%` 指单轨效率，双轨是 `61.58% → 56.83%`。C2 同时改 `θ` 与 `q/p` 方差，单独降低 `q/p` 的证据来自前一轮。`0.283 → 1.076` 是已选中的同真值候选的新息标准差。创新协方差涵盖已拟合且被选择的 measurement，数值接近 1 不能独立标定先验 covariance。[逐级结果](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/independent_summary.json)
- **C3。** `0.473662 T` 能校正冻结场下、前向三站且满足闭合选择的总体 seed 标度。五类无噪声闭合的总体中位数差在约 `0.012%` 内。在 `110051` 的 `80≤r_max<100 mm` 子集，场闭合样本 150 条，中位数 **−0.159%**，16–84% 范围 **[−1.02%, +1.50%]**。`110033` 同孔径中位数 **+0.096%**，16–84% 范围 **[−0.80%, +1.79%]**。全孔径精度、其他 field IOV、IFT 与反向应用需各自验证。[场闭合独立分层](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/independent_summary.json)
- **C4、C5。** 保证非有限分数比较器满足严格弱序并让 bad score 变为有限惩罚，是明确改进。`C5` 修复代码在根号无效时回退既存方向，没有报告所称的 `99.9% r` 夹取。局部切线采用拟合圆与段位置，线性斜率拟合输入是段中心和 `±30 mm` 伪点。修复后 `110051` good-seed 同真值从 711 增至 2,740，串扰降至 **8.01%**。无 good-seed gate、purity cohort 的最终 `w0` 动量中心在 C0 已约为零。C5 优势主要体现在 seed 方向、样本产额、串扰和效率。[局部切线实现](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:596) · [冻结五样本汇总](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/summary.json)
- **H1。** 冻结样本上电子与 μ 子假设的动量中心及效率变化小，这一结果可复现。报告对机制的解释错误。对 Athena 24.0.41 的 `libActsCore.so` 反汇编确认该路径调用 `computeEnergyLossBethe`、`computeEnergyLossLandauSigmaQOverP` 和依赖 PDG 的 `computeMultipleScatteringTheta0`。换假设改变散射模型与质量相关量，本次 CKF / KF 路径没有启用额外的电子辐射能损或辐射协方差。因此 H1 不能回答启用电子辐射过程模型之后尾部会否改善。pull 变窄还混合了误差、选择与参考真值差异。[实际 ACTS 源码](https://github.com/acts-project/acts/blob/v32.0.2/Core/src/Propagator/detail/PointwiseMaterialInteraction.cpp)
- **S3。** 冻结任务将 `SeedMaxReducedChi2` 从默认 35 设为 `10⁹`，保留原 `>6` 重叠淘汰，排序变化本身足以改变最终候选池。对同一 19 文件的 36,812 条 fiducial truth，`C5→S3` 单轨效率 **78.99% → 78.52%**，按事件配对的差 **−0.465 pp**（95% CI `[−0.672,−0.257] pp`）。17,853 个双 truth 事件的配对效率 **67.59% → 67.11%**，差 **−0.476 pp**（95% CI `[−0.762,−0.190] pp`）。同时串扰 **8.01% → 7.51%**，而 GitHub 切选内低 purity 比例 **24.48% → 25.11%**。其他样本的单轨效率变化与零相容。这是串扰、候选产额与效率间的实测取舍。[排序和淘汰](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:334) · [配对区间](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/independent_summary.json)

### 材料、九层图与证据层级

两次冻结 JSON 与包内 gzip 解压逐字节相同，均为 11,866,364 bytes。重查 baseline smoke 的 `writer.npz`，43 条 tracks、731 条 measurement states、59 条 seeds 全字段相等，包含 NaN 位置。1,131 是审计轨迹状态数量，731 才是两边 writer 都保存的 measurement 状态。372 次 native lookup 涉及 9 个 surfaces。先前材料能损、散射开关的证据覆盖 1,000 事件控制样本和当时保留的同路径轨迹，支持该机制样本的负尾对这些开关不敏感。材料建模、其他几何分区和电子辐射无全量排除证据。[native smoke](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/smoke_validation.json) · [开关对照](2026-11-16_stage04_latest_momentum_origin.md#第一文件的机制控制)

逐 shard 用 uproot 对照两轮的 18 幅九层图，486 个 TH1 的 `fEntries` 与 **49,572** 个含 flow 分箱一致。用参考宏的 `t_eQOP`、100 bins、`[−1.2,1.2]`，两轮都附加了 purity / 可选 good-seed selection。表格里 First S2 使用 production truth，分别保留有限值分母。逐 bin 一致性核验的是已冻结宏输出的合并与图表，generation 口径仍需看 writer。比如 `C5` 无 gate 的 S2 L0 filtered 峰仍在约 **−0.828**，S1 L0 的 P/S 虽同在零附近的最高 bin，50,148 个状态的两条柱形分布相差半 L1 计数 3,791。报告中的“完全重合”超过了图中证据。[ROOT 核对](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/independent_root_checks.json)

## 重建源码发现

| 优先级 | 问题与影响 | 证据范围 |
| --- | --- | --- |
| 高 | **ACTS bound → Trk curvilinear → refit plane 缺少 Jacobian。** `CreateTrkTrackTool` 只转换 `q/p` 单位。对沿 `+x_A` 的轨迹，Trk curvilinear `U,V` 近似 `(+y_A,+z_A)`，CKF target 是 `(+y_A,−z_A)`，在一般 wafer 面还有 stereo 旋转。当前输出沿用原 `loc1` 与其他参数的相关项。refit 又把原矩阵赋给新平面，使用旧位置中心横坐标并提前 10 mm，缺少平移应有的角度项和传播协方差。默认路径另外先用 `SeedCovarianceScale=10` 放大对角元，再把整矩阵乘以 10，实际为对角 ×100、非对角 ×10。 | 源码和 native Trk frame 定义已核对。对效率与 pull 的影响未测量。[输出转换](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CreateTrkTrackTool.cxx:105) · [refit](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/KalmanFitterTool.cxx:362) |
| 高 | **C4 让 score 有限，病态 seed 仍可进入 CKF。** CircleFit 对近共线输入的 determinant 缺少有效性判断。`getChi2` 设 `1e30` 后，选择循环仍可选中该 seed，`fakeFit` 的半径及动量没有成套验证，`get_params` 还直接用 `charge / momentum` 和 `effDir.z()` 分母。零 NaN score 与初值全有限是两件事。`C4` 本身还保留 `getY` 展开多项式的大数相消。 | 代码可达，当前样本失败频率未归因到此缺陷。[圆拟合](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CircleFit.cxx:150) · [seed 选择](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:360) · [初值构造](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:699) |
| 高 | **不同顺位候选使用不同的最小 hit 数规则。** `MinNumberMeasurements` 存在于配置却未在选择处使用。CKF 先无条件接受队首，再丢弃其他 `nMeasurements<=12` 候选和共享超过 6 个 cluster 的候选。会接受首条短轨、丢掉独立的次条 12-hit 轨，同时依赖候选池顺序。`chi2` tie-breaker 也缺少有限性保护。 | 源码确认，尚未用同源事件量化各类短轨的效率和假轨率。[CKF 选择](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx:265) · [配置属性](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.h:103) |
| 中 | **原生 refit 在类型检查前解引用。** `getMeasurementsFromTrack` 先调用 `dynamic_cast` 结果的 `prepRawData()`，再判断指针。混合 measurement 输入或无 PRD 的输入可触发空指针。当前由 CKF2 自己构造的普通 SCT 测量通常符合预期。 | 路径和触发条件已确认。[refit 测量转换](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/KalmanFitterTool.cxx:553) |
| 中 | **跨事件静态缓存与隐式数值状态。** `s_spacePointMap` 从不清空，保存上一个事件的指针。`s_indexMap` 与 CKF 的 `TrackInfo::nClusters` 也是进程静态变量。已验证的串行重建主要透过 cluster source links 走主路径，旧 space point 指针未用于这条拟合。启用跨实例并发或以后复用这些指针时会影响正确性。`SeedStripFieldFit` 的轴向积分场缓存也只在首次有效条件时生成。 | 代码确认了生命周期，未测量跨 IOV 或并发影响。[seed 缓存](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:228) · [磁场缓存](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:70) · [CKF 候选位图](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.h:51) |
| 中 | **生产 4/3 站与反向使用相同三站调参。** 正常入口还启动带 IFT 和反向 CKF。圆拟合优先四段，仍用三站等效场常数。反向找迹后 refit 默认按正向传播。seed tool 的 `initialSurface()` 仍是原生 Z 法向，但 CKF2 会换成旋转面，因此直接影响其他使用该接口的调用者。 | 源码确认配置差异，多模式逐 event 物理性能未纳入两次 campaign。[生产入口](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Control/CalypsoExample/Reconstruction/scripts/faser_reco.py:300) · [seed 面](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:444) · [CKF 面](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx:204) |

### 生产入口实跑

使用 Athena 24.0.41、latest 原生 geometry 和冻结 `software_fix_v1` 的库与 Python overlay，通过标准 `faser_reco.py --isMC -g TI12MC04 -c OFLCOND-FASER-05 -n 20 -r independentAudit` 重建 `110051` 的 `file01 [0,20)`。进程映射确认实际加载库。20 个事件处理完成，退出状态 0，wall time 101.74 s，峰值 RSS 约 1.28 GiB。正反向各运行带 IFT 与去 IFT 两个模式。每个模式分别看到 50 个 seeds，正向选中 42 条轨迹，反向选中 40 条。两个反向模式各出现一条 4-measurement 候选，refit 明确拒绝并使用 CKF fallback。普通日志中的 “selected and re-fitted” 计数包含这些 fallback。[运行日志](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/production_fixed.log) · [实际加载库](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/production_loaded_libraries.txt)

正向模式各出现 6 次 seed 找迹失败，反向各 8 次，日志中的错误为达到 5,000 步上限。启动时 geometry 另有 5 条 `LayerArrayCreator` overlap 错误，位置在 `x_A ≈ −1800..−1609 mm` 的上游区域。同样五条错误见原 `C5` worker 日志。20 事件 smoke 确认了路径可运行和故障存在。几何 overlap 与具体失败 seed 的关系、IFT 有覆盖情况下的效率变化、反向动量偏移未测量。[原 worker log](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/runs/C5_local_tan/110051/file01/tracking.log.gz)

## 复现与产物

```bash
audit_results=/eos/home-x/xcheng/FASER/geoRotation_test_results
OPENBLAS_NUM_THREADS=1 python3 analysis/independent_latest_audit.py \
  --results "$audit_results" \
  --output "$audit_results/stage04_latest_validation/20261119_independent_audit"
python3 -m unittest discover -s tests -v
git diff --check
```

115 项现有单测通过。重画的首片 S2 图用固定分箱、原始计数，并保留 underflow / overflow。五页 PDF 附经过逐 bin 复核的原宏基线、C5 条件与无条件对照及场闭合图。[复算脚本](../analysis/independent_latest_audit.py) · [bin counts](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/station2_bin_counts.npz) · [逐图 PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/station2_independent.pdf) · [完整图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261119_independent_audit/report.pdf)

## 下一步计划

优先用同一 CKF 候选及 hit 序列完成 bound → curvilinear → bound 的 covariance Jacobian 闭合，比较现有和修正初始面在同一批轨迹上的 refit 成功率、动量和多参数 pull。独立检查共线 seed 的有限性和 CKF 最小 measurement 规则，按轨迹顺位、station 覆盖及共享 hits 多重数算效率。电子辐射量化需要固定 hit 和 truth 身份并记录 G4 局部动量或带过程标签的真值，再以同面 residual 约束模型。
