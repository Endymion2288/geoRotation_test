# S04：latest Calypso 动量残差左偏的源码审计、18 晶圆演化与 C++ 单因素因果解耦

## 阶段性结论与核心实测结果

在 `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest`（commit `7fd40df8cd3b60e81b4095a4f84310c4a54700a0`）及其包内旋转材料图 `Tracking/Acts/FaserActsGeometry/data/material-maps-alma9-rotated-rowmajor.json.gz`（原样解压 `11,866,364` bytes）下，对 DSID `110051` 全量 20 个 RDO 文件（`20,000` 源事件，6 组受控实验共 `120` 个 HTCondor 作业、`120,000` 重建事件）开展 C++ 源码级审计与单因素实验。旧 `good-seed` 测试中观测到的三个“左偏/负尾”现象——**seed 残差的 `−1.0` 尖峰**、**Station 2 `predicted` 的 `−0.92` 强左峰**、以及 **`good-seed` 条件下 `smoothed` 约 `−11.7%` 的左移**——由完全不同的源码机制产生；**Station 2 `predicted` 左偏不是 `smoothed` 左移的原因**：

1. **材料图核对结果（已排除材料查表致偏）**：`geoRotation_calypso_latest` 自带的 `material-maps-alma9-rotated-rowmajor.json.gz` 含 22 个 surface 条目（volumes `2/25/27/29` 的 12 个有源 approach 面采用行主序旋转矩阵 `[0,0,1, 0,1,0, -1,0,0]`，精确对应 `(x_A, y_A, z_A) = (z_F, y_F, −x_F)`；另含 10 个被动体空条目），与历史 AFS `material-maps-alma9-rotated.json`（12 条目、转置旋转序列及 2 mm 平移差）不同。372 次 native `Acts::BinnedSurfaceMaterial` 查表的局部坐标、bin 索引和 slab 厚度与解压 JSON 完全一致；关闭 CKF 能损或散射对同路径残差中位数的改变小于 `2×10⁻⁵`。
2. **全量 `52,338` 条 seed 的残差分布完全分解（`CircleFitTrackSeedTool.cxx`）**：
   - **干净单粒子低辐射 3-station seed（`N = 12,012`，占有效 truth seed 的 `23.85%`）无 `−1.0` 尖峰，呈窄峰但中心右移至中位数 `+17.20%`（无噪声场闭合 `+16.12%`）**：原因是 `CircleFitTrackSeedTool.h:92` 与 `.cxx:344` 在 `Seed::fakeFit(double B = 0.55)` 中硬编码均匀磁场 `0.55 T`，而真实 `FaserFieldTable_v2.root` 跨三站（`z_F ≈ 41, 1241, 2441 mm`，含无场间隙）的圆拟合有效场为 **`0.473662 T`**（`0.55 / 0.473662 − 1 = +16.12%`）。9 点 `±30 mm` 人造点（`.cxx:270-272`）相对 3 个站中心圆拟合（`.cxx:331`）的相对动量差中位数仅 `+5.70×10⁻⁵`（均值 `+0.093%`）。
   - **全体 seed 中 `35.22%` 堆积在 `< −0.5`（众数 bin `[−1.008, −0.984]`）的来源**：**`59.76%` 来自 `A' → e⁺e⁻` 正负电子跨轨迹混合组合**（`N = 20,404`，占有效 seed `40.51%`），**`32.74%` 来自次级粒子/非纯段/2-station 回退组合**（`N = 9,118`，占 `18.10%`），**`7.27%` 来自 `>5%` 硬轫致辐射单粒子轨迹**（`N = 8,830`，占 `17.53%`），干净低辐射单粒子仅占该负尾的 **`0.23%`**。跨轨迹假组合大量存活源于四处源码实现：
     1. `CircleFitTrackSeedTool.cxx:217-233` 对各站全部 `GhostBusters` 段做笛卡尔积；
     2. `.cxx:365-376` 的 `getChi2()` 在 `y(z)` 上对 3 个站中心恒有 0 自由度精确圆解，跨 `e⁺/e⁻` 毫米级大矢高组合给出小半径 `r`（`p_seed ≪ p_truth`，残差 `→ −1`）却仍获低 `χ²`；
     3. `.cxx:162` 去重条件 `((p.clusterSet & selected.clusterSet).count() > 6)` 使用严格大于 `6`，允许被选中的 seed 完整共享同一站的 **6 个 cluster**（全样本 `52,338` 条选中 seed 中有 **`21,681` 条（`41.42%`）** 与前序已选 seed 恰好共享 6 个 cluster，其中 `57.54%` 落在 `< −0.5`）；
     4. `.cxx:356` 的 `getY(z)` 未保护负根号判别式，当 `|z − cx| > r` 时产生 `NaN` `chi2`（全样本原始组合 `675,151` 条中 `1,971` 条 `NaN`，最终选中 `52,338` 条中仍有 **`750` 条 `NaN` `chi2`**、**`608` 条负判别式**），破坏 `.cxx:147-152` `std::list::sort` 的严格弱序。
   - **seed 方向与旋转平面协方差缺陷**：
     1. `.cxx:302` 取 `direction = positions[1] − positions[0]`（跨第一磁铁的 S1→S2 弦线而非 seed 面切线），造成电荷反号的系统性 `φ_A` 偏置（`e⁻` 中位数 **`+0.107 mrad`**，`e⁺` 中位数 **`−0.105 mrad`**；改用圆切线后分别为 `+0.005 mrad` 与 `+0.003 mrad`，标准差由 `0.960 mrad` 降至 `0.063 mrad`）。
     2. `.cxx:397-408` 将旋转版 seed `PlaneSurface` 局部轴设为 `eBoundLoc0 = +y_F`（精细弯转方向，`σ = 0.019 mm`）、`eBoundLoc1 = +x_F`（粗立体角方向，`σ = 0.9 mm`），但 `CKF2Config.py:62-63` 沿用未旋转版的反序赋值 `covLoc0 = 100 × 0.9² = 81.0 mm²`、`covLoc1 = 100 × 0.019² = 0.0361 mm²`。单因素修正该轴互换（`fix_cov_loc_swap`）后，`good-seed` 下 CKF 跨轨迹串扰率（seed 与最终轨迹 majority 粒子不一致）由 **`23.52%`（`223/948`）降至 `6.30%`（`54/857`）**，同真值有效轨迹数由 **`725` 增至 `803`（`+10.76%`）**。
3. **Station 2 `predicted` 强左峰（`−0.92`）的逐晶圆面根源**：
   - **前因 1（S1 末层提前拟合发散）**：`CKF2Config.py:56` 将 `initialVarianceInflation[4]` 从原生 `beforeGeo` 的 `1,000` 放大 100 倍至 **`100,000`**（`covQOverP = 10⁻³ GeV⁻²`，`σ(q/p) = 0.0316 GeV⁻¹ ≈ 33 |q/p_truth|`）。前 4 个晶圆面（`S1L0_s0..S1L1_s1`）仅约束 `(loc0, loc1, φ, θ)`，`q/p` 保持 seed 初值（`predicted < −0.1` 为 `0%–0.7%`）；到第 5–6 晶圆面（`S1L2_s0, S1L2_s1`）5 参数满秩后，`GainMatrixUpdater` 利用 S1 内部仅 `120 mm` 臂长的边缘场微小矢高强行更新 `q/p`，使 `S1L2_s1 filtered`（即 **`S2L0_s0 predicted`**，跨磁铁传播相对动量变化中位数仅 `−4.96×10⁻⁷`）产生 **`36.0%` 的 `< −0.1` 左尾与 `16.2%` 错电荷**。单因素恢复 `qop_1000` 后，`S2L0_s0 predicted < −0.1` 立即从 **`36.0%` 降至 `0.6%`**（`qop_1000_loc_fixed` 为 **`0.5%`**，错电荷 **`0.0%`**）。
   - **前因 2（`S2L0_s0` 首次更新放大 + 同层双晶圆合并 + `1/|q/p|` 折叠）**：经 Dipole 1（`~1.2 m`）传播至第 7 晶圆面 `S2L0_s0` 时，宽 `q/p` 先验映射为 `σ_inn = 3.094 mm`（`corr(q/p, φ) = −0.9977`，卡尔曼增益 `K_qop = −0.01009 GeV⁻¹/mm`），S1 短臂长斜率噪声产生 `std(r) = 0.604 mm` 新息并被更新器归入磁偏转 `Δ(q/p) = K r`，使有符号曲率误差标准差从 `1.10` 暴涨至 **`8.49`**（临时错电荷 **`43.1%`**），动量残差中位数跌至 **`−73.18%`**，随即作为第 8 晶圆面 **`S2L0_s1 predicted`（中位数 `−73.63%`）** 进入九层图。`RootTrajectoryStatesWriterTool` 与 `ckf_residual_perLayer.C` 将同层 `_s0` 与 `_s1` 合并绘制，且对称展宽的 `q/p = q/p_truth(1 + ε)` 经 `p = 1/|q/p|` 非线性变换（雅可比 `(1 + δ_p)⁻²`）把 `ε > 0` 与错号 `ε < −2` 两侧全部折叠进 `(−1, 0)`，在 `−0.92` 形成尖峰；直至第 13 晶圆面 **`S3L0_s0`** 形成 `2.4 m` 三站长基线后，错电荷骤降至 **`0.7%`**。
4. **Station 2 `predicted` 与 `smoothed` 的因果解耦，及 `smoothed` 左移的两项精确来源**：
   - **因果独立性证明**：`CKF2.cxx:277` 对选中路径独立运行逆向 `Acts::GainMatrixSmoother`，起点 `S3L2_s1` 已具完整三站精度（`σ(q/p) ≈ 0.13 |q/p_truth|`，初始宽先验信息权重 `< 10⁻⁴`）。在 `qop_1000_loc_fixed` 中，`S2L0_s0 predicted` 左尾已从 `36.0%` 消除至 `0.5%`，`S2L0_s0 filtered` 中位数从 `−73.18%` 升至 `−5.63%`，但同路径配对 `First S2 smoothed` 残差仅移动 **`+0.057%`**（19 个独立文件总体中位数仍为 **`−11.08%`**，对比 baseline **`−11.67%`**）。因此 **S2 `predicted` 左偏并不导致 `smoothed` 左偏**。
   - **来源一（占 `−9.2%`，条件选择偏差）**：在 C++ `CircleFitTrackSeedTool` 中单因素设置 `SeedBFieldTesla = 0.473662`（仅由 `file 0 [0,100)` 确定性场闭合标定，**严格排除 `file 0` 在其余 19 个独立文件、`19,000` 事件上评估**）：
     - 若固定原 baseline `good-seed` 已选中的同路径轨迹（`692/702` 条配对），`First S2 smoothed` 中位数差值仅 **`−0.001%`**（`−11.71% → −11.64%`），证明改 `SeedBFieldTesla` 不改变给定 hits 上的 CKF 拟合值；
     - 若对各自重建出的 C++ seed 施加**完全未修改的原始 `good-seed` 门限**（`|p_seed/p_truth − 1| < 0.1`、同号、`|Δθ| < 10⁻⁴ rad`），19 个独立文件的同真值通过轨迹数从 **`702` 翻倍增至 `1,465`（`seed_b_04737`）和 `1,666`（`combined_phys_fix`）**，`First S2 smoothed` 中位数从 **`−11.67%` 移至 `−0.49%` 与 `−0.50%`**（`< −0.1` 比例从 **`61.82%` 降至 `4.78%` 与 `3.78%`**）。原 `B = 0.55 T` 使未辐射正常轨迹的 `p_seed` 落在 `+16.1%` 而被 `±10%` 窗口剔除，专门筛出了三站矢高偏大（含磁铁前 `~10%–15%` 辐射或正曲率涨落）的子样本。
   - **来源二（占 `−2.4%`，无门限总体的电子真实轫致辐射 vs 产生点真值定义）**：不加任何 seed 门限时，全 20 文件同真值轨迹（`N = 17,125` 有限状态）相对产生点 HepMC 动量 `p_prod` 的 `First S2 smoothed` 中位数为 **`−2.41%`**。按 HepMC 同 PDG 延续顶点分层：
     - **无辐射电子（全径迹 `p_min / p_prod ≥ 0.995`，`N = 6,624`）相对 `p_prod` 的中位数为 `+1.65%`，不存在负偏**；
     - **轻度辐射（`0.95 ≤ p_min / p_prod < 0.995`，`N = 2,795`）中位数为 `+0.51%`**；
     - **硬轫致辐射（`p_min / p_prod < 0.95`，`N = 7,376`）相对 `p_prod` 中位数为 `−16.99%`，而改用进入 S2 前局部真值上界 `p_pre_S2` 作分母后中位数立即回到 `−0.52%`**（S2 前损失 `>10%` 的 `3,566` 条轨迹从相对 `p_prod` 的 **`−37.11%`** 回到相对 `p_pre_S2` 的 **`+0.33%`**）。
   - **CKF `smoothed` 与后续 `KalmanFitterTool1` refit 对比**：`good-seed` 同真值 725 条轨迹中 refit 成功率 `99.86%`，CKF `First S2 smoothed` 中位数 `−11.71%`，最终 `Trk::Track` 中位数 `−11.39%`，逐轨迹相对动量差 `(p_final − p_ckf)/p_ckf` 中位数仅 **`+0.033%`**（无门限 `17,126` 条中位数 **`+0.012%`**），排除最终 refit 引入额外左偏。

---

## 阶段目标

针对用户在 `workbook/2026-11-13_stage04_good_seed_comparison.md` 与 `workbook/2026-11-16_stage04_latest_momentum_origin.md` 基础上提出的深入源码核查要求，在保持 `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest` 原始目录干净的前提下完成四项任务：

1. 核实 `geoRotation_calypso_latest` 实际使用的旋转材料图文件与 native `BinnedSurfaceMaterial` 解码行为，排除材料图错配；
2. 深入 `SegmentFitAlg`、`GhostBusters`、`CircleFitTrackSeedTool`、`RootSeedWriteTool` 源码，完整量化全量 `52,338` 条 seed 残差分布（尤其是 `35.2%` 的 `−1.0` 尖峰、`+16%` 干净峰、cluster 共享、方向弦线近似、`NaN` `chi2` 及旋转平面 `covLoc0/covLoc1` 轴互换）；
3. 沿正向 CKF 全部 **18 个独立 SCT 晶圆面**（`S1L0_s0` 至 `S3L2_s1`）逐面追踪 `predicted → filtered → smoothed`、新息标准差 `σ_inn`、卡尔曼增益 `K_qop`、错电荷率及 `p = 1/|q/p|` 折叠效应，查明 Station 2 `predicted` 左偏的精确起点与机制；
4. 在隔离的 `geoRotation_calypso_latest` 编译副本中实现真实的 C++ `SeedBFieldTesla` 与旋转协方差单因素开关，在 DSID `110051` 全 20 文件（及排除标定文件后的 19 个独立文件）上检验 Station 2 `predicted` 与 `smoothed` 的因果关系，并拆解无门限总体的剩余 `−2.41%` 偏移与最终 refit 影响。

---

## 工作内容与实验设计

### 1. 环境、隔离源码树、材料图与不变性验证

结果根目录为 `/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1`（下文简称 `R2`）。

| 项目 | 实际配置与路径 |
| --- | --- |
| 基准源码（保持 clean） | `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest`（HEAD `7fd40df8cd3b60e81b4095a4f84310c4a54700a0`） |
| 隔离诊断 worktree | `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_momentum_v2`（同一 commit，稀疏检出 `FaserActsKalmanFilter`） |
| 仓库可复用补丁 | [0049_latest_seed_and_cov_audit.patch](../patches/0049_latest_seed_and_cov_audit.patch)、[LatestSeedAndStateAudit.h](../patches/LatestSeedAndStateAudit.h) |
| 编译产物与配置 | [libFaserActsKalmanFilter.so](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/software_deep_v1/lib/libFaserActsKalmanFilter.so)（Athena `24.0.41`，`x86_64-el9-gcc13-opt`，ACTS `v32.0.2`，链接 `geoRotation_calypso_latest/run/lib`） |
| 旋转材料图来源 | `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsGeometry/data/material-maps-alma9-rotated-rowmajor.json.gz` |
| 冻结解压材料文件 | [material-latest-rowmajor.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/material-latest-rowmajor.json)（`11,866,364` bytes，22 个 surface 条目） |
| 输入样本 | `/eos/experiment/faser/data0/sim/mc24/foresee/110051/rdo/master_dev/` 全部 20 个 RDO 文件，每文件 `[0, 1000)`，共 `20,000` 源事件（以 `file_index + source_entry` 为唯一事件键） |
| 坐标系与单位 | 原始 FASER 系 `F: (x_F, y_F, z_F)`；旋转 ACTS 系 `A: (x_A, y_A, z_A) = (z_F, y_F, −x_F)`；动量 `GeV`，曲率 `q/p` `GeV⁻¹`，位置 `mm`，角度 `rad/mrad`，磁场 `T` |
| 基准筛选与九层宏 | candidate-isolated smoothing 开启；GitHub 原宏 `ckf_residual_perLayer.C`（commit `33ca2281b961935456dfd4dc9fc06235e24b551e`）三项筛选 `chi2Sum/NDF < 20 && rMax_smt < 90 && NDF > 9`；诊断附加 `strict final SDO purity > 0.7`；`good-seed` 门限为 `|p_seed/p_truth − 1| < 0.1`、`qop_seed · qop_truth > 0`、`|\theta_seed − \theta_truth| < 10⁻⁴ rad`（不含 `φ` 筛选） |

**补丁内容与 20 事件逐位闭合检查**（[smoke_validation.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/smoke_validation.json)）：
- 在 `CircleFitTrackSeedTool.{h,cxx}` 中增加可选属性 `SeedBFieldTesla`（默认 `0.55`，行为与原生完全相同）及只读 seed 取证写出（保存每个选中 seed 的 3 站 cluster 多数粒子、纯度、9 点 `fakeFit` 与 3 站中心 `fit` 半径对比、弦线/段方向/圆切线三种角度、与前序选中 seed 的最大共享 cluster 数、`getY` 负根号判别式计数及原始组合总数）；
- 在 `CKF2.{h,cxx}` 与 `TrackFinderFunction.{h,cxx}` 中记录真实 `iseed` 祖先索引、全部 18 个晶圆面的 `predicted/filtered/smoothed` 六维参数及协方差对角/相关系数、新息与卡尔曼增益，以及 native `BinnedSurfaceMaterial` 查表记录；
- 在 `file 0 [0, 20)` 上对比未打补丁的 `geoRotation_calypso_latest/run/lib`（`native_smoke`）与默认参数下的 `software_deep_v1/lib`（`audit_smoke`）：59 条 seeds、43 条 writer tracks、1,131 个轨迹状态的全部标量、向量、偏移与 `NaN` 位置 **100% 逐位相等（`exact_equal_nan: true`）**；372 次 native 材料查表全部精确闭合；731 次测量更新满足 `|(q_f − q_p) − K r| ≤ 1.39×10⁻¹⁷ GeV⁻¹`。补丁通过 `git apply --check -R` 反向验证。

### 2. 六组全量 20 文件单因素与组合控制矩阵

所有 6 组均在相同 20 个 RDO 文件（每组 `20,000` 源事件，共 `120` 个作业）上运行完整 Calypso 重建（含 `SegmentFitAlg → GhostBusters → CircleFitTrackSeedTool → CKF2 → KalmanFitterTool1`）：

| 组名 | `SeedBFieldTesla` [T] | seed 面 `(covLoc0, covLoc1)` [mm²] | `covQOverP` inflation (`covQOverP` [GeV⁻²]) | 评估文件范围 | 设计目的 |
| --- | ---: | --- | --- | --- | --- |
| `baseline` | `0.55`（原生） | `(81.0, 0.0361)`（原生反序） | `100,000`（`10⁻³`） | 全 20 文件 & 独立 19 文件 | `geoRotation_calypso_latest` 原生基线 |
| `fix_cov_loc_swap` | `0.55` | **`(0.0361, 81.0)`（修正轴对应）** | `100,000`（`10⁻³`） | 全 20 文件 & 独立 19 文件 | 单因素检验旋转 seed 平面 `loc0=y_F, loc1=x_F` 协方差反序的影响 |
| `qop_1000_only` | `0.55` | `(81.0, 0.0361)` | **`1,000`（`10⁻⁵`，对齐 `beforeGeo`）** | 全 20 文件 & 独立 19 文件 | 单因素检验 `100,000` 超宽曲率方差对 S1L2 与 S2 的放大作用 |
| `qop_1000_loc_fixed` | `0.55` | **`(0.0361, 81.0)`** | **`1,000`（`10⁻⁵`）** | 全 20 文件 & 独立 19 文件 | 同时消除协方差轴反序与 S1L2 曲率发散，检验 S2 `predicted` 恢复后 `smoothed` 是否仍左偏 |
| `seed_b_04737` | **`0.473662`** | `(81.0, 0.0361)` | `100,000`（`10⁻³`） | **严格独立 19 文件（排除 `file 0`）** | 单因素在 C++ seed 中使用有效场积分常数，保留宽先验与 S2 发散，检验原始 `good-seed` 下 `smoothed` 是否归零 |
| `combined_phys_fix` | **`0.473662`** | **`(0.0361, 81.0)`** | **`1,000`（`10⁻⁵`）** | **严格独立 19 文件（排除 `file 0`）** | 同时修正 seed 有效场、局部轴协方差与曲率膨胀，检验全链路联合表现 |

---

## 详细数据与源码机制分析

### 1. 源码审计一：全量 `52,338` 条 seed 残差分布的精确构成

阅读 [CircleFitTrackSeedTool.cxx](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx) 与 [RootSeedWriteTool.cxx](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/RootSeedWriteTool.cxx) 定位以下实现机制：

- **站间笛卡尔积与 0 自由度圆拟合**（[CircleFitTrackSeedTool.cxx:217-233](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx#L217-L233)、[338-376](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx#L338-L376)）：`DSID 110051` 为 `A' → e⁺e⁻` 双轻子事件，全样本平均每事件在 S1/S2/S3 产生 `(3.75, 4.82, 3.99)` 个原始段、`GhostBusters` 后仍保留 `(2.74, 3.40, 3.20)` 个段（`20,000` 事件中有 `18,179` 个事件至少一站含 `≥3` 个段），共生成 **`675,151` 个原始 seed 组合**。3 个站中心在弯转面 `(z_F, y_F)` 恒能唯一确定一个精确圆，跨 `e⁺/e⁻` 组合因正负电子向相反方向偏转而形成毫米级假矢高，拟合出极小半径 `r`（`p_seed = 0.3 × 0.55 × r × 10⁻³ ≪ p_truth`，故 `(p_seed − p_truth)/p_truth → −1`），但 9 个 `fakePositions`（每站 `±30 mm` 共线点）相对该圆的 `getChi2()` 依然很小。
- **去重阈值 `> 6` 放行整站共享**（[CircleFitTrackSeedTool.cxx:161-163](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx#L161-L163)）：每个 3 层双面 SCT 站恰好含 `6` 个 cluster。代码写为 `((p.clusterSet & selected.clusterSet).count() > 6)`，导致两条共享完整 1 个站（6 个 cluster）的组合**不会被去重剔除**。全样本 `52,338` 条选中 seed 中，与前序已选 seed 共享 0 个 cluster 的为 `25,997` 条，共享 `1–5` 个的为 `4,660` 条，**恰好共享整站 6 个 cluster 的高达 `21,681` 条（`41.42%`）**，其中 **`57.54%` 落在 `< −0.5` 负尾**。
- **`getY()` 负判别式导致 `NaN` `chi2` 破坏排序**（[CircleFitTrackSeedTool.cxx:147-152](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx#L147-L152)、[355-358](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx#L355-L358)）：`std::sqrt(-cx*cx + 2*cx*z + r*r - z*z)` 在 `±30 mm` 外推点超出小圆水平跨度 `|z − cx| > r` 时返回 `NaN`。在 C++ `std::list::sort` 比较器中 `NaN < x` 与 `x < NaN` 均为 `false`，违反严格弱序，使 **`750` 条 `NaN` `chi2` 组合（含 `608` 条负根号判别式）** 直接混入最终选中的 `52,338` 条 seed。
- **固定 `B = 0.55 T` 与跨磁铁弦线方向**（[CircleFitTrackSeedTool.h:92](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.h#L92)、[.cxx:302](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx#L302)、[.cxx:344](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx#L344)）：下表给出全 20 文件有效 truth seed（`50,364` 条；另有 `1,974` 条因 `RootSeedWriteTool` 无匹配主粒子而无 truth）按物理来源的互斥分解（完整数值见 [summary.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/summary.json)）：

| Seed 互斥物理类别（全 20 文件，ACTS/FASER 标量 `p` 相对 `RootSeedWriteTool` 传播真值） | 数量 | 占有效 seed 比例 | 中位数 `(p_seed−p_truth)/p_truth` | 16% — 84% 分位数 | `< −0.5` 比例 | **占全样本 `< −0.5` 尖峰的份额** | 众数 bin |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| 全部有效 truth seed | 50,364 | 100.00% | −1.31% | −89.84% — +39.36% | 35.22% | 100.00% | `[−1.008, −0.984]` |
| 1. 纯单粒子 3 站、低辐射（`p_min/p_prod ≥ 0.95`） | 12,012 | 23.85% | **+17.20%** | +6.05% — +38.79% | **0.34%** | **0.23%** | **`[+0.144, +0.168]`** |
| 2. 纯单粒子 3 站、硬轫致辐射（`p_min/p_prod < 0.95`） | 8,830 | 17.53% | −2.15% | −46.50% — +20.24% | 14.61% | 7.27% | `[+0.072, +0.096]` |
| 3. 主 `e⁺/e⁻` 正负电子跨轨迹混合 seed | 20,404 | **40.51%** | −56.64% | −93.51% — +55.03% | 51.95% | **59.76%** | **`[−1.008, −0.984]`** |
| 4. 次级粒子 / 段纯度不足 / 2 站回退 seed | 9,118 | 18.10% | −75.84% | −98.31% — +78.95% | 63.68% | **32.74%** | **`[−1.008, −0.984]`** |

对上述 `12,012` 条干净单粒子低辐射 3 站 seed 的角度估计器对比（ACTS 系 `φ_A, θ_A` 相对 seed 面传播真值，单位 `mrad`）：

| 方向估计量 | 样本 | 中位数 [mrad] | 标准差 [mrad] | 16% — 84% 分位数 [mrad] |
| --- | ---: | ---: | ---: | --- |
| 原生 S1→S2 弦线 `Δφ_A`（`e⁻` 子集） | 7,632 | **+0.107** | 0.969 | +0.053 — +0.332 |
| 原生 S1→S2 弦线 `Δφ_A`（`e⁺` 子集） | 4,380 | **−0.105** | 0.801 | −0.387 — −0.046 |
| 圆拟合在 seed 面切线 `Δφ_A`（`e⁻` 子集） | 7,632 | **+0.005** | **0.068** | −0.026 — +0.037 |
| 圆拟合在 seed 面切线 `Δφ_A`（`e⁺` 子集） | 4,380 | **+0.003** | **0.055** | −0.024 — +0.032 |
| 原生 S1→S2 弦线 `Δθ_A`（全体干净 seed） | 12,012 | +0.001 | 0.651 | −0.597 — +0.606 |
| 3 站 `x(z)` 直线拟合 `Δθ_A`（全体干净 seed） | 12,012 | +0.001 | **0.350** | −0.312 — +0.314 |

![Seed population forensics](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/seed_population_forensics.png)

---

### 2. 源码审计二：旋转 seed 平面 `covLoc0/covLoc1` 轴互换与全部 18 晶圆面演化

在 [CircleFitTrackSeedTool.cxx:391-415](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx#L391-L415) 中，全局旋转矩阵 `trafo` 将 `(x_F, y_F, z_F)` 映射为 `(z_F, y_F, −x_F)`，同时构造旋转后的 seed `PlaneSurface`：
- `rotation.col(0) = (0, 1, 0)_A = +y_F` → **`eBoundLoc0` 对应精细弯转方向 `y_F`（strip pitch 分辨率 `~0.019 mm`）**；
- `rotation.col(1) = (0, 0, −1)_A = +x_F` → **`eBoundLoc1` 对应粗立体角方向 `x_F`（stereo 分辨率 `~0.9 mm`）**。

然而 [CKF2Config.py:46-66](file:///eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/python/CKF2Config.py#L46-L66) 直接保留了未旋转 `master` 分支（`loc0 = x_F, loc1 = y_F`）的交叉赋值，并把曲率方差膨胀因子从 `1,000` 改为 `100,000`：

```python
sigma_loc0 = 1.9e-2   # fine pitch (y_F)
sigma_loc1 = 9e-1     # coarse stereo (x_F)
initial_variance_inflation = [100, 100, 100, 10000, 100000]
track_seed_tool.covLoc0 = initial_variance_inflation[0] * sigma_loc1 * sigma_loc1  # 81.0 mm^2 -> 赋给了精细轴 y_F!
track_seed_tool.covLoc1 = initial_variance_inflation[1] * sigma_loc0 * sigma_loc0  # 0.0361 mm^2 -> 赋给了粗轴 x_F!
track_seed_tool.covQOverP = initial_variance_inflation[4] * sigma_qop * sigma_qop  # 1e-3 GeV^-2 (sigma = 0.0316 GeV^-1)
```

这导致两个直接后果：（1）seed 在粗立体角轴 `x_F` 上先验标准差仅 `0.19 mm`（比真实段误差小约 `4.7` 倍），而在弯转轴 `y_F` 上高达 `9.0 mm`，在双径迹紧邻事件中容易捕获伴随轻子测量点；（2）`covQOverP = 10⁻³ GeV⁻²` 在第 5 个晶圆面（`S1L2_s0`）满秩时立即触发单站短基线曲率拟合。

下表列出全 20 文件、`good-seed` 且同 seed/track 真值的轨迹在全部 **18 个独立晶圆面**（按正向 CKF 实际访问顺序）的演化，对比 `baseline`（`N=725`）与 `qop_1000_loc_fixed`（`N=826`）：

| 序号 | 晶圆面 | `baseline` 有效面数 | `baseline` `prt` 中位数 | `baseline` `prt < −0.1` | `baseline` `flt` 中位数 | `baseline` `flt` 错电荷 | `baseline` `σ_inn` 中位数 [mm] | `baseline` `pull` RMS | `qop_1000_loc_fixed` `prt < −0.1` | `qop_1000_loc_fixed` `flt` 错电荷 | `qop_1000_loc_fixed` `smt` 中位数 |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `S1L0_s0` | 721 | +2.74% | 0.0% | +2.74% | 0.0% | 9.584 | 0.005 | 0.0% | 0.0% | −11.12% |
| 2 | `S1L0_s1` | 705 | +2.58% | 0.0% | +2.56% | 0.0% | 0.276 | 0.096 | 0.0% | 0.0% | −11.21% |
| 3 | `S1L1_s0` | 721 | +2.74% | 0.0% | +2.47% | 0.0% | 1.250 | 0.706 | 0.0% | 0.0% | −11.10% |
| 4 | `S1L1_s1` | 695 | +2.47% | 0.7% | +2.55% | 0.0% | 0.057 | 0.530 | 0.0% | 0.0% | −11.03% |
| 5 | **`S1L2_s0`** | 723 | +2.55% | 1.5% | **+6.65%** | **8.3%** | 0.049 | 0.834 | **0.0%** | **0.0%** | −11.13% |
| 6 | **`S1L2_s1`** | 699 | +7.47% | **26.0%** | **+6.13%** | **16.2%** | 0.046 | 0.728 | **0.3%** | **0.0%** | −11.23% |
| 7 | **`S2L0_s0`** | 719 | **+4.80%** | **36.0%** | **−73.18%** | **43.1%** | **3.094** | 0.189 | **0.5%** | 27.1% | −11.16% |
| 8 | **`S2L0_s1`** | 700 | **−73.63%** | **78.3%** | −69.83% | 39.9% | 0.456 | 0.532 | 48.7% | 27.6% | −11.23% |
| 9 | `S2L1_s0` | 722 | −69.62% | 78.0% | −63.69% | 38.6% | 0.035 | 1.161 | 56.9% | 27.0% | −11.08% |
| 10 | `S2L1_s1` | 692 | −64.50% | 76.3% | −53.68% | 35.7% | 0.035 | 1.101 | 55.9% | 23.2% | −11.05% |
| 11 | `S2L2_s0` | 724 | −53.89% | 72.5% | −38.91% | 40.7% | 0.033 | 0.908 | 52.2% | 22.0% | −11.13% |
| 12 | `S2L2_s1` | 699 | −38.53% | 62.1% | +8.27% | 38.9% | 0.032 | 0.806 | 51.7% | 13.9% | −11.23% |
| 13 | **`S3L0_s0`** | 721 | +5.19% | 47.2% | **−10.14%** | **0.7%** | 0.587 | 1.254 | 41.0% | **1.1%** | −11.16% |
| 14 | `S3L0_s1` | 691 | −10.13% | 50.4% | −11.21% | 0.4% | 0.050 | 0.647 | 49.3% | 0.8% | −11.23% |
| 15 | `S3L1_s0` | 723 | −11.38% | 55.2% | −11.32% | 0.6% | 0.031 | 1.062 | 52.2% | 1.0% | −11.10% |
| 16 | `S3L1_s1` | 680 | −11.33% | 56.0% | −11.86% | 0.3% | 0.031 | 1.064 | 54.8% | 0.6% | −11.07% |
| 17 | `S3L2_s0` | 720 | −11.82% | 61.3% | −11.58% | 0.6% | 0.028 | 0.993 | 56.1% | 1.0% | −11.12% |
| 18 | `S3L2_s1` | 687 | −11.53% | 59.4% | **−11.65%** | **0.6%** | 0.028 | 0.930 | 57.0% | **0.5%** | −11.08% |

**物理与数值机制解析**：
1. **S1 末层（晶圆 5–6）**：恢复 `qop_1000` 后，刚进入 Station 2 的 `S2L0_s0 predicted` 左尾从 `36.0%` 直接消失至 `0.5%`，证实入口左尾 100% 源自 `100,000` 倍曲率膨胀因子导致的 S1 短臂长过早拟合，与站间磁铁传播或材料无关。
2. **S2 首层双面合并与曲率折叠（晶圆 7–8）**：即使在 `qop_1000_loc_fixed` 下，仅靠 S1+S2 两站测量时，S1 站内 `120 mm` 基线的出射角分辨率（`σ_φ ≈ 0.2 mrad`）与 `~1 TeV` 轨迹跨 `1.2 m` 磁铁的偏转角（`~0.24 mrad`）同量级，因此两站前向滤波在 `S2L0_s1..S2L2_s1` 仍具 `σ(q/p)/|q/p_truth| ~ 2–3` 的固有曲率不确定度（`baseline` 更达 `8.49`）。对称的宽曲率分布经 `p = 1/|q/p|` 折叠后必偏向负动量残差侧；而九层图将 `S2L0_s0`（更新前）与 `S2L0_s1`（首次 S2 更新后）画入同一条绿色 `predicted` 曲线，形成视觉上的“刚进第二站 predicted 即有 `−0.92` 尖峰”。
3. **S3 首层锁定与逆向平滑（晶圆 13–18）**：一旦到达 `S3L0_s0`，三站 `2.4 m` 磁铁矢高无需依赖单站内部斜率即可唯一锁定曲率，错电荷率在单次更新内从 `38.9%`（或 `13.9%`）骤降至 `0.7%`（`1.1%`），`S3L2_s1` 末端 `filtered = smoothed` 稳定在约 `−11.6%`，并由 `GainMatrixSmoother` 无损回传至全径迹所有 18 个晶圆面（各面 `smoothed` 中位数恒为 `−11.7%`）。

![Eighteen sensors evolution](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/eighteen_sensors_evolution.png)

![Station 2 wafer separation and curvature folding](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/station2_wafer_and_folding.png)

---

### 3. 独立 19 文件（`19,000` 源事件，排除标定文件 `file 0`）C++ 单因素与联合验证

`0.473662 T` 仅由 `file 0 [0, 100)` 的无噪声确定性场传播闭合得出；因此所有涉及 `SeedBFieldTesla = 0.473662` 的评估均严格在其余 **19 个独立文件（`file01..file19`，`19,000` 源事件）** 上统计：

#### A. 各自重建轨迹上的原始 `good-seed` 门限与无门限总体对比（19 个独立文件）

| 实验组（19 独立文件，`19,000` 事件） | writer tracks | GitHub + `purity>0.7` | 同 seed/track 真值轨迹 | 原 `good-seed` 全部 / 同真值轨迹 | `good-seed` 首次 S2 `prt` 中位数 (`<−0.1`) | `good-seed` 首次 S2 `flt` 中位数 (`<−0.1`) | **`good-seed` 首次 S2 `smt` 中位数 (`16% ~ 84%`)** | **`good-seed` `smt < −0.1`** | 无门限同真值 `smt` 中位数 (`<−0.1`) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| `baseline` | 40,094 | 24,119 | 16,269 | 919 / **702** | +4.81% (35.9%) | −74.01% (78.2%) | **−11.67%** (−18.36% ~ −6.62%) | **61.82%** | −2.44% (35.0%) |
| `fix_cov_loc_swap` | 38,480 | 22,643 | **19,586** | 827 / **775** | +7.38% (32.6%) | −77.05% (79.4%) | **−11.56%** (−18.17% ~ −6.43%) | 61.03% | −1.96% (33.6%) |
| `qop_1000_only` | 39,967 | 24,357 | **21,390** | 885 / **833** | +2.85% (**0.6%**) | −16.84% (52.3%) | **−11.20%** (−17.98% ~ −6.22%) | 58.58% | −1.71% (33.0%) |
| `qop_1000_loc_fixed` | 39,351 | 23,185 | **20,839** | 839 / **799** | +2.97% (**0.5%**) | **−5.63%** (48.2%) | **−11.08%** (−17.94% ~ −5.97%) | 57.57% | −1.60% (33.0%) |
| **`seed_b_04737`** | 40,140 | 24,130 | 16,417 | 1,939 / **1,465** | +2.09% (32.0%) | −46.24% (67.3%) | **−0.49%** (**−6.39% ~ +5.27%**) | **4.78%** | −2.42% (35.0%) |
| **`combined_phys_fix`** | 39,565 | 23,197 | **20,921** | 1,747 / **1,666** | **−0.48%** (**0.6%**) | **−1.20%** (39.4%) | **−0.50%** (**−6.36% ~ +5.66%**) | **3.78%** | −1.71% (33.2%) |

#### B. 固定同一组候选路径（`file + entry + seed_index + barcode + ordered hits`）的配对差值与事件聚类 95% 置信区间（19 个独立文件）

固定 **原 `baseline` 已经通过 `good-seed` 的 702 条轨迹**（不重新切换门限选样，检验修改参数本身是否改变同一串 hits 上的拟合结果）：

| 控制组 vs `baseline`（固定原 `good-seed` 候选） | 完全同 hits 配对数（保留率） | 首次 S2 `prt` 逐轨迹差中位数 | 首次 S2 `prt` `Δ(<−0.1)` [95% CI] | 首次 S2 `flt` 逐轨迹差中位数 | **首次 S2 `smt` 逐轨迹差中位数** | **首次 S2 `smt` `Δ(<−0.1)` [95% CI]** |
| --- | ---: | ---: | --- | ---: | ---: | --- |
| `fix_cov_loc_swap` | 628 / 702 (89.46%) | +0.318% | −2.55 pp `[−4.00, −1.10]` | +0.179% | **+0.005%** | +0.16 pp `[−1.13, +1.45]` |
| `qop_1000_only` | 674 / 702 (96.01%) | −2.644% | **−35.01 pp `[−38.60, −31.43]`** | +24.699% | **+0.028%** | −2.08 pp `[−3.15, −1.00]` |
| `qop_1000_loc_fixed` | 607 / 702 (86.47%) | −4.193% | **−34.60 pp `[−38.36, −30.83]`** | +27.481% | **+0.058%** | −2.14 pp `[−3.61, −0.67]` |
| `seed_b_04737` | 692 / 702 (98.58%) | −12.361% | +12.28 pp `[+9.66, +14.90]` | −0.003% | **−0.001%** | **0.00 pp `[−0.57, +0.57]`** |
| `combined_phys_fix` | 603 / 702 (85.90%) | −18.385% | +26.37 pp `[+21.29, +31.45]` | +23.558% | **+0.033%** | −1.99 pp `[−3.43, −0.55]` |

表 A 与表 B 的对照构成严密的双向因果闭合：
- 在表 B 中，`seed_b_04737` 对固定 hits 的 `smoothed` 动量逐轨迹改变量仅为 **`−0.001%`**（在无门限 `16,031 / 16,227 = 98.79%` 同路径轨迹上改变量为 **`0.000%`**），证明 CKF 最终动量完全由三站测量点决定，不受 seed `q/p` 初值牵引；
- 在表 A 中，同一 `seed_b_04737` 在原始未改动的 `good-seed` 门限下却使 `smoothed` 中位数从 **`−11.67%`（`N=702`）变为 `−0.49%`（`N=1,465`）**，证明旧 `−11.7%` 左移纯粹是因为 `B = 0.55 T` 把正常轨迹的 `p_seed` 推高 `+16.12%` 后，`±10%` 真值窗口发生共享测量条件选样偏移。

---

### 4. 无门限总体剩余 `−2.41%` 偏移的电子轫致辐射归因与 Refit 对比（全 20 文件）

在 `RootTrajectoryStatesWriterTool.cxx:426` 中，writer 保存的 `t_eQOP` 来自 `SiHit->particleLink()->momentum()`，即关联 HepMC 粒子的产生顶点动量 `p_prod`（并非各硅片处的瞬时 G4 动量），且 CKF 全程使用 `Acts::ParticleHypothesis::muon()`（仅含电离能损，不含电子轫致辐射模型）。利用全 20 文件提取的 HepMC 同 PDG 延续顶点链（`p_min` 为径迹全长最小延续动量，`p_pre_S2` 为进入 S2 前 `z_F < 1200 mm` 的局部真值动量上界），对无 seed 门限的全部 `17,126` 条同真值轨迹（`17,125` 条有限 `s2_first_smt`）分层统计：

| 无门限同真值电子/正电子子样本（全 20 文件） | 轨迹数 | 相对产生点 `(p_smt − p_prod)/p_prod` 中位数 (`16% ~ 84%`) | 相对产生点 `< −0.1` / `< −0.5` | **相对 S2 前局部真值上界 `(p_smt − p_pre_S2)/p_pre_S2` 中位数 (`16% ~ 84%`)** |
| --- | ---: | --- | ---: | --- |
| 全部无门限同真值轨迹 | 17,125 | **−2.41%** (−34.41% ~ +17.63%) | 34.83% / 10.88% | **+0.68%** (−12.77% ~ +22.29%) |
| 无轫致辐射（`p_min / p_prod ≥ 0.995`） | 6,624 | **+1.65%** (−9.49% ~ +27.63%) | 15.35% / 2.69% | **+1.65%** (−9.47% ~ +27.74%) |
| 轻度辐射（`0.95 ≤ p_min / p_prod < 0.995`） | 2,795 | **+0.51%** (−9.22% ~ +23.63%) | 15.24% / 2.97% | **+1.53%** (−8.20% ~ +24.76%) |
| 硬轫致辐射（`p_min / p_prod < 0.95`） | 7,376 | **−16.99%** (−61.38% ~ +5.92%) | 60.40% / 21.68% | **−0.52%** (−18.90% ~ +16.39%) |
| 其中 S2 前硬辐射（`p_pre_S2 / p_prod < 0.90`） | 3,566 | **−37.11%** (−77.75% ~ −11.06%) | 85.47% / 36.65% | **+0.33%** (−10.92% ~ +14.31%) |

同时对比 CKF 选中轨迹的 `First S2 smoothed` 与后续 `KalmanFitterTool1` 最终输出的 `Trk::Track` 动量：
- `good-seed` 同真值（`N = 725`，refit 成功率 `99.86%`）：CKF 中位数 **`−11.71%`** vs 最终 refit **`−11.39%`**，逐轨迹 `(p_final − p_ckf)/p_ckf` 中位数 **`+0.033%`**（16%–84% 为 `−0.15% ~ +0.63%`）；
- 无门限同真值（`N = 17,126`，refit 成功率 `99.95%`，电荷一致率 `99.90%`）：CKF 中位数 **`−2.41%`** vs 最终 refit **`−2.28%`**，逐轨迹相对差中位数 **`+0.012%`**。

![Smoothed bias causal decomposition](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/smoothed_bias_decomposition.png)

---

### 5. 原宏规范九层残差分布对照（`ckf_residual_perLayer.C`）

以下四图严格采用用户参考脚本 `ckf_residual_perLayer.C`（`chi2Sum/NDF<20 && rMax_smt<90 && NDF>9`，叠加 `purity>0.7` 与原始未改动的 `good-seed` 门限，100 bins，`[−1.2, 1.2]`，raw `Entries`，`P`=绿 / `F`=红 / `S`=蓝，native `t_eQOP` 真值）。注意原宏绘制全部通过 `good-seed` 的轨迹（含跨轨迹串扰候选，如 `baseline` 共 `948` 条轨迹、每层双面约 `1,865–1,887` raw entries；修正协方差后串扰大幅减少，同真值有效轨迹反而上升）：

1. **`baseline`（全 20 文件，`20,000` 事件）**：`S1L2` 已出现红/绿展宽，`S2L0–S2L2` 绿色 `predicted` 与红色 `filtered` 呈强 `−0.92` 左峰，蓝色 `smoothed` 整体左偏在 `−12%`；
2. **`qop_1000_loc_fixed`（全 20 文件，`20,000` 事件）**：`S1L0–S1L2` 红/绿完全收窄在真值附近，`S2` 的 `−0.92` 假峰大幅压低，但因 `SeedBFieldTesla` 仍为 `0.55 T`，蓝色 `smoothed` 依然停留在 `−11%`；
3. **`seed_b_04737`（独立 19 文件，`19,000` 事件，排除 `file 0`）**：仍保留 `100,000` 曲率膨胀，故 `S2` 红/绿仍有强左峰，但蓝色 `smoothed` 在原始 `good-seed` 门限下已居中归零（`−0.49%`），直接直观展示 **`S2 predicted` 左峰存在时 `smoothed` 依然无偏**；
4. **`combined_phys_fix`（独立 19 文件，`19,000` 事件，排除 `file 0`）**：`S1` 发散消除、`S2` 尖峰压低、蓝色 `smoothed` 居中归零且有效轨迹产额翻倍（raw entries 达 `~3,470`/层）。

![Baseline nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/110051_baseline_nine_layers.png)

![Fix covLoc swap and qop 1000 nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/110051_qop_1000_loc_fixed_nine_layers.png)

![Seed B 0.4737 T held-out 19 files nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/110051_seed_b_04737_heldout19_nine_layers.png)

![Combined physics fix held-out 19 files nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/110051_combined_phys_fix_heldout19_nine_layers.png)

其余单因素九层图见 [110051_fix_cov_loc_swap_nine_layers.pdf](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/110051_fix_cov_loc_swap_nine_layers.pdf) 与 [110051_qop_1000_only_nine_layers.pdf](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/110051_qop_1000_only_nine_layers.pdf)；全部 10 页矢量图汇编见 [report.pdf](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/report.pdf)。

---

## 执行记录、验证与产物清单

- **HTCondor 调度选择与作业完成记录**：
  - 提交前查询 `condor_status -schedd` 与 `condor_q -totals`，排除曾记录 EOS 兼容异常的 `bigbird21.cern.ch`，选择已验证兼容且普通 `bigbird` 中 `running+idle` 最低的 `bigbird20.cern.ch`：
    1. 全量重建矩阵 `latest_deep_campaign_v1_20261118`（cluster **`1584662`**，6 组 × 20 文件 = `120` 作业，`1 CPU / 3800 MB / 8 GB scratch / workday`）：**`120/120` 全部 `exit_code == 0`**，提交记录见 [submission.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/latest_deep_campaign_v1_20261118_worker/submission.json)；
    2. 全量归约与出图初版 `latest_deep_full_study_v1_20261118`（cluster **`1584675`**，`exit_code == 0`）：目视核查发现 `eighteen_sensors_evolution` 中间子图旧标题文字与图例高度略近、两张九层图图例标签略长，随即在 `analysis/latest_deep_study.py` 修正；
    3. 全量归约与出图定稿版 `latest_deep_full_study_v2_20261118`（cluster **`1584685`**，`1 CPU / 5500 MB / 4 GB / microcentury`，`exit_code == 0`）：提交记录见 [submission.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/latest_deep_full_study_v2_20261118_worker/submission.json)。
  - 大型中间 ROOT/xAOD 全部在 worker 本地 scratch 内完成匹配与直方图归约并当场清理，EOS 仅保存紧凑 NPZ 与 27 个原生 `TH1` 计数。
- **图表逐 bin 闭合与回归单测**：
  - 6 组九层图共 `162` 个 `TH1`、含 underflow/overflow 共 `16,524` 个分箱与 worker 原宏输出 100% 精确一致（`figure_validation.json`），10 页 [report.pdf](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/report.pdf) 全部通过目视检查；
  - 仓库单元测试 `python3 -m unittest discover -s tests -v`（新增 [test_latest_deep_momentum.py](../tests/test_latest_deep_momentum.py) 覆盖旋转 seed 平面 `loc0=y_F / loc1=x_F` 映射、曲率折叠、`>6` 共享 cluster 规则、弦线 vs 圆切线 `φ` 偏置及独立文件排除 `file 0` 检查）共 **108 项全部通过**。
- **代码与关键产物路径**：
  - 仓库脚本与配置：[latest_deep_tracking.py](../configs/latest_deep_tracking.py)、[prepare_latest_deep.py](../scripts/prepare_latest_deep.py)、[latest_deep_tasks.py](../scripts/latest_deep_tasks.py)、[run_latest_deep.py](../scripts/run_latest_deep.py)、[latest_deep_reduce.py](../analysis/latest_deep_reduce.py)、[latest_deep_study.py](../analysis/latest_deep_study.py)；
  - 数据与图索引：[summary.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/summary.json)、[figure_validation.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/figure_validation.json)、[bin_counts.npz](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/bin_counts.npz)、[reference_bin_counts.npz](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/reference_bin_counts.npz)、[histogram_metadata.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2/histogram_metadata.json)。

复现命令（通过 HTCondor `scripts/submit_momentum.py` 执行）：

```bash
python3 analysis/latest_deep_study.py \
  --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1 \
  --chains-dir /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/truth_vertices \
  --output /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261118_v1/full_study_v2
```

---

## 事实与假设区分

- **已证实事实**：
  1. `geoRotation_calypso_latest` 正常加载 `material-maps-alma9-rotated-rowmajor.json.gz`，native 查表闭合，材料不是本轮 `S2 predicted` 或 `smoothed` 左偏的原因；
  2. 全量 seed 残差的 `−1.0` 尖峰 `92.50%` 来自双径迹跨组合（`59.76%`）与次级/非纯段（`32.74%`），由站间笛卡尔积、3 站 0 自由度圆拟合、去重 `> 6` 放行整站 6-cluster 共享及 `getY()` 负根号 `NaN` 排序共同造成；干净低辐射 3 站 seed 因 `B = 0.55 T` 硬编码高估约 `+16%`；
  3. `CKF2Config.py` 存在旋转 seed 平面 `covLoc0/covLoc1` 轴互换（`y_F` 赋 `81 mm²`、`x_F` 赋 `0.0361 mm²`）及 `covQOverP` 100 倍额外膨胀（`100,000` vs `1,000`）；后者导致 `S1L2` 短基线发散并在 `S2L0_s0` 更新后结合双晶圆合并与 `1/|q/p|` 折叠形成 `−0.92` 左峰；
  4. `S2 predicted` 左偏与 `smoothed` 左移无因果传递；`good-seed` 下 `−11.7%` 的 `smoothed` 偏移由 `B = 0.55 T` 引起的 `−9.2%` 共享测量选样偏移与电子真实轫致辐射相对产生点真值的 `−2.4%` 物理能损共同组成；无辐射电子（`p_min/p_prod ≥ 0.995`）相对产生点的中位数为 `+1.65%`。
- **待进一步验证的边界（不作过度外推）**：
  1. 本次 RDO 未保存每个硅片处的瞬时 Geant4 step 动量；HepMC 延续顶点链给出严格的辐射位置区间上下界，精确单硅片局部残差与电子 GSF（Gaussian Sum Filter）拟合收益需在带 step-truth 的模拟样本中量化；
  2. `0.473662 T` 是 3-station forward 几何下的有效圆拟合常数，若推广到 4-station（含 Station 0）或 backward 配置，应改为直接由真实场积分计算的 propagation-aware seed 估计器而非单一全局常数。

---

## 下一步计划

1. **生产级 Seed 重构与全拓扑验收（S07）**：
   - **具体任务**：在隔离 Calypso 分支中：（a）修正 `CKF2Config.py` 的旋转平面 `covLoc0/covLoc1` 对应关系；（b）在 `CircleFitTrackSeedTool` 中保护 `getY()` 负判别式、将去重阈值收紧以避免无约束整站共享、用圆切线 + 3 站 `x(z)` 拟合替代跨磁铁弦线方向，并用场积分替代硬编码 `0.55 T`；
   - **所需输入**：DSID `110051`（高动量 `e⁺e⁻`）、`110033`（低动量/不同寿命 `A'`）及 `100116/100117`（单 `μ±`，无轫致辐射且无双径迹组合本底）；
   - **验收标准**：在**不加任何 truth seed 门限**的独立样本上，同时报告双径迹重建效率、假组合率、跨轨迹串扰率、各层测量 pull RMS、电荷误判率与最终 `Trk::Track` 动量分辨率，要求干净单粒子动量中心偏差 `< 1%` 且双轻子配对效率不低于当前基线。
