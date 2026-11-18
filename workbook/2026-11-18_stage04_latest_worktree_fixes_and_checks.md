# S04：latest Calypso 新建修复 worktree 五项逐级验证、粒子假设对照与 Seed 组合/拟合检查（五类样本 600k 事件）

## 阶段性结论与核心实测结果

基于干净的 `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest`（commit `7fd40df8cd3b60e81b4095a4f84310c4a54700a0`，包内旋转材料图 `material-maps-alma9-rotated-rowmajor.json.gz`）新建独立 Git worktree `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1`（分支 `fix/latest-seed-and-momentum-v1`，commit `e47b85fd74ecb785dcc0ece6b8f61022075b2e7d`），默认启用五项核心源码修复（`C1` 修正旋转平面 `covLoc0/covLoc1` 轴互换、`C2` 恢复原版方差膨胀因子 `[100, 100, 100, 100, 1000]`、`C3` 标定 `SeedBFieldTesla = 0.473662 T`、`C4` 修复 `getY()` 负根号 `NaN` 与严格弱序排序、`C5` 将 seed 初始方向从跨磁铁 `S1→S2` 弦线改为局部圆切线 + 三站 `x_F(z_F)` 线性拟合），并保留两类受控检查开关（`ParticleHypothesis::muon()` vs `electron()`、以及 18-cluster 双磁铁积分场拟合与去重/预筛选检查）。在 **5 类物理样本**——DSID `110051`（`A' → e⁺e⁻`, `ε = 10⁻⁴`, 20 文件 `20,000` 事件，排除标定文件后独立 19 文件 `19,000` 事件）、`110033`（`A' → e⁺e⁻`, `ε = 10⁻⁵`, 5 文件 `10,000` 事件）、`100116`（`PG μ⁺` 2D FLUKA, `10,000` 事件）、`100117`（`PG μ⁻` 2D FLUKA, `10,000` 事件）、`100130`（`PG K_S` FASERnu 末端强子簇射与衰变 `π±/K±/p`, `10,000` 事件）——上完成 **10 组配置 × 5 类样本 = `330` 个双核 HTCondor 作业、共 `600,000` 重建事件**（`330/330` 全部 `exit_code == 0`），得到以下实测结论：

1. **`SeedBFieldTesla = 0.473662 T` 的跨粒子种类与全磁铁孔径普适性（回答常数是否适用于 muon 及其他样本）**：
   - 在真实 Calypso `FaserFieldTable_v2.root` 三维磁场与完整 `CircleFit` 确定性传播闭合（5 类样本各 `[0, 250)` 事件，覆盖 `p` 中位数 `32.3–920.6 GeV`、入射斜率 `90%` 分位数 `0.36–26.9 mrad`、横向孔径 `r_max < 100 mm`）中，5 类样本的圆拟合等效磁场 `B_eff` 中位数分别为 **`0.473661 T`（`110051` `e±`）**、**`0.473663 T`（`110033` `e±`）**、**`0.473673 T`（`100116` `μ⁺`）**、**`0.473670 T`（`100117` `μ⁻`）**、**`0.473717 T`（`100130` `K_S` 强子）**，跨粒子种类、正负电荷与 75 倍角度跨度的最大中位数差异仅 **`0.000056 T`（`0.012%`）**；采用 `0.473662 T` 时确定性闭合动量残差中位数全部落在 **`−0.012%` 至 `+0.0002%`**（原 `0.55 T` 在 5 类样本上一致高估 **`+16.10%` 至 `+16.12%`**）。
   - 在全量 RDO 重建的干净 3-station seed（相对进入 Station 1 前局部真值 `p_st1`）中，将 `0.55 T`（`C2`）改为 `0.473662 T`（`C3`）使残差中位数在全部 5 类样本上同步归零：`110051` 独立 19 文件 **`+16.35% → +0.20%`**、`110033` **`+15.39% → −0.63%`**、`100116` `μ⁺` **`+13.99% → −1.83%`**、`100117` `μ⁻` **`+14.06% → −1.77%`**、`100130` `K_S` 强子 **`+15.22% → −0.77%`**（其中 `100116/100117` 的 `~ −1.8%` 对应 `z_F = −3990 mm` 产生后穿越 FASERnu 钨靶末段及 `15 cm` 铅屏蔽时、在最后一次离散 HepMC 顶点之后约 `0.8–1.1 GeV` 的连续电离 `dE/dx` 真实能损）。因此 **`0.473662 T` 由三站几何位置与双偶极铁积分场唯一决定，严格适用于电子、正负 μ 子与带电强子全部样本**。
2. **五项核心修复逐级对照（`C0 → C1 → C2 → C3 → C4 → C5`）证明每一步均不可省略**：
   - **`C0 → C1`（修正旋转 seed 平面 `covLoc0 = 0.0361 mm², covLoc1 = 81.0 mm²`）**：消除弯转精细轴 `y_F` 被错赋 `9 mm` 先验导致的跨径迹抢点，CKF 跨粒子串扰率在 `110051`（独立 19f）从 **`32.54%` 骤降至 `13.53%`**，在 `110033` 从 **`13.76%` 降至 `5.41%`**，在 `100130` 强子样本从 **`2.42%` 降至 `1.26%`**。但单独启用 `C1` 时因 `[10000, 100000]` 超宽角度/曲率方差仍在，`S1L2` 发散在紧缩 `y_F` 门限下反而使双径迹事件效率从 `74.90%` 暂降至 `71.57%`，证明 `C1` 必须与 `C2` 配合。
   - **`C1 → C2`（恢复原版方差膨胀因子 `[100, 100, 100, 100, 1000]`）**：抑制 `S1L2` 短臂长曲率发散与 `S2L0` 卡尔曼增益爆炸，**Station 2 首块晶圆测量新息 pull RMS 从失真低估恢复到物理单位高斯附近**（`110051`: `0.283 → 1.076`；`100116` `μ⁺`: `0.207 → 0.918`；`100130` 强子: `0.353 → 0.878`），`110051` 同真值轨迹 `S2` 首片 `filtered` `< −0.5` 负尾从 **`62.91%` 压低至 `34.79%`**（中位数 `−71.05% → −27.18%`），串扰率进一步降至 **`9.37%`**（`110033`: `2.52%`，`100130`: `0.22%`），且 fiducial 重建效率跃升至 **`79.05%`**（超越 `C0` 的 `74.90%`，双径迹配对效率 **`61.58% → 67.52%`**）。
   - **`C2 → C3`（标定 `SeedBFieldTesla = 0.473662 T`）**：消除全部 5 类样本 seed 动量的 `+16.1%` 标度偏置与九层图 `S1 predicted/filtered` 的 `+0.16` 右偏峰；在未改动的原始 `good-seed` 门限下，通过门限的有效同真值轨迹数大幅增加且 `S2 smoothed` 左移彻底消失（`110051` 独立 19f: `N = 917 → 1,893`，中位数 **`−11.15% → −0.51%`**；`110033`: `N = 376 → 1,390`，**`−9.69% → −0.79%`**；`100130`: `N = 32 → 141`，**`−2.74% → −0.58%`**）。
   - **`C3 → C4`（修复 `getY()` 负根号判别式 `NaN` 与 `std::list::sort` 严格弱序崩溃风险）**：在全部 5 类样本中将原始与最终选中 seed 的 `NaN` `getY()/chi2` 数量**从数千条全部清零**（`110051` 20f 原始 `1,971 → 0`、选中 `750 → 0`；`110033`: `1,804 → 0`；`100116` `μ⁺`: **`3,088 → 0`**；`100117` `μ⁻`: **`3,042 → 0`**；`100130` `K_S`: **`7,158 → 0`**），消除 C++ 未定义排序行为隐患而不损失正常轨迹效率。
   - **`C4 → C5`（将 `S1→S2` 弦线改为 seed 平面局部圆切线 + 三站 `x_F(z_F)` 线性拟合）**：彻底消除跨 Dipole 1 弦线造成的正负电荷反号弯转角偏置（`110051` 独立 19f 干净 seed 的 `Δφ_A` 中位数由 `e⁻: +0.114 mrad / e⁺: −0.113 mrad` 归零至 **`+0.001 / +0.007 mrad`**），并在低动量大偏转样本上将 seed 弯转角分辨率 `σ_68(Δφ_A)` 提升 **`5.4–19.6` 倍**（`110051`: `0.185 → 0.034 mrad`；`110033`: **`1.546 → 0.093 mrad`**；`100116` `μ⁺`: **`1.357 → 0.112 mrad`**；`100117` `μ⁻`: **`1.400 → 0.105 mrad`**；`100130` 强子: **`3.485 → 0.178 mrad`**），同时将非弯转角分辨率 `σ_68(Δθ_A)` 在全部 5 类样本上一致提升约 **`1.7–1.9` 倍**（`~0.62–0.67 mrad → ~0.32–0.39 mrad`），使原始 `good-seed` 有效轨迹产额再次跃升 **`+45%` 至 `+110%`**（`110051` 独立 19f: `1,894 → 2,740`，累计较 `C0` 增长 **`3.85` 倍**；`100116` `μ⁺`: `439 → 917`；`100117` `μ⁻`: `475 → 934`），并将 `110051` 串扰率进一步降至 **`8.01%`**、双径迹配对效率升至 **`67.59%`**。
3. **检查一结论（相同事件、相同物质效应下 `ParticleHypothesis::muon()` [`C5`] vs `electron()` [`H1`]）**：
   - **动量偏差、动量分辨率与重建效率几乎不变**：在全部 5 类样本中，最终 refit 动量残差中位数变化 **`≤ 0.02%`**（`110051`: `−0.00%` vs `−0.01%`，硬轫致辐射子集 `p_min/p_st1 < 0.95` 为 `−9.12%` vs `−9.16%`；`110033`: `−0.91%` vs `−0.92%`；`100116` `μ⁺`: `−1.69%` vs `−1.71%`），68% 半宽分辨率变化 **`≤ 0.03%`**，fiducial 轨迹效率变化 **`≤ 0.06%`**，post-refit `χ²/ndf` 中位数变化 **`≤ 0.02`**。根本原因是 Calypso 的 CKF 与 `KalmanFitter` 均为基于 `Acts::MaterialInteractor` 的**单高斯卡尔曼滤波器**而非多分量 Bethe-Heitler 高斯和滤波器（`Acts::GaussianSumFitter`, GSF）：切为 `electron()` 仅在每层硅片加入确定性平均辐射能损与单高斯方差项，无法拟合单次大角度/大能损硬轫致辐射折线。
   - **协方差 Pull 表现出显著粒子种类分化**：在电子样本上，`electron()` 补充的辐射过程噪声使首块晶圆 `w0` 的 `q/p` pull `σ_68` 更贴近单位高斯 `1.0`（`110033` 全体 `1.94 → 1.68`、无硬辐射子集 **`1.24 → 1.08`**；`110051` 全体 `1.85 → 1.78`、无硬辐射子集 `1.42 → 1.38`）；但在**非电子样本（`100116/100117` μ 子与 `100130` 强子）上，强行使用 `electron()` 会注入非物理的韧致辐射过程方差**，使 `μ±` 的 `φ` pull `σ_68` 被过度压缩至 `0.99 < 1.0` 并虚增约 `25%` 曲率误差棒。因此**不应在通用重建中全局替换为 `electron()`**；电子专用通道可仅作协方差标定选配，而真正改善电子硬轫致辐射动量尾需引入 `Acts::GaussianSumFitter`。
4. **检查二结论（3 站圆拟合判别力不足 vs 18-cluster 双磁铁积分场拟合 [`S3`] 与笛卡尔积预筛选/强化去重 [`S1, S2`]）**：
   - **纯算法替换 `S3_strip_only`（用 18-cluster 5 参数双磁铁真实积分场最小二乘拟合 `ndf ≈ 13` 替代 3 点圆拟合 `ndf = 1` 进行排序与初值估计，不删减组合）**：在全部 5 类样本上保持 fiducial 效率不变（`110033`: `86.85%` vs `86.87%`；`100116`: `99.50%`；`100117`: `99.53%`；`100130`: `23.38%`，配对效率 `6.65%` vs `6.59%`；`110051`: `78.52%`），同时取得 **10 组配置中最低的 CKF 跨径迹串扰率**（`110051`: **`7.51%`**；`110033`: **`2.17%`**）、**最高的 `good-seed` 产额**（`110051`: **`2,792` 条**，中位数 **`−0.32%`**；`110033`: **`2,167` 条**）及略优的最终动量分辨率（`110051` `σ_68`: `21.43%` vs `21.52%`；`100130`: `6.31%` vs `6.48%`）。但只要沿用逐对 `> 6` cluster 去重且不设 `χ²` 上限，排在真实轨迹之后的残余拉丁方跨径迹假组合仍会被继续选出（`110051` 选中 seed 仍为 `50,108` 条）。
   - **叠加预筛选与强化去重 `S1_dedup_prefilter` / `S2_strip_field_fit`（`x_F` 预筛选 + 已选并集 `> 6` 与 `S2/S3 ≥ 4` 共享否决 + `χ²/ndf ≤ 35`）**：在单径迹 `μ±` 样本上**零效率损失**（`99.53% / 99.58%`），在多径迹样本上大幅压减选中 seed 总数（`110051`: `−31.9%`；`110033`: `−48.4%`；`100130`: `−50.9%`，跨正负强子假组合从 `10.8%` 降至 **`0.3%`**），使 `110051` 全 seed 的 `< −0.5` 假峰占比从 `36.5%` 降至 **`21.6%–22.8%`**，并压低最终 GitHub 筛选后的假轨迹比例（`24.48% → 22.57%`）；但在双电子样本上因硬轫致辐射折线轨迹超出 `χ²/ndf ≤ 35` 以及 `> 1 TeV` 极小张角电子对在 Station 2 仍合并共享 cluster，造成 **`−2.6 pp`（`110051`）至 `−5.5 pp`（`110033`）** 的电子效率折损。因此生产上宜保留宽松前端 seed 候选池配合 `18-cluster` 积分场排序（`S3`），而将严格共享 cluster 去重放在 CKF 之后的轨迹级歧义消解（`AmbiguityResolution`）阶段执行。

---

## 阶段目标

在保持 `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest` 原始目录完全不动的前提下完成四项工作：

1. 基于 `geoRotation_calypso_latest` 创建新 worktree，修复前序源码审计定位的 5 个问题（旋转 seed 局部协方差轴颠倒、初始方差膨胀过大、seed 磁场常数硬编码 `0.55 T`、`getY()` 负根号 `NaN`、seed 初始方向采用 `S1→S2` 弦线），检验 `0.473662 T` 在 `μ±` 与强子等样本上的普适性，并通过逐级累加实验（`C0 → C1 → C2 → C3 → C4 → C5`）量化每一项修复对动量及其他动力学量分布的独立贡献；
2. 在完全相同事件与相同旋转材料图设置下对比 `ParticleHypothesis::muon()` 与 `ParticleHypothesis::electron()` 的动量偏置、动量分辨率、找迹效率、`χ²/ndf` 及首片晶圆 5 参数协方差 pull；
3. 对 seed 站间笛卡尔积跨径迹假组合、去重规则及 3 站圆拟合判别力不足做改进探索与对照检查（区分纯 18-cluster 双磁铁积分场拟合 `S3` 与叠加预筛选/强化去重 `S1, S2`）；
4. 将全部对比推广至暗光子（硬谱 `110051`、软谱 `110033`）、正负 μ 子（`100116`、`100117`）及 `K_S` 强子（`100130`）共 5 类 MC24 样本，全程采用 C++ 二进制零 JSON 审计、本地 scratch 多进程并发重建与多线程归约加速处理。

---

## 工作内容与工程实现

### 1. 新建 Worktree、源码修改与编译验证

结果根目录为 `/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1`（下文简称 `R3`）。

| 项目 | 实际配置与路径 |
| --- | --- |
| 基准源码（保持 clean） | `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest`（HEAD `7fd40df8cd3b60e81b4095a4f84310c4a54700a0`） |
| 新建修复 worktree | `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1`（分支 `fix/latest-seed-and-momentum-v1`，commit `e47b85fd74ecb785dcc0ece6b8f61022075b2e7d`） |
| 仓库可复用补丁与组件 | [0050_latest_seed_fixes_and_checks.patch](../patches/0050_latest_seed_fixes_and_checks.patch)、[LatestFastValidationAudit.h](../patches/LatestFastValidationAudit.h)、[MultiSampleSeedFieldClosure.cxx](../patches/MultiSampleSeedFieldClosure.cxx) |
| 编译库产物 | [libFaserActsKalmanFilter.so](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/software_fix_v1/lib/libFaserActsKalmanFilter.so)、[libMultiSampleSeedFieldClosure.so](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/software_field_closure/libMultiSampleSeedFieldClosure.so)（Athena `24.0.41`, `x86_64-el9-gcc13-opt`, 14 核并行编译） |
| 旋转材料图 | [material-latest-rowmajor.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/material-latest-rowmajor.json)（原样解压自 `geoRotation_calypso_latest` 包内 `material-maps-alma9-rotated-rowmajor.json.gz`，`11,866,364` bytes） |
| 坐标系与参考面约定 | 原始 FASER 系 `F: (x_F, y_F, z_F)`；旋转 ACTS 系 `A: (x_A, y_A, z_A) = (z_F, y_F, −x_F)`；seed 参考面为 `targetSurface`（`PlaneSurface` at `x_A = z_F = −10 mm`, `eBoundLoc0 = +y_F`, `eBoundLoc1 = +x_F`）；首测量参考面 `w0` 为轨迹在 Station 1 的第一块有源 SCT 晶圆 `Acts::Surface` |

**五项默认修复与两项可选检查的具体源码改动**（均位于 `Tracking/Acts/FaserActsKalmanFilter`）：
1. **修复 1 + 修复 2**（`python/CKF2Config.py`）：修正旋转平面局部测量标准差顺序为 `[0.019, 0.9, 0.033, 0.00029, 0.0001]`（使 `covLoc0` 对应精细弯转轴 `y_F`、`covLoc1` 对应粗立体角轴 `x_F`），并将默认 `initialVarianceInflation` 恢复为原版 `[100, 100, 100, 100, 1000]`（即 `covLoc0 = 0.0361 mm², covLoc1 = 81.0 mm², covPhi = 0.1089, covTheta = 8.41×10⁻⁶, covQOverP = 1.0×10⁻⁵ GeV⁻²`）；保留 `legacySwappedCovLoc` 与 `initialVarianceInflation` 参数用于 `C0–C2` 回退对照。
2. **修复 3**（`CircleFitTrackSeedTool.{h,cxx}`）：将默认 `SeedBFieldTesla` 设为 `0.473662`（支持运行期传入 `0.55` 作 `C0–C2` 对照）。
3. **修复 4**（`CircleFitTrackSeedTool.{h,cxx}`）：`FixGetYNaN = true` 时在 `getY(z)` 中对判别式 `rad = r² − (z − cx)² ≤ 0` 保护平方根并在 `getChi2()` 中施加圆外几何距离惩罚 `((|z − cx| − r)/σ)²`，同时在 `std::list::sort` 比较器中显式将非有限分数排至末尾以保证严格弱序。
4. **修复 5**（`CircleFitTrackSeedTool.{h,cxx}`）：`UseLocalTangent = true` 时在 `get_params()` 中以 seed 参考面 `z_F = origin` 处的局部圆切线斜率 `dy_F/dz_F = −(z_eval − cx)/(y_eval − cy)`（当 `|z_F − cx| ≥ r` 时平滑夹取至 `99.9% r`）结合全部站中心的 `LinearFit` 斜率 `dx_F/dz_F = c1` 构造局部切向单位矢量，替代原跨磁铁 `S1→S2` 弦线。
5. **可选检查 1（默认关闭，`AbsPdgHypothesis = 13`）**：在 `CircleFitTrackSeedTool`、`CKF2`、`TrackFinderFunction`、`CreateTrkTrackTool` 与 `KalmanFitterTool` 全链路贯通 `AbsPdgHypothesis`（`13` → `ParticleHypothesis::muon()`，`11` → `electron()`，`211` → `pion()`），并在 `KalmanFitterTool` 中同时保存 `targetSurface` 与首片晶圆 `w0` 的拟合参数及协方差对角元，实现与 CKF 和局部 `FaserSiHit` 真值的严格同面（`w0`）pull 对比。同时修正 `RootSeedWriteTool.cxx` 原 `FIXME` 中对 `π±/K±/p` 正负强子电荷符号取反的缺陷。
6. **可选检查 2（默认关闭）**：在 `CircleFitTrackSeedTool` 中实现基于真实 `FASERMagneticFieldWrapper` 轴向积分场表 `(J_y, I_y, J_z, I_z)` 的 18-cluster 一维条测量 5 参数 `LDLT` 线性最小二乘拟合（`SeedStripFieldFit`，单 seed 耗时 `< 1.5 µs`，`ndf ≈ 13`）、站间粗 `x_F` 预筛选（`SeedPreFilter`）及防 `S2/S3` 复用的并集去重（`SeedImprovedDedup`）。
7. **计算加速与逐位回归闭合**：实现 C++ 原生二进制审计头 `LatestFastValidationAudit.h`（直接输出紧凑 `float64` 矩阵，省去慢速 `xAODTruthCnvAlg` 与逐状态 JSON 序列化）；在 `C0_baseline` 参数下与未打补丁的 `geoRotation_calypso_latest/run/lib`（`native_smoke`）在 `110051` 前 20 事件上比对 `track_states.root` 的 `chi2Sum`、`NDF` 与全部轨迹状态 `eQOP_prt/flt/smt`，验证 **100% 逐位相等**。单作业双核（`ProcessPoolExecutor` + 本地 `$TMPDIR`）重建并归约 `5,000` 个 muon 事件总耗时仅 **`36.3 s`（`138 events/s/job`）**。

### 2. 五类冻结样本与十组实验配置

| 样本 DSID | 物理过程与产生位置 | 冻结文件数与事件范围 | 源事件数 | 物理特征与检验目的 |
| --- | --- | --- | ---: | --- |
| `110051`（及 `heldout19`） | `A' → e⁺e⁻`（`m = 10 MeV, ε = 10⁻⁴`，衰变体真空产生） | `master_dev` 全 20 文件 `[0,1000)`；标定评估严格排除 `file00` 用后 19 文件 | `20,000`（独立 `19,000`） | `p` 中位数 `921 GeV`、极小张角（`0.36 mrad`）共线双电子 + 强轫致辐射 + 高跨径迹组合本底 |
| `110033` | `A' → e⁺e⁻`（`m = 10 MeV, ε = 10⁻⁵`，衰变体真空产生） | `master_dev` 前 5 文件 `00001..00005` 各 `[0,2000)` | `10,000` | `p` 中位数 `91 GeV`、较大张角（`1.80 mrad`）软谱双电子，检验低动量弦线偏置与电子假设 |
| `100116` | `PG μ⁺` 2D FLUKA（`z_F = −3990 mm`，FASERnu 钨靶上游产生） | `s0012` 文件 `00000-00009` `[0,10000)` | `10,000` | `p` 中位数 `67 GeV`、大斜率（`13.1 mrad`）、满孔径单正 μ 子，无双径迹本底且无电子硬轫致辐射 |
| `100117` | `PG μ⁻` 2D FLUKA（`z_F = −3990 mm`，FASERnu 钨靶上游产生） | `s0012` 文件 `00000-00009` `[0,10000)` | `10,000` | `p` 中位数 `85 GeV`、大斜率（`12.3 mrad`）、满孔径单负 μ 子，与 `100116` 检验正负电荷对称性 |
| `100130` | `PG K_S`（`s0013`, `z_F ≈ −1960 mm` FASERnu 钨靶末端注入） | `s0013` 文件 `00000-00009` `[0,10000)` | `10,000` | `p` 中位数 `32 GeV`、极大斜率（`26.9 mrad`），含 `K_S → π⁺π⁻` 衰变与钨靶强子簇射（`π±, K±, p`）高多重数环境 |

十组受控配置定义如下（每组均重建上述 `60,000` 事件，独立评估计 `59,000` 事件）：
- **`C0_baseline`**：`geoRotation_calypso_latest` 原生设置（反序 `covLoc`、膨胀 `[100,100,100,10000,100000]`、`B = 0.55 T`、无 `NaN` 保护、`S1→S2` 弦线、`muon()` 假设）；
- **`C1_cov_loc`**：`C0` + 修正旋转平面 `covLoc0/covLoc1` 轴对应；
- **`C2_orig_var`**：`C1` + 恢复原版方差膨胀 `[100, 100, 100, 100, 1000]`；
- **`C3_seed_b`**：`C2` + 标定 `SeedBFieldTesla = 0.473662 T`；
- **`C4_no_nan`**：`C3` + `getY()` 负根号保护与严格弱序排序；
- **`C5_local_tan`（新建 worktree 默认配置）**：`C4` + 局部圆切线与三站线性 `x_F(z_F)` 初始方向；
- **`H1_electron_hypo`（检查 1）**：`C5` 基础上将全链路粒子假设改为 `ParticleHypothesis::electron()`（旋转材料图、散射与能损开关保持完全一致）；
- **`S3_strip_only`（检查 2a）**：`C5` 基础上仅将 3 点圆拟合替换为 18-cluster 双磁铁积分场拟合用于候选排序与 seed 初值，保留全部笛卡尔积组合与原 `> 6` 去重；
- **`S1_dedup_prefilter`（检查 2b）**：`C5`（3 点圆拟合）基础上启用站间 `x_F` 预筛选 + 已选并集 `> 6` 与 `S2/S3 ≥ 4` 共享否决 + `χ²/ndf ≤ 35`；
- **`S2_strip_field_fit`（检查 2c）**：同时启用 18-cluster 双磁铁积分场拟合 + `S1` 的预筛选与强化去重。

---

## 详细实测数据与物理分析

### 1. `SeedBFieldTesla = 0.473662 T` 的五类样本普适性验证

在无散射/测量噪声干扰的真实三维场 `FaserActsExtrapolationTool` 外推 + 原生 `CircleFit` 闭合测试（每样本 `[0, 250)` 事件，要求三站均在磁铁孔径 `r < 100 mm` 内）与全量 RDO 重建干净 3-station seed（单粒子同 barcode、纯度 `≥ 0.7`、`p_pre_S2 / p_st1 ≥ 0.95`）上对比 `0.55 T` 与 `0.473662 T`：

| 样本 DSID（粒子） | 闭合轨迹数 | `p` 中位数 [GeV] | `r_max` 90% 分位 [mm] | 斜率 90% 分位 [mrad] | 确定性 `B_eff` 中位数 `[p16, p84]` [T] | 闭合残差中位数（`0.55 T` → `0.473662 T`） | 全量 RDO 干净 3 站 seed `(p_seed − p_st1)/p_st1` 中位数（`C2: 0.55 T` → `C3: 0.4737 T`） |
| --- | ---: | ---: | ---: | ---: | --- | --- | --- |
| `110051`（`A' → e⁺e⁻`；RDO 取独立 19f） | 492 | 920.6 | 93.4 | 0.36 | `0.473661 [0.473337, 0.474472]` | `+16.117% → +0.0002%` | **`+16.35% → +0.20%`**（`σ_68`: `17.13% → 14.76%`） |
| `110033`（`A' → e⁺e⁻`, 10k） | 442 | 91.0 | 93.5 | 1.80 | `0.473663 [0.472919, 0.474042]` | `+16.116% → −0.0003%` | **`+15.39% → −0.63%`**（`σ_68`: `4.40% → 3.79%`） |
| `100116`（`PG μ⁺` 2D, 10k） | 165 | 66.6 | 92.6 | 13.13 | `0.473673 [0.473600, 0.474006]` | `+16.114% → −0.0023%` | **`+13.99% → −1.83%`**（`σ_68`: `6.39% → 5.50%`） |
| `100117`（`PG μ⁻` 2D, 10k） | 153 | 84.5 | 92.5 | 12.34 | `0.473670 [0.472639, 0.473801]` | `+16.115% → −0.0016%` | **`+14.06% → −1.77%`**（`σ_68`: `6.41% → 5.52%`） |
| `100130`（`PG K_S` 强子, 10k） | 370 | 32.3 | 94.5 | 26.86 | `0.473717 [0.473609, 0.474089]` | `+16.103% → −0.0116%` | **`+15.22% → −0.77%`**（`σ_68`: `2.46% → 2.12%`） |

**关于 `100116/100117` μ 子上游钨/铅能损的物理核实**：`100116` 与 `100117` 在 FASERnu 钨靶上游 `z_F = −3990 mm` 产生，穿越 `1.1 m` 钨靶与 `15 cm` 铅屏蔽后才进入 Station 1。实测最终重建动量相对产生点 `p_prod` 的绝对差 `med(p_rec − p_prod)` 在 `[15, 40)`, `[40, 80)`, `[80, 160)`, `[160, 400) GeV` 分箱中分别为 `−3.64, −4.14, −5.04, −6.86 GeV`（`−13.3%` 至 `−3.3%`）；追踪 `TruthEvent` 离散延续顶点链得到 `p_st1`（`91.5%` 的 μ 子含离散顶点）后，残余 `med(p_rec − p_st1)` 在各能区收敛至恒定的 **`−0.80` 至 `−1.32 GeV`**（`[160, 400) GeV` 相对偏差仅 **`−0.35%`**），精确对应最后一次离散顶点之后剩余钨板与铅屏蔽中的连续 Bethe-Bloch 电离能损（Geant4 连续电离不额外生成 HepMC 顶点）。

![Seed B_eff universality across 5 species](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/seed_b_multisample_closure.png)

---

### 2. 五项修复逐级消融（`C0 → C1 → C2 → C3 → C4 → C5`）：动量与动力学量演化

#### 表 A：`110051`（独立 19 文件，`19,000` 事件；`NaN` 计数含全 20 文件）逐级累加指标

（坐标系：旋转 ACTS 系 `A`；角度单位 `mrad`；动量残差无量纲百分数；轨迹质量筛选：GitHub 三项切选 + `purity ≥ 0.7`）

| 步骤配置 | 20f 原始 / 选中 `NaN` seed | 干净 seed `δp_st1` 中位数 | 干净 seed `Δφ_A` 中位数 `(e⁻ / e⁺)` [mrad] | `σ_68(Δφ_A)` / `σ_68(Δθ_A)` [mrad] | S2 首片新息 pull RMS | 同真值 S2 首片 `flt` 中位数（`< −0.5` 占比） | CKF 跨粒子串扰率 | Fiducial 单轨 / 双轨配对效率 | 原始 `good-seed` 同真值数 `N`（S2 `smt` 中位数） | 最终 refit `δp_st1` 中位数（`σ_68`） |
| --- | ---: | ---: | --- | --- | ---: | --- | ---: | --- | --- | --- |
| `C0_baseline` | `1,971 / 750` | `+16.35%` | `+0.114 / −0.113` | `0.185 / 0.618` | `0.272` | `−68.41%`（`60.54%`） | `32.54%` | `74.90% / 61.58%` | `711`（**`−11.71%`**） | `−0.01%`（`21.96%`） |
| `C1_cov_loc` | `1,971 / 750` | `+16.35%` | `+0.114 / −0.113` | `0.185 / 0.618` | `0.283` | `−71.05%`（`62.91%`） | **`13.53%`** | `71.57% / 56.83%` | `780`（`−11.53%`） | `+0.07%`（`22.71%`） |
| `C2_orig_var` | `1,971 / 750` | `+16.35%` | `+0.114 / −0.113` | `0.185 / 0.618` | **`1.076`** | **`−27.18%`（`34.79%`）** | **`9.37%`** | **`79.05% / 67.52%`** | `917`（`−11.15%`） | `+0.04%`（`21.48%`） |
| `C3_seed_b` | `1,971 / 750` | **`+0.20%`** | `+0.114 / −0.113` | `0.185 / 0.618` | `1.058` | `−30.90%`（`36.27%`） | `9.11%` | `78.93% / 67.45%` | **`1,893`（`−0.51%`）** | `+0.00%`（`21.52%`） |
| `C4_no_nan` | **`0 / 0`** | `+0.20%` | `+0.114 / −0.113` | `0.185 / 0.618` | `1.058` | `−30.91%`（`36.28%`） | `9.10%` | `78.92% / 67.46%` | `1,894`（`−0.51%`） | `−0.00%`（`21.53%`） |
| `C5_local_tan` | **`0 / 0`** | `+0.20%` | **`+0.001 / +0.007`** | **`0.034 / 0.322`** | `1.058` | `−31.27%`（`36.58%`） | **`8.01%`** | **`78.99% / 67.59%`** | **`2,740`（`−0.41%`）** | `−0.00%`（`21.52%`） |

#### 表 B：其余四类样本（`110033`, `100116`, `100117`, `100130`，各 `10,000` 事件）逐级关键改善汇总

| 样本 DSID | 原始 `NaN` seed（`C3 → C4`） | Seed `σ_68(Δφ_A)` [mrad]（`C4 → C5`） | Seed `σ_68(Δθ_A)` [mrad]（`C4 → C5`） | S2 首片新息 pull RMS（`C0 → C2 → C5`） | CKF 跨粒子串扰率（`C0 → C1 → C2 → C5`） | Fiducial CKF 效率（`C0 → C1 → C2 → C5`） | 原始 `good-seed` 轨迹数与 S2 `smt` 中位数（`C0 → C3 → C5`） |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `110033`（`e⁺e⁻` 软谱） | **`1,804 → 0`** | **`1.546 → 0.093`**（`16.6×`） | **`0.625 → 0.329`** | `0.268 → 1.015 → 1.003` | `13.76% → 5.41% → 2.52% → 2.31%` | `85.05% → 82.33% → 86.32% → 86.87%`（配对 `75.53% → 79.12%`） | `342 (−10.05%) → 1,390 (−0.79%) → 2,101 (−0.70%)` |
| `100116`（`PG μ⁺`） | **`3,088 → 0`** | **`1.357 → 0.112`**（`12.1×`） | **`0.668 → 0.393`** | `0.208 → 0.918 → 0.916` | `0.00%`（单径迹） | `99.45% → 99.48% → 99.50% → 99.52%` | `383 (−8.77%) → 439 (−5.34%) → 917 (−5.29% vs p_prod)` |
| `100117`（`PG μ⁻`） | **`3,042 → 0`** | **`1.400 → 0.105`**（`13.3×`） | **`0.663 → 0.386`** | `0.209 → 0.922 → 0.919` | `0.00%`（单径迹） | `99.41% → 99.55% → 99.55% → 99.55%` | `351 (−8.81%) → 475 (−5.06%) → 934 (−4.88% vs p_prod)` |
| `100130`（`PG K_S` 强子） | **`7,158 → 0`** | **`3.485 → 0.178`**（`19.6×`） | **`0.628 → 0.335`** | `0.351 → 0.878 → 0.884` | `2.42% → 1.26% → 0.22% → 0.22%` | `23.73% → 23.17% → 23.63% → 23.39%`（高多重数钨靶簇射） | `30 (−5.67%) → 141 (−0.58%) → 204 (−0.65%)` |

> **低动量样本上 `C5`（局部圆切线）效应放大的物理机制**：跨第一磁铁（`Δz ≈ 1200 mm`）弦线与 Station 1 局部切线的夹角偏差满足 `Δ(dy_F/dz_F) ≈ L_eff / (2 R) = 0.15 B_eff L_eff / p`，严格与动量 `p` 成反比。在 `921 GeV` 的 `110051` 上该弦线误差为 `±0.11 mrad`，而在 `32–91 GeV` 的 `110033`、`100116/100117` 与 `100130` 上放大至 **`1.4–3.5 mrad`**，因此 `C5` 在软谱电子、μ 子与强子样本上带来了 **12 至 20 倍**的初始角度精度提升。

![Step-by-step momentum impact C0 to C5](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/ablation_step_by_step_momentum.png)

![Observable-specific necessity of C1 C2 C4 C5](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/ablation_kinematics_necessity.png)

---

### 3. 原宏规范九层残差分布对照（`ckf_residual_perLayer.C`）

以下分布严格采用参考仓库 `ckf_residual_perLayer.C`（`chi2Sum/NDF < 20 && rMax_smt < 90 && NDF > 9` + `purity > 0.7`，`100` bins，`[−1.2, 1.2]`，绿 `prt` / 红 `flt` / 蓝 `smt`，`REFERENCE_SCOPE_LABEL` 标明独立文件与事件覆盖）：

1. **无 seed 真值门限（`110051` 独立 19 文件，`19,000` 事件）从 `C0` → `C2` → `C5`**：
   - `C0` 中 `S1 L0–L2` 的 `prt/flt` 峰因 `B = 0.55 T` 右偏在 `+0.16` 且 `S1 L2` 已向左发散，`S2 L0` 出现高度 `3,000`/bin 的 `−0.92` 增益爆炸红峰；
   - `C2` 恢复 `[100..1000]` 方差后，`S2 L0` 的 `−0.92` 尖峰高度下降近 3 倍，`S1 L2` 发散消失，仅余 `+0.16` seed 标度偏移；
   - `C5` 完成全部 5 项修复后，**`S1 L0–L2` 的绿 `prt`、红 `flt`、蓝 `smt` 三线完全重合居中于 `0.00`**，`S2 L0` 的绿 `prt` 主峰对准 `0.00`。
2. **原始 `good-seed` 门限下（`110051` 独立 19 文件，`19,000` 事件）`C0` vs `C5`**：
   - `C0` 中蓝线 `smt` 在全部 9 层左移至 `−0.12`（峰值约 `235`/bin），`S2 L0` 红/绿线堆积在 `−0.92`；
   - `C5` 中**全部 9 层的蓝线 `smt` 对称居中于 `0.00`（峰值跃升至 `1,000`/bin，超 `4` 倍），`S1` 与 `S3` 三线完全重合，`S2 L0` 的 `−0.92` 假峰彻底消失**。

![C0 ungated 19 files nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/110051_C0_baseline_19f_ungated_nine_layers.png)

![C2 ungated 19 files nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/110051_C2_orig_var_19f_ungated_nine_layers.png)

![C5 ungated 19 files nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/110051_C5_local_tan_19f_ungated_nine_layers.png)

![C0 good-seed 19 files nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/110051_C0_baseline_19f_goodseed_nine_layers.png)

![C5 good-seed 19 files nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/110051_C5_local_tan_19f_goodseed_nine_layers.png)

其余 7 幅九层参考分布图（全 20 文件 `C0`、独立 19 文件 `C1`、`C3`、`C4`、`H1_electron`、`S3_strip_only`、`S2_strip_fit`）均收录于 18 页汇总图集 [report.pdf](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/report.pdf)。

---

### 4. 探索性检查一：相同物质效应下 `ParticleHypothesis::muon()` vs `ParticleHypothesis::electron()`

在完全相同的源事件、完全相同的 `material-latest-rowmajor.json` 及 `CkfMultipleScattering = true, CkfEnergyLoss = true` 设置下对比 `C5_local_tan`（`muon()`）与 `H1_electron_hypo`（`electron()`），在第一块有源测量晶圆 `w0`（旋转 ACTS 系 `A`）上严格匹配 `FaserSiHit` 局部真值与 `p_st1` 评估 5 参数协方差 pull：

| 样本 DSID | 粒子假设 | Fiducial CKF 效率（双轨配对效率） | Refit 成功率 | 最终 `δp_st1` 中位数（全体 / 无辐射 / 硬辐射 `<0.95`） | 最终 `σ_68(δp_st1)` | Post-refit `χ²/ndf` 中位数 | 首片晶圆 `w0` Refit Pull `σ_68`: `(loc0, loc1, φ, θ, q/p_全体, q/p_无辐射)` |
| --- | --- | --- | ---: | --- | ---: | ---: | --- |
| `110051_heldout19` | `muon()` (`C5`) | `78.99%`（`67.59%`） | `99.87%` | `−0.00% / +2.31% / −9.12%` | `21.52%` | `1.066` | `(1.60, 1.07, 1.75, 1.11, 1.85, 1.42)` |
| `110051_heldout19` | `electron()` (`H1`) | `78.93%`（`67.51%`） | `99.87%` | `−0.01% / +2.30% / −9.16%` | `21.51%` | `1.065` | `(1.59, 1.07, 1.69, 1.11, 1.78, 1.38)` |
| `110033`（`e⁺e⁻`） | `muon()` (`C5`) | `86.87%`（`79.12%`） | `99.64%` | `−0.91% / +0.06% / −8.38%` | `6.93%` | `1.056` | `(1.27, 1.07, 1.52, 1.09, 1.94, 1.24)` |
| `110033`（`e⁺e⁻`） | `electron()` (`H1`) | `86.87%`（`79.15%`） | `99.66%` | `−0.92% / +0.07% / −8.40%` | `6.96%` | `1.054` | **`(1.27, 1.07, 1.30, 1.07, 1.68, 1.08)`** |
| `100116`（`PG μ⁺`） | `muon()` (`C5`) | `99.52%`（—） | `98.63%` | `−1.69% / −1.69% / −13.01%` | `3.40%` | `0.989` | `(1.06, 1.01, 1.18, 1.04, 2.12, 2.12)` |
| `100116`（`PG μ⁺`） | `electron()` (`H1`) | `99.52%`（—） | `98.65%` | `−1.71% / −1.71% / −12.48%` | `3.38%` | `0.984` | `(1.06, 1.01, 0.99, 1.03, 1.67, 1.66)` |
| `100117`（`PG μ⁻`） | `muon()` (`C5`) | `99.55%`（—） | `98.69%` | `−1.54% / −1.54% / −8.00%` | `3.47%` | `0.977` | `(1.05, 1.02, 1.17, 1.03, 1.99, 1.99)` |
| `100117`（`PG μ⁻`） | `electron()` (`H1`) | `99.55%`（—） | `98.69%` | `−1.55% / −1.55% / −8.35%` | `3.48%` | `0.973` | `(1.05, 1.02, 0.99, 1.02, 1.59, 1.59)` |
| `100130`（`PG K_S`） | `muon()` (`C5`) | `23.39%`（`6.59%`） | `88.00%` | `−0.80% / −0.80% / −3.29%` | `6.48%` | `1.024` | `(1.04, 1.07, 1.17, 1.02, 1.30, 1.30)` |
| `100130`（`PG K_S`） | `electron()` (`H1`) | `23.38%`（`6.59%`） | `88.04%` | `−0.78% / −0.78% / −3.32%` | `6.48%` | `1.003` | `(1.04, 1.07, 1.06, 1.00, 1.12, 1.12)` |

> **注**：`100116/100117` 的 `q/p` pull 是相对 `p_st1`（含最后一次离散顶点后约 `0.9 GeV` 未记录连续钨/铅电离能损）计算的，因而存在整体正偏移（`100116` `muon()` pull 中位数 `+0.94`）；位置与角度 pull `(loc0, loc1, φ, θ)` 直接使用 `w0` 上的局部 `FaserSiHit` 真值，在 `muon()` 假设下已高度接近 `1.0`（`1.05, 1.01, 1.17, 1.03`）。

![Particle hypothesis comparison muon vs electron](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/particle_hypothesis_comparison.png)

---

### 5. 探索性检查二：18-cluster 双磁铁积分场拟合（`S3`）与笛卡尔积预筛选/强化去重（`S1, S2`）

在 `110051_heldout19` 上，18-cluster 双磁铁积分场拟合的约化卡方 `χ²_strip/ndf`（`ndf ≈ 13`）对干净单粒子 3-station seed 的中位数仅 **`0.98`**（仅 `1.80%` 超过 `35`，对比原圆拟合 `3.68%` 超过 `35`），而对跨正负电子混合假组合显著右移至高卡方长尾（图 5 左）。下表对比标准组合圆拟合（`C5`）、纯 18-cluster 积分场拟合（`S3_strip_only`）、预筛选 + 强化去重（`S1_dedup_prefilter`）及两者叠加（`S2_strip_field_fit`）：

| 样本 DSID | 配置 | 选中 Seed 数 | 跨正负主粒子假 Seed 占比 | 全 Seed `< −0.5` 占比 | Fiducial Seed → CKF 效率（双轨配对效率） | GitHub 切选内假轨迹率 | CKF 跨粒子串扰率 | 最终 `δp_st1` 中位数（`σ_68`） | `good-seed` 同真值轨迹数 `N` |
| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: |
| `110051_heldout19` | `C5_local_tan` | `49,749` | `39.9%` | `36.5%` | `82.76% → 78.99%`（`67.59%`） | `24.48%` | `8.01%` | `−0.00%`（`21.52%`） | `2,740` |
| `110051_heldout19` | `S3_strip_only` | `50,108` | `41.1%` | `36.4%` | `81.53% → 78.52%`（`67.11%`） | `25.11%` | **`7.51%`** | `−0.04%`（**`21.43%`**） | **`2,792`** |
| `110051_heldout19` | `S1_dedup_prefilter` | **`33,638`** | **`32.1%`** | **`21.6%`** | `74.65% → 76.36%`（`62.90%`） | **`22.57%`** | `12.03%` | `−0.09%`（`21.91%`） | `2,606` |
| `110051_heldout19` | `S2_strip_field_fit` | **`33,884`** | `33.5%` | `22.8%` | `73.54% → 75.87%`（`62.31%`） | `23.19%` | `11.68%` | `−0.08%`（`21.77%`） | `2,595` |
| `110033`（`e⁺e⁻`） | `C5_local_tan` | `31,201` | `26.1%` | `28.9%` | `89.53% → 86.87%`（`79.12%`） | `6.97%` | `2.31%` | `−0.91%`（`6.93%`） | `2,101` |
| `110033`（`e⁺e⁻`） | `S3_strip_only` | `31,275` | `26.4%` | `28.8%` | **`89.36% → 86.85%`**（`78.94%`） | `7.09%` | **`2.17%`** | **`−0.90%`（`6.90%`）** | **`2,167`** |
| `110033`（`e⁺e⁻`） | `S1_dedup_prefilter` | **`15,954`** | **`13.0%`** | **`12.8%`** | `80.59% → 81.28%`（`67.72%`） | **`6.43%`** | `3.70%` | `−0.91%`（`7.08%`） | `1,950` |
| `110033`（`e⁺e⁻`） | `S2_strip_field_fit` | **`16,098`** | `13.5%` | `13.3%` | `80.69% → 81.42%`（`67.86%`） | `6.44%` | `3.65%` | `−0.91%`（`7.09%`） | `1,976` |
| `100116`（`PG μ⁺`） | `C5` / `S3` / `S2` | `9,230 / 9,232 / 9,176` | `0.0%` | `5.2%` | **`99.52% / 99.50% / 99.53%`** | `0.00%` | `0.00%` | `−1.69%`（`3.40% / 3.40% / 3.52%`） | `917 / 922 / 874` |
| `100117`（`PG μ⁻`） | `C5` / `S3` / `S2` | `9,178 / 9,177 / 9,148` | `0.0%` | `5.6%` | **`99.55% / 99.53% / 99.58%`** | `0.04%` | `0.00%` | `−1.54%`（`3.47% / 3.47% / 3.56%`） | `934 / 915 / 870` |
| `100130`（`PG K_S`） | `C5` / `S3` / `S2` | `11,436 / 11,416 / 5,607` | **`10.8% / 10.8% / 0.3%`** | `45.1% / 45.0% / 23.0%` | `23.39% / 23.38% / 23.05%` | `20.70% / 20.48% / 20.17%` | `0.22% / 0.22% / 0.23%` | `−0.80% / −0.80% / −0.87%`（`6.48% / 6.31% / 6.33%`） | `204 / 199 / 187` |

![Exploratory Check 2 seeding combinatorics and 18-cluster strip fit](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/seeding_combinatorics_checks.png)

![Cross-sample generalization summary across all 5 samples](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/multisample_summary_comparison.png)

---

## 执行记录、资源验证与产物索引

- **HTCondor 调度记录（全部选择已验证兼容且普通 `bigbird` 中负载最低的 `bigbird20.cern.ch`）**：
  1. 主矩阵 9 组 × 33 分片 = `297` 个双核作业 `latest_val_campaign_v1_20261118`（Cluster **`1584737`**，`2 CPUs / 6400 MB / 8 GB scratch`，`540,000` 事件）：**`297/297` 全部 `exit_code == 0`**，清单见 [submission.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/latest_val_campaign_v1_20261118_worker/submission.json)；
  2. 补充纯积分场拟合对照 `S3_strip_only` 共 `33` 个双核作业 `latest_val_s3_only_20261118`（Cluster **`1584791`**，`2 CPUs / 6400 MB / 8 GB scratch / espresso`，`60,000` 事件）：**`33/33` 全部 `exit_code == 0`**，清单见 [submission.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/latest_val_s3_only_20261118_worker/submission.json)；
  3. 全样本 10 组 × 5 样本多线程汇总与 18 页图表导出 `latest_val_study_v3_20261118`（Cluster **`1584793`**，`4 CPUs / 14000 MB / 6 GB / espresso`）：**`exit_code == 0`**，清单见 [submission.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/latest_val_study_v3_20261118_worker/submission.json)。
- **图表逐 bin 闭合与回归单测**：
  - 12 幅九层参考分布图共 `324` 个 `TH1`、含 underflow/overflow 共 `33,048` 个分箱与 worker 原宏输出 100% 逐 bin 闭合（[figure_validation.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/figure_validation.json)），18 页 [report.pdf](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/report.pdf) 全部通过目视核验；
  - 仓库单元测试 `python3 -m unittest discover -s tests -v`（含新增 [test_latest_validation_fixes.py](../tests/test_latest_validation_fixes.py)，覆盖旋转平面 `covLoc0/covLoc1`、`getY()` 负判别式与严格弱序、局部圆切线消除正负电荷弦线偏置、强子/轻子 PDG 电荷映射、18-cluster 双磁铁积分场线性最小二乘闭合、并集/`S23` 去重保护合并 S1 双轨及 `K_S` 中性初态 Geant4 次级衰变链条码解析）共 **115 项全部通过**。
- **关键代码与数据产物路径**：
  - 仓库脚本与配置：[latest_validation_tracking.py](../configs/latest_validation_tracking.py)、[multi_sample_field_closure.py](../configs/multi_sample_field_closure.py)、[run_validation_shard.py](../scripts/run_validation_shard.py)、[prepare_validation_campaign.py](../scripts/prepare_validation_campaign.py)、[latest_validation_study.py](../analysis/latest_validation_study.py)；
  - 统计与图集索引：[summary.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/summary.json)、[field_closure_summary.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/field_closure_summary.json)、[bin_counts.npz](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/bin_counts.npz)、[reference_bin_counts.npz](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/reference_bin_counts.npz)、[report.pdf](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/full_validation_v3/report.pdf)。

---

## 事实与假设区分

- **已证实事实**：
  1. `SeedBFieldTesla = 0.473662 T` 在电子（`110051, 110033`）、正负 μ 子（`100116, 100117`）与 `K_S` 强子（`100130`）五类样本及 `r < 100 mm` 全孔径内确定性闭合偏差均在 `0.012%` 以内，全量 RDO 干净 seed 标度一致对齐；
  2. 五项默认修复各司其职且缺一不可：`C1` 压低跨径迹串扰、`C2` 恢复 S2 测量新息 pull 至 `~1.0` 并消除 `-0.92` 增益爆炸尖峰及恢复双轨效率、`C3` 消除 `+16.1%` seed 动量标度偏置与 `good-seed` 条件左移、`C4` 清零全部五类样本数千条 `NaN` seed、`C5` 消除正负电荷反号弦线弯转角偏置并将软谱样本初始 `φ_A` 精度提升 12–20 倍；
  3. 单高斯 CKF/`KalmanFitter` 下将 `ParticleHypothesis` 切为 `electron()` 不改变硬轫致辐射动量偏置与效率，仅增大过程方差使电子 `q/p` pull 更接近 `1.0`，同时对 μ 子/强子造成非物理误差棒膨胀；
  4. 纯 18-cluster 双磁铁积分场拟合（`S3`）在不损失效率的前提下进一步降低串扰并提升好 seed 产额，而前端硬 `χ²/ndf ≤ 35` 与 `S2` 共享否决（`S1, S2`）虽削减 32%–51% 假组合却折损 2.6%–5.5% 轫致辐射/高能合并双电子效率。
- **待进一步验证的假设**：
  1. 若在 Calypso 中接入 `Acts::GaussianSumFitter`（Bethe-Heitler 多高斯分量拟合器）专门对电子候选重拟合，预期可在不改变 μ 子重建的前提下实质压缩 `110051/110033` 的硬轫致辐射负尾；
  2. 将 `S3_strip_only` 前端排序与 CKF 后的 `Acts::AmbiguityResolutionAlg`（按完整轨迹测量数、孔洞数与轨迹 `χ²/ndf` 全局消歧）组合，预期可在保持 `C5/S3` 双电子配对效率的同时剔除残余跨径迹假轨迹。

---

## 下一步计划

1. **电子通道 `Acts::GaussianSumFitter`（GSF）重拟合原型评估**：
   - **具体任务**：在隔离 worktree 中为电子假设轨迹配置 ACTS GSF 组件（`BetheHeitlerApprox`，6–12 分量），保持前端 `C5` CKF 找迹不变，仅对比单高斯 `KalmanFitter` 与 `GSF` 在 `110051` 和 `110033` 硬轫致辐射子集（`p_min / p_st1 < 0.95`）上的动量中位数、`σ_68` 与 pull 标定；
   - **所需输入**：已冻结的 `110051`（独立 19 文件）与 `110033`（5 文件）RDO 样本及 `material-latest-rowmajor.json`；
   - **验收标准**：在相同 CKF 找到的轨迹集合上量化 GSF 相对单高斯 `KalmanFitter` 的动量尾改善幅度、单事件 CPU 开销与对非辐射轨迹分辨率的影响。
2. **4-station（含 IFT）与 Backward 模式下的积分场 Seed 及轨迹级消歧验证**：
   - **具体任务**：测试 `SeedStripFieldFit`（`S3`）在开启 Station 0（IFT）及反向寻迹（`BackwardPropagation = True`）时的普适性（自动按积分场替代固定 `0.473662 T` 常数），并在 CKF 后引入轨迹级共享测量消歧以替代前端硬删 seed；
   - **所需输入**：`110051`、`110033`、`100116/100117` 与 `100130` 样本；
   - **验收标准**：双轻子配对效率不低于 `C5`（`110051 ≥ 67.6%`, `110033 ≥ 79.1%`），同时将 GitHub 切选内假轨迹率压低至 `S1` 水平以下。
