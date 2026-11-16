# S04：latest Calypso 的 good-seed 动量左偏与第二站状态溯源

## 阶段目标

在用户指定的 `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest` 中，检验旧 good-seed 测试的 seed、Station 2 predicted 和 smoothed 左偏是否构成同一个因果链。逐项核对真实 Calypso 的 seed、传播、measurement update、候选归属、smoothing 和材料查表。

结果根目录：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1`，下文记为 `R`。本轮使用 latest 的实际重建及匹配库；控制实验也在该版本的原生组件中运行。

## 阶段性结论

**第二站九层图的强左峰主要是中间曲率分布过宽，经 `p=1/abs(q/p)` 变换形成的低动量峰；需要区分本站第一 sensor 的 predicted 和后续 sensor 的 predicted。** 同 seed/track truth 的 725 条轨迹中，首次 S2 predicted 的 residual 中位数是 **+4.72%**，首次 S2 filtered 变为 **−73.17%**，下一 sensor predicted 继承该结果。九层图把同层两个 sensor 合并，因而绿色 predicted 中包含首次更新之后的强左峰。

**good-seed 子样本中约 −12% 的 smoothed 中心，主要由 seed 动量筛选造成。** 固定同一个全量重建结果、同一 seed/track truth，不加 seed gate 时中位数为 **−2.41%**；只加 seed 动量 ±10% 条件就变为 **−11.21%**；只筛电荷和 theta 时为 **−1.68%**；三条全加为 **−11.71%**。seed 和最终拟合共享测量，seed 动量筛选不是独立的质量检验。

**原生 seed-model closure 找到与上述选择偏移一致的具体来源：** 真实完整磁场下传播的200个电子，用原生 CircleFit 与固定0.55 T公式反推动量，中位数高估16.12%，对应圆拟合有效场0.47366 T。seed近似与完整场拟合的响应不同；将前者筛选到truth附近，会偏向后者的低动量子集。

**独立文件验证支持该机制：** 用第一文件100事件测得的1.16116响应系数，仅校准其余19文件的 seed-p 分析窗口，保持所有重建参数不变。S2 smoothed 中位数由原 gate 的 −11.674% 变为 −0.510%；因此主要左移可以在不改变 predicted/filter/smoothing 计算的情况下由选样解释。

这些机制需要分开解释。真实 Calypso 的先验控制显示，抑制前向宽尾不必然消除 smoothed 的条件偏移。不能据 good-seed 图直接认定全体轨迹具有 −12% 的动量标定错误，也不能把 truth 辅助窄先验当成 production 修复。

## 软件、材料、样本与选择

| 项目 | 实际设置 |
| --- | --- |
| latest commit | `7fd40df8cd3b60e81b4095a4f84310c4a54700a0` |
| 原始源码 | `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest`，保留 clean |
| 诊断 worktree | `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_momentum_diag`，同 commit |
| 重建库 | `R/software_v4/lib/libFaserActsKalmanFilter.so` |
| Python config | `R/software_v4/python/FaserActsKalmanFilter/`，来自同一 latest 源码 |
| 依赖 | latest `run/` 中的 geometry、RootWriter、TrackingTruthUtils、TrackerSegmentFit；Athena 24.0.41 / EL9 gcc13；ACTS v32.0.2 |
| 材料来源 | latest `Tracking/Acts/FaserActsGeometry/data/material-maps-alma9-rotated-rowmajor.json.gz` |
| 实际加载材料 | `R/material-latest-rowmajor.json`，上述文件原样解压，11,866,364 bytes |
| 输入 | 110051 `mc24/foresee/110051/rdo/master_dev/` 全20文件，每文件 `[0,1000)`，共20,000源 entries |
| geometry / conditions | `FASERNU-04` / `OFLCOND-FASER-05` |
| 重建 | 3station forward；mask `[0,1,2,3,4,5]`；candidate-isolated smoothing 开 |
| 原宏选择 | `chi2Sum/NDF<20 && rMax_smt<90 && NDF>9` |
| 本次额外选择 | strict final SDO purity `>0.7`；原始 ancestor seed 的三条 good-seed 条件 |
| 九层主图 truth | native writer `t_eQOP`，即 SiHit link 对应的 HepMC particle momentum |
| 操作分解 truth | 同 writer majority particle 的 production momentum，单位 GeV；曲率单位 GeV⁻¹ |
| frame | ACTS A：`(x_A,y_A,z_A)=(z_F,y_F,-x_F)` |

good-seed 条件保持旧测试定义：

```text
abs(p_seed / p_seed_truth - 1) < 0.1
qop_seed * qop_seed_truth > 0
abs(theta_seed - theta_seed_truth) < 1e-4 rad
```

没有增加 phi cut。truth-mean 控制仍使用**原始 seed 的 gate**，避免“注入真值后所有 seed 自动过关”的选样变化；实际 CKF 输入均值另行记录。seed/track 多数粒子不同的轨迹单列，不用重建电荷或动量参与主 truth matching。

最新材料含22个 surface entries，旧 AFS JSON 含12个。共有 entries 的 raw cells 相同，但原 AFS 文件的 rotation 序列与 latest row-major 序列不同，translation 还相差2 mm；不能将两个 JSON 当成同一份输入。本轮未使用旧 AFS 文件。实际 ACTS decoder 及 Calypso surface 的 matrix、local coordinate、bin indices、raw thickness/material 参数已在372次 native 查询、9个 surfaces 上与冻结 latest JSON 核对；见 [smoke_validation.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/smoke_validation.json) 和 [材料差异](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/material_vs_afs_detail.json)。该 native 验证范围是20-event smoke，不宣称遍历全部22个 entries。

原宏来自用户提供的 [test_pResolution](https://github.com/zsm2024/FASER_tracking/tree/main/test_pResolution)，冻结 commit `33ca2281b961935456dfd4dc9fc06235e24b551e`。复用旧 good-seed 阶段记录过的 adapted macro，保留实际启用公式、原 cuts、100 bins、`[-1.2,1.2]`、raw Entries、P绿/F红/S蓝和九层布局；本次额外 purity/seed gate 明示于图。ROOT TTreeFormula 的除零、NaN、flow 行为保持原样，不把公式约成更改这些行为的表达式。

## 新增诊断与重建不变性

[0048_latest_momentum_audit.patch](../patches/0048_latest_momentum_audit.patch) 与 [LatestMomentumAudit.h](../patches/LatestMomentumAudit.h) 为完整 tracking 诊断改动。补丁在 finder 前后记录 seed、selected path、P/F/S 的六参数及 covariance、measurement/projector、material 查表和 hit 身份；默认关闭的控制可替换 seed 均值或切换原生 CKF 的 energy loss、scattering、step size。

原 CKF 的 `nSeed` 未递增，`trackGroup=nSeed-1` 恒为 UINT32_MAX；诊断改用 finder 循环的 `iseed` 保存真实祖先。该字段不参与本配置的 ranking/fitting，alignment writer 关闭。isolated smoothing 使用 latest 自带实现。所有控制保留真实 Calypso 的 hit finding、candidate selection 和 native ACTS updater/smoother，没有另写 tracking。

完整 FaserActsKalmanFilter package 按 latest CMake 实际存在的13个 translation units 编译，使用 Athena 24.0.41 匹配编译/链接配方，依赖重定向 latest `run/`。构建及完整配方见 `R/software_v4/software.json`、`compile_recipes.json`、`link_recipe.json`。旧 CMake 中三个已不存在的 source patterns 不产生编译单元，已列入 manifest。

同一20-event输入的 native baseline 与 audit baseline，59个 seeds、43条 writer tracks 的全部保存 scalars、P/F/S `q/p` 与误差、seed fields、offsets **逐值精确一致，NaN 一致**。731次可用 measurement update 满足 `Δ(q/p)=K r`，最大绝对误差 `1.39e-17 GeV⁻¹`。补丁与实际编译源码逐项一致，reverse apply 检查通过；见 `R/source_validation.json`。全量归约额外核对每个 audit measurement state 与 ROOT writer 的 geometry id、P/F/S 数值和顺序。

## 详细数据：第二站第一 sensor 与下一 sensor

全量 baseline 输出42,151条 writer tracks；good-seed 选中948条，其中223条 seed majority 与 writer track majority 不同。下面只用**相同 majority 的725条**，每条轨迹取 forward path 首个 S2 measurement，分母是 track，不是九层图的 measurement-state Entries。

| 状态 | momentum residual 中位数 | residual <−0.1 | 临时错电荷 |
| --- | ---: | ---: | ---: |
| 原始 seed | +2.741% | 0% | 0% |
| 最后 S1 filtered | +4.716% | 35.862% | 16.138% |
| 首次 S2 predicted | +4.716% | 35.862% | 16.138% |
| 首次 S2 filtered | −73.166% | 78.207% | 43.034% |
| 下一 S2 sensor predicted | −73.166% | 78.207% | 43.034% |
| 首次 S2 smoothed | −11.706% | 61.793% | 0.828% |

首次 S2 predicted 的 signed curvature error `qop/qop_truth−1` 标准差为1.1033，首次 filtered 增至8.4859；其中位数从 −0.0990 变为 −0.3067。**这是两侧曲率显著展宽，并非所有曲率向同一个方向平移。** 对很大正/负曲率，`p=1/abs(qop)` 都给出小动量，从而把两侧折叠到 residual 接近 −1 的区域。

包含所有948条 gate 通过轨迹时，912条的首次及下一次 S2 measurement 均位于 layer0。首次 filtered 到下一 predicted 的 relative-p 变化中位数和95%分位均为0；少数材料状态产生微小变化。九层 S2L0 predicted 因混合这两个 sensor，形状不能代表“刚进入 S2 时的预测”。

same-truth 中有700条的两个measurement均在S2L0，单列为 summary 的 `same_truth_and_both_S2_layer0`，避免少量缺失前部measurement的轨迹混淆层号。全948条首次 S2 的 `corr(q/p,phi)` 中位数 −0.997656、gain中位数 −0.0100916 GeV⁻¹/mm；innovation sigma中位数3.094 mm，实际 residual标准差0.754 mm，pull标准差0.241。选中updates的均值恒等式最大误差 `6.94e-18 GeV⁻¹`，covariance恒等式相对误差最大 `1.69e-12`；全部保存states未检出小于−1e−8的correlation-matrix eigenvalue。这些证据支持宽且高度相关的先验对中间更新的放大作用，没有显示本轮 updater 的代数实现错误。

入口 predicted 的左尾已经在 S1 更新产生：725条中的260条首次 S2 predicted `<−0.1`，其首次越界为 S1L2 update 249条、S1L1 update 11条；97条 `<−0.5` 均首次在 S1L2 update 越界。这是逐轨迹历史定位，不把“首次越界”标签本身当成完整因果证明。

从最后 S1 filtered 到首次 S2 predicted，全948条的 relative-p 变化中位数为 `−4.96e-7`，不是10%量级；累计 transport 的 signed 曲率变化远小于累计 measurement updates。因而本轮观测不支持“站间传播/能损直接把入口动量降低约10%”的解释。

![S2 sensors separated](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/full_controls_v12/station2_sensors.png)

## 详细数据：seed gate 产生的条件偏移

原生全部seed population共52,338条，50,364条有有限动量residual。其中位数 −1.312%，16%/84%分位为 −89.836%/+39.356%，`residual<−.5` 占35.216%。这不是单个+16%的窄峰：native-model closure解释的是估计器对无噪声完整场轨迹的系统响应，不能代替全部seed的错误组合、measurement/segment误差和辐射尾部分解。本轮主结论针对最终选中轨迹的good-seed条件中心，整个seed population的宽尾尚未逐类完成归因。

下表来自**同一个全20文件 baseline**，固定 native reconstruction，不重新重建，仅切换原始 ancestor seed 的分析条件；均要求 seed majority = writer track majority、GitHub cuts、final purity >0.7，并有 S1/S2 measurement。

| seed gate | tracks | 首次 S2 smoothed residual 中位数 | 原宏分箱的最高 bin |
| --- | ---: | ---: | --- |
| 无 | 17,126 | −2.406% | [−0.024,0] |
| 只选 p 在 ±10% | 3,039 | −11.208% | [−0.096,−0.072] |
| 只选同电荷 | 15,826 | −2.286% | [−0.024,0] |
| 只选 theta | 3,635 | −1.759% | [−0.024,0] |
| p + 同电荷 | 3,011 | −11.221% | [−0.096,−0.072] |
| 同电荷 + theta | 3,492 | −1.684% | [−0.024,0] |
| 原始三条 | 725 | −11.706% | [−0.120,−0.096] |

seed residual 与 smoothed residual 的二维图显示二者强相关，且二者的响应中心不同。seed 由相同站的 segment/hits 构造，再从这些 seeds 中选择“接近 truth”的子集，会选择特定测量涨落；后续拟合并不独立于这一选择。去掉 seed-p 条件后仍有约 −2%的中位偏移，其余来源不能从此表单独判定。

![Seed gate selection](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/full_controls_v12/seed_gate_selection.png)

## 详细数据：同路径原生 Calypso 控制

配对键为 **input file + source entry + 原 seed index + writer majority barcode + ordered measurement hits/stations**。baseline 限制为与每个控制完全相同的输入覆盖。表中的两个中位数是各自分布的中位数，不能拿二者之差充当逐轨迹差值的中位数；后者另存 summary 的 `difference`。候选新增/丢失导致的配对覆盖也单列。

### 全20文件控制

每组均20,000源事件。所有重建与归约任务exit0，逐文件核对相同input、完整 `[0,1000)`、实际加载库、材料与isolation。表中同路径配对允许 seed majority 与 writer majority 不同，但严格固定同一 writer majority；两者相同的子集另外计数。

| 控制 | writer tracks | good-seed tracks | 同路径配对 | 其中同 seed/track truth | 首次 S2 P 中位数：baseline→control | 首次 S2 F 中位数：baseline→control | 首次 S2 S 中位数：baseline→control |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| truth q/p，C不变 | 41,973 | 932 | 905 | 696 | +.04451→+.02172 | −.73014→−.72999 | −.10802→−.10784 |
| truth full mean，C不变 | 40,632 | 788 | 574 | 527 | +.05888→+.02235 | −.72706→−.73605 | −.11049→−.10976 |
| q/p inflation100 | 40,760 | 920 | 694 | 679 | +.05486→+.02882 | −.73274→+.05325 | −.11483→−.10977 |
| q/p inflation1 | 37,995 | 885 | 631 | 630 | +.06141→+.02782 | −.76503→+.02746 | −.11513→−.03581 |

truth-q/p 控制的逐轨迹 S residual差中位数为 `+3.88e-8`；truth-full-mean为 `+4.59e-5`。在same-seed/track-truth子集中，前者S中位数 `−.11751→−.11729`，后者 `−.11513→−.11439`。**即使修正seed完整mean，只要保留原宽covariance与原gate，主要S偏移和S2首次更新的宽尾仍在。**

inflation100 的同路径首次P左尾 `<−.1` 为 `247/694→6/694`，率差 −34.73 pp，按source-event聚类的95%区间 `[−38.32,−31.13] pp`；首次F为 `546/694→185/694`。S左尾仅从 `421/694→400/694`，率差 −3.03 pp，95%区间 `[−4.29,−1.76] pp`，S中位数仍约−11%。这确认“显著压低P/F宽尾，却没有移走S条件中心”的现象在全量上成立。

inflation1的同路径首次P与首次F左尾均为 `0/631`，S左尾 `385/631→120/631`，但重建writer tracks与配对覆盖也明显改变。不能只选其好看的条件分布而忽略效率、候选选择与未筛选总体。

### 第一文件的机制控制

以下均为第一文件1000源事件，不标为全量；baseline 29条 gate 通过轨迹。

| 控制 | 控制选中 tracks | 同 seed/hits/truth 配对 | 首次 S2 P 中位数：baseline→control | 首次 S2 F 中位数：baseline→control | 首次 S2 S 中位数：baseline→control |
| --- | ---: | ---: | ---: | ---: | ---: |
| 仅注入 truth q/p，C 不变 | 30 | 28 | −.05150→−.04963 | −.62616→−.62610 | −.12070→−.12073 |
| 注入完整 truth mean，C 不变 | 25 | 19 | −.05356→−.08078 | −.62854→−.72993 | −.10078→−.10026 |
| q/p inflation 100000→100 | 31 | 24 | −.05150→+.03423 | −.62616→+.07258 | −.12704→−.12421 |
| q/p inflation 100000→1 | 28 | 20 | −.01681→+.03485 | −.62616→+.03237 | −.12704→−.03215 |
| truth q/p + inflation1 | 28 | 20 | −.01681→−.00000475 | −.62616→+.00027 | −.12704→−.04899 |
| 关闭 CKF energy loss | 29 | 29 | −.05356→−.05355 | −.62854→−.62854 | −.11609→−.11609 |
| 关闭 CKF scattering | 30 | 29 | −.05356→−.05356 | −.62854→−.62854 | −.11609→−.12878 |
| CKF max step 10 mm | 23 | 23 | −.05356→−.05356 | −.62854→−.62859 | −.12877→−.12555 |

`inflation100` 将配对首次 S2 P `<−0.1` 从9/24降为0/24，首次 F 从19/24降为6/24，但 S 仍接近 −12%；这直接区分前向宽尾与 smoothed 条件偏移。`inflation1` 的首次 P/F 左尾均为0/20，S 从11/20降为1/20；该设置加强了已知 good-seed 信息，不能按分布变窄直接用于无 truth 的 production。

关闭 energy loss 的逐轨迹 residual 差中位数，P 为 `+1.90e-6`、S 为 `−6.58e-7`，尾部计数不变。关闭 scattering 的 paired S 差中位数仅 `−1.08e-5`，尽管两组分布中位数相差约0.013；这是分布排序与少量事件变化，不能表述成每条轨迹都移动了1.3%。10 mm step 控制的 paired S 差中位数 `+1.14e-5`，主要中心偏移仍在。以上材料/步长控制只覆盖该1000-event机制样本。

本样本首次 S2 predicted 的 innovation sigma 中位数3.079 mm，而实际 innovation residual 标准差0.513 mm，pull 标准差0.167；`corr(q/p,phi)` 中位数 −0.997655，gain 中位数 `−0.0101384 GeV⁻¹/mm`。宽且高度相关的 predicted covariance 允许一个亚毫米 residual 驱动很大的曲率改变量。选中 measurement updates 的均值恒等式最大误差 `8.67e-19 GeV⁻¹`，covariance 恒等式相对误差最大 `6.96e-13`。这与 updater 算错矩阵乘法的解释不符。

## truth 参考与真实辐射能损

native `t_eQOP` 与同 writer majority 的 production `q/p`，705,440个共同有效 states 的相对曲率差 RMS 为 `2.53e-8`，量级来自 writer 精度；不能用这两种定义的差别解释本轮10%偏移。native truth 缺失的额外覆盖另行统计，不把 raw Entries 等同于有限 residual 分母。

`SiHit->particleLink()->momentum()` 是关联 HepMC particle 的动量，**不是 sensor 处精确 G4 momentum**。seed writer 的 truth 也是 production particle 经 native extrapolator 传播到 seed 面的结果；该工具默认不加能损和 scattering，且用 muon hypothesis，不能称为精确局部 electron truth。

因此另用 Athena 读取全部20文件的 HepMC 同PDG continuation vertices，构造 forward、动量单调的上下界；前后无后续 vertex 时 lower=0，不进行插值。对 same-truth good-seed 轨迹，整个 tracking 跨度的 continuation lower/upper≥0.99 有157条；首次 S2 smoothed 的**真实局部 residual 下/上界中位数**分别为 −0.093298 / −0.092588。由此可排除“所有条件左偏都是真实大辐射能损”这一解释，同时不能据这157条推断每一条轨迹的精确局部动量。完整 bounds 与未解类别见 summary 的 `recorded_truth_bounds`。

## 真实 Calypso seed-model closure

为区分 measurement noise、真实辐射与 seed 圆拟合近似，使用 latest 的真实 `FaserActsExtrapolationTool`，在同一完整 geometry、field、conditions 和冻结材料中，把110051第一文件 `[0,100)` 的 production electrons 传播到 detector manager 的三个 station centers。选取向前、production 位于station1之前、p>10 GeV的电子；200个全部成功。关闭工具的能损/scattering，使用与原 seed truth 一致的 muon hypothesis；每个 station 用传播得到的真实位置、方向产生原来 ±30 mm points，直接调用 latest 编译库的 **`CircleFit::circleFit`**，随后采用原来 `.3 × .55 × R × .001` 公式。

这里调用了真实 Calypso 原有传播和 circlefit，没有编写传播器、场模型或 circlefit。此项是无measurement的确定性 seed-model closure，**不是重建样本，也不包含 G4 随机散射、digitization 或 segment fit 误差**；100-event覆盖与全量20,000-event重建分别标注。

| quantity | median | 16%–84% quantiles |
| --- | ---: | --- |
| native circle p / input p −1 | +0.161165 | +0.158340 — +0.161590 |
| `p_input / (.3 R)` 的有效场 | 0.473662 T | 0.473489 — 0.474817 T |

0.55 T 是 seed 估计器中的固定常数，不是对完整非均匀场积分的计算。本轮 closure 证实该估计器的响应中心约高16%，与全量二维 seed/smoothed 相关结构一致。作为量级检查，same-truth good-seed 的 seed residual 中位数+0.02741，除以closure比例1.16116得到约 −0.1152，接近实测 S 中位数 −0.1171；这是两个统计中位数的量级比较，不是对每条轨迹的恒等式，也没有据此硬编码新的常数场。

进一步冻结这个由file0前100事件测得的比例，对**其余19文件**进行独立选样检验，file0完全排除。仍要求same seed/track truth、GitHub cuts与final purity>0.7；校准gate保留原电荷和theta条件，只将动量窗口替换为 `abs((p_seed/p_seed_truth)/1.16116468−1)<0.1`。tracking与全部states均不变。

| 独立19文件选择 | tracks | S2 smoothed residual中位数 | `<−.1` / `>+.1` |
| --- | ---: | ---: | ---: |
| 无seed gate | 16,227，有限S为16,226 | −.02444 | 34.95% / 21.87% |
| 原seed gate | 702 | −.11674 | 61.82% / 0.57% |
| native-model校准后的seed gate | 1,468 | **−.00510** | **4.70% / 4.90%** |

原窗口选择的S明显左移，按独立native模型响应调整窗口后中心与左右10%尾部接近对称。此项直接支持**seed近似响应差 + 共享测量 + truth-p条件选择**产生主要S偏移的机制；它不是已经修改seed或修复production tracking的结果。

![Independent seed model gate validation](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/full_controls_v12/seed_model_gate_validation.png)

观察函数保存在 Calypso诊断 worktree 的 `Tracking/Acts/FaserActsKalmanFilter/src/LatestSeedFieldClosure.cxx`，仓库版本为 [LatestSeedFieldClosure.cxx](../patches/LatestSeedFieldClosure.cxx)；编译入口为 [build_latest_seed_field_closure.py](../scripts/build_latest_seed_field_closure.py)，Athena入口为 [latest_seed_field_closure.py](../configs/latest_seed_field_closure.py)。实际加载 helper 为 `R/software_field_closure_v3/libLatestSeedFieldClosure.so`，共享 same `software_v4/libFaserActsKalmanFilter.so` 和 latest geometry库。source、compile/link配方和loaded libraries保存在 `R/software_field_closure_v3/build.json`、[closure.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/seed_field_closure_v3/closure.json)。

早期 ROOT JIT 观察程序先出现 GeoPrimitives/Eigen plugin include错误，调整顺序后仍在JIT构造ACTS surface时崩溃，两次均exit129且无可用物理结果。最终采用与重建package相同配方原生编译后exit0，FPE auditor三类异常均为0。失败输出保留在 `R/seed_field_closure/`、`seed_field_closure_v2/`，有效结果仅用v3。

![Native seed model closure](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/full_controls_v12/seed_model_closure.png)

## 源码路径与物理机制

下表路径相对 `/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/`，行号对应未加 audit 的 latest 源码。

| 环节 | 真实源码 | 本轮核对内容 |
| --- | --- | --- |
| seed covariance | [Tracking/Acts/FaserActsKalmanFilter/python/CKF2Config.py:49](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/python/CKF2Config.py:49) | inflation `[100,100,100,10000,100000]`；Cqq=0.001 GeV⁻²、σq=0.0316228 GeV⁻¹；Cphi=.1089 rad²、σphi=.33 rad |
| seed measurements | [Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:96](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:96) | 读取 cluster covariance 但 measurement V 实际硬编码为 `.0231² mm²`；本轮保留该原生定义 |
| seed segment/fake points | [CircleFitTrackSeedTool.cxx:235](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:235) | segment fit position/momentum；position 及沿 segment direction 的 ±30 mm points |
| seed direction | [CircleFitTrackSeedTool.cxx:302](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:302) | 前两站 position 的 chord，不是 seed 面精确局部 tangent |
| seed momentum | [CircleFitTrackSeedTool.cxx:338](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:338)；[CircleFitTrackSeedTool.h:92](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.h:92) | 原生 CircleFit 半径乘 `.3 × B × .001`，默认固定 B=.55 T |
| seed reference/mean | [CircleFitTrackSeedTool.cxx:378](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx:378)；[CKF2.cxx:197](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx:197) | 旋转 seed/target surfaces；CKF2 重建实际 target 面，不能把 seedtool 旧 initialSurface 成员单独认定为本次原因 |
| native finder | [TrackFinderFunction.cxx:14](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/TrackFinderFunction.cxx:14) | native GainMatrixUpdater、EigenStepper、GainMatrixSmoother；控制改变同一个 CKF 的 options |
| candidate smoothing | [CKF2.cxx:277](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx:277) | 复制选中路径为 private states 后重算 smoothing，避免候选共享状态串扰 |
| writer stage | [CKF2.cxx:356](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx:356) | 九层 writer 输出 selected CKF states；不是 refit 后的 states；final purity 是 refit/final 轨迹的 cut |
| field rotation | [Tracking/Acts/FaserActsGeometry/FaserActsGeometry/FASERMagneticFieldWrapper.h:45](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsGeometry/FaserActsGeometry/FASERMagneticFieldWrapper.h:45) | inverse rotation 查真实 field，再旋转 B/gradient；本轮未发现明显方向/单位错误 |
| seed truth | [RootSeedWriteTool.cxx:959](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/FaserActsKalmanFilter/src/RootSeedWriteTool.cxx:959) | production particle → native extrapolator → seed reference surface |
| state truth | [Tracking/Acts/RootWriterTools/src/RootTrajectoryStatesWriterTool.cxx:426](/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/Tracking/Acts/RootWriterTools/src/RootTrajectoryStatesWriterTool.cxx:426) | SiHit link 对应 HepMC momentum |

实际 ACTS headers 位于 `/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/include/Acts/`。`TrackFinding/CombinatorialKalmanFilter.hpp:697` 起顺序为 transport covariance → pre-material → predicted → measurement update → 用 filtered 更新 stepper → post-material。`TrackFitting/GainMatrixUpdater.hpp` 实现 `q_f=q_p+K r`；`TrackFitting/GainMatrixSmoother.hpp:81` 起从末端 F=S 逆向结合未来测量。smoother 有后续两站的信息，前向 P/F 没有这些信息；三者的统计约束不同。

本轮已观测到 seed 宽 prior、S1 更新产生入口尾、S2 强曲率更新、gate 选择偏差，以及固定0.55 T seed近似在真实非均匀场中的约+16%响应。不能据此判定磁场 wrapper 或 field scale 有 bug。实际运行使用 `FaserFieldTable_v2.root`，COOL field scale 为1。

## 图表、执行和验收

![Latest full sample nine layers](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/full_controls_v12/110051_latest_good_seed.png)

全量及控制图的输入分箱与原宏计数逐 bin 核对，含 underflow、overflow 和 raw Entries；英文图保存 PNG、逐图 PDF、多页 PDF，附 `reference_bin_counts.npz`、`bin_counts.npz`、`histogram_metadata.json`、`figure_validation.json`。baseline 图和 controls 图各自标注输入覆盖，不把第一文件机制测试标为全量。

最终全量图集14页，135个reference histograms、13,770个含flow的bins全部与冻结计数精确一致；第一文件机制图集16页。14张全量图和第一文件9组九层图均已目视检查，正文相对/绝对链接全部可解析。全量最终summary与先前v11统计逐项完全相等，绘图修正没有改变物理统计。`git diff --check`通过。

[全量与四组控制图集 PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/full_controls_v12/report.pdf) · [第一文件全部机制控制图集 PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/pilot_controls_v10/report.pdf) · [全量统计 summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/full_controls_v12/summary.json) · [图索引与原宏逐bin验收](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/full_controls_v12/figure_validation.json) · [最终源码一致性](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/source_validation_final.json)。

入口代码：[prepare_latest_momentum.py](../scripts/prepare_latest_momentum.py)、[build_latest_momentum.py](../scripts/build_latest_momentum.py)、[latest_momentum_tasks.py](../scripts/latest_momentum_tasks.py)、[run_latest_momentum.py](../scripts/run_latest_momentum.py)、[latest_momentum.py](../analysis/latest_momentum.py)、[latest_momentum_study.py](../analysis/latest_momentum_study.py)、[latest_momentum_validation.py](../analysis/latest_momentum_validation.py)。

重建、编译、全量统计和绘图均使用 HTCondor。每次提交前保存实际 schedd 广告及队列，限制为已验证兼容的 `bigbird20.cern.ch`，排除已记录 EOS 兼容问题的 bigbird21；完整负载、选择时间、提交文件、cluster、资源见各 `*_worker/submission.json`。每文件1 CPU、3500 MB、6 GB scratch、workday；大型 audit/ROOT/xAOD 只留 worker scratch，归约后保留紧凑 NPZ、histograms 和压缩日志。输入身份使用 file/entry，重复生产 event number 不误删。

baseline 原20文件 cluster `1583106`；旧 reducer 失败的file02/07/14/15/18经终止后重跑，替代关系在 `R/baseline_replacements.json`，最终只纳入每文件一个成功范围。初版 PyROOT reducer 因 JIT memory 失败改为只读必要 branches 的 uproot-v5，版本与替换记录在 `R/reducer_update_v5.json`；未改变 reconstruction。仍在运行的作业没有因慢或 idle 而重复提交。

truth均值控制全量由file0控制与extension clusters `1583116`、`1583124`、`1583126` 组成；prior1、prior100的files1–19分别为 `1583118`、`1583120`，加此前成功的file0。native closure编译/运行cluster `1583132` 成功。第一文件最终图集cluster `1583133`；全量图集最终重跑cluster `1583135`。前一全量绘图 `1583134` 在ROOT TFileMerger合并输入canvas时exit134，统计summary已经写完；改为只读取并累加27个TH1、核对逐文件Entries及102个bins后重导出，不重复重建。第一文件早期v8也出现单输入TFileMerger的空输出，终止失败后以v10成功导出；两类故障均属于绘图汇总，不用失败目录中的PDF作为最终产物。

最终分析复现命令：

```bash
python3 analysis/latest_momentum_study.py \
  --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1 \
  --output /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_momentum/20261116_v1/reproduce_new_directory \
  --groups baseline truth_qop truth_mean qop_prior_100 qop_prior_1
```

通过相同批量 analysis task 入口运行，不在交互节点做全量归约。小型回归 `python3 -m unittest discover -s tests -v` 已通过102项，日志 `R/unittest_final_v2.log`；本任务覆盖 file/entry identity、被动状态计账、两个S2 sensor的区分、非正定 covariance 检出、原 seed gate 不含 phi、原公式除零行为和校准文件不得进入独立gate验证。

## 下一步计划与生产修改边界

1. seed quality 评估同时报告无seed-p gate、p-gate 条件分布与同 seed/track majority 子集，避免共享测量的条件选择冒充全体标定偏差。
2. 以真实 Calypso propagation-aware seed mean/covariance 为下一轮生产候选；窄先验必须在未按 truth 选中的独立数据上检验效率、假轨迹、tail、charge、pull coverage。当前扫描只证明机制，不提供通用最优 inflation。
3. 对无 seed-p gate 时剩余约 −2%的中心、真实局部 electron 能损与强辐射 continuation 做独立检验；精确局部 G4 momentum 仍未在本次 RDO 中直接测量。
4. 本轮没有将固定 field wrapper、measurement variance 或 production inflation 修改为默认修复。任何后续修改先固定候选/hits/truth，再用完整无 truth-gate 的样本验收，另行提交 Calypso patch。
