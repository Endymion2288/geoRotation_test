# S06：真实材料噪声记录与方差计算精度

**相同物理输入下，原生材料噪声在两个坐标表达中的约10⁻⁸差异，来自φ方差和θ方差采用不同的乘法精度。384条真实原生记录的θ方差全部精确符合float平方，φ方差符合double平方至3.96×10⁻¹⁶；仅对记录项统一平方精度的算术诊断，使共同方向噪声tensor的差异降到8×10⁻¹⁶以内。本次没有改变材料公式或重新生成修正后的track covariance，完整posterior继续验收。**

## 阶段目标

延续[完整输入原生transport](2026-11-14_stage06_bound_transport.md)及[内部CKF步长控制](2026-11-14_stage06_ckf_step_control.md)，直接读取实际Calypso材料相互作用计算出的方差，区分原生噪声的坐标差和继承输入/传播差异。新增观察首先必须精确重放既有transport结果及原trace字段。

结果目录 `N=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/material_noise_v1`。权威分析为`N/analysis_v2`，英文图为`N/figures_v1/results`。

## 工作内容与实际软件

实际Calypso源码为：

```text
/eos/home-x/xcheng/FASER/geoRotation_calypso_bound_transport_v1/
  Tracking/Acts/FaserActsKalmanFilter/src/CalypsoSlopeDiagnosticsAlg.cxx
```

[0044补丁](../patches/0044_native_recorded_material_noise.patch)在既有material CSV末尾追加三列，直接写`RecordedMaterial::MaterialInteraction`的`hit.sigmaPhi2`、`hit.sigmaTheta2`、`hit.sigmaQoP2`，分别为rad²、rad²和GeV⁻²。没有改动propagator、interactor、update或smoother。原编译副本与修改前源码精确相同，新实际源码、冻结副本及编译副本逐字节一致，increment reverse-check通过，见`N/auxiliary_v1/source_validation.json`。

仅重编诊断plugin，构建为`N/software_probe`。两版`overlay_beforeGeo` / `overlay_afterGeo`继续使用上一阶段真实geometry、field和extrapolation库；全部运行按实际resolved路径检查加载库。环境为Athena24.0.41、EL9/GCC13、FASERNU-04 / OFLCOND-FASER-05。

这些记录来自已有的`propagationSteps()`调用；bound `propagate()`仍提供原完整mean/C输出。本次没有在内部CKF的每个材料操作中另加observer，也没有把trace的噪声tensor冒称完整bound transport的逐步Q/J实测。

## 输入与覆盖

原seed来自三个DSID各首RDO的entries `[0,20)`，isolated smoothing ON。当前transport重放取其中四条已定位路径：100044 entry8/component6、100116 entry10/component14、100117 entry1/component24及entry8/component4。

每条路径的相邻完整states分别使用各版本原生输入及固定before物理输入，在两种坐标中调用实际Calypso。包括reference显式10 m、100 mm、10 mm步长上限和material OFF；这里的10 m是extrapolator属性，区别于CKF未显式设置的原生默认上限。

| DSID | 每输入模式、每步长的links | paired material interactions |
| --- | ---: | ---: |
| 100044 | 27 | 8 |
| 100116 | 23 | 8 |
| 100117 | 52 | 16 |

24个运行共有816对输入links；材料开启的18组记录612对trace序列，全部一致。其中192对有scattering interaction，另420对没有材料交互，仍保留在覆盖分母。材料OFF的6组原配置未请求trace，其204对links的噪声为**未记录**；完整transport arrays仍通过重放。没有将无trace文件当作零噪声实测。

## 详细数据：观察重放验证

[完整gate](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/material_noise_v1/legacy_gate_v2/validation.json)核对24个运行：

- 576个transport array fields逐值相同。
- 18组material trace的384行全部原字段逐字符串相同，calls CSV逐字节相同。
- 三列新增方差全部finite且非负；输入、材料、options与对应旧任务一致。
- 实际加载库及编译源码副本正确。

首次gate 1582787因脚本误要求材料OFF组也有trace而exit1；没有重跑tracking。修正为按冻结的`--trace-material`选项检查存在性，并将未记录与空trace区分，新gate1582789 exit0。四项定向测试覆盖旋转不变的各向同性tensor、float平方导致的方向变化、零scattering记录，以及未请求trace和缺失必需trace的区别。

## 详细数据：相同物理输入的方向噪声

从native方向角的导数，将记录的`diag(sigmaPhi2,sigmaTheta2)`转换成FASER frame中3×3 unit-direction噪声tensor。比较指标为tensor最大元素差除以两侧平均`varianceTheta`，无量纲。`same-double-square`只用记录项恢复native float θ₀，再将θ方差也以double乘法求平方；不是一次新的track fit。

| DSID | reference原生Q最大差异指标 | same-double-square算术诊断 | 10 mm原生Q最大差异指标 | 10 mm算术诊断 |
| --- | ---: | ---: | ---: | ---: |
| 100044 | 3.001×10⁻⁸ | 4.803×10⁻¹⁶ | 3.001×10⁻⁸ | 7.922×10⁻¹⁶ |
| 100116 | 4.991×10⁻⁸ | 5.495×10⁻¹⁶ | 4.991×10⁻⁸ | 4.081×10⁻¹⁶ |
| 100117 | 3.622×10⁻⁸ | 6.723×10⁻¹⁶ | 3.622×10⁻⁸ | 6.564×10⁻¹⁶ |

相同输入的θ方差、q/p方差在两版逐值相同；差异来自微小非各向同性在物理空间的方向随坐标表达改变。对各自继承的native输入，100117的最大Q指标仍为1.823×10⁻⁶，算术诊断仍为1.809×10⁻⁶；输入方向本身已有差异，该诊断不会将其消除。

![Native material noise coordinates](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/material_noise_v1/figures_v1/results/material_noise_coordinates.png)

## 详细数据：原生精度关系与源码解释

384条recorded interactions为192对×两版，均有scattering。全部记录满足：

| 检查 | 最大绝对relative差 |
| --- | ---: |
| recorded θ方差 vs float(θ₀×θ₀) | 0，逐值精确 |
| recorded φ方差×sin²θ vs double(θ₀)² | 3.9585×10⁻¹⁶ |
| recorded各向同性偏离幅度 | 5.4714×10⁻⁸ |

运行使用的ACTS头文件声明32.0.2，`CommitHash="GITDIR-NOTFOUND-dirty"`，不能从它证明精确runtime commit。对应版本的[官方源码](https://github.com/acts-project/acts/blob/v32.0.2/Core/src/Propagator/detail/PointwiseMaterialInteraction.cpp)为：

```cpp
const float theta0 = computeMultipleScatteringTheta0(...);
const auto sigmaPhi = theta0 * (dir.norm() / VectorHelpers::perp(dir));
variancePhi = sigmaPhi * sigmaPhi;
varianceTheta = theta0 * theta0;
```

`sigmaPhi`经double表达式提升后平方；`theta0*theta0`先按float乘法舍入，再存进double方差。若两者统一精度，方向噪声为切平面上的各向同性tensor；原生float舍入留下约10⁻⁸的不等方差。该源码解释由实际recorded方差的逐值关系支持，不依赖将公开tag当作runtime精确commit，也不使用错误的`1−dirz²`消去猜想。

![Native material variance arithmetic](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/material_noise_v1/figures_v1/results/material_noise_arithmetic.png)

## 执行、验证与产物

实现为[任务冻结](../scripts/prepare_charge_material_noise.py)、[原生记录分析](../analysis/charge_material_noise.py)、[英文图表](../analysis/charge_material_noise_figures.py)及[来源/资源归约](../analysis/charge_material_noise_artifacts.py)。所有批量运行、归约、图表及报告组装经`scripts/submit_momentum.py`提交HTCondor，每次重新查询schedd并按兼容普通bigbird的running+idle最低选择bigbird20，排除bigbird21；查询及选择依据随`submission.json`保存。

构建cluster1582774；before/after各12个运行1582782/1582783；成功gate1582789；权威分析1582790；绘图/资源1582791；图表验证1582792。重建、归约各1CPU，构建4CPU；24个worker CPU合计308.94 s，最大RSS964,304 KiB。native Athena CPU没有单独计时，记为**未测量**。

复现入口，输出使用新的版本目录：

```bash
python3 scripts/prepare_charge_material_noise.py --root <root> \
  --phase analysis-tasks --analysis-version v2
python3 scripts/submit_momentum.py --tasks <material_noise>/analyze_tasks_v2.json \
  --library <material_noise>/overlay_afterGeo/lib --name <unique> \
  --stage S06-material-noise --cpus 1 --memory-mb 2300 --disk-gb 3 \
  --flavour microcentury --exclude-schedd bigbird21.cern.ch
```

主要产物为[完整统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/material_noise_v1/analysis_v2/summary.json)、[覆盖验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/material_noise_v1/analysis_v2/validation.json)、[三页英文图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/material_noise_v1/figures_v1/results/material_noise_report.pdf)。逐记录及trace序列随结果保存。

三张PNG已按原始分辨率目视检查；绘图数值从paired records独立重算一致，逐图与汇编PDF的归一化vector内容、page boxes及rotation一致。`py_compile`、四项定向测试及`git diff --check`通过。[142页S06汇报PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v20/charge_stage06_report.pdf)由cluster1582793生成、exit0；其前139页的索引、decoded vector streams和page boxes精确继承report_v19，新增本次三页，见[继承验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v20/lineage_validation.json)。

## 阶段性结论

事实：原生噪声中的约10⁻⁸坐标差有实际方差精度关系支持；固定输入的算术诊断闭合到double舍入尺度。继承输入差异仍保留。所有原有transport数值与trace字段精确不变。

未测量：修改原生variance计算后完整CKF/refit的输出、完整有限先验的likelihood posterior及生产misID变化。局部噪声精度不是本批59例afterGeo生产错号的已证明根因，本项没有关闭S06 Q5。

## 下一步计划

在真实Calypso中固定共同物理measurement集合，直接从原生predicted/measurement/C计算完整innovation likelihood；保留所有先验分量、失败、不同路径和fallback，不将不同observation序列的似然混成同一posterior。随后执行width/分量数收敛及共同物理分布比较，对剩余可见差异继续定位。
