# S06：原生 fixed-hit fitter 的观测覆盖与导航缺项

**相同60个源事件、共同有限先验和冻结cluster观测下，实际Calypso KalmanFitter每版完成1,311/1,760条路径，比CKF finder增加14条；完整posterior仍为16/55个有先验事件。其余449个分量保留为未测量likelihood，未归一化幸存子集。27个原生失败均为达到5,000 steps上限；另外422条成功返回的路径缺少部分观测。源码和定向日志提示external measurement surface的pre/post-step边界检查不一致，单因素因果验证尚未完成。S06 Q5继续开放。**

## 阶段目标

承接[共同观测likelihood实验](2026-11-14_stage06_fixed_observation_posterior.md)，检验缺失的463个分量能否通过Calypso既有fixed-source-link fitter补齐，区分组合finder/selector与导航限制。要求默认关闭精确重放、逐分量核对实际传入的hits、保留所有失败与短路径，并为下一项控制提供实际native返回值和导航证据。

结果目录为`U=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_fixed_fitter_v2`。权威归约为`U/analysis_v1`，资源为`U/auxiliary_v1`，图为`U/figures_v1/results`；定向重放为`N=U/navigation_log_v1`。

## 工作内容与冻结配置

输入为100044、100116、100117各自第一RDO文件的entries `[0,20)`，共60个源事件。20/17/18个事件有原生seed，五个无seed事件保留。完整RDO路径、material JSON、conditions、mask与选定source entry见`U/tasks/*.json`，沿用前项的同一输入。事件身份为输入文件和entry。

先验仍为`matched_n32_f256_v1`：每个有先验事件32个等质量Gaussian分量，内部C为原C/256。before完整mean/C/surface不变，after沿用解析seed covariance Jacobian；有限Gaussian非线性映射仍是一阶近似。冻结观测文件与前项逐byte相同，拟合输出没有参与观测选择。

新增默认关闭的`ComponentFixedHitFitter`只在有限分量诊断中启用，从**已经过滤的`slAccessor.container`**构造输入，调用实际`TrackFitterFunction`/ACTS `KalmanFitter`。普通geometry Navigator、实际磁场/material、calibrator、measurement update、smoother与后续原生conversion/refit均保留；未使用独立tracking实现或DirectNavigator。reverse-filtering阈值设0，isolated smoothing开启，内部propagator步长上限10 mm。before使用既有regular update/smoothing chart，refit diagonal scale为1、后续完整C×10保留。配置中的selector仍为chi2Max=10³⁰、nMax=1，但fixed fitter不调用CKF measurement selector。

| 版本 | 实际Calypso源码 | 实际加载库 |
| --- | --- | --- |
| beforeGeo | `/eos/home-x/xcheng/FASER/geoRotation_calypso_fixed_fitter_v2_beforeGeo`，master `1eb885ae`加冻结诊断patch | `U/overlay_beforeGeo/lib` |
| afterGeo | `/eos/home-x/xcheng/FASER/geoRotation_calypso_fixed_fitter_v2_afterGeo`，geoRotation `88c457c1`加冻结诊断patch | `U/overlay_afterGeo/lib` |

相对既有CKF-step before与smoother-audit after基底的当前增量为[0046 before](../patches/0046_native_component_fixed_fitter_beforeGeo.patch)和[0046 after](../patches/0046_native_component_fixed_fitter_afterGeo.patch)。完整tracked diff、added headers和reverse-apply验证在`U/source_validation.json`；各构建的`software.json`记录实际源码。默认路径使用新库而不启用fitter开关，逐字段与前项fixed-observation CKF运行比较。

每个ON输入分量新增`component_fixed_observations`审计，保存实际count及hit IDs。六组ON运行的3,520个分量均已验证：实际输入集合与冻结观测完全相同、IDs无重复、每个分量恰好一次记录。完整seed mean/C/surface、实际加载库、chart/refit属性和measurement q/p更新恒等式继续通过既有validator。

### 前一实现错误与更正

`common_prior_fixed_fitter_v1`错误地从未过滤的原始`sourceLinks`构造fitter输入，未实现冻结观测实验。其ON覆盖率、likelihood和posterior无效；原目录及归约目录均有`RESULTS_INVALID_FOR_FIXED_OBSERVATION.md`，没有接受其图或加入汇报。0045仅保留为这一错误实现的历史记录，0046不以0045为应用基底。

v1的OFF重放证据仍单独有效；本报告全部ON结果来自新worktree、新构建和新输出v2。实际hit集合审计是本次增加的运行门槛，避免仅凭CLI路径推断fitter使用了正确输入。

## 详细数据：完整分母和solver迁移

下表分别对beforeGeo、afterGeo成立，两版计数相同。refit的分母为全部selected，不仅限于完整观测路径。

| DSID | 输入分量 | selected | CKF完整路径 | fitter完整路径 | fitter短/不同路径 | 未selected | terminal-only | refit成功/fallback | 完整posterior/有先验事件 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 100044 | 640 | 614 | 460 | 468 | 146 | 26 | 15 | 458/156 | 4/20 |
| 100116 | 544 | 543 | 390 | 394 | 149 | 1 | 0 | 513/30 | 5/17 |
| 100117 | 576 | 576 | 447 | 449 | 127 | 0 | 0 | 427/149 | 7/18 |
| 合计 | 1,760 | 1,733 | 1,297 | 1,311 | 422 | 27 | 15 | 1,398/335 | 16/55 |

terminal-only是短路径的子类，不另加到缺项总数。CKF→fitter的完整→不完整迁移为0；不完整→完整为8、4、2，共14。没有新增完整事件，39个事件的posterior仍为`None`，449个缺项分量的likelihood未测量，未指定为0。

![Native fitter coverage](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_fixed_fitter_v2/figures_v1/results/posterior_solver_coverage.png)

所有实际measurement值、R、projector和首个sensor surface继续通过跨分量/版本检查。完整分量的最大`|ΔlogL|`分别为3.446×10⁻⁶、4.342×10⁻⁶、6.842×10⁻⁶。16个完整事件的最大q/p posterior mean版本差仍为8.748×10⁻¹² GeV⁻¹、最大差/σ为6.383×10⁻⁸。图中相近的版本结果不补足缺项，也不证明原始Gaussian先验的采样收敛；当前ESS约1–1.46。

![Native likelihood comparison](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_fixed_fitter_v2/figures_v1/results/posterior_likelihood.png)

likelihood使用实际prediction的`r=y−Hm`和`S=HCHᵀ+R`，保留quadratic、logdet及维数常数项。完整六维mixture保留within/between C；refit重复使用同一观测，其weighted output单列，不另乘第二次likelihood。具体定义沿用前项报告。

## 详细数据：剩余路径和native导航重放

对全部3,520条solver记录归约后的下表，同样分别对每个版本成立。类别按冻结有序观测序列定义，缺前缀/后缀不是物理原因标签。

| DSID | 完整 | 仅缺前缀 | 仅缺后缀 | 中间缺项 | 两端均缺 | 空路径/native错误 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 100044 | 468 | 68 | 22 | 55 | 1 | 26 |
| 100116 | 394 | 80 | 6 | 63 | 0 | 1 |
| 100117 | 449 | 52 | 35 | 38 | 2 | 0 |
| 合计 | 1,311 | 200 | 63 | 156 | 3 | 27 |

27个原生错误都记录为category=`PropagatorError`、message=`Propagation reached the configured maximum number of steps`，没有按整数码猜测类别。其余422个fitter成功返回并不表示已使用完整观测：当前ACTS在已有measurement且navigation break时也可完成smoothing并返回成功。

随后在100044中固定四个分量，分别重放before/after，共8个作业。每个作业仍读取同一文件的20个源事件，只注入所列一个分量；观测文件、输入mean/C/surface、材料和实际库保持原样，只将CKF ACTS logger设为VERBOSE。八个分量共160条native component audit与各自完整运行逐值精确相同，故日志和减少独立分量调用没有改变这些拟合。selected states、conversion/refit及失败terminal均包括在相等检查中。

| 源entry/component | 例子 | 输入/实际measurement数 | 实际记录 |
| --- | --- | ---: | --- |
| 0/0 | 完整路径 | 12/12 | 正常完成 |
| 0/16 | 缺早期hits | 12/8 | vol27 layer2/4的部分surface未进入measurement update |
| 1/4 | 原生失败 | 6/0 | 日志确认配置5,000 steps、实际5,000 steps；末次forward位置z_F约49,786.9 mm |
| 17/0 | terminal-only | 5/1 | 仅vol29/layer2/sensitive16被更新，随后navigator未继续整套观测而进入smoothing |

两版日志在F-frame的上述末点与序列一致。日志位置由ACTS有限小数位打印，只作路径定位；完整精度版本比较使用audit，不用这些打印值声明10⁻¹²级闭合。`N/analysis_v1/native_log_paths.npz`保存实际打印路径，`summary.json`保存native markers、实际geometry IDs、终止返回值和八组逐值检查。

### 可检验的导航机制

匹配的CVMFS ACTS头文件位于`/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/include/Acts/`：

- `TrackFitting/KalmanFitter.hpp:373`插入external measurement surfaces；`:415`允许已有measurement后的navigation break触发finalization。
- `Propagator/Navigator.hpp:643`在`targetSurfaces()`中对external surface关闭boundary check；`:560`的`surfaceStatus()`却仍传入`BoundaryCheck(true)`。
- `Propagator/detail/SteppingHelper.hpp:42`用传入的boundary check计算实际surface intersection。

entry0/component16的日志出现同一surface在post-step被判NOT reachable，随后pre-step忽略边界时已为ON SURFACE，接着切换layer；该surface没有进入measurement update。entry17/component0也出现下一measurement surface未交付后转入navigation break。这与pre/post边界处理不一致相符，但**尚未测量只改变这一检查能恢复多少分量**。27个step-limit错误也不由这个小样本直接归为同一原因。

下一项在真实Calypso内部仅对本次固定输入的measurement surfaces统一pre/post boundary语义，保持普通Navigator、实际surface/material/field及原生Stepper计算，并先做默认关闭精确重放。需重新测量全部1,760分量，而非只展示恢复的个例。

## 执行、验证与产物

两版构建1582810/1582811、12组native重建1582812/1582813均exit0。pipeline1582814顺序执行默认gate、likelihood归约和solver coverage，exit0；资源1582815、绘图1582816、独立验证1582819均exit0。默认六组187,160条audit、76,166条smoother operation记录、348个observer数组字段精确重放。8项既有posterior测试通过，新增/修改Python文件通过`py_compile`，`git diff --check`通过。

定向导航重放before1582817、after1582818各四个作业均exit0，归约1582820 exit0。准备脚本首轮因未创建tasks父目录失败，未提交任何作业；修正后保留`navigation_log_v1_prepare_failed`并重新冻结成功。未因观察超时重启任何运行。

全部提交使用`scripts/submit_momentum.py`查询当时schedd广告与队列，选择兼容且running+idle最低的bigbird20，排除bigbird21。每次选择依据与资源在`*_worker/submission.json`和schedd记录中。主重建申请1 CPU、3,500 MB/4 GB scratch；pipeline1 CPU、4,000 MB/3 GB；绘图/资源/验证1 CPU、2,500 MB/3 GB；导航重放及归约1 CPU、3,000 MB/3 GB。12组Athena CPU合计313.32 s、worker CPU合计438.60 s、最大worker RSS为1,031,008 KiB；不含独立build、归约与导航作业。

准备、执行入口：

```bash
python3 scripts/prepare_charge_fixed_fitter.py --root <stage06-root> --version v2 --phase source
# 提交并完成build_tasks后，冻结runs、analysis、pipeline、figures，再按依赖提交。
python3 scripts/submit_momentum.py --tasks <U/pipeline_tasks_v1.json> --library <U/overlay_afterGeo/lib> --name <unique-name> --stage S06-fixed-fitter --cpus 1 --memory-mb 4000 --disk-gb 3 --flavour microcentury --exclude-schedd bigbird21.cern.ch
python3 scripts/prepare_charge_fitter_navigation.py --input <U> --phase runs
# 两版replay均成功后：
python3 scripts/prepare_charge_fitter_navigation.py --input <U> --phase analysis
```

准备入口写入新目录，不覆盖历史输出。pipeline与单独gate/analyze/coverage tasks是替代入口，不同时提交。

五张图的输入、likelihood、事件结果和横坐标逐值检查通过；analytic density的独立重算差为0。五张PNG完整，逐图PDF与五页汇编的normalized vector content及page boxes一致，五张图均目视检查英文标题、表格、轴、图例与q/p指数标记。证据为`U/figures_v1/results/validation.json`和`visual_review.json`。

最终产物：

- [新增五页PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_fixed_fitter_v2/figures_v1/results/posterior_report.pdf)与[图索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_fixed_fitter_v2/figures_v1/results/figure_index.json)。
- [全部分量、likelihood和posterior summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_fixed_fitter_v2/analysis_v1/summary.json)、[solver coverage](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_fixed_fitter_v2/analysis_v1/solver_coverage.json)、[导航重放归约](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_fixed_fitter_v2/navigation_log_v1/analysis_v1/summary.json)。
- [151页S06汇报PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v22/charge_stage06_report.pdf)，汇编1582821 exit0。前146页index、decoded streams与boxes精确继承report_v21，新增147–151页；数值、lineage及目视证据在report_v22相应JSON中。

## 阶段性结论

事实：切换到实际native fixed-hit fitter仅补回14个完整分量，没有补齐任何新的完整事件；失败与短路径都不能由放宽selector自动消除。3,520个实际fitter输入集合正确，默认OFF完整重放通过，8个定向日志重放保留原计算值。两版完整有限混合仍相近。

待检验：native Navigator的external-surface边界语义不一致与部分漏测日志相符；其全部发生率及单因素恢复量未测量。step-limit错误的首个失败导航决定仍需定位。

未测量：其余39个事件的完整posterior、当前frozen observations下另外三组width/分量数收敛、边界控制对全量生产59例wrong的净改善。当前有限混合结果不证明完整非线性Gaussian先验等价。

## 下一步计划

1. 增加默认关闭、仅作用于固定measurement surface IDs的native boundary控制；保持原生Stepper算法和geometry Navigator，OFF必须精确重放U的全部输出，再测全部分量的有序观测覆盖、likelihood与终止类别。
2. 对仍缺项的分量，按原生geometry/volume状态、首个遗漏surface及step-limit类别继续最小控制。缺失likelihood保持未测量，不用零或幸存子集归一化代替。
3. 覆盖补齐后运行既有另外三组有限width/分量数，测event evidence、ESS、完整六维moments、q/p边际及charge probability的收敛，继续Q5成因验收与生产回归。
