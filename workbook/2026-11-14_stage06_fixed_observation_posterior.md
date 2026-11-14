# S06：共同观测下的原生有限先验 likelihood 与 posterior

**固定同一实际 cluster 观测后，两个版本各有1,297/1,760个先验分量完成整套观测；16/55个有先验的事件可以对全部32个分量归一化。它们的q/p posterior均值最大版本差为8.748×10⁻¹² GeV⁻¹、最大标准化差6.383×10⁻⁸。其余39个事件的463个缺失分量继续保留，不归一化幸存子集。有效分量数约1–1.46，分量数/宽度收敛尚未测量，S06 Q5继续开放。**

## 阶段目标

承接[内部CKF步长控制](2026-11-14_stage06_ckf_step_control.md)和[原生材料噪声记录](2026-11-14_stage06_material_noise.md)，将原先逐条首选candidate的版本比较推进到共同观测下的有限先验likelihood加权。验收包括所有输入分量的路径覆盖、同一measurement的数值与误差、归一化项、分量内部和分量之间的covariance，以及未完成路径和refit fallback。

结果目录为`P=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_posterior_v1`。权威分析为`P/analysis_v2`、资源与六维比较为`P/auxiliary_v2`，最终图为`P/figures_v3/results`。

## 工作内容与实际配置

使用100044、100116、100117各自原共同先验诊断输入的第一文件、entries `[0,20)`，共60个源事件；具体RDO、材料图、conditions、mask及完整seed component文件保存在`P/tasks/*.json`和每个输出的`source_task.json`。20/17/18个事件有原生seed，其余0/3/2个事件完整保留。事件身份使用输入文件和entry。

本次复用`matched_n32_f256_v1`：每个有先验的事件32个等质量分量，分量内部C为原C/256，均值云与内部C一起保持原F-frame的一阶、二阶矩。before组件均值和完整C保持原样；after的均值使用实际物理旋转，局部C使用既有解析Jacobian。有限宽度Gaussian的非线性坐标变换仍是一阶近似，不能由总mean/C相同推出完整有限分布相同。

实际库及源码为：

| 版本 | 实际库目录 | 实际Calypso源码 |
| --- | --- | --- |
| beforeGeo | `common_prior_ckf_chart_v1/ckf_step_control_v1/overlay_beforeGeo/lib` | `/eos/home-x/xcheng/FASER/geoRotation_calypso_ckf_step_beforeGeo`，master `1eb885ae` 加冻结诊断patch |
| afterGeo | `common_prior_ckf_chart_v1/smoothing_operations_v1/overlay_afterGeo/lib` | `/eos/home-x/xcheng/FASER/geoRotation_calypso_smoothing_audit_afterGeo`，geoRotation `88c457c1` 加冻结诊断patch |

表中相对路径均在`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`下。每个运行实际加载的四个ACTS/observer库都与声明overlay逐项一致，见`P/default_gate_v2/validation.json`和`P/analysis_v2/summary.json`中的runtime。

全部运行开启isolated smoothing，CKF内部step上限10 mm；before的CKF update/smoothing和refit使用既有regular angular chart；refit diagonal scale=1，后续原生完整C×10保留。转换、refit、失败后的converted fallback仍为Calypso原生实现。本次没有修改C++源码，也没有重新编译。

新增的Python CLI只设置既有Gaudi属性：`--ckf-chi2-max`和`--ckf-nmax`，同时记录请求值与实际值。默认重放不设置二者，实际为chi2Max=25、nMax=10；fixed-observation控制显式设置chi2Max=10³⁰、nMax=1，并读取同一份`SourceLinkSelectionFile`。这是观测路径控制，不是生产cuts。

### 共同观测的定义

每个源事件从原`common_prior_ckf_v1/smoke/afterGeo/<DSID>/audit`中取一条原生selected track，按measurement数降序、chi2升序、candidate ID升序固定。该组样本中，有seed的事件都恰好有一条reference selected。其cluster ID集合在有限分量拟合前冻结，`P/observations.json`同时保留实际reverse linked顺序及forward顺序；selector文件为`P/observations/<DSID>_hits.txt`。

有限分量的拟合结果不参与观测选择。reference track可以包含混合truth hits，本项是坐标及观测一致性控制，不是pure-hit误判率。无seed事件保留空观测。

分析从`acts_previous_index`逐边检查实际状态链，反转audit的reverse traversal得到物理forward顺序。只有恰好一条selected路径、完整有序measurement IDs相同、无outlier且innovation有效的分量才具有本项完整likelihood。多候选不事后挑首选；短路径、不同序列、未选中及数值失败均保留为缺项。所有完整路径的measurement值、R和projector在所有分量、两个版本间逐值完全相同，首个sensor的实际surface在物理旋转后通过2×10⁻⁹ mm容差核对。

### Likelihood与有限混合

每个measurement的原生记录提供`predicted`、`predicted_cov`、`measurement`、`measurement_cov`和`projector`。只做这些真实输出的统计归约：

```text
r = y − H m_pred
S = H C_pred Hᵀ + R
log L_i = −1/2 Σ_k [r_kᵀ S_k⁻¹ r_k + log det S_k + d_k log(2π)]
```

measurement坐标单位mm，R和S为mm²。两个版本的同一观测使用相同单位与维数。使用Cholesky solve计算quadratic form和logdet；非正定innovation不加floor补成有效likelihood。没有只使用`exp(−chi2/2)`，也没有使用native累计float chi2替代逐measurement的完整密度。

分量先验质量为1/32；全部32个分量都有完整likelihood时，以log-sum-exp计算event evidence和posterior weights。后验位置/角度/q/p/time的mean及完整六维C在同一首sensor、共同regular chart中归约：

```text
C_total = Σ_i w_i C_i + Σ_i w_i (m_i − m_total)(m_i − m_total)ᵀ
```

保存分量内部C和分量之间的spread，未将32个相关分量视为32个独立事件。q/p边际还给出native Gaussian混合的wrong-charge probability，它不是生产misID率。native refit重新使用相同hits并inflate C，其weighted output单列，未额外乘一次refit likelihood。

## 详细数据：完整分母与缺失路径

以下每行分别对beforeGeo、afterGeo成立；两版计数相同。

| DSID | 输入分量 | selected | 完整观测分量 | 不同/短观测路径 | 未selected | terminal-only | refit成功 / fallback | 完整posterior事件 / 有先验事件 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 100044 | 640 | 614 | 460 | 154 | 26 | 15 | 454 / 160 | 4 / 20 |
| 100116 | 544 | 543 | 390 | 153 | 1 | 0 | 512 / 31 | 5 / 17 |
| 100117 | 576 | 576 | 447 | 129 | 0 | 0 | 425 / 151 | 7 / 18 |
| 合计 | 1,760 | 1,733 | 1,297 | 436 | 27 | 15 | 1,391 / 342 | 16 / 55 |

terminal-only是不同/短路径的子类，不另加进缺失总数。refit成功/fallback分母为全部selected，包含观测不完整路径。两版总共保存3,520条输入分量记录和120条源事件记录；五个无先验事件在每个版本各保留一次。

对39个不完整事件，`posterior=None`，逐事件记录`missing_prior_mass`。即使一个事件31/32个分量成功，也没有用31个重新归一化代替原先验。未完成分量的likelihood是未测量，不指定为零。

![Complete coverage](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_posterior_v1/figures_v3/results/posterior_coverage.png)

## 详细数据：完整分量与完整事件的版本差

所有完整观测分量的`|ΔlogL|`最大值如下；本表不能据以归一化含缺项的事件。

| DSID | 完整分量对 | median | p95 | maximum |
| --- | ---: | ---: | ---: | ---: |
| 100044 | 460 | 3.613×10⁻⁹ | 4.165×10⁻⁷ | 3.446×10⁻⁶ |
| 100116 | 390 | 6.304×10⁻⁸ | 1.259×10⁻⁶ | 4.342×10⁻⁶ |
| 100117 | 447 | 5.590×10⁻⁸ | 1.300×10⁻⁶ | 6.842×10⁻⁶ |

![Native likelihood differences](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_posterior_v1/figures_v3/results/posterior_likelihood.png)

只对全部分量完成的事件计算以下完整有限混合比较。q/p取同一首sensor；σ为两版posterior variance平均值的平方根。TV为同一分量ID的权重total variation。

| DSID | 完整事件对 | 最大q/p均值绝对差 / GeV⁻¹ | 最大q/p均值差/σ | 最大log evidence绝对差 | max weight TV | 最大wrong-charge probability差 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 100044 | 4 | 9.087×10⁻¹³ | 4.340×10⁻⁹ | 4.887×10⁻⁸ | 2.961×10⁻⁹ | 1.752×10⁻¹⁴ |
| 100116 | 5 | 3.530×10⁻¹² | 3.256×10⁻⁸ | 8.751×10⁻⁸ | 4.878×10⁻¹⁶ | 7.424×10⁻¹²⁶ |
| 100117 | 7 | 8.748×10⁻¹² | 6.383×10⁻⁸ | 2.703×10⁻⁷ | 7.380×10⁻¹⁰ | 4.149×10⁻¹¹ |

六维mean最大逐维标准化差依次为2.367×10⁻⁸、3.256×10⁻⁸、6.383×10⁻⁸；完整六维C最大scaled差依次为1.735×10⁻⁷、2.206×10⁻⁸、7.190×10⁻⁸。scaled采用两版对应方差的平均尺度，与原逐state C诊断一致。逐事件完整矩阵、分量权重和原生Gaussian保存在`analysis_v2`；六维差异在`auxiliary_v2/full_posterior_comparisons.json`。

![Complete finite mixtures](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_posterior_v1/figures_v3/results/posterior_comparison.png)

完整事件的afterGeo有效分量数`1/Σw²`：100044为1.000175–1.457724，100116为1–1.000000008，100117为1–1.177956。大量事件几乎只有一个分量承载后验质量。这个结果说明当前有限云的采样收敛必须继续检验，不能将相近的before/after结果当成对原始完整Gaussian先验的收敛证明。

![Signed qop marginal densities](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_posterior_v1/figures_v3/results/posterior_qop.png)

## 验证、执行与产物

默认六组重放与归档10 mm CKF控制逐条比较：58,160条完整audit、26,464条smoother operation记录及348个observer数组字段完全一致。每个输入分量的完整mean/C/surface、native terminal marker、实际加载库、component refit及chart配置通过既有validator；原生measurement q/p update identity同时通过。

首次gate分析cluster1582796失败，原因是误把C++属性默认chi2Max=15当成实际运行值；本批`CKF2Config`实际覆盖为25。修正运行值检查后，将分析冻结到新`analysis_code_v2`，cluster1582797通过；没有重做或更改任何tracking输出。首次准备阶段的失败文件保留在相邻`common_prior_posterior_v1_prepare_failed`，其对应第一次before提交未创建Condor job。

12个native jobs为before cluster1582794、after1582795，各六组，全部exit0；完整归约1582798、绘图/资源1582799成功。每次提交均通过`scripts/submit_momentum.py`重新查询schedd广告与所选队列，选择兼容且running+idle最低的bigbird20，排除bigbird21。各任务申请1 CPU，native 3,500 MB/4 GB scratch，分析3,500 MB/3 GB，绘图2,500 MB/3 GB，线程数=1；实际分配以Condor日志为准。

12组native Athena CPU合计205.14 s，worker CPU合计275.42 s，最大worker RSS=1,031,616 KiB；均从真实`tracking_resource.txt`和`resource.txt`累计，见`P/auxiliary_v2/resource_summary.json`。这些数值不包含独立分析/导图作业。

八项有针对性的测试覆盖measurement likelihood的determinant/常数项、log-sum-exp、小权重稳定性、between-component C、缺失分量不归一化、多候选不挑首条、观测顺序和非SPD innovation不自动修补，全部通过：

```bash
python3 -m unittest discover -s tests -p test_charge_posterior.py -v
python3 -m py_compile analysis/charge_posterior.py analysis/charge_posterior_figures.py analysis/charge_posterior_artifacts.py scripts/prepare_charge_posterior.py configs/slope_tracking.py
```

可执行准备与分析入口为：

```bash
python3 scripts/prepare_charge_posterior.py --root <stage06-root> --phase runs
python3 scripts/prepare_charge_posterior.py --root <stage06-root> --phase analysis --analysis-version <new-version>
python3 scripts/prepare_charge_posterior.py --root <stage06-root> --phase figures --analysis-version <analysis-version> --figure-version <new-version>
python3 scripts/submit_momentum.py --tasks <frozen-tasks.json> --library <matching-overlay/lib> --name <unique-name> --stage S06-posterior --cpus 1 --memory-mb 3500 --disk-gb 3 --flavour microcentury --exclude-schedd bigbird21.cern.ch
```

准备在新输出目录执行；default gate通过后才提交归约，plot完成后再执行validate和报告汇编。原生重建、matching、likelihood归约及图表均由冻结任务保存，不长期保留xAOD或完整ROOT。

最终产物包括[四页新增PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_posterior_v1/figures_v3/results/posterior_report.pdf)、[图索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_posterior_v1/figures_v3/results/figure_index.json)、[完整分母与逐分量记录](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_posterior_v1/analysis_v2/summary.json)。

导图初次独立验证cluster1582800失败：不同CPU上的`np.exp`重新计算analytic density存在末位浮点差。没有放宽原生记录、事件结果、likelihood或横坐标的逐值检查；只对重新计算的analytic density采用`rtol=5×10⁻¹³, atol=10⁻³⁰⁰`。修正q/p刻度重叠及floor标签后的最终绘图cluster1582801、冻结的独立validator v4 cluster1582802均exit0；实测最大density相对差为4.352×10⁻¹⁶、绝对差为4.547×10⁻¹³。四张PNG完整性及逐图PDF与四页汇编的normalized vector content和page boxes检查通过。comparison、qop以原分辨率目视检查；coverage、likelihood与已目视检查的v2 PNG逐byte相同，审查证据为`P/figures_v3/results/visual_review.json`，验证为同目录`validation.json`。

报告汇编cluster1582804 exit0，[S06完整汇报PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v21/charge_stage06_report.pdf)共146页。前142页的index、decoded vector streams和page boxes与report_v20精确一致，新增上述四页；继承检查见`report_v21/lineage_validation.json`，目视审查继承记录见[visual_review.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v21/visual_review.json)。

## 阶段性结论

事实：在共同实际观测、相同有限物理分量和10 mm内部步长控制下，完整分量likelihood及16个完整事件的有限混合posterior版本差很小。默认路径完整重放一致。放宽原生selector后仍有463个分量缺少完整观测，故本次没有覆盖其余39个事件的完整posterior。

假设与待查：缺失路径可能涉及导航、surface接受度或fitter终止；当前仅凭短序列/未selected不足以分配这些原因。小ESS提示有限云的后验分辨率不足，width/分量数扫描需要在覆盖补齐后测量。原始有限Gaussian的非线性坐标变换误差同样需要收敛检查。

未测量：全部55个有先验事件的完整likelihood posterior、当前固定观测下的其他三组width/分量数、更多分量的收敛、这些控制对生产59例错号的净改善。原生refit的额外数据复用/inflation继续作为独立行为分析，未将其命名为同一Bayesian posterior。

## 下一步计划

1. 在真实Calypso中对同一源cluster集合调用既有fixed-source-link `KalmanFitter`，保留实际geometry navigator、material、calibrator和regular chart；先做default OFF完整重放，再量化全部分量的观测覆盖，判别finder/selector与导航终止的作用。
2. 对仍缺失的分量记录最早缺少的surface、原生返回码、目标与边界条件，固定完整mean/C/surface继续最小控制。禁止将这些分量直接设为零likelihood或只对幸存者归一化。
3. 覆盖补齐后，扩展既有四组width/分量数，对event evidence、ESS、完整六维moment、q/p边际与charge probability检查收敛；必要时增加native分量数。达到共同有限分布及完整posterior的验收后，再完成Q5归因与独立生产回归。
