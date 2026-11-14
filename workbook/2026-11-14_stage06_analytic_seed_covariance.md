# S06：只替换 seed covariance 坐标 Jacobian 的实际重建

**解析 Jacobian 将输入物理C差异中位数从约4×10⁻⁹降至约10⁻¹⁴，并显著减小部分 smoothed 差异；4,056个选中分量的有序hits、refit成功状态和最终电荷符号全部不变。100117最差smoothed direction-C仍约1.757×10⁻³，最终|Δq/p|最大仍约9.115×10⁻⁶ GeV⁻¹。因此输入数值微分解释部分差异，不能解释剩余全部差异。**

后续[原生 smoother 逐步审计](2026-11-14_stage06_smoothing_operations.md)已完成：默认和observer各六组输出逐值重放，13,232对实际后向运算均未触发PSD修正。剩余差异已在raw covariance及输入中存在，继续定位前向transport与后向继承。

## 阶段目标

延续[四组宽度扫描](2026-11-14_stage06_ckf_chart_width_scan.md)，检验其剩余版本间差异是否由输入准备中固定步长central difference引起。固定所有物理中心、before C、surface、有限分量设置和实际Calypso配置，只替换after输入的covariance Jacobian。

结果根目录 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`；本次 `A=C/common_prior_ckf_chart_v1/analytic_covariance_v1`，参考 `W=C/common_prior_ckf_chart_v1/width_scan_v1`。

## 工作内容与输入

100044、100116、100117各首个RDO文件entries `[0,20)`，共60个不同源事件；四组仍为n32 C/256、n32 C/1024、n64 C/1024、n64 C/4096。共10,560个相关输入分量，不是独立事件或生产误判率分母。

复用W的12个before运行，新增12个after实际Calypso重建。before为master `1eb885ae22299fb39bb352bead6005d41b0f953e`加已冻结0035–0038补丁的 `common_prior_ckf_chart_v1/overlay_beforeGeo/lib`；after为geoRotation `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`加已记录诊断的 `common_prior_component_refit_v1/overlay_afterGeo/lib`。本次没有修改或重编tracking源码。before的CKF/refit联合regular chart、after原生chart、isolated smoothing ON、refit额外对角scale=1及完整C×10均沿用W。

[prepare_charge_component_tasks.py](../scripts/prepare_charge_component_tasks.py)新增可选 `--covariance-jacobian analytic`，默认仍为历史 `finite_difference`。本次[准备程序](../scripts/prepare_charge_analytic_covariance.py)直接读取W实际运行引用的冻结文本，保留全部非C tokens逐字相同，只替换after的36个C字段。没有重新抽样或重新生成中心。每个before输入C均验证为对角；ACTS对负theta的等价角规范化因此不改变其物理C。任意非零相关矩阵的坐标变换另有独立回归检查。

若R为F→A旋转，T为surface transform，局部位置Jacobian为 `(T_A⁻¹)[:2,:3] R T_F[:3,:2]`；方向2×2块为 `∂(phi_A,theta_A)/∂n_A · R · ∂n_F/∂(phi_F,theta_F)`；q/p和time块严格为1。after C=`J C_F Jᵀ`，仅作对称化和正定性检查。它是中心处的解析微分，不将有限Gaussian的非线性pushforward当作另一个精确Gaussian。

## 详细数据：实际注入的物理C

将真实seed audit的bound mean/C表达为共同FASER frame的position、unit direction、q/p、time。指标为 `max |Cb_ij−Ca_ij|/sqrt(D_i D_j)`，`D_i=(|Cb_ii|+|Ca_ii|)/2`，无量纲。n32 C/256的完整8D物理C结果为：

| DSID | 输入分量 | 中位数：central difference→analytic | 最大值：central difference→analytic |
| --- | ---: | --- | --- |
| 100044 | 640 | 4.120×10⁻⁹ → 2.249×10⁻¹⁴ | 1.119×10⁻⁷ → 3.426×10⁻¹¹ |
| 100116 | 544 | 3.968×10⁻⁹ → 1.072×10⁻¹⁴ | 2.013×10⁻⁷ → 6.056×10⁻¹¹ |
| 100117 | 576 | 4.865×10⁻⁹ → 8.280×10⁻¹⁵ | 2.375×10⁻⁷ → 5.875×10⁻¹¹ |

全部四组的旧最大值为1.6665×10⁻⁶，新最大值6.0561×10⁻¹¹；新中位数范围8.082×10⁻¹⁵至2.249×10⁻¹⁴。完整12组的max/p95/median与逐分量指标见 `A/comparison_v1/summary.json`、`seed_input_metrics.jsonl.gz`。实际after seed记录除C及从C复制的 `variance_qop` 外全部字段逐值一致；两者的 `variance_qop==C[4,4]` 另行断言。

![Seed physical covariance](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analytic_covariance_v1/figures_v1/results/analytic_covariance_seed.png)

## 详细数据：完整链路与相同states

12组after的central-difference/analytic输出均保持相同selected成员和有序hits、相同refit/fallback状态、相同CKF/conversion/final sign。与固定before比较，仍为4,056个共同selected、19,440个配对measurement states、368个双方成功refit和3,688个双方fallback。两版最终错号分量各640；各stage版本间sign分歧均0。没有以成功子组替换全部分母。

下表为n32 C/256；每格为同一完整ordered-hit集合中smoothed direction-C指标的旧→新结果，q/p为绝对版本间差异。

| DSID | direction-C中位数 | direction-C最大值 | smoothed最大 abs(Δq/p) / GeV⁻¹ |
| --- | --- | --- | --- |
| 100044 | 5.185×10⁻⁹ → 1.729×10⁻¹⁰ | 6.282×10⁻⁴ → 5.539×10⁻⁶ | 5.336×10⁻⁸ → 1.427×10⁻⁹ |
| 100116 | 5.646×10⁻⁹ → 9.440×10⁻¹⁰ | 6.845×10⁻⁷ → 3.554×10⁻⁸ | 9.165×10⁻¹⁰ → 8.100×10⁻¹² |
| 100117 | 8.015×10⁻⁹ → 4.800×10⁻⁹ | 1.757226×10⁻³ → 1.757227×10⁻³ | 2.630×10⁻⁸ → 2.105×10⁻⁸ |

全部设置的最终最大|Δq/p|为9.115269031×10⁻⁶→9.115269031×10⁻⁶ GeV⁻¹，主要极端差异仍存在。所有P/F/S指标由旧、新same-hit集合的明确交集计算；本批交集保留全部19,440个states。没有删除outliers或用中位数代表最大值。

![Actual chain outcomes](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analytic_covariance_v1/figures_v1/results/analytic_covariance_outcomes.png)

![Smoothed physical differences](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analytic_covariance_v1/figures_v1/results/analytic_covariance_states.png)

100117的最差例仍为n32 C/256、entry1/component24、hit10183852618321756160、path73.00156 mm。该位置P/F/S的方向差norm分别1.556×10⁻¹³、4.138×10⁻¹²、1.505×10⁻⁶；direction-C分别5.365×10⁻⁹、8.010×10⁻⁹、1.757×10⁻³。这里S比P/F大，支持继续追踪真实后向递推；仅由该位置的三项快照还不能确定哪个下游操作首先引入差异。完整native states、C、measurement与Jacobian在 `A/analysis_v1/extreme_native_states.json`。

## 执行、验证与产物

12次重建cluster1582690全部exit0，CPU user+system合计168.75 s，最大RSS1,034,216 KiB；每job申请1 CPU/3500MB/3GB scratch。统计、对照、导图及验证的顺序任务1582691全部exit0；申请1 CPU/5000MB/3GB。运行实际库、完整seed mean/C/surface及20-event覆盖由24组audit核对，框架仍为Athena24.0.41。

六项focused tests通过，包括相关6×6 C、独立surface局部basis/origin、phi分支边界、负theta/等价角、theta=10⁻⁸近极点、q/p/time恒等块、默认数值路径不变及只改after C。单位方向C有受约束径向零模，极近pole测试以direction block的trace作为向量方差尺度，避免用单个接近零的对角项放大纯表示roundoff；真实图仍保留原有逐分量指标，未改变定义。

三张PNG按原始分辨率逐张检查通过；绘图独立从逐seed及共同state记录重算全部指标，从配对分量重算sign/hit计数。单页/多页PDF的归一化vector content、page boxes、rotation一致。汇编1582693 exit0，127页PDF的前124页index、decoded vector content及page boxes与上一版精确相同，新增3页继承目视记录。

提交前均重新查询兼容schedd负载和队列，选择bigbird20，排除bigbird21；UTC、选择依据、cluster和资源见各worker的 `submission.json`。冻结任务为 `A/afterGeo_tasks.json`、`analysis_figure_sequence_tasks_v1.json`、`report_tasks_v1.json`。统计入口为 [charge_analytic_covariance.py](../analysis/charge_analytic_covariance.py)，绘图为 [charge_analytic_covariance_figures.py](../analysis/charge_analytic_covariance_figures.py)。

- [权威比较统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analytic_covariance_v1/comparison_v1/summary.json)、[原生完整链路统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analytic_covariance_v1/analysis_v1/summary.json)、[实际输入与运行验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analytic_covariance_v1/analysis_v1/validation.json)。
- [三页英文图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analytic_covariance_v1/figures_v1/results/analytic_covariance_report.pdf)、[绘图验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analytic_covariance_v1/figures_v1/results/validation.json)、[目视记录](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analytic_covariance_v1/figures_v1/results/visual_review.json)。
- [127页S06图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v15/charge_stage06_report.pdf)、[前124页衔接检查](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v15/lineage_validation.json)。

## 阶段性结论与下一步

事实：输入数值Jacobian确实引入可测的物理C差异；只替换其解析微分，能降低部分实际smoothed差异。所有选中hits/refit/sign不变，μ−大差异及最终q/p极端差异保留。有限Gaussian内的非线性分布映射、likelihood加权posterior聚合及完整收敛仍未验收；Q5继续。

源码提供下一项可判别实验：ACTS32.0.2的[GainMatrixSmoother.cpp](https://github.com/acts-project/acts/blob/v32.0.2/Core/src/TrackFitting/GainMatrixSmoother.cpp)先显式求predicted C的逆，再对非正定smoothed C调用修正。实际安装的 `GainMatrixSmoother.hpp`、`covariance_helper.hpp` 与该tag逐字节相同，下载来源及原文在 `A/source_inspection/source_validation.json`。这些是运算入口证据，尚未实测剩余差异由求逆或PSD修正导致。

下一步在真实Calypso后向执行顺序中记录gain、regular-chart C/J、更新前后mean/C、PSD修正前后矩阵和condition，先做observer逐值不改结果的gate，再追踪首个差异增长操作。需要实际单因素重建来区分输入精度、数值病态和非线性chart贡献，不能把当前最大值直接命名为Calypso缺陷或生产误判率问题。
