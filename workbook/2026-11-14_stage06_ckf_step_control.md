# S06：真实 CKF 内部传播步长控制

**完整 Calypso CKF 的 μ− 最大 smoothed direction-C 差异指标从默认的 1.940×10⁻³ 降至 10 mm 上限下的 1.317×10⁻⁷。三个样本的首选存在性、有序 hits、完整 state 链、CKF/final sign 和 refit 成败均不随本批步长控制改变。较大的版本间分叉对内部传播步长敏感；剩余材料噪声与完整有限先验 posterior 继续检查。**

## 阶段目标

延续[完整输入原生 transport 控制](2026-11-14_stage06_bound_transport.md)，直接改变实际 CKF 内部的传播步长上限，检查此前约 134→1194 mm 长段的 mean/C 分叉是否在完整 finder、measurement update、候选独立 smoothing、conversion/refit 中响应。默认关闭的新配置必须重现原有全部输出。

结果目录 `S=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1`。权威分析为 `S/analysis_v1`；三页图为 `S/figures_v2/results`。

## 工作内容与实际软件

- beforeGeo 实际源码为 `/eos/home-x/xcheng/FASER/geoRotation_calypso_ckf_step_beforeGeo`，基底 master `1eb885ae22299fb39bb352bead6005d41b0f953e`，继承此前已验证的原生 component、regular chart、refit 和 smoothing observer 控制。
- 新增 [0043 补丁](../patches/0043_before_native_ckf_step_control.patch)只增加默认值 0 的 `MaxStepSize` 属性；值大于 0 时设置实际 `Acts::PropagatorPlainOptions::maxStepSize`，单位 mm。没有替换 propagator、update 或 smoother。
- beforeGeo 重编 28 个 KalmanFilter translation units，复用已匹配的 geometry 库。构建位于 `S/software_beforeGeo`，运行使用 `S/overlay_beforeGeo`。
- afterGeo 复用 `smoothing_operations_v1/overlay_afterGeo`，其真实 CKF 已有同义属性，本次不重编。
- 新 CLI 为 `configs/slope_tracking.py --ckf-max-step-mm {100,10}`；不传参数保持原生默认，同时记录 requested/actual Gaudi property。这里的默认不等同于前一独立 extrapolator 实验的显式 10 m 上限。refit 和 observer 的传播步长保持原设置。

实际执行环境仍为 Athena 24.0.41、EL9/GCC13、FASERNU-04 / OFLCOND-FASER-05；before/after 的 native geometry、场与正确 row-major 材料设置继承上一阶段。每个运行实际加载的库均按 resolved 路径检查。完整 tracked diff、附加 headers 与反向 apply 检查保存在 `S/source_validation.json` 和 `S/auxiliary_v1/source_validation.json`。

## 输入、选择与覆盖

100044、100116、100117 各首 RDO 的 entries `[0,20)`，共 60 个独立源事件。每个已有原生 seed 使用 n32、within-C/256 的完整物理 mean/C/surface 分量；after seed-C 使用已验收的解析 Jacobian。before 的 CKF/refit 使用 regular angular charts，两版 isolated smoothing 均 ON。refit 对角 scale=1，原有完整 C×10 保留。

每个版本比较 default、100 mm、10 mm，共 18 个真实重建。三种设置之间 seed 输入逐值一致；本次没有重新抽样。按原生 selected candidates 的 measurements 降序、χ²升序、candidate ID 取第一条用于配对；全部 native candidate 数及 finder errors 另存于 runtime checks。没有加 GitHub residual 图的质量 cuts，也没有用 truth charge/p 选择 candidate。

| DSID | 每版本、每设置输入分量 | 首选分量 | 双方成功 refit | 双方 fallback | terminal-only 首选 | 完整配对 states |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 100044 | 640 | 311 | 37 | 274 | 7 | 4,157 |
| 100116 | 544 | 254 | 24 | 230 | 0 | 3,580 |
| 100117 | 576 | 360 | 69 | 291 | 2 | 6,420 |
| 合计 | 1,760 | 925 | 130 | 795 | 9 | 14,157 |

表中计数在三个步长设置中相同。所有版本间比较及同版本 default→100/10 mm 比较均无单侧新增/丢失首选、不同 ordered hits、不同 full-state signature 或 selected/final sign 分歧；invalid sign 为 0。未选中的 835 个输入分量仍在输入分母中。两版均错的诊断分量会保留，这些 cloud counts 不作为生产样本 misID。

21 组比较共有 99,099 对 states，包含 passive/material states 及 terminal-only 路径。真实 linked-state 链按 geometry ID、hit 和 measurement/outlier/hole flags 配对；不是按列表位置将不同拓扑强行配成同一状态。

## 详细数据：默认设置重放

新增属性不设置上限时，六组完整重放通过：

- 58,160 条 decoded audit records 精确相同，含 truth、seed、found、selected、conversion、refit、失败与 component end。
- 26,464 条 decoded smoothing-operation records 精确相同。
- 348 个紧凑 observer array fields 逐值相同。
- 每个源事件、完整输入分量、实际 runtime chart/refit/isolation 和真实库路径另由原生 component validator 检查。

权威检查为 [default gate](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1/default_gate_v1/validation.json)。

## 详细数据：完整链路的版本间差异

方向 mean 使用 FASER frame 的 unit direction。方向 C 转到共同 `(phi_A, latitude_A)` chart，保留完整 6×6 C 及交叉项。表中的 scaled 指标是方向 2×2 block 元素差除以两侧平均 diagonal variance 构成的尺度，再取最大值；它是无量纲矩阵比较量，不是所有矩阵元素共同的相对误差。

| DSID | 设置 | 最大 smoothed unit-direction 差 | 最大 smoothed direction-C 指标 |
| --- | --- | ---: | ---: |
| 100044 | default | 2.038×10⁻⁸ | 7.615×10⁻⁶ |
| 100044 | 100 mm | 2.966×10⁻¹⁰ | 1.103×10⁻⁷ |
| 100044 | 10 mm | 6.563×10⁻¹¹ | 3.875×10⁻⁸ |
| 100116 | default | 1.379×10⁻¹⁰ | 3.552×10⁻⁸ |
| 100116 | 100 mm | 1.379×10⁻¹⁰ | 3.553×10⁻⁸ |
| 100116 | 10 mm | 1.379×10⁻¹⁰ | 3.554×10⁻⁸ |
| 100117 | default | 1.505×10⁻⁶ | 1.940×10⁻³ |
| 100117 | 100 mm | 9.453×10⁻⁷ | 1.064×10⁻³ |
| 100117 | 10 mm | 2.475×10⁻¹⁰ | 1.317×10⁻⁷ |

100117 的 smoothed direction-C 指标中位数从 2.761×10⁻⁹ 变为 2.592×10⁻⁹，p95 从 1.748×10⁻⁷ 变为 2.246×10⁻⁸。步长主要压低极端尾部，没有使所有 covariance 差异归零。各 P/F/S 的 max/p95/median、绝对 C 差与 q/p 差保存在 [summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1/analysis_v1/summary.json)。

首选 CKF 与最终 track 的最大绝对 q/p 差另从完整 audit 独立归约；下表单位为 GeV⁻¹，分母为每组全部共同首选，含 refit fallback。

| DSID | 首选分量数 | selected：default→10 mm | final：default→10 mm |
| --- | ---: | --- | --- |
| 100044 | 311 | 1.42347×10⁻⁹→1.22593×10⁻¹² | 1.31555×10⁻⁶→1.31677×10⁻⁶ |
| 100116 | 254 | 7.87955×10⁻¹²→4.54970×10⁻¹² | 2.35024×10⁻⁸→2.24847×10⁻⁸ |
| 100117 | 360 | 2.09031×10⁻⁸→1.81999×10⁻¹¹ | 9.11527×10⁻⁶→9.11668×10⁻⁶ |

selected 数值差在 μ− 中显著缩小，但 final 最大差没有同时消失；refit 的输入、数值路径及有限分布仍需检验。本次只改变 CKF 内部上限，不能将其称为完整 refit 数值收敛。全部21组 selected/final q/p 比较见 [qop_comparisons.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1/auxiliary_v1/qop_comparisons.json)。

![Native CKF step distributions](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1/figures_v2/results/ckf_step_distributions.png)

主例 100117 entry1/component24 的完整链路图如下。默认与 100 mm 在长段后仍出现较大 predicted mean/C 分叉，10 mm 将该分叉压低；后向 smoothing 随之变化。该结果来自每次重新运行整条真实 CKF 链，区别于上一阶段每段独立从旧归档 state 重启的 extrapolator 控制。

![Native CKF extreme path](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1/figures_v2/results/ckf_step_extreme.png)

## 验证、运行与产物

实现为 [prepare_charge_ckf_step_control.py](../scripts/prepare_charge_ckf_step_control.py)、[charge_ckf_step_control.py](../analysis/charge_ckf_step_control.py)、[图表](../analysis/charge_ckf_step_figures.py)及[最终参数/资源归约](../analysis/charge_ckf_step_artifacts.py)。四项针对性测试覆盖单侧 selection、同 hits 但 passive topology 不同、terminal-only 和 invalid sign 的独立计数。

所有批量重建、归约和导图通过 `scripts/submit_momentum.py` 提交 HTCondor。每次提交重新查询 schedd 广告和队列，选择兼容且 running+idle 最少的 bigbird20，排除已有 EOS 兼容记录的 bigbird21；请求 CPU 与 BLAS/OpenMP 线程数一致。构建申请 4 CPU/8500 MB；18 个重建各 1 CPU/3500 MB；完整分析 1 CPU/6500 MB。实际资源保存于 `S/auxiliary_v1/resource_summary.json`。

18 个真实 Athena 进程 CPU 合计287.55 s，完整 worker CPU 合计354.95 s，最大RSS 1,038,264 KiB；逐作业分别读取 `tracking_resource.txt` 和 `resource.txt`，不将CPU合计当作并行wall time。

已成功的 clusters：构建 1582758；after/before 重建 1582759/1582763；默认 gate 1582767；完整分析 1582769；初版图 1582773；修正版图 1582781；v2图表数值/PDF验证1582786；q/p、资源和来源归约1582785。初版覆盖表文字重叠，v2 调整表格边界与说明位置后重新导出三页，不改变数据。正式使用 v2 图；三张PNG已按原始分辨率目视检查，所有绘图数值从paired states重算一致，逐图与合并PDF的归一化vector内容及page boxes一致。四项定向unit tests通过。

复现入口示例，已冻结任务名与输出目录不可覆盖，重跑使用新版本：

```bash
python3 scripts/prepare_charge_ckf_step_control.py \
  --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1 \
  --phase analysis-tasks

python3 scripts/submit_momentum.py \
  --tasks /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1/analyze_tasks_v1.json \
  --library /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/overlay_afterGeo/lib \
  --name <unique_task_name> --stage S06-CKF-step \
  --cpus 1 --memory-mb 6500 --disk-gb 4 --flavour microcentury \
  --exclude-schedd bigbird21.cern.ch
```

主产物为[三页英文图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1/figures_v2/results/ckf_step_report.pdf)、[完整覆盖表](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1/figures_v2/results/ckf_step_coverage.png)、[分析验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/ckf_step_control_v1/analysis_v1/validation.json)。逐图 PNG/PDF、绘图数值、冻结任务和软件来源均随结果保存。

[139页S06汇报PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v19/charge_stage06_report.pdf)由cluster1582788生成、exit0；其前136页的索引、decoded vector streams和page boxes精确继承report_v18，新增三页均为已检查的v2图，见[继承验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v19/lineage_validation.json)。新文件的`py_compile`、四项定向测试、本阶段更新文档的104个本地链接及`git diff --check`通过。

## 阶段性结论

事实：默认关闭新增属性精确重放原结果；明确限制 CKF 内部步长能够大幅压低 μ− 极端状态的版本间 mean/C 差异。候选、hits、sign 与 refit 成败在本批不变。较大的分叉来自对传播数值路径的敏感性；本次没有将步长效应解释为几何、场或材料 raw slab 本身的不同。

未测量：10 mm 控制在完整有限先验 likelihood posterior、额外宽度/分量数及生产大样本上的收敛与联合效果。约 10⁻⁸–10⁻⁷ 剩余 C 差异的原生材料方差记录正在独立控制中检验。本批结果不替代完整 posterior 或生产修复验收，S06 Q5 保持开放。

## 下一步计划

1. 记录真实 material interactor 的 `sigmaPhi2`、`sigmaTheta2`、`sigmaQoP2`；首先要求全部既有 propagation arrays 和旧 CSV fields 精确重放，再检查共同物理方向上的噪声 tensor。
2. 使用数值已受控的真实 Calypso 输出计算有限 Gaussian 的完整分布及 likelihood 加权结果，保留未选中、失败、分支和 fallback 的分母，测量 width/分量数收敛。
3. 对可用于正式配置的最小变化，按 Q5/S07 的独立文件、μ± 与 A′、方向/站数协议验收 joint momentum/angle/charge/efficiency/CPU；本批小分量云不替代该验收。
