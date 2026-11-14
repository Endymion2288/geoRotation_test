# S06：原生 smoother 逐步审计与 PSD 修正检验

**实际 Calypso 的 13,232 对后向 smoothing 运算中，两版均没有触发 PSD 修正，raw C 与 validated C 逐值相同。因此本批剩余 covariance 差异在 PSD 检查前已经存在。原先 μ− 的 S1 极端例可沿后向执行顺序追到下游输入与递推；完整 posterior 闭合继续留在 Q5。**

## 阶段目标

延续[解析 seed covariance 单因素实验](2026-11-14_stage06_analytic_seed_covariance.md)，直接记录真实 ACTS smoother 的操作数、修正前后 covariance 和执行顺序，判断 PSD 修正是否造成剩余差异。验收包括默认关闭和开启记录的完整输出重放，以及共同 ordered-hit 路径的全部后向步骤覆盖。

结果目录 `O=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1`。权威统计为 `O/analysis_v2`，三页英文图为 `O/figures_v3/results`。

本次 isolated smoothing 始终开启。这里研究的 covariance 差异与 [S04 共享 smoothed slots 覆写](2026-11-12_stage04_smoothing_bump_cause.md)已分开控制。

## 工作内容与输入

100044、100116、100117 各首个 RDO 文件 entries `[0,20)`，共 60 个不同源事件。固定已有 n32、within-C/256 的有限分量，复用相同 before 输入与 after 解析 Jacobian 输入；分别有 640、544、576 个物理分量，共 1,760 个。输入文件、样本和任务原文在 `O/manifest.json`、`O/tasks/`、各运行的 `source_task.json`。分量具有相关性，不作为生产误判率的独立分母。

两版使用以下真实 Calypso worktree，均在 Athena 24.0.41 / ACTS 32.0.2 中构建：

| 配置 | 源码根目录 | branch 基准 commit | 材料 |
| --- | --- | --- | --- |
| beforeGeo | `/eos/home-x/xcheng/FASER/geoRotation_calypso_smoothing_audit_beforeGeo` | master `1eb885ae22299fb39bb352bead6005d41b0f953e` | `stage02_closure/material-current-source.json` |
| afterGeo | `/eos/home-x/xcheng/FASER/geoRotation_calypso_smoothing_audit_afterGeo` | geoRotation `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` | `stage06_charge/20261113_isolated_v1/physical_prior_material_bins_v4/material-rowmajor-rotated.json` |

before 沿用 CKF/refit regular update+smoothing；after 沿用原生旋转坐标。两版 refit 额外对角 scale=1、完整 C×10 均保持前次设置。几何库分别复用匹配的 `common_prior_ckf_chart_v1/software_beforeGeo` 和 `common_prior_component_refit_v1/software_afterGeo`，只重编 KalmanFilter package。完整 tracked diff、附加 headers、增量 patch 正向/反向检查和库信息在 `O/source_validation.json`、`O/*_complete_source.patch`、`O/software_*/software.json`。

新增只读接口位于两棵 worktree 的 `Tracking/Acts/FaserActsKalmanFilter/src/CKF2.h`、`CKF2.cxx`、`CandidateSmoothingAudit.h`。`SmoothingOperationFile` 与 `SmoothingOperationSelectionFile` 默认均为空。选择文件只控制记录对象，不限制重建输入。调用实际 `Acts::GainMatrixSmoother::calculate`；regular 分支保留既有方向变换与 periodic phi 处理。实际框架继续完成 finder、transport、update、conversion 和 refit。

记录当前 filtered mean/C、下一状态 predicted/smoothed mean/C、实际存储的 transport J、raw/validated C、native 输出、geometry/hit、state 索引及 operation 序号。gain 字段明确命名为 `reconstructed_gain`，由相同操作数及原公式重算；不声称直接读取了 ACTS 私有局部变量。raw C 在冻结 ACTS 32.0.2 的 destination accessor 第二次访问时捕获，包含零矩阵情形；全部记录的捕获标记均通过核对。

上游实际运算为 `G=CF·Jᵀ·CP_next⁻¹`、`xS=xF+G·(xS_next−xP_next)`、`CS_raw=CF+G·(CS_next−CP_next)·Gᵀ`，随后调用 native covariance validator。版本来源及原文见[前阶段源码检查](2026-11-14_stage06_analytic_seed_covariance.md)。原生计算结果供给保存的 track states。

## 详细数据：观察开关逐值重放

先运行默认 OFF 的六组，再比较 observer ON 的六组。两种模式均与前阶段冻结基线逐值一致：完整 audit stream，包括 found/seed/selected/all states、conversion/refit/final、truth 和错误/end 标记；每组另核对 58 个 observer 数组字段。每种模式共比较 58,160 条完整 audit 记录。实际加载的四个库、seed mean/C/surface、20-event 覆盖及 refit 控制也通过原有验证。

- [默认 OFF gate](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/gate_off_v1/validation.json)。
- [observer ON gate](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/gate_on_v1/validation.json)。

| DSID | 每版选中分量 | 每版 selected states | 只有终端 state 的分量 | 有后向运算的路径 | 配对后向运算 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 100044 | 311 | 4,157 | 7 | 304 | 3,846 |
| 100116 | 254 | 3,580 | 0 | 254 | 3,326 |
| 100117 | 360 | 6,420 | 2 | 358 | 6,060 |
| 合计 | 925 | 14,157 | 9 | 916 | 13,232 |

每条路径的终端 state 取 `smoothed=filtered`，后向运算数等于 states 数减一；因此 14,157−925=13,232。九个单 state 分量保留在 selected 分母，无需产生后向步骤。所有执行序号连续、前后 state 索引衔接、共同 geometry/hit/measurement flags 及 inherited mean 均核对通过。

![Native audit coverage](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/figures_v3/results/smoothing_operations_counts.png)

## 详细数据：PSD 检查前后的 covariance

before 操作数已经位于 regular `(phi_A, latitude_A)` chart；after 原生 polar chart 以 `latitude_A=pi/2−theta_A` 线性重表达，`D=diag(1,1,1,−1,1,1)`，C/J/gain 均在对应两端应用 D。局部 surface 坐标、q/p 和 time 保留原定义。

差异指标为 `max |Cb_ij−Ca_ij|/sqrt(d_i d_j)`，其中 `d_i=(|Cb_ii|+|Ca_ii|)/2`，取 common chart 的方向 2×2 子块，无量纲。它与前阶段 Cartesian unit-direction 3×3 C 指标定义不同；数值大小不直接作为前后改善量。

| DSID | raw direction-C 差异中位数 | raw 最大值 | validated 最大值 | before / after 实际 PSD 修正次数 |
| --- | ---: | ---: | ---: | ---: |
| 100044 | 1.1609×10⁻¹⁰ | 7.6147×10⁻⁶ | 7.6147×10⁻⁶ | 0 / 0 |
| 100116 | 6.4287×10⁻¹⁰ | 3.5521×10⁻⁸ | 3.5521×10⁻⁸ | 0 / 0 |
| 100117 | 2.7367×10⁻⁹ | 1.94035×10⁻³ | 1.94035×10⁻³ | 0 / 0 |

两版各 13,232 个 raw C 均通过 native PSD 判定，validated C 与 raw C 逐值相同。由此排除**本批记录的 private smoothing 运算**中 PSD 修正造成版本差异的解释。没有将此结果外推为所有 Calypso 样本都不触发该修正。

![Raw and validated covariance differences](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/figures_v3/results/smoothing_operations_covariance.png)

## 详细数据：μ− 极端例的真实后向顺序

最大 common-chart direction-C 差异仍在 100117、entry 1、component 24 的路径上。以下 operation 从下游向上游实际执行，含非 measurement states；位置为原 audit 的 path length，并非 station z。

| operation | path / mm | measurement | 当前 filtered-C 差异 | inherited next-S-C 差异 | raw output-S-C 差异 |
| --- | ---: | --- | ---: | ---: | ---: |
| 0 | 2461.039 | 是 | 8.5046×10⁻⁶ | 8.4443×10⁻⁶ | 8.4447×10⁻⁶ |
| 11 | 1256.921 | 否 | 1.6194×10⁻⁴ | 1.1562×10⁻⁵ | 1.6019×10⁻⁴ |
| 16 | 1193.919 | 否 | 5.4082×10⁻⁶ | 1.8239×10⁻⁴ | 1.94035×10⁻³ |
| 17 | 133.934 | 否 | 7.1930×10⁻⁹ | 1.94035×10⁻³ | 1.57975×10⁻³ |
| 19 | 73.002 | 是 | 6.8454×10⁻⁹ | 1.57428×10⁻³ | 1.57426×10⁻³ |

operation 19 就是前报告 hit `10183852618321756160` 的位置。其 next-predicted C 的对角标准化 condition number 约 3.195，局部 filtered C 差异约 7×10⁻⁹；传入的 next-S 已有约 1.574×10⁻³ 差异。该步骤的 smoothed q/p 差异为 7.25521×10⁻⁹ GeV⁻¹。因此不能仅由 S1 当层 P/F 很接近而 S 不同，断言问题始于该层 update 或局部求逆。

阈值 10⁻⁴ 下，该路径的首个 backward output-C 超阈值位置为 operation 11，但其当前 filtered C 已超阈值。operation 0 也已有非零前向差异。记录支持继续追踪前向 C 和后向输入传播；尚未证明显式求逆是原因。上述指标使用各自矩阵的方差归一化，指标升高本身不能替代绝对误差或因果分解。

![Actual backward execution order](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/figures_v3/results/smoothing_operations_extreme.png)

完整矩阵与路径在 `O/analysis_v2/extreme_operation_paths.json`；所有配对步骤在 `paired_operations.jsonl.gz`，按阈值的首个执行位置在 `first_backward_growth.jsonl.gz`。阈值只作定位，全部后向运算进入统计。

## 执行、验证与产物

构建 clusters 1582695/1582696 均 exit 0，各编译 28 个真实 KalmanFilter translation units。before/after 的编译 CPU user+system 分别 1151.97/1138.88 s，wall 分别 6:10.11/6:07.71，`time -v` 最大 RSS 分别 2,285,328/2,608,852 KiB；每 job 申请 4 CPU、10 GB memory、3 GB scratch，编译并行度为 4。

默认 OFF 重建 clusters 1582698/1582699、ON 重建 1582701/1582702 共 12 个 jobs 全部 exit 0；实际 Calypso 重建 CPU 合计 193.45 s，最大 RSS 1,040,628 KiB。每 job 申请 1 CPU、3500 MB、3 GB scratch，Athena 单 event slot，BLAS/OpenMP 线程为 1。

默认 gate 1582700 成功；ON gate 与首次统计在 1582704 成功。绘图首次因单 state 分母检查停止，补齐统计后 1582705 的数值分析成功；随后绘图读取 Counter 省略的零计数项时停止。修正为显式按零读取，1582707 的导图及数值/PDF 验证均成功。两次停止均在绘图核验，所有重建、mean/C、hits、cuts 和 gate 保持不变。权威版本为 `analysis_v2` / `figures_v3`。

三张 PNG 均按原始分辨率目视检查通过；图中统计由逐 operation 数据重算，single/multipage PDF 的规范化 vector content、page boxes 和 rotation 一致。六项 focused finite-prior 坐标/输入 tests 通过，[日志](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/unit_tests.log)。所有提交前重新查询兼容 schedd 负载与队列，选择 bigbird20，排除 bigbird21；UTC、选择依据、实际 cluster 和请求资源在各 `*_worker/submission.json`。

实际准备入口：

```bash
python3 scripts/prepare_charge_smoothing_audit.py --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1 --phase source
python3 scripts/prepare_charge_smoothing_audit.py --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1 --phase tasks
```

构建后 `--phase overlay --variant beforeGeo` / `afterGeo` 组装匹配四库。执行任务依次保存在 `O/*_build_tasks_v1.json`、`O/*_off_tasks.json`、`O/gate_off_tasks_v1.json`、`O/*_on_tasks.json`、`O/observer_analysis_sequence_tasks_v1.json`、`O/analysis_figure_sequence_tasks_v2.json`、`O/figure_sequence_tasks_v3.json`；均由 [submit_momentum.py](../scripts/submit_momentum.py) 提交。新增统计与绘图代码为 [charge_smoothing_operations.py](../analysis/charge_smoothing_operations.py)、[charge_smoothing_operations_figures.py](../analysis/charge_smoothing_operations_figures.py)。

- [权威统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/analysis_v2/summary.json)、[执行配对验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/analysis_v2/validation.json)。
- [三页英文图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/figures_v3/results/smoothing_operations_report.pdf)、[导图验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/figures_v3/results/validation.json)、[目视检查](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/smoothing_operations_v1/figures_v3/results/visual_review.json)。
- [130页 S06 汇报图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v16/charge_stage06_report.pdf)、[前127页衔接验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v16/lineage_validation.json)。汇编 cluster 1582708 exit 0，前127页的 index、decoded vector content 和 page boxes 与上一版精确相同，新增三页继承已通过的目视记录。

## 阶段性结论

事实：新增记录的默认与开启两种模式均完整逐值重放；本批全部 PSD 修正计数为零，差异已经存在于 raw covariance 计算及其输入。S1 极端例包含从下游 next-S 继承的差异，部分更早的 filtered C 已经不同。当前观测排除了本批 PSD 修正解释，尚未区分前向 chart/transport 数值路径、非线性 Gaussian 表示及 gain 运算各自的贡献。

有限分量 sign 一致和观察开关重放，均不替代完整物理 prior/posterior 的等价检验。既有生产主率与参考分布继续使用前阶段结果；本次没有新增全样本生产率测量，Q5 保持开放。

## 下一步计划

1. 在实际 Calypso 前向路径中定位第一处明显的物理 C 分岔，记录 surface 前后 mean/C、transport J 与过程噪声。固定同一完整物理输入到同一目标面，检查旋转后的 native covariance transport，分离继承输入、chart Jacobian 和 material 项。验收需要可逐 surface 复查的操作数及共同物理 C。
2. 对 operation 11、16、17、19 进行 mean/C 差异的操作数分解；求逆/solve 或 chart 的单因素干预在真实 Calypso 中执行，以相同输入重放和输出迁移检验解释。condition number 单项不作为原因结论。
3. 继续有限 Gaussian 的分布映射、likelihood 加权 posterior 及 width/分量数收敛；保留所有无选中、terminal-only、refit/fallback 类别，再按 Q5/S07 的留出样本和方向/站数协议验收。
