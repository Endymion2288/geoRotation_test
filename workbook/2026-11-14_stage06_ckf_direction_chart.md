# S06：实际 CKF update/smoothing 的方向参数化控制

**固定同一组有限先验输入，beforeGeo 的 regular update 使39个 CKF hit 路径差异全部消失；联合 regular smoothing 后，本批925个首选分量的最终电荷符号与 afterGeo 全部一致。同 hit 共同集合的 smoothed direction-C 差异中位数降至约5–7×10⁻⁹。两例原有 path/fallback 分歧变成两版均 fallback、均错号；该结果定位参数化效应，不代表生产误判率修复。**

## 阶段目标

继续 [Q5](../roadmap/04_charge_misidentification.md)，分解[实际 refit 输入报告](2026-11-14_stage06_refit_input_covariance.md)仍存的上游 CKF mean/C/path 差异。保持 isolated smoothing 开启，分别改变真实 CKF update 与 smoother 的方向 chart；检查首选存在性、hits、物理 mean/C、conversion、refit/fallback 与最终 sign。

结果根目录 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`；本次 `K=C/common_prior_ckf_chart_v1`，上游 `I=C/common_prior_refit_input_v1`。

## 工作内容与固定条件

输入为100044、100116、100117各首个 RDO 文件 entries `[0,20)`，共60个源事件，文件身份与上游报告相同。固定 `matched_n32_f256_v1`：每个原生 seed 32个中心，分量内 covariance 为原 C/256，保留相关项；分别640、544、576个输入分量。这些是相关的诊断分量，不是1760个独立物理事件。无原生 seed 的事件保留空输出。

beforeGeo 实际源码 `/eos/home-x/xcheng/FASER/geoRotation_calypso_ckf_chart_beforeGeo` 基于 master `1eb885ae22299fb39bb352bead6005d41b0f953e` 及上游0035–0037。新增[0038增量补丁](../patches/0038_native_ckf_direction_chart.patch)将已经验证的 refit chart adapters 原样抽到 `Tracking/Acts/FaserActsKalmanFilter/src/RegularDirectionKalman.h`，在真实 `TrackFinderFunction.cxx` 的 updater/smoother delegate 与 `CKF2.cxx` 的候选独立 smoother 中复用。`RegularCkfUpdateChart`、`RegularCkfSmoothingChart` 默认均关闭。方向代数使用远离 beam 极点的 chart，实际 geometry、field、material、local measurements 与传播保持本版原生实现；记录及后续传播前将 bound mean/C 还原为本版 chart。

构建为 `K/software_beforeGeo`，实际四个运行库来自 `K/overlay_beforeGeo/lib`。afterGeo 使用上游实际旋转源码 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` 及 `common_prior_component_refit_v1/overlay_afterGeo/lib`。完整源码、新增 headers、tracked diff、reverse-check 见 `K/source_validation.json`；0038在实际审计基底的 apply-check 见 `K/patch_validation.json`。

四种 before CKF 设置为 native、regular update、regular smoothing、both。**所有 before 控制均固定 refit regular update+smoothing 开启；两版 refit 均固定 `SeedCovarianceScale=1`，完整 C×10 保留。** after 参考固定为 `I/scale1_v1/.../afterGeo_native`。因此本表的 native 指 CKF chart；不要与上游 native-refit/default-scale10 的表混用。

## 详细数据：观察门槛与候选覆盖

先执行3个默认关闭重放：输入 seed、selected、全部 states、conversion、final、primary truth，以及真实 refit input/state 记录均与上游逐字段精确相同。合计1760个输入、925个 selected candidates、14,168条 selected state records；另核对真实 refit C 的复制、单位和 inflation 契约。证据为 `K/gate_v1/validation.json`。随后9个干预重建均 exit0，12组运行的实际库、开关及完整输入逐项验证通过。

下表每格依次为 native / update / smoothing / both；after 参考不变。首选依据 measurement 数降序、chi2升序、candidate ID升序，无 Good95 gate；全部 selected candidates 另行保存在原生审计中，本批计数与首选数相同。

| DSID | before/after selected | 相同有序 CKF hits | before 成功 refit | after 成功 refit | 最终版本间 sign 分歧 |
| --- | ---: | --- | --- | ---: | --- |
| 100044 μ+ | 311 / 311 | 306 / 311 / 306 / 311 | 38 / 37 / 38 / 37 | 37 | 1 / 1 / 1 / 0 |
| 100116 μ+ | 254 / 254 | 253 / 254 / 253 / 254 | 24 / 24 / 24 / 24 | 24 | 0 / 0 / 0 / 0 |
| 100117 μ− | 360 / 360 | 327 / 360 / 327 / 360 | 70 / 69 / 70 / 69 | 69 | 1 / 0 / 1 / 0 |

所有设置 CKF selected sign 与 after 均相同；更新 chart 改变的是上游方向/协方差及部分 hit 选择，未把这些 CKF 分量的错误电荷恢复为正确。CKF 错号数分别89、5、45，各开关不变。最终 before 错号 native→both 为88→89、5→5、44→45，after分别89、5、45。首选存在性没有迁移；39个有序路径差异为5+1+33，update 和 both 均全部消除。

![Native CKF chart outcomes](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/figures_v2/results/ckf_chart_outcomes.png)

## 详细数据：固定共同 states 的物理 mean/C

使用真实 audit surface，将 bound mean/C 解析表达为共同 FASER frame 的 position、unit direction、q/p、time；没有另行传播或拟合。direction-C 指标为 `max_ij |Cb_ij−Ca_ij|/sqrt(D_i D_j)`，`D_i=(|Cb_ii|+|Ca_ii|)/2`，仅取方向3×3子块，无量纲。

固定四种 CKF 控制中均与 after 使用相同有序 hits 的交集：100044为306个分量、1547个 measurement states；100116为253/1154；100117为327/1903。表中数值是该固定交集的中位数，不能把各设置变化后的样本混作同一分母。

| DSID | smoothed direction-C：native→update→smoothing→both | smoothed abs(Δq/p)：native→both / GeV⁻¹ |
| --- | --- | --- |
| 100044 | 0.135612 → 0.00164951 → 0.119760 → 5.14313×10⁻⁹ | 3.29315×10⁻⁷ → 6.59889×10⁻¹⁵ |
| 100116 | 0.0800195 → 0.00300257 → 0.0509335 → 5.64170×10⁻⁹ | 6.89175×10⁻⁸ → 8.38657×10⁻¹⁵ |
| 100117 | 0.188127 → 0.0203368 → 0.161356 → 7.33594×10⁻⁹ | 1.41466×10⁻⁶ → 8.25993×10⁻¹⁴ |

update 已将 filtered direction-C 中位数降到4.75×10⁻⁹、5.15×10⁻⁹、7.99×10⁻⁹；若 smoother 继续用原 near-pole chart，smoothed 中位数重新放大到上表的 update 列。只改 smoother 则保留 native 的 predicted/filtered 差异。这给出两个运算各自的直接因果证据。

联合控制没有达到每个浮点值完全相等。固定共同集合上，smoothed direction-C 的95th percentile分别5.08×10⁻⁸、3.08×10⁻⁸、4.06×10⁻⁶，最大值分别6.28×10⁻⁴、6.84×10⁻⁷、1.76×10⁻³；smoothed |Δq/p|最大为5.34×10⁻⁸、9.17×10⁻¹⁰、2.63×10⁻⁸ GeV⁻¹。全部 selected 的最终 |Δq/p|最大为1.32×10⁻⁶、2.35×10⁻⁸、9.12×10⁻⁶ GeV⁻¹。保留这些极端值，后续定位其首个增长操作，不以中位数代表全部状态闭合。

![Common physical states](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/figures_v2/results/ckf_chart_common_states.png)

## 详细数据：原四例与新增 update-only 分歧

原两条 same-hit 例的最终 q/p（GeV⁻¹）：100044 entry1/component23，before native CKF −0.002111883597，both −0.002043966797，after −0.002043966563，三者均错；entry3/component19分别 +0.001268384172、+0.001333094422、+0.001333094463，均正确。

原 path/fallback 两例在 update/both 后均与 after 使用相同 hits，并且两版均 fallback：

| DSID / entry / component | native CKF：before 最终 q/p | both：before / after 最终 q/p / GeV⁻¹ | 结果 |
| --- | ---: | --- | --- |
| 100044 / 12 / 2 | +0.000836452747 | −0.001581931348212 / −0.001581931348210 | 两版均错 |
| 100117 / 0 / 25 | −0.001072306801 | +0.000272719884131 / +0.000272719884131 | 两版均错 |

native CKF 的 before 曾有额外 hits 并成功 refit，之后恢复了正确 sign；update改变了实际路径，使该额外恢复不再发生。这里的版本间一致性不能当作质量改善。

只改 update 还出现100044 entry10/component19的新最终 sign 分歧：same hits、双方 refit 成功，before +0.007968821131、after −0.002644735480 GeV⁻¹。联合 smoother 后这例也不再分歧。该例纳入全分量计数；原四例图不代表全部 update-only 异常。

![Original four sign cases](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/figures_v2/results/ckf_chart_original_cases.png)

## 执行、验证与产物

准备入口 [prepare_charge_ckf_chart.py](../scripts/prepare_charge_ckf_chart.py)，统计 [charge_ckf_chart.py](../analysis/charge_ckf_chart.py)，绘图 [charge_ckf_chart_figures.py](../analysis/charge_ckf_chart_figures.py)。构建1582674、3个精确重放1582675、gate归约1582677、9个干预1582678、统计/首轮导图1582680均 exit0。构建4 CPU、28 translation units；重建每job申请1 CPU/3500MB/3GB scratch，12个实际 Athena 进程 user+system CPU合计190.48s、最大RSS1,030,548KiB，见 `K/reconstruction_resources.json`。提交前按实际普通兼容schedd的 running+idle 最低选择bigbird20，排除bigbird21；选择时间、负载及submit配置随各任务保存。

首轮图例遮挡两个计数，保留未验收记录；v2增加图例空间并区分median circles/p95 triangles。v2导图1582681 exit0，三页PNG已按原始分辨率逐张检查，所有计数、坐标、图例、fallback与hit标签可读。绘图独立从paired rows重算counts、共同集合及max/median/p95；单页与多页PDF的归一化vector content和page boxes一致，记录见下列验证文件。

汇编1582685 exit0，产出121页 S06 PDF；前118页的figure index、decoded vector content和page boxes与上一版精确相同，新增3页继承上述原始尺寸目视验收。新增Python入口编译及`git diff --check`通过。

- [权威统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analysis_v1/summary.json)、[7040条分量比较](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analysis_v1/paired_components.json)、[运行验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/analysis_v1/validation.json)。
- [3页英文图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/figures_v2/results/ckf_chart_report.pdf)、[绘图值及索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/figures_v2/results/figure_index.json)、[数值/PDF验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/figures_v2/results/validation.json)、[目视记录](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/figures_v2/results/visual_review.json)。
- [121页 S06 汇报图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v13/charge_stage06_report.pdf)、[前118页衔接验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v13/lineage_validation.json)。

实际提交使用冻结的 `K/gate_tasks.json`、`control_tasks.json`、`analysis_figure_sequence_tasks_v1.json`、`figure_sequence_tasks_v2.json`；例如：

```bash
python3 scripts/submit_momentum.py --tasks <K>/control_tasks.json \
  --library <K>/overlay_beforeGeo/lib --name <unique-name> \
  --stage S06 --cpus 1 --memory-mb 3500 --disk-gb 3 \
  --exclude-schedd bigbird21.cern.ch
```

## 阶段性结论

事实：在本批固定有限输入下，CKF update chart解释全部39条路径差异；update与smoothing联合控制大幅降低共同state的物理mean/C差异，并消除最终sign分歧。refit的额外对角inflation、CKF update chart、smoothing chart是已被实际单因素控制区分的三个来源。native版本间不一致不能全部归为物理旋转本身。

尚未验收：有限分量内Gaussian的非线性映射仍为一阶近似，未作likelihood加权posterior聚合或收敛证明；极端state数值差异的首次增长及完整宽度/分量数稳健性待查。生产μ±误判率、原生完整先验及效率不由这60事件的诊断分量替代。S05已有的角度成因交付保留，S06 Q5继续。

## 下一步计划

四组检查已经完成，结果及后续输入精度控制见[宽度扫描报告](2026-11-14_stage06_ckf_chart_width_scan.md)。

以相同真实Calypso联合chart及refit scale1，对既有四组宽度/分量数作完整before/after检查，复用本批n32 C/256，新增另外三组共18个重建。逐组报告选择、hits、sign、refit/fallback与物理state差异的max/p95/median，并对最差状态定位首个增长操作；不能仅用无sign分歧作为posterior等价证据。冻结入口 [prepare_charge_ckf_chart_scan.py](../scripts/prepare_charge_ckf_chart_scan.py)，任务保存在 `K/width_scan_v1`。
