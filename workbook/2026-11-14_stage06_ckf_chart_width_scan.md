# S06：联合 CKF chart 控制的宽度与分量数检查

**四组有限先验设置的4,056个共同首选分量全部保持相同有序 hits，CKF、conversion、final 的版本间 sign 分歧均为0。smoothed direction-C 差异中位数约为5×10⁻⁹至10⁻⁸，但最大值仍达1.757×10⁻³，最终 |Δq/p|最大为9.115×10⁻⁶ GeV⁻¹。368个双方成功 refit 与3,688个双方 fallback 分开计数；该批结果不证明完整 posterior 收敛。**

## 阶段目标

延续[实际 CKF chart 四模式控制](2026-11-14_stage06_ckf_direction_chart.md)，检查联合 regular update/smoothing 在既有四组分量数和内宽度下的稳健性。保留全部 states、极端值及真实 native mean/C/measurement/Jacobian，找出沿路径首先出现较大差异的测量位置，为下一项单因素控制提供输入。

结果根目录 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`，本次 `W=C/common_prior_ckf_chart_v1/width_scan_v1`。

## 工作内容

使用100044、100116、100117各首文件 entries `[0,20)`，共60个不同源事件。四组为 n32 C/256、n32 C/1024、n64 C/1024、n64 C/4096；C为原始 beforeGeo seed covariance。组间分量相关，不能把10,560个输入分量当作独立事件或生产率分母。

复用已完成的n32 C/256两版六次运行，新增其他三组的18次真实 Calypso 重建。beforeGeo 使用 `common_prior_ckf_chart_v1/overlay_beforeGeo/lib`，CKF与refit均开启 regular update/smoothing；afterGeo 使用 `common_prior_component_refit_v1/overlay_afterGeo/lib` 的原生旋转实现。两版都开启 isolated smoothing，并固定 refit 额外对角 scale=1，保留完整 C×10。材料使用已经由实际 ACTS decoder 验证的对应图。全部中心、covariance、surface、loaded libraries及20-event覆盖由逐运行 audit 核对。

准备入口为 [prepare_charge_ckf_chart_scan.py](../scripts/prepare_charge_ckf_chart_scan.py)，统计为 [charge_ckf_chart_scan.py](../analysis/charge_ckf_chart_scan.py)，绘图为 [charge_ckf_chart_scan_figures.py](../analysis/charge_ckf_chart_scan_figures.py)。本次仍使用历史 central-difference covariance Jacobian；没有重写旧输入。

## 详细数据：全部输入与选中覆盖

| 分量数 / C除数 | DSID | 输入 | before / after selected | same ordered hits | 两版成功 refit |
| --- | --- | ---: | ---: | ---: | ---: |
| 32 / 256 | 100044 | 640 | 311 / 311 | 311 | 37 / 37 |
| 32 / 256 | 100116 | 544 | 254 / 254 | 254 | 24 / 24 |
| 32 / 256 | 100117 | 576 | 360 / 360 | 360 | 69 / 69 |
| 32 / 1024 | 100044 | 640 | 209 / 209 | 209 | 10 / 10 |
| 32 / 1024 | 100116 | 544 | 186 / 186 | 186 | 8 / 8 |
| 32 / 1024 | 100117 | 576 | 292 / 292 | 292 | 38 / 38 |
| 64 / 1024 | 100044 | 1280 | 425 / 425 | 425 | 29 / 29 |
| 64 / 1024 | 100116 | 1088 | 369 / 369 | 369 | 14 / 14 |
| 64 / 1024 | 100117 | 1152 | 584 / 584 | 584 | 74 / 74 |
| 64 / 4096 | 100044 | 1280 | 303 / 303 | 303 | 13 / 13 |
| 64 / 4096 | 100116 | 1088 | 262 / 262 | 262 | 5 / 5 |
| 64 / 4096 | 100117 | 1152 | 501 / 501 | 501 | 47 / 47 |

总计10,560个输入、4,056个共同首选、19,440个配对 measurement states。选择存在性、ordered CKF/final hits、refit成功状态逐分量一致。两版各640个最终错号分量，CKF→final sign无迁移；这些只是有限先验诊断计数。成功refit合计368，fallback合计3,688；不能将其合并称为4,056次成功refit。

![Width scan coverage](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/width_scan_v1/figures_v1/results/ckf_chart_scan_coverage.png)

## 详细数据：物理 state 差异

每个设置内使用全部 same-ordered-hit 集合。共同 FASER frame 的 position、unit direction、q/p、time 与 C 只作坐标表达，没有额外传播。direction-C 指标沿用上游定义，`max |Cb_ij−Ca_ij|/sqrt(D_i D_j)`，`D_i=(|Cb_ii|+|Ca_ii|)/2`，仅取方向3×3子块。下表各项为无量纲 smoothed direction-C 的中位数 / 最大值。

| 分量数 / C除数 | 100044 | 100116 | 100117 |
| --- | --- | --- | --- |
| 32 / 256 | 5.185×10⁻⁹ / 6.282×10⁻⁴ | 5.646×10⁻⁹ / 6.845×10⁻⁷ | 8.015×10⁻⁹ / 1.757×10⁻³ |
| 32 / 1024 | 5.333×10⁻⁹ / 1.133×10⁻⁶ | 5.816×10⁻⁹ / 4.991×10⁻⁸ | 9.868×10⁻⁹ / 2.522×10⁻⁴ |
| 64 / 1024 | 4.889×10⁻⁹ / 7.593×10⁻⁵ | 4.764×10⁻⁹ / 8.153×10⁻⁸ | 9.025×10⁻⁹ / 6.318×10⁻⁴ |
| 64 / 4096 | 5.409×10⁻⁹ / 2.112×10⁻⁵ | 5.878×10⁻⁹ / 6.712×10⁻⁸ | 1.006×10⁻⁸ / 5.442×10⁻⁴ |

该分母不同于上游四个chart控制的固定交集，因此n32 C/256的中位数不要求与上游交集表相同。各组全部 P/F/S 的max、p95、median保存在 `W/analysis_v2/summary.json`。最终 |Δq/p|全局最大值为9.115269×10⁻⁶ GeV⁻¹，并未随分量内宽度单调消失。宽度改变也明显改变selected及refit成功数，当前尚未计算likelihood加权的mixture posterior。

最大 smoothed direction-C 对应100117、n32 C/256、entry1/component24、hit10183852618321756160，路径73.00156 mm。before seed角为 `[1.014481689765, 0.006538064019]` rad，该hit的P/F/S theta分别0.006540536315、0.006012723285、0.005577846488 rad；after seed角为 `[0.005552191574, 1.574248800157]` rad。完整native参数、C、surface、measurement和transport Jacobian随 `extreme_native_states.json` 保存，避免只保留汇总指标。

![All-state residuals](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/width_scan_v1/figures_v1/results/ckf_chart_scan_residuals.png)

## 详细数据：沿路径的首次超阈位置

描述阈值为方向差norm>10⁻⁸、|Δq/p|>10⁻⁹ GeV⁻¹、direction-C>10⁻⁵；不是接受或删选条件。按measurement path length搜索首个位置，再在该位置依次查看P/F/S。**这不是算法执行顺序中的首个因果运算**：该处smoothed值可能包含下游来源。

100117各设置的首位置 P/F/S 计数如下；未超过阈值的分量仍进入全状态统计。

| 分量数 / C除数 | 选中分量 | direction | q/p | direction-C |
| --- | ---: | --- | --- | --- |
| 32 / 256 | 360 | 11 / 0 / 9 | 1 / 1 / 10 | 4 / 2 / 7 |
| 32 / 1024 | 292 | 5 / 0 / 8 | 0 / 3 / 4 | 3 / 2 / 4 |
| 64 / 1024 | 584 | 21 / 0 / 13 | 0 / 0 / 14 | 11 / 5 / 7 |
| 64 / 4096 | 501 | 19 / 0 / 8 | 0 / 0 / 4 | 6 / 2 / 5 |

![First measurement locations](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/width_scan_v1/figures_v1/results/ckf_chart_scan_first_growth.png)

## 执行、验证与产物

before新增9次重建cluster1582682，after新增9次1582684，全部exit0。新增18个Athena进程CPU user+system合计279.10 s，最大RSS1,042,040 KiB，详见 `W/reconstruction_resources.json`。每job申请1 CPU/3500MB/3GB scratch，线程数=1。

首个分析1582686在读取重建前因冻结目录漏带 `charge_ckf_chart.py`、`charge_refit_input.py`而导入失败。补齐完整依赖，冻结目录 `--help` 通过后重新提交1582687，exit0；没有重复重建。失败与修复见 `W/analysis_dependency_fix.json`。权威分析目录为 `analysis_v2`。

绘图1582688 exit0，独立从paired rows重算所有分组计数及P/F/S物理指标，plot values完全一致；三个PNG已按原始尺寸逐张检查，标题、legend、表格与footer清楚，无裁切或重叠。单页/多页PDF的归一化vector content和page boxes一致。汇编任务1582689将三页追加到前121页，产出124页；前121页的index、decoded vector content和page boxes逐页精确保留，新增三页继承原始尺寸目视记录。

全部提交都由 `submit_momentum.py` 重新查询兼容schedd的running+idle，选择bigbird20，排除bigbird21；查询、队列、UTC时间、cluster与资源写入各worker的 `submission.json`。复现任务为 `W/beforeGeo_tasks.json`、`afterGeo_tasks.json`、`analysis_tasks_v2.json`、`figure_sequence_tasks_v1.json`、`report_tasks_v1.json`；其中analysis实际冻结命令以对应task JSON为准。

- [权威统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/width_scan_v1/analysis_v2/summary.json)、[全部分量](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/width_scan_v1/analysis_v2/paired_components.json)、[输入与运行验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/width_scan_v1/analysis_v2/validation.json)。
- [三页英文图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/width_scan_v1/figures_v1/results/ckf_chart_scan_report.pdf)、[绘图验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/width_scan_v1/figures_v1/results/validation.json)、[目视记录](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/width_scan_v1/figures_v1/results/visual_review.json)。
- [124页S06图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v14/charge_stage06_report.pdf)、[衔接核对](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v14/lineage_validation.json)。

## 阶段性结论

事实：联合chart控制对四组设置均消除本批版本间选择、路径和sign差异，保留非零的极端state差异。分量内宽度会改变绝对接受和fallback，不能把sign一致或缩小内宽度当作原始完整Gaussian posterior等价。

待检验：输入准备中的F→A covariance Jacobian仍是固定步长central difference，原输入的物理C差异可达约10⁻⁷。它是否贡献当前差异，需保持全部中心、before C与surface不变，只以解析Jacobian替换after C进行实际Calypso单因素重建；目前未将其作为已证明原因。

## 下一步计划

解析Jacobian单因素已经完成，见[实际重建与结果](2026-11-14_stage06_analytic_seed_covariance.md)。输入C精度贡献得到分离，剩余极端差异继续按真实后向执行顺序诊断。

新增可选解析Jacobian，默认保留历史实现。检查完整position/direction/qop/time covariance及近极点与角度等价表示，冻结只改after C的新输入。复用现有before运行，比较after两种输入下的ordered hits、P/F/S、refit/fallback和最终sign；保留全部max/p95/median与极端状态。Q5保持未完成，posterior聚合、纯旋转证据与生产修复的验收分别记录。
