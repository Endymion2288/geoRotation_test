# S06：固定完整 seed 和 cluster IDs 的测量坐标、过程噪声控制

**792 次实际 Calypso 重建已完成。六条原先纯错号轨迹在相同完整 seed、相同 admitted cluster IDs 下，SDO strip 中心控制恢复五条 selected 正确；仅关闭 CKF scattering 恢复一条。51839 的剩余错号不随 CKF scattering/energy-loss 开关消失，SiHit centroid 控制则在同一路径上恢复正确。**

这组结果将此前“坐标控制同时改变了候选路径”的解释进一步分开：六条轨迹的全部12种设置都保留原始 ordered measurement hits。它是按结果选取的成因实验，不用于估计生产样本的误判率改善。

## 阶段目标

继续 `problem.pdf` p4 / roadmap Q4–Q5，判断六条 primary-only 错号的坐标与材料过程项作用，并进一步核查上一轮 SDO 控制后仍错号的 100116 / source entry 51839。保持 isolated smoothing ON，所有 finding、update、smoothing、conversion/refit 均由 Calypso 实际组件执行。

验收要求为完整 seed mean、6×6 covariance、reference surface 和原始 cluster IDs 固定；逐步核对 measurement 值与 V；分别计数 sign、质量失败和 ordered path 变化。共同物理先验的完整 CKF 比较另行推进，未由本实验替代。

## 工作内容与输入

结果目录 `P=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/fixed_hit_process_v1`。

- 两类输入为 DSID100116/100117 的原始 `s0012` RDO，具体文件、source entry 和完整命令在 `P/tasks/`、`P/resolved_tasks_v1.json` 与逐作业 `source_task.json`。六条 pure wrong 为 μ+ entries 51839、224394，μ− entries 227965、390966、439421、458033；另有原先冻结的正常对照各30条，共66个源事件。
- 原始祖先 seed 与 hits 来自 `ckf_measurement_v1` 保存的实际 Calypso audit。每次只运行该 seed，并只接纳其原 selected 路径的 cluster IDs。finder 仍可拒绝 measurement，实际完整路径逐项比较，不能由 admitted 集合相同推断路径相同。
- 三类坐标为 native X、仅使用有 primary SDO 支持的原生 strip 中心均值、同 primary SiHit centroid；保留 localY、原 Lorentz shift、IDs、width/error。seed 在 cluster 坐标改变后显式恢复原始 mean/C/surface。
- 四类 CKF 过程设置为 material on、scattering off、energy loss off、两项 off。refit 和 observer 保留各自 native 设置。三坐标×四过程×66事件＝792次，不混入 refit material 的变化。
- Calypso worktree 为 `/eos/home-x/xcheng/FASER/geoRotation_calypso_charge_fixed_hits`，基底 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`；Athena24.0.41、FASERNU-04 / OFLCOND-FASER-05。53个 translation units 使用4CPU重编；tracking、observer、cluster 三组实际加载路径在每次 `run.json` 中核对。材料为本批已验证的 `material_encoding_control_v1/corrected_afs_map.json`。

增量 [0033](../patches/0033_fixed_hit_full_seed_input.patch) 提供默认关闭的完整 seed、运行指定 seed、cluster-ID 接纳与 CKF 材料开关。配置入口为 [slope_tracking.py](../configs/slope_tracking.py)，归约为 [charge_process_noise.py](../analysis/charge_process_noise.py)。未运行独立 tracking 或简化几何模型。

12个 identity/default 门槛作业各有58个保存字段与原基线精确一致，完整 mean/C/surface 逐值相同。全矩阵核对13,540次实际 measurement updates：使用的值与对应 cluster 审计相同，native CKF 的 V 始终为 `(0.0231 mm)^2`；flags、实际库和 admitted IDs 逐作业验证。

## 详细数据

下表为六条原 pure wrong 的 selected 正确数；括号内是质量失败数，符号正确但质量失败仍不计 selected 正确。

| 测量坐标 | material on | scattering off | energy loss off | 两项 off |
| --- | ---: | ---: | ---: | ---: |
| Native X | 0/6 | 1/6 | 0/6 | 1/6 |
| SDO strip X | 5/6 | 5/6 | 5/6 | 5/6 |
| SiHit centroid X | 5/6（1） | 5/6（1） | 5/6（1） | 5/6（1） |

全部六条、全部12设置的 ordered hits 与原始路径一致。centroid 的六条 final sign 均正确；439421 因半径质量条件失败，不能计作重建恢复。

60条 normal 在 native/SDO 两类坐标、全部过程设置均为60/60 selected 正确、60/60完整路径相同。centroid 则各为59/60 selected 正确、57/60路径相同，final sign 均正确；100116 / 109670 仅剩6 hits、3 layers、station bits=2、ndof=1，失败于三站/层数/measurement 数条件。该结果保留了真值坐标控制可能造成的效率代价。

![Fixed seed and hit-ID outcomes](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/fixed_hit_process_v1/figures_v3/results/fixed_hit_process_outcomes.png)

### 51839：过程项不消除剩余错号

μ+ 的正确 q/p 应为正。单位为 GeV⁻¹，全部使用同一完整 ordered path：

| 坐标 | material on 的 final q/p | scattering off 的 final q/p | selected 结局 |
| --- | ---: | ---: | --- |
| Native X | −0.00153412433 | −0.00153412712 | wrong / wrong |
| SDO strip X | −0.00002263651 | −0.00002263786 | wrong / wrong |
| SiHit centroid X | +0.00188565176 | +0.00188564662 | correct / correct |

energy loss off 和两项 off 也不改变上述结局。SDO 控制将 q/p 推近零，却仍在错误一侧；centroid 在原 hits 上恢复。其 reported covariance 已在图中给出，不能将一个近零、低显著性的错号当作已证明的系统偏置。

事实支持该轨迹对剩余 measured-coordinate/strip 分辨率效应敏感；本轮未用关闭物理过程的方法恢复它。原 SiHit centroid 与数字化 strip 中心之间还包含真实电荷输运、共享和量化，不能仅凭该对照将差值全部命名为量化误差，或声称 centroid 是可用于数据的生产修复。

![51839 predicted filtered and smoothed states](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/fixed_hit_process_v1/figures_v3/results/fixed_hit_51839_states.png)

### 439421：scattering 影响的是下游恢复能力

Native X 下 material on 的 final q/p 为 `+0.00259627373 GeV⁻¹`，scattering off 为 `−0.00199416795 GeV⁻¹`。两者都在 S2 L0、`z_F=1211.415 mm` 从正确负号穿到正号，关闭 scattering 没有避免首次错号。

S3 L0、`z_F=2401.415 mm` 的同一次 measurement update 则出现以下差别：

| 实际 update 量 | material on | scattering off |
| --- | ---: | ---: |
| predicted q/p / GeV⁻¹ | +0.00995489122 | +0.00985467543 |
| residual / mm | 2.25687804 | 2.33463863 |
| innovation variance / mm² | 0.50656200 | 0.34680479 |
| q/p gain / (GeV⁻¹/mm) | −0.00316380738 | −0.00528030703 |
| Δq/p / GeV⁻¹ | −0.00714032737 | −0.01232760876 |
| filtered q/p / GeV⁻¹ | +0.00281456385 | −0.00247293334 |

因此，该事件的材料过程项通过 covariance/gain 改变下游纠正能力。这是同 seed、同 hits 上的实际运算证据，不能据此推广为 scattering 在整个样本中产生了所有错号。SDO 坐标控制在保留材料时也恢复该条 selected 正确。

## 执行、验证与图表

全部规模运行与归约经 HTCondor；每次提交比较 schedd 广告及队列，选兼容的最低 running+idle 的 bigbird20，排除已有 EOS 兼容问题的 bigbird21。申请及实际资源、UTC、选择依据保留在 `P/*worker/submission.json` 与 Condor logs。完整矩阵计入的792次重建 CPU 合计9,196.23 s，最大 RSS 1,026,952 KiB；不将 CPU 之和当作并行实际历时。

初轮16个失败在独立 `recovery_v1` 原样重放；另两条先确认是 CKF 前 metadata 停顿，明确移除对应 live jobs 后在 `recovery_v2` 重放，保留原始观察与 removal closure。`resolved_tasks_v1.json` 为792个唯一 event/mode 的成功尝试，每个只有输出路径改变；未重复计入重试。归约断言覆盖全部66事件×12设置。

四页英文图的数值与冻结 summary 精确相同，PNG完整解码，单页与合并 PDF 的 vector graphics、page boxes 和 rotation 通过验证。四张最终 PNG 已逐张全尺寸目视检查；记录为 `P/figures_v3/results/visual_review.json`。早期 PDF 校验遇到 Matplotlib 局部 image/resource 别名差别，v3 解析实际 resource 后通过，不改变任何物理统计。

- [四页机制图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/fixed_hit_process_v1/figures_v3/results/fixed_hit_process_report.pdf)：结局矩阵、51839逐层、六轨迹 q/p及误差、路径与sign计数。
- [冻结统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/fixed_hit_process_v1/analysis_v2/full/summary.json)、[逐事件与逐update数值](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/fixed_hit_process_v1/analysis_v2/full/events.json)、[导图验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/fixed_hit_process_v1/figures_v3/results/validation.json)。

已有60项仓库单测通过，新增4项输入先验检查通过；包括同一矩阵的完整 covariance 保持、精确旋转的物理均值、有限宽度误差随收缩减小和非法分量数。三项 worker restart 检查包含在已有60项中，验证相同完成任务可复用、不同任务不能覆盖及原始失败返回值保留。

## 阶段性结论

事实：在完整 seed 与相同 ordered hits 上，测量坐标足以恢复大部分这组 pure wrong；scattering 对一条轨迹的下游 sign 恢复具有实际影响，不能解释全部六条。51839 的剩余错号不由本轮 CKF 材料过程开关移除；其测量坐标控制、近零曲率与误差保留为后续研究依据。60条正常对照和质量失败均已计数。

这些控制解释特定错号的敏感环节；truth 辅助坐标和关闭 scattering 均未作为生产配置推广。正负总体误判率、输入谱及历史 Alma 比较继续引用[全量基线](2026-11-13_stage06_charge_baseline.md)与[谱/模式控制](2026-11-13_stage06_charge_controls.md)。

## 下一步计划

1. 用真实 Calypso 的逐分量 CKF hook 完成共同物理先验检查，固定相同物理 seed 分布、正确材料查询及相同候选路径，分别测量 finite width 与分量数收敛。prototype 的“移动均值后仍保留完整 C”会展宽先验，已更正为分量间/分量内 covariance 分配；不能用此前 prototype 的分布宣称先验闭合。
2. 保留所有空 seed、finder failure、未通过选择及候选迁移，逐项给实际更新与后验差别；共同先验的完整 CKF 验收仍未完成。
3. 在 S06 总报告合并本次四页图与输入、cutflow、Alma参考图，不以这组 outcome-selected 样本给生产修复率。S05 的 isolated ON 全量角度结论继续由[独立全量复核](2026-11-13_stage05_isolated_validation.md)支持。
