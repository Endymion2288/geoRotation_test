# geoRotation_test

研究旋转版 Calypso 中的 track 动量低估、`tan(theta_y)` 异常结构及 μ+/μ− 电荷误判差异。

- [完整研究路线图](roadmap/README.md)：执行顺序、阶段验收和三个研究分支。
- [材料与源码核查](roadmap/00_materials_and_source.md)：三份报告的证据、样本路径和当前实现中的优先检查点。
- [公共实验协议](roadmap/01_common_protocol.md)：坐标与 truth 定义、诊断字段、对照配置、统计和批量执行。
- [阶段 00 报告](workbook/2026-11-11_stage00_material_review.md)：本次实际核查与验证。
- [阶段 01 报告](workbook/2026-11-11_stage01_inputs_and_definitions.md)：输入 metadata、50 事件完整性、物理定义和局部 truth 记录。
- [阶段 02 报告](workbook/2026-11-11_stage02_rotation_and_audit_closure.md)：实际 sensor/field/material/接口检查、观察开关配对与来源审计。
- [阶段 03 报告](workbook/2026-11-11_stage03_baseline.md)：50,000 事件基线、最终候选路径、动量尾、角度定义与 μ± 电荷误判统计。
- [阶段 04 报告](workbook/2026-11-12_stage04_momentum.md)：材料坐标根因、逐运算审计、Calypso内T0–T4、50k独立配对、seed与辐射后KF残余。
- [Calypso源码与结果来源](workbook/2026-11-12_stage04_calypso_provenance.md)：材料坐标问题的详细源码路径、实际任务/库/patch来源，以及独立toy退出验收的记录。
- [两份旋转材料图的比较](workbook/2026-11-12_stage04_material_maps.md)：AFS rotated12与strict current22的共同材料bin、空surface和法向2 mm平移差异。
- [S04 动量补图与原脚本核对](workbook/2026-11-12_stage04_reference_figures.md)：材料图旋转前后seed/九层CKF/final分布、good-seed条件图；[完整汇报PDF](../geoRotation_test_results/stage04_momentum/reference_figures/figures/momentum_material_comparison.pdf)。
- [各阶段参考图清单](roadmap/05_reference_figures.md)：S04–S07对应的参考页、图组和PNG/PDF交付要求。
- [Calypso patches 与构建](patches/README.md)；[冻结定义](configs/definitions.json)；[诊断字段](configs/fields.md)。
- [仓库协作规范](AGENTS.md)。

研究起点：`/eos/home-x/xcheng/FASER/geoRotation_calypso`，commit `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。现有构建为 Athena 24.0.41 / EL9；所有运行产物保存在 `/eos/home-x/xcheng/FASER/geoRotation_test_results`。

S00–S04已完成并分别提交阶段报告。S04证实本配置current22材料查表坐标与旋转几何不一致；严格旋转同一map，在五类独立样本上均降低动量负尾。110051的production p ratio<0.5从34.151%降至10.543%，Good95效率提高4.018个百分点。

完整truth seed和固定真hits仍可复现错误材料图的低估；正确材料下，剩余seed混合与电子辐射后KF偏差已量化并定位。S05角度结构、S06谱平衡后的μ±误判差异和S07完整修复验收继续按路线图执行。
