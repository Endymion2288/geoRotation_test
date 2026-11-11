# geoRotation_test

研究旋转版 Calypso 中的 track 动量低估、`tan(theta_y)` 异常结构及 μ+/μ− 电荷误判差异。

- [完整研究路线图](roadmap/README.md)：执行顺序、阶段验收和三个研究分支。
- [材料与源码核查](roadmap/00_materials_and_source.md)：三份报告的证据、样本路径和当前实现中的优先检查点。
- [公共实验协议](roadmap/01_common_protocol.md)：坐标与 truth 定义、诊断字段、对照配置、统计和批量执行。
- [阶段 00 报告](workbook/2026-11-11_stage00_material_review.md)：本次实际核查与验证。
- [阶段 01 报告](workbook/2026-11-11_stage01_inputs_and_definitions.md)：输入 metadata、50 事件完整性、物理定义和局部 truth 记录。
- [阶段 02 报告](workbook/2026-11-11_stage02_rotation_and_audit_closure.md)：实际 sensor/field/material/接口检查、观察开关配对与来源审计。
- [Calypso patches 与构建](patches/README.md)；[冻结定义](configs/definitions.json)；[诊断字段](configs/fields.md)。
- [仓库协作规范](AGENTS.md)。

研究起点：`/eos/home-x/xcheng/FASER/geoRotation_calypso`，commit `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。现有构建为 Athena 24.0.41 / EL9；所有运行产物保存在 `/eos/home-x/xcheng/FASER/geoRotation_test_results`。

S00–S02 的材料、输入和诊断验证已记录；S03 的 50,000 事件 pilot 正在 HTCondor 运行。三个物理问题的总体根因按后续 S04–S06 实验验证。
