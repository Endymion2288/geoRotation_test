# geoRotation_test

研究旋转版 Calypso 中的 track 动量低估、`tan(theta_y)` 异常结构及 μ+/μ− 电荷误判差异。

- [完整研究路线图](roadmap/README.md)：执行顺序、阶段验收和三个研究分支。
- [材料与源码核查](roadmap/00_materials_and_source.md)：三份报告的证据、样本路径和当前实现中的优先检查点。
- [公共实验协议](roadmap/01_common_protocol.md)：坐标与 truth 定义、诊断字段、对照配置、统计和批量执行。
- [阶段 00 报告](workbook/2026-11-11_stage00_material_review.md)：本次实际核查与验证。
- [仓库协作规范](AGENTS.md)。

研究起点：`/eos/home-x/xcheng/FASER/geoRotation_calypso`，commit `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。现有构建为 Athena 24.0.41 / EL9；所有运行产物保存在 `/eos/home-x/xcheng/FASER/geoRotation_test_results`。

本次完成材料研究与实验设计；三个物理问题的当前版本发生率和根因留待路线图规定的实验测量。
