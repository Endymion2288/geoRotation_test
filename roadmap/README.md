# 几何旋转后的 tracking 研究路线图

本研究要给出三个可复核的答案：动量低估首次出现在哪个运算或候选选择环节、角度分布的剩余结构来自物理混合还是实现问题、μ− 比 μ+ 更易误判电荷的差异在旋转后是否存在及其来源。完成标准包括事件级因果证据、总体发生率及其统计区间、最小修复或经验证的物理解释，以及独立样本回归。

研究起点为本地旋转版Calypso `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。S00–S04已交付：S03建立五类各10k事件的原版基线，S04完成动量逐运算研究、独立50k母事件的材料配对验证及Calypso内T0–T4控制。此前独立ACTS toy退出验收，T5在Calypso实际几何/组件内的完整统计closure未测量，列入S07；[源码与软件来源](../workbook/2026-11-12_stage04_calypso_provenance.md)列明实际路径与patch。

实测入口：[S01](../workbook/2026-11-11_stage01_inputs_and_definitions.md)、[S02](../workbook/2026-11-11_stage02_rotation_and_audit_closure.md)、[S03](../workbook/2026-11-11_stage03_baseline.md)、[S04](../workbook/2026-11-12_stage04_momentum.md)。S04证实current22材料与旋转几何的坐标不一致，主要通过scattering covariance/gain制造大负尾；严格旋转同一map后，独立110051尾率从34.151%降至10.543%，Good95效率提高4.018个百分点。接口chart/refit起点问题、seed混合和辐射后KF残余分别保留受控实验入口。补充G4退出134已定位修复，同sensor truth关联通过2400事件、81153 steps逐值验证。

110051完整输入的最新动量图见[全20文件专项报告](../workbook/2026-11-12_stage04_full110051_pfs.md)：20,000源事件，材料坐标修正前后final尾率33.858%→10.577%。该报告补齐原10文件逐层图，给出predicted/filtered/smoothed的信息范围、同state尾部迁移、曲率不确定度及7页图集；描述性全样本与上述独立验收样本分别报告。

下表保留S03的current22历史基线；strict与独立样本的完整配对表见S04。

| DSID | Good95 匹配粒子 | production p ratio <0.5 / % | charge misID / % |
| --- | --- | --- | --- |
| 100044 | 6965 | 11.357 | 5.298 |
| 100116 μ+ | 4868 | 1.541 | 0.103 |
| 100117 μ− | 5047 | 2.180 | 0.674 |
| 110033 | 12325 | 12.252 | 2.150 |
| 110051 | 12460 | 33.563 | 13.604 |

以上使用 fiducial+H18、SDO purity≥0.7、Good95；完整 cutflow、区间与匹配歧义检查见 S03。μ−−μ+ 总误判率差为 0.571 个百分点，bootstrap95 为 [0.333,0.829] 个百分点；输入谱不同，尚未隔离算法的电荷不对称或旋转净效应。角度的物理定义已恢复，剩余峰的显著性未测量。

## 阅读与执行入口

| 文档 | 内容 |
| --- | --- |
| [00 材料与源码](00_materials_and_source.md) | PDF 页码、历史数值、输入样本、源码检查点及尚缺信息 |
| [01 公共协议](01_common_protocol.md) | frame/units/truth、运行矩阵、诊断数据、统计与 HTCondor |
| [02 动量研究](02_momentum_bias.md) | segment→seed→CKF→conversion→refit→output 的逐环节审计与消融 |
| [03 角度研究](03_slope_structure.md) | 物理方向恢复、峰的量化、条件分布与数值稳定性检查 |
| [04 电荷研究](04_charge_misidentification.md) | Alma 9 复现、配对比较、谱匹配、符号首次翻转与对称性实验 |
| [05 参考图与汇报产物](05_reference_figures.md) | 参考PDF页码到各阶段图清单、共同绘图规范、PNG/PDF与阶段完成标准 |

## 目前证据决定的优先级

1. **S05在current/strict配对上统一物理面。** S04已证明材料能显著改变选择与curvature；先将tanθ_y和truth传播到同一st1/VetoNu面，再分解reco-p窗口、hit pattern和candidate迁移。量化峰的位置、宽度、占比及统计显著性，不能用production方向residual直接充当角分辨率。
2. **S06先平衡谱与接受度。** 独立strict样本μ+、μ−各剩1条wrong；事件数尚不足以判定残余正负差异。先比较共同truth p/位置/方向支持域与有效样本量，再按目标精度追加。材料配对不代替完整G0/G1的旋转净效应。
3. **S07逐项处理seed和辐射模型。** 110051 strict有3020条最终祖先seed tail，其中2346条含other-majority segment。固定真hits的电子仍有23条guaranteed local tail，11条在update、12条在smoothing出现同sensor可证实的穿越，并有显著沿程能损。使用纯segment组合/同点方向与field integral、G4 bremsstrahlung控制及piecewise q/p/GSF分别检验。
4. **保留接口、pruning与IFT实验。** S04全量100k观察配对的final p/charge/Good95均一致，少量rmax改变；转换mean闭合不代表covariance已修复。4-station finding会使用IFT，但既有refit将IFT判outlier；完整四站refit须显式控制该策略。候选pruning损失也有独立计数。
5. **验证使用新的保留输入。** S04的独立50k已用于材料验收，后续修复使用新的文件/ranges或随机种子；分别报告效率、fakes、failure、angle和charge。历史truth界待定的类别继续保留发生率和输入定位。

## 阶段、依赖与验收

各阶段规模为预算起点，按有效 truth 粒子数、误判数和测量精度调整。凡达到批量规模的 reconstruction/simulation/扫描均走 HTCondor；交互节点只做少量事件与配置调试。

| 阶段 | 依赖与问题 | 实施与起始规模 | 必交付结果与验收 |
| --- | --- | --- | --- |
| S00 材料研究 | 输入报告与本地软件 | 阅读 3 份 PDF，检查源码、样本路径与现成构建 | 本路线图、补全 AGENTS、中文 workbook；已完成，验证见阶段报告 |
| S01 输入与定义冻结 | S00；历史图究竟比较了什么 | 已完成 metadata/schema 与每类10事件检查，冻结 Good Track、fiducial、frame、units、truth surface | manifest、诊断入口和字段字典已交付，50事件计数闭合，见S01报告；当时的补充G4退出134已在S04修复并验证 |
| S02 旋转与诊断闭合 | S01；观察量与 frame 是否可靠 | 已完成 sensor/field/material/transport/interface 与观察开关配对；保留原版并提供默认关闭的实验修正 | 来源唯一、观察开关前后输出一致；chart/refit问题与material差异已量化，见S02报告；总体物理影响留待S04–S06 |
| S03 原版复现与基线 | S01；与 S02 可交错做原始量检查 | 已完成五类各10k事件，50/50 shards exit0；紧凑审计归档后完整复算 | particle表、图、cutflow、区间、32条代表track与资源已交付，见S03报告；历史selection未恢复、同面truth与完整G0未测量均已记录 |
| S04 动量因果定位 | S02–S03 | current/strict全量逐层、Calypso内T0–T4、独立50k配对、backward/4-station/真实merge已完成；T5完整Calypso统计closure列入S07 | 材料主因与独立改善、seed类别和23条电子local tail运算证据、完整剩余类/资源/测试；[S04报告](../workbook/2026-11-12_stage04_momentum.md) |
| S05 角度结构定位 | S02–S03；共享 S04 的误差与材料观察 | [角度实验 A1–A5](03_slope_structure.md)，truth 谱与 detector/hit pattern 条件分解 | 峰位置、宽度、占比、稳定性和对应物理/程序来源；剩余结构在 forward model 中可复现或被最小修复消除 |
| S06 电荷与正负差异 | S02–S03；共享 S04 的 curvature trace | [电荷实验 Q1–Q5](04_charge_misidentification.md)；先 pilot，再扩至每电荷至少约 100k 有效选中粒子或所需精度 | before/after 配对表、误判率/效率/差值区间、谱匹配结果、首次符号翻转与解释；无显著差异时给可排除范围 |
| S07 修复与独立验收 | S04–S06 各分支形成可检验解释 | 每次一个修改；保留未调参文件作验证；覆盖μ±、A′两样本、forward/backward、3/4 station；补齐实际Calypso内T5统计closure | 物理指标、效率、失败率与资源回归；中文最终报告、配置、最小patch、测试和各阶段commits |

S04、S05、S06 在公共定义与诊断闭合后可以独立推进，结果互相反馈。某一分支暂未找到原因时，记录具体未解类别并继续其他分支；少量离群事件不足以否定整条研究路线。

```mermaid
flowchart TD
    S0["S00: Materials and source review"] --> S1["S01: Inputs and definitions"]
    S1 --> S2["S02: Frame and diagnostics closure"]
    S1 --> S3["S03: Baseline reproduction"]
    S2 --> M["S04: Momentum mechanism"]
    S3 --> M
    S2 --> A["S05: Slope structure"]
    S3 --> A
    S2 --> Q["S06: Charge asymmetry"]
    S3 --> Q
    M --> R["S07: Minimal fixes and validation"]
    A --> R
    Q --> R
```

## 结论成立所需的证据

| 结论类型 | 所需证据 |
| --- | --- |
| 程序缺陷导致某类异常 | 固定输入与候选，在明确函数/更新前后复现偏差；最小修改恢复期望关系；独立事件上对应异常率下降且效率与失败率可解释 |
| 材料、能损或谱混合解释结构 | 使用同一位置的 truth、匹配的输入谱和完整选择，forward simulation/transport 能量化重现结构；关闭相应效应的诊断实验按预期改变结果 |
| 旋转后电荷差异仍存在 | μ± 在共同支持域内的误判率差区间及分箱趋势；同一事件 before/after 的配对比较；重建效率与选择迁移一并给出 |
| 在目标精度内未见问题 | 预先定义指标、区间与等效性范围；给剩余效应的上限或可排除范围；不能用“没有显著性”代替“完全不存在” |

最终比较至少提供 `f(p_reco/p_truth<0.5)`、momentum core median/68% width、angle residual/peak fraction、μ± charge misID、truth reconstruction efficiency、duplicates/fakes、refit failure/fallback 及 CPU/I/O。主效应按事件聚类估计不确定度；逐层与多候选记录不作为独立统计样本。

每阶段同时交付[参考图清单](05_reference_figures.md)规定的分布图，使用参考材料对应的物理量和面板布局，展示before/after及必要的legacy/正确定义对照。PNG、逐图PDF、完整汇报PDF和bin counts写入results，主图嵌入workbook并链接图集。

## 下一轮落地任务

1. S05读取S04的paired particles与compact traces，冻结同一物理面和p窗口，分别看truth、seed、selected、refit/final的角度峰及选择迁移；使用A1–A5的显著性与forward-model验收。
2. S06使用correct material后的μ±共同支持域，检查首次sign flip与谱平衡后的率差；记录既有样本的精度，再确定追加规模。Alma历史selection与完整G0/G1仍单独控制。
3. S07从已定位的other-segment seed、真实radiative loss后的update/smoothing、chart/refit covariance、IFT策略分别做单因素实验。Gaussian mixture能损的文献可行性和典型事件见S04，未把它当作已验证修复。
4. 所有新批量按S04实测CPU/RSS/I/O配置资源，提交前重新比较兼容schedd；所有源码、数据身份与observer开关继续冻结。

每阶段先写 `workbook/YYYY-MM-DD_stageNN_topic.md`，列明事实、假设、配置、统计和产物路径，再提交代码/配置与报告。根因检验使用预先保留的验证事件，最终不以在调参样本上得到更好看的直方图作为验收。
