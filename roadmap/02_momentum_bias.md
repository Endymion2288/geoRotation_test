# 动量低估：逐运算与候选来源的因果研究

目标是分别解释 seed 自带的负尾和后续拟合新产生的负尾，并给出每类 track 首次异常的函数、surface、执行操作与总体占比。研究同时检查真实材料能损、truth 定义与程序实现，不能把低 p 相对生产点的 residual 自动归为重建偏差。

依赖：[材料与源码](00_materials_and_source.md)、[公共实验协议](01_common_protocol.md)。主样本为 A′ 110051、110033，控制为 100044 与 100116/100117 μ±。先 3-station forward，再 backward、4-station 和实际合并。

## M1：先复现分布并分清 truth 与绘图问题

1. 在原版配置保留 generator-referenced residual，复现报告 p24 的 seed 与 p16–18 的 per-layer predicted/filtered/smoothed 图，同时记录 event/track/state、单位与 truth 来源。
2. 逐事件核对 `p=1/|κ|`、Cartesian norm、truth charge/barcode 和 reference surface。检查无匹配 truth、0/负占位 p、q/p=0、invalid state、histogram under/overflow；有效正 p 的 `δp≥−1`，若小于 −1 则先查定义、单位或无效数据。
3. 增加同面 truth residual；将生产点→sensor 的真实 ΔE 与 reco−local-truth 分别画出。states writer 当前用 HepMC p，需要 S01 确认真值可用性；没有局部 p 的真实样本只报告 generator-relative 观察，另用 G4 audit 样本做物理归因。
4. 按正确/错误 charge、truth purity、e+/e−、μ+/μ−、station/hit pattern、生产点、opening angle/track separation 和 p_truth 分组；报告 tail 占比与 matched efficiency。A′ 同一事件两条 track 以 event 聚类估计误差。
5. 同时审计绘图读取的文件、tree、branch、筛选和复用对象，确认 before/after 统计相同是预期还是对象/文件混用。当前 circle fitting 仍在 F frame 做，旋转只改参数传递等环节，seed 标量 p 不随纯旋转改变本身可以是合理结果。

**交付**：原始定义/修正定义两套基线、选样 cutflow、局部 truth 可用性表、正常/tail/wrong-charge/未匹配的代表事件清单。验收为能将每一张图追溯到具体候选与 truth reference。

## M2：拆解 segment 和 seed 的初始估计

跟踪 `SegmentFitAlg → GhostBusters → CircleFitTrackSeedTool::Segment/Seed/fakeFit/get_params`：

| 要检查的运算 | 保存与比较 | 支持或排除的机制 |
| --- | --- | --- |
| cluster/space point→segment | strip/stereo side、local coordinate/C、truth 局部方向、segment slopes/C、shared hit 与 segment χ² | 前级测量、错误 stereo 组合、ghost removal 或方向估计偏差 |
| segment 组合 | 参与 station、真实 barcode、各 segment center/direction、fallback 层级与被剔除原因 | A′ 的 cross-track 混合、短 lever arm、候选先验与组合选择 |
| fake points | 每 segment 的 center 和 ±30 mm 方向点、重复数、与真实 space points 的偏离 | 人造点相关性、错误 segment direction、缺失 station 与拟合权重 |
| circle fitting | 实际 y–z 输入、radius/center、迭代停止、DET、近共线条件数、残差、getY 分支/radicand | 高 p 的小曲率不稳定、分支选择、p 接近零的异常回退 |
| radius→p/charge | 固定 B=0.55 T、单位、circle.cy 符号、p_perp 与 p 的关系；独立 sagitta/field-integral 估计 | 常量 B 近似、field nonuniformity、radius 单位/符号与 total/perpendicular p 混淆 |
| direction 与初始 plane | 两 station 连线/单 segment direction、投影 origin、F→A、局部坐标与 reference surface | 正反向 plane、θ/φ、position extrapolation 与坐标误用 |
| fallback 与排序 | 单 station 默认 p=9999999、charge=+1；stations 优先及 normalized χ² | charge asymmetry 与错误初值是否集中于 fallback；区分实际路径与未调用 fit() |
| covariance | 完整 C、单位、eigenvalues/correlations、inflation 分量 | 宽先验、局部基变换、κ–angle correlations 与后续 gain 的影响 |

使用相同物理几何下的已知直线/helix 点、带有限测量误差的点和真实 segment 做分级 closure；对 circle 方法加独立几何/数值解，只用于判别稳定性。单因素比较 fake points、真实 space points、线段夹角与真实场积分。变化前冻结 seed 组合与排序，变化后再放开选择，区分数值变化和候选替换。

高 p 的 curvature 接近零时同时检查 signed κ，不能用非线性 `1/|κ|` 的长尾外观直接评价估计器无偏性。对比两电荷的 seed charge confusion 与错误率，交给 Q 分支共同定位。

**交付**：seed tail 的机制分类、各类分子分母与曲率证据；每类给一组固定 hit 的可复现输入和替代估计结果。若 seed 均值正确而 CKF 出现新 tail，继续 M3–M5。

## M3：使用受控 seed 与固定 hit 分开初值和 CKF 问题

对同一 truth particle、初始 plane、event 与候选 lineage 设计以下阶梯。替换前后均保留有限且合理的 covariance；“零偏均值”不等于将 covariance 设为零。

| 配置 | 改变什么 | 回答什么 |
| --- | --- | --- |
| T0 | 正常重建 seed | 当前生产基线 |
| T1 | 仅筛选报告条件：`abs(δp)<0.1`、同 sign、`abs(δtheta)<10⁻⁴` | 重现 good-seed 现象；保留剩余 position/direction/φ/C 误差与条件选择偏差 |
| T2 | 在同一 seed 与 surface 只替换 signed κ 为局部 truth | 后续问题能否单由 momentum/charge 初值解释 |
| T3 | 替换完整 position、direction、κ/time 的 truth mean，先保留原 C，再使用经验证的同 frame C | 未替换参数、C 及坐标转换分别是否影响 CKF；truth/time 不可得项明确标注 |
| T4 | T3 + 固定正确 truth-associated hit 序列，用非组合 KF | 排除错误 hit association、branch/pruning/selection 后是否仍有偏差 |
| T5 | T4 + 理想几何/场、与输入生成一致的无材料/材料模型、已知测量噪声 | 验证 propagation、update、smoothing、units/C 的数值与统计 closure |

T3 使用局部 truth 并正确变换到初始 plane；若真实样本无法提供，先在有充分 truth 的诊断 simulation 上执行，再联系实际问题。T4 固定 hit 需检查 shared/mixed clusters 和二次粒子，记录纯度与 ambiguity，不能把所有与 barcode 有任意 deposit 的 cluster 当作纯测量。

现有 `TruthBasedInitialParameterTool` 不能直接视为精确 truth 注入：当前实现选 μ、默认 charge=+1、主动随机 smearing、使用 production vertex/p、`sigmaP` 取 vertex norm，且 covariance 填写有可疑项。应审计并实现独立、确定性、支持 e±/μ± 的诊断注入，记录替换后的实际 mean/C 与 reference plane。

**判别**：T2 仍异常而 T3 正常，优先查其他 seed 分量/C；T3 CKF 异常而 T4 正常，优先查 hit association/候选；T4 仍异常，查 transport/update/smoothing；理想 T5 失败则先解决闭合。各比较只在同一母样本上进行，另报因 T1 筛选引入的样本差异。

## M4：逐层定位 prediction、update 和 smoothing 的首次变化

对每个 seed/candidate 在**执行顺序**记录：

```text
seed
  -> propagation/material prediction at surface i
  -> measurement selection and branching
  -> filtered update at surface i
  -> next surface ...
  -> backward smoothing / target propagation
  -> candidate selection
  -> ACTS-to-Trk -> refit seed conversion
  -> refit prediction/update/smoothing
  -> final candidate / extrapolation / Ntuple
```

同时存 physical z_F 与 path length；backward finding 和 smoothing 的运算顺序与 detector layer 顺序不同。“最前层已坏”不代表错误发生在该层的 forward update。

| 变化 | 应记录的量 | 首要因果检查 |
| --- | --- | --- |
| seed→第一个 predicted | position/direction/κ/C，target plane，穿越 surface/B/material，F/Q 与 ΔE | frame/plane、过度材料、积分方向、错误几何导航和 seed 到首层传播 |
| previous filtered→next predicted | Δκ、path、B/gradient、material thickness、mean energy loss/straggling、transport J | 无材料应保持动量模长；材料开启时与同面真值/预期 ΔE 比较，检查重复穿层/材料 |
| predicted→filtered | innovation r、S、K、χ²、measurement dim、H、hit barcode、C 的 κ–position/angle correlations | 在本地 bound convention 下验证 `Δκ≈K_κ r`，记录 normalize 等额外非线性操作；错误测量还是错误 gain |
| filtered→smoothed | smoother gain、邻层 covariance/J、状态连接、首个引入大差值的执行操作 | smoothing 从后续错误层传播、reference surface/free-to-bound、矩阵条件数 |
| target fitted parameters | 返回 surface、smoothed/fitted state 对应关系、p/κ/C | output target 与预期相同；不是比较不同位置造成的能损差 |

先为所有事件聚合每一步的 `Δκ`、`Δp/p` 与 failure 分类，完整 r/S/K/F/Q 只保留按已冻结规则抽取的异常和匹配正常事件。对每种主 tail 类别给“正常进入该操作、异常离开该操作”的事件表。

单因素实验包括：CKF 的 material resolving 与 ACTS 实际能损/散射选项分别控制；精度/step size/maxSteps 收敛；真实场对理想场；measurement 原 C 对固定 C；seed C 的映射/inflation；outlier/measurement χ²、最大 branch 数；固定 hits 后的 smoothing 开关。仅 `resolveMaterial=False` 不自动代表所有材料效应已关闭，须记录实际选项与累计积分。

关闭 reconstruction 材料并不能让已包含 G4 能损的 RDO 成为“无材料 truth”。严格无材料 closure 使用匹配生成模型的 toy/诊断 simulation；实际 RDO 消融只测模型敏感性。

## M5：转换、refit 和选择审计

1. **ACTS→Trk→ACTS**：同一 state 保存位置、Cartesian momentum、charge、plane/curvilinear local basis 与完整 C。分别检查 1000 倍单位转换、local coordinate Jacobian、front/back 参数排序、附加 target-state 被标 Hole 的影响，以及 origin 重新向前移动 10 mm 时是否真的传播了 mean/C。
2. **refit**：对同一 CKF winner 比较 bypass/refit 两条诊断路径；记录 measurement IDs/顺序、每次更新、实际 covariance scaling、outliers、material、reverse filtering、成功与 fallback。当前两次 covariance inflation 分开开关，先对固定 hits 做 closure，再评估真实选择。
3. **粒子假设**：当前正常 seed、CKF 转换与 refit 多处固定 muon hypothesis。先验证 e± 与 μ± 在对应能区的传播/材料影响，再做相同 hits 下 electron/muon 假设的单因素比较。electron bremsstrahlung 与错误 particle mass 必须分开；已有 RDO 的 stochastic radiation 不能通过更换一个 hypothesis 完全补偿。
4. **候选选择**：在 shared-hit pruning、nMeasurements/χ² 排序、ghost removal、refit fallback 和 merge 前后保留全部候选 p/κ、hit purity 和 lineage；量化“好候选被淘汰”和“同候选参数变坏”。一个 candidate winner 改变不归类为一次 KF momentum update。
5. **最终输出**：对比 raw ACTS、Trk、final collection、Ntuple norm(p) 与各参考面 p；核对外推失败、raw frame 与 physical frame。将图使用的前级 CKF states 和最终 refit 状态明确区分。

**交付**：每个最终 tail track 的最后正常 stage、首个异常操作或 candidate switch、refit 状态、真值匹配和机制标签；按互斥主类别归类，未解释类单列。

## M6：验证、扩展与停止条件

根据发现的机制选最小修复，固定一项修改后重放全部诊断事件并在预留文件上批量验证。初始 toy p 网格可取 10、50、100、250、1000 GeV，两电荷及小斜率/近 beam 方向；只运行能区与当前问题相关的子集。真实样本覆盖 e± 两个 A′ 寿命/衰变拓扑、μ±、不同 station 模式；新 simulation 批量走 HTCondor。

验收输出：

- seed 与各后续 stage 的 median、68% width、`r_p<0.8/0.5/0.2`、wrong sign 与 failure rate，给事件数和区间。
- 每类 tail 的首次异常阶段占比、配对修复差值、未解释占比；truth-relative 与 generator-relative 分开。
- reconstruction efficiency、duplicates/fakes、track separation 依赖、angle structure 与 charge misID，避免通过丢弃难 track 获得更好的 momentum 图。
- 资源、有效输出 schema、完整复现命令、最小 Calypso patch 和有意义的 closure/regression checks。

若局部 truth 与经验证的能损模型已解释负尾，应将其报告为相应物理/定义效应并量化剩余偏差。若仍有无法解释的类，给其发生率、典型运算与下一项能区分假设的实验；不以少数事件或“峰变小”关闭分支。

## S04实施结果与后续入口

M1–M6阶段交付已完成，实际样本、所有控制、分子分母、配对区间、失败记录与复现入口见[S04报告](../workbook/2026-11-12_stage04_momentum.md)。严格旋转current22材料已通过独立样本验证；完整truth seed与纯hits控制证明错误材料图能在初值正确时产生大负尾。全量逐运算和局部truth界将seed、update、smoothing/refit及候选选择分别量化。

剩余研究由S07继续：同truth segment/方向与field-integral seed、真实辐射后的KF能损模型、measurement C/branch策略及chart/refit/IFT控制。23条固定真hits的电子local tail已有同sensor前后运算和真实沿程p证据，GSF文献与具体区分实验已列入报告；历史truth界待定类也保留发生率。S05/S06分别验收角度峰和谱匹配电荷差异。
