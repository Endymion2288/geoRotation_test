# 电荷误判与 μ+/μ− 差异

目标是重新测量旋转后 μ± 的 charge misID，区分旋转的改善、仍在的谱/接受度差异和实现偏差，并定位错误符号首先出现在 seed、CKF、conversion、refit 还是候选选择。

历史参照为 Alma 9 报告 p17–18：μ+ old/new 为 0.101%/0.197%，μ− 为 0.217%/0.359%；old-correct/new-wrong 为 0.122%/0.223%。这些表值缺少原始分母和误差，不能据此认定旋转后的误判率或统计显著性。

## Q1：严格定义分母、匹配与历史复现

1. 使用 100116（μ+、PDG=−13）、100117（μ−、PDG=+13）的 `s0012` RDO 和对应 `s0012-dev` / `s0012-r0019` PHYS 做历史配对。100044 的单 μ+ 5 mrad 样本作为额外控制，不能替代 Alma 的 μ± 2D FLUKA 对照。
2. 获取历史 Good Track 定义，精确复现 fiducial truth at vertex/三 station、`r_atMaxRadius<100 mm`、exactly one Good Track。与当前 Ntuple 判据分别给 cutflow，避免把历史 r<100 和当前 r<95、可配置质量阈值混合。
3. 主 matching 采用 hit/SDO/barcode，与 reco sign/p 无关；检查原分析按 PDG 与 charge 匹配的含义。若旧方法过滤了 wrong sign，保留旧结果复现并增加修正分析，报告两个口径差异。
4. 用 truth charge 与 signed κ 一致性检查 PDG 号、q/p 单位、track->charge()、最终 branch 和 forward/backward。κ=0、缺拟合参数、unmatched 与非法 charge 单列，不能计作“正确”。
5. 独立报告两个分母：

```text
conditional_misID = N(selected, valid matched track, q_reco*q_truth<0)
                  / N(selected, valid matched track)
truth_correct_eff = N(fiducial truth with selected correct matched track)
                  / N(fiducial truth)
truth_wrong_eff   = N(fiducial truth with selected wrong matched track)
                  / N(fiducial truth)
```

exactly-one-Good-Track 的历史条件与 inclusive matched-particle 条件分开。真值粒子找到多个候选时记录 duplicates/conflicts，并按固定、与 charge 无关的 matching/质量规则选 candidate。

**交付**：历史选择的可执行定义、样本交集/覆盖、匹配纯度、分子分母与 68%/95% 区间。历史数字无法精确复现时，给具体差异及待查输入，不以“版本不同”代替排查。

## Q2：测量旋转效应与正负差异

对同一 RDO 运行公共协议的 G0/G1/G2，先 3-station forward，再 backward、4-station 和实际 merge。历史 old/new comparison 与本次 before/after rotation 是两条不同的比较轴；不能将 historical dev 直接称作当前未旋转控制。

每种 charge、每个配置报告：

- global/分箱 `f_mis`、correct/wrong/unreconstructed/fake/duplicate/invalid 比例、refit failure/fallback、Good Track 通过率。
- 对同一 truth 的 before/after 状态转移表：correct→correct、correct→wrong、wrong→correct、wrong→wrong，另有失去/新增 track 与未匹配。
- `Δ_rotation=f_mis,after−f_mis,before`，按共同 truth 事件配对估计区间；conditional samples 的选择迁移另外展示。
- `Δ_charge=f_mis(μ−)−f_mis(μ+)` 与 ratio，给各电荷分母和区间。μ± 两个独立样本不按 event number 强行配对。

先按 p_truth（例如 <50、50–100、100–250、250–500、>500 GeV）、true position/slope 和 curvature significance 分组，按样本统计调整冻结边界。reco-p 分箱保留作与历史图对比，但不作为唯一 kinematic 轴。关注 high-p 小曲率区，也检查是否存在低 reco-p 的 wrong-sign 尾与 M 分支交集。

pilot 约 10k events/电荷仅用于配置、吞吐和粗分类；主测量起点约每电荷 100k **有效选中粒子**。以误判数、分箱有效 N 和目标差值精度调整追加作业；零计数给置信上限。置信区间包含事件内相关性，reweighting 后按 effective N 重新评估精度。

## Q3：排除 μ± 输入谱与接受度的混杂

FLUKA 正负 μ 的 p、位置、方向和 detector 照射分布可以不同；同一误差模型对两个不同谱积分也能产生不同总误判率。

| 实验 | 实施 | 解释 |
| --- | --- | --- |
| 原始谱对比 | p_truth、vertex/st1 x/y、true slopes、fiducial/hit pattern、field/material path、station 模式的 μ± 分布 | 明确差异是否已经存在于输入或接受度 |
| 共同支持域 | 选择两电荷都有统计的 truth kinematic 区域；比较分箱率 | 避免在只出现一个 charge 的区域推断算法不对称 |
| matching/reweighting | 对 truth p、x/y、slopes 等分层匹配或重加权；报告 balance、weight tail、sumw/sumw2 和 effective N | 测试相同输入谱下的 residual charge difference；不使用 reco p/sign 构造权重 |
| 选择迁移 | 全 truth 母样本、各自通过 Good Track、共同条件的 conditional 结果 | 区分 reconstruction/quality/几何接受度与电荷估计本身 |
| 理想对称 gun | 窄 p/方向、对称 x/y 与 slopes、两电荷相同生成条件 | 去除 FLUKA 谱差异后的算法检查；需要的 simulation 通过 HTCondor |

detector 命中模式可作为机制分组；若它是 charge bending/接受度造成的中介变量，分别报告“仅 truth 输入谱校正”和“再按 hit pattern 条件化”的结果，不能用后者抹掉待解释的几何效应。

匹配后差异缩小支持输入谱/接受度贡献；仍有差异则以分组和不确定度判断算法/几何机制。两种情况都继续定位错误事件本身的首次符号翻转。

## Q4：找到第一次错误符号和引发它的运算

共享 M 分支的 lineage 与 per-state trace，按执行顺序记录：

```text
truth -> segment/seed sign -> CKF predicted -> filtered -> smoothed
      -> candidate winner -> ACTS/Trk conversion
      -> refit -> final/merge -> Ntuple charge
```

每步保存 signed κ、σ_κ、`z_κ=κ/σ_κ`、truth purity、B/material、r/S/K/C、candidate ID、fallback 与 rejection reason。区分“同一 candidate 的 κ 穿零”“换成另一 candidate”“仅输出 sign 错误”；记录正确→错误以及错误→正确的修复路径。

| 优先机制 | 具体源码/实验 | 需要的证据 |
| --- | --- | --- |
| charge convention | PDG μ+=−13、μ−=13，field wrapper、beam direction、plane normals 与 propagation direction | 已知 q/B 的 helix 弯曲和所有 stage sign 一致；纯 proper rotation 不反号 |
| seed prior/fallback | `CircleFitTrackSeedTool` 的 circle.cy 符号、单 station charge=+1/p sentinel、fallback 频率 | 两电荷按 seed topology 的 confusion matrix；相同 hit 条件下先验与 CKF 纠正率 |
| hit/branch | mixed hits、outlier、holes、缺站、最大 branch/χ² 和 shared-hit pruning | 固定正确 hits 后错误是否消失；好 sign 候选是否被排序/淘汰 |
| measurement update | `Δκ≈K_κ r`、错误 hit、κ–angle/local covariance、极小曲率 C、smoothing | 固定 seed/hits，某次 update 或 smoother 使 sign 首次翻转；对称控制和最小修复 |
| conversion/refit | GeV/MeV、local/curvilinear C、初始 plane、两次 inflation、reverse filtering、material | 无传播 conversion sign 应精确保留；refit bypass 与逐步 trace 判别 |
| output/merge/matching | final collection 来源、Ntuple charge、truth association、重复候选 | raw fitted sign 对最终 branch，识别误标/换候选，不把 matching 错误当 charge fit 错误 |

优先重放历史 old-correct/new-wrong 的事件；报告 p17 可读到 run 100117 / event 1033，先用完整 production/event identity 核实再定位 entry，不能直接用 `--skip 1033` 代替物理 event ID。统计样本必须同时包含随机正常事件和所有定义内的错判类别，不能只解释一张 display。

## Q5：用对称性、独立样本和回归完成归因

1. **纯坐标旋转**：G0/G1 相同物理输入/先验，仅重新表达坐标；验证 sign 不变。出现差异时定位 parameterization、geometry/field 或 covariance 数值路径。
2. **q 与 B 同时反号**：在理想 transport/已知 hit 的条件下，Lorentz 项 qB 保持不变，同一轨迹的期望 κ sign 反号、绝对曲率与拟合质量相同。field provider、seed sign convention 与粒子假设一起一致变更；只改 q 或只改 B 会改变实际弯曲，不能称作同轨迹对照。
3. **空间对称性**：在 geometry 与 field 有相应对称性的位置比较镜像轨迹。若用反射 S，B 是轴矢量，应按 `B'=det(S) S B` 变换并一起处理 surfaces/material；不对真实非对称 apparatus 假定精确镜像相等。
4. **nominal/ideal geometry 与匹配 simulation**：区分真实 field/acceptance/alignment、散射能损与代码效应；带随机过程的样本只要求在统计误差与已知 charge 物理差异内一致。
5. **修复回归**：每次一项最小变化，在预留 μ± 文件与 A′ e± 上测 paired change；覆盖 forward/backward/3/4 station、不同 p 和 detector 区域。报告 charge 改善与 momentum/angle/efficiency/failure/CPU 的联合变化。

最终报告分别回答：旋转相对一致未旋转控制改善多少；μ−−μ+ 总差有多少来自谱与接受度；相同输入条件下还剩多少差异；各类错误首次产生在哪个操作。若不再观察到差异，给预先设定精度内的区间/上限；若仍存在但原因未完结，给未解释类别、比例和下一项可判别实验。

最终图组按Alma报告p16–18、p23–25画μ±的signed q/p、wrong-sign尾部和误判率，按p19–22、p27画truth/reco动量相关与宽范围动量谱；增加带区间的truth-p分箱misID、共同支持域/谱平衡前后对比。历史old/new、G0/G1旋转和current/strict材料比较分别标清，按[参考图规范](05_reference_figures.md)交付PNG、逐图PDF和汇报PDF。
