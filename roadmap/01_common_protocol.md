# 公共实验协议与诊断数据

本协议用于让三个分支共享同一组物理定义、输入和候选来源。以下实验、字段和目标容差为待实现设计；已完成的环境检查另见阶段 00 报告。

## 物理坐标、参数与单位

使用 F 表示原始 FASER frame，A 表示旋转 ACTS frame：

```text
R = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
r_A = R r_F ; p_A = R p_F ; r_F = R^T r_A ; p_F = R^T p_A
(x_A, y_A, z_A) = (z_F, y_F, -x_F)
(px_F, py_F, pz_F) = (-pz_A, py_A, px_A)
s_y = tan(theta_y) = py_F / pz_F = py_A / px_A
theta_y = atan2(py_F, pz_F)   [forward physical direction]
r_transverse = hypot(x_F, y_F) = hypot(z_A, y_A)
kappa = q / |p| ; p = sqrt(px^2 + py^2 + pz^2)
delta_p = (p_reco - p_truth_at_same_surface) / p_truth_at_same_surface
```

`det(R)=+1`；旋转不改变电荷、`|p|`、物理方向和传输轨迹。ACTS 的 spherical theta 在 A 中相对 A-z 定义，不能当作探测器的 θy；ACTS 的 pT 也不能直接当作相对物理 beam 的 pT。Backward propagation 是沿轨迹反向积分，不自动意味着物理粒子动量反向或电荷反号。

ACTS 通常使用 mm、GeV、GeV⁻¹；Trk/HepMC 常使用 mm、MeV、MeV⁻¹。S01 从实际对象验证单位，输出表统一使用 mm、GeV、GeV⁻¹、rad，并保留原 branch 单位。`kappa_MeV = kappa_GeV/1000`，covariance 的 q/p 行列均要转换，方差有 10⁶ 比例。磁场 wrapper 有 kT→T 转换，必须与 field service 的返回单位一起验证。

协方差遵循完整参数映射 `C_new = J C_old J^T`；J 同时包含 surface local basis、角度与曲率的参数化及单位转换。不能把完整旋转简化为角度对调，也不能把 plane covariance 直接当作 curvilinear covariance。在旧 frame 接近 θ=0 时，优先在 Cartesian/free 参数或物理方向残差空间进行比较。

## Truth 与选择条件

每个残差明确标注 truth 的参考位置：production vertex、seed reference plane、单个 sensor surface、固定 st1 plane 或 VetoNu plane。不得将生成点 p 和局部 SiHit 方向拼合后命名为“同层真实参数”。

S01 逐类核查 RDO/HITS/PHYS 的 truth 信息。若 SiHit 只提供位置、方向、能量沉积和 HepMC link，不能据此唯一恢复该处 p。先报告 generator-referenced 与 same-surface 两种量的可用性；需要补充时在小型 G4 验证样本记录 pre/post-step 四动量、process、barcode 和穿面信息，或使用实际已有的等价 truth 记录。独立传播 production truth 只能作为模型预测，不能冒充经历实际 stochastic loss 的局部 truth。

主匹配使用 hit/SDO 的 barcode 贡献，保留 matched-hit purity、runner-up 贡献、parent/daughter、shared hits、unmatched 与 duplicates。预先固定工作点，例如多数匹配纯度≥0.7，并用 0.5/0.9 做系统比较；这些是拟采用值，S01 核对当前 hit 计数定义后冻结。匹配中不使用 reco charge 或 p 是否接近 truth。

样本选择分别定义：

- **truth 母样本**：稳定的 signal particle/μ、vertex 与三 tracker station 的 fiducial、实际 truth hit 要求；是否≥18 hits 按报告定义单列。以 truth 母样本测 reconstruction efficiency。
- **Alma 历史选择**：fiducial truth、`r_atMaxRadius<100 mm`、exactly one Good Track；取得原 Good Track 定义后精确复现。不将当前 Ntuple 的 r<95 mm 条件替代历史定义。
- **诊断选择**：所有候选、固定 best matched candidate、同一 truth 在多个配置中均通过选择的交集，以及各自通过选择的样本同时保留。
- **角度窗口**：全谱、历史 `Track_pz0` 窗口、物理 `50<p_reco<250 GeV`，并增加相同 `p_truth` 窗口；A′ 另复现 50–100 GeV。重建 p 切选的迁移必须单独量化。

报告中的 station 平均 truth position 另作历史量；精确残差统一映射到同一 sensor/plane，不能跨参考面计算位置、角度和能损残差。

## 对照软件与最小运行矩阵

先保留原版作为不可变的测量基线，再增加诊断补丁。建立 on/off 对照的优先方案是相同共同父版 `c95170f2` 加相同诊断，旋转版对应当前两次旋转修改；记录其余源码差异。需要在同一实现中开关时，将 geometry、alignment context、field/gradient、seed、reference surfaces、refit、writers 和 material bins 全部纳入统一配置，并做闭合验证。

| 组别 | Frame | 初始 covariance | 用途 |
| --- | --- | --- | --- |
| G0 | F | 旧版原参数 | 重现旧版行为 |
| G1 | A | 将 G0 的同一物理误差分布按完整 J 映射 | 分离 frame/数值效应 |
| G2 | A | 当前旋转版 tuned 参数 | 当前待测软件基线 |
| G3（需要时） | F | 将 G2 的物理误差反映射 | 检验调优与 frame 的交互 |

G0/G2 是“版本差异”；G0/G1 才检验保持物理先验一致时的坐标效应。非线性与 θ≈0 时用 sampled free-state covariance/方向空间验证映射，不在奇异点强行比较 ill-defined phi。

先用 **3-station forward** 一条路径完善诊断，再以相同 inputs 扩展 backward、4-station 和实际 merge。四 station 不默认代表所有 IFT layers 都启用：当前 CKF header 的 maskedLayers 默认是 `[0,1,2]`，需按实际实例配置记录；3-station 为 `[0,1,2,3,4,5]`。所有运行显式记录 masked layers，而不只记录“3/4 station”标签。

truth seed 实验、材料开关、measurement covariance、B 模型、refit 开关等仅在单个分支实验中变化，避免无目标的全参数笛卡尔积。普通数据对照保持相同 RDO、geometry/conditions、磁场、校准、随机种子与 event range。

## 必须先完成的闭合检查

| 检查 | 输入与判断 | 初始目标；在 S02 冻结最终阈值 |
| --- | --- | --- |
| 纯坐标变换 | 已知轴向/小斜率/两电荷 free states，RᵀR 与 forward/inverse | 位置差<10⁻⁸ mm、方向差<10⁻¹²、相对 p 差<10⁻¹²；charge 完全一致 |
| surfaces 与材料 | 对应 sensor 的中心/法向/local basis、sensor ID、binning/gluing；沿同一物理 ray 的 surface 顺序、X/X0、路径与能损积分 | 变换后中心/方向满足上项；无重复或漏 surface；积分差与容差一同报告 |
| field/gradient | 对应点验证 B_A(Rr)=RB_F(r)，gradient 与 finite difference；含边缘与三轴分量 | 先在平滑区域建立插值误差预算，FD 相对差目标 10⁻⁴；近零分量用绝对误差 |
| ACTS↔Trk | 同一 physical state 和 surface 的 mean、charge、完整 covariance round trip；analytic J 对 finite difference/MC | 无传播时 p 相对差目标<10⁻⁸；covariance 相对差目标<10⁻⁶，含 off-diagonal、对称性与正定性 |
| 传播 | 零场直线、均匀场 helix、真实场无材料、材料逐项打开 | 无材料相对 p 漂移目标<10⁻⁶，并对步长/精度收敛；material 开启时与预期 ΔE 和 Q 对应 |
| 观察是否影响算法 | 同事件 diagnostics off/on；仅修复输出定义/索引前后 | physics candidate、hit 集、拟合结果、计数一致；浮点顺序变化给可解释容差 |
| 事件状态 | 某事件独立运行、连续运行、倒序/分片、不同 CKF 实例与资源设置 | 同一 event 的 seed/hit/fit 结果一致；用于定位 static maps、未重置变量和跨实例状态 |

阈值是测试设计而非本次测得精度。全链 frame 对照若不一致，需要量化并定位数值差异；不能通过宽松容差掩盖研究对象。

## 诊断输出设计

每条记录包含稳定 key：`production_id, input_file_id, run, event, truth_barcode, collection, seed_id, candidate_id, parent_candidate_id, stage, state_id, physical_sensor_id`。输入文件与 entry 只作重放定位，配对物理事件时用 event/truth identity。A′ 的 e+、e− 与同一 parent 保留独立 key。

`trackGroup` 在当前实现中不能直接作为可靠 seed key。先在 finder 返回处用实际 loop index 建立 lineage，并检查 container copy、sort、refit 和 merge 后的归属。分支 splitting/pruning 保留父子关系，不按不同阶段的 vector index 强行 join。

| 数据级别 | 必存内容 | 大小控制 |
| --- | --- | --- |
| Run manifest | 软件 commit/patch、安装路径、ACTS/Athena 版本、flags、输入/event ranges、geometry/conditions/alignment/material/field、seed、threads、schedd/资源、schema version | 每 run 一份；worker 只记录 shard 差异，合并去重 |
| Event/particle | truth species/charge、生产与各参考面的有效 truth、fiducial/cut bits、输入/匹配/失败计数、event weight | 对完整母样本保留紧凑行 |
| Seed/candidate | segment/hit IDs、circle 数值与 fallback 原因、拟合/选择分数、p/κ/方向、covariance、候选来源及 rejection reason | pilot 保留全部；大样本保留标量与分类计数 |
| Per-state | surface transform、位置/动量、6D bound/free mean/C、F/Q/H、path/material/B、predicted/filtered/smoothed 的 κ 与角度、innovation/S/K、χ²、hit truth/purity、hole/outlier | 代表正常与各类异常事件保留完整矩阵；其余只聚合 step 指标 |
| Conversion/refit/output | mean/C 转换前后、reference surface、measurement 集与顺序、refit success/error/fallback、final collection/merge 来源、Ntuple 外推成败及 raw/physical frame | 与上游 candidate ID 明确关联 |

每个量有 frame、units、surface、validity/status。只绘制 `hasPredicted/Filtered/Smoothed` 为真的状态；missing、999999、NaN、q/p=0 与失败分别计数。最终 compact 表同时存 `p`、signed `κ`、方向与 charge，便于交叉检查。

现有 writers 可复用，但必须校验内容与单位。开启普通 seed/states/summary 需同时 `noDiagnostics=False`、`seedWriter=True`、`StatesWriter=True`、`SummaryWriter=True`；refit 的 `KalmanFitterTool` 及其 writer 也须独立打开，当前 CKF2Config 内硬编码的 `noDiagnostics=True` 需要覆盖。多实例为 tools、writers 和文件配置唯一名字；新配置作为 S01 的实际交付，本仓库当前尚未实现该运行入口。

## 异常定义与统计

采用可比较的初始指标：

- 动量：`r_p=p_reco/p_truth`、`δp`、`κ` residual/pull；主 tail 定义 `r_p<0.5`，另报 <0.8、<0.2、wrong sign、invalid 和 histogram under/overflow。先区分 production 与 same-surface truth。
- 角度：`s_y`、`theta_y`、同面 truth residual/pull；峰位置、面积、宽度、零值或离散值占比，以及条件分布的稳定性。主峰窗口按独立 pilot 与 angular resolution 冻结。
- 电荷：`f_mis=N(q_reco*q_truth<0)/N(valid matched selected tracks)`；另给按 truth 母样本的 correct/wrong/unreconstructed 比例，避免质量 cuts 通过丢掉难事件降低误判率。

率给分子/分母和 68%/95% Wilson 或 Clopper–Pearson 区间；零误判给上限。加权结果同时给原始计数、sumw/sumw2 与 effective N，重采样按 event/文件，不能套用未加权二项区间。多 track/多层按 event 聚类 bootstrap；paired comparisons 用事件配对转移表和 paired bootstrap/McNemar（适用时）处理相关性。μ+ 与 μ− 两份样本间的比较独立于 before/after 同事件配对。

报告 median、central 68% interval 与 tail mass，不以全分布 Gaussian fit 的 σ概括强尾部分布。角度 mixture 的多分量拟合只作描述；bootstrap、bin-width/offset 变化、truth forward model 和预留样本决定结构是否稳定。

pilot 先每类约 10k 事件用于发现路径与估算吞吐；charge 约 0.2%–0.36% 时，每电荷 100k **有效选中**粒子的理想独立统计误差约 0.014–0.019 个百分点（规划估算）。细分箱、reweighting 和稀有类别按有效 N 追加。主比较与精度目标在批量前冻结；探索后结论用独立事件验证。

## 计算、存储与执行

各阶段分布图按[参考图与汇报产物](05_reference_figures.md)执行：冻结对应参考页、样本、选择、bin edges、轴范围和归一化，输出before/after主图及PDF汇编，并把图中entries/tail与数值summary逐项核对。

所有运行输出位于 `/eos/home-x/xcheng/FASER/geoRotation_test_results`。S00 核查产物为 `stage00_material_review/`；后续目录为 `stageNN_<topic>/<run_id>/`。Git 保存配置、schema、分析和报告；EOS 保留紧凑表、合并统计、英文图表与少量重放事件。

S01 的单事件模板在确认 metadata 后选定 geometry；示意流程只使用现有 CLI：

```bash
source /eos/home-x/xcheng/FASER/geoRotation_calypso/setup_environment.sh
mkdir -p /eos/home-x/xcheng/FASER/geoRotation_test_results/stage01_inputs/smoke
cd /eos/home-x/xcheng/FASER/geoRotation_test_results/stage01_inputs/smoke
python3 /eos/home-x/xcheng/FASER/geoRotation_calypso/Control/CalypsoExample/Reconstruction/scripts/faser_reco.py --help > reco_help.txt 2>&1
```

实际 reconstruction 命令在 metadata 核对和诊断配置落地后写入 workbook，不把旧报告带 `-o` 的命令直接复制为当前可运行入口。

HTCondor 使用 `~/condor/geoRotation_test/<stage或任务名>/` 保存轻量 submit/wrapper。提交前 `condor_status -schedd` 和逐候选 `condor_q -name ... -totals`，记录时间、负载、任务兼容性与选择理由，使用所选 schedd。worker 在 scratch 内 setup 相同构建，读取固定 shard，重建并就地提取紧凑字段，仅传出必要结果；日志在 results。检查 EOS/CVMFS/凭据、EL9、wall time/memory/disk 条件。

开始时因当前 seed static maps 和 writer 状态尚未闭合，使用每进程单 event slot、多独立进程分片；与申请 CPU 对齐，并限制 BLAS/OpenMP/ROOT 线程，防止嵌套超额。线程安全通过 S02 后再比较 Athena 多线程与多进程吞吐。ACTS 默认 CPU 工作流不申请无用途 GPU。

每次分片使用固定输出名与唯一 run/shard ID，失败重试按输入范围去重。合并时核对输入覆盖、重复 key、各 stage 计数、sumw/sumw2、失败清单与 cutflow；保留重试来源。先完成合并与抽样读回，再清理中间 xAOD 和完整 states。
