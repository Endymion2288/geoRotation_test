# 阶段 02：旋转、接口与诊断闭合

日期：2026-11-11。192 个 sensor 的旋转位置和 basis 与 GeoModel 完全一致；无材料 ACTS 传播的动量模长保持不变。10 个 μ+ 事件的原版、审计版与 ROOT diagnostics 开关对照，最终 track 数和记录参数完全相同。实际接口测试确认原版 plane→curvilinear covariance 缺少 chart 变换；实验修正已通过有限差分，但其对总体异常率的影响属于后续因果研究。

## 阶段目标

验证诊断能可靠追踪 seed、候选、CKF states、conversion、refit 与 fallback；闭合 frame、units、实际 sensor、field/gradient、material 和 covariance。保留原版物理基线，用独立开关测试接口修正与 event maps，不将这些实验结果混入原版发生率。

## 工作内容

`patches/0002-tracking-audit.patch` 以 Calypso `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` 为基底，增加 C++ opt-in JSONL observer。每行包含 schema=2、run/lumi/event、collection、stage、seed/candidate/state；同一物理事件的文件与 entry 由运行 manifest 补全。seed 使用实际循环下标，finder 返回时赋予 candidate ID，选择和复制保留它；refit 通过 `fitAudited` 接受显式 lineage，原 `fit` ABI 保留。

完整模式写出 predicted/filtered/smoothed mean/C、transport Jacobian、measurement/H/R、sensor transform 和 SiHit 方向；紧凑模式保留 stage mean、q/p variance、hit/SDO 来源、拟合质量、refit/fallback 和 truth。SDO 匹配不使用电荷或动量。闭合检查调用实际 ACTS/Trk converter、refit 参数接口、field wrapper 与 ACTS propagator；不是另写一个替代拟合器验证自身。

`CorrectCovariance`、`CorrectInterfaces`、`ResetEventMaps` 均默认关闭。接口实验将 plane→curvilinear 的 covariance 变换为 `J C Jᵀ`，反向使用同点投影 Jacobian；refit 初始 plane 放在同一物理 x。原有两次 covariance inflation 留作独立实验，未在此处重新调参。新增接口的首版存在 Eigen block 引用临时 rotation matrix 的问题，造成 covariance 往返相对误差约 1.416；改为拥有数据的 3×2 矩阵后，冻结源码重新构建并完成下述验证。失败版本未用于 pilot。

编译只更新修改的 tracking target；Configurables、confdb 和安装同步，日志置于 results。最终构建日志为 `build_frozen.log`、`build_final.log` 和对应 install 日志；配置生成方法见 `patches/README.md`。最初构建过程中修改头文件曾造成 schema 不一致，已通过冻结后强制重编 CKF、converter、fitter 消除。批量使用独立软件快照，不读取正在编辑的配置或组件库。

## 详细数据

### 诊断完整性与观察开关

样本为 100044 首文件 entries 0–9，FASERNU-04 / OFLCOND-FASER-05，3-station forward，masked `[0..5]`，历史 rotated material。原版对照加载修改前保存的组件库。

| 配对 / 检查 | 实测结果 |
| --- | --- |
| 原库 vs 新审计库，最终 observer | 10 events / 10 tracks；全部记录字段 normalized difference=0 |
| ROOT diagnostics off vs on | 10 events / 10 tracks；difference=0 |
| 原版连续运行 vs ResetEventMaps=true | 10 events / 10 tracks；difference=0 |
| 异常 event 1 独立 vs 连续 / reset | 各 1 event / 1 track；difference=0 |
| 正常 event 2 独立 vs 连续 / reset | 各 1 event / 1 track；difference=0 |
| entries 5–9 独立分片 vs 连续运行 | 5 events / 5 tracks；difference=0，77 个诊断 key 唯一 |
| A′ 110051 的原库 vs audit | 10 events / 21 tracks；difference=0；19967 个诊断 key 唯一 |
| schema=2 audit 唯一 key | 1488 行，0 重复 |
| seed / found / selected / converted / final | 10 / 22 / 10 / 10 / 10 |
| refit_input / refit / final fallback | 9 / 9 / 1 |
| found / selected / refit state 行数 | 307 / 248 / 231 |
| SiHit / truth 记录数 | 471 / 80，后者包含 continuation 与次级粒子，不等于 signal 分母 |

所有选中/final candidate 均能返回原始 found 和 seed。S01 的 observer 尚无 lumi 列，原库比较仅对这一个固定输入文件显式使用 `--legacy-observer`；C++ audit 的 lumi/schema 不作缺省处理。独立与连续结果相同只支持上述事件范围的 event 一致性，未测量多 slot / 多线程的线程安全。

### 几何、场和传播

| 实际测试 | 实测值 |
| --- | --- |
| 192 sensors 的最大 center / basis 差 | 0 mm / 0 |
| 24 个点的 `B_A−R B_F` | 0 T |
| gradient 与 `R G_F Rᵀ` | 0 T/mm |
| 24 个 grid 点 + 24 个偏移点，各 3 个差分步长 | 共 144 trials；每组 3 个零梯度 trial，用绝对误差表示 |
| 偏移点 analytic gradient vs native double interpolation gradient | 最大相对误差 1.13443×10⁻⁷ |
| 偏移点 native double central FD | 最大相对误差 5.56317×10⁻⁹，最大绝对误差 4.62241×10⁻¹³ T/mm |
| grid 点 native double 最接近单侧 FD | 最大相对误差 6.19642×10⁻⁹ |
| float cache vs native double field | 最大差 7.08730×10⁻⁸ T |
| zero / uniform / real field，μ±，step 100 / 10 mm | 12/12 propagation 成功；最大相对 p 漂移=0 |
| zero field，2 m 直线终点误差 | 6.82123×10⁻¹² mm |
| 0.55 T uniform field，2 m helix 解析终点误差 | 3.40091×10⁻¹¹ mm |
| real field，100 / 10 mm 步长终点差，q=−1 / +1 | 6.91939×10⁻⁶ / 3.91913×10⁻⁶ mm |

测试使用实际 conditions 的 `FaserFieldTable_v2.root`、scale=1，路径为 `/cvmfs/faser.cern.ch/repo/sw/software/22.0/faser/offline/ReleaseData/v20/MagneticFieldMaps/FaserFieldTable_v2.root`。field 查询点 F/mm 为 `x∈{−50,0,50},y=17,z∈{−1700,−1000,0,500,1000,1500,2000,2400}`；偏移点加 `(0.371,0.293,0.419) mm`，差分步长为 0.1、0.01、0.001 mm。

初次对 float cache 作 h=0.01 mm 中心差分时，相对差最大 0.3833。native `BFieldCache` 的插值权重、field corners 与 arithmetic 使用 float，小 h 的 field 可因舍入相等；map 的分段三线性插值在 grid 边界的左右导数也不同。随后调用同一 map 的 `BFieldMesh::getB` double 接口作独立检查，内部 central FD 和边界单侧 FD 得到上表结果。因此这次差异有数值解释，未发现梯度旋转错误。该检查未量化 float 场插值对 track tail 的总体影响。

### conversion 与 refit 起点

合成初态使用实际 beam-normal plane，p=150 GeV、q=±1、phi=−0.003/0/0.004 rad、theta=π/2−0.002 rad，非对角正定 covariance；两电荷共 6 trials。有限差分 chart 使用 loc=10⁻³ mm、angle=10⁻⁵ rad、q/p=10⁻⁶ GeV⁻¹。

| 指标 | 原版 | 接口修正实验 |
| --- | --- | --- |
| mean position 往返误差 | 0 mm | 0 mm |
| mean momentum 往返相对误差 | 1.89479×10⁻¹⁶ | 1.89479×10⁻¹⁶ |
| covariance chart 有限差分 normalized 最大差 | 0.105495 | 2.02218×10⁻¹¹ |
| 除去配置 inflation 的 C 往返相对误差 | 1.37064×10⁻¹⁶ | 5.48258×10⁻¹⁶ |

原版 C 往返闭合只说明两个接口重复复制同一矩阵；有限差分显示第一次复制与输出 curvilinear chart 不一致。修正需同时满足 chart FD 与往返，不能只检查 C 往返。检查范围是前五个 tracking 参数；time 不进入 Trk 5-parameter fit。

实际 10 事件中，selected→converted 的最大 p 相对变化和 q/p variance 相对变化均为 2.22045×10⁻¹⁶，说明这次 conversion 的低 p 不是 mean 的单位变换造成。原版 9 个实际 refit_input 均比 converted 起点平移 10 mm，没有伴随 mean propagation；接口实验将其作为独立可切换问题。它对三个总体问题的影响尚未测量。

### 材料与数值预算

当前默认未旋转 map 有 22 surfaces；历史 rotated map 有 12 surfaces，部分中心还相差 2 mm。因此不能把两者对比等同于纯旋转。`scripts/rotate_material.py` 从冻结的当前 map 生成 22-surface 旋转版本，每 surface 4096 points，bin data 和 geometry IDs 完全相同，变换前后 local-coordinate 最大差=0 mm。

2026-11-12逐项补查：AFS旋转图与strict共同12项的86400个材料bin完全相同；strict额外10项的118个bin全为零厚度空材料。2 mm差沿binUtility局部法向，不改变共同项的binX/binY查询，不能将12/22解释为缺失10层非零材料。具体数据与已有Calypso小样本差异见[专项比较](2026-11-12_stage04_material_maps.md)。

实际 Navigator+MaterialInteractor 用 p=150 GeV、2.4 m path、max step=20 mm，在历史 map 和一致旋转 22-surface map 分别测试 energy loss / scattering 的四组合、两电荷；各配置均 8/8 成功。这条路径两种 map 得到下表相同结果，不代表所有轨迹的材料分布相同。

| q | X/X0 | interactions / unique surfaces | loss off 的 Δp | loss on 的 Δp |
| --- | --- | --- | --- | --- |
| −1 | 0.0823869081 | 7 / 7 | −6.694719 keV | −7.471543 MeV |
| +1 | 0.0822111079 | 7 / 7 | −6.694719 keV | −7.473543 MeV |

在 deterministic ACTS MaterialInteractor 中，开启 scattering 的 mean p 与同 energyLoss 配置一致；散射角 covariance 预算本阶段未测量。loss on 的累计 sigma(q/p)² 为 1.23621×10⁻¹⁶ / 1.23580×10⁻¹⁶ GeV⁻²，loss off 为 0。

`PointwiseMaterialInteraction` 的 momentum/mass/qOverP 为 float；即使 Eloss=0，也计算 `hypot(m,p)` 后再开方还原 p。按该表达式预测一次 drift=−6.694725556 keV，与实测差 6.22435×10⁻⁶ keV。`MaterialInteraction.deltaP` 相对每次 rounded p 记录，7 次相加为 −46.863035 keV，而最终 p 只保留一次同样 rounding 后的值。所有组合的总 Δp 与 ΣdeltaP 最大差 40.168316 keV；按实际 7 次 interaction、p 与 q/p 各半个 float ULP 累加的界为 90.076588 keV，实测在界内。不能将该日志差额解释为丢失了几十 GeV 的能量。

### 正常与异常事件的范围

100044 event 1：production truth p=419.987668 GeV，seed0 p=828.622830 GeV；selected candidate4 p=6.954013 GeV，6 个 measurement clusters 中 4 个贡献来自 signal barcode10001，purity=2/3，未达到冻结的 0.7 matching 阈值；conversion 保持 p，refit 因 measurement 数不足返回 fallback，final 比值 0.016558。它也不通过 Good95。该例支持 CKF/候选路径已有低 p 输出，**不能作为合格 signal track tail 的个例或发生率证据**。

100044 event 2：production p=293.267169 GeV，seed0=344.813521、selected candidate0=296.504163、ACTS refit=295.228344、Trk final=295.223848 GeV；18/18 signal clusters，final 通过 Good95，比值 1.006672。refit stage 和 final 的参考位置不同，不能把二者 p 差直接当作一次无传播 conversion 的误差。

### 命令、产物与验证

初始化环境后在各 results 子目录执行：

```bash
python3 .../configs/tracking_runner.py --dsid 100044 --events 10 --audit --states --closure
python3 .../configs/tracking_runner.py --dsid 100044 --events 10 --audit --states --diagnostics
python3 .../configs/tracking_runner.py --dsid 100044 --events 10 --audit --states --closure --correct-interfaces --material .../stage02_closure/material-current-rotated.json
python3 .../scripts/check_field.py --output .../stage02_closure/field_interpolation_final.json
python3 analysis/check_closure.py --run <run_dir> --paired <paired_dir> --output <checks.json>
python3 -m unittest discover -s tests -v
```

results 根目录为 `/eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/`。正式小样本目录为 `sealed_original_100044/`、`sealed_diagnostics_100044/`、`sealed_corrected_100044/`；原库在 `original_lib/`。汇总为 `original_audit_checks.json`、`corrected_checks.json`、`diagnostics_pair.json`、`event_independence.json`、`small_event_summary.json`、`field_interpolation_final.json` 和 `material_rotation.json`。运行目录保存日志、resource、manifest 和完整 states。带旧 schema 或错误接口实验的早期目录只作调试记录，不计入验收数值。

Python 测试覆盖电荷/动量独立 matching、旋转 slope、未重建 truth 分母和 status=1 筛选，2/2 通过；Python AST、bash syntax 与 `git diff --check` 通过。外部 patch 对当前源码的 reverse apply check 通过。

## 阶段性结论

**事实**：旋转、GeV/MeV、实际 surface/field、无材料 transport 和诊断配对在上述范围闭合。材料漂移与导数差分中的小量已得到源码和数值证据支持的解释。原版 covariance chart 与 refit 起点存在可重现的接口问题，已交付可关闭的实验修正；没有据此声称已解释总体 momentum tail、剩余 slope peaks 或 charge asymmetry。

**事实**：原版首批异常 event1 不满足主 purity 和 Good95，说明后续必须同时报告 cutflow、匹配失败和 fallback。静态 event maps 未在这 10 事件及代表事件对照中改变输出；多线程安全未验证。

**待检验**：原版 C/plane 接口是否造成合格 track 的偏差、seed/CKF 的主要低估类别、不同 material 的全样本影响和正负电荷差异，交由 S03 基线及 S04–S06 单因素实验回答。

## 下一步计划

S03 已冻结 50 个 1000-event shard，总计 50,000 events，通过 HTCondor 执行原版物理开关关闭的 baseline。首轮 cluster1573723 在 setup 阶段 exit127，修正 worker 对 CVMFS 初始化的 shell 模式后，单事件 cluster1573728 exit0；正式重提 cluster1573730。完成后核对覆盖和来源、合并紧凑表，报告 selection、tail、slope、misID、区间及 CPU/I/O；各 root-cause 分支继续使用稳定 lineage 重放代表事件。
