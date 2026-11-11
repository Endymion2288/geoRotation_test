# 阶段 01：输入、定义与可运行诊断配置

日期：2026-11-11。五类样本均完成各 10 事件的完整性重建，共 50 事件、70 个 generated signal 粒子、65 条最终 track。重建作业均正常返回；无 track 的事件保留在分母中。sensor 处精确 truth momentum 的补充记录已取得，补充模拟在程序结束时存在单独的退出故障，详见下表。

## 阶段目标

冻结输入文件与 event ranges、production geometry/conditions、坐标和单位、truth 与 Good Track 定义；交付能观测 seed/CKF/refit 的小规模入口。输入验收要求每类 10–100 事件、处理事件数与预期一致、零 track 事件不丢失。

## 工作内容

`configs/samples.json` 固定 5 类样本首批 10,000 事件的具体文件。`inspect_inputs.py` 读取每类首文件完整 POOL metadata 和 CollectionTree schema；`tracking_runner.py` 仅运行 clusterization、space points、segment fit、GhostBusters、3-station forward CKF、refit 和 truth conversion，省去 waveform/calo 与长期 xAOD。

入口分别配置 CKF 与 refit writers，refit writer 使用唯一名字和输出文件；普通 `--diagnostics` 的 seed/summary、`--states` 的完整状态及独立最终 track observer 均有明确开关。字段和统计定义集中在 `configs/definitions.json`、`configs/fields.md`。

## 详细数据

### metadata 与 schema

共同输入根目录 `/eos/experiment/faser/data0/sim/mc24/`。完整文件名和文件大小见 `configs/samples.json`，完整原始 metadata 保存于 results。

| DSID | 粒子/来源 | RDO tag | 首文件 entries | pilot 输入文件数 | 生产 Athena | 顶层 branches |
| --- | --- | --- | --- | --- | --- | --- |
| 100044 | μ+、5 mrad FLUKA-E | s0013 | 500000 | 1 | 22.0.49 | 16 |
| 100116 | μ+、2D FLUKA | s0012 | 100000 | 1 | 22.0.49 | 16 |
| 100117 | μ−、2D FLUKA | s0012 | 100000 | 1 | 22.0.49 | 16 |
| 110033 | A′→e+e−，10 MeV、ε=10⁻⁵ | master_dev | 2000 | 5 | 24.0.41 | 20 |
| 110051 | A′→e+e−，10 MeV、ε=10⁻⁴ | master_dev | 1000 | 10 | 24.0.41 | 20 |

五类均为 **GeoFaser=FASERNU-04、IOVDbGlobalTag=OFLCOND-FASER-04**。μ 样本使用 McEventCollection_p5 / EventAuxInfo_v2，A′ 使用 p6 / v3；共同提供 `TruthEvent`、`BeamTruthEvent`、`SCT_Hits`、`SCT_SDO_Map`、`SCT_RDOs`。这排除了把 TI12MC03 自动映射到 FASERNU-03 当作生产配置的做法。

### 50 事件完整性重建

此组沿用历史旋转重建配置：FASERNU-04 / OFLCOND-FASER-05、历史 rotated material map、masked layers `[0,1,2,3,4,5]`。它与生产 conditions 的差异显式保留；S03 的当前入口默认材料基线与该组分开。

| DSID | 事件 | generated signal | SiHits 总数 | 最终 track | 零 track 事件 | exit code | 总 wall time/s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 100044 | 10 | 10 | 471 | 10 | 0 | 0 | 179.906 |
| 100116 | 10 | 10 | 528 | 9 | 1 | 0 | 126.986 |
| 100117 | 10 | 10 | 720 | 10 | 0 | 0 | 106.648 |
| 110033 | 10 | 20 | 3820 | 15 | 1 | 0 | 103.484 |
| 110051 | 10 | 20 | 1378 | 21 | 0 | 0 | 123.719 |

wall time 包括 POOL metadata、geometry/conditions、Python/ROOT dictionary 初始化，不能将 `10/wall time` 当作稳定事件吞吐。A′ 的 track 数可低于或高于 signal 粒子数，需依赖 SDO matching 解释漏重建、额外候选和次级粒子。

100116 的第一次配置调试在独立 observer 调用不存在的 `TruthParticle.p()` 时失败，已改为从 px/py/pz 计算模长，并完整重跑 10 事件。失败日志保留为 `original_100116/`，成功输出为 `original_100116_v2/`。没有把失败作业的部分输出计入本表。

### 冻结定义与 truth 可用性

- proper rotation `r_A=(z_F,y_F,-x_F)`，动量同样旋转，charge 不变；物理 `tan(theta_y)=py_A/px_A`，角度为其 atan。legacy `atan(py_A/pz_A)` 单独复现。
- ACTS mean 的动量使用 GeV、q/p 使用 GeV⁻¹；Trk/生产 truth 使用 MeV、q/p 使用 MeV⁻¹；位置 mm、角度 rad。每次跨接口显式转换 C 的 q/p 行列。
- generated signal 使用 status=1、`0<barcode<200000` 和目标 PDG；G4 continuation vertices 另存。fiducial 使用 vertex 半径<100 mm、station 1–3 各自有效 SiHit 平均位置半径<100 mm；H18 要求这三站至少 18 SiHits。
- 当前 public Ntuple 源码的 Good95：χ²/ndof<25、≥7 layers、ndof+5≥12、拟合 state 最大半径<95 mm。Alma PDF 的 r<100 / exactly one Good Track 与未恢复的 Twiki selection 作为历史复现差异保留，不能称为已精确复现。
- 主 matching 使用 cluster RDO 对应的 SDO barcode，purity≥0.7，charge 与动量不参与匹配。多贡献和未匹配分别保留，best candidate 按 hit 数/chi2/稳定 ID 选取。
- `FaserSiHit` 保存 start/end、edep、time、particle link，**不保存 sensor 处动量大小**。现有 `FiducialParticleTool::getMomentum()` 使用局部 hit 方向乘 HepMC production p。RDO 的部分 continuation vertices 确实保存沿途 p snapshot，但它们不是每个 sensor 的精确 p。

已落实的补充入口为 `configs/local_truth_sim.py` 和 FaserSctSensorSD 的可选 sidecar：记录同一步的 pre/post G4 momentum、F-frame position、G4 event/track/parent、barcode 和 sensor 索引。后续若对历史样本补充，需要使用其 BeamTruthEvent、随机状态及生产配置重新模拟；新模拟单独标识，不声称还原了历史随机相互作用。

### 同一步局部动量补充验证

单事件 gun 使用 μ+、100000 MeV、`z=-10 mm, theta=0.002 rad, phi=0`、seed=80217。原生 `G4FaserAlgConfig_Test.py` 入口分别启用与关闭记录器；patch 保存于 `patches/0001-sct-local-truth.patch`。

| 实测项目 | 结果 |
| --- | --- |
| sidecar 行数 | 59 |
| primary G4 steps / unique sensor 索引 | 29 / 25 |
| primary 覆盖 station | 1、2、3 |
| primary pre-step p 范围 | 99995.898060–99999.857745 MeV |
| 启用 / 关闭记录器的 HITS entries | 1 / 1 |
| 两个 HITS 的顶层 branches / ROOT zombie | 19 / false |
| 启用 / 关闭记录器的进程 exit code | 134 / 134 |

两种配置均处理事件、完成 HITS 写入和 Athena finalize，随后报告 `free(): invalid pointer`。改用原生 exec 入口后仍相同，关闭 sidecar 的对照也相同；退出故障与记录器的启用开关无关，根因尚未定位。本阶段证据支持局部 p 字段能实际写出、HITS 可打开并有一个事件；补充模拟尚未通过正常进程退出的验收。未将它作为历史样本的局部真值，后续大批量重模拟前须解决退出故障。

执行命令为 `python3 .../configs/local_truth_sim.py` 和 `python3 .../configs/local_truth_sim.py --no-sidecar`，分别在 results 的 `stage01_inputs/local_truth_sim_v2/` 与 `local_truth_no_sidecar/` 中运行。精确数值汇总在 `stage01_inputs/local_truth_validation.json`；pre/post CSV、配置、日志和小型 HITS 位于对应目录。

### 执行与产物

在已初始化的 Calypso 环境执行：

```bash
python3 scripts/inspect_inputs.py --samples configs/samples.json --output /eos/home-x/xcheng/FASER/geoRotation_test_results/stage01_inputs/input_metadata.json
python3 /eos/home-x/xcheng/FASER/geoRotation_test/configs/tracking_runner.py --dsid 100116 --events 10 --diagnostics
```

每个重建在独立 results 子目录执行，用 `/usr/bin/time -v -o resource.txt` 记录资源。原版组件库保存在 `stage02_closure/original_lib/libFaserActsKalmanFilter.so`；其余完整性作业通过 LD_LIBRARY_PATH 优先加载此库，避免诊断编译覆盖基线。

产物根目录 `/eos/home-x/xcheng/FASER/geoRotation_test_results/stage01_inputs/`：`input_metadata.json`、`integrity_summary.json`、`metadata.log`、`original_<DSID>/`；100116 成功目录带 `_v2`。ROOT diagnostics 只对 100116 小样本打开，其他样本保存紧凑 event 表。

## 阶段性结论

**事实**：五类输入在现成旋转版构建上可读取和重建；50 个输入事件均有对应输出记录，两个无 track 事件保留。production geometry/conditions 已读实值。三种 material 配置不是同一文件的简单别名，须独立比较。

**事实**：已在少量 μ 事件中再次看到低 p 输出，例如 100044 event 1 的唯一 track `p_reco/p_production_truth≈0.0166`。此处尚未通过 SDO 审计确认其完整来源，不用该个例估算总体 tail 率或指定根因。

**待检验**：production→tracker 的真实能损、seed 初值、CKF hit 选择、conversion/refit 和输出定义均可能参与偏差。S02/S03 将通过稳定 lineage 和配对测试区分。

## 下一步计划

S02 完成实际 sensor/field/gradient/ACTS propagation、conversion covariance 的数值闭合，比较 diagnostics 开关及独立/连续事件。S03 通过 HTCondor 执行冻结的 50 个 shard、共 50,000 事件，测量三问题的基线与资源消耗。
