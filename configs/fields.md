# 诊断字段字典

本字典对应 `tracking_runner.py` 和 Calypso audit patch。定义及 cuts 的机器入口为 [definitions.json](definitions.json)。每个输出目录的 `run.json` 指定输入文件、entry 范围、conditions、material、安装路径和开关；同一目录的记录继承这些信息。

| 字段/记录 | 定义、单位与参考面 |
| --- | --- |
| `run,lumi,event,collection` | Athena EventContext 的生产标识和输出 collection；本批 RDO 在同一文件内会复用 run/lumi/event，不能独立用作事件身份 |
| `dsid,file_index,source_entry` | 冻结的 manifest 文件与 entry 构成稳定事件身份；粒子再加 barcode，诊断再加 collection/stage/seed/candidate/state；bootstrap 按此事件身份分组 |
| `stage,seed,candidate,state` | 运算阶段、实际 seed loop index、finder 返回时赋予的候选 ID、反向 state traversal 序号；`-1` 表示不适用 |
| `acts_state_index,acts_previous_index` | ACTS container 内的 state 与 predecessor；需连同 event、seed 和 stage 使用 |
| `event` 记录 | 每个输入事件恰好一条；seed、cluster、measurement 数及 `event_weight` |
| `truth` 记录 | HepMC barcode、PDG、status、signal 标志、生产 vertex F/mm、momentum F/MeV、各 station 的 SiHit 数/平均位置/平均方向；平均位置不是指定 sensor 面 |
| `truth_hit` | 只在完整 states 模式保存；实际 SiHit 中点 F/mm、方向 F、物理 sensor ID、station/layer；没有 sensor 处的精确动量大小 |
| `seed,seed_circle` | seed mean/C 与 cluster IDs；circle 中心/半径 mm、chi2、station 数；一 station 的默认 p/charge 可由 station 数识别 |
| `found,selected` | 所有返回候选与 greedy selection 后候选；hit 集、层数、chi2/ndof、holes/outliers、ACTS reference-surface mean/C；copy 保留 seed/candidate ID |
| `converted` | CreateTrkTrackTool 输出的 Trk 首个参数，A/mm、A/GeV；其 chart 是 curvilinear，`bound_trk`/`cov_trk` 的 q/p 单位是 1/MeV |
| `refit_input,refit` | 实际送入 refit 的 plane 参数与 covariance inflation 后 C、refit 输出；mean/C 参考各自记录中的 surface |
| `final` | refit 成功后的最终 TrackCollection，或 conversion fallback；`refit=1/0` 区分两种来源；参数是 Trk 首个参数所在面 |
| `bound` | ACTS `[loc0/mm,loc1/mm,phi/rad,theta/rad,qop/(1/GeV),time]`；`cov` 为对应完整 6×6 矩阵 |
| `variance_qop` | audit 中统一为 1/GeV²；从 Trk C 转出时乘 10⁶；不能与原 ROOT writer 的 raw q/p residual 混用 |
| `hits,matches,layers` | 原始 cluster compact IDs、各 barcode 的贡献 cluster 数、独立 station/layer 数；多个 truth 对同一 cluster 的贡献分别保留，和不必等于 measurement 数 |
| `*_state` | predicted/filtered/smoothed 及其完整 C、surface transform、Jacobian、measurement/H/R；只写实际存在的 component，缺省不是零 |
| `closure_*` | 合成状态、真实 sensor/field、实际 ACTS stepping 的闭合检查；独立于 physics event/candidate 的统计母样本 |
| `events.jsonl` | 独立 Python observer 的输入 truth/SiHit 计数和最终 Trk 参数；S01/S02 交叉核查使用，批量 `--compact-only` 省去此重复输出 |

JSON 非有限数写 `null`，失败另存 `finder_error`、`refit=0` 和作业 exit code。不能以缺行代替零 track 事件。生产 truth 的 continuation vertices 保留原 barcode；补充模拟的 `sensor_truth.csv` 使用 G4 event/track/parent/barcode 和物理 sensor 索引，与独立 run manifest 一起解释。

S03 compact particle 表保留所有 signal truth，包括未重建者；独立的各 stage best match 用于发生率，逐 stage 的首次数值变化分析必须通过原始 seed/candidate lineage 追踪同一条路径。
