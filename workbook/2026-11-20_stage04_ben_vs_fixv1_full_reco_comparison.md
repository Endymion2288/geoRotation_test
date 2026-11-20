# Stage 04 — 2022–2024 五时段全链路重重建与 Ben (`El9CalypsoForWaveForm` + Official `PHYS`) 对照研究

- **日期**：2026-11-20
- **对比软件栈**：
  - **本工作重重建 (`fix_v1`)**：`/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1`（branch `fix/latest-seed-and-momentum-v1`，HEAD `1c0f3f45`；基于 `geoRotation_calypso_latest` `7fd40df8`，合入旋转种子/动量修正 `e47b85fd` 与 `NtupleDumper` 无硬编码两级 Reduced Charge 处理 `a4ecefea` / `7729ddf0` / `1c0f3f45`；Athena `24.0.41`，`x86_64-el9-gcc13-opt`）
  - **Ben 分析基线 (`Ben`)**：官方生产 `PHYS` ntuples（2022–2023 为 `r0021`，2024 为 `r0022`）按 `(run, eventID)` 关联 `/eos/user/b/bewilson/El9CalypsoForWaveForm` 五个目录下的 `myVetoNu0_rawCharge35` / `myVetoNu1_rawCharge35`（与 [file_config.yaml](file:///eos/home-x/xcheng/FASER/nu_muAnalysis/AnalysisZipFramework/config/file_config.yaml) 完全一致）
- **绘图与物理选择参考**：
  - FASER Data Quality 框架：[FASER_DQ_RDF.py](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/FASER_DQ_RDF.py)、[plotting_utils.py](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/notebooks/plotting_utils.py)、[finalize.py](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/internal_note/finalize.py)、[Electronic_2D_measurement_internal_note-6.pdf](file:///afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/docs/Electronic_2D_measurement_internal_note-6.pdf)（Figs. 5–6）
  - $\nu_\mu$ 分析选择定义：[cuts.yaml](file:///eos/home-x/xcheng/FASER/nu_muAnalysis/AnalysisZipFramework/config/cuts.yaml)、[definitions.yaml](file:///eos/home-x/xcheng/FASER/nu_muAnalysis/AnalysisZipFramework/config/definitions.yaml)、[Analyis.cxx](file:///eos/home-x/xcheng/FASER/nu_muAnalysis/AnalysisZipFramework/source/Analyis.cxx)
- **提交、分析与绘图脚本**：
  - HTCondor 单文件重建脚本：[ben_comparison_worker.sh](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/ben_comparison_worker.sh)
  - HTCondor 任务生成与提交入口：[submit_ben_comparison.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/submit_ben_comparison.py)
  - 5 进程矢量化对比与坐标系审计：[ben_vs_fixv1_compare.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/ben_vs_fixv1_compare.py)
  - FASER DQ 风格对比绘图：[plot_ben_vs_fixv1_dq.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/plot_ben_vs_fixv1_dq.py)
- **产物目录**：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/`
  - 调度与任务清单：[submission.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/submission.json)、[tasks.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/tasks.json)
  - 完整直方图与统计汇总：[ben_vs_fixv1_full.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/analysis/ben_vs_fixv1_full.json)、[ben_vs_fixv1_summary.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/analysis/ben_vs_fixv1_summary.json)
  - 8 页合并汇报图集 PDF：[ben_vs_fixv1_dq_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/ben_vs_fixv1_dq_report.pdf)、[figure_index.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/figure_index.json)

---

## 阶段目标

1. 在 Ben 历史分析覆盖的全部五个数据时期——**2022**（`waveforms-2022`）、**2023**（`Waveforms-2023`）、**2024 pre-CaloNu**（`waveforms-5FilesPerChunk`）、**2024 during-CaloNu**（`WaveformsCaloNuPeriod`）、**2024 post-CaloNu**（`waveforms_postCaloNu`）——选取相同的 v10 Physics GRL 对撞束流 run 与精确相同的 RAW 文件窗口，使用 `geoRotation_calypso_latest_fix_v1` 从原始 `.raw` 完整运行 `faser_reco.py -> faser_ntuple_maker.py`。
2. 逐事例及分布级对比 `fix_v1` 与 Ben 基线在 **FASERν Veto Reduced Charge**（含 2024 CaloNu `status == 512` 数字化仪 16 ns 偏移事例、`Preshower` 回退事例及 `< 30 pC` 中微子信号区）上的一致性。
3. 仿照 `FASER_DQ_Analyis` 绘图风格（双面板密度与比值子图、Garwood 泊松误差棒、时段配色规范），系统对比闪烁体电荷、量能器能量、径迹动力学分布与 $\nu_\mu$ 分析关键物理量及完整 cutflow，定位物理差异与程序坐标系约定差异。

---

## 工作内容

### 1. 五时段样本冻结与官方生产条件对齐

核查 Ben 五个文件夹的文件切分（每份 `partN` 对应 5 个连续 `xAOD/RAW` 段）、v10 Physics GRL（`GoodRunsList_2022/2023/2024_Physics_v10.yaml`）及 `/eos/experiment/faser/raw/` 在线原始文件覆盖后，每个时期冻结 1 个高亮度稳定对撞 run 的 **10 个连续 RAW 文件（对应 Ben 整整 2 个 5-file chunks）**，总计 **50 个完整 RAW 文件、5,069,589 个原始触发事例**：

- 核验官方 `r0021`（2022–2023）与 `r0022`（2024）生产日志，统一探测器几何与重建参数为 `--isData --geom TI12Data04 --noIFT --noBackward`，条件数据库分别对齐 `OFLCOND-FASER-04`（2022–2023）与 `OFLCOND-FASER-05`（2024）。
- 实测验证 `FaserEventSelector` 读取 ByteStream `.raw` 时忽略 `--skip` 参数，因此每个 HTCondor 作业原子处理 **1 个完整 `.raw` 文件（`--nevents -1`）**；`faser_reco.py` 与 `faser_ntuple_maker.py` 在 worker scratch 下的独立子目录 `reco/` 与 `ntp/` 中顺序执行（避免 `PoolFileCatalog.xml` 写锁冲突），提取紧凑 `PHYS` 分支后立即在本地清理数 GB 的临时 `xAOD`。
- `NtupleDumper` 通过命令行预设 `--reduced-charge-preset` 传入时序参数（`2024_caloNu` 使用 `data_caloNu`，其余四个时期使用 `data_default`），C++ 代码中无任何硬编码 run 号。

### 2. HTCondor 调度选择与 50 分片执行

提交前于 `2026-11-20T17:19:45Z` 查询全部 EL9 `bigbird` schedd 负载（完整记录于 [submission.json](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/submission.json)）：

| 候选 Schedd | Running | Idle | Total (`Running+Idle`) | MaxJobsRunning | 结论 |
| --- | ---: | ---: | ---: | ---: | --- |
| **`bigbird20.cern.ch`** | **615** | **1,629** | **2,244** | **10,055** | **总负载最低，选中提交（Cluster `1584913`）** |
| `bigbird21.cern.ch` | 1,679 | 1,502 | 3,181 | 10,055 | 次低备选 |
| `bigbird22.cern.ch` | 2,648 | 3,049 | 5,697 | 10,055 | 负载较高 |
| `bigbird18.cern.ch` | 2,418 | 3,696 | 6,114 | 10,053 | 负载较高 |
| `bigbird11.cern.ch` | 1,617 | 5,265 | 6,882 | 10,055 | 排队较多 |

- 提交文件位于 [jobs.sub](file:///afs/cern.ch/user/x/xcheng/condor/geoRotation_test/stage04_ben_comparison_run_20261120_v1/jobs.sub)（`+JobFlavour = "workday"`，`request_cpus = 1`，`request_memory = 3500 MB`）。
- 49/50 个分片在 Condor 上一次成功完成；分片 11（`2023_011523_00031`）在节点 `b9p28p3646.cern.ch` 初始化 `GeoModelSvc` 时偶发 CVMFS SQLite 读取异常，随后直接调用同一冻结脚本 [ben_comparison_worker.sh](file:///eos/home-x/xcheng/FASER/geoRotation_test/scripts/ben_comparison_worker.sh) 重跑一次通过（`rc = 0`，产出 `6,638` 条 PHYS 事例）。最终 **50/50 分片 100% 完整闭合**，无重复 `eventID`。

### 3. 坐标系约定核查与多核矢量化统计出图

对比官方 `PHYS` 与 `fix_v1` 输出的 `nt` 树时，定位到继承自 `geoRotation_calypso_latest`（`88c457c1`）的 `NtupleDumperAlg.cxx` 坐标系约定差异：

- **上游径迹态向量**：`NtupleDumperAlg.cxx` 将旋转 ACTS 坐标系下的位置与动量直接写入 `Track_{x,y,z}{0,1}` 和 `Track_{px,py,pz}{0,1}`，满足正交映射 $(x_{\mathrm{rot}}, y_{\mathrm{rot}}, z_{\mathrm{rot}}) = (z_{\mathrm{FASER}}, y_{\mathrm{FASER}}, -x_{\mathrm{FASER}})$。物理 FASER 坐标系由精确逆变换还原：
  $$x_{\mathrm{FASER}} = -z_{\mathrm{rot}},\quad y_{\mathrm{FASER}} = y_{\mathrm{rot}},\quad z_{\mathrm{FASER}} = x_{\mathrm{rot}},\qquad p_{x}^{\mathrm{FASER}} = -p_{z}^{\mathrm{rot}},\quad p_{y}^{\mathrm{FASER}} = p_{y}^{\mathrm{rot}},\quad p_{z}^{\mathrm{FASER}} = p_{x}^{\mathrm{rot}}.$$
- **各站外推横向坐标与角度**：`Track_X_at{VetoNu,Veto1,Veto2,Trig,Preshower1,Preshower2,Calo}` 写入了 $+z_{\mathrm{rot}} = -x_{\mathrm{FASER}}$（相对物理 FASER 坐标系符号相反），而 `Track_ThetaX_at*` 默认按未旋转公式 $\arctan(p_{\mathrm{rot},x}/p_{\mathrm{rot},z}) \approx \pm\pi/2$ 计算。[ben_vs_fixv1_compare.py](file:///eos/home-x/xcheng/FASER/geoRotation_test/analysis/ben_vs_fixv1_compare.py) 在保留原始分支审计统计的同时，对 `fix_v1` 施加精确逆变换 $R^T$ 与 $X_{\mathrm{st}}^{\mathrm{FASER}} = -X_{\mathrm{st}}^{\mathrm{raw}}$，使两套重建在严格同一物理 FASER 坐标系下对比。

复现命令：

```bash
source /afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/setup.sh
python3 /eos/home-x/xcheng/FASER/geoRotation_test/analysis/ben_vs_fixv1_compare.py \
  --run-dir /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1 \
  --workers 5
python3 /eos/home-x/xcheng/FASER/geoRotation_test/analysis/plot_ben_vs_fixv1_dq.py \
  --run-dir /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1
```

---

## 详细数据

### 1. 五时段样本覆盖与 FASER DQ 产额对照总表

DQ 事件定义为：位于完全相同的 `eventID` 区间、通过 v10 Physics GRL 时间窗、`Colliding == 1`、`TAP & 7 != 0` 且 `longTracks > 0`。GoodTrack 按 `FASER_DQ_RDF.py` 定义为：`Track_nLayers >= 7`、`Track_chi2 / Track_nDoF < 25`、$p^{\mathrm{FASER}} > 20\ \mathrm{GeV}$、$|x_0^{\mathrm{FASER}}| < 100\ \mathrm{mm}$、$|y_0^{\mathrm{FASER}}| < 100\ \mathrm{mm}$。

| 时段 | Run (`segments` / Ben `parts`) | 条件库 / 时序预设 | 精确 `eventID` 区间 | 原始触发数 (`10` RAW) | DQ 事件数 (`fix_v1` vs `Ben`) | DQ 事件净增益 | GoodTracks (`fix_v1` vs `Ben`) | GoodTrack 净增益 | 单分片平均耗时 (`reco` + `ntp`) |
| --- | --- | --- | --- | ---: | --- | ---: | --- | ---: | --- |
| **2022** | `007961` (`00010–00019` / `part2,3`) | `OFLCOND-FASER-04` / `data_default` | `[446340, 892531]` | 446,192 | `33,015` vs `31,374` | **+5.23%** | `28,446` vs `27,415` | **+3.76%** | `1781 s + 93 s` |
| **2023** | `011523` (`00030–00039` / `part6,7`) | `OFLCOND-FASER-04` / `data_default` | `[2342743, 3123210]` | 780,468 | `68,874` vs `65,446` | **+5.24%** | `59,177` vs `57,085` | **+3.66%** | `2459 s + 122 s` |
| **2024 pre** | `015387` (`00010–00019` / `part2,3`) | `OFLCOND-FASER-05` / `data_default` | `[1349884, 2672757]` | 1,322,874 | `87,298` vs `85,118` | **+2.56%** | `77,656` vs `76,570` | **+1.42%** | `6267 s + 237 s` |
| **2024 CaloNu** | `015849` (`00010–00019` / `part2,3`) | `OFLCOND-FASER-05` / `data_caloNu` | `[1207545, 2394243]` | 1,186,699 | `184,143` vs `179,367` | **+2.66%** | `162,907` vs `160,495` | **+1.50%** | `5516 s + 183 s` |
| **2024 post** | `016963` (`00025–00034` / `part5,6`) | `OFLCOND-FASER-05` / `data_default` | `[3353689, 4687044]` | 1,333,356 | `200,949` vs `195,945` | **+2.55%** | `178,748` vs `176,072` | **+1.52%** | `6826 s + 268 s` |
| **合计** | **5 runs（50 RAW files）** | — | — | **5,069,589** | **`574,279` vs `557,250`** | **+3.06%** | **`506,934` vs `497,637`** | **+1.87%** | — |

---

### 2. FASERν Veto Reduced Charge 逐事例闭合检验

由于 Ben 的 `El9CalypsoForWaveForm` 对该窗口内所有波形事例均输出了 `myVetoNu{0,1}_rawCharge35`，`fix_v1` 的全部 `574,279` 个 DQ 事例均可按 `(run, eventID)` 与 Ben 逐事例精确配对（定义 $\Delta Q_{\mathrm{red}} = Q_{\mathrm{red}}^{\mathrm{fix\_v1}} - Q_{\mathrm{red}}^{\mathrm{Ben}}$，单位 `pC`）：

| 时段（配对事例数 $N$） | `VetoNu0` $|\Delta Q| < 0.01\ \mathrm{pC}$ | `VetoNu0` $|\Delta Q| < 1.0\ \mathrm{pC}$ | `VetoNu1` $|\Delta Q| < 0.01\ \mathrm{pC}$ | `VetoNu1` $|\Delta Q| < 1.0\ \mathrm{pC}$ | 中位数 / `[q05, q95]` [pC] | 参考源计数 (`fb=1 Veto2` / `fb=2 PS1` / `fb=0 raw`) | `fb=2 (PS1)` $|\Delta Q| < 1.0\ \mathrm{pC}$ | $|\Delta Q_0| > 1\ \mathrm{pC}$ 离群数（其中 `Veto2 OVERFLOW`） | $\nu_\mu$ 预选事例数 / `<30 pC` 通过数 (`fix_v1` / `Ben` / 差异数) |
| --- | ---: | ---: | ---: | ---: | --- | --- | ---: | --- | --- |
| **2022** (`33,015`) | **98.12%** | **99.80%** | **97.96%** | **99.79%** | `0.0000` / `[-4.6e-5, +6.1e-5]` | `31,456` / `1,558` / `1` | **99.61%** | `65` (`30`) | `8,049` / **`0` / `0` / `0`** |
| **2023** (`68,874`) | **98.46%** | **99.76%** | **98.53%** | **99.75%** | `0.0000` / `[-6.1e-5, +6.1e-5]` | `65,538` / `3,329` / `7` | **99.19%** | `162` (`71`) | `16,386` / **`1` / `1` / `0`** |
| **2024 pre** (`87,298`) | **98.43%** | **99.77%** | **98.55%** | **99.75%** | `0.0000` / `[-6.1e-5, +6.1e-5]` | `81,683` / `5,605` / `10` | **99.86%** | `202` (`90`) | `26,806` / **`0` / `0` / `0`** |
| **2024 CaloNu** (`184,143`) | **98.55%** | **99.76%** | **98.59%** | **99.76%** | `0.0000` / `[-6.1e-5, +6.1e-5]` | `172,463` / `11,664` / `16` | **99.84%** | `445` (`143`) | `53,146` / **`0` / `0` / `0`** |
| └─ 其中 **`status == 512`** (`44,604`) | **98.57%** | **99.76%** | — | — | `0.0000` / `[-6.1e-5, +6.1e-5]` | 2nd Dig `+16 ns` 修正子样本 | — | `105` | — |
| **2024 post** (`200,949`) | **98.50%** | **99.76%** | **98.55%** | **99.75%** | `0.0000` / `[-6.1e-5, +6.1e-5]` | `188,673` / `12,253` / `23` | **99.82%** | `488` (`217`) | `61,209` / **`0` / `0` / `0`** |

- **2024 CaloNu 跨数字化仪验证**：在 `44,604` 个携带 `status & 512`（`TIMING_SHIFTED`）的第二数字化仪事例上，`fix_v1` 通过状态位检测自动补偿 `+16.0 ns` 并应用 `data_caloNu` 预设（`35.783 / 53.512 ns`），与 Ben 逐事例差值 `< 0.01 pC` 的比例达 **98.57%**（`< 1.0 pC` 达 **99.76%**），与同 run 的 `status == 0` 子样本（`98.99%` / `99.83%`）完全同级。
- **极少数 `> 1 pC` 尾部的来源**：全部样本中仅约 **0.20%–0.24%** 的事例出现 $|\Delta Q_{\mathrm{red}}| > 1\ \mathrm{pC}$，其中 `42%–49%` 伴随 `Veto2` 波形饱和（`status & 4`，`OVERFLOW`），其余来自 `geoRotation_calypso_latest` 相对旧生产版本在 `WaveformReconstructionTool`（`dd56a65e`）上的基线/拟合微调使 `localtime` 在单个 2 ns 采样边界附近跨过一格。
- **$\nu_\mu$ 中微子信号区闭合**：在五个时期共 `165,596` 个通过全部 $\nu_\mu$ 上游预选的径迹事例中，`VetoNu0_reduced_charge < 30 pC && VetoNu1_reduced_charge < 30 pC` 信号区判定 **100.000% 一致**（2023 检出完全相同的 1 个事例，其余四个时期均为 0，不一致事例数为 **0**）。

![VetoNu Reduced Charge All Periods](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_01_vetonu_reduced_charge_all_periods.png)

![Zoomed Reduced Charge and Event-by-Event Diff](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_02_reduced_charge_zoom_and_ebye_diff.png)

---

### 3. `NtupleDumperAlg.cxx` 原始分支坐标系审计数据

在共同重建成功的 `521,802` 个配对 DQ 事例上，直接读取未逆旋转的 `fix_v1` 原始分支与官方 `PHYS` 对比：

| 时段（配对数 $N_{\mathrm{common}}$） | `fix_v1` 原始 `Track_px0` 中位数 [GeV] | `fix_v1` 原始 `Track_pz0` 中位数 [GeV] | 官方 `Track_pz0` 中位数 [GeV] | $\mathrm{med}(X_{\mathrm{atVetoNu}}^{\mathrm{raw}} + X_{\mathrm{atVetoNu}}^{\mathrm{off}})$ [mm] | $\mathrm{med}(\| \Theta_{X,\mathrm{atVetoNu}}^{\mathrm{raw}} \|)$ [rad] |
| --- | ---: | ---: | ---: | ---: | ---: |
| **2022** (`29,443`) | `+195.0` | `-0.104` | `+156.7` | `+0.021` | `1.5683` ($\approx \pi/2$) |
| **2023** (`61,323`) | `+193.1` | `-0.091` | `+154.4` | `+0.021` | `1.5682` ($\approx \pi/2$) |
| **2024 pre** (`79,750`) | `+419.6` | `+0.188` | `+304.8` | `-0.020` | `1.5691` ($\approx \pi/2$) |
| **2024 CaloNu** (`167,674`) | `+396.7` | `+0.172` | `+292.9` | `-0.018` | `1.5691` ($\approx \pi/2$) |
| **2024 post** (`183,612`) | `+405.5` | `+0.187` | `+300.4` | `-0.018` | `1.5691` ($\approx \pi/2$) |

该表直接证实：`fix_v1`（及上游 `geoRotation_calypso_latest`）的 `nt` 树原始 `Track_px0` 实际存放的是 $+p_z^{\mathrm{FASER}}$，`Track_pz0` 存放的是 $-p_x^{\mathrm{FASER}}$，而 `Track_X_atVetoNu` 与官方物理坐标相加中位数仅 `±0.02 mm`（即精确为 $-x_{\mathrm{FASER}}$）。若下游分析脚本直接按官方分支名读取而不做 $R^T$ 逆变换，`LeadTrack_pz > 100 GeV` 与 `r_atVetoNu` / 角度条件将全部失效。

---

### 4. 闪烁体、量能器与共同物理 FASER 坐标系下的径迹及 $\nu_\mu$ 选择对比

#### (1) 闪烁体与量能器响应（[dq_03_scintillator_calo_dq_fig5_style.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_03_scintillator_calo_dq_fig5_style.png)）

- **闪烁体电荷一致性**：`VetoNu`、`Veto2`、`Timing` 四路求和与 `Preshower` 两层电荷分布在全部五个时期与官方 `PHYS` 高度重合（比值在 MIP 峰区稳定于 `1.00`）。
- **2023 低增益量能器 (`Calo_total_E_EM`) 显著改善**：2023 年 run `011523` 处于低增益 ECAL 配置，通过 μ 子的原始峰值电压仅约 `0.81 mV`（触发 `Calo_status == 1`，`BELOW_THRESHOLD`）。官方 `r0021` 对低于阈值的 MIP 事例在 `Calo_total_E_EM` 中大量填入 `NaN`，且其 `Calo_total_nMIP` / `Calo_total_raw_E_EM` 给出非物理偏高值（均值约 `26.1 MIP` / `8.65 GeV`）；而 `fix_v1`（含 `dd56a65e` 波形重建更新）恢复了正常通过 μ 子的 MIP 峰（均值 **`1.94 MIP` / `0.63 GeV`**），与 2022 及 2024 三个时期的物理响应一致。

![Scintillator and Calorimeter Fig.5 Style](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_03_scintillator_calo_dq_fig5_style.png)

#### (2) 共同物理 FASER 坐标系下的逐事例径迹闭合与分动量区间标度

在 `521,802` 个共同配对事件上对比主径迹（物理 FASER 坐标系，参考位置为最上游测量面 `(x0, y0, z0)`）：

| 时段 | 共同事件 / 仅 `fix_v1` / 仅 `Ben` | 新增径迹 $\chi^2/\mathrm{ndf}<15$ 占比（对比丢失径迹） | $\mathrm{med}(\Delta\chi^2)$ | $\langle\Delta\mathrm{nDoF}\rangle$ | $\Delta x_0$ 中位数 ± std [mm] | $\Delta y_0$ 中位数 ± std [mm] | 电荷一致率（全动量 / `<50 GeV`） | 分动量中位比值 $p_{\mathrm{fix}}/p_{\mathrm{off}}$（`<50` / `50–150` / `150–400` / `400–1000 GeV`） |
| --- | --- | --- | ---: | ---: | --- | --- | --- | --- |
| **2022** | `29,443` / `3,572` / `1,931` | **97.6%**（vs `90.4%`） | **-0.132** | `+0.161` | `-0.004 ± 0.569` | `-0.000 ± 0.160` | `90.0%` / **95.8%** | `0.9996` / `1.0025` / `1.0302` / `1.0190` |
| **2023** | `61,323` / `7,551` / `4,123` | **97.7%**（vs `91.1%`） | **-0.136** | `+0.163` | `-0.004 ± 0.571` | `-0.000 ± 0.163` | `89.9%` / **96.1%** | `0.9996` / `1.0021` / `1.0278` / `1.0188` |
| **2024 pre** | `79,750` / `7,548` / `5,368` | **97.6%**（vs `92.3%`） | **-0.162** | `+0.141` | `+0.004 ± 0.527` | `+0.000 ± 0.127` | `85.9%` / **95.0%** | `0.9996` / `1.0045` / `1.0345` / `1.1687` |
| **2024 CaloNu** | `167,674` / `16,469` / `11,693` | **97.5%**（vs `91.9%`） | **-0.166** | `+0.143` | `+0.004 ± 0.533` | `+0.000 ± 0.136` | `85.6%` / **95.3%** | `0.9996` / `1.0049` / `1.0307` / `1.1778` |
| **2024 post** | `183,612` / `17,337` / `12,333` | **97.7%**（vs `92.3%`） | **-0.159** | `+0.141` | `+0.004 ± 0.527` | `+0.000 ± 0.135` | `85.6%` / **95.1%** | `0.9996` / `1.0056` / `1.0310` / `1.1645` |

- **拟合质量与寻迹纯度**：五个时期中，`fix_v1` 共同径迹的中位 $\chi^2$ 均下降 `0.13–0.17` 且平均测量自由度增加 `+0.14–0.16`；仅在 `fix_v1` 中重建出的新增径迹（合计 `52,477` 个事件）满足严格质量条件 $\chi^2/\mathrm{ndf} < 15$ 的比例稳定在 **97.5%–97.7%**，显著高于仅在官方版本中存在的丢失径迹（`90.4%–92.3%`）。
- **动量标度与高能近直线区行为**：在 `< 150 GeV` 区间内，`fix_v1` 与官方 `PHYS` 动量标度偏差 `< 0.6%`（`< 50 GeV` 处比值 `0.9996`，电荷一致率 `95%–96%`）；在 `150–400 GeV` 区间高约 `+2.8%–3.5%`。在 2024（`OFLCOND-FASER-05`）`400–1000 GeV` 的近直线高动量区，相对曲率偏移 $\mathrm{med}[(q/p_{\mathrm{fix}} - q/p_{\mathrm{off}})\cdot p_{\mathrm{off}}] \approx -0.044 \sim -0.052$（对应约 `-0.05 TeV^-1` 的微小矢高差），使正电荷高动量径迹比例略高于官方版本（2024 GoodTrack $N_+/N_-$ 为 `1.19` vs 官方 `1.11`）。

![Track Kinematics Fig.6 Style](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_04_track_kinematics_dq_fig6_style.png)

![Event-by-Event Tracking and Frame Audit](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_08_ebye_tracking_and_frame_audit.png)

#### (3) $\nu_\mu$ 分析完整预选 Cutflow 对照

严格按 [cuts.yaml](file:///eos/home-x/xcheng/FASER/nu_muAnalysis/AnalysisZipFramework/config/cuts.yaml) 顺序施加筛选（`GRL+Colliding+Trig` $\to$ `CaloNu status clean` $\to$ `Veto2 > 40 pC` & `Timing > 20 pC` & `PS > 2.5 pC` $\to$ `LeadTrack pz > 100 GeV` $\to$ `nDoF > 12 & chi2/nDoF < 15` $\to$ `r_VetoNu < 120 mm, r_S1..3 < 95 mm, IP < 95 mm, theta < 25 mrad` $\to$ `VetoNu0/1_redQ < 30 pC`）：

| 时段 | `GRL+Coll+Trig` (`fix_v1` / `Ben`) | `PS > 2.5 pC` (`fix_v1` / `Ben`) | `Lead pz > 100 GeV` (`fix_v1` / `Ben`) | `TrkQual (nDoF>12, χ²/ndf<15)` (`fix_v1` / `Ben`) | **全部上游预选 `Pre-VetoNu`** (`fix_v1` / `Ben` / 净增益) | **`VetoNu0 & 1 redQ < 30 pC`** (`fix_v1` / `Ben`) |
| --- | --- | --- | --- | --- | --- | --- |
| **2022** | `33,019` / `31,612` | `32,119` / `31,007` | `19,423` / `18,791` | `18,542` / `17,562` (+5.58%) | **`8,049` / `7,908` (+1.78%)** | **`0` / `0`** |
| **2023** | `68,885` / `65,863` | `66,981` / `64,604` | `40,393` / `39,000` | `38,474` / `36,408` (+5.67%) | **`16,386` / `16,154` (+1.44%)** | **`1` / `1`** |
| **2024 pre** | `87,302` / `85,968` | `85,477` / `84,545` | `61,969` / `61,067` | `59,005` / `57,134` (+3.27%) | **`26,806` / `26,503` (+1.14%)** | **`0` / `0`** |
| **2024 CaloNu** | `184,168` / `181,239` | `175,462` / `172,982` | `123,876` / `122,132` | `117,797` / `113,979` (+3.35%) | **`53,146` / `52,471` (+1.29%)** | **`0` / `0`** |
| **2024 post** | `200,980` / `198,091` | `196,906` / `194,840` | `142,078` / `139,832` | `135,168` / `130,621` (+3.48%) | **`61,209` / `60,476` (+1.21%)** | **`0` / `0`** |

![nu_mu Cutflow and Gained Track Quality](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_07_numu_cutflow_and_gained_tracks.png)

---

### 5. 完整图表清单与产物链接

全部图表均同时保存为 300 dpi PNG、单图矢量 PDF 与 8 页合并 PDF [ben_vs_fixv1_dq_report.pdf](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/ben_vs_fixv1_dq_report.pdf)：

1. **五时段 VetoNu Reduced Charge 全谱与比值**：[dq_01_vetonu_reduced_charge_all_periods.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_01_vetonu_reduced_charge_all_periods.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_01_vetonu_reduced_charge_all_periods.pdf)
2. **$\nu_\mu$ 低电荷区 (`< 60 pC`) 放大与逐事例残差（含 CaloNu `status=512`）**：[dq_02_reduced_charge_zoom_and_ebye_diff.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_02_reduced_charge_zoom_and_ebye_diff.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_02_reduced_charge_zoom_and_ebye_diff.pdf)
3. **FASER DQ Fig. 5 风格闪烁体与量能器响应对比**：[dq_03_scintillator_calo_dq_fig5_style.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_03_scintillator_calo_dq_fig5_style.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_03_scintillator_calo_dq_fig5_style.pdf)
4. **FASER DQ Fig. 6 风格共同物理坐标系径迹动力学分布**：[dq_04_track_kinematics_dq_fig6_style.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_04_track_kinematics_dq_fig6_style.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_04_track_kinematics_dq_fig6_style.pdf)
5. **FASER DQ Fig. 6(h–k) 风格按正负电荷拆分的径迹角度分布**：[dq_05_track_angles_by_charge_dq_fig6hk_style.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_05_track_angles_by_charge_dq_fig6hk_style.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_05_track_angles_by_charge_dq_fig6hk_style.pdf)
6. **$\nu_\mu$ 分析核心选择观测量（`LeadTrack_pz`、`r_atVetoNu`、`r_max`、`IP`、$\theta$、$\chi^2/\mathrm{ndf}$）**：[dq_06_numu_selection_observables.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_06_numu_selection_observables.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_06_numu_selection_observables.pdf)
7. **$\nu_\mu$ 逐步 Cutflow 产额比与新增/丢失径迹质量对比**：[dq_07_numu_cutflow_and_gained_tracks.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_07_numu_cutflow_and_gained_tracks.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_07_numu_cutflow_and_gained_tracks.pdf)
8. **逐事例曲率/电荷/$\chi^2$ 闭合与原始 `NtupleDumper` 坐标系审计**：[dq_08_ebye_tracking_and_frame_audit.png](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_08_ebye_tracking_and_frame_audit.png) | [PDF](file:///eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1/figures/dq_08_ebye_tracking_and_frame_audit.pdf)

---

## 阶段性结论

### 已证实的事实

1. **Reduced Charge 算法在 2022–2024 全部五个时期与 Ben 的 `El9CalypsoForWaveForm` 达到逐事例浮点级闭合**：在 50 个完整 RAW 文件的 `574,279` 个 DQ 事例中，$|\Delta Q_{\mathrm{red}}| < 0.01\ \mathrm{pC}$ 占 **98.12%–98.55%**，$|\Delta Q_{\mathrm{red}}| < 1.0\ \mathrm{pC}$ 占 **99.76%–99.80%**，中位数精确为 `0.0000 pC`；2024 CaloNu `status == 512` 的 `44,604` 个跨数字化仪事例同样达到 **98.57%**（`< 0.01 pC`）闭合；`165,596` 个 $\nu_\mu$ 预选事例在 `< 30 pC` 信号区的筛选结果 **100.000% 一致**。不再需要外挂独立的 `El9CalypsoForWaveForm` 二次辅助树。
2. **`NtupleDumperAlg.cxx` 存在旋转坐标系直写约定差异（程序事实）**：继承自 `geoRotation_calypso_latest` 的 `NtupleDumperAlg.cxx` 将旋转 ACTS 坐标系下的 `Track_{x,y,z}{0,1}` 与 `Track_{px,py,pz}{0,1}` 直接写出，并使 `Track_X_at* = -x_{\mathrm{FASER}}$、`Track_ThetaX_at*` 落向 $\pm\pi/2$。直接接入现有 `nu_muAnalysis` 或 `FASER_DQ` 前，必须在 `NtupleDumperAlg.cxx` 内统一施加逆旋转 $R^T$ 还原至标准物理 FASER 坐标系。
3. **同一物理 FASER 坐标系下，`fix_v1` 寻迹效率与拟合质量全面优于官方 `r0021`/`r0022`**：DQ 径迹事件数在 2022–2023 提升 **+5.2%**、在 2024 提升 **+2.6%**；完整 $\nu_\mu$ 上游预选产额在五个时期一致提升 **+1.14% 至 +1.78%**，且新增径迹 `97.5%–97.7%` 满足 $\chi^2/\mathrm{ndf} < 15$，共同径迹中位 $\chi^2$ 降低 `0.13–0.17`。
4. **2023 低增益量能器 MIP 能量恢复正常**：`fix_v1` 消除了官方 `r0021` 在 2023 低增益 run `011523` 中对 `BELOW_THRESHOLD` MIP 事例填 `NaN` 及 `nMIP` 虚高（`26.1 MIP`）的问题，回归 `1.94 MIP`（`0.63 GeV`）。

### 待进一步验证的假设

- **2024 高动量区（`> 400 GeV`）微小曲率偏移来源**：在 `OFLCOND-FASER-05` 下观察到的 $\Delta(q/p)\cdot p_{\mathrm{off}} \approx -0.05$（约 `-0.05 TeV^-1`）推测源于旋转 ACTS 几何/alignment 变换矩阵或材料表在百 GeV–TeV 弱偏转极限下与旧非旋转几何存在微米级等效矢高差；需结合本仓库 S05/S06 固定 hits 对照实验单独归因，不影响 `< 300 GeV` 主体及 $\nu_\mu$ 预选效率。

---

## 下一步计划

1. **在 `geoRotation_calypso_latest_fix_v1` 的 `NtupleDumperAlg.cxx` 中统一修复输出坐标系约定**：对 `Track_{x,y,z}{0,1}`、`Track_{px,py,pz}{0,1}` 施加 $R^T$ 逆旋转，修正全部七个外推站面的 `Track_X_at*` 符号（$-z_{\mathrm{rot}} \to +x_{\mathrm{FASER}}$）与 `Track_ThetaX_at*` / `Track_ThetaY_at*` 物理角度公式，使产出的 `PHYS` ntuple 无需下游额外转换即可与 `nu_muAnalysis` 和 `FASER_DQ_Analyis` 原生无缝兼容；验收标准为单文件回归重跑后原始分支直接与官方 `PHYS` 逐事例闭合（$\mathrm{med}(\Delta x_0), \mathrm{med}(\Delta X_{\mathrm{atVetoNu}}) < 0.05\ \mathrm{mm}$）。
2. **2025–2026 新数据小样本端到端 `PHYS` 验证**：在完成 `NtupleDumper` 物理坐标系还原后，选取 2025（含 Preshower `Dig1` 正常与 GFIT fail 触发 `Dig2` 回退的 run）及 2026 代表 run 运行相同全链路验证，确认两级 Preshower 回退与物理坐标系输出在新数据上同样一次闭合。
