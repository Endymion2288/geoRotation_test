# S06：共同有限先验下的原生 CKF、conversion 与 refit

**四组有限先验控制中，beforeGeo/afterGeo 的 CKF 候选存在性与首选候选的 charge sign 全部一致；conversion 保留 q/p 和 sign。n32、分量内 C/256 的 refit 后出现四个版本间符号差异，其中两例同 hits 且双方 refit 成功，两例同时伴随路径差别和 afterGeo fallback。其余三组没有最终符号差异。CKF q/p 数值、部分 hit 路径及 refit 数值仍有差异，完整共同物理先验的全链路归因继续执行。**

后续[实际 refit 输入与对角 inflation 单因素检查](2026-11-14_stage06_refit_input_covariance.md)已完成15个重建与精确 gate：只取消额外对角放大，两条 same-hit 符号分歧消失，一条两版都错、一条都对；两条不同路径/after fallback 的分歧保留。下文“待检验”和下一步保留为本阶段提出的问题，实测结论以上述后续报告为准。

## 阶段目标

推进 [Q5 纯坐标旋转检查](../roadmap/04_charge_misidentification.md)：在真实 Calypso 中将同一个物理 seed 分布分别表达于 FASER / rotated ACTS 坐标，沿 actual finder、measurement update、候选独立 smoother、ACTS→Trk conversion 和 native refit 比较结果。将进入 CKF 的共同先验、后续 refit 重新构造的先验以及 fallback 分开检查。

isolated smoothing 始终开启。这里的 component 是相关的诊断分量；物理事件数为三个 DSID 各20个，共60个源事件。component 比例只描述本项控制，生产误判率继续使用 μ± 各500k 的[全量基线](2026-11-13_stage06_charge_baseline.md)。

## 工作内容与输入

结果根目录为 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`，本轮全链路目录为 `P=C/common_prior_component_refit_v1`，CKF-only 输入与对照为 `C/common_prior_components_v4`。

使用各 DSID 首个 RDO 的 entries `[0,20)`，保留原生空 seed 事件：

| DSID | 输入目录与文件 | 有 seed 的源事件 / 总源事件 |
| --- | --- | ---: |
| 100044 μ+ | `/eos/experiment/faser/data0/sim/mc24/particle_gun/100044/rdo/s0013/FaserMC-MC24_PG_mupl_fasernu_5mrad_flukaE-100044-00000-00099-s0013-RDO.root` | 20/20 |
| 100116 μ+ | `/eos/experiment/faser/data0/sim/mc24/particle_gun/100116/rdo/s0012/FaserMC-MC24_PG_mupl_fasernu_2D_fluka-100116-00000-00009-s0012-RDO.root` | 17/20 |
| 100117 μ− | `/eos/experiment/faser/data0/sim/mc24/particle_gun/100117/rdo/s0012/FaserMC-MC24_PG_mumi_fasernu_2D_fluka-100117-00000-00009-s0012-RDO.root` | 18/20 |

材料采用已由 native decoder 验证的 corrected AFS map 及其 FASER 表达；路径逐任务冻结。几何、磁场、measurement 和 transport 仍由本版 Calypso 执行。

### 有限先验的构造与适用范围

原始 beforeGeo Gaussian 的完整均值和6×6 covariance 为 `m,C`。以反号配对、白化后的有限点云分配分量中心，使中心总体均值为 `m`、中心间 covariance 为 `(1−f)C`；每个分量内部 covariance 为 `fC`。因此总均值/C 保持原始值，没有同时移动中心并保留完整 C 的重复展宽。

每个分量均值按物理位置和方向精确旋转；分量内 C 在其自身均值处按 Jacobian 转换。有限 Gaussian 的非线性坐标映射仍是一阶近似，另外记录每分量12个内部 sigma points 的方向差。该指标是已检查点上的误差，不是整个 Gaussian 支持域的上界。

| 设置 | 每 seed 分量数 | 分量内 covariance | 总 C 的最大尺度化误差 | 内部检查点最大方向 norm 差 |
| --- | ---: | --- | ---: | ---: |
| `matched_n32_f256_v1` | 32 | C/256 | 2.403×10⁻¹⁴ | 4.67788×10⁻⁵ |
| `matched_n32_f1024_v1` | 32 | C/1024 | 2.672×10⁻¹⁴ | 1.16981×10⁻⁵ |
| `matched_n64_f1024_v1` | 64 | C/1024 | 1.403×10⁻¹⁴ | 1.18668×10⁻⁵ |
| `matched_n64_f4096_v1` | 64 | C/4096 | 1.019×10⁻¹⁴ | 2.96694×10⁻⁶ |

相同分量数的 width 对照使用同一标准化点云，中心按保持总 C 所需的因子重缩放。n64 的 C/1024→C/4096 中心因子为1.0003665017531884，均值关系最大误差4.26×10⁻¹⁴，分量 C 比例关系精确通过。

各分量独立调用实际 Calypso finder 和本候选 private smoother，不在不同分量之间进行 shared-hit winner 竞争，也没有用诊断输出实现 likelihood 加权的 Gaussian-mixture tracking。选中轨迹按 measurements 降序、χ²升序和 candidate ID 取第一条用于配对；本项没有额外 Good95 gate。失败、空 seed、无候选与 refit fallback 均保留。

### 实际源码、构建与开关

| 用途 | 实际 Calypso worktree |
| --- | --- |
| beforeGeo CKF component hook | `/eos/home-x/xcheng/FASER/geoRotation_calypso_components_beforeGeo` |
| afterGeo CKF component hook | `/eos/home-x/xcheng/FASER/geoRotation_calypso_components_afterGeo` |
| beforeGeo native conversion/refit | `/eos/home-x/xcheng/FASER/geoRotation_calypso_component_refit_beforeGeo` |
| afterGeo native conversion/refit | `/eos/home-x/xcheng/FASER/geoRotation_calypso_component_refit_afterGeo` |

beforeGeo commit 为 `1eb885ae22299fb39bb352bead6005d41b0f953e`，afterGeo 为 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`，均使用 Athena24.0.41、FASERNU-04 / OFLCOND-FASER-05。component hook 增量为 [0035 before](../patches/0035_native_ckf_components_beforeGeo.patch) / [0035 after](../patches/0035_native_ckf_components_afterGeo.patch)，全链路增量为 [0036 before](../patches/0036_component_native_refit_beforeGeo.patch) / [0036 after](../patches/0036_component_native_refit_afterGeo.patch)。这些 opt-in 开关默认关闭。

`ComponentRefit` 使用本版 `CreateTrkTrackTool` 和 `KalmanFitterTool1`；beforeGeo 另设联合 regular update/smoothing chart。`P/overlay_beforeGeo` / `overlay_afterGeo` 包含匹配 observer、geometry、CKF 及 Python metadata，逐运行核对实际加载的四个 ACTS/observer 库。

`P/*_complete_source.patch` 同时包含 tracked diff 和新增头文件，`source_additions_*` 保存新增文件。`P/source_validation.json` 核对 tracked diff、每个新增文件及完整 patch 的 reverse-check，补足仅保存 `git diff` 会漏掉未跟踪头文件的来源缺项。

## 详细数据

四组共24个 CKF-only jobs、36个 full-chain jobs，60个作业均重建成功。10,560个物理输入分量在两种坐标中配对，四类比较共42,240条记录。以下计数使用 signal、status=1 的 primary μ truth。

**truth 归约更正：**早期 exploratory component 分析取同一 event 的最后一条 truth，可能被 secondary 覆盖。当前只允许唯一 primary μ，权威统计为 `P/analysis_v3`，正式六页图从该目录生成。早期 component correct/wrong 计数已在结果目录用 `STATISTICS_SUPERSEDED.md` 标明替代关系。此更正不改变几何间 sign migration 数；生产基线使用独立 primary matcher。

### CKF 与最终 sign

表中 selected 数在两版完全相同；CKF wrong 数亦逐 component 相同。final wrong 按 before native / before regular / after native 排列。

| 设置 | DSID | 输入分量 | selected | CKF wrong | final wrong：before / regular / after | same selected hits |
| --- | --- | ---: | ---: | ---: | --- | ---: |
| n32 C/256 | 100044 | 640 | 311 | 89 | 87 / 84 / 86 | 306 |
| n32 C/256 | 100116 | 544 | 254 | 5 | 5 / 5 / 5 | 253 |
| n32 C/256 | 100117 | 576 | 360 | 45 | 45 / 46 / 46 | 327 |
| n32 C/1024 | 100044 | 640 | 209 | 59 | 59 / 59 / 59 | 204 |
| n32 C/1024 | 100116 | 544 | 186 | 1 | 1 / 1 / 1 | 186 |
| n32 C/1024 | 100117 | 576 | 292 | 41 | 41 / 41 / 41 | 256 |
| n64 C/1024 | 100044 | 1,280 | 425 | 134 | 134 / 134 / 134 | 417 |
| n64 C/1024 | 100116 | 1,088 | 369 | 12 | 12 / 12 / 12 | 369 |
| n64 C/1024 | 100117 | 1,152 | 584 | 81 | 81 / 81 / 81 | 536 |
| n64 C/4096 | 100044 | 1,280 | 303 | 96 | 96 / 96 / 96 | 293 |
| n64 C/4096 | 100116 | 1,088 | 262 | 7 | 7 / 7 / 7 | 262 |
| n64 C/4096 | 100117 | 1,152 | 501 | 70 | 70 / 70 / 70 | 449 |

![Native full-chain component outcomes](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_component_refit_v1/figures_v3/results/component_full_chain_outcomes.png)

conversion 对实际每条 selected candidate 保留 q/p 与 sign；增加 conversion/refit 后，原 CKF seed、selected 及全部 selected-state 记录与 CKF-only 对照精确一致。

### 实际 refit 成功与 fallback

| 设置 | DSID | selected | refit 成功 before / after | 双方成功 |
| --- | --- | ---: | --- | ---: |
| n32 C/256 | 100044 | 311 | 38 / 37 | 37 |
| n32 C/256 | 100116 | 254 | 24 / 24 | 24 |
| n32 C/256 | 100117 | 360 | 70 / 69 | 67 |
| n32 C/1024 | 100044 | 209 | 11 / 10 | 10 |
| n32 C/1024 | 100116 | 186 | 8 / 8 | 8 |
| n32 C/1024 | 100117 | 292 | 38 / 38 | 36 |
| n64 C/1024 | 100044 | 425 | 30 / 29 | 29 |
| n64 C/1024 | 100116 | 369 | 14 / 14 | 14 |
| n64 C/1024 | 100117 | 584 | 76 / 74 | 73 |
| n64 C/4096 | 100044 | 303 | 14 / 13 | 13 |
| n64 C/4096 | 100116 | 262 | 5 / 5 | 5 |
| n64 C/4096 | 100117 | 501 | 48 / 47 | 45 |

总分母4,056个首选 component 中，native before/after 成功376/368个，双方成功361个。大量低 measurement / outlier 候选回退到 converted CKF；因此 final q/p 在大部分分量上接近，不能用来证明 refit 的坐标闭合。CDF 另列双方 refit 成功的子组。

n32 C/256 中，100044 CKF 的 |Δq/p| 中位数为2.19443×10⁻⁹ GeV⁻¹、p95=1.42578×10⁻⁵ GeV⁻¹；final native 的 p95 增为5.01439×10⁻⁴ GeV⁻¹。在双方 refit 成功的37个 regular-before/after 对中，|Δq/p| 中位数9.46624×10⁻⁵ GeV⁻¹、p95=1.08678×10⁻³ GeV⁻¹，数值差尚未消失。

### 四个最终符号差异的定位

均位于 n32 C/256；source entry 以本次首个文件中的 entry 定义：

| DSID | entry / component | same CKF hits | refit 成功 before / after | final q/p before / after，GeV⁻¹ |
| --- | --- | --- | --- | --- |
| 100044 | 1 / 23 | 是 | 是 / 是 | −0.000652921 / +0.000460927 |
| 100044 | 3 / 19 | 是 | 是 / 是 | −0.000339920 / +0.001104541 |
| 100044 | 12 / 2 | 否 | 是 / 否 | +0.001457856 / −0.001581931 |
| 100117 | 0 / 25 | 否 | 是 / 否 | −0.002178604 / +0.000272720 |

before native→before regular 也有4个 sign migrations；regular before→after 剩2个，均在100044。单独更换 refit chart 尚未闭合整个链路，需要实际 refit 输入 C、mean 和 surface 的配对证据。

### width 与分量数敏感性

n64 的分量内 C/1024→C/4096 时，三个 DSID 的 selected 数分别425→303、369→262、584→501。虽然内部方向映射误差随 width 收缩，candidate 接受与 refit 成功集合仍变化。当前没有已收敛的 likelihood-weighted posterior 积分，不能由上述 component wrong 比例估计生产误判率。

![Sampling and covariance split sensitivity](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_component_refit_v1/figures_v3/results/component_sampling_sensitivity.png)

## 执行、验证与图表

全部批量重建、统计和导图经 HTCondor。提交前查询实际 schedd 广告/队列，选择兼容且 running+idle 最低的 bigbird20，排除已记录 EOS 兼容问题的 bigbird21；选择依据、UTC、资源和 cluster ID 均保存在 worker 的 `submission.json`。

full-chain build cluster1581748 两版均exit0。完整重建包括1581750/51、1581753/54、1581759/60，补充输入归约为1581752/58；精确任务清单在 `P/*tasks.json` 和各 job 的 `source_task.json`，CKF-only 对照清单在 `common_prior_components_v4`。重建申请1 CPU，构建申请4 CPU、12000 MB memory，线程数按 worker 配置固定。60个实际 Athena 进程的 CPU 合计987.30 s、最大 RSS1,041,876 KiB，逐运行 GNU time 数值保存在 `P/reconstruction_resources.json`；CPU 合计不当作并行 wall time。

canonical 统计 cluster1581761 的统计与绘图子步骤成功，初次 PDF 校验因局部 font alias 未通过；该序列退出1，未据此重跑 tracking。其 GNU time 记录 wall339.08 s、CPU311.13 s、最大 RSS1,139,196 KiB。新版六图导出1581775 成功，校验子步骤因冻结目录漏拷贝 validator helper 失败；补齐 helper 后单独校验1581780 exit0，无物理数据改动。

`analysis_v3/validation.json` 对60个作业核对完整 seed mean/C/surface、空 seed、terminal records、实际加载库、chart/refit 开关及 native update 的 q/p gain 恒等式。36组 full-chain 的 seed/CKF 状态与原记录精确相同，conversion sign 无变化。绘图从原始 paired rows 独立重算 stage counts 与 q/p max/median/p95；6页 PDF 与逐图 PDF 的归一化 vector operations、font glyph outlines/encoding/widths、page boxes 和 rotation 一致。六张 PNG 已逐张原始分辨率目视检查，坐标、图例和文字可读。

复现入口示例，先冻结配置和匹配库再提交：

```bash
python3 scripts/prepare_charge_component_refit.py \
  --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1 \
  --phase tasks --series matched_n32_f256_v1 matched_n64_f1024_v1
python3 analysis/charge_component_refit.py \
  --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1 \
  --full-chain /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_component_refit_v1 \
  --series matched_n32_f256_v1 matched_n32_f1024_v1 matched_n64_f1024_v1 matched_n64_f4096_v1 \
  --output <new-results-directory>
```

实际批量命令使用 `P/analysis_sequence_tasks_v3.json`，经 `scripts/submit_momentum.py --stage S06` 执行。已有输出目录不可覆盖；代码、配置和逐命令日志与结果一同冻结。

- [6页 component 图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_component_refit_v1/figures_v3/results/component_refit_report.pdf)：charge、q/p配对、CDF、逐stage sign、fallback及敏感性。
- [权威统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_component_refit_v1/analysis_v3/summary.json)、[逐component配对](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_component_refit_v1/analysis_v3/paired_components.json)、[数值验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_component_refit_v1/analysis_v3/validation.json)。
- [图表验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_component_refit_v1/figures_v3/results/validation.json)、[目视记录](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_component_refit_v1/figures_v3/results/visual_review.json)。
- [115页 S06 汇报 PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v11/charge_stage06_report.pdf)保留原109页，新增本轮六页；[页面衔接核对](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v11/lineage_validation.json)验证前109页的index、vector内容及page boxes精确不变，继承其已有目视记录。组装cluster1582082 exit0。

## 阶段性结论

事实：共同有限输入的四组 CKF 首选存在性/sign 一致，conversion 保留符号；数值差及部分路径差仍在，n32 C/256 的最终 sign 差异集中到 refit/fallback。后验统计对分量内 width 有可见敏感性。新的 truth gate 更正了 exploratory component 计数，权威图与统计已统一到 primary μ。

待检验原因：两版 refit 都在 `getParametersFromTrack` 中单独放大 C 的对角项，随后整体乘10；对角放大在不同 angular chart 中通常不表示同一物理先验。native conversion/refit 也重建参考面和局部 covariance basis。源码提示这些候选原因，尚需实际输入与单因素结果区分。

## 下一步计划

1. beforeGeo 增加只读 `refit_input` / `refit_state` 审计，首先验证原 seed、CKF、conversion、final 全字段不变；与 afterGeo 的原审计比较物理 mean/C/surface。
2. 在同一 n32 C/256 输入中，将两版 `SeedCovarianceScale` 从10变为1，保留完整 C×10，逐项核对只有对角 inflation 改变；单列同 hits、双方 refit 成功和 fallback，检验四个 sign 差异是否保留。
3. 如仍有差异，对这两条 same-hit、双方成功的路径追踪 refit 第一次 predicted/update/smoothing 分叉，必要时固定共同 refit prior。完成后更新 Q5 审计；S06 保持 active。
