# S06：实际 refit 输入与对角 covariance 放大的单因素检查

后续[实际 CKF chart 四模式控制](2026-11-14_stage06_ckf_direction_chart.md)已定位本批39条路径差异，并联合降低物理state差异；宽度/分量数稳健性及极端状态仍在继续检查。

**两版 refit 都先将 covariance 对角项额外乘10，再将完整矩阵乘10。只取消前一步，两例“相同 CKF hits、双方 refit 成功”的版本间最终符号差异均消失；另两例不同路径、afterGeo fallback 的差异保留。第一例变成两版都错，第二例两版都对，因此该控制证明 refit 对物理先验和 chart 的敏感性，不代表生产误判率改善。**

## 阶段目标

继续 [Q5](../roadmap/04_charge_misidentification.md)，检验[共同有限先验全链路报告](2026-11-14_stage06_common_prior_full_chain.md)中四个最终 sign 差异是否来自实际 refit 重新构造的 mean、covariance、surface 和对角 inflation。isolated smoothing 全程开启；tracking、conversion、refit、transport 均由实际 Calypso 执行。

结果根目录 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`；本次 `I=C/common_prior_refit_input_v1`，上游全链路 `P=C/common_prior_component_refit_v1`。

## 工作内容与固定条件

使用100044、100116、100117各首个 RDO 文件 entries `[0,20)`，共60个源事件；精确文件名见上游报告。有限先验固定为 `matched_n32_f256_v1`：每个实际原生 seed 32个分量，分量内 covariance 为 C/256，保留所有相关项；源事件无 seed 时保留空输出。三个样本分别640、544、576个输入分量，首选 selected 分别311、254、360个。这些分量相关，不能当作独立物理事件或生产误判率分母。

实际 beforeGeo 源码为 `/eos/home-x/xcheng/FASER/geoRotation_calypso_component_refit_input_beforeGeo`，基底 `1eb885ae22299fb39bb352bead6005d41b0f953e`；afterGeo 继续使用 `/eos/home-x/xcheng/FASER/geoRotation_calypso_component_refit_afterGeo`，基底 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。beforeGeo 的[0037增量](../patches/0037_component_refit_input_beforeGeo.patch)保留原 `fit()` ABI，仅为 component 入口增加 `fitAudited()`，记录真正传入 native fitter 的 mean/C/surface 及实际 refit states。构建输出 `I/software_beforeGeo`，组合库 `I/overlay_beforeGeo`；afterGeo 组合库沿用 `P/overlay_afterGeo`。完整 patch、新增头文件和 reverse-check 见 `I/source_validation.json`。

本次15个实际重建作业分两步：

- `I/replay_v1`：beforeGeo 三样本×native/regular refit，共6 jobs；新增观察记录与原 full-chain 逐字段比较。
- `I/scale1_v1`：beforeGeo 两种 refit chart 与 afterGeo native，共9 jobs；`SeedCovarianceScale=10→1`，其余输入和设置保持不变。

两版实际 `KalmanFitterTool.cxx` 均在 `getParametersFromTrack()` 对 `cov(i,i)` 乘 `SeedCovarianceScale`，并在 `fit()` 中对完整 C 再乘10；beforeGeo 后一步硬编码，afterGeo 用默认值10的 `RefitCovarianceMultiplier`。**完整 C×10 不是两个版本之间的差异。** 两版还将 refit plane 向上游移10 mm，并从 Trk covariance 构造 ACTS bound C；本实验保留这些操作。

## 验证与物理量定义

`I/analysis_v1/validation.json` 的6个 observer gates 证明新增记录前后 seed、selected、全部 selected states、conversion、final 和 primary truth 全部精确相同。9个 scale gates 证明 CKF 与 conversion 全字段不变；refit 输入 bound mean 和 surface 精确相同，C 的非对角项不变，仅对角项除10。实际输入还逐项核对 Trk 5×5 covariance 的复制、MeV/GeV 单位转换、对角 inflation 和完整×10；最大绝对矩阵误差约6.78×10⁻²¹。

记录的 bound mean/C 通过解析 Jacobian 表达为共同 FASER frame 的 `(position[3], unit direction[3], q/p, time)`，无独立轨迹传播或拟合。矩阵差异指标为 `max_ij |Cb_ij−Ca_ij| / sqrt(D_i D_j)`，其中 `D_i=(|Cb_ii|+|Ca_ii|)/2`；该尺度化指标无量纲。direction 指标只取3×3方向子块。它不是 covariance 元素的普通相对百分比。

原输入 seed 经两个 Jacobian 表达后的尺度化 C 差异最大为1.12×10⁻⁷、2.01×10⁻⁷、2.37×10⁻⁷；它与有限分量在原 chart 中保持总均值/C 的约10⁻¹⁴ moment 误差是不同指标。分量内部的有限 Gaussian 非线性映射仍为一阶近似，当前结果不证明完整 posterior 收敛。

## 详细数据：实际 refit 物理 covariance

下表固定 same CKF hits、两版 refit 成功、两种 inflation 设置都满足的相同分量；N不随开关改变。数值为上述差异指标的中位数：

| DSID | N | 完整物理 C：scale10→scale1 | direction C：scale10→scale1 |
| --- | ---: | --- | --- |
| 100044 μ+ | 36 | 0.792992→0.286272 | 0.752848→0.277691 |
| 100116 μ+ | 24 | 0.737490→0.146619 | 0.737490→0.146619 |
| 100117 μ− | 62 | 0.282936→0.138297 | 0.269528→0.138297 |

![Paired native refit covariance](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/figures_v1/results/refit_input_covariance.png)

缩小对角 inflation 明显降低物理 C 差异，但不使其为零；CKF 的上游 mean/C 和部分路径原先已经不同，后续还存在 covariance basis 重建。两者的进一步因果分离未测量。

## 详细数据：两条 same-hit 符号差异

下表是 native beforeGeo / native afterGeo 的实际最终 q/p，单位 GeV⁻¹；100044 truth charge 为正：

| source entry / component | 原 scale10：before / after | scale1：before / after | scale1 结果 |
| --- | --- | --- | --- |
| 1 / 23 | −0.000652921137 / +0.000460927059 | −0.002286313868 / −0.002043966563 | 两版都错 |
| 3 / 19 | −0.000339919599 / +0.001104540612 | +0.000930912047 / +0.001333094463 | 两版都对 |

两例同 hits 且所有设置双方 refit 均成功。原 direction-C 差异分别0.760346、0.625081；scale1后分别0.0719740、0.177369。输入物理 mean 的版本间差别不随 scale 改变：

| entry / component | position norm 差 / mm | direction norm 差 | q/p 差 / GeV⁻¹ |
| --- | ---: | ---: | ---: |
| 1 / 23 | 0.001056223 | 2.433318×10⁻⁵ | 3.543280×10⁻⁵ |
| 3 / 19 | 0.000357845 | 2.973204×10⁻⁵ | 1.286523×10⁻⁵ |

![Actual signed q/p controls](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/figures_v1/results/refit_input_sign_cases.png)

另外两例100044 entry12/component2及100117 entry0/component25均在原始 CKF 中已经使用不同 hits；before refit 成功、after fallback，scale1 后仍有最终符号差异。不得将这些 fallback 当作双方 refit 比较。

三个样本的 native 版本间 sign differences 为 `3/0/1→1/0/1`；before regular-refit 对 after native 则为 `2/0/0→1/0/1`，并非联合 chart 与 scale1 消除了全部差异。native before/after 的 refit 成功数分别38/37、24/24、70/69，在 scale10/scale1 间不变；同样保留大量 converted-CKF fallback。

![Sign migration and fallback counts](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/figures_v1/results/refit_input_migrations_coverage.png)

## 执行、验证与产物

源码准备及运行任务由 [prepare_charge_component_refit_audit.py](../scripts/prepare_charge_component_refit_audit.py)冻结，统计由 [charge_refit_input.py](../analysis/charge_refit_input.py)执行，英文图组由 [charge_refit_input_figures.py](../analysis/charge_refit_input_figures.py)生成。批量通过 HTCondor，按当前兼容普通 schedd 的 running+idle 最低选择 bigbird20，排除有 EOS 兼容问题记录的 bigbird21；负载、时间、资源、submit 与 cluster 均保存在 worker manifest。

build1581774、after scale1 的1581777、before observer replay1581778、before scale1 的1581779均 exit0；统计1582083 exit0，归约3700条配对记录。15个实际 Athena 进程 CPU user+system 合计238.80 s，最大 RSS1,036,620 KiB，详见 `I/reconstruction_resources.json`；此合计不是并行 wall time。

导图及数值/PDF核对1582672 exit0。3张PNG已逐张原始分辨率目视检查，坐标、图例、英文说明及 fallback 标注可读，无裁切和文字重叠。逐图 PDF 与3页 PDF 的归一化 vector operations、font glyph outlines 和page boxes一致。图表程序独立从 paired rows 重新计算全部 cohort counts 与矩阵指标的max/median/p95，保存用于绘图的精简数据。67项单测通过，其中新增3项检验完整物理 C 的旋转、整体 inflation 的协变性及对角 inflation 的非协变性；日志 `C/verification_components_v1/unit_tests_with_refit_input.log`。

复现入口需使用新输出路径，实际批量命令、匹配库和任务以 `I/*tasks*.json` 为准：

```bash
python3 scripts/prepare_charge_component_refit_audit.py --root <root> \
  --phase tasks --run-name <new-replay-name>
python3 scripts/prepare_charge_component_refit_audit.py --root <root> \
  --phase tasks --run-name <new-scale1-name> --variants beforeGeo afterGeo \
  --refit-covariance-scale 1
python3 analysis/charge_refit_input.py --root <root> --replay <replay> \
  --scale-control <scale1> --output <new-analysis>
```

- [3页 refit 输入图组](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/figures_v1/results/refit_input_report.pdf)。
- [权威统计](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/analysis_v1/summary.json)、[3700条配对](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/analysis_v1/paired_refit_inputs.json)、[数值验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/analysis_v1/validation.json)。
- [绘图数据与索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/figures_v1/results/figure_index.json)、[导图验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/figures_v1/results/validation.json)、[目视记录](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_refit_input_v1/figures_v1/results/visual_review.json)。
- [118页 S06 汇报图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v12/charge_stage06_report.pdf)，组装1582673 exit0；前115页的index、vector content与page boxes精确不变，新增3页继承本次原始尺寸目视验收，详见[衔接验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v12/lineage_validation.json)。另核对上一版PDF实际渲染的112/114页，CDF及monospace计数表可读。

## 阶段性结论

事实：对角 inflation 是本批两条 same-hit refit sign 分歧的可控来源；统一整体×10在两版相同。取消额外对角 inflation 后，差异缩小、两例不再分歧，但最终正确性不一致改善。尚存上游 CKF mean/C/path 差及两例 path/fallback sign 差异，完整物理 posterior 等价未验收。

## 下一步计划

在实际 beforeGeo CKF 中复用已验证的 regular direction update/smoothing adapter，分别运行 OFF、update-only、smoothing-only、both；保留相同 component 输入、材料、measurement 和 refit。先用 OFF 重放精确验证新构建未改变原记录，再比较首个 physical mean/C 分叉、同 hit 与不同路径以及最终 sign/fallback。若剩余差异进入 refit，固定共同物理 refit prior 分离输入与 refit 运算；生产修复仍须独立 μ±/A′回归。
