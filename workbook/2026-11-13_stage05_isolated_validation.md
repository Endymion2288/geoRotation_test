# S05：isolated smoothing 开启后的全量复核

**100044 每版完整500,000事件已用 isolated smoothing ON 重建。与旧 OFF 基线相比，最终动量、电荷、方向、传播结果和全部选样成员逐条相同；仅 beforeGeo 73个、afterGeo 268个 rmax 元素变化。st1/VetoNu 的全部比较直方图逐 bin 相同。** 因此，S05已定位的最终角度结构及其主要来源在本次要求的 ON 版本中仍成立。

## 阶段目标

落实用户要求的“使用开启 isolated smoothing 的 geoRotation Calypso”研究 `problem.pdf` p4 的角度问题。核对开关对最终物理量与选择迁移的影响，以 ON 的实际重建重新生成 `geometryRotation.pdf` p19–20 所要求的参考图。

本报告是[角度来源结论](2026-11-13_stage05_angle_structure_conclusions.md)与[A1–A5验收](2026-11-12_stage05_acceptance.md)的补充；之前的完整实验记录保留各自实际开关，不追溯修改其 provenance。

## 工作内容

结果根目录为 `/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`，下文记为 `C`。阶段标签对应 S05，目录置于本次 ON 批次下，避免混淆旧冻结批次。

输入为 DSID100044 的单个完整500k RDO：

```text
/eos/experiment/faser/data0/sim/mc24/particle_gun/100044/rdo/s0013/FaserMC-MC24_PG_mupl_fasernu_5mrad_flukaE-100044-00000-00099-s0013-RDO.root
```

beforeGeo 基底为 master `1eb885ae22299fb39bb352bead6005d41b0f953e`，afterGeo 为旋转分支 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。全部 tracking/propagation 由匹配的 Calypso 实际组件执行；两版均开启 `IsolateCandidateSmoothing=True`。Athena24.0.41、FASERNU-04、OFLCOND-FASER-05，材料分别为 source/current map 与用户指定 AFS rotated JSON。原生 seed covariance 差别继续保留，因此版本差别仍不称为相同物理先验的纯旋转效应。

ON 数据为 `stage05_slopes/isolated_100044_v1`，OFF 为 `stage05_slopes/main_100044_v1`。每版100 shards，每 shard5000个连续 entries。`analysis/slope_isolation.py` 按同一文件、entry、track ordinal，比较events/truth/tracks/segments全部56个公共字段，并直接比较 all/reco-p-window/Good95 的通过集合。

## 详细数据

| 项目 | beforeGeo OFF / ON | afterGeo OFF / ON |
| --- | ---: | ---: |
| 源事件 |500,000 / 500,000|500,000 / 500,000|
| final tracks |484,845 / 484,845|484,931 / 484,931|
| 三站 long tracks |436,859 / 436,859|432,124 / 432,124|
| reco p 50–250 GeV |176,702 / 176,702|157,132 / 157,132|
| Good95 |356,951 / 356,951|353,016 / 353,016|
| 三种选择逐条成员变化 |0 / 0 / 0|0 / 0 / 0|
| `tracks/rmax_mm` 变化元素 |73|268|
| 其余公共字段变化元素 |0|0|

未将“数量相同”当作“成员相同”；实际 membership 对比全为0。rmax 的少数变化未跨越本次Good95阈值。最终 q/p 的新增观察字段不在旧 OFF 文件的公共字段中；已有最终 momentum 与 charge 逐元素完全相同。

物理定义沿用 S05：`sy=py_F/pz_F=py_A/px_A`，st1固定`z_F=0 mm`、VetoNu固定`z_F=-3112 mm`；production p 与 SiHit chord truth model 区别标注。峰窗口复用冻结的 `stage05_slopes/main_analysis_v4/results/peak_windows.json`，未重新按 ON 图调窗。

完整逐字段、逐选择、逐bin证据：

- [beforeGeo验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/s05_validation_v2/beforeGeo/summary.json)。
- [afterGeo验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/s05_validation_v2/afterGeo/summary.json)。
- [ON参考图索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/s05_isolated_reference/figures/index.json)与[57页PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/s05_isolated_reference/figures/slope_reference_report.pdf)。包括physical/legacy、同面truth/residual、条件分组与response。

![Isolated smoothing ON: st1 physical slope reference](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/s05_isolated_reference/figures/100044_st1_physical_reference.png)

![Isolated smoothing ON: VetoNu physical slope reference](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/s05_isolated_reference/figures/100044_VetoNu_physical_reference.png)

归约 cluster1574107 的2/2任务exit0，ON参考绘图1574108 exit0；全量重建的逐shard状态由上述输入任务及coverage记录闭合。全部57张PNG的contact sheet已目视检查，主要st1/VetoNu参考图另行检查，未见布局截断。bin counts和entries随参考图及验证summary保存。

复现入口：

```bash
python3 scripts/validate_slope_isolation.py --help
python3 analysis/slope_isolation.py --help
```

实际冻结命令见 `C/s05_validation_v2_tasks.json`、`C/s05_isolated_reference_tasks.json`，源码快照和提交时的schedd比较保存在对应worker目录。批量运行使用HTCondor，匹配CPU资源与单线程重建；未在交互节点重做完整样本。

## 阶段性结论

事实：candidate-isolated smoothing修正的是逐候选states；本样本中最终角度、动量、电荷和选择集合没有变化。原S05的 near-pole refit update/smoothing 来源、afterGeo strip/response结构及VetoNu散射解释均与本次ON最终量一致。新的ON参考图可用于当前汇报。

事实：这项ON/OFF验证没有重新测量原S05各个诊断控制的所有中间states。原控制的配置和直接因果证据继续引用原报告；最终物理量不变的全量验证不替代新修复的独立效率/pull验收。

## 下一步计划

S06继续用两版ON库完成μ±全量率、谱/接受度控制及逐运算错判定位。S05完整生产修复仍按S07验收新的保留输入、效率、误差校准与资源；既有联合chart的约2个百分点独立gun效率损失及pull余量保留在结论中。
