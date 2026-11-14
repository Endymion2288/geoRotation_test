# S06：μ±与A′的联合refit回归

**冻结的regular update＋smoother控制在两类μ和两类A′的全部六种模式中减少错号计数。110051三站forward的misID从1212/13589降至297/13648，correct-charge效率增加4.9722个百分点，95%配对区间[4.6127,5.3316]；相对production p的半动量负尾36.7577%→9.8916%，角残差68%绝对分位171.609→41.221 μrad。四站μ样本的变化较小，完整分母、选择迁移、fallback与CPU同时保留。**

## 阶段目标

完成[Q5修复回归](../roadmap/04_charge_misidentification.md)中已冻结的beforeGeo联合refit控制，覆盖μ±、A′、forward/backward/merge、3/4站；同时核对charge、momentum、angle、efficiency、fallback和计算资源。本次不修改afterGeo，也不将原生版本差解释为纯坐标效应。S06的完整共同物理先验与剩余measurement机制继续研究。

结果根目录`C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`。

## 工作内容与冻结输入

实际使用Calypso master `1eb885ae`的原生几何和tracking，加[0022/0027诊断与默认关闭的refit控制](../patches/README.md)。`native`与`both`固定相同RDO、几何、field、material和CKF isolated smoothing ON，只切换既有`RegularRefitUpdateChart`与`RegularRefitSmoothingChart`。新任务通过`C/overlay_refit_before/lib`加载`C/software_refit_before`的三份tracking/geometry库及`C/observer`；复用native任务的原库路径保存在各run manifest，逐job核对。

| DSID | 输入file indices | 每文件entries | 源事件数 |
| --- | --- | --- | ---: |
| 100116 μ+ | 3–4 | [0,20000) | 40000 |
| 100117 μ− | 3–4 | [0,20000) | 40000 |
| 110033 A′、ε=10⁻⁵ | 5–9 | [0,2000) | 10000 |
| 110051 A′、ε=10⁻⁴ | 10–19 | [0,1000) | 10000 |

每类运行`3forward/3backward/3merge/4forward/4backward/4merge`的native/both，共48组、累计1.2M次重建事件。96个已验收μ-native shards复用，336个任务新重建，合计432个有效任务。新重建720k事件、复用480k事件。这是冻结补丁的保留范围回归；此前已查看native结果，μ三站forward也已有报告，不称为新的blind discovery。110033这五文件未进入S05输入；110051这十文件仅代表本次回归范围，不标为全20文件。

发现merge配置此前只向主CKF的refit设置控制，而第二个backward CKF没有设置。现由[configure_refit](../configs/slope_tracking.py)对实际两个方向分别应用，并保存`actual_refit_controls_by_direction`。60个smoke任务全部成功；24组默认58字段逐值相同、48组实际开关与任务要求一致，再运行全量。

matching、purity严格`>0.7`、Good95及fiducial沿用S06 charge-independent选择；每truth先按固定质量规则排序，再判quality。A′选HepMC status=1、barcode在(0,200000)、abs(PDG)=11的电子。残差动量参考HepMC production p；角度参考Calypso将production-p/首个ST1 SiHit chord模型传播至共同`z_F=0 mm`的方向，**不是同sensor局部G4动量或方向**。电子同源两个daughter的相关性按源事件聚类保留。

## 详细数据：各模式错号与分母

数值为`wrong / selected matched truth`。F/B/M分别表示forward/backward/真实merge。

| DSID | 模式 | native | regular update＋smoother |
| --- | --- | ---: | ---: |
| 100116 | 3F | 35/20084 | 2/20076 |
| 100116 | 3B | 48/20328 | 6/20341 |
| 100116 | 3M | 23/20656 | 3/20667 |
| 100116 | 4F | 13/18406 | 11/18444 |
| 100116 | 4B | 22/19940 | 12/19994 |
| 100116 | 4M | 11/20254 | 9/20300 |
| 100117 | 3F | 68/20756 | 4/20747 |
| 100117 | 3B | 95/21079 | 4/21107 |
| 100117 | 3M | 54/21387 | 4/21391 |
| 100117 | 4F | 38/19146 | 31/19154 |
| 100117 | 4B | 51/20657 | 20/20710 |
| 100117 | 4M | 20/20993 | 12/21006 |
| 110033 | 3F | 87/12730 | 41/12704 |
| 110033 | 3B | 123/12563 | 64/12608 |
| 110033 | 3M | 75/13207 | 38/13247 |
| 110033 | 4F | 87/12730 | 41/12704 |
| 110033 | 4B | 125/12555 | 68/12600 |
| 110033 | 4M | 76/13205 | 37/13245 |
| 110051 | 3F | 1212/13589 | 297/13648 |
| 110051 | 3B | 1319/12563 | 333/12740 |
| 110051 | 3M | 1194/14021 | 307/14205 |
| 110051 | 4F | 1212/13589 | 297/13648 |
| 110051 | 4B | 1325/12522 | 345/12701 |
| 110051 | 4M | 1193/14020 | 304/14210 |

各组完整truth表和segments表逐值相同；7×7选择迁移、truth-p分箱、ST1四象限、truth charge分组、重复候选及未匹配/质量失败/非法曲率均存于[全量summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/joint_regression_full_analysis_v2/summary.json)。各模式在其自身共同母样本中配对；不同模式的重复运行不合并为独立统计。

三站forward的主要配对差值如下，单位为百分点，Δ=both−native；区间是95% source-event cluster delta区间。

| DSID | misID Δ [95%区间] | correct效率 Δ [95%区间] | 共同selected / 同candidate且ordered hits相同 |
| --- | --- | --- | ---: |
| 100116 | −0.1643 [−0.2203,−0.1083] | +0.1070 [+0.0334,+0.1806] | 20058 / 19985 |
| 100117 | −0.3083 [−0.3861,−0.2306] | +0.2270 [+0.1380,+0.3159] | 20724 / 20652 |
| 110033 | −0.3607 [−0.5016,−0.2198] | +0.1189 [−0.0425,+0.2803] | 12655 / 12534 |
| 110051 | −6.7428 [−7.2232,−6.2625] | +4.9722 [+4.6127,+5.3316] | 13502 / 13261 |

同candidate且ordered hits相同的子集另保存动量和角度分布；总体改善也包含quality/ranking迁移，不能将上表当作固定track的纯参数差。误差为0的经验重采样不再给[0,0]：使用注明条件为observed exposures的保守any-event边界上限。

## 详细数据：动量、角度与fallback

以下为三站forward，箭头为native→both。半动量负尾比例以各自selected为分母，角度68为`|Δtanθ_y|`的68%分位。

| DSID | p_reco/p_prod<0.5 / % | 角度68 / μrad | 角度绝对残差>1 mrad / % |
| --- | --- | --- | --- |
| 100116 | 1.1004→0.3835 | 75.646→48.386 | 0.80163→0.09962 |
| 100117 | 1.7489→0.4290 | 76.651→43.400 | 0.90094→0.08194 |
| 110033 | 11.1076→6.6200 | 150.149→62.420 | 2.27808→1.03117 |
| 110051 | 36.7577→9.8916 | 171.609→41.221 | 3.56170→0.64478 |

四站μ的方向参考含ST1模型与上游信息的差别，不能直接与三站的分位大小比较。固定四站forward，μ+角度68为624.151→625.260 μrad，大尾20.52048%→20.56495%；μ−为526.499→524.619 μrad、大尾18.28580%→18.24684%。这些小变化不支持“所有模式的每个指标都改善”。四站μ−correct效率增量+0.0619 pp的95%区间[−0.0072,+0.1310]包含0。

逐流的refit counters与最终tracks逐项闭合；真实merge的输出按hits、q/p和upstream参数与两个单方向来源匹配。3F、3B、4F、4B四次处理的全部fallback次数，在native/both之间分别相同：μ+为26290/147937、μ−24887/147570、110033为2277/69857、110051为912/88346；这些分母是重复配置中的refit calls，不是独立粒子。对应selected中的fallback次数也分别相同，为15、7、0、2；selected分母随选择变化。该控制并未通过减少全体fallback次数获得上述改善。

![110051 joint metrics](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v11/joint/110051_joint_metrics.png)

![110051 momentum response](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v11/joint/110051_relative_p.png)

[18页完整图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v11/joint/joint_regression_report.pdf)含四类样本各自的联合指标、momentum、angle、signed q/p，以及配对效率/misID和CPU。动量采用参考宏100 bins、[−1.2,1.2]及raw Entries；该图是最终track的六模式比较，不标为原宏逐sensor九层图。每个panel记录finite N和flows，图头注明文件及entry范围。

## 执行、资源与验证

全部批量处理经HTCondor，每次提交重新比较兼容schedd负载并查询所选队列，均选bigbird20，排除bigbird15/21；时间、选择依据和资源随各`*_worker/submission.json`记录。全量原cluster1578686共336 jobs，325成功、11在进入事件前导入awkward失败。worker将Python3.9额外user-site显式追加在release PYTHONPATH之后；20-event两组环境控制1578692证明58字段精确、release NumPy1.23.5/SciPy1.10.1优先，追加awkward2.8.12/uproot5.6.9。

仅恢复已终止的11个失败任务，cluster1578693全部exit0；原失败输出保留，恢复结果通过effective manifest替换，不覆盖或重跑325项成功结果。`joint_regression_effective_v1`的48组432任务逐range覆盖无遗漏/重叠、实际两方向开关正确。全量归约1578704、最终绘图1578707、独立图表验证1578709均exit0。

观测process CPU范围如下，已含初始化，运行worker不同，不能直接解释为相同CPU上的算法耗时比。

| DSID | native / ms每源事件 | both / ms每源事件 | 两组最大RSS / MiB |
| --- | --- | --- | --- |
| 100116 | 8.879–13.791 | 13.912–19.612 | 990.80 / 994.86 |
| 100117 | 8.695–14.544 | 13.957–17.659 | 989.73 / 987.58 |
| 110033 | 184.533–217.582 | 140.558–220.961 | 1094.60 / 1094.25 |
| 110051 | 965.983–1341.987 | 1192.831–1310.823 | 1258.39 / 1257.76 |

分析实现：[charge_joint_regression.py](../analysis/charge_joint_regression.py)；任务准备、范围/配置验收、失败恢复分别为[prepare](../scripts/prepare_charge_joint_regression.py)、[validate](../scripts/validate_charge_joint_regression.py)、[recover](../scripts/recover_charge_joint_regression.py)。只读覆盖验证命令为：

```bash
python3 scripts/validate_charge_joint_regression.py \
  --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1 \
  --phase full
```

实际分析输入为`C/joint_regression_effective_v1`，冻结执行代码为`C/joint_regression_analysis_v2_code/analysis`。图代码为`C/figures_v11/code/analysis`；[validation.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/figures_v11/joint/validation.json)证明224项summary对照、144个直方图含flows/sumw2、18页PDF vector/boxes与单图一致，PNG全部解码。18张contact sheet与110051两张主图已目视检查。54个unit tests通过，记录在`C/verification_v8/unit_tests.log`；新增三项joint测试覆盖相关daughter的误差、零计数边界和配对统计。refit双方向配置由实际smoke和全量run manifest检查。

## 阶段性结论

事实：该联合refit控制在冻结范围、全部模式中减少错号计数；三站μ与两类A′的动量负尾和角度指标同步改善，110051效应最大。correct效率总体增加，但个别区间包含0。全体fallback calls未减少，四站μ的角度变化较小，CPU需同硬件复测才能给算法开销比。

解释：结果支持beforeGeo near-pole refit的更新与smoothing参数化为先前额外错号的重要来源；该解释已有逐运算与独立full-G4证据，本轮补齐模式/样本回归。afterGeo剩余错号的measurement、branch和seed机制不由此控制解决；完整Gaussian先验的非线性坐标闭合及材料查询仍需独立验证。

## 下一步计划

1. native bin-query与共同有限物理先验transport已完成：旧JSON解码转置导致非闭合，修正后7,680点达到浮点精度，见[材料decoder报告](2026-11-14_stage06_physical_prior_material_decoder.md)。其真实CKF影响另作全量控制。
2. 在匹配材料及物理先验后继续CKF、measurement及模型误差控制，报告剩余类别发生率，不将transport sign保持等同于全部CKF posterior闭合。
3. 汇入S06总图集，完成Q1–Q5逐项验收；S06保持未完成。
