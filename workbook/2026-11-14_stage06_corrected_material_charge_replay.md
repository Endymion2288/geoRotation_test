# S06：修正材料编码后的119事件错号机制复核

**修正材料图后，原59条错号的机制分类保持：56条最后在CKF measurement update进入持续到final的错误符号，3条从seed保留错误符号。56次最后穿零的hit也逐一相同，其中31个没有SDO、25个只有primary支持。全部119个事件的完整seed、selected ordered hits及共同measurement坐标相同；最终q/p数值全部改变。因此此前逐运算归因没有被这项材料修正推翻，但不能将两种材料下的数值trace视为相同。**

## 阶段目标

在[150万事件全量材料编码对照](2026-11-14_stage06_material_encoding_reconstruction.md)确认μ±仍是相同59个wrong之后，重新检查真实Calypso的seed、CKF逐state update、native target fit、conversion和refit。验收要求包括独立复现corrected full shard、候选/seed身份、逐hit支持以及矩阵恒等式，不由最终sign相同推断内部机制相同。

根目录`C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`，本轮根目录`M=C/material_encoding_control_v1`。

## 工作内容与冻结配置

目标集合完全沿用`C/charge_replay_v1/targets.json`：100116 μ+的27个wrong、100117 μ−的32个wrong，另有原先冻结的每电荷30个正常控制。正常控制按原文件分层的modulus-997规则选择；未按本轮outcome重挑事件。所有输入以RDO路径与entry识别。

只将原详细replay task的材料路径换为`M/corrected_afs_map.json`，输出换到`M/replay_v1`，逐数组精确参考换为`M/<DSID>/full/shards/afterGeo/sNNNN`。runner、worker、实际full-audit library、input/range及options均复用原始冻结版本：

```text
C/overlay_full_audit_afterGeo/lib
C/charge_replay_v1_code/scripts/run_charge_replay.py
--no-seed-writer --actual-covariance --isolate-smoothing --stage-audit --charge-audit
```

source task保留旧target的q/p作为provenance，归约不使用它强制新的q/p。每个encoding内均从该次observer的charge-independent最佳matched track出发，通过完整ordered hit IDs及实际final q/p唯一定位audit final record，再用其seed/candidate定位selected trace。分类序列按实际前向P/F顺序，之后是native CKF target reference、conversion和final refit；候选隔离smoothed states与native target fitted mean分开。

本轮没有重建独立tracking，也没有用新定义的选样代替原misID分母。119事件是已知outcome的机制样本，比例不解释为总体误判率。

## 详细数据：精确重放与身份

两条wrong smoke先与各自corrected full shard的58个数组字段精确闭合，再提交剩余117条。全部119条真实重建exit0、58字段逐值闭合；event、truth、tracks、segments都在比较范围内。每次运行的实际audit/isolation/charge-audit开关为ON。

| 集合 | 事件 | 完整seed相同 | selected ordered hits相同 | 首次wrong stage相同 | 最后进入final wrong的stage相同 | final q/p数值改变 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| μ+ wrong | 27 | 27 | 27 | 27 | 27 | 27 |
| μ− wrong | 32 | 32 | 32 | 32 | 32 | 32 |
| μ+ normal | 30 | 30 | 30 | 30 | 30 | 30 |
| μ− normal | 30 | 30 | 30 | 30 | 30 | 30 |

完整seed检查包括六维bound mean、6×6 covariance、surface、seed ID及hit vector顺序，均逐值一致。两版各1816次selected measurement update的共同measured coordinates全部相同；原RDO及cluster配置不变。per-hit SDO来自原compact measurement audit，在当前相同hit路径上全部1816次都有覆盖；`not_in_legacy_measurement_audit`为0，没有把未审计hit当作无SDO。

两个normal集合仍全部correct。部分normal在中途曾短暂进入错误符号后恢复，因此它们的`last_update_into_wrong`存在不表示final wrong；主机制计数始终按最终wrong集合及其最后持续穿零定义。

## 详细数据：最后进入错误符号的操作

| charge | wrong总数 | seed保留 | CKF measurement update | 最后穿零hit相同 | 穿零hit无SDO | 穿零hit仅primary支持 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| μ+ | 27 | 1 | 26 | 26/26 | 15 | 11 |
| μ− | 32 | 2 | 30 | 30/30 | 16 | 14 |
| 合计 | 59 | 3 | 56 | 56/56 | 31 | 25 |

legacy/corrected两套计数相同。本集合未出现最后由conversion或refit引入的新增错误符号；不将其推广到其他beforeGeo或direction/station配置。conversion在每条新trace上独立核对native q/p到MeV⁻¹ branch的单位转换和符号一致。

实际记录的`r,S,K,C`重新计算满足`Δ(q/p)=(K r)_(q/p)`：legacy/corrected各1816次update的最大绝对恒等式误差分别为`1.96041×10⁻¹⁴`与`4.22566×10⁻¹⁴ GeV⁻¹`。该恒等式闭合证明使用了实际记录的更新运算，不能证明误差模型或hit假设正确。primary-only的SDO支持也不能单独证明该轨迹假设正确；后者需要既有hit/centroid/fixed-seed控制和进一步局部轨迹证据。

## 原六条pure wrong在新材料下的状态

下表是此前primary-SDO过滤后仍错、并已做centroid/fixed-seed控制的冻结六条。当前表只改变材料编码，没有重新施加centroid或hit过滤，单位GeV⁻¹。

| DSID/sample entry | legacy final q/p | corrected final q/p | 完整seed及ordered hits |
| --- | ---: | ---: | --- |
| 100116/51839 | −0.001534125092 | −0.001534124330 | 相同 |
| 100116/224394 | −0.001469602076 | −0.001469602131 | 相同 |
| 100117/227965 | +0.000264475596 | +0.000264435713 | 相同 |
| 100117/390966 | +0.000407198602 | +0.000407252880 | 相同 |
| 100117/439421 | +0.002789309002 | +0.002596273726 | 相同 |
| 100117/458033 | +0.000327301885 | +0.000327301812 | 相同 |

439421的q/p幅值变化较其他五条明显，但仍错号且没有换seed/hits。材料会改变数值响应，不能因总体率接近就宣称内部q/p不受影响。centroid固定seed使五条selected恢复、另一条仅rmax失败的既有证据见[固定seed报告](2026-11-14_stage06_fixed_seed_alignment.md)；该旧centroid实验未在本轮新map下重新测量。

## 执行、验证与产物

smoke重建1578758、gate1578760、剩余117重建1578761、full validator1578763全部exit0。归约smoke1578762及full1578764均exit0，完整119条候选定位与逐运算检查通过，无重跑。每次提交重新比较兼容schedd负载，排除bigbird15/21，选择bigbird20；实际资源、库入口和source task随`*_worker/submission.json`及每事件`run.json`保留。请求1 CPU，tracking3500MB/4GB disk，归约3500MB/3GB disk，线程数与申请一致。

实现：[任务冻结与58字段门槛](../scripts/prepare_material_encoding_replay.py)、[候选定位与逐运算归约](../analysis/material_encoding_replay.py)、[归约任务冻结](../scripts/prepare_material_replay_analysis.py)、[机制图与数据核对](../analysis/material_encoding_replay_figures.py)、[绘图任务冻结](../scripts/prepare_material_replay_figures.py)。

```bash
python3 scripts/prepare_material_encoding_replay.py --root "$C" --phase prepare
python3 scripts/prepare_material_replay_analysis.py --root "$C" --revision <new-analysis>
python3 scripts/prepare_material_replay_figures.py --root "$C" --analysis <validated-analysis> --revision <new-figures>
python3 scripts/submit_momentum.py --tasks <frozen-tasks.json> \
  --library "$C/overlay_full_audit_afterGeo/lib" --name <unique-name> --stage S06-material-replay \
  --cpus 1 --memory-mb 3500 --disk-gb 4 \
  --exclude-schedd bigbird15.cern.ch --exclude-schedd bigbird21.cern.ch
```

[119事件完整验收](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/replay_v1/full_validation.json) · [机制summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/replay_analysis_v1/full/summary.json) · [逐事件r/S/K/C与身份](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/replay_analysis_v1/full/events.json)。

![Corrected-material charge origins](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/replay_figures_v1/results/material_replay_origins.png)

绘图1578765、独立validator1578767与汇编1578770均exit0。两张PNG已全尺寸目视检查；逐图与两页PDF的vector内容、page boxes及rotation精确一致，plot_values与119事件归约逐项相等。新增两页合入此前96页，形成[98页完整S06图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v7/charge_stage06_report.pdf)，汇编保持98页原始vector内容及PNG完整性，并继承全部来源图的目视记录。[两页机制PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/replay_figures_v1/results/material_replay_report.pdf) · [新图验收](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/replay_figures_v1/results/validation.json) · [汇编验收](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v7/validation.json)。

## 阶段性结论与下一步计划

事实：材料编码修正没有改变本批59个wrong的seed/measurement来源分类、同一最后穿零hit或原119个所选hit路径；q/p数值及更新响应会变化。修正map下的错号仍需从seed假设、measurement关联/位置及其误差模型解释。31次无SDO与25次primary-only是穿零hit的观测分类，不是将59条按唯一物理根因作穷尽分配。

下一步以修正map推进完整共同物理先验的真实CKF控制，以及六条pure事件的measurement error、strip量化及局部G4轨迹拆分。transport闭合、实际Kalman恒等式和总体率变化很小均不能代替这两项。μ±/A′留出模式的联合refit回归已在[独立回归报告](2026-11-14_stage06_joint_refit_regression.md)交付；S06继续保持未完成。
