# S06：完整前向状态的 covariance 分叉定位

**原先在 smoothing 中看到的 μ− 极端差异，已有可定位的前向来源：entry 1、component 24 在约134→1194 mm的原生传播段，physical direction差从4.79×10⁻¹²增至2.86×10⁻⁶。全部14,157对states的raw material slab相同。后续update的归一化C指标增加不等于绝对误差等比例增加；完整固定输入传播和posterior检验继续执行。**

**位置诊断更正已完成：`analysis_v2`将被动surface的loc0/loc1统一按平面x/y转换，该假设不成立。权威更正为`analysis_v3`，使用实际Calypso的Plane/Disc类型覆盖8,986对states，其余5,171对位置未测量；旧归档相关位置指标撤回。方向结论不变，逐项重算差异及真实传播结果见[完整C transport报告](2026-11-14_stage06_bound_transport.md)。下表index7位置已按typed转换更新。**

## 阶段目标

延续[真实smoother审计](2026-11-14_stage06_smoothing_operations.md)，将检查范围从measurement扩展到全部保存的被动材料面，沿ACTS真实state links定位前向mean/C差异，区分继承输入、保存的transport Jacobian和其余covariance增量。验收包括完整路径覆盖、共同surface检查、相同定义下的数值/分布图及可复核的极端路径。

结果目录`F=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/forward_covariance_v1`；权威更正统计为`F/analysis_v3`，图为`F/figures_v2/results`。原`analysis_v2`及`figures_v1`保留作历史记录。

## 工作内容

读取前阶段真实Calypso `smoothing_operations_v1/runs/<DSID>/<variant>/on/prior_audit.jsonl.gz`，不重新生成tracking状态。源码、matching libraries、输入与observer的完整重放gate沿用前报告；本次新增脚本只做坐标表达、逐状态配对、矩阵记账和绘图。

100044、100116、100117各首个RDO entries `[0,20)`，共60个独立源事件；n32、within-C/256，after使用解析seed Jacobian；isolated smoothing始终开启。每个有限分量独立调用实际Calypso finder/refit，这些分量相关，不是新的生产误判率分母。各分量按measurement数、chi2、candidate排序取首选，未选中/终端单state类别全部保留在计数中。

前向顺序直接反转实际`trackStatesReversed()`，逐项检查`acts_previous_index`连接。部分路径有signed path length递减，不能据此重新排序或丢弃；三样本分别有6、4、19条负path-step，均保留。归约v1错误地要求path length单调而中止；v2改用真实links，不改输入或拟合值。

每对surface完整4×4变换验证`S_after=R·S_before`，包括原生local basis。共同bound表达为`(local0,local1,phi_A,latitude_A,q/p,time)`：before用与`RefitDirectionChart.h`一致的解析Jacobian；after用`latitude_A=pi/2−theta_A`及对应符号矩阵。方向2×2 C指标与Cartesian unit-direction 3×3 C指标分别保存。

## 详细数据：覆盖与整体分布

| DSID | 输入分量 | 两版共同首选路径 | 单state路径 | measurement states | 被动states | 合计states | raw slab不同 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
|100044|640|311|7|1,584|2,573|4,157|0|
|100116|544|254|0|1,156|2,424|3,580|0|
|100117|576|360|2|2,158|4,262|6,420|0|
|合计|1,760|925|9|4,898|9,259|14,157|0|

每个输入分量均在两版同时选中或同时未选中，无one-sided selection；925条共同首选路径全部same ordered hits及完整state signature。5,081对states带surface material，其余均无材料查询字段。raw slab比较的是surface material查询返回的完整X0厚度和mm厚度，不把相同raw slab等同于相同全部材料作用；实际path correction及前后材料分配仍可能不同。

归一化指标为`max |Cb_ij−Ca_ij|/sqrt(d_i d_j)`，`d_i=(|Cb_ii|+|Ca_ii|)/2`，方向取共同chart的2×2子块，无量纲。

| DSID | predicted中位数 / 最大值 | filtered中位数 / 最大值 | smoothed中位数 / 最大值 |
| --- | --- | --- | --- |
|100044|5.603×10⁻¹¹ / 5.344×10⁻⁶|5.734×10⁻¹¹ / 5.346×10⁻⁶|1.125×10⁻¹⁰ / 7.615×10⁻⁶|
|100116|4.919×10⁻¹⁰ / 3.745×10⁻⁸|5.025×10⁻¹⁰ / 3.745×10⁻⁸|6.537×10⁻¹⁰ / 3.552×10⁻⁸|
|100117|2.614×10⁻⁹ / 4.331×10⁻⁴|2.548×10⁻⁹ / 4.331×10⁻⁴|2.761×10⁻⁹ / 1.940×10⁻³|

![Complete native-state covariance distributions](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/forward_covariance_v1/figures_v1/results/forward_covariance_distributions.png)

## 详细数据：前向分叉与covariance记账

定义保存的矩阵余额`E=C_pred−J C_prev Jᵀ`。每端均将前一filtered C、当前predicted C和J两端解析变换到共同chart，再作以下恒等式：

```text
ΔC_pred = Jb (Ca−Cb) Jbᵀ
        + [Ja Ca Jaᵀ − Jb Ca Jbᵀ]
        + (Ea−Eb)
```

下标b/a为before/after。第一项为前一C差异的收缩，第二项为保存J变化的收缩，第三项为剩余矩阵余额。它们是对原生输出的代数记账，不是重新运行tracker的反事实实验。E不能称为直接观测的native Q，因为材料更新和其他没有编码在保存J中的操作均可能进入余额。全部14,157个恒等式的最大标准化残差为2.1532×10⁻¹⁵。

100117、entry1、component24关键位置如下，path为原native保存的signed path length：

| forward index | path / mm | measurement | predicted方向C差异 | filtered方向C差异 | physical predicted方向差 | predicted位置差 / mm |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
|0|3.888|否|1.726×10⁻¹⁵|1.726×10⁻¹⁵|2.321×10⁻¹⁶|0|
|1|10.000|是|2.336×10⁻⁸|2.336×10⁻⁸|1.230×10⁻¹⁶|1.421×10⁻¹⁴|
|7|133.934|否|7.193×10⁻⁹|7.193×10⁻⁹|4.791×10⁻¹²|7.597×10⁻¹⁰|
|8|1193.919|否|5.408×10⁻⁶|5.408×10⁻⁶|2.858×10⁻⁶|1.528×10⁻³|
|9|1208.862|是|5.365×10⁻⁶|1.463×10⁻⁴|2.859×10⁻⁶|1.570×10⁻³|

index8的source surface为geometryID`1873497582425080064`、z_F=137.401 mm，target为`1945555176731443200`、z_F=1197.3575 mm。index8余额差异仅1.013×10⁻¹⁵，保存J变化收缩为5.415×10⁻⁶，前一C收缩为7.179×10⁻⁹。这把可判别的下一实验定位为该段固定完整mean/C/surface的原生transport及步长控制。当前两版输入mean/C仍有极小继承差异，因此不能直接把全部差值认定为相同物理输入下的旋转缺陷。

首个材料面之后index1已经有约2.34×10⁻⁸的C差异，主要在余额项；该位置mean仍近机器精度一致，支持另测material前后分配及噪声坐标表达。

index9的measurement update后归一化C差异约增加27.3倍，但方向C最大元素的绝对差只从4.2602×10⁻¹²增至4.5244×10⁻¹² rad²，Frobenius差只从6.5528×10⁻¹²增至6.8377×10⁻¹² rad²。不能把归一化指标的增长说成等倍绝对误差放大。后向index8 smoothed归一化差为1.940×10⁻³，但其绝对差1.0767×10⁻¹² rad²反而小于本层predicted绝对差。

![Extreme path with absolute and scaled covariance differences](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/forward_covariance_v1/figures_v1/results/forward_covariance_extreme.png)

## 验证、运行与产物

实现为[charge_forward_covariance.py](../analysis/charge_forward_covariance.py)、[绘图程序](../analysis/charge_forward_covariance_figures.py)、[任务冻结入口](../scripts/prepare_charge_forward_covariance.py)。新增5项坐标/矩阵记账/state-link检查、原有6项有限prior检查和1项PDF resource alias检查通过。新PDF检查修正了跨行`/A2\ngs`资源别名的解析；数值图和PDF均未修改，单页/汇编内容在解析同一资源后相等。

HTCondor归约v2 cluster1582714；绘图1582716；验证解析修正后1582718；汇编1582726，最终均exit0。提交前比较兼容schedd负载和队列，均选择bigbird20，排除已有EOS兼容记录的bigbird21。首轮1582712因path-length单调断言退出；首轮绘图中的导图成功，验证因跨行资源alias未解析失败，修正后只重做验证。归约CPU59.50 s、最大RSS429,852 KiB；绘图CPU9.75 s、最大RSS304,088 KiB。各任务原文、请求资源和运行日志在F中。

```bash
python3 scripts/prepare_charge_forward_covariance.py \
  --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1 --revision 2
```

冻结tasks为`F/analysis_tasks_v2.json`、`figure_sequence_tasks_v1.json`和`figure_validation_tasks_v2.json`，经`scripts/submit_momentum.py`提交。权威[summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/forward_covariance_v1/analysis_v2/summary.json)、[配对验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/forward_covariance_v1/analysis_v2/validation.json)、[三页图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/forward_covariance_v1/figures_v1/results/forward_covariance_report.pdf)与[导图验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/forward_covariance_v1/figures_v1/results/validation.json)均保留。三张PNG已按原始分辨率目视检查。

[133页S06汇编](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v17/charge_stage06_report.pdf)继承report_v16的130页；前缀索引、decoded vector stream和page boxes逐页精确一致，验证见[lineage_validation.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v17/lineage_validation.json)。原生产参考分布继续保留，本次追加三页逐状态诊断。

## 阶段性结论

事实：raw material slab在本批paired states相同；μ−极端mean/C差异在前向传播已出现，随后被update和smoothing继承。共同chart指标、Cartesian物理C和绝对差分别保存，不能把指标变化混为同一效应。多数states的残差很小，但极端路径保留，未凭中位数或sign一致宣布Q5闭合。

假设：长段积分步长/数值路径和首材料面噪声表达分别可能解释上述两个增长点。固定完整输入的实际Calypso传播尚待测量，不能以本次代数记账确认其因果。

## 下一步计划

1. 以真实source/target surfaces和完整mean/C运行Calypso extrapolator；固定before物理输入与两端各自native输入分别比较，保留失败类。
2. 材料ON/OFF、原步长/100 mm/10 mm分别比较；先完成扩展probe的旧输入逐值重放gate，再执行新模式。
3. 按结果记录原生J、材料过程项或更细传播步骤，继续有限Gaussian分布映射、likelihood加权posterior及宽度/分量数收敛；S06 Q5保持开放。
