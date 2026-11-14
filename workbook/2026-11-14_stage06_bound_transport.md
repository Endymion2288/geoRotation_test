# S06：完整 covariance 的原生 surface-to-surface 传播

**真实Calypso传播支持两项区分：固定相同物理mean/C后，四条极端路径的局部传播mean在浮点精度内一致；保留两版各自的微小输入差异时，约134→1194 mm长段对步长明显敏感。固定输入仍有约10⁻⁸的方向C差异，关闭材料后降到约10⁻¹⁴。内部CKF逐段传播及完整posterior尚未闭合，S06 Q5继续开放。**

## 阶段目标

延续[前向状态审计](2026-11-14_stage06_forward_covariance.md)，在真实Calypso source/target surfaces之间注入完整6D mean和6×6 covariance，区分继承输入、材料与积分步长的作用。同时更正将被动面loc0/loc1一律视为平面x/y的位置诊断。

结果目录`T=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/bound_transport_v1`；权威传播分析为`T/analysis_v2`，图为`T/figures_v1/results`。更正后的完整前向分析为相邻`forward_covariance_v1/analysis_v3`，图为其`figures_v2/results`。

## 工作内容与软件

传播调用实际`IFaserActsExtrapolationTool::propagate()`，材料调用另由同工具的`propagationSteps()`记录。没有另写propagator、measurement update或smoother。

| 项目 | 实际来源 |
| --- | --- |
| Probe源码 | `/eos/home-x/xcheng/FASER/geoRotation_calypso_bound_transport_v1/Tracking/Acts/FaserActsKalmanFilter/src/CalypsoSlopeDiagnosticsAlg.cxx` |
| Probe基底 | `geoRotation_calypso_material_bins_v4`中的既有真实Calypso诊断组件，HEAD `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` |
| 新增补丁 | `0040_full_covariance_native_transport_probe.patch` → `0041_native_passive_surface_transport_lookup.patch` → `0042_native_transport_surface_types.patch` |
| 实际传播实现 | `geoRotation_calypso{,_main}/Tracking/Acts/FaserActsGeometry/src/FaserActsExtrapolationTool.cxx` |
| beforeGeo几何/场库 | `stage05_geometry_comparison/software_before/lib/libFaserActsGeometry{,Lib}.so`，master `1eb885ae22299fb39bb352bead6005d41b0f953e` |
| afterGeo几何/场库 | `stage06_charge/20261113_isolated_v1/physical_prior_transport_v2/software_afterGeo/lib/libFaserActsGeometry{,Lib}.so` |
| Probe构建 | `T/software_probe_v3`；运行overlay为`T/overlay_{beforeGeo,afterGeo}_v2` |
| 环境 | Athena24.0.41，EL9/GCC13，FASERNU-04 / OFLCOND-FASER-05 |
| 材料 | before为`stage02_closure/material-current-source.json`；after为`physical_prior_material_bins_v4/material-rowmajor-rotated.json` |

新kind2输入每行52列，包含source geometry ID、完整bound mean/C和target geometry ID。source/target均从真实TrackingGeometry取出；ACTS的`findSurface()`敏感面索引不足以查找被动面，因此沿原生volumes、layers、approach/boundary surfaces遍历。未用人工plane代替Disc。ROOT新增完整C、两端surface矩阵、surface type与输出bound参数，原生global position/momentum同时保存用于独立交叉核对。

两版实际`/proc/self/maps`中的三个相关库逐个检查resolved路径。新probe对旧14条输入/版的18个数组字段精确重放，材料CSV与调用记录逐字节相同，见[legacy gate](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/bound_transport_v1/legacy_gate_v3/validation.json)。三份增量patch顺序apply后与实际及编译冻结源码逐字节一致，见`T/patch_chain_validation.json`。

## 输入、选择与控制

输入来自isolated smoothing ON、n32、within-C/256、解析seed Jacobian的原生CKF归档。源研究为100044、100116、100117各首RDO entries `[0,20)`，共60个独立源事件；本次选此前方向C极端的四条路径：100044 entry8/component6、100116 entry10/component14、100117 entry1/component24与entry8/component4。

对每条相邻native state link，从source filtered mean/C传播到target真实surface。两种输入模式：

- `native_inputs`：两版各自已保存的source mean/C。
- `fixed_before_input`：两版使用相同物理mean和完整tangent Gaussian。先验证完整surface basis满足整体旋转，再原样保留local coordinates；仅解析旋转global angles及完整C，包括交叉项。不把径向local coordinates转换成平面坐标。

每段独立从归档source重新开始；不会把本次较短步长的输出接入下一段。该实验不是内部CKF传播段的精确重放，也不用于测量生产误判率。

| DSID | 原生links | 两种输入合计/配置 | 四个控制合计paired transports |
| --- | ---: | ---: | ---: |
|100044|27|54|216|
|100116|23|46|184|
|100117|52|104|416|
|合计|102|204|816|

控制为extrapolator默认10 m步长上限、scattering/energy loss同时关闭、材料ON且步长上限100 mm、材料ON且上限10 mm。`MaxStepSize`在此工具的单位为m，真实代码使用`m_maxStepSize * 1_m`。这不是声明内部CKF也采用10 m上限。

24个实际传播jobs均成功，816对全部两侧valid；失败计数为0而不是删去失败样本。完整输入C逐元素与原生读回值相同；固定输入的physical position、direction、q/p及common-chart完整C通过核对。

## 详细数据：固定完整输入的局部传播

下表为每组fixed-input links的最大差；direction-C为共同`(phi_A,latitude_A)`下的归一化2×2差异，无量纲。

| DSID | 默认最大unit-direction差 | 默认最大position差/mm | 默认最大direction-C指标 | 材料OFF最大direction-C指标 |
| --- | ---: | ---: | ---: | ---: |
|100044|3.432×10⁻¹⁶|9.237×10⁻¹⁴|2.998×10⁻⁸|2.201×10⁻¹⁵|
|100116|2.819×10⁻¹⁶|7.848×10⁻¹⁴|2.202×10⁻⁸|2.278×10⁻¹⁵|
|100117|2.968×10⁻¹⁶|8.882×10⁻¹⁴|3.349×10⁻⁸|2.511×10⁻¹⁴|

100/10 mm控制中fixed-input mean仍在相同浮点量级；材料ON的最大C差异基本不随步长改变。本批残余C差异与材料运算相关，但尚未直接测量native过程噪声matrix或每个材料操作的完整Jacobian，不能由此指定某一材料公式为最终原因。

![Full input transport controls](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/bound_transport_v1/figures_v1/results/bound_transport_fixed_inputs.png)

## 详细数据：μ−长段的步长敏感性

主例为100117 entry1/component24、forward index7→8，source geometry ID `1873497582425080064`为原生Disc(type2)，target `1945555176731443200`为Plane(type4)。路径长度约133.934→1193.919 mm，surface中心z_F约137.401→1197.3575 mm。

各自native输入的source unit-direction差为4.791×10⁻¹²，位置差7.597×10⁻¹⁰ mm；这两个输入在四种控制中保持不变。

| 控制 | native-input输出direction差 | 输出position差/mm | 输出direction-C指标 | fixed-input输出direction差 |
| --- | ---: | ---: | ---: | ---: |
|默认10 m|4.783×10⁻⁷|6.130×10⁻⁵|7.817×10⁻⁷|1.270×10⁻¹⁶|
|材料OFF、10 m|4.783×10⁻⁷|6.130×10⁻⁵|7.874×10⁻⁷|1.270×10⁻¹⁶|
|材料ON、100 mm|4.889×10⁻⁹|2.882×10⁻⁷|9.162×10⁻⁹|7.980×10⁻¹⁷|
|材料ON、10 mm|8.127×10⁻¹²|6.470×10⁻⁹|7.037×10⁻⁹|2.225×10⁻¹⁶|

长段的差异对步长有明确响应，关闭材料并未消除这项mean差。这支持检验传播积分/数值路径对微小输入差异的敏感性；固定完整输入时两版mean精确到浮点量级，不能据各自native输入的差异直接认定同物理输入的旋转不等价。

该表的默认输出方向差4.783×10⁻⁷不同于原内部CKF段的2.858×10⁻⁶；独立extrapolator与内部CKF的action/material边界、步长或调用方式不完全相同。没有用局部实验替代原CKF端点，也未证明将步长改小就能修复完整track。后续links仍从各自旧归档state重启，因此图中更后面继承的大差异不会因本段的局部控制自动消失。

![Extreme path independent links](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/bound_transport_v1/figures_v1/results/bound_transport_extreme.png)

## 被动surface位置诊断更正

此前完整前向审计错误地将所有loc0/loc1作为plane x/y。Disc实际使用r/phi，例如本例source `[101.30695,0.6203785]`应先转为`(r cosφ,r sinφ)`；单有4×4 transform无法判断local coordinates的含义。

本次保存的154个两版native surface元数据均为Plane或Disc，并将typed local-to-global位置和原生保存global position交叉核对。该catalog用于重新归约全部14,157对前向states：8,986对有已核实类型，5,171对位置为未测量；未知类型未默认为Plane。四条极端路径全部位置已覆盖，主例index7位置差从旧算法6.951×10⁻¹⁰更正为7.597×10⁻¹⁰ mm。

更正前后identifiers、flags、raw material统计完全相同；方向C归一化指标最大绝对变化6.098×10⁻¹⁴、物理方向C指标1.001×10⁻¹³、direction差指标4.342×10⁻¹⁶，q/p差逐值相同。这是重算的舍入量级，未改变原方向分叉结论。逐项验证见[位置更正比较](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/forward_covariance_v1/surface_correction_validation.json)。旧位置指标由新`analysis_v3`取代，不再引用旧`analysis_v2`的被动位置作物理结论。

## 验证、运行与产物

实现为[原生transport任务入口](../scripts/prepare_charge_bound_transport.py)、[分析](../analysis/charge_bound_transport.py)、[坐标诊断](../analysis/bound_surface_coordinates.py)、[绘图](../analysis/charge_bound_transport_figures.py)和[更正验证](../analysis/validate_forward_surface_correction.py)。16项针对性检查通过，覆盖Disc坐标、完整交叉C旋转、未知类型、真实state links、有限prior及PDF resource aliases；Python编译、文档链接和`git diff --check`另行核对。

所有批量重建、归约、绘图和汇编经HTCondor。每次提交重新比较schedd广告与队列，选择兼容且running+idle最少的bigbird20，排除已有EOS兼容记录的bigbird21。每job申请1 CPU，BLAS/OpenMP线程=1；传播申请3000 MB/4 GB scratch。24个传播jobs合计user+system CPU322.15 s，最大RSS957,392 KiB；记录见`T/resource_summary.json`。

新增surface metadata构建1582735；旧重放1582736/1582737、gate1582741；smoke1582739/1582740；余下22个传播1582742/1582743；分析1582745；typed前向更正1582747；新三页图1582748；前向更正/图验证1582751。最终均exit0。更早第一版full-C任务因敏感面索引不能查被动面而失败；引入真实被动面遍历后重新运行，失败日志保留。准备更正图集时曾遇到已占用版本目录，提交前退出，随后使用新目录；未覆盖历史产物。

实际入口示例：

```bash
python3 scripts/prepare_charge_bound_transport.py \
  --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1 \
  --phase analysis-tasks

python3 scripts/submit_momentum.py \
  --tasks /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/bound_transport_v1/analysis_tasks_v2.json \
  --library /eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/bound_transport_v1/overlay_afterGeo_v2/lib \
  --name s06_bound_transport_analysis_v2_20261114 --stage S06-bound-transport \
  --cpus 1 --memory-mb 2300 --disk-gb 2 --flavour microcentury \
  --exclude-schedd bigbird21.cern.ch
```

这些命令对应已存在的冻结目录与唯一任务名；复跑使用新版本目录和新任务名。输入及24个任务原文在`T/manifest.json`、`T/inputs_v2`、`T/tasks`及两个`*_tasks_v2.json`中。

主产物：[逐case与汇总](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/bound_transport_v1/analysis_v2/summary.json)、[816对验证](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/bound_transport_v1/analysis_v2/validation.json)、[三页transport图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/bound_transport_v1/figures_v1/results/bound_transport_report.pdf)、[更正后的前向三页图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/common_prior_ckf_chart_v1/forward_covariance_v1/figures_v2/results/forward_covariance_report.pdf)。六张PNG按原始分辨率目视检查，逐图/汇编PDF的normalized graphics与page boxes一致。

[136页完整S06报告](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v18/charge_stage06_report.pdf)由cluster1582754生成；继承report_v16的130页，索引、decoded vector streams及page boxes全部精确相同，另列三页更正前向图和三页transport图。该报告更新report_v17，保留原Alma参考分布和全量生产图；验证见[报告继承检查](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v18/lineage_validation.json)。

## 阶段性结论

事实：实际native surfaces与完整输入下的局部传播mean闭合；原生微小输入差异在长段传播中存在明显步长敏感性；fixed-input方向C的残余与材料操作相关。被动面位置错误已透明更正，未知类型不再产生位置数字。

尚未测量：内部CKF同一执行链的100/10 mm步长控制、native transport J及材料前后过程项、完整有限先验的likelihood加权posterior与分量数/宽度收敛。本次不能据mean闭合、sign一致或零失败宣布S06 Q5完成。

## 下一步计划

1. 在真实CKF内部设置明确的步长上限，保持完整物理seed、hits、regular chart及isolated smoothing ON；默认关闭先做全输出精确重放，再比较100/10 mm对极端mean/C、候选选择、refit与sign的影响。
2. 记录实际native材料covariance操作和transport Jacobian，区分material ON的约10⁻⁸ residual来源；不将`P−J C Jᵀ`余额直接当native Q。
3. 将已验证的坐标与数值控制用于有限Gaussian mixture的完整物理分布、likelihood加权posterior和独立回归，保留失败及未选择分量的完整分母。
