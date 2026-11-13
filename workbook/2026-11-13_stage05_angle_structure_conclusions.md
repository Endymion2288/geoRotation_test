# S05：角度结构来源与验收结论

S05角度结构定位已完成：主要来源由生产样本和独立 full-G4 控制验证，290页统一图集及A1–A5交付核对通过。beforeGeo额外近零结构来自near-pole方向参数化下的refit update/smoothing；联合regular chart控制在RDO与独立μ±中移除主要增强并改善角误差核心和大尾。afterGeo主要结构具有strip控制及预留样本response证据。完整实验记录见[阶段报告](2026-11-12_stage05_slope_structure.md)，要求见[角度路线图](../roadmap/03_slope_structure.md)及[验收审计](2026-11-12_stage05_acceptance.md)。

## 阶段目标

解释整体坐标旋转前后的角度结构，确定首次引入异常的实际 Calypso 环节，用固定输入的最小控制和独立模拟验证来源，并量化剩余结构。主图对应 `problem.pdf` p4、`geometryRotation.pdf` p19–20、p23；同时报告效率、动量、电荷及外推失败。

## 输入、物理定义与研究范围

结果根目录为 `/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes`，下文记为 `S`。

主样本为 100044 μ+ 的完整 500,000 events/版本，RDO 路径为：

```text
/eos/experiment/faser/data0/sim/mc24/particle_gun/100044/rdo/s0013/FaserMC-MC24_PG_mupl_fasernu_5mrad_flukaE-100044-00000-00099-s0013-RDO.root
```

entries `[0,250000)` 用于 discovery，`[250000,500000)` 用于 validation；峰窗口及响应分组在 validation 前冻结。辅助 100116/100117 各10k，110033 全5文件10k、110051 全20文件20k。所有跟踪、传播和模拟来自实际 Calypso、完整 geometry/field/material 及 full G4/digitization。

beforeGeo 为 master `1eb885ae22299fb39bb352bead6005d41b0f953e`，afterGeo 为旋转分支 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`，同用 Athena24.0.41、ACTS32.0.2、FASERNU-04/OFLCOND-FASER-05。两个原生分支的 seed covariance 不同，因此版本总差值不作为相同物理先验下的纯旋转因果效应。beforeGeo 材料为 `stage02_closure/material-current-source.json`，afterGeo 为用户指定的 AFS `material-maps-alma9-rotated.json`。

统一物理量为 `sy=py_F/pz_F=py_A/px_A`。st1 参考面为 `z_F=0 mm`；VetoNu 为 `z_F=-3112 mm`，位于两块板之间。RDO sensor truth 为 SiHit chord 方向；fixed-plane truth model 将 production p 大小附在首个 st1 SiHit chord 后调用 Calypso 传播，不冒充局部 G4 momentum。新 near-pole gun 从 st1 生成，st1 truth 直接采用生成方向。

## 工作内容与主要来源

| 观察结构 | 首次出现位置及直接证据 | 结论的适用范围 |
| --- | --- | --- |
| 历史角度及 pz 选样失真 | 实际 rotated Ntuple 存 A-frame momentum，旧 `atan(py/pz)` 的分母不是束流方向；物理公式及 writer 开关逐字段闭合 | 输出解释已定位；正确主图使用同一物理面和 physical sy |
| 旋转前 beam 附近额外窄峰 | 固定 final 候选的 native refit F→S 增强；periodic φ、canonical angles、空间 covariance basis 和 input plane 均不足以解释；regular 方向 chart 的前后向联合控制显著减少增强 | 来源定位到 near-pole 参数化下的 refit update/smoothing；须同步检查残差质量，不能单独以消峰宣称修复 |
| segment 周期及 afterGeo st1 两肩 | 实际 pitch=0.08mm、lever arm=63mm、stereo≈±0.02rad；相同完整 hits 的独立40k locX 控制改变最终谐波与冻结肩部计数 | measured strip coordinate 有可重复贡献；不是把两肩全部面积归因于单一离散机制 |
| VetoNu 的 mrad 宽分布 | 实际靶材料约250.31X0；full-G4 μ± 100/250GeV 的宽度符合约1/p标度 | 靶散射是主要宽度来源；该宽度不是 tracker 本征角分辨率 |
| A′ 双track 的窄中心 | 全110033/110051 按 decay z、opening angle、st1 truth 间距分组；小间距组 truth 本身较窄 | 输入拓扑及接受度有贡献，不将其用作单μ样本 refit 增强的原因 |
| CKF 逐层 smoothed states 异常 | 候选隔离后 smoother 恒等式恢复，所测 final 方向不变 | 解释逐 state 的共享 slot 污染，不解释所测 final 近零峰 |

## near-pole 数值原因与最小控制

beforeGeo 束流接近 native polar θ=0。微小物理方向变化可以对应很大的 φ 变化；在这个 chart 中对均值、协方差和方向差分作线性 Gaussian 运算，不能由简单的 φ 周期差分或 θ 反射规范化保证稳定。

实际 Calypso worktree `geoRotation_calypso_slope_chart_update` 的[完整补丁0017](../patches/0017_refit_regular_update_chart.patch)将 update 与 smoothing 中的方向改写为远离该极点的 longitude/latitude。变换使用 `(x_A,y_A,z_A)=(z_F,y_F,-x_F)`，latitude 以0为中心；local measurement、geometry、field、材料及原生 transport 保持相同。测量更新调用实际 `Acts::GainMatrixUpdater`；RTS 调用实际 `Acts::GainMatrixSmoother::calculate`，同时变换完整6×6 covariance及 transport Jacobian。两个开关独立、默认关闭，构成 native / smoothing-only / update-only / both 四个控制。

1887组实际 native writer 角度的 header 检查给出方向往返最大误差6.94×10⁻¹⁸、逆 Jacobian 误差2.18×10⁻¹¹、中心差分 scaled误差3.60×10⁻¹¹。实际53个 translation units 的编译1574002成功。四模式20-event smoke1574003验证全部 native 保存字段复现基线，smoothing-only复现前一版本，实际加载的三份 Calypso 库及 refit outcome 身份均闭合。

10k RDO 为 entries `[0,5000)` 与 `[250000,255000)`；下表统计后半段 validation，先固定原生候选及原生 p=50–250GeV，要求 primary μ、purity≥0.7及三站 long。近零窗口固定为 `|sy|<100µrad`，残差参考为同面 truth model。

| 指标 | native | smoothing-only | update-only | both |
| --- | ---: | ---: | ---: | ---: |
| 固定 native-p cohort 的近零 /1758 |109|43|84|61|
| 同 cohort truth 近零 /1758 |58|58|58|58|
| 固定 truth-p cohort 的近零 /1595 |56|32|52|31|
| 同 truth-p cohort truth 近零 /1595 |30|30|30|30|
| st1 residual RMS /mrad |0.46667|1.49993|1.87396|0.32907|
| 绝对角误差68%分位 /µrad |108.71|100.06|100.67|41.76|
| `|residual|>1mrad` /1758 |38|122|126|8|
| `|residual|>3mrad` /1758 |5|26|31|2|

仅改一侧会增大残差尾部，联合修改同时改善近零结构、核心和大尾。both−truth 的近零比例差为 **+0.1706 [−0.4077,+0.7928]个百分点**；truth-p 为 **+0.0627 [−0.4392,+0.5608]个百分点**。区间采用1000次 source-event Poisson bootstrap，对应 `chart_update_analysis_v1`；质量图中的 empirical bootstrap 区间另有定义，不能混用。这些是所列 cohort 和窗口的局部余量，不外推为全500k或任意峰宽的上限。

全部9707条候选身份及majority不变，38998条 seed/CKF states逐值不变。update-only/both 的完整 final hits 有76条变化；在相同完整hits的1744条 native-p 子组中，近零仍为108→61，truth58，residual RMS为0.46156→0.28220mrad。变化因此不能由丢弃难候选或更换 truth 解释。actual refit success/fallback 为9110/597，联合控制前后相同；固定1758条中为1753/5，原109条和控制后61条近零均 refit 成功。

![Production near-pole control](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/chart_update_analysis_v1/regular_both/production_muplus_st1_regular_both_refit.png)

独立 μ± gun 使用 seeds50981/50982，各4000事件、100GeV、st1生成点、x/y均匀±30mm、theta_x/y均匀±0.5mrad。所有模式共用同一次 full G4/digitization 输入，实际 geometry 与 tracking library maps 均已验证；新环境中的 native/smoothing-only 全部保存字段精确复现旧有效 beforeGeo 控制。st1参考为精确生成方向，固定native候选和native-p窗口的独立结果为：

| 独立样本 | 分母 | native近零 | smoothing-only | update-only | both | 同面truth |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| μ+ |3906|1299|686|1610|765|783|
| μ− |3912|1321|660|1647|753|749|

以下区间统一采用refit报告中的1000次source-event Poisson bootstrap。both−truth 的剩余近零差为μ+ **−0.4608 [−1.1779,+0.2753]个百分点**、μ− **+0.1022 [−0.6480,+0.8544]个百分点**；两者均与0相容。same-hit子组分别3670/3713条，native→both近零为1229→722、1268→722，对应truth730/713；配对剩余差区间同样包含0。全部majority不变。独立两电荷因此重复验证了方向chart效应，不依赖RDO的传播truth model或改变hits的子群。

独立μ+ RMS为0.29319→0.10404mrad，绝对误差68%分位191.14→31.27µrad，超过1mrad的误差47→7；μ−分别为0.28731→0.09288mrad、185.14→30.88µrad、40→13。两电荷均3965/4000 refit成功、35 fallback，联合控制不改变该计数，全部近零候选refit成功。旧smoothing-only在独立样本中形成的近零欠量，在联合控制后与0相容；update-only则使近零堆积更强，说明只切换单个环节并不能恢复一致的方向估计。

![Independent muon control](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/chart_gun_analysis_v2/regular_both/muplus_st1_regular_both_refit.png)

## 效率、动量、电荷及误差范围

10k配对控制的效率口径为 primary μ、purity≥0.7、三站 long，truth 去重；动量左尾定义为 `p_reco<0.5*p_production`，与全量 fiducial+H18/Good95效率的分母不同。

| 指标 | native | smoothing-only | update-only / both |
| --- | ---: | ---: | ---: |
| 效率 |8680/10000|8680/10000|8674/10000|
| wrong charge |317/8680|317/8680|55/8674|
| momentum left tail |773/8680|773/8680|135/8674|

这些是单个10k控制的伴随指标，不是完整 production 修复验收。固定 native-p 的 angular pull RMS在联合控制后仍为2.03，truth-p为1.90；误差未校准为标准normal。完整500k原生版本比较的 fiducial H18/Good95效率为84.8189%→83.8805%，charge misID为4.29765%→0.816043%，动量左尾为9.47544%→1.82637%；各自source-event区间及分母在主报告。版本比较中的效率下降和原生先验差异均需保留。

独立gun中，联合控制的pure-long效率从μ+3920/4000降至3840/4000、μ−3935/4000降至3858/4000，即下降2.000/1.925个百分点；这是实际副作用。wrong charge分别从1/3920到0/3840、1/3935到1/3858；动量左尾从13/3920到0/3840、20/3935到1/3858，三个参考面的外推失败均为0。全部4000候选仍保留，角度主比较固定原native cohort，same-hit子组也改善，因而没有用效率下降筛出消峰结果。独立pull RMS为1.373/1.331；该补丁用于本阶段成因验证，仍需生产效率与误差校准研究。

## afterGeo 剩余结构及可约束范围

独立40k同粒子、同完整hits、固定 native-p 的11511对中，locX truth-coordinate控制使几何谐波从0.027504降至0.006630，配对变化−0.020874 [−0.026757,−0.004612]。−0.71/+0.59mrad冻结肩窗的计数分别358→329、260→234，占比变化为−0.2519 [−0.4256,−0.0653]及−0.2259 [−0.3929,−0.0602]个百分点。该控制保留真实 digitization、噪声与测量误差，不能称为 digitization off。

五种预定义 response 分组均由 discovery 学习并用于 validation；共同 covered 范围和全部 uncovered reco分开保存，24组 response flow 重放闭合。加入分组系统项和未覆盖贡献后，两肩冻结窗的保守额外正过量余量分别为全 validation选样的 **0.4504% / 0.6048%**。它支持谱、接受度、测量response混合能重现主要结构；不是对任意条件、任意窄峰的普适排除，也不把 Gaussian mixture拟合本身当原因证明。

## 运行、验证与交付

本轮 build1574002、smoke1574003、RDO production1574005和analysis1574007均成功。独立gun1574004的μ+成功；μ−首次在事件处理前因EOS共享父目录并发mkdir失败，保留原错误后仅重提失败任务1574006，已成功。三模式独立gun归约1574010、RDO/gun质量归约1574008/1574011均成功。提交端预建共享父目录，各worker仍创建自己的全新输出目录。每次提交重新比较兼容schedd及队列，选择最低running+idle的bigbird20，排除已记录EOS兼容问题的bigbird21；具体资源与选择证据在各submission.json及Condor日志。

`S/loader_audit_v1/validation.json`给出旧nested driver混用master geometry与rotated tracking的20-event实证：相同生成、SiHit、segments及seed covariance下，清理跨版本库/Python路径后0→20 tracks。修正后的完整gun每电荷afterGeo均为4000 tracks；旧afterGeo空参考不作为物理证据。主500k两版分别以本版overlay启动，不属于这种nested fallback。

本轮RDO六张refit图的108个histograms、独立gun十八张refit图的324个histograms、RDO/gun十八张质量图的498个histograms已分别通过counts/sumw2/flow及原始数组重放；regular gun的旧beforeGeo数值也全部精确复现。质量图显示完整±3mrad残差、对数计数和图外entries；目视发现初版图外计数与图例重叠后，仅从冻结counts重导展示，不重复bootstrap或tracking。40项仓库测试通过；0013–0017与实际worktree对应关系验证见 `S/patch_validation_v2/validation.json`，其中0013+0014共同组成basis worktree的完整diff。

最终质量图 `S/chart_diagnostics_display_v1` 的六报告、498个histograms、18图全部由原始Calypso数组独立重放；494991个原数值字段精确复现。新独立gun的st1质量图和主对照已目视，均显示完整残差范围及图外计数。对应脚本为 [slope_chart_diagnostics.py](../analysis/slope_chart_diagnostics.py) 与 [validate_slope_chart_diagnostics.py](../scripts/validate_slope_chart_diagnostics.py)。

统一交付在 `S/final_report_v1`，冻结catalog包含52组报告、75项证据，最终PDF为290页、4.18MB。汇编1574013 exit0，wall61.59s、CPU24.26s、最大RSS81076KiB，申请1CPU/3500MB/3GB scratch。工具使用已冻结pypdf5.9.0，先在worker scratch复制输入并关闭PDF，再复制至EOS；PNG逐张完整解码。实际任务与选择schedd的证据在 `tasks.json` 和 `s05_final_report_v1_worker/submission.json`。

- [290页完整图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/final_report_v1/results/stage05_analysis_report.pdf)：第1–57页为100044参考主图；第264–266页为生产RDO联合chart控制；第267–272页为独立两电荷联合控制；第273–290页为三模式RDO/gun残差与pull质量。
- [图组页码与来源索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/final_report_v1/results/index.md)及[逐图PNG/PDF、参考与证据索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/final_report_v1/results/index.json)。bin counts/sumw2/flows保存在索引所列各组summary及验证记录中。

汇编与源报告的290页decoded stream、text、page boxes、rotation，以及52组bookmarks/index均由 [validate_slope_final_artifacts.py](../scripts/validate_slope_final_artifacts.py)独立核对通过；290张PNG完整解码，290个单页PDF与源报告的绘图指令及实际resource定义一致。首轮1574014检查发现Matplotlib单页与PdfPages的字体、透明度resource编号不同；具体失败页的5369条绘图指令在解析实际resource定义后完全一致。1574015保留汇编对源报告的逐字节检查，单页对原报告另核对resource语义、绘图指令、文字与page boxes，exit0；wall246.89s、CPU193.54s、最大RSS100800KiB，图表及统计未变动。完整结果及逐图汇编页码在 [validation.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_slopes/final_report_v1/validation.json)。已渲染并目视汇编第4页主图及第267页独立μ+控制，另已目视最终μ+质量图和μ−主控制图。

复现入口为 `S/final_report_v1/tasks.json` 及 `validation_tasks_v3.json`，分别冻结汇编和最终验证的完整命令、代码与catalog；新运行使用新的output/job目录。实际提交入口均为 `scripts/submit_momentum.py`，`--stage S05 --cpus 1 --memory-mb 3500 --disk-gb 3 --flavour microcentury --exclude-schedd bigbird21.cern.ch`，相应name为 `s05_final_report_v1`、`s05_final_report_validation_v2`。40项已有单测通过，最终18个Python文件编译、本地文档链接及 `git diff --check` 已检查；阶段报告与相关代码/配置同次Conventional Commit交付。

## 阶段性结论与下一步

事实：输出frame错误、measured strip coordinate周期、靶散射及输入拓扑均有实际Calypso证据。beforeGeo的额外近零结构由refit near-pole方向参数化引入；联合regular update/smoothing在RDO与独立两电荷gun、固定原生选样及相同完整hits子组中均移除主要增强，所测剩余近零差与0相容，同时改善角误差核心和尾部。

来源定位、所列窗口余量及A1–A5交付验收已完成。后续以新保留输入验收正式生产补丁的全样本指标、效率副作用处理与pull误差校准；这三项尚未测量。
