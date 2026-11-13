# S06：isolated smoothing 下的全量电荷误判与路径统计

**原生3-station forward配置、两电荷各500,000个相同RDO源事件、两版均开启isolated smoothing：μ+误判率从0.163695%降至0.010923%，μ−从0.314528%降至0.012519%。旋转后μ−−μ+差为+0.001596个百分点，95%区间[−0.004386,+0.007578]个百分点。** 在此精度下未见Alma历史量级的正负差异；这些是整体原生版本比较，seed covariance也不同。

全量路径记录中，afterGeo最终59例错误的最后一次“正确→错误”发生在seed的3例、CKF measurement update的56例；conversion与refit没有新生的最终错误。后续[谱、模式与逐运算控制](2026-11-13_stage06_charge_controls.md)已完成production-only平衡、100k方向/IFT/merge、119事件精确重放及全500k matched seed covariance：56次翻转中31次发生于无SDO关联hit、25次于仅primary支持hit。

## 阶段目标

回答 `problem.pdf` p4 的μ±电荷误判问题，执行[Q1–Q5路线图](../roadmap/04_charge_misidentification.md)：恢复历史口径、测量整体旋转前后及正负差异、分离输入谱和接受度、定位符号变化，随后用物理先验与对称性控制完成归因。本报告记录已完成的主测量及全量compact lineage，不将尚在执行的机制/对称性实验计作完成。

## 软件、样本与执行

结果根目录 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`。

| 样本 | 电荷 / PDG | RDO目录 | 输入 |
| --- | --- | --- | --- |
|100116|μ+ / −13|`/eos/experiment/faser/data0/sim/mc24/particle_gun/100116/rdo/s0012/`|全部5文件，各100,000事件|
|100117|μ− / +13|`/eos/experiment/faser/data0/sim/mc24/particle_gun/100117/rdo/s0012/`|全部5文件，各100,000事件|

每电荷before/after各500k，共2百万次事件重建。事件身份为冻结的输入文件及entry；生产run/event保留但不用于去重。discovery为文件0/1/2，held-out为3/4。

before基底master `1eb885ae22299fb39bb352bead6005d41b0f953e`，after基底旋转版`88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。源码worktrees分别为 `/eos/home-x/xcheng/FASER/geoRotation_calypso_charge_before`、`.../geoRotation_calypso_charge_after`；冻结库为 `C/software_beforeGeo`、`C/software_afterGeo`，与只读observer组成`C/overlay_*`。Athena24.0.41、ACTS32.0.2、FASERNU-04、OFLCOND-FASER-05。before材料为source/current map，after为用户指定AFS rotated JSON。所有重建、refit、传播均调用实际Calypso。

[patch0021](../patches/0021_charge_diagnostics_after.patch)、[patch0022](../patches/0022_charge_diagnostics_before.patch)分别相对两个原生基底，开启默认isolated smoothing并添加只读seed/candidate/update/conversion/refit记录。保留原生有问题的`trackGroup`值，真实ancestor seed用独立索引记录，避免观察改动影响原重建。两个补丁不能串联应用。[patch0023](../patches/0023_charge_observer_qop.patch)只给已有observer增加q/p与variance字段；ACTS单位GeV⁻¹，Trk q/p×1000、variance×10⁶。

6组20-event观察ON/OFF对照的58数组全部精确一致，after还复现干净基线，实际加载四份库逐项核对：`C/charge_smoke/validation.json`。S05完整500k的ON/OFF最终量和membership验证见[专项报告](2026-11-13_stage05_isolated_validation.md)。

主after cluster1574103、before1574104各200任务；before μ− shard0066首次在事件处理前读取CVMFS GaudiConfig2数据库时发生`Errno 5`。失败输出归档后，以原冻结代码和range重提1574106，exit0；没有把失败尝试加入分母。每个成功shard5000事件，400个最终shards覆盖全部输入。主归约1574109成功。提交端逐次比较schedd负载，选择兼容且running+idle最少的bigbird20，排除已验证EOS不兼容的bigbird21。

## 匹配、cuts与分母

使用primary μ：正确PDG、status1、`0<barcode<200000`；逐source entry恰好一个。truth fiducial要求生产vertex及ST1/2/3平均SiHit位置的F-frame横向半径均<100mm。主分母不强制H18。

SDO majority matching与重建charge/p无关；多个匹配候选先按majority hits数、总hits数、chi2、ordinal排序，**再**应用质量选择，不以次佳通过者替代最佳失败者。主样本要求final purity>0.7、三站long、ndof>0、chi2/ndof<25、layers≥7、ndof+5≥12、rmax<95mm。r<100、H18、exactly-one、无purity和held-out均独立报告。

非法/零κ、unmatched、duplicates、低purity、quality failed分别计数；signed κ与final charge逐条闭合。两个分母分别为：

```text
conditional_misID = selected wrong / selected valid matched
truth correct/wrong efficiency = fiducial truth with selected correct/wrong / fiducial truth
```

每个源事件一个primary，未把逐层states作为独立样本。率给精确binomial68/95%区间；旋转差使用同truth的ratio差delta method，保留配对covariance。全部选择、表与直方图在`configs/charge_study.json`、`analysis/charge_tables.py`和`analysis/charge_analysis.py`实现。

## 主测量

| 电荷 / 配置 | wrong / selected | conditional misID / % | 95%区间 / % |
| --- | ---: | ---: | --- |
|μ+ before|412 / 251687|0.163695|[0.148280,0.180276]|
|μ+ after|27 / 247182|0.010923|[0.007199,0.015892]|
|μ− before|817 / 259754|0.314528|[0.293358,0.336819]|
|μ− after|32 / 255612|0.012519|[0.008563,0.017673]|

after−before的配对差为μ+ **−0.152772 [−0.168713,−0.136832]个百分点**、μ− **−0.302009 [−0.323686,−0.280333]个百分点**。before正负差为+0.150833 [0.124128,0.177538]个百分点；after差与0相容，不能表述为证明两电荷完全相同。

| 电荷 / 配置 | fiducial truth | correct / fiducial / % | wrong / fiducial / % |
| --- | ---: | ---: | ---: |
|μ+ before|293601|85.583837|0.140326|
|μ+ after|293601|84.180572|0.009196|
|μ− before|303839|85.221779|0.268892|
|μ− after|303839|84.116917|0.010532|

误判率下降伴随约1.1–1.4个百分点correct efficiency下降，二者均保留。r<100的conditional结果为μ+513/272605→40/268009、μ−966/281510→50/277270；这些不与主r95分母混用。exactly-one在r95下仅μ−after少1例wrong：31/255582；其余wrong数与inclusive相同。

完整源覆盖、7×7状态迁移、truth-p分箱、无purity/H18/duplicates及失败分类见[主summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/main_analysis_v1/summary.json)。同RDO的production momentum、position、sensor truth direction与station truth position在before/after逐元素相同。

![Signed q/p, both charges and native geometry versions](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/main_analysis_v1/figures/signed_qop.png)

![Wrong-sign q/p distributions](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/main_analysis_v1/figures/wrong_signed_qop.png)

## 入口条件谱平衡

最初joint-bin重加权在held-out的sy SMD仍为−0.2328，未达到谱平衡验收，不能据此归因；该探索性图不纳入正式主结论。

后续truth-only logistic overlap模型仅用discovery文件训练，输入为log production p、st1 truth x/y、sensor truth sx/sy。包含标准化linear、squares、全部pairwise项、20/50/80% hinges，L2=10⁻⁶；等charge先验。μ+权重1−e(x)、μ−权重e(x)，共同支持`0.01≤e≤0.99`。模型优化成功596 iterations，参数在`C/overlap_v1/model.json`冻结。

预留文件的五项`|SMD|max=0.019243`、weighted marginal `KSmax=0.042688`，均通过预先记录的0.05/0.05阈值。truth有效N为μ+71195.9、μ−77026.6；selected有效N为after μ+59983.7、μ−64593.6，sumw/sumw2和weight tails随结果保存。

入口条件平衡后的held-out μ−−μ+差：before **+0.021769 [−0.053414,+0.096951]个百分点**，after **+0.002284 [−0.009223,+0.013791]个百分点**。全文件after差为−0.004029 [−0.012045,+0.003988]个百分点。区间条件于冻结的discovery模型，未包含训练模型不确定度。

这是**入口条件**对照：st1位置与sensor方向已包含真实上游场弯曲和接受度，不能全部归为generator输入谱。另行完成的production p/x/y/sx/sy、全部primary母样本控制在held-out达到最大|SMD|=0.005970、KS=0.006875；before正负差为−0.026085 [−0.067004,+0.014833]个百分点，after为+0.001575 [−0.007485,+0.010635]个百分点，详见[控制报告](2026-11-13_stage06_charge_controls.md)。结果入口：[overlap summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/overlap_v1/summary.json)。方法参考Li/Morgan/Zaslavsky，[DOI:10.1080/01621459.2016.1260466](https://doi.org/10.1080/01621459.2016.1260466)；采用overlap权重及报告有效N，不将加权曲线平滑本身作为物理解释。

## 全量lineage：错误符号在哪一步形成

全400shards由实际Calypso的`charge_lineage.npz`归约，不重跑tracking。after1574114、before1574115各200任务成功，合并1574278成功。全部selected/wrong数与主summary精确闭合，final q/p与observer逐条相同；conversion最大相对差<6.4×10⁻¹⁶且signed κ保持，排除本测量中的conversion sign错误。

**下表为最终wrong路径最后一次进入错误符号的阶段，不等于第一次曾出现错误。** 中途错误可被纠正，之后再翻转，两个标签分别保存。

| 样本 | seed保留至最终 | CKF measurement update | refit | conversion / 新native reference错误 |
| --- | ---: | ---: | ---: | ---: |
|μ+ before|3|125|284|0|
|μ+ after|1|26|0|0|
|μ− before|9|258|550|0|
|μ− after|2|30|0|0|

after首次记录错误分别为μ+ seed3/update24、μ− seed5/update27。最终59条中56条最后在CKF update进入错误、3条从seed保留。before有834/1229条最后在refit进入错误，关联S05已证明的near-pole refit问题仍需在本μ±样本中以详细trace和最小控制测量，不能仅据stage标签替代机制证据。

隔离后的逐state S重算不替换finder已缓存的native reference mean；conversion使用后者。分类按seed→forward P/F→native reference→conversion→refit/final，另存isolated S，避免把未参与后续conversion的S当成因果输入。

after两电荷各5个wrong事件存在primary-supported、正确符号的native候选，且这些候选被pruning。逐事件核查已确认这些备选全部没有majority tie，淘汰为nMeasurements≤12和/或shared hits>6，完整候选见[pruning表](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/replay_analysis_v1/pruning.json)。before对应正确候选存在数为306/573，其中被pruned42/70。此为机制线索，不把“存在另一候选”直接称为可无损修复。

主selected refit：μ+before251687、after247182全部成功；μ−before259753成功/1fallback，after255612全部成功。after错判中μ+25/27、μ−25/32的`|κ|/σ_κ>3`，其中23/21超过5；剩余错误不能统称为高动量下统计分辨不出符号。

证据：[全量lineage summary](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/lineage_summary_v1/summary.json)及同目录`finals/candidates`紧凑表、`replay_events.json`。详细audit的20-event验证两电荷58数组逐值相同，update恒等式最大绝对误差<9.2×10⁻¹⁵ GeV⁻¹；全部59wrong与60normal的冻结重放输入见`C/charge_replay_v1/targets.json`。

## Alma历史复现与具体差异

Alma p2对应 **old=r0019、new=dev**。v1反向配对已标SUPERSEDED；正式使用`C/history_analysis_v2`，1574105成功。

当前可执行25/7/12质量阈值、r<100、truth fiducial及exactly-one Good下：

| 样本 | dev new | r0019 old |
| --- | ---: | ---: |
|μ+|424/217967=0.194525%，400k源|263/269683=0.097522%，500k源|
|μ−|970/281265=0.344871%，500k源|564/278490=0.202521%，500k源|

μ+ dev缺file1，不能以两个全量分母强行配对。共同400k文件的old-correct/new-wrong为266/217967=0.122037%，μ−500k为611/281265=0.217233%；分别接近Alma的0.122%/0.223%。

历史完整cuts尚未全部恢复：25/7/12来自当前可执行源码，Alma明确的是r<100和exactly-one；旧`t_truthHitRatio`分母为ndof+5，当前purity用实际hits；primary truth key受MinMomentum门槛影响，未知fiducial不能当false；Alma p3的PDG/charge matching含义仍需源码证实。因此这项为具体定义下的复现和差异清单，不宣称精确重现所有历史数字。原生before/after与历史dev/r0019是两条比较轴。

## 验证与产物

全部48项现有单测通过，包含无charge偏置的ranking、empty truth、weighted有效N、配对ratio covariance和错误路径分类；日志`C/verification_v2/unit_tests.log`。主图按Alma p16–27的q/p、wrong-sign、truth/reco p、2D相关及truth-p分箱率组织，英文PNG、单页vector PDF、bin counts与索引随图保存。原15页v1图集的末尾三页joint-bin平衡图属于未通过平衡验收的探索控制；正式图集使用accepted overlap结果。

## 阶段性结论与后续验收

事实：在本次完整样本和明确cuts下，两版ON的主误判率显著降低，after正负差异在上述区间内与0相容；原生before有refit新生错误，after剩余主要在CKF update。入口条件谱差异显著，accepted overlap后正负差异区间包含0。

已完成的production-only、100k/电荷/模式、119事件r/S/K/C及500k matched-prior实测见[控制报告](2026-11-13_stage06_charge_controls.md)。native组合CKF hit过滤、q/B全重放及100k before refit chart电荷控制均已完成；独立full-G4 ON测量与模式fallback精确计数见[独立控制报告](2026-11-13_stage06_independent_charge_controls.md)。[测量/seed与联合角度回归](2026-11-13_stage06_measurement_seed_controls.md)进一步显示：59错号中centroid恢复23条selected正确，production seed q/p-only无selected恢复；六条pure wrong的centroid结果为5正确、1质量失败，固定完整seed后的测量归因继续执行。before updater-only虽修电荷却增加角大尾，both才联合改善。prior仅匹配一阶完整covariance，不证明有限宽度非线性physical prior或全链纯旋转闭合；Q5后续回归保持未完成。

下一步以这些冻结控制完成Q4/Q5：任何“算法不对称”“纯旋转改善”或“已修复”的结论须有相应独立证据；本报告不缩小Q1–Q5的原目标。
