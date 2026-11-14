# S05/S06：材料 rotation 编码修正的150万事件全量重建对照

**仅修正 AFS 材料图12项 rotation 的序列化编码后，μ+仍为原27个错号事件、μ−仍为原32个错号事件，没有 selected correct↔wrong 迁移。100044错号292→291，包含7次correct→wrong、7次wrong→correct及1次wrong→quality_failed。冻结角峰窗口与角误差核心变化很小。材料编码错误需要修正，但不是这批μ±剩余错号及100044肩峰的主要来源。**

## 阶段目标

在[原生材料decoder诊断](2026-11-14_stage06_physical_prior_material_decoder.md)之后，测量真实CKF/refit对编码修正的响应；固定输入、实际库、runner、seed配置和选样，同时检查角度、动量、电荷、效率、候选迁移及truth-model变化。所有tracking使用真实geoRotation Calypso，candidate-isolated smoothing始终开启。

根目录 `C=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1`；本轮目录 `M=C/material_encoding_control_v1`。

## 工作内容与输入配置

原文件为 `/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json`。冻结副本为`M/legacy_afs_map.json`；新文件`M/corrected_afs_map.json`只转置原12个rotation数组。逐项验证其余所有JSON值完全相同，包括geometry IDs、translations、bin axes/data、surface存在性及metadata；原路径未修改。编码的runtime根因和真实材料查询闭合见上述decoder报告。本轮是afterGeo内部的材料编码单因素，不是整体beforeGeo/afterGeo。

| 样本 | 完整输入与事件 | 原始结果 | 实际库入口 |
| --- | --- | --- | --- |
| 100116 μ+ | 全5文件，各100k，共500k | `C/main_100116_v1` | `C/overlay_afterGeo/lib` |
| 100117 μ− | 全5文件，各100k，共500k | `C/main_100117_v1` | 同上 |
| 100044 μ+ | 原单文件的全部500k | `stage05_slopes/isolated_100044_v1` | `C/software_after_clean/lib` |

每组100个5k-event ranges，共300个完整重建任务；冻结路径、entry范围及实际加载的4个library由`M/manifest.json`和`M/full_validation.json`核对。重用原始runner和worker；未改变measurement、CKF/refit、质量条件或covariance。每组legacy/corrected各20-event smoke共6任务：legacy与原完整shard前20事件58/58/56个数组字段逐值相同，corrected的原truth/segments与事件身份相同。

charge主选择沿用[S06基线](2026-11-13_stage06_charge_baseline.md)：与reco charge/p无关的matching、purity>0.7及Good95。角度图保持S05原选择，并在共同物理面比较`p_y,F/p_z,F`：st1 `z_F=0 mm`、VetoNu `z_F=−3112 mm`；p窗口为50–250 GeV。动量参考是HepMC production p。角度residual主对照固定legacy production-p/首ST1-SiHit chord的原生Calypso truth model；actual local G4 momentum不在本轮测量中。

100044的旧S05 observer未单存q/p，采用同一个upstream参数的`charge/|p_cartesian|`。μ±有原生q/p，独立核对此表达式的最大相对差约2.22×10⁻¹⁶；没有为100044虚构q/p variance。第一次归约v1因此缺列而失败，修正后新冻结analysis_v2，smoke及full均成功；原失败目录保留。

## 详细数据：电荷、效率和动量

| DSID | wrong/selected：legacy→corrected | misID：legacy→corrected / % | p_reco<0.5 p_prod计数 | correct/fiducial：legacy→corrected | 68%绝对角误差：legacy→corrected / μrad |
| --- | --- | --- | --- | --- | --- |
| 100116 | 27/247182→27/247181 | 0.01092313→0.01092317 | 964→961 | 247155/293601→247154/293601 | 47.832011→47.844861 |
| 100117 | 32/255612→32/255610 | 0.01251897→0.01251907 | 1145→1145 | 255580/303839→255578/303839 | 43.931735→43.949190 |
| 100044 | 292/352812→291/352813 | 0.08276362→0.08247995 | 2211→2211 | 352520/421268→352522/421268 | 27.806062→27.800175 |

误判率分母为selected valid matched truth；correct efficiency分母为fiducial truth。区间采用source-event cluster delta method。correct efficiency的配对变化分别为：

- 100116：−0.0003406个百分点，95% `[−0.0027475,+0.0020663]`。
- 100117：−0.0006582个百分点，95% `[−0.0028928,+0.0015763]`。
- 100044：+0.0004748个百分点，95% `[−0.0018976,+0.0028471]`。

修正后的μ−−μ+原始misID差为+0.0015959个百分点，独立样本95%区间`[−0.0043862,+0.0075780]`。沿用原冻结的production-only overlap weights后，μ+和μ−的weighted misID分别为0.0123154%和0.0108352%；差为−0.0014802个百分点，条件于冻结权重模型的95%区间`[−0.0075543,+0.0045938]`。held-out结果及分箱原始分子分母见`M/figures_v1/results/charge_differences.json`及full summary。

μ±同一59个wrong仍全部wrong；其q/p数值均发生变化，不能据此声称原内部trace逐值不变。两样本的selected correct↔wrong迁移均为0/500k源事件；独立源事件假设下，事件级单侧95%概率上限均为`1−0.05^(1/500000)=5.99145×10⁻⁶/event`，见`M/sign_migration_upper_limits.json`。此上限只针对所定义的selected sign迁移，quality边界变化另计。分子完全相同使配对misID比值差的delta区间很窄，该区间不用于断言总体分布完全不受材料影响。

μ+有7次correct→quality_failed、6次反向迁移；μ−为7次与5次。100044除了7+7次sign交换，还损失1个原wrong至quality_failed；wrong净少1不能解释成仅改变1条轨迹。相同final ordinal且ordered hit IDs相同的共同selected数分别为247150、255585、352787；这一身份定义不额外证明independent seed lineage相同。

![Material encoding joint metrics](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/figures_v1/results/material_joint_metrics.png)

## 详细数据：100044冻结角峰窗口

使用S05独立预留区间`sample_entry≥250000`，保持`stage05_slopes/main_analysis_v4/results/peak_windows.json`中的原窗口；未按corrected图重新调窗。reco-p选择的有限states数为78428→78409。

| 物理面 | 冻结中心 / rad近似 | legacy窗口计数 | corrected窗口计数 |
| --- | ---: | ---: | ---: |
| st1 | −0.00189 | 11541 | 11549 |
| st1 | −0.00071 | 2552 | 2548 |
| st1 | +0.00059 | 1731 | 1733 |
| VetoNu | −0.00037 | 3916 | 3900 |
| VetoNu | +0.00043 | 11892 | 11877 |

窗口按dimensionless tan(theta_y)定义，表中rad仅指小角近似。峰fraction、sideband-subtracted fraction及各自事件bootstrap95区间完整保存在`M/analysis_v2/full/100044/summary.json`。本轮没有由曲线目视重合推断余量为零。

三样本的st1/VetoNu truth-model momentum、slope及validity均逐值不变，因此主对照的角度变化来自reconstruction或选择。角度图同时保留all/reco-p/truth-p三类分布，不把charge纯样本的分母混用于S05全长track图。

![100044 material-encoding angle control](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/figures_v1/results/100044_material_angle_physical.png)

## 执行、验证与产物

重建clusters1578722–1578724全部300 jobs exit0，无恢复重跑；full validator1578752 exit0。analysis_v1 smoke1578751中的100044缺q/p而失败，v2 smoke1578753及full1578754各3组均exit0。绘图1578755、逐图验收1578757均exit0；127项绘制数据与冻结summary精确相同，直方图finite/flows、sumw/sumw2闭合，12页vector PDF与逐图页面的graphics、page boxes及rotation一致。12张PNG已通过两张contact整体目视检查，并全尺寸检查角度、动量、联合指标、迁移及truth-p misID主图；记录在图目录的`validation.json`。

每次提交重新读取schedd广告与queue，排除实际不兼容EOS提交路径的bigbird15/21，选择兼容负载最低的bigbird20；查询时间、选择依据、CPU/内存及磁盘申请保存在各`*_worker/submission.json`。tracking为单event slot、计算线程数与1CPU申请一致。CPU ms/source-event的legacy→corrected描述统计为100044 17.47392→17.37220、μ+10.36500→11.01562、μ−10.87456→13.84804；不同worker条件下不据此报告严格加速比。完整RSS与每job原始资源见full summary。

本轮7个Python文件编译、3项配对/零计数/相关子粒子统计回归、64个本地文档链接及`git diff --check`通过，记录在`C/verification_v11`。图集增加的链接随本阶段提交再次核对。

主要实现：[材料全量分析](../analysis/material_encoding_reconstruction.py)、[绘图与逐值验收](../analysis/material_encoding_figures.py)、[归约任务冻结](../scripts/prepare_material_encoding_analysis.py)、[绘图任务冻结](../scripts/prepare_material_encoding_figures.py)。

```bash
python3 scripts/prepare_material_encoding_analysis.py --root "$C" --revision <new-analysis>
python3 scripts/prepare_material_encoding_figures.py --root "$C" --analysis "$M/analysis_v2" --revision <new-figures>
python3 scripts/submit_momentum.py --tasks <frozen-tasks.json> --library <matching-lib> \
  --name <unique-name> --stage S06-material --cpus 1 --memory-mb 4000 --disk-gb 4 \
  --exclude-schedd bigbird15.cern.ch --exclude-schedd bigbird21.cern.ch
```

[12页材料控制图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/figures_v1/results/material_encoding_report.pdf) · [逐图索引与bin counts](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/figures_v1/results/figure_index.json) · [完整覆盖验收](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/material_encoding_control_v1/full_validation.json)。

[96页S06完整汇报图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage06_charge/20261113_isolated_v1/report_v6/charge_stage06_report.pdf)已由cluster1578759组装、exit0；全部页vector content、MediaBox、CropBox和rotation与来源单页精确相同，96张来源PNG完整性通过。前84页的既有目视记录与本轮12页目视记录均在`report_v6/validation.json`引用，未用重新栅格化替代vector页。

## 阶段性结论与下一步计划

事实：修正材料编码不改变本批μ± selected错号成员，100044的净变化很小且存在双向sign迁移；角峰窗口及核心变化较小，truth model精确不变。该结果排除了“这一个编码修正能消除当前主要剩余结构”的解释，未把历史错误材料图重新定义为正确配置。

继续使用修正map对原59wrong+60normal做真实Calypso详细audit，要求58字段与各自corrected full shard精确闭合，再检验旧measurement/seed归因是否保持。完整共同有限先验的CKF非线性响应仍未测量；确定性transport闭合不代替该项。S05主要来源的独立验证保留，S06继续推进Q5及剩余measurement机制，不将本轮材料对照作为全目标完成。
