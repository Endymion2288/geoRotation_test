# 阶段 03：旋转版原始物理配置的 50,000 事件基线

日期：2026-11-11。五类样本各 10,000 events、共 50,000 events 和 70,000 generated signal 粒子均完成重建；50/50 shard exit0，覆盖检查无缺失。当前旋转版的合格匹配 track 仍存在低动量尾与电荷误判；μ−、μ+ 的边缘误判率不同。角度先恢复物理分量并记录选择迁移，峰的统计显著性及来源留待 S05。

## 阶段目标

在冻结输入、frame、truth、matching 和 Good Track 定义下，取得可重放的原版 baseline：紧凑 particle/track 表、stage 来源、cutflow、tail/misID/efficiency 和区间、原始/物理角度图，以及 CPU/I/O 预算。S02 的物理实验开关全部关闭，避免把接口修正效果混入原版结果。

## 工作内容

使用 `configs/tracking_runner.py` 的 3-station forward 入口；FASERNU-04 / OFLCOND-FASER-05，masked layers `[0,1,2,3,4,5]`，一进程一个串行 event slot。输入文件与 ranges 见 `configs/samples.json` 和 `pilot_current_v2/tasks.json`。五类输入生产 conditions 均为 OFLCOND-FASER-04，本阶段保留历史重建使用的 -05，未把二者混同。

primary material 为当前 `faser_reco.py` 默认的未旋转 alma9 map 的冻结副本 `stage02_closure/material-current-source.json`。`CorrectInterfaces=false`、`CorrectCovariance=false`、`ResetEventMaps=false`。启用紧凑 C++ audit，关闭 ROOT diagnostics 和完整 states，不保存批量 xAOD。完整 Calypso 起点与 audit patch 见 S02；运行库和小型 worker 配置冻结于 `stage03_baseline/software_v2/`。该快照记录的主仓库 HEAD 是冻结时的 15b9c22，精确运行文件另外保存，不能仅按那个 commit 重建 worker。

## 详细数据

### 执行、身份与计数闭合

50 jobs × 1000 events，CPU1、memory3500MB、disk4GB、EL9、microcentury。提交前比较普通 bigbird schedd；bigbird20 的 running/idle/capacity 为 1/2/10055，总队列最低且 EL9/EOS/CVMFS probe 已成功。选择记录和完整候选负载位于 submission manifest；正式 cluster=1573730。

首轮 cluster1573723 的 50 jobs 在软件 setup 阶段 exit127，未处理事件；原因是 worker 在 CVMFS 初始化时启用了 shell errexit。关闭 setup 期间的 errexit/nounset 后，单事件 cluster1573728 exit0，随后重新提交全部 pilot。失败尝试不进入物理分母。

合并发现同一输入 RDO 也会复用 run/lumi/event。100044 source entry0 与 entry5000 均标为 run100044/lumi0/event0，但 p_truth 分别为239.484334与783.080658 GeV，vertex也不同。稳定身份采用 DSID+file_index+source_entry，barcode在事件内定位粒子；stage/seed/candidate/state在事件内定位运算。合并检查输入 ranges不重叠、每条源事件恰有一条header、所有stage身份一致和lineage唯一。bootstrap按源事件分组，A′两子粒子一起抽样。

| DSID | events | signal truth | fiducial+H18 truth | final tracks | 零 track events | 文件内复用生产 ID 的额外 entries |
| --- | --- | --- | --- | --- | --- | --- |
| 100044 | 10000 | 10000 | 8274 | 9687 | 324 | 5000 |
| 100116 | 10000 | 10000 | 5794 | 9104 | 897 | 0 |
| 100117 | 10000 | 10000 | 5975 | 9068 | 934 | 0 |
| 110033 | 10000 | 20000 | 16737 | 19556 | 360 | 5000 |
| 110051 | 10000 | 20000 | 19344 | 22755 | 438 | 9000 |

每个样本的 event sumw=sumw2=10000；当前使用的 nominal event weights 全为1，计数率与加权率相同。所有源事件保留，包括零 track、未匹配与未重建 truth。

| DSID | seeds | found candidates | selected | finder errors | refit fallback | 未通过 signal 主匹配的 final | signal 重复 track 超额 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 100044 | 9916 | 15728 | 9687 | 4 | 723 | 72 | 0 |
| 100116 | 9230 | 13652 | 9104 | 0 | 1374 | 64 | 0 |
| 100117 | 9178 | 13075 | 9068 | 1 | 1265 | 70 | 0 |
| 110033 | 31203 | 530373 | 19556 | 617 | 487 | 4947 | 56 |
| 110051 | 26030 | 847114 | 22755 | 1806 | 315 | 8866 | 73 |

finder errors按seed计数，fallback按最终track计数，不是作业失败率。未通过signal主匹配的track包含次级粒子、低纯度或无匹配，不能一律称为fake。五类样本的final无效动量计数均为0。

### 动量尾、效率与电荷误判

主表：status1、目标PDG、`0<barcode<200000`，fiducial+H18，SDO majority purity≥0.7，并要求最终Good95。Good95与fiducial完整定义见S01；动量参照是HepMC production p。每个truth使用与charge/p无关的固定最佳匹配规则。表中区间为Wilson95；事件bootstrap95另存summary，A′相关性以事件bootstrap为主解释。

| DSID | Good95 n | p_reco/p_truth<0.5 | tail % [Wilson95] | wrong charge | misID % [Wilson95] | Good95效率 % [event bootstrap95] |
| --- | --- | --- | --- | --- | --- | --- |
| 100044 | 6965 | 791 | 11.357 [10.633, 12.123] | 369 | 5.298 [4.796, 5.849] | 84.179 [83.309, 84.989] |
| 100116 | 4868 | 75 | 1.541 [1.231, 1.927] | 5 | 0.103 [0.044, 0.240] | 84.018 [83.086, 84.933] |
| 100117 | 5047 | 110 | 2.180 [1.812, 2.620] | 34 | 0.674 [0.482, 0.940] | 84.469 [83.564, 85.439] |
| 110033 | 12325 | 1510 | 12.252 [11.684, 12.842] | 265 | 2.150 [1.909, 2.421] | 73.639 [72.883, 74.338] |
| 110051 | 12460 | 4182 | 33.563 [32.739, 34.397] | 1695 | 13.604 [13.013, 14.217] | 64.413 [63.647, 65.204] |

matching歧义单独检查：SDO的一个cluster可同时接受多个truth贡献，故贡献数之和不要求等于hit数。定义次多数barcode也贡献≥70% hits为歧义，名义结果仍按固定majority/barcode规则归属。五类主Good95样本的歧义数依次为0/0/0/0/9；110051排除9条后，n=12451、tail=4173（33.515380%）、wrong=1690（13.573207%），其余主样本不变。未作Good95选择的all_matched中，110033/110051分别有43/132条歧义，后续不能将这些track作为清洁truth-hit因果个例。

| DSID | p ratio Q16 / median / Q84 | p ratio≥2（图外） | wrong∩tail |
| --- | --- | --- | --- |
| 100044 | 0.59278 / 0.92623 / 1.18978 | 302 | 147 |
| 100116 | 0.76789 / 0.91782 / 1.10296 | 26 | 1 |
| 100117 | 0.75103 / 0.91964 / 1.10943 | 58 | 2 |
| 110033 | 0.58345 / 0.93260 / 1.14641 | 203 | 128 |
| 110051 | 0.22493 / 0.78105 / 1.24140 | 837 | 884 |

这些tail率不是同sensor的分辨率或纯算法故障率：production→tracker的真实能损和电子bremsstrahlung也包含在p ratio中。未测量每个sensor处精确p；S01明确了补充G4记录的可用字段和退出故障。

### 最终候选路径上的首次已观测异常

将最终合格track追溯到其原始seed与candidate，不使用每个stage独立best match拼接。下表为最终tail/wrong样本中首次观测到相应条件的stage；found是CKF返回后，refit是ACTS refit返回后。某一早期异常可能被暂时修复后再出现，故此表不直接指定最终异常的根因。

| DSID | final tail：seed / found / refit | final wrong：seed / found / refit |
| --- | --- | --- |
| 100044 | 33 / 726 / 32 | 8 / 354 / 7 |
| 100116 | 9 / 59 / 7 | 0 / 5 / 0 |
| 100117 | 21 / 84 / 5 | 0 / 30 / 4 |
| 110033 | 841 / 621 / 48 | 158 / 106 / 1 |
| 110051 | 2065 / 2014 / 103 | 977 / 684 / 34 |

100044的791条合格tail中，726条在CKF返回时首次出现低比值，33条在seed已出现、32条在refit首次出现；不能只研究seed或只看最终refit。μ−的34条wrong中，30条首次在found、4条首次在refit；下一步要在同一candidate的predicted/filtered/smoothed中定位具体更新。

`analysis/select_cases.py`已选出32条可重放track，保存为`merged_current/representative_cases.json`，覆盖五类样本中实际存在的normal、seed/found/refit首次tail与wrong类别。要求fiducial+H18+Good95且final matching无上述歧义；每个DSID/category固定取源顺序首条，normal另要求0.9<p_ratio<1.1、正确charge及50<p_truth<250GeV。这是调试清单，不是随机统计子样本。下列6条均来自file_index=0：

| DSID / 类别 | source entry | barcode | seed/candidate | p_truth / GeV | p_final / GeV | final purity |
| --- | --- | --- | --- | --- | --- | --- |
| 100116 normal | 2 | 10001 | 0/1 | 51.792585 | 48.749897 | 1.000000 |
| 100044 seed tail | 670 | 10001 | 0/0 | 336.625313 | 121.850057 | 1.000000 |
| 100044 found tail | 21 | 10001 | 0/0 | 1708.984924 | 262.020357 | 1.000000 |
| 100044 refit tail | 99 | 10001 | 0/0 | 484.126982 | 241.492984 | 1.000000 |
| 100117 found wrong | 1645 | 10001 | 0/0 | 114.341886 | 17.300857 | 0.923077 |
| 100117 refit wrong | 504 | 10001 | 0/0 | 424.364045 | 2710.199954 | 1.000000 |

entry99的found/refit p_ratio为0.503721/0.498831，属于阈值附近跨越；entry504的charge则在refit从−1变+1，p_ratio从0.506090变6.386512。这两例的异常幅度与含义不同，S04分别检查，不将“跨越0.5阈值”等同于大幅动量损失。此阶段未对32条个例重新运行完整states。

### Alma比较与μ±差异

100116为μ+、100117为μ−。主定义下总率差为0.570956个百分点，事件bootstrap95为[0.333173,0.829263]个百分点；两独立样本的Fisher双侧p=2.65096×10⁻⁶。这是两个输入谱的边缘总率差，未隔离输入p/position/slope、接受度和选择，也未测量同输入的完整旋转前后净效应。

另提供Alma radius100 proxy：保留fiducial，使用当前质量条件、rmax<100mm、exactly one Good100，取消额外H18；仍采用独立SDO匹配。μ+为8/5358=0.149309%，μ−为34/5534=0.614384%，差0.465074个百分点，bootstrap95 [0.244078,0.701767]。Alma报告new tracking是0.197%/0.359%；原Twiki/selection尚未恢复，geometry/material/conditions/代码也有明确差异，所以本表不是历史数值的精确复现。

| 物理p_truth / GeV | μ+ wrong/n | μ− wrong/n |
| --- | --- | --- |
| 0–50 | 0/1363 | 0/1264 |
| 50–100 | 0/1662 | 2/1359 |
| 100–250 | 1/1459 | 4/1630 |
| 250–500 | 1/331 | 14/604 |
| ≥500 | 3/53 | 14/190 |

各p分箱Wilson区间在summary中；bin内位置与方向谱也需匹配，不能仅凭这张分箱表归因。

### 角度原始定义、物理定义与选择迁移

物理`s_y=py_A/px_A=py_F/pz_F`；角度为atan(s_y)。legacy分量公式`atan(py_A/pz_A)`只用于复现未反变换分量的输出误用，不等同于完整旋转前软件的物理角度。所有表来自同一Good95主样本。

| DSID | physical s_y Q16 / median / Q84，单位10⁻³ | physical abs(theta_y)≥10mrad | legacy abs(angle)≥10mrad | exact s_y=0 | reco-p 50–250GeV n / truth-p同窗n |
| --- | --- | --- | --- | --- | --- |
| 100044 | -3.4845 / -0.7450 / 0.7206 | 189 | 6924 | 0 | 2908 / 2493 |
| 100116 | -8.4133 / -2.2104 / 1.1416 | 661 | 4852 | 0 | 2941 / 3121 |
| 100117 | -0.8281 / 1.8889 / 7.8400 | 616 | 5018 | 0 | 2887 / 2989 |
| 110033 | -1.7626 / -0.0503 / 1.7189 | 336 | 12294 | 0 | 6299 / 6310 |
| 110051 | -0.4190 / -0.0855 / 0.2803 | 26 | 12353 | 0 | 3374 / 1401 |

100044的reco-p窗口保留2908粒子，truth-p同窗为2493；110051分别3374与1401，显示选样迁移显著改变角度的组成。没有exact zero积累，但这不能排除其他离散值或窄峰。

五张英文图每张含stage p ratio、raw/physical/production-truth角度和reco/truth momentum窗口。stage图各自使用独立best match、分母可不同，只用于分布比较；上面的lineage表才用于追踪同一candidate。production truth在磁场前，μ±方向差可含正常弯曲；它与首个reco state也不在同一物理面。当前角度图不作为同面angular resolution或非物理多峰显著性的证据，S05需传播至共同st1/VetoNu plane并作条件分解、bin/KDE/Dip检查。

![100044 baseline](../../geoRotation_test_results/stage03_baseline/merged_current/100044_baseline.png)

其余图：[100116](../../geoRotation_test_results/stage03_baseline/merged_current/100116_baseline.png)、[100117](../../geoRotation_test_results/stage03_baseline/merged_current/100117_baseline.png)、[110033](../../geoRotation_test_results/stage03_baseline/merged_current/110033_baseline.png)、[110051](../../geoRotation_test_results/stage03_baseline/merged_current/110051_baseline.png)。

### 相同10事件的material控制

110051首文件entries0–9、相同软件与全部其他开关，仅切换材料。cluster1573736/1573733/1573734对应current/historical/consistent，提交前bigbird20队列总量分别15/13，低于bigbird21的17/22；EL9/EOS/CVMFS已验证。

| map | surfaces | wall/s（含初始化） | final tracks | matched signal | 主Good95粒子 | wrong | tail<0.5 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| current，未旋转 | 22 | 32.791 | 25 | 14 | 11 | 2 | 2 |
| historical，已旋转 | 12 | 35.078 | 21 | 14 | 11 | 0 | 2 |
| consistent，当前map严格旋转 | 22 | 35.432 | 21 | 14 | 11 | 0 | 2 |

小样本支持材料配置会改变候选与电荷结果，需要后续扩大同输入配对；2个wrong变0不能代表全样本修复率。110051批量耗时更高，但这10事件没有显示current map更慢，不能把完整样本的高CPU直接归因于map旋转。

### 资源、压缩产物与重放

| DSID | primary Σ(user+system) /s | CPU s/event（含每shard初始化） | max RSS/MiB | audit/MB（十进制压缩） |
| --- | --- | --- | --- | --- |
| 100044 | 271.37 | 0.027137 | 1005.77 | 18.277 |
| 100116 | 176.06 | 0.017606 | 982.22 | 14.039 |
| 100117 | 200.79 | 0.020079 | 1322.70 | 14.153 |
| 110033 | 1539.43 | 0.153943 | 1486.50 | 127.459 |
| 110051 | 8849.67 | 0.884967 | 1203.29 | 190.213 |

primary合计CPU=11037.32s（3.066h），不含失败setup、probe、控制和分析；所有作业在1h预算内。最终归档保留5个压缩audit、共364.141212MB，附member offset/source ranges索引和一个job-record归档；长期xAOD为0。`scripts/archive_pilot.py`按gzip member直接合并，校验字节范围后清理原shard重复文件。50k完整验证后清理了旧40k部分分析的7个重复文件（15119991bytes）。分析可直接读索引归档，或恢复原shards重放。

最后一次从`audit_current/`完整复算exit0，wall=149.67s、user/system=141.18/2.25s、maxRSS=923072kB（901.438MiB）。结果为events=50000、missing_or_failed_shards=[]，每个源事件的stage与最终lineage唯一性断言全部通过。去掉新增的`matching_ambiguous/unambiguous`统计字段后，`summary.json`与归档前保存的`summary_before_matching_qa.json`严格相等，包含全部名义cutflow、发生率、分位数和bootstrap区间。

```bash
python3 scripts/prepare_condor.py --output <pilot_dir> --submit-dir ~/condor/geoRotation_test/<task> --mode current
condor_submit -name bigbird20.cern.ch <submit_dir>/jobs.sub
python3 analysis/pilot.py --input <results>/stage03_baseline/audit_current --output <results>/stage03_baseline/merged_current
python3 analysis/select_cases.py --particles <merged>/particles.jsonl.gz --output <merged>/representative_cases.json
python3 scripts/archive_pilot.py --restore-from <results>/stage03_baseline/audit_current --output <scratch>/restored_shards
python3 -m unittest discover -s tests -v
```

实际运行入口使用冻结software_v2/code；上述prepare命令展示可复用入口，精确任务为归档index中的50个ranges。合并的每条particle包含未重建truth、物理与legacy量、最终候选lineage、文件与entry。所有cutflow和区间在summary.json；figure溢出、missing shard和matching失败均显式计数。

验证：`python3 -m unittest discover -s tests -v`的4项检查全部通过，覆盖charge/p独立matching、F/A角度与未重建分母、重复生产ID保留，以及最终candidate祖先与独立best-stage的区别。五张英文图完成视图检查；Python AST、shell语法、Markdown本地链接/表格/代码块与`git diff --check`的最终检查记录见`merged_current/final_validation.json`。

产物根目录：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage03_baseline/`。主交付为`merged_current/`、`audit_current/`、`software_v2/`、`material_controls.json`和`material_controls_results.json`；提交记录保留于`pilot_current_v2/`及归档。AFS提交文件位于`/afs/cern.ch/user/x/xcheng/condor/geoRotation_test/stage03_pilot_current_v2/`。

## 阶段性结论

**事实**：50,000源事件全部覆盖，70,000 signal truth均进入particle母表。当前旋转版的Good95主样本中，100044动量尾11.357%，μ+/μ−误判0.1027%/0.6737%；A′两样本也有明显尾部。最终候选路径显示seed、CKF返回与refit都需要检查，主要新增环节随样本而不同。

**事实**：物理角度与raw分量比差异很大，reco-p选样迁移可量化；当前证据尚未给出剩余多峰的非物理显著性。未复现完整历史Good Track代码，也未运行完整未旋转G0、backward、4-station或merged配置。

**待检验**：μ±总率差来自哪些输入谱/接受度或算法环节；同面真值下的p/angle误差；material/接口修正能改善多少合格track异常。小样本控制和首次出现stage是后续因果实验的输入。

## 下一步计划

1. S04按representative_cases重放合格seed/found/refit首次异常，保存同一lineage完整states；truth参数注入与truth hits实验分开，逐次核对propagation、measurement update、smoothing和candidate替换。验收为首次数值异常、对应操作及独立样本修复效果。
2. S05统一st1/VetoNu物理面，比较SiHit方向与production transport，固定truth-p/charge/module/hit pattern，检验binning/KDE/峰窗口稳定性；不以Gaussian分量数判定根因。
3. S06先对μ±的truth p、位置与方向做共同支持域比较，再扩充有效选中粒子与冻结bin统计；定位found/refit符号翻转并执行一致G0/G1及q/B对称性控制。
4. 同输入批量控制当前/严格旋转material、chart Jacobian、refit plane，单次只改一项；保留原版和独立验证文件。G4补充模拟先解决退出134，再生成精确同sensor p truth。
