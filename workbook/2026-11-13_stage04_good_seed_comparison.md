# S04：Purity 基准上新增 good-seed 条件的几何前后对比

结果目录：`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2`，下文记为 `R`。

## 阶段目标

以用户确认的 **candidate-isolated smoothing 开启、strict final purity > 0.7、native truth** 为后续动量画图基准。在此基础上增加 seed 初值接近 MC truth 的条件，对比 beforeGeo/afterGeo，并测量 afterGeo smoothed residual 的负尾和左右不对称。

参考 [problem.pdf](../docs/problem.pdf) 第3页左下角的三条 seed limit；本页图像保存为 `R/evidence/problem_page3.png`。保留 [GitHub test_pResolution](https://github.com/zsm2024/FASER_tracking/tree/main/test_pResolution) commit `33ca2281b961935456dfd4dc9fc06235e24b551e` 原生的 residual、cuts、100 bins、[-1.2,1.2]、raw Entries、颜色和九层布局。

## 工作内容与冻结设置

两组完全沿用[上一轮 Purity 基准](2026-11-12_stage04_native_purity_comparison.md)的材料、covariance、mask、事件分组与 native writer truth。beforeGeo 基底为本地 master `1eb885ae22299fb39bb352bead6005d41b0f953e`；afterGeo 为旋转分支 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。实际源码分别在 `geoRotation_calypso_smoothing_before`、`geoRotation_calypso_smoothing_after`；相同诊断 writer 补丁也同步至主诊断源码树 `geoRotation_calypso`，本次运行使用冻结的两套匹配库。

| 设置 | beforeGeo | afterGeo |
| --- | --- | --- |
| initial variance inflation | `[100,100,100,100,1000]` | `[100,100,100,10000,100000]` |
| material JSON | `stage02_closure/material-current-source.json` | `/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json` |
| isolated smoothing | 开 | 开 |
| final purity | 严格 >0.7 | 严格 >0.7 |
| audit / reset-event-maps | 关 / 关 | 关 / 关 |
| 新增 seed limit | 同下述三条 | 同下述三条 |

共同使用 3station forward、mask `[0,1,2,3,4,5]`、FASERNU-04、OFLCOND-FASER-05。GitHub cuts 仍为 `chi2Sum/NDF<20 && rMax_smt<90 && NDF>9`；其上增加 `final_purity>0.7 && seed_good`。两套原生 covariance 不同，比较仍是整体版本比较。

新增 seed 条件为：

```text
abs(p_seed / p_seed_truth - 1) < 0.1
qOverP_seed * qOverP_seed_truth > 0
abs(theta_seed - theta_seed_truth) < 0.0001 rad
```

三条必须同时满足。动量 magnitude 不将错误电荷符号混入 residual；电荷符号单独判断。实现用 `abs(abs(t_QoverP/track_QoverP)-1)<0.1`，两者均为 GeV⁻¹。

seed 参数来自真实 Calypso `RootSeedWriterTool` 输出的 `track_QoverP`、`track_theta`；truth 来自该 writer 的 `t_QoverP`、`t_theta`，即 seed majority HepMC particle 经实际 Calypso extrapolator 传播至 seed reference surface 的参数。beforeGeo 使用 native FASER/ACTS frame，afterGeo 使用旋转 ACTS frame，分别比较同 frame、同 surface 的 theta；这与直接对 production vertex momentum 求 theta 的旧补图定义存在参考面的差别。theta为各自native ACTS坐标系的polar angle，旋转前后该阈值约束的方向分量不同，因此不能把通过样本差异当成同一物理角度窗口。该 seed 条件只决定样本是否通过，最终九层图继续使用原 writer 的 native `t_eQOP`（SiHit-linked HepMC momentum），未替换 truth 分支。

### 精确 seed 归属与源码修正

原 `RootTrajectoryStatesWriterTool` 未输出 seed id；按 event 和 barcode 查找任意好 seed 会误接同一粒子的其他 seed。本次增加只读 `seed_nr`，从实际 CKF `trackGroup` 读取，按 **input file + source entry + seed_nr = seed writer seedIndex** 联结，event number 只作顺序核对，不作为独立事件身份。

检查发现两个原生 CKF 均声明 `unsigned int nSeed=0`，随后写 `seedNumber(track)=nSeed-1`，但 `nSeed` 从未递增，故标签恒为 `UINT32_MAX`。修正为当前 finder 循环的 `iseed`；该标签在本次配置中不参与 fitting、candidate ranking 或 selection，alignment writer 关闭。

可复用补丁为 [0018_writer_seed_lineage.patch](../patches/0018_writer_seed_lineage.patch) 和 [0019_native_ckf_seed_index.patch](../patches/0019_native_ckf_seed_index.patch)，另保留已确认的 [0010 isolated smoothing](../patches/0010_candidate_isolated_smoothing.patch)。主诊断树的 audit patch 已使用 `iseed`，只同步0018，不重复应用0019。两个原生基底重建 matching tracking package，几何库复用未变构建；最终 `R/seed_software_before_v2` 与 `R/seed_software_after_v2` 的 `software.json` 保存完整累计 diff、commit、28个 translation units、共享依赖和重用几何来源。补丁顺序应用后8项文件比较全部相等：`R/patch_validation/validation.json`。

实现入口为 [reference_tracking.py](../configs/reference_tracking.py)、[native_purity.py](../analysis/native_purity.py)、[run_reference_task.py](../scripts/run_reference_task.py)。`--seed-unbiased` 只启用 seed 诊断，worker 在 scratch 中应用 gate 并输出 compact histograms/逐事件计数，不保存全量 states 或 xAOD。

## 详细数据：smoke 与重建不变性

两个修正后的20-event smoke均成功，clusters 为 afterGeo `1574079`、beforeGeo `1574080`。新增标签和 seed writer 开启时，仍用旧 Purity 宏产生控制 histogram，和上一轮相同输入 entries `[0,20)` 核对：**54个 histogram、5508个 bins（含 flow）、全部 raw Entries 逐项相等**。见 `R/smoke_validation.json`。同 barcode 多 seed 的组在每组20-event中均为9组，因此精确 seed 编号关联有实际必要。

| 配置 | source events | writer tracks | GitHub cuts | Purity 基准通过 | seed总数 | p cut通过 | 同电荷通过 | theta cut通过 | seed三条全过 | 最终 good-seed轨迹 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| beforeGeo | 20 | 44 | 30 | 25 | 59 | 4 | 42 | 5 | 0 | 0 |
| afterGeo | 20 | 42 | 30 | 23 | 59 | 4 | 42 | 7 | 1 | 1 |

20-event只用于实现检查；物理结论使用全量。43项回归检查通过，见 `R/unit_tests.log`，包含 strict purity 边界、seed电荷与动量分离、theta rad阈值、invalid输入，以及重复event number / 同barcode多seed的真实ROOT联结与friend落盘；`git diff --check`通过。

早期 smoke 失败及修正：ROOT读取 seed file 后 current directory改变，必须显式恢复输出文件再写 friend tree；初版 writer暴露了 `UINT32_MAX` 的旧来源标签。重新提交使用独立输出目录，未混入全量。失败输出和日志保留于 `R/smoke*`、`R/lineage_smoke/`；最终依据为 `R/lineage_smoke_v2/`。

## 详细数据：全量执行

冻结输入为110051的全部20个RDO文件，每文件1000个`CollectionTree` entries、完整范围`[0,1000)`，每组20,000个独立源事件。worker逐文件检查实际entries，归约检查file/entry唯一键、event ranges、加载库、covariance、isolation与seed诊断开关。输入清单及配置见`R/inputs.json`。

HTCondor全量clusters：afterGeo `1574081`，beforeGeo `1574082`，各20个作业，1 CPU、3500 MB、12 GB scratch、workday，OMP/BLAS/MKL线程为1。每次提交前查询schedd和队列，选择兼容且最空闲的`bigbird20.cern.ch`，排除已记录EOS兼容问题的bigbird21；查询时间、候选负载、队列、submit文件、cluster与资源均保存在相应`*_worker/submission.json`及load/queue文件。交互节点`lxplus901.cern.ch`，可见Tesla T4，本次CPU tracking未使用GPU，见`R/interactive_environment.json`。

40个最终纳入的任务全部成功，逐文件均为1000 events、范围`[0,1000)`；两组分别有20,000个唯一file/entry键且一一配对，无遗漏、重叠或重复计入。writer tracks、GitHub cuts和Purity基准通过数逐文件与上一轮一致。所有54个histogram的raw Entries、finite计数和flow独立闭合，见`R/reductions/summary.json`；bin edges、counts、sumw/sumw2见`bin_counts.json`。

| 配置 | writer tracks | GitHub cuts | Purity基准通过 | 再加seed条件 | 相对Purity基准保留率 | smoothed raw states | 有限smoothed residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| beforeGeo | 44,238 | 36,178 | 26,530 | 635 | 2.3935% | 11,193 | 11,193 |
| afterGeo | 41,940 | 33,597 | 25,306 | 951 | 3.7580% | 16,866 | 16,862 |

两组全部smoothed invalid reconstruction均为0；afterGeo有4个缺失native truth的states，按原宏保留在raw Entries和flow中。通过seed gate的轨迹中，native writer majority与final majority不同的两组均为0。

每组初始seed均为52,338个；p cut通过6,317、同电荷38,866、p与电荷同时通过6,023。theta cut分别通过5,584/7,761，三条同时通过1,033/1,200；不是所有好seed都得到最终输出轨迹。每组seed writer有1,974条无法形成有效momentum比较、1,973条无法形成有效theta比较；这类seed不通过三条共同条件。每组共有8,237个同source entry/majority barcode下的多seed组，已按真实seedIndex消除联结歧义。

最终通过的轨迹中，seed majority与final majority不同的分别为110/635（17.32%）、228/951（23.97%）。这是实际seed来源联结后观测到的truth归属变化，当前未额外加同barcode条件；因此“好seed”并不保证CKF最终保留的是该seed原本对应的粒子。该计数不把电荷或重建动量用于主truth匹配。

### smoothed负尾与seed筛选作用

分母为各层有限native residual states；负尾为`residual<-0.5`，更远左尾为`residual<-0.8`。率差和95%区间按同一20,000个源事件配对，事件内多条track/state共同聚类计算，方法为source-event cluster ratio delta method。

| 层 | beforeGeo <−0.5 / % | afterGeo <−0.5 / % | after−before / pp | 95%区间 / pp |
| --- | ---: | ---: | ---: | --- |
| S1 L0 | 26.0176 | 2.5899 | -23.4277 | [-26.9395, -19.9159] |
| S1 L1 | 26.3911 | 2.5067 | -23.8844 | [-27.4022, -20.3666] |
| S1 L2 | 26.2321 | 2.6553 | -23.5768 | [-27.1025, -20.0511] |
| S2 L0 | 26.0835 | 2.4573 | -23.6262 | [-27.1404, -20.1120] |
| S2 L1 | 26.1809 | 2.4155 | -23.7655 | [-27.2727, -20.2582] |
| S2 L2 | 25.5387 | 2.6525 | -22.8862 | [-26.3723, -19.4001] |
| S3 L0 | 25.9289 | 2.3479 | -23.5810 | [-27.1423, -20.0197] |
| S3 L1 | 25.5390 | 2.0529 | -23.4860 | [-27.0252, -19.9468] |
| S3 L2 | 25.3247 | 2.5174 | -22.8073 | [-26.3135, -19.3011] |

S1 L0相对于上一轮Purity基准的单项变化：

| 配置 | 基准 <−0.5 / % | 新增seed后 / % | 变化 / pp | 95%区间 / pp | 基准 <−0.8 / % | 新增seed后 / % |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| beforeGeo | 36.5994 | 26.0176 | -10.5819 | [-13.9455, -7.2182] | 12.1204 | 8.2203 |
| afterGeo | 10.5236 | 2.5899 | -7.9338 | [-8.9798, -6.8877] | 4.0172 | 1.1099 |

afterGeo的<−0.5负尾相对减少75.39%，<−0.8由4.0172%降至1.1099%。这一结果测量的是新增truth seed选择的总体效应；未更改seed均值、hits或Kalman更新。

### 左右是否对称

以S1 L0为例，同一组、同一有限states分母，左侧定义`residual<-a`，右侧定义`residual>a`。差值为左减右，区间仍按源事件聚类；相同事件内左右计数的相关性包含在误差中。

| 配置 | a | 左侧 / % | 右侧 / % | 左−右 / pp | 95%区间 / pp |
| --- | ---: | ---: | ---: | ---: | --- |
| beforeGeo | 0.1 | 65.2035 | 13.4078 | 51.7957 | [46.1961, 57.3953] |
| beforeGeo | 0.2 | 42.6177 | 11.4924 | 31.1253 | [25.9409, 36.3097] |
| beforeGeo | 0.5 | 26.0176 | 9.4174 | 16.6002 | [12.2181, 20.9822] |
| afterGeo | 0.1 | 54.0698 | 6.2368 | 47.8330 | [44.1069, 51.5591] |
| afterGeo | 0.2 | 12.7907 | 4.1755 | 8.6152 | [6.0476, 11.1828] |
| afterGeo | 0.5 | 2.5899 | 1.5328 | 1.0571 | [-0.2346, 2.3488] |

**事实：新增seed条件显著降低afterGeo远左尾；但负偏尚未消失。** 共同0.024宽分箱下，S1L0 smoothed主峰最高bin由Purity基准的`[-0.024,0]`移至beforeGeo的`[-0.120,-0.096]`、afterGeo的`[-0.096,-0.072]`；这是histogram modal bin，不是拟合峰位，见`R/reductions/mode_check.json`。 S1L0的±0.5远尾差为1.0571 pp，95%区间[−0.2346,2.3488] pp，当前不能确认该阈值下有非零左右差；±0.2和±0.1窗口外的左侧仍分别比右侧高8.6152与47.8330 pp，区间均不含0。因此不能把远尾变少解释为完整分布已左右对称。

**尚未检验的原因假设：** good seed之后仍可能在组合跟踪中更换truth归属；也可能存在固定正确hits路径下的更新/传播偏差。当前新增seed筛选本身不能区分这两者。theta条件为两套native frame下的polar angle，也不构成严格相同物理角度窗口的旋转控制。

### 全量作业与慢作业处理

最终全部文件成功。原afterGeo首文件`1574081.0`所在`b9p09p5581.cern.ch`持续计算但进度明显慢于其余作业；另在重新查询schedd后提交相同冻结输入/配置的备用`1574086.0`，排除该worker。备用任务在`R/straggler_replay/afterGeo/f00_e0000`成功，cutflow与旧基准一致；最终`afterGeo_tasks.json`只纳入这一个首文件结果，原慢作业移除。记录见`R/straggler_resolution.json`，不把重复处理当成新增输入事件。备用提交最初错误地用`-af Requirements`读取了没有TARGET时求值的false；改用`-long -attributes Requirements`恢复raw ClassAd后正常匹配，该更正也保存在submission.json。

最终两个版本的重建总wall分别约8.14/8.68 worker小时，不含构建、失败smoke与未采用的慢作业。执行副本与本仓库源码、两套库的累计source diff共12项核对均相等，见`R/code_validation.json`。

### 正式图集

主图保留原宏raw Entries、P=green/F=red/S=blue、共同分箱、可比纵轴和九层布局。geometry overlay中实线为beforeGeo、虚线为afterGeo：

![Full good-seed geometry comparison](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/110051_geometry_overlay.png)

新增seed条件对afterGeo的形状影响，灰线为已确认Purity基准、红线为再加good seed；两者均按各自raw Entries归一化：

![Full afterGeo seed selection effect](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/110051_afterGeo_seed_effect.png)

镜像图在相同绝对residual下比较负侧和正侧：

![Full afterGeo mirrored residual tails](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/110051_afterGeo_smoothed_mirrored_tails.png)

[完整十页PDF](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/native_purity_comparison.pdf) · [beforeGeo九层图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/110051_beforeGeo_native_nine_layers.png) · [afterGeo九层图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/110051_afterGeo_native_nine_layers.png) · [两组smoothed左尾](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/110051_smoothed_left_tail.png) · [左尾log图](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/110051_smoothed_left_tail_log.png) · [seed初值分布](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/110051_seed_residuals.png) · [beforeGeo筛选作用](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/110051_beforeGeo_seed_effect.png) · [图集索引](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/figures/figure_index.md)。每图均保存PNG与单页PDF。

导图作业`1574087`退出0，合并后54个逐层N标注与冻结raw Entries全部相等。10个单页PDF和合并PDF逐页content、media/crop boxes与rotation全部一致，见`R/figures/figure_validation.json`。全部PNG与十页PDF渲染已目视检查，未见截断/图例遮挡，预览见`R/png_review_contact.png`与`R/pdf_review_contact.png`。

### 执行与复现入口

prepare入口在原基准之上加`--seed-unbiased`，保留原macro diff；实际runner、worker、分析副本都冻结于`R/code/`或`*_worker/`。重建所用库必须为带0018/0019标签修正、并保持0010 isolated smoothing的匹配版本。两个版本均使用`build_native_calypso.py --source <native_worktree> --source-patch <R/before或after_seed_lineage_v2.patch> --reuse-geometry <R/seed_software_before或after> --jobs 4`，实际完整命令在`R/before_seed_build_v2_tasks.json`及`after_seed_build_v2_tasks.json`。

全量冻结任务为`R/beforeGeo_tasks.json`、`R/afterGeo_tasks.json`；后者首文件指向唯一采用的备用结果。直方图和统计重画入口为：

```bash
python3 R/code/analysis/native_purity_figures.py --root R --phase all
```

这里`R`需展开为本页顶部实际绝对路径，执行环境为Athena24.0.41；本次该命令通过`R/all_analysis_tasks.json`在HTCondor中完成，非交互式大规模处理。全部submission清单在`R/all_submissions.json`。



## 阶段性结论

1. 后续基准固定为isolated smoothing开启、strict final purity > 0.7、native truth。新增seed条件单列为诊断cohort。
2. 三条seed条件进一步压低负尾。全量S1 L0的<−0.5比例：beforeGeo从36.5994%降至26.0176%，afterGeo从10.5236%降至2.5899%，配对区间均支持下降。
3. afterGeo远尾左右差的95%区间包含0，但靠近主峰的负偏仍显著，不能称完整分布已经对称。good seed也不能保证最终粒子归属不变：通过样本中afterGeo有228/951条seed/final majority变化。
4. 精确seed编号来源问题已修正，新增诊断与标签不改变基准histogram；全量cutflow与上一轮逐文件相同。

## 下一步计划

固定本次实际seed与候选路径，分开检查seed/final majority保持一致与变化的贡献；对相同truth且正确hits的子集检查CKF逐层负偏是否仍在。验收继续使用本次native residual、相同分箱与源事件区间估计，不把额外truth选择当作动量更新修复。
