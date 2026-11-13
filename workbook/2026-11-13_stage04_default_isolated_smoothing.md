# S04：将 candidate-isolated smoothing 集成为 CKF2 默认行为

## 阶段目标

将已验证的 candidate-isolated smoothing 从研究开关整理为真实 Calypso 的默认功能，避免输出候选的 smoothed parameters 和 covariance 继承其他 CKF 分支的值。为 `sazhang/calypso:geoRotation_rebasedOfficialMaster` 准备最小、可审查的修复。

## 工作内容

目标分支通过 `git ls-remote` 和 fetch 核对，其 HEAD 为 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`，与原 geoRotation 工作目录基底一致。独立 worktree 为 `/eos/home-x/xcheng/FASER/geoRotation_calypso_smoothing_fix`，分支为 `fix/ckf2-isolated-smoothing`。功能修复 commit 为 `1c06c7b65135a88e5ebd864050b10a638b152c98`；验证说明 commit 为 `25147d8028e2e7439700a1395a38884d640b7581`，只补充英文文档，未修改已验证的实现。补丁 0020 包含这两个可用 `git am` 导入的提交。

源码定位：`Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx` 在候选筛选后的 `destProxy.copyFrom(..., true)` 后重新调用 `Acts::GainMatrixSmoother`。`src/CKF2.h` 中 `IsolateCandidateSmoothing` 默认 `true`；显式设为 `false` 可复现旧 shared-state 输出。重算失败返回 `StatusCode::FAILURE`。

ACTS 来源为 AthenaExternals 24.0.41 的 `include/Acts`。`TrackFinding/CombinatorialKalmanFilter.hpp` 按候选依次调用 smoother，并在当次处理时生成参考面参数；`TrackFitting/GainMatrixSmoother.hpp` 沿 `previous` 回溯并原地写入 smoothed mean/covariance；`EventData/TrackProxy.hpp::copyFrom` 在 `copyTrackStates=true` 时复制完整路径。共享祖先的输出值可能已被后续分支改写，因此需要在深拷贝之后重算。

该修复保持候选排序、共享 hit 去重、seed covariance、predicted/filtered states 和 CKF 已生成的 reference-surface parameters 的原流程。它修复 selected track 的逐层 smoothed states；后续 Kalman refit 仍使用原实现。这里定位的是共享 states 的所有权和结果保存问题，不将其推广为所有 ACTS 版本的 smoother 算法错误。

MR 只有五个文件：CKF2 的 `.cxx/.h`、package CMake、真实 ACTS 回归测试 `test/CandidateIsolatedSmoothing_test.cxx`、英文功能说明 `doc/candidate_smoothing.md`。没有加入 purity/good-seed 筛选、truth 处理、variance inflation 调参或 audit 诊断补丁。可复用提交补丁见 [0020_default_candidate_isolated_smoothing.patch](../patches/0020_default_candidate_isolated_smoothing.patch)。

## 仓库整理

原目录 `/eos/home-x/xcheng/FASER/geoRotation_calypso` 的研究改动保留。修改前的 tracked diff、七个未跟踪源码/辅助文件及状态清单保存于结果目录中的 `research_before.patch`、`research_untracked_before.tar.gz`、`research_before.json`。`build/`、`run/` 不归档也不加入 Git。

原目录同步了默认开启、源码注释、回归测试和功能说明。其本地 `package_filters.txt` 增加 `-run/.*` 和 `-build/.*`：原配置重建时会扫描 `run/src` 安装副本并报重复 package；排除后可配置。Git 本地 `info/exclude` 排除 build/run，`README.local.md` 记录修复分支与研究改动的关系。这些目录整理不属于 MR。

研究入口 `configs/reference_tracking.py` 增加显式 `--shared-smoothing` 对照入口；不传开关则继承组件默认值，`--isolate-smoothing` 仍显式开启。reference、tracking 和 slope 入口记录组件实际生效的 `IsolateCandidateSmoothing`，避免把 CLI 未传参数误记为关闭。冻结的历史软件和代码副本保留原行为。

## 详细数据与验证

结果目录：`R=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_default/20261113_mr`。实际主机 `lxplus901.cern.ch`，GPU 为 Tesla T4；本次编译及 ACTS 重建使用 CPU。环境为 Athena 24.0.41、`x86_64-el9-gcc13-opt`。

干净 MR 源码按现有 Calypso compile/link recipes 编译 28 个 tracking translation units，两个 geometry 库复用相同基底的冻结 afterGeo 构建，其余依赖复用原 `run/`。`R/software_clean/software.json`、`compile_recipes.json`、`link_recipes.json` 和 `mr_source.patch` 记录实际来源；构建后核对提交树与构建 patch 一致。Gaudi metadata 从新库重新生成，组件检查得到 `IsolateCandidateSmoothing=True`。

新增 C++ 回归使用真实 `VectorMultiTrajectory`、`TrackContainer` 和 `GainMatrixSmoother`。构造两个共享两层祖先、末端参数和 covariance 不同的候选，先证明参数及完整 covariance 均发生覆盖，再深拷贝各自路径重算，验证解析 RTS 解、两种候选顺序、先保存候选不变、源树不变、predicted/filtered/Jacobian 不变、参考面参数不变及单 state 端点。直接按匹配 compiler recipe 编译和运行通过；CMake 生成的测试 target 也已构建，正式 CTest **1/1 通过**，测试耗时 0.60 s。日志为 `R/unit_compile.log`、`R/unit_run.log`、`R/ctest_target_build.log`、`R/ctest.log`。研究仓库原有 **43 项 unittest 全部通过**，见 `R/unittest.log`。

真实 Calypso 对照固定 DSID 110051 第一文件 `FaserMC-MC24_FS_Aee_010MeV_1Em4-110051-00000-00009-master_dev-RDO.root` 的 entries `[0,20)`，共一个文件、20 个源 entries。三种模式为新默认、显式 ON 和显式 OFF；统一使用旋转 ACTS geometry、native writer truth、mask `0,1,2,3,4,5`、variance inflation `[100,100,100,10000,100000]`，其他重建设置固定。三次重建在 worker scratch 的临时目录中运行，完成后只保留 compact snapshots、统计、日志和图表。本次是小样本软件回归，不标记为全量。

| 检查 | 实际结果 |
| --- | --- |
| 真实 Calypso 覆盖 | 1 个 RDO 文件，entries `[0,20)`，20 个源事件 |
| writer 输出 | 42 条轨迹，718 个 measurement states |
| 新默认 / 显式 ON | 103 个快照数组逐元素完全相等，包括参数、误差、hit 标识、path length 和 track metadata |
| 新默认 / 显式 OFF | 所有非 smoothed 数组相等；156 个 states 的 smoothed q/p 改变 |
| 与此前冻结 ON 对照 | 27 张 native P/F/S 直方图，每张全部 102 bins（含 flow）、raw Entries 和 sumw2 相等 |
| 正式 CTest | 1/1 通过；安装前复核再次通过（0.58 s） |
| 原目录组件构建 | `FaserActsKalmanFilter` 构建完成；测试 target 单独调用构建通过 |
| 安装后默认配置 | `run/lib/libFaserActsKalmanFilter.so` 加载成功，Gaudi property 为 `True` |
| 研究仓库回归 | 最后一次 43 项 tests 全通过 |

详细结果为 [validation.json](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_default/20261113_mr/smoke/validation.json)。`writer_snapshot.npz` 保留选定数组；`all_histogram_counts.json` 保留三模式全部 P/F/S 直方图的 `sumw/sumw2`、flow 和 raw Entries。`bin_counts.json` 与 `figure_index.json` 保存归一化叠图 counts 和图表索引。

![New default versus explicit isolated ON](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_default/20261113_mr/smoke/figures/default_vs_explicit_on.png)

![Native source-style P/F/S: legacy OFF versus new default](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_default/20261113_mr/smoke/figures/default_vs_legacy_PFS.png)

完整 [3 页 PDF 图集](/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_default/20261113_mr/smoke/figures/smoothing_default_validation.pdf)，逐图 PNG/PDF 保存在同目录。主图已目视核对，图中 Entries 与冻结 counts 相符。

验证脚本适配原宏时先后遇到输入名必须包含 `3station`、保存时必须保留 `h_component_station_layer` 键名两个问题；均已修正，失败日志保存在 `R/smoke_failed_macro_name.log`、`R/smoke_failed_histogram_key.log`。这是回归脚本的文件/直方图命名错误，没有更改 CKF2 修复。最终三组检查完全通过，临时 ROOT/full states 已清理。

复现真实重建检查：先加载环境，并把 `R/software_clean/lib` 和 `R/software_clean/python` 放在动态库/Python 搜索路径首位，同时设置 `GEOROTATION_LIBRARY_DIR`，然后执行：

```bash
python3 scripts/validate_smoothing_default.py \
  --reference-task /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/afterGeo_smoke_tasks.json \
  --software "$R/software_clean" --output "$R/smoke_reproduction" \
  --previous-on /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/smoke/afterGeo
```

本次直方图用于软件回归，采用原 GitHub cuts，不额外加 purity gate，以便逐 bin 核对上一轮冻结 ON smoke。用户确认的后续物理画图基准仍为 **isolated smoothing + strict final purity >0.7 + native truth**，其全样本结果见 [good-seed 比较报告](2026-11-13_stage04_good_seed_comparison.md)。

## 阶段性结论

真实 ACTS state-container 回归支持共享祖先 smoothed 状态覆盖的定位；深拷贝后重新 smoothing 同时恢复候选各自的参数与 covariance，并消除处理顺序对保存结果的影响。真实 Calypso 小样本证明新默认精确复现此前显式 isolated ON 的输出；本次集成没有再引入额外 seed、purity 或 variance 设置变化。原日常 `run/` 的库及 metadata 均已安装更新，后续不传 smoothing 开关即采用隔离方法。

## 合入目标开发分支

源分支为 `xcheng/calypso:fix/ckf2-isolated-smoothing`；目标为 `sazhang/calypso:geoRotation_rebasedOfficialMaster`。源分支已发布，包含一个功能提交和一个仅文档的验证提交。不要把原研究目录中的全部未提交改动一起推送。

已创建 [Draft MR !1](https://gitlab.cern.ch/sazhang/calypso/-/merge_requests/1)。通过 GitLab API 核对 source project ID `241372`、target project ID `201326`、两侧分支、HEAD `25147d80` 和上述五个 changed files；GitLab 返回 `can_be_merged`，当前阻止合并的状态是 `draft_status`。这是 GitLab 的合并检查结果，不代表目标项目已经跑过 CI。API 核对结果保存为 `R/merge_request.json`，发布日志为 `R/push.log`、`R/merge_request_creation.log`。

当前 shell 的 SSH 认证不可用；CERN Kerberos HTTPS 已验证能读取个人 fork。发布命令为：

```bash
cd /eos/home-x/xcheng/FASER/geoRotation_calypso_smoothing_fix
git -c remote.origin.pushurl=https://:@gitlab.cern.ch:8443/xcheng/calypso.git \
  push --set-upstream origin fix/ckf2-isolated-smoothing
```

已有 MR 中的 Source 为 `xcheng/calypso` 的 `fix/ckf2-isolated-smoothing`，Target 为 `sazhang/calypso` 的 `geoRotation_rebasedOfficialMaster`；Changes 已核对只有上述五个文件。英文标题、正文位于 `R/merge_request_title.txt` 和 `R/merge_request_description.md`。开发者 review 和所需 CI 检查通过后，点击 **Mark as ready**，再由有目标分支权限的维护者点击 **Merge**。本次停在 Draft MR，实际合并由维护者执行。

如果目标分支在 review 期间更新，仍在干净修复目录执行 `git fetch sazhang geoRotation_rebasedOfficialMaster`，再 `git rebase FETCH_HEAD`，解决冲突并重跑 CTest/小样本回归后更新源分支。维护者也可以单独 cherry-pick `1c06c7b65135a88e5ebd864050b10a638b152c98`；无需合并研究工作目录。

## 下一步计划

由目标分支维护者审阅、运行其环境中的 package CTest，并将修复合入目标开发分支。后续动量图沿用用户确认的 baseline；若升级 ACTS 或改变 CKF 状态容器，重跑共享祖先回归和同一样本的默认/显式 ON 数值检查。
