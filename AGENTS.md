# Repository Guidelines

## 项目结构

- `docs/`：输入参考材料，保留原始 PDF；用户提到的 geoRotation 报告实际文件名为 `geometryRotation.pdf`。
- `roadmap/README.md`：研究总入口、阶段依赖和完成标准；`roadmap/00_materials_and_source.md` 保存材料证据和源码定位，其他文件分别规定公共实验协议与三个研究分支。
- `workbook/`：中文阶段报告，包含实际执行、数据和结论。计划写入 roadmap，测量结果写入 workbook。
- `configs/`：冻结样本、物理定义、字段字典与 Athena 入口；`scripts/`：metadata、material map 和 Condor 工具；`analysis/`：合并、统计与英文图表；`tests/`：matching、frame、分母等回归检查；`patches/`：外部 Calypso 的可复用补丁。
- Calypso 源码与构建位于仓库外。诊断配置、分析和提交脚本保存在本仓库；修改 Calypso 时记录其 commit、patch 和构建路径，避免将本仓库结果与另一套安装混用。

## 运行环境与软件

- 当前工作节点为 `lxplus-gpu`，终端可使用 GPU。交互式调试和小规模验证可直接使用 GPU，大批量数据处理必须通过 HTCondor 提交。
- CVMFS 已提供所需软件包，Calypso 等工具通常已有编译产物。优先复用与实际配置匹配的现成构建，只有源码变化时才重新构建。
- 已核对的环境为 Athena `24.0.41`、`x86_64-el9-gcc13-opt`；使用 `/eos/home-x/xcheng/FASER/geoRotation_calypso/setup_environment.sh` 和对应 `run/`。软件版本、实际 hostname、GPU 可用性和资源使用随运行记录；ACTS 重建的 CPU 作业按 CPU 吞吐量配置。
- Kubeflow 服务也可用；需要使用而用户尚未登录时，提醒用户登录。

## 数据、产物与任务目录

| 用途 | 路径 |
| --- | --- |
| 输入样本 | `/eos/experiment/faser/data0/sim/mc24/`；已确认的 DSID、RDO/PHYS 和版本见 `roadmap/00_materials_and_source.md`，执行前冻结具体文件与 event range |
| 旋转版 Calypso | `/eos/home-x/xcheng/FASER/geoRotation_calypso`；研究起点 HEAD `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` |
| 已有构建与安装 | 上述源码目录下 `build/`、`run/`；源码未变时复用 |
| Geant4 模拟源码，探测器响应细节以此为准 | 上述 Calypso 下 `Simulation/`、`Tracker/TrackerG4/FaserSCT_G4_SD/`；digitization 见 `Tracker/TrackerDigitization/FaserSCT_Digitization/`，实际样本的生产 release/conditions 另行核对 |
| 所有运行输出，包括日志、图表和模型 | `/eos/home-x/xcheng/FASER/geoRotation_test_results` |
| HTCondor 提交文件 | `~/condor/geoRotation_test/<stage或任务名>/` |

home 空间有限，提交文件应尽可能小，作业日志和数据输出写入上述 results 目录。阶段报告和代码、配置保存在仓库中。

EOS 仅保留尽量小的最终结果，避免保存 xAOD 等中间文件；合并或去重批量作业产生的 JSON 等重复输出。大型数据、checkpoints 和冗长 notebook 输出不纳入 Git。

结果采用 `<results>/stageNN_<topic>/<run_id>/`，内含 manifest、日志、紧凑统计和图表；每个 shard 使用独立目录，合并时检查重复事件键和覆盖范围。当前样本在同一 RDO 文件内也会复用 run/lumi/event；以冻结的输入文件与 entry 作为事件身份，保留生产标识，不据它误删独立事件。批量作业将临时 xAOD 和完整 states 放在 worker scratch，在同一作业内提取所需字段后清理；需复现的少量异常事件保留紧凑诊断表及输入定位。

## 开发与验证命令

以下命令对应现有软件。初始化与帮助检查：

```bash
cd /eos/home-x/xcheng/FASER/geoRotation_test
georotation_results=/eos/home-x/xcheng/FASER/geoRotation_test_results
mkdir -p "$georotation_results/manual_checks"
source /eos/home-x/xcheng/FASER/geoRotation_calypso/setup_environment.sh > "$georotation_results/manual_checks/environment.log" 2>&1
python3 /eos/home-x/xcheng/FASER/geoRotation_calypso/Control/CalypsoExample/Reconstruction/scripts/faser_reco.py --help > "$georotation_results/manual_checks/reco_help.txt" 2>&1
git diff --check
```

源码修改后才按所申请的 CPU 数重建，构建日志放 results：

```bash
bash /eos/home-x/xcheng/FASER/geoRotation_calypso/build_local.sh 8 > "$georotation_results/manual_checks/build.log" 2>&1
```

`verify_local.sh` 提供库加载、组件配置和两个 GeoModel 单测，可按需运行；它的 CTest 还会在外部 build 目录写 `Testing/`。完整重建样本的配置、计数和物理检查按 `roadmap/01_common_protocol.md` 实施。普通 `faser_reco.py` 没有 `-o` 和 diagnostics CLI；本仓库显式入口为：

```bash
python3 scripts/inspect_inputs.py --samples configs/samples.json --output "$georotation_results/manual_checks/input_metadata.json"
cd "$georotation_results/manual_checks"
python3 /eos/home-x/xcheng/FASER/geoRotation_test/configs/tracking_runner.py --dsid 100116 --events 10 --diagnostics
```

`--audit --states --closure` 需要本仓库记录的 Calypso audit patch 及匹配构建；`--correct-interfaces` 和 `--reset-event-maps` 是独立受控实验，不作为未经说明的默认修复。Python observer 的 `events.jsonl` 用于小样本交叉核对，批量用 `--audit --compact-only` 避免重复输出。回归检查在本仓库根目录运行 `python3 -m unittest discover -s tests -v`。

S04 的批量入口为 `scripts/submit_momentum.py`，接受冻结 task JSON、library 目录和唯一任务名，提交前重新选择 schedd。`--states --selected-states-only` 配合 worker 的 `audit_reduce.py` 只保留选中路径的 compact traces；`--isolate-smoothing` 会改变少量 rmax 元数据，不能预设所有字段不变。sensor CSV 必须用 patch0005 的 Athena event/source_entry，G4 event ID 在当前调用中恒为1。实际配置、50k独立配对和复现入口见 S04 workbook；大型统计也通过 Condor 执行。

提交前读取当前 schedd 广告与队列，按实际可见 schedd 进一步查询：

```bash
condor_status -schedd -af Name Machine TotalRunningJobs TotalIdleJobs MaxJobsRunning
condor_q -name <schedd> -totals
```

结合 EL9、EOS/CVMFS、凭据、资源限制及 idle/running/capacity 选择兼容且负载最低的 schedd。将查询时间、候选负载、选择依据、submit 文件和 cluster ID 写入该阶段 manifest；随后使用 `condor_submit -name <schedd> <submit_file>`。`scripts/prepare_condor.py` 只生成 ranges 与 submit 文件，实际提交及 schedd 选择必须显式执行和记录。

## 编码规范与计算效率

- Python 使用 4 空格缩进、`pathlib`、明确的函数与字段名；C++ 遵循所改 Calypso package 的现有风格，单位优先使用 ACTS/CLHEP 常量。配置项与分析阈值集中管理，避免在多份脚本中复制。
- 每个参数、表和图声明坐标系、参考 surface、单位和 stage。区分原始 FASER frame 与旋转 ACTS frame，以及 GeV/MeV、`q/p` 和 `p`、`py/pz` 和 `atan(py/pz)`。
- 保留稳定的 event/particle/seed/candidate/state 标识；truth matching 使用 hit/SDO/barcode，电荷与重建动量不参与主匹配条件。失败、未匹配、重复候选和 refit fallback 分别计数。
- 新增诊断保持 reconstruction 行为一致；先验证观察开关不改变事件和候选结果，再做单因素消融。几何旋转的 on/off 必须覆盖几何、磁场、seed、refit、输出和 covariance；不能只切换 `geometryRotation` 成员。
- 所有重建结果、物理结论和验收图均基于 `/eos/home-x/xcheng/FASER/geoRotation_calypso` 的实际组件及其匹配构建。新增诊断、truth seed、固定 hits 和消融在该源码树内实现并保存 patch；配置脚本调用 Calypso 组件，仓库分析脚本负责 matching、统计和绘图。禁止使用另写的独立 tracking、简化 detector/field 或 standalone ACTS toy 代替 Calypso 结果。冻结的 `.so` 是该 Calypso 构建的运行副本，记录来源、commit/patch、配置和库目录；实际源码的本地改动必须如实列明。
- 对坐标、单位、covariance Jacobian、候选归属和已定位缺陷增加有意义的小型回归检查。文档调整使用链接检查和 `git diff --check`；物理结论使用明确样本上的比较与区间估计。
- ROOT 分析只读取必要 branch，优先 ROOT RDataFrame 或批量数组计算，避免逐事件 Python 对象循环；保存 `sumw`、`sumw2` 和原始计数。

代码与 HTCondor 脚本均须考虑多核、多线程、I/O 和文件读写效率。尽可能并行处理，线程数与作业申请的 CPU 资源匹配；减少重复读取和小文件写入，按任务选用 C++、Rust、ROOT RDataFrame 等合适语言或软件包加速。除非确有必要，避免 SHA 校验和大段防御性编程带来的额外 I/O 与运行开销。

提交 HTCondor 作业前，检查可用 schedd 的排队与负载情况，选择兼容任务且最空闲的 schedd，记录选择依据以减少排队时间。

## 物理分析与研究方法

每一步都应推进一个有数据证据支持的明确物理结论，同时检查物理原因与程序原因。少数离群事件不足以终止分析分支；原因未知时明确记录待查问题，不将其作为结论。报告说明样本、选择条件、系统效应和统计依据，结论的确定程度须与证据相符。

需要新思路时保持发散思考，主动检索网络论文，兼顾高能物理、计算机科学和机器学习方法；记录引用及其对当前分析的可行性判断。

## 阶段报告、Commit 与 Pull Request

每个阶段完成后，先在仓库根目录 `workbook/` 写中文 Markdown 报告，再将报告与相关代码或配置一起提交；每阶段至少一次 Git commit。推荐文件名为 `YYYY-MM-DD_stageNN_topic.md`，至少包含：

- **阶段目标**：待回答的问题及验收目标。
- **工作内容**：实现、配置和分析步骤。
- **详细数据**：实际测得的数值及单位、输入与配置、执行命令、验证结果、产物路径；适用时引用数据表和图表。未测量项目明确标注“未测量”，不得用推测替代结果。
- **阶段性结论**：证据支持的结果，分别表述事实与假设。
- **下一步计划**：具体任务、所需输入和验收标准。

每阶段最终产物必须包含与参考材料对应的分布图，具体图清单与规范见 `roadmap/05_reference_figures.md`。before/after使用共同分箱和可比较的轴/归一化，标明样本、选择、frame、truth参考与entries；历史定义有误时同时标出legacy与正确物理定义。图表保存PNG、逐图PDF和多页汇报PDF，bin counts及产物索引随结果保存；workbook嵌入主图并链接完整图集。阶段完成前核对图中统计与冻结summary，并目视检查导出效果。

后续动量图优先参考用户提供的 `https://github.com/zsm2024/FASER_tracking/tree/main/test_pResolution` 原代码，保留可用的公式、分箱、颜色和九层布局，并明确必要的truth/选择差异。每张图声明输入文件数、源事件数及覆盖范围；小样本图不能标为全量。全量重画必须逐文件核对entries和event ranges，worker内完成matching及histogram归约，避免长期保存大规模逐state JSON。

后续提交统一遵循 Angular Conventional Commits：`<type>(<scope>): <subject>`，例如 `feat(gravnet): add vertex distance features`、`docs(workbook): record stage01 validation`。Pull Request 说明改动、意义、验证命令和样本选择，关联相关 issue；模型或可视化变更附相应图表。

## 沟通与文档

默认使用清晰自然的中文，用户指定其他语言时按要求使用；技术术语保留 English，以表达清楚为准。所有图表使用英文，包括标题、坐标轴、图例和标注。

先给结论或实际结果，再简洁说明改动、意义和验证。事实与假设分别表述；失败时说明具体失败步骤、已知原因和下一步，按需要解释语法与基本概念。所有文档禁止防御性表达，避免套话、泛泛警告、空泛限定语、无实质内容的章节和无关扩展。
