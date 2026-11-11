# 阶段 00：参考材料研究、源码核查与实验路线图

日期：2026-11-11。已完成三份报告的材料研究、当前旋转版 Calypso 的关键源码核查、样本路径盘点、现有构建加载验证和研究路线图，补全了 `AGENTS.md`。三个物理问题的当前发生率和根因未测量；后续按 [路线图](../roadmap/README.md) 执行。

## 阶段目标

将动量低估、角度分布结构、μ± charge misID 三个问题转化为可执行且能检验因果关系的实验，明确已有证据、具体软件环节、输入样本、统计口径、诊断字段与验收条件。验收要求为路线图覆盖每个可能改变动量/charge 的操作与候选选择，并给出独立的角度定义和峰形验证方案。

## 工作内容

1. 读取三份 PDF 的全部文本内容，渲染并人工查看 problem 全部页及 geometry/Alma 的关键示意图、直方图和电荷表。problem 为图像型 PDF，pdftotext 无法提供正文，改用 pdftoppm 查看。
2. 核查 geometry rotation 的源代码范围：geometry/context/alignment/binning、field/gradient、segment/seed、CKF、ACTS↔Trk、Kalman refit、truth/state writers、Ntuple 与外推；记录 [12 项优先检查点](../roadmap/00_materials_and_source.md)。
3. 检查 EOS 中 5 类样本的 PHYS/RDO 路径和部分文件名，明确 100116 为 μ+、100117 为 μ−，区分 historical s0012、s0013 与 A′ master_dev。
4. 在现成安装上加载 4 个库，配置 4 个 Athena 组件并执行 reconstruction CLI 帮助检查；未重新编译 Calypso，未运行事件重建或提交 HTCondor。
5. 在线核对 ACTS、Kalman track fitting 与多峰检验三项方法参考的题名/摘要或书目信息，将可行用途写入材料文档。
6. 写入 roadmap 总览、公共协议和三个研究分支；补全 AGENTS 项目结构、软件/样本/Geant4 路径、有效命令、数据规范和计算要求；增加 README 导航及运行产物忽略规则。原始三份 PDF 作为研究依据纳入本阶段提交。

## 详细数据与验证

### 材料与历史数值

| 项目 | 实际核查结果 |
| --- | --- |
| problem.pdf | 4 页；重点第 3、4 页 |
| geometryRotation.pdf | 24 页；用户称 geoRotation 的报告实际文件名为此 |
| Alma 9 Validation Overview.pdf | 27 页 |
| 总页数 | 55 页 |
| report good-seed 条件 | `abs(δp)<0.1`、同电荷、`abs(δtheta)<0.0001`；最后一项单位/定义待旧配置核对 |
| 历史 μ+ misID（old/new） | 0.101% / 0.197% |
| 历史 μ− misID（old/new） | 0.217% / 0.359% |
| 历史新 tracking μ−−μ+ 差 | 0.162 个百分点，由表值相减；统计误差未测量 |
| 当前版本 tail/peak/misID/efficiency | 未测量 |
| 各输入文件 event 数、schema、生产 geometry/conditions、局部 p truth | 未测量 |

报告中的动量/角度图未给足原始计数、配对和误差；完整转录和参考页见 [材料证据](../roadmap/00_materials_and_source.md)。

### 软件、硬件与输入

| 项目 | 实际结果 |
| --- | --- |
| Calypso HEAD | `88c457c1e4bf4154e67b96114ee8e7020d9b1a97` |
| 两次旋转修改前的父版本 | `c95170f28813c273c47f8a521aa78da50871a204` |
| faser-common | `1770d15df9311f6535c09909edf9ad64a5779e32` |
| Athena / platform | 24.0.41 / x86_64-el9-gcc13-opt |
| 加载安装 | `/eos/home-x/xcheng/FASER/geoRotation_calypso/run` |
| hostname | `lxplus901.cern.ch` |
| GPU | NVIDIA Tesla T4，nvidia-smi 报告总显存 15360 MiB；本次未使用 GPU 计算 |
| 已确认样本 | 100044、100116、100117 的 PHYS/RDO 目录；110033、110051 的 RDO 目录 |
| 任务资源与吞吐 | 未测量；未提交批量任务，未进行 schedd 选择 |

### 执行命令与结果

PDF 信息与图像读取使用：

```bash
pdfinfo docs/problem.pdf
pdfinfo docs/geometryRotation.pdf
pdfinfo 'docs/Alma 9 Validation Overview.pdf'
pdftotext -layout docs/geometryRotation.pdf -
pdftoppm -f 3 -l 3 -scale-to 1700 -png -singlefile docs/problem.pdf /eos/home-x/xcheng/FASER/geoRotation_test_results/stage00_material_review/figures/problem_page3
```

对每份 PDF 文本按页提取到 results，并对关键页作渲染。当前 Poppler 21.01.0 的 `pdfinfo` 不支持 `-url`，因此没有从该命令取得报告的 Twiki 链接；历史 Good Track 全定义列为 S01 待核对项。

软件验证在 results 目录下执行，使用下列核心步骤（完整输出见对应日志）：

```bash
cd /eos/home-x/xcheng/FASER/geoRotation_test_results/stage00_material_review
source /eos/home-x/xcheng/FASER/geoRotation_calypso/setup_environment.sh > environment.log 2>&1
PYTHONDONTWRITEBYTECODE=1 python3 - <<'PY' > build_imports.log 2>&1
import ctypes
import os
from pathlib import Path
from AthenaConfiguration.ComponentFactory import CompFactory

print("Athena:", os.environ.get("AtlasVersion"))
print("Calypso:", os.environ.get("Calypso_DIR"))
for name in ("GeoModelFaserUtilities", "FaserActsGeometry", "FaserActsKalmanFilter", "NtupleDumper"):
    path = Path("/eos/home-x/xcheng/FASER/geoRotation_calypso/run/lib") / ("lib" + name + ".so")
    ctypes.CDLL(str(path), mode=os.RTLD_GLOBAL)
    print("Loaded:", path)
for name in ("FaserActsTrackingGeometrySvc", "CKF2", "KalmanFitterTool", "NtupleDumperAlg"):
    print("Configured:", CompFactory.getComp(name)("RoadmapImportCheck").getType())
PY
PYTHONDONTWRITEBYTECODE=1 python3 /eos/home-x/xcheng/FASER/geoRotation_calypso/Control/CalypsoExample/Reconstruction/scripts/faser_reco.py --help > faser_reco_help.txt 2>&1
```

结果：4/4 库加载成功、4/4 组件配置成功、`--help` exit code 0。确认 CLI 无 `-o`；普通入口默认关闭 diagnostics。库可加载与组件可配置支持复用现有构建，事件级 physics closure 尚未执行。

另执行 `git rev-parse HEAD`、`git show --stat`、定向 `rg`/源码读取、EOS 目录元数据读取和 `nvidia-smi --query-gpu=name,memory.total --format=csv,noheader`。文档交付检查包括本地链接、表格列数、代码块闭合、遗留占位项和 `git diff --check`；最终检查结果记录在 results 的 `documentation_checks.txt`。

### 产物路径

| 产物 | 路径 |
| --- | --- |
| 研究路线图 | `roadmap/README.md` 与 `00_materials_and_source.md`、`01_common_protocol.md`、`02_momentum_bias.md`、`03_slope_structure.md`、`04_charge_misidentification.md` |
| 本阶段报告与规范 | 本文件、`AGENTS.md`、`README.md`、`.gitignore` |
| 软件与硬件起点 | `/eos/home-x/xcheng/FASER/geoRotation_test_results/stage00_material_review/software_snapshot.json`、`hardware.txt` |
| 输入路径盘点 | 同目录 `input_inventory.json`、`rdo_inventory.json` |
| 软件验证 | 同目录 `environment.log`、`build_imports.log`、`faser_reco_help.txt` |
| 材料提取/视觉核对 | 同目录 `pdf_text/`、`figures/`；只保留少量核查图，其余可从原 PDF 再生成 |
| 外部参考元数据 | 同目录 `references.json` |

## 阶段性结论

**事实**：seed 及 CKF 的动量低估、旋转后角度改善但剩余结构、Alma 9 的 μ± 误判率差均有历史材料依据。现成旋转版软件已安装且关键库和组件可加载；输入样本路径存在。源码表明最终 track 还有 conversion/refit/fallback，现有 CKF states 不是全链终态；旋转提交也改变 covariance。

**事实**：源码中存在必须先核验的诊断/接口问题，包括 seed lineage 赋值、raw/physical frame、truth p 来源及 q/p 单位、plane/curvilinear covariance、两步 refit covariance inflation 和单 station 默认正电荷。它们与三个总体物理现象的因果关系尚未测量。

**待验证假设**：动量尾可能同时包含初值、hit association、材料、covariance/拟合与选择效应；角度结构可能包含输入谱/分辨率混合与数值/输出问题；μ± 差异可能包含输入谱、接受度、默认先验和后续 sign flip。路线图为每类假设配置了可区分的实验。

## 下一步计划

1. S01 冻结 5 类样本的具体文件、schema、production geometry/conditions、Good Track 和 truth 可用性；验收为 manifest、计数口径和可运行的小规模诊断配置。
2. S02 验证 frame/单位/field/material/conversion-C，修复诊断 lineage，开启 refit 观测；验收为已知 state 和固定事件闭合、观察不改变重建、连续/分片事件一致。
3. S03 在 HTCondor 执行 pilot，测吞吐、内存、输出大小、有效选中率及三问题的当前发生率，随后进入 M/A/Q 三分支受控实验。

每阶段按本仓库 Conventional Commits 约定提交配置、代码和中文报告；提交不包含 EOS 运行产物。
