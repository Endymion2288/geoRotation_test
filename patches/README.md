# Calypso 外部修改与构建

基底为 `/eos/home-x/xcheng/FASER/geoRotation_calypso` 的 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`，Athena 24.0.41 / x86_64-el9-gcc13-opt。本仓库不复制整套 Calypso。

补丁保留外部源码上下文和删除行的原有空白，`.gitattributes` 为这些文件关闭行尾空白检查；补丁中的新增源码行另行确认没有行尾空白。

- `0001-sct-local-truth.patch`：环境变量 `FASER_SCT_LOCAL_TRUTH_FILE` 启用 G4 SCT pre/post momentum CSV。使用单进程、单 G4 worker；退出134已在S04定位修复，跨事件关联需同时应用0005。
- `0002-tracking-audit.patch`：event/truth/seed/candidate/CKF/conversion/refit/final JSONL、完整 states 和实际 transport closure；修正 seed lineage 的索引。`AuditFile` 为空时不写记录。
- `CorrectCovariance` / `CorrectInterfaces`：实验开关，分别转换 plane↔curvilinear covariance，并将 refit 初始 plane 放在相同物理点；默认 false。
- `ResetEventMaps`：实验开关，在 event 开始清空 seed 静态 maps；默认 false。
- `0003-momentum-controls.patch`：在 `0002` 后应用；增加独立 smoother 快照、seed 输入审计、实际材料/粒子假设/step/C 控制、确定性 truth mean 注入及纯 hits 的非组合 KF。物理控制默认关闭，energy loss/scattering 保持原 true。S04 对全部名义事件核对诊断控制是否保持 reconstruction 输出一致。
- `0004-diagnostic-generation.patch`：ParticleGun HepMC3 显式使用 MeV/mm；Pythia8Decayer 隐藏静态链接库符号，避免其 G4ios finalizer 经符号抢占调用另一库的 finalizer，导致第二次 delete。与 `0001` 配合用于补充 G4；不修改历史 RDO。
- `0005-sensor-event-identity.patch`：在 `0001` 后应用，将 sensor CSV 的主键改为 Athena event number/source entry，同时保留 G4 event ID。当前 serial G4 调用的 G4 event ID 恒为1，不能用于跨事件真值关联。S04 用相同随机种子重新生成2400事件的sidecar，逐步验证物理数值与原模拟一致后再关联。

保留原 `KalmanFitterTool::fit` ABI，新审计使用 `fitAudited` 显式传递 event/seed/candidate 来源。批量原版基线关闭全部物理实验开关；诊断只改变来源标签和输出。观察开关配对与 covariance 数值检查见 S02 报告，不能把实验修正的结果当作原版结果。

在匹配基底上先 `git apply --check <patch>`，再应用。初始化 `setup_environment.sh --build` 后，只构建修改的组件及其配置：

```bash
cmake --build build --target FaserActsKalmanFilter --parallel 3
bash /eos/home-x/xcheng/FASER/geoRotation_test/scripts/install_tracking_metadata.sh
```

SD patch 对应 `FaserGeant4` target 和 `build/Simulation/FaserGeant4/cmake_install.cmake`。命令的 stdout/stderr 重定向到 results。修改源码/头文件结束后再开始编译；共享库、Python Configurables 和 confdb 必须来自同一冻结版本。

本轮 SD event metadata 增量编译触发 CMake 全量扫描时，已安装的 `run/src` 被原 package filter 再次选中，报重复 package。实际使用既有 `build.make` 编译一个 SD object，随后按既有 `link.txt` 重链 `FaserGeant4` 与 `FaserGeant4Lib`，最后执行该 package 的 install script；具体命令与日志见 S04。没有把失败的全量 configure 当成成功构建。

`genconf` 必须接收新 build plugin 的绝对路径。只传库 basename 时，已初始化的 `run/` 环境曾优先载入旧库，生成不含新属性的 metadata；`install_tracking_metadata.sh` 已显式绑定 build 库并检查关键属性。S04 首轮失败 smoke 的配置未进入事件分母。

`0004` 的 Pythia 修正已通过现有 `link.txt` 加 `-Wl,--exclude-libs,ALL` 重新链接并安装；实测完整命令及日志在 `stage04_momentum/g4_debug/relink_pythia_command.json`、`relink_pythia.log`。正常重建使用更新后的 CMake 定义；本轮未把未经执行的全项目 reconfigure 当作验证结果。

材料最小修复是使用与 A frame 一致的 current22 map，保持全部 geometry IDs 与 bin data。运行 [rotate_material.py](../scripts/rotate_material.py) 后，在配置中显式设置 `flags.TrackingGeometry.MaterialSource` 或通过 `tracking_runner.py --material <map>` 指定。材料 JSON 保存在 results，不进入 Git；历史 rotated12 map 不作为严格旋转 current22 的替代品。

`0002 → 0003` 与 `0004` 已在仅含基底源码的小目录顺序应用，检查所得 12 个 tracking 文件、2 个生成文件与实际源码逐字节一致。记录在 `stage04_momentum/patch_validation/validation.json`；`0003` 是增量 patch，不应与冻结软件的累计 `tracking.patch` 重复应用。

最终 `0001 → 0002 → 0003 → 0004 → 0005` 顺序应用后，全部15个文件与实际Calypso源码逐字节一致，见 `stage04_momentum/patch_validation_final/validation.json`。controls_v5冻结生成库与sensor event metadata；使用同随机种子的2400事件、81153个step验证物理数值未变。

`scripts/freeze_software.py` 将小型 worker 配置、精确 tracking patch、编译库和 provenance 复制到独立 results 目录。Condor 读取冻结副本并优先加载其库，运行期间不依赖正在编辑的仓库代码。原组件库另存于 `stage02_closure/original_lib/`；未使用新属性的原版对照通过 `LD_LIBRARY_PATH` 加载它。
