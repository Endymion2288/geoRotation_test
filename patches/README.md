# Calypso 外部修改与构建

基底为 `/eos/home-x/xcheng/FASER/geoRotation_calypso` 的 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`，Athena 24.0.41 / x86_64-el9-gcc13-opt。本仓库不复制整套 Calypso。

补丁保留外部源码上下文和删除行的原有空白，`.gitattributes` 为这些文件关闭行尾空白检查；补丁中的新增源码行另行确认没有行尾空白。

S06 isolated-smoothing批次的增量与基底：

- `0021_charge_diagnostics_after.patch`：相对旋转基底`88c457c1`，默认isolated smoothing ON；只读seed/candidate、P/F/S、native cached reference、conversion/refit和pruning记录。真实seed索引另存，不修正原生trackGroup，保留物理基线。
- `0022_charge_diagnostics_before.patch`：相对master`1eb885ae`的对应诊断与默认ON。0021/0022是两套替代基底补丁，不能串联。
- `0023_charge_observer_qop.patch`：相对已有S05 observer的q/p及variance字段增量；不包含此前transport probe。与tracking库组成四库overlay，6组20-event ON/OFF的58数组精确一致，见`stage06_charge/20261113_isolated_v1/charge_smoke/validation.json`。
- `0024_matched_seed_prior_after.patch`：在0021之后应用，仅修改CircleFitTrackSeedTool；默认关闭的`MatchFaserSeedCovariance`以实际两套seed surface basis和方向Jacobian映射完整6×6 covariance。`SeedPriorAuditFile`保存实际输入C、J和输出C供验证。基底prior数值必须使用beforeGeo配置；它只匹配seed先验，不代表全链纯旋转已闭合。实际编译与后续数值验证状态见[S06报告](../workbook/2026-11-13_stage06_charge_baseline.md)。
- `0025_charge_hit_and_polarity_controls.patch`：在S04 `software_controls_v3/tracking.patch`之后应用，保留native组合CKF，默认关闭地过滤无指定primary SDO支持的source links；seed mean/C不变。另有seed q/p及完整C交叉项反号、CKF/refit/extrapolator的field/gradient一起反号。与0021不是串接关系。
- `0026_charge_conjugated_truth_transport.patch`：相对0023后的observer增量，只将原输入production truth用于model transport的charge反号；reco charge已由CKF反号，不重复乘。v2的三DSID smoke默认58字段精确、q/B路径和完整C闭合，见`stage06_charge/20261113_isolated_v1/hitcontrol_smoke_v2/validation.json`。
- `0027_charge_refit_chart_before.patch`：0022之后的refit五文件增量，来自S05已测量的native regular update/smoother chart；默认关闭，不包含CKF替换。实际beforeGeo build及两电荷20-event默认58字段、全部CKF行闭合，见`refit_control_smoke_v1/validation.json`。
- `0028-nominal-alignment-store-control.patch`：真实`FaserSCT_AlignCondAlg`默认关闭的nominal store控制与sensor transform只读审计；nominal模式不读alignment folder、不施加conditions displacement，仍生成实际SCT conditions store供detector elements和ACTS消费。相对`88c457c1`、仅conditions package；该package在before/after基底源码相同。实际11个translation units编译成功，6组20事件默认58字段精确，192个sensor实际/default矩阵差为0。匹配G4/digitization/reco链通过，default/nominal六组58字段精确相同，见[固定seed与alignment报告](../workbook/2026-11-14_stage06_fixed_seed_alignment.md)。

S06全量主测量及历史Alma配对的cuts、分母、软件快照与未完成实验同见该报告；不将探索控制视为正式生产修复。

- `0001-sct-local-truth.patch`：环境变量 `FASER_SCT_LOCAL_TRUTH_FILE` 启用 G4 SCT pre/post momentum CSV。使用单进程、单 G4 worker；退出134已在S04定位修复，跨事件关联需同时应用0005。
- `0002-tracking-audit.patch`：event/truth/seed/candidate/CKF/conversion/refit/final JSONL、完整 states 和实际 transport closure；修正 seed lineage 的索引。`AuditFile` 为空时不写记录。
- `CorrectCovariance` / `CorrectInterfaces`：实验开关，分别转换 plane↔curvilinear covariance，并将 refit 初始 plane 放在相同物理点；默认 false。
- `ResetEventMaps`：实验开关，在 event 开始清空 seed 静态 maps；默认 false。
- `0003-momentum-controls.patch`：在 `0002` 后应用；增加独立 smoother 快照、seed 输入审计、实际材料/粒子假设/step/C 控制、确定性 truth mean 注入及纯 hits 的非组合 KF。物理控制默认关闭，energy loss/scattering 保持原 true。S04 对全部名义事件核对诊断控制是否保持 reconstruction 输出一致。
- `0004-diagnostic-generation.patch`：ParticleGun HepMC3 显式使用 MeV/mm；Pythia8Decayer 隐藏静态链接库符号，避免其 G4ios finalizer 经符号抢占调用另一库的 finalizer，导致第二次 delete。与 `0001` 配合用于补充 G4；不修改历史 RDO。
- `0005-sensor-event-identity.patch`：在 `0001` 后应用，将 sensor CSV 的主键改为 Athena event number/source entry，同时保留 G4 event ID。当前 serial G4 调用的 G4 event ID 恒为1，不能用于跨事件真值关联。S04 用相同随机种子重新生成2400事件的sidecar，逐步验证物理数值与原模拟一致后再关联。
- `0006_slope_diagnostics.patch`：实际 Calypso 的只读 final/segment/truth observer 和 InputStates transport 组件；raw/physical frame 并存，所有 fixed-plane propagation 用真实 field、geometry 和 extrapolator。独立编译为诊断 plugin，与原生 before/after tracking 库组合，不替换 tracking。
- `0007_rotated_ntuple_angles.patch`：NtupleDumper 的 `RotatedPhysicalAngles` opt-in，默认 false；保存正确物理 VetoNu angles。20+500 events 的真实 writer 配对和逐字未改 GitHub 宏核对见 S05。
- `0008_master_ntuple_waveform_accessor.patch`：仅 master NtupleDumper 的 waveform accessor 编译兼容；保持 before tracking 源码与原生 frame。
- `0009_truth_cluster_position.patch`：真实 TrackerClusterMakerTool 的默认关闭 `UseTruthClusterPosition`，仅替换 local measured strip coordinate，保留原误差/宽度/噪声及无 truth fallback。依赖同 release 的 HitManagement/TrackerSimEvent；构建入口 `scripts/build_cluster_control.py`。是诊断 ideal centroid 控制，不是关闭完整 digitization 或正式修复。

S05 的 0006–0009 各自以其记录基底及适用 package 应用，不能把 before 专用 0008 当作 after tracking patch 串接。`stage05_slopes/patch_validation_v1/validation.json` 核对四份 patch 与实际源码一致，软件来源、观察闭合及测量限制见 [S05 报告](../workbook/2026-11-12_stage05_slope_structure.md)。

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

问题对应的Calypso源码路径、默认材料入口和实际作业来源见[源码与结果来源报告](../workbook/2026-11-12_stage04_calypso_provenance.md)。当前`faser_reco.py`第196行的默认材料路径尚未修改，已经验证的是显式flag替换；诊断及修复继续在同一Calypso源码树和匹配构建中进行。

`0002 → 0003` 与 `0004` 已在仅含基底源码的小目录顺序应用，检查所得 12 个 tracking 文件、2 个生成文件与实际源码逐字节一致。记录在 `stage04_momentum/patch_validation/validation.json`；`0003` 是增量 patch，不应与冻结软件的累计 `tracking.patch` 重复应用。

最终 `0001 → 0002 → 0003 → 0004 → 0005` 顺序应用后，全部15个文件与实际Calypso源码逐字节一致，见 `stage04_momentum/patch_validation_final/validation.json`。controls_v5冻结生成库与sensor event metadata；使用同随机种子的2400事件、81153个step验证物理数值未变。

`scripts/freeze_software.py` 将小型 worker 配置、精确 tracking patch、编译库和 provenance 复制到独立 results 目录。Condor 读取冻结副本并优先加载其库，运行期间不依赖正在编辑的仓库代码。原组件库另存于 `stage02_closure/original_lib/`；未使用新属性的原版对照通过 `LD_LIBRARY_PATH` 加载它。

- `0018_writer_seed_lineage.patch`：只读诊断，为 `RootTrajectoryStatesWriterTool` 增加 `seed_nr` 分支，从实际 CKF `trackGroup` 读取 seed 编号；不含该列的 refit writer 输出 `UINT32_MAX`。同一补丁适用于 before/after 原生基底，不更改状态、选择或拟合。与 seed writer 的 `seedIndex` 按 input file + source entry 联结，避免重复 event number 或同 barcode 多 seed 造成错误归属。
- `0019_native_ckf_seed_index.patch`：修正两个原生基底 CKF 的来源标签 `nSeed-1`（`nSeed` 恒为0）为真实 `iseed`。该标签在本次禁用 alignment writer 的配置中不参与 tracking；与 0018 配合导出精确 seed 归属。主诊断源码树的 audit patch 已使用 `iseed`，不重复应用 0019。
