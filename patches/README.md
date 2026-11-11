# Calypso 外部修改与构建

基底为 `/eos/home-x/xcheng/FASER/geoRotation_calypso` 的 `88c457c1e4bf4154e67b96114ee8e7020d9b1a97`，Athena 24.0.41 / x86_64-el9-gcc13-opt。本仓库不复制整套 Calypso。

- `0001-sct-local-truth.patch`：环境变量 `FASER_SCT_LOCAL_TRUTH_FILE` 启用 G4 SCT pre/post momentum CSV。仅单进程、单 G4 worker 验证；当前模拟退出故障见 S01 报告。
- `0002-tracking-audit.patch`：event/truth/seed/candidate/CKF/conversion/refit/final JSONL、完整 states 和实际 transport closure；修正 seed lineage 的索引。`AuditFile` 为空时不写记录。
- `CorrectCovariance` / `CorrectInterfaces`：实验开关，分别转换 plane↔curvilinear covariance，并将 refit 初始 plane 放在相同物理点；默认 false。
- `ResetEventMaps`：实验开关，在 event 开始清空 seed 静态 maps；默认 false。

保留原 `KalmanFitterTool::fit` ABI，新审计使用 `fitAudited` 显式传递 event/seed/candidate 来源。批量原版基线关闭全部物理实验开关；诊断只改变来源标签和输出。观察开关配对与 covariance 数值检查见 S02 报告，不能把实验修正的结果当作原版结果。

在匹配基底上先 `git apply --check <patch>`，再应用。初始化 `setup_environment.sh --build` 后，只构建修改的组件及其配置：

```bash
cmake --build build --target FaserActsKalmanFilter --parallel 3
cmake --build build --target FaserActsKalmanFilterConfigurables --parallel 3
cmake --build build --target CalypsoConfdb2Merge/fast
cmake -DCMAKE_INSTALL_PREFIX=/eos/home-x/xcheng/FASER/geoRotation_calypso/run -P build/Tracking/Acts/FaserActsKalmanFilter/cmake_install.cmake
cp -p build/x86_64-el9-gcc13-opt/lib/Calypso.confdb2 run/lib/
```

SD patch 对应 `FaserGeant4` target 和 `build/Simulation/FaserGeant4/cmake_install.cmake`。命令的 stdout/stderr 重定向到 results。修改源码/头文件结束后再开始编译；共享库、Python Configurables 和 confdb 必须来自同一冻结版本。

`scripts/freeze_software.py` 将小型 worker 配置、精确 tracking patch、编译库和 provenance 复制到独立 results 目录。Condor 读取冻结副本并优先加载其库，运行期间不依赖正在编辑的仓库代码。原组件库另存于 `stage02_closure/original_lib/`；未使用新属性的原版对照通过 `LD_LIBRARY_PATH` 加载它。
