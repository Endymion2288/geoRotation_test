# S04补充：材料坐标问题的Calypso源码位置与结果来源

日期：2026-11-12。Calypso根目录为`/eos/home-x/xcheng/FASER/geoRotation_calypso`，基底commit为`88c457c1e4bf4154e67b96114ee8e7020d9b1a97`。结果根目录为`/eos/home-x/xcheng/FASER/geoRotation_test_results`，下文记为`$O`。

## 阶段目标

定位材料图坐标与旋转几何不一致的实际Calypso调用链；核对已交付动量图和独立配对的重建来源；执行用户“全部结果基于Calypso，不使用另写的独立或简化tracking版本”的要求。

## 结论与源码路径

已验证的问题位于**旋转后的tracking geometry与仍按原坐标查表的material map之间的接口**。源码已经旋转sensor和layer，但重建入口仍选取原材料JSON，geometry service直接加载JSON中的binning transform，没有根据`geometryRotation`同步旋转它。严格旋转同一JSON后，实际Calypso重建的动量负尾明显减少。这里不涉及将动量矢量旋转后改变模长；同一输入的seed动量不变。

以下行号对应本次核对时的实际源码。后续编辑可能改变行号，函数名和文件路径是定位入口。

| 环节 | Calypso详细路径与行号 | 实际行为 |
| --- | --- | --- |
| 旋转开关 | [Tracking/Acts/FaserActsGeometry/FaserActsGeometry/FaserActsGeometryContext.h](../../geoRotation_calypso/Tracking/Acts/FaserActsGeometry/FaserActsGeometry/FaserActsGeometryContext.h)，第21行 | `geometryRotation{true}`，默认使用旋转后的ACTS frame |
| Sensor旋转 | [Tracking/Acts/FaserActsGeometry/src/FaserActsDetectorElement.cxx](../../geoRotation_calypso/Tracking/Acts/FaserActsGeometry/src/FaserActsDetectorElement.cxx)，第89–108、130–151行 | `transform()`和`storeTransform()`分别为construction与alignment store应用`R`；坐标为`(x_A,y_A,z_A)=(z_F,y_F,-x_F)` |
| Layer导航分箱 | [Tracking/Acts/FaserActsGeometry/src/CuboidVolumeBuilder.cxx](../../geoRotation_calypso/Tracking/Acts/FaserActsGeometry/src/CuboidVolumeBuilder.cxx)，`buildLayer()`第60–86行 | 沿束流的layer导航从`binZ`切到`binX`。它与外部材料JSON中的局部分箱transform是两件事 |
| 默认材料选择 | [Control/CalypsoExample/Reconstruction/scripts/faser_reco.py](../../geoRotation_calypso/Control/CalypsoExample/Reconstruction/scripts/faser_reco.py)，第192–196行 | `useCKF`时仍指定`/cvmfs/faser.cern.ch/repo/sw/database/DBRelease/current/acts/material-maps-alma9.json` |
| Python配置传递 | [Tracking/Acts/FaserActsGeometry/python/ActsGeometryConfig.py](../../geoRotation_calypso/Tracking/Acts/FaserActsGeometry/python/ActsGeometryConfig.py)，`ActsTrackingGeometrySvcCfg()`第8–17行 | `flags.TrackingGeometry.MaterialSource`被交给`MaterialMapInputFile`，启用`UseMaterialMap` |
| JSON实际加载 | [Tracking/Acts/FaserActsGeometry/src/FaserActsTrackingGeometrySvc.cxx](../../geoRotation_calypso/Tracking/Acts/FaserActsGeometry/src/FaserActsTrackingGeometrySvc.cxx)，第88–101、118–128行 | 构造`Acts::JsonMaterialDecorator`并交给`TrackingGeometryBuilder`。这段没有将原JSON的`binUtility.transform`与旋转几何同步变换 |
| 几何服务接口 | [Tracking/Acts/FaserActsGeometry/src/FaserActsTrackingGeometryTool.cxx](../../geoRotation_calypso/Tracking/Acts/FaserActsGeometry/src/FaserActsTrackingGeometryTool.cxx)，第37–40行 | `trackingGeometry()`返回上述service建立的实际ACTS几何 |
| CKF入口 | [Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx](../../geoRotation_calypso/Tracking/Acts/FaserActsKalmanFilter/src/CKF2.cxx)，第102–111、325–342行 | 使用该geometry建立finder/fitter；普通路径调用`m_finder`，固定truth hits的诊断分支调用同组件内`m_fitter` |
| Prediction/update/smoothing | [Tracking/Acts/FaserActsKalmanFilter/src/TrackFinderFunction.cxx](../../geoRotation_calypso/Tracking/Acts/FaserActsKalmanFilter/src/TrackFinderFunction.cxx)，第15–44、49–81行 | 实际`Acts::CombinatorialKalmanFilter`、`GainMatrixUpdater`和`GainMatrixSmoother`，材料散射及能损选项传入ACTS；不是仓库内重写的KF |
| Refit入口 | [Tracking/Acts/FaserActsKalmanFilter/src/KalmanFitterTool.cxx](../../geoRotation_calypso/Tracking/Acts/FaserActsKalmanFilter/src/KalmanFitterTool.cxx)，第39行；[TrackFitterFunction.cxx](../../geoRotation_calypso/Tracking/Acts/FaserActsKalmanFilter/src/TrackFitterFunction.cxx)，第27–64行 | Refit也使用同一tracking geometry，实际执行ACTS KalmanFitter及其update/smoother |

Calypso所用ACTS依赖也有明确路径，当前头文件根目录为：

```text
/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/include/Acts/
```

其中`Utilities/BinUtility.hpp`第153–155行用`m_itransform * position`将global position送入材料分箱。Calypso加载的JSON若保留F-frame transform，而传入位置已是A frame，就会查到不同bin。`Propagator/MaterialInteractor.hpp`第88–106行计算材料相互作用、transport covariance并更新state；随后实际CKF的measurement gain/smoothing会受到影响。这些是Calypso链接的ACTS实现，不是独立toy。

## 已验证的修复落点

本次采用的修复是**配置同一Calypso加载坐标一致的材料JSON**，没有把独立程序的拟合结果写回Calypso。修复入口就是`faser_reco.py`第196行对应的`flags.TrackingGeometry.MaterialSource`；诊断入口[tracking_runner.py](../configs/tracking_runner.py)第80行设置相同flag，再从Calypso导入`FaserGeometryCfg`、clusterization、space points、`SegmentFitAlgCfg`、`GhostBustersCfg`和`CKF2Cfg`。

| 配置 | 实际材料文件 |
| --- | --- |
| current | `$O/stage02_closure/material-current-source.json`，原current材料的冻结副本 |
| strict | `$O/stage02_closure/material-current-rotated.json`，同一份材料的严格坐标旋转 |

[rotate_material.py](../scripts/rotate_material.py)只预处理JSON数据：对每个`binUtility`的local→global旋转/平移执行`U_A=R U_F`、`t_A=R t_F`，保持geometry IDs、材料bin data及局部分箱边界不变。其自身不传播或拟合track。22个surface entries、每surface 4096点的局部坐标一致性检查最大差为0 mm，见[material_rotation.json](../../geoRotation_test_results/stage02_closure/material_rotation.json)。

例如一个approach surface的材料transform原平移为`(0,0,-1885.1075)` mm，strict中为`(-1885.1075,0,0)` mm，并具有相应`R`旋转；bin data不变。

`geoRotation_calypso/Control/CalypsoExample/Reconstruction/scripts/faser_reco.py`的默认材料路径本次仍未改写。现有修复结果来自显式材料flag；直接使用该脚本默认值不会自动获得strict材料修正。服务端自动转换或默认材料发布尚未实现，不能把已验证的配置替换描述成已经合入默认软件的修复。

## 工作内容与实际运行来源

新增[audit_calypso_provenance.py](../scripts/audit_calypso_provenance.py)，读取已有任务、`run.json`、软件manifest及压缩日志，不运行新的reconstruction。实际核对如下：

| 图或统计 | 实际重建来源 | 本次逐作业核对 |
| --- | --- | --- |
| 独立样本final分布与34.151%→10.543%等主配对数值 | `$O/stage03_baseline/software_v2/lib/libFaserActsKalmanFilter.so` | 100个作业，current/strict各50000事件，50个配对ranges；全部成功 |
| Seed、九层P/F/S、good-seed、CKF/final诊断图 | `$O/stage04_momentum/software_controls_v3/lib/libFaserActsKalmanFilter.so` | 100个作业，current/strict各50000事件，50个配对ranges；全部成功 |

两份`.so`的来源均记录为`/eos/home-x/xcheng/FASER/geoRotation_calypso/run/lib/libFaserActsKalmanFilter.so`，使用同一Calypso基底及各自`tracking.patch`。实际环境为Athena24.0.41、EL9/GCC13；作业日志报告ACTS32.0.2。冻结副本用于防止长批次运行期间源码或构建更新，并非另维护一套tracking实现。软件manifest与patch分别见：

- [S03冻结软件manifest](../../geoRotation_test_results/stage03_baseline/software_v2/software.json)、[tracking.patch](../../geoRotation_test_results/stage03_baseline/software_v2/tracking.patch)。
- [S04逐层冻结软件manifest](../../geoRotation_test_results/stage04_momentum/software_controls_v3/software.json)、[tracking.patch](../../geoRotation_test_results/stage04_momentum/software_controls_v3/tracking.patch)。

200个作业的task/run material路径、源输入、event range、Calypso安装位置、library override及实际Calypso initialization/finalization日志全部一致。同组配对只改变材料路径，runner和其他options相同。检查记录见[reconstruction_audit.json](../../geoRotation_test_results/stage04_momentum/calypso_provenance/reconstruction_audit.json)。历史作业没有保存进程`/proc/maps`；库来源证据是worker的`LD_LIBRARY_PATH`、`run.json`与冻结软件manifest，不将它表述成历史动态链接映射的直接测量。

来源命令：

```bash
python3 scripts/audit_calypso_provenance.py \
  --root /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum \
  --output /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum/calypso_provenance/reconstruction_audit.json
```

## 本地修改与结论边界

该Calypso工作树有本地诊断和受控实验修改，不是只有HEAD commit的干净源码。本仓库保存[0001–0005补丁及构建说明](../patches/README.md)；应用到该基底后15个文件与实际源码逐字节一致，见[patch_validation_final](../../geoRotation_test_results/stage04_momentum/patch_validation_final/validation.json)。每批使用对应冻结patch，而不是将当前所有补丁回溯当作每批已启用的配置。

逐层图使用`--states --selected-states-only --reset-event-maps --isolate-smoothing`。candidate独立smoothing在Calypso的`CKF2.cxx`内执行实际ACTS smoother，解决共享ancestor诊断slot被覆写的问题。原始名义重建与该诊断构建已对100000次事件处理逐字段配对：current的70170条final、strict的66370条final，其动量、电荷和Good95判定均一致；converted的rmax分别8832/6625条、final的rmax分别37/36条改变。因此不能把诊断构建称作逐字段完全未改动的原版。详细差异见[current](../../geoRotation_test_results/stage04_momentum/observer_integrity/current/summary.json)、[strict](../../geoRotation_test_results/stage04_momentum/observer_integrity/strict/summary.json)。

已有2400事件的truth-seed/固定hits控制确实调用Calypso的ParticleGun、G4、digitization与实际KF，但生成时仅开启`FaserSCT`与`Dipole`。它们是注明子探测器范围的Calypso诊断实验，不是完整探测器效率样本；不替代官方RDO上的独立材料配对。后续新生成与验收使用匹配的完整Calypso几何及实际组件，保留每个受控开关。

此前独立的`configs/acts_momentum_toy.cxx`使用均匀0.55 T和自建18个strip planes；虽然调用ACTS，它不属于Calypso实际重建。本次将该30000-fit结果移出S04物理证据、资源主表和T5验收，保留代码和旧产物作为历史记录。离线circle替代解和理想点只作数学诊断，也不构成Calypso性能或修复结果。已交付的31页动量图、官方RDO材料配对和真实merge结果不读取toy输出。

## 验证与下一步

事实：本次200个实际重建作业的来源核对通过，没有运行新track fit；材料配对和图集中使用的是Calypso组件输出。独立toy退出后，材料主结论仍有官方RDO配对、同seed/hits逐层记录和实际Calypso材料/scattering消融证据。已完成的物理数值不因这次来源整理改变。

未测量：T5在完整Calypso几何/生成/重建配置中的完整统计closure。此项列入S07，不能用以前的standalone toy补足。S05–S07的重建、传播和算法替换均在Calypso源码/配置中执行；本仓库保留配置、patch、matching、统计、绘图和报告。AGENTS与roadmap已同步这一要求。

本次新增脚本的实际执行、Python语法/帮助检查、Markdown本地链接和`git diff --check`结果汇总在[validation.json](../../geoRotation_test_results/stage04_momentum/calypso_provenance/validation.json)。
