# tan(theta_y)：剩余峰形的来源与验证

目标是判断旋转后剩余尖峰、肩部或多峰结构是否超出 truth 谱、磁场弯曲、材料散射、测量离散性与选样能够解释的范围，并定位首次出现该结构的 reconstruction/output 环节。“看起来像多个 Gaussian”仅是观察，不作为程序问题或物理合理性的结论。

主复现样本为 100044 μ+；用 100116/100117 μ± 检验 charge/输入谱，110033/110051 检验双 track 拓扑。st1 与 VetoNu 分开分析，并按 [公共协议](01_common_protocol.md) 对齐 frame、surface、单位与 selection。

## A1：恢复正确物理角度并核对原图

1. 同时保存 raw A-frame 和 physical F-frame 的 mean/direction。由 F 的 `s_y=py/pz` 与 `theta_y=atan2(py,pz)` 独立计算；在 A 中使用 `py_A/px_A`。用已知 +z_F、±小 y-slope 的 free states 和一条真实 track 检查符号、范围与反变换。
2. 沿实际 track 检查 seed、CKF、conversion、refit、Ntuple 的 frame。当前 Ntuple 直接存 rotated momentum，VetoNu 外推角度仍计算 [1]/[2]；先证明具体运行走到该路径，再用正确物理定义与旧 branch 分别画图。
3. 明确 `Track_py0/pz0` 是 track 参数序列上游 state 的量，不保证等于报告标题的固定 st1 plane；将所有配置传播/关联到同一物理 st1 plane，VetoNu 则使用明确且与样本 geometry 一致的 plane。检查 fixed −3112 mm 与具体 detector version 的位置关系。
4. 每张图保存精确表达式、branch 单位、binning、entry 数、truth 母样本、成功外推条件及 collection；历史 `atan(py0/pz0)` 与历史 pz 切选另作原样复现。
5. 记录失败占位值 999999、NaN、zero direction、p_parallel≈0、分母符号与所有 validity bits；按源数据计数，不能先裁剪画图区间后将剩余样本当总体。

**验收**：同一 mean 用 F/A 两种公式得到一致的物理角度；同一 physical surface 的 Ntuple 与重算量闭合。若仅 writer 修复改变图，保留原 reconstruction 结果并量化此前失真的范围，随后继续检查真实剩余结构。

## A2：量化“尖峰/多峰”，排除绘图和谱混合

| 检查 | 方法与输出 | 判别意义 |
| --- | --- | --- |
| 分箱稳定性 | 原 bin 宽、1/2 与 2 倍宽、半 bin 平移；ECDF、连续 KDE 的 bandwidth 扫描 | 仅依赖 histogram edge 的尖峰属于绘图/离散效应线索 |
| 精确值与离散值 | double 原值、zero/near-zero 占比、相邻 unique spacing；查 rounding、clamp、float cast 和初始化默认值 | 区分真正窄分布与一批相同数值/未更新参数 |
| 选择迁移 | 全谱；物理 p_reco 50–250 GeV；同样 p_truth 窗口；历史 pz 窗口；A′ 50–100 GeV | 区分 momentum bias、frame 错误和 reconstructed-p 选样制造的结构 |
| truth 与 residual | 同面 truth 的 s_y、reco 的 s_y、Δs_y 和其 pull；truth/reco 二维图 | 原始 truth/磁弯结构不能直接归为 reconstruction peak |
| 多峰统计 | 对冻结条件样本作 Dip test/局部过量检验，按 event bootstrap；扫描点与窗口记录 | 描述结构稳定性；在预留事件上确认，不由一个 p-value 判定非物理 |
| 峰参数 | 峰中心/宽度/面积及区间、峰区与 sidebands 事件列表 | 便于比较 stage、charge、momentum 与 detector 条件 |

峰窗口按 pilot 的独立参考或物理 angular resolution 规定，在后续 before/after 与 fix 比较中保持相同。多 Gaussian 拟合可辅助分量描述，但分量数和拟合优度不能证明 detector 中存在对应机制。

## A3：对条件分布逐项分解

每次优先按一个物理因素分组，再分析必要的二因素交互；低统计分组给区间或合并相近条件。

- **kinematics**：p_truth、charge、true slope、生产位置、st1 真值位置/方向、A′ decay z/opening angle/两 track separation。FLUKA 谱与接受度本身可使边缘分布偏离单 Gaussian。
- **材料与磁场**：同一轨迹的累计 X/X0、路径、磁场积分、局部 B、energy loss；scattering 的宽度随 p 改变，不同能量与路径的混合可能形成窄核加宽肩。
- **measurement/geometry**：station、module/eta/phi/side、stereo pattern、hit 数、lever arm、strip cluster width、缺层/holes/outliers、cross-module hits；检查峰中心是否与 module 边界或局部 basis 相关。
- **algorithm**：seed station/fallback、segment direction、circle 条件数、seed covariance、forward/backward、IFT mask、candidate selection、refit success/fallback 和 merge source。
- **问题交集**：momentum-tail/normal、correct/wrong charge、high/low truth purity、shared hits 与重复候选。统计各分组对总峰区的占比，不只挑选漂亮事件。

对“多个子群叠加”形成可以验证的模型：

```text
f_reco(s) = sum_c integral
            f_truth(s_t, p, c) * efficiency(s_t, p, c)
            * response(s | s_t, p, c) ds_t dp
```

c 为预先定义的 detector/algorithm 条件；response 从独立控制样本或 matched simulation 中得到，允许非 Gaussian。用单个条件的 residual 学到 response，再按实测权重重构整体结构，验证峰位置、面积和宽度。用同一多峰图拟合出的 Gaussian mixture 再解释该图不算独立证据。

**判别**：条件样本 residual 平滑且经验 response 能在预留事件重建总结构，支持谱/分辨率混合；固定 truth、p 和 detector 条件后仍有异常离散/偏置，进入 A4 的逐 stage 试验。

## A4：寻找结构首次形成的环节

对峰区事件和在 p_truth/position/charge/hit pattern 匹配的 sideband/control 事件，保存：

```text
truth -> SiHit/digitization/cluster -> segment -> seed
      -> CKF predicted / filtered / smoothed
      -> ACTS-to-Trk -> refit -> Ntuple -> st1/VetoNu propagation
```

同一 surface 比较各 stage 的 s_y 及 Δs_y；逐层误差与 M 分支共享 r/S/K/F/Q/C。若峰只在拟合后出现，查看导致方向更新的 measurement 和 covariance；若只在 VetoNu 出现，查外推的 B、material、surface 和 frame。改变状态顺序或 reference surface 造成的图形变化要与物理传播效应分开。

受控实验按证据排序：

1. **输出定义**：只修正 frame/atan/单位/无效值处理，重建候选保持不变，量化对历史分布的影响。
2. **方向数值稳定性**：近 beam 方向扫过零、小正负 slopes、不同 p 与两电荷；比较 A/F 下 reconstructed free directions。检查 θ≈0 时的 Jacobian、φ wrap、normalizePhiTheta、κ–angle correlations 与 covariance。
3. **measurement 离散性**：验证实际 strip pitch、stereo angle 和 lever arm。估算 `Δs≈pitch/lever_arm` 的间距，比较峰间距随 lever arm、stereo side 与半 strip 初始位置移动的变化；再用 digitization on/off 的 uniform-direction gun 验证。
4. **材料和场**：使用匹配 toy/simulation 做无材料、散射、能损逐项实验；固定 p_truth、position/方向与 local hit 序列比较。对实际 RDO 的 reconstruction-only 开关只作敏感性研究。
5. **选择与 refit**：固定 candidate，比较 seed C、measurement C、smoothing、refit、outlier 和 forward/backward/merge 的贡献；峰若对应默认值/fallback，给其覆盖率与触发原因。
6. **几何和 conditions**：检查 construction/alignment transform、module local basis、material bins 与 field gradients，比较同一物理 hit/residual map，避免通过不匹配 production geometry 的配置“平滑”曲线。

需要新 simulation 时，在同一几何条件用窄 p、窄 slope、均匀 transverse position 与对称 slopes 的 μ±/e± gun 作为独立响应测量；批量任务走 HTCondor。峰型不符合单 Gaussian 不作为要求更换模型或 smoothing 参数的理由。

## A5：验收与最终表述

输出每个主要结构的来源表：stage、physical surface、峰位/宽度/占比/区间、对应 truth/detector/algorithm 条件、最小消融/修复效果，以及未解释余量。图表至少包括 raw 对重算角度、truth/reco/residual、条件分布重组、s_y 对 p_truth/position/module、逐 stage peak fraction。

程序问题须由固定输入的最小修改和独立样本验证；物理混合解释须量化重现结构；如果在目标精度内未见额外尖峰，给可排除的额外峰占比/宽度范围。同步报告 reconstruction efficiency、momentum tail、charge misID 和外推失败率，确认改善来自参数解释或算法修复，而非难事件丢失。
