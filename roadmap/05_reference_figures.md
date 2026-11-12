# 各阶段参考图与汇报产物

每个研究阶段在数值和运算证据之外，交付与输入报告对应的分布图。参考页码均为从1开始的PDF物理页码。图用于直接比较异常的形状、位置和占比；统计结论仍使用该阶段冻结的分母、事件配对和不确定度。

## 参考页与阶段交付

| 阶段 | 参考图 | 必须交付的图 |
| --- | --- | --- |
| S04 动量 | `problem.pdf` p3；`geometryRotation.pdf` p16–18、p24 | seed相对动量残差；3 stations×3 layers的predicted/filtered/smoothed残差；good-seed条件下的九层残差；selected CKF与最终refit/fallback分布；独立样本的修正前后叠图。保留相同seed/hits的条件对照，标明与总体样本的区别 |
| S05 角度 | `problem.pdf` p4；`geometryRotation.pdf` p19–20、p23 | st1与VetoNu的`tan(theta_y)`分布，before/after并排及叠图；全谱和参考材料的动量切片；truth、seed、selected CKF、refit/final的分布或同面residual；峰区放大图和分组贡献。旧`atan(py/pz)`/pz切片另标legacy，正确frame、物理p窗口和同面truth构成主图 |
| S06 电荷 | `Alma 9 Validation Overview.pdf` p16–18、p23–25；p19–22、p27 | μ+与μ−分别的signed q/p及truth residual/pull，correct/wrong-sign尾部和事件比例；misID随truth p的分箱曲线及区间；truth/reco动量相关图和宽范围动量谱；原始谱与共同支持域/谱平衡后的对比。保留旧/new、几何旋转G0/G1、材料current/strict三种比较轴的准确标签 |
| S07 修复验收 | 对应S04–S06的全部主图 | 每项修复在预留验证输入上重新生成相同图组；提供原配置、单因素修正及最终配置的共同分箱对比，同时给效率、失败/fallback、duplicates和unmatched-to-signal。最终PDF汇编覆盖三个研究问题 |

S04历史补图入口为 `analysis/momentum_reference_figures.py`，其compact traces使用额外的fiducial/Good95/SDO选择和production truth，记录见[S04补图报告](../workbook/2026-11-12_stage04_reference_figures.md)。这些图用于相应诊断，不能作为GitHub原宏或整体几何旋转前后的复现。S05–S07按本表在阶段开始时冻结各自的图清单。

用户补充的[原绘图脚本目录](https://github.com/zsm2024/FASER_tracking/tree/main/test_pResolution)已在S04补图中按commit核对。阶段开始时必须区分源码实际启用的公式/cuts、注释中的备用配置和PDF文字；核对writer的truth/frame/单位及逐层取值方式后，再声明与参考图的一致程度。

后续同类动量主图严格使用该目录的实际公式、cuts、分箱、颜色和九层布局。`ckf_residual_perLayer.C`冻结commit为`33ca2281b961935456dfd4dc9fc06235e24b551e`，cut为`chi2Sum/NDF<20 && rMax_smt<90 && NDF>9`，100 bins、[-1.2,1.2]、raw Entries，P=green/F=red/S=blue。I/O、标签、无用越界绘图块和共同纵轴范围的必要适配保留patch；worker调用该宏的TTree::Draw生成TH1，合并后调用相同绘图块，不以另一套分析选择代替。

本次第17页对照使用110051全部20文件、20,000事件；beforeGeo为Calypso本地master，afterGeo为旋转分支。材料图bug已解决，AFS rotated JSON与strict rotated JSON另作单因素对照。当前两个分支的原生seed covariance不同，整体版本比较与协方差匹配后的纯旋转比较分别命名。第19页的100044 μ+500k角度研究留到后续阶段；保留图中的全谱/pz窗口legacy及正确物理p窗口定义。

truth来源检查在同一份native writer ROOT上运行原`t_eQOP`与相同majority barcode的production q/p两次；另画共同有效states对照。两种参考都不是sensor局部G4 p；NaN、production-only/link-only states、逐bin差异和有限分母单列。实现为`configs/reference_tracking.py`、`scripts/prepare_reference_macro.py`、`scripts/run_reference_task.py`与`analysis/reference_figures.py`。全20文件的beforeGeo、afterGeo、strict主图及truth对照已交付，[复核报告](../workbook/2026-11-12_stage04_native_geometry_truth.md)链接14页完整图集、bin counts和覆盖验证；原宏的TTreeFormula除零/NaN行为与有限residual统计分别记录。

用户指定开启候选独立 smoothing 后的全量主图见[全110051 smoothing/covariance 对照](../workbook/2026-11-12_stage04_smoothing_covariance.md)：beforeGeo、afterGeo，以及 afterGeo 恢复原 inflation 三组各20文件、20,000事件；保留其余原生设置和原宏 cuts。16页图集同时给出 P/F/S 九层、几何/variance 叠图和 smoothed 左尾 log 图。相同数值 inflation 与经过 Jacobian 匹配的相同物理 covariance 分别表述。

小样本发现图与全输入图明确区分；全输入版本逐文件检查CollectionTree entries及无重叠、无遗漏的entry ranges，图上列出源事件数。predicted/filtered/smoothed诊断另保留同state的tail迁移计数、末端filtered=smoothed检查及covariance变化，区分前向信息不足、更新新生tail与后向修复。批量state在worker scratch归约，只保存histogram和紧凑诊断统计。

## 每张比较图的约定

1. 标明DSID/粒子、源事件范围、所比较的配置和stage，声明坐标系、物理参考面与truth来源。`before/after`必须说明改的是几何、材料还是其他具体配置。
2. 参考图对应的横轴量、分组和面板布局必须保留。历史坐标/单位/选择存在问题时，给明确标注的legacy复现及正确物理量的并排对照；不能把两者混成同一个分布。没有历史输入时记录未恢复的条件。
3. before/after使用共同bin edges、横轴范围和归一化；对应面板使用可比较的纵轴范围。主图保留参考材料中的线性分布，必要时增加log或峰区放大图。
4. 图例或图注给有限entries、分母和tail/peak/misID计数；归一化包含的underflow/overflow、invalid以及丢失的tracks另存统计表，不通过裁剪横轴重归一化掩盖尾部。多层states是分布entries，率的统计单位仍是源事件/粒子。
5. 总体样本使用各配置自己的冻结选择，配对图使用共同truth并说明相同seed/hits等附加条件；good-seed筛选图列清p、charge、角度阈值及发生位置。失败/pruning留在相应效率的母分母，不能用最终成功轨迹重新定义good-seed效率。
6. 样本量、预留验证、效率和选择迁移与图一起报告。标题、坐标轴、图例和标注全部使用English。

## 文件与阶段完成标准

- results目录保存逐图PNG与可直接插入幻灯片的PDF，以及一份多页PDF汇编；主要曲线保留vector输出。
- 保存bin edges、原始bin counts、sumw/sumw2、underflow/overflow/invalid、selection和输入来源，提供可复现绘图命令；图索引记录每张图对应的参考报告/页码和产物路径。
- 中文workbook嵌入该阶段主要before/after图并链接完整PDF；README或roadmap提供可发现的入口。只列指标表不构成完整的可视化交付。
- 完成前检查图中文字、面板顺序、legend、轴范围和导出裁切，核对图中entries/tail与冻结summary一致；批量统计按仓库规范通过HTCondor执行。

用户追问 smoothing 是否实际生效的检查见[GitHub cuts 下逐参数复核](../workbook/2026-11-12_stage04_github_smoothing_recheck.md)：保留 GitHub cuts，补充全量 OFF/ON 九层叠图和旧 audit 实现到 writer 的实际参数检查。
