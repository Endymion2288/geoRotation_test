# S04补充：AFS旋转材料图与strict current22的逐项比较

日期：2026-11-12。结果目录为`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum/material_map_comparison`。

## 阶段目标

核对用户给出的`/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json`，与S04实际使用的`/eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/material-current-rotated.json`有何差异，区分文件结构、材料数值和实际局部分箱坐标。

## 结论

两份JSON存在结构差异，但共同12个surface的86400个材料bin完全相同，旋转矩阵和分箱定义也完全相同。strict额外保留10个surface，共118个bin，全部为`material=null`且`thickness=0`。共同surface的transform translation相差2 mm，方向沿材料面的法向，不改变当前局部`binX/binY`查询。

因此，不能把12与22的surface数量差解释为AFS文件缺少10层非零材料。AFS文件已经包含正确旋转；S04已报告的尾部改善对应未旋转current22→strict current22，并非AFS旋转图→strict的改善。strict保留原current22的所有IDs和transform来源，是研究纯坐标旋转时的受控配置。

## 工作内容与详细数据

新增[compare_material_maps.py](../scripts/compare_material_maps.py)，解析原JSON，按geometry identity关联条目；精确比较所有材料cell、metadata、binning与transform。脚本只检查Calypso输入数据，不产生重建结果。

| 项目 | AFS旋转图 | S04 strict旋转图 | 核对结果 |
| --- | --- | --- | --- |
| 文件字节数 | 45045313 | 11866365 | strict使用紧凑JSON格式，大小不代表材料多少 |
| Surface entries | 12 | 22 | 共同12项，仅strict多10项 |
| Volume entries | 0 | 0 | 相同 |
| JSON hierarchy header | format-version=0 | format-version=0 | 完全相同 |
| 共同surface geometry IDs | volumes 2/25/27/29、layers 2/4/6、approach=1 | 同左 | 完全相同 |
| 共同材料grid | 每surface 120行×60列，共86400 bins | 同左 | 全部`material`数组及`thickness`精确相同 |
| 共同metadata | `mapMaterial`、`mappingType`、`type` | 同左 | 完全相同 |
| 共同局部分箱 | binX 60 bins、binY 120 bins | 同左 | 轴、范围、open/type均相同 |
| 共同旋转矩阵 | F→A：`(z_F,y_F,-x_F)` | 同左 | 完全相同 |
| 共同平移 | 下表 | 下表 | 仅global A-frame x差±2 mm |

共同surface的平移只有A-frame x分量非零，单位mm；下表的差为AFS−strict。这是材料bin utility的transform，不是一次修改实际tracking geometry surface位置。

| ACTS volume | layer | AFS x_A / mm | strict x_A / mm | 差 / mm |
| --- | --- | --- | --- | --- |
| 2 | 2 | −1883.1075 | −1885.1075 | +2 |
| 2 | 4 | −1851.6075 | −1853.6075 | +2 |
| 2 | 6 | −1820.1075 | −1822.1075 | +2 |
| 25 | 2 | 7.3575 | 9.3575 | −2 |
| 25 | 4 | 38.8575 | 40.8575 | −2 |
| 25 | 6 | 70.3575 | 72.3575 | −2 |
| 27 | 2 | 1197.3575 | 1199.3575 | −2 |
| 27 | 4 | 1228.8575 | 1230.8575 | −2 |
| 27 | 6 | 1260.3575 | 1262.3575 | −2 |
| 29 | 2 | 2387.3575 | 2389.3575 | −2 |
| 29 | 4 | 2418.8575 | 2420.8575 | −2 |
| 29 | 6 | 2450.3575 | 2452.3575 | −2 |

strict额外保留的10项均为`sensitive=1`，具体如下：

| volume | layer | bins | 非空或非零厚度bins |
| --- | --- | --- | --- |
| 3、4、5、6、7 | 1 | 各4 | 0 |
| 9、10 | 2 | 各4 | 0 |
| 8、26、28 | 2 | 各30 | 0 |
| 合计10项 | — | 118 | 0 |

## 查表语义与Calypso结果范围

实际ACTS的`BinUtility.hpp`第153–161行先计算`bPosition = m_itransform * position`，再按各binning axis查表。该头文件位于：

```text
/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/include/Acts/Utilities/BinUtility.hpp
```

两图共同surface都使用`u=-z_A`、`v=y_A`、`w=x_A-t_x`的局部坐标，分箱轴只有u/v。±2 mm的`t_x`差仅改变w，因此对任意有限global position，共同材料的bin选择相同。比较脚本逐surface检查`Rᵀ(t_strict−t_AFS)`，其局部x/y分量均精确为0。

零材料条目与没有material object仍是不同的配置表示。完整tracking还包含surface导航、seed顺序与branch/pruning，不能从上述材料数据一致直接推出逐事件输出完全一致。

已有S03实际Calypso小样本使用这两个材料路径，对同一110051文件entries0–9运行，保留了日志与审计输出；本次重新读取这些结果，没有新运行reconstruction。两次run均为同一Calypso安装、同一`software_v2`库及相同输入/flags。重新核对结果：

| 观察量 | historical/AFS路径 | consistent/strict路径 |
| --- | --- | --- |
| 源事件 | 10 | 10 |
| Seed candidates | 34 | 34 |
| Found candidates | 807 | 863 |
| Final tracks | 21 | 21 |
| Final hit集合共同项 | 21 | 21 |
| 相同seed且momentum逐值相同 | 18 | 18 |

逐事件seed mean的多重集合全部相同；source entry9的seed编号与候选结果不同。按相同final hit集合配对，另外3条final来自不同seed，p分别为385.112112→385.012935 GeV、9794.391064→36287.990978 GeV、7.686140→7.691162 GeV；所有21条配对charge相同。不能将这3条的差异直接归因于2 mm偏移或额外空surface，具体候选执行原因未定位。

历史任务记录了AFS路径，但没有在run目录单独保留当时的材料JSON快照；上述是历史Calypso作业对照，当前两份文件的逐项结论以本次直接读取为准。这里不提供两份旋转图的总体性能差值。

## 验证与产物

```bash
python3 scripts/compare_material_maps.py \
  --first /afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json \
  --second /eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/material-current-rotated.json \
  --output /eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum/material_map_comparison/comparison.json
```

[完整逐项比较](../../geoRotation_test_results/stage04_momentum/material_map_comparison/comparison.json)保存来源路径、size/mtime、每个共同及额外条目的材料/transform检查；[既有Calypso小样本核对](../../geoRotation_test_results/stage04_momentum/material_map_comparison/historical_calypso_check.json)保存run来源、counts和21条配对。[验证记录](../../geoRotation_test_results/stage04_momentum/material_map_comparison/validation.json)汇总语法、帮助、本地链接与diff检查。没有修改两份输入材料图，也没有新增tracking实现。

## 下一步计划

保留strict current22作为当前研究的纯旋转控制，并在材料配置记录中明确AFS已旋转图与未旋转current的区别。如果需要量化AFS图与strict的全链性能差异，在实际Calypso内固定物理seed及输入顺序，分别控制空条目和法向平移，再做同输入验证；验收同时检查candidate迁移、p/charge、效率和导航记录。此项本次未测量，不将已有current→strict改善移用于该比较。
