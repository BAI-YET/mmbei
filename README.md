# MathorCup 2026赛道B：四问代码、方法与PNG成果（v2.2）

本仓库用于妈妈杯代码图片上传保存留档，以确保能够正确回退版本。本文整合数据准备、四问作答、模型原理、运行方法、图题图注与全部图像索引。当前成果由实际输入、分析计算和数值实验生成；参数的观测、文献、标定和实验设计角色分别记录。

## 1. 当前计算回答了什么

| 问题 | 计算结果 |
|---|---|
| Q1关联特征 | 62项候选、56项有效；按机制合并相关变量后，五项代表通道为水汽VPD、近地层热力、地表能量余量、短波辐射、土壤水分。原始逐变量排序同时保留。 |
| Q1极端过程 | 事件前七日VPD有29小时达到匹配历史q95，其中3小时超过匹配历史记录；2m气温有15小时达到q95，最长连续6小时。 |
| Q2临界条件 | 非零24h降雨q90为10.20mm，源区7日正积温q90为39.05℃·d；新增30日热量、辐射能量和冻融通道。 |
| Q2本次能量通道 | 正能量余量7日当量为192.72mm，历史q90为189.46mm，达到历史第92.3百分位；现已进入主动态指数。 |
| Q3预测预警 | 6h条件超限频率0.853；按0.25阈值，预测预警时间2026-08-26 09:00:00，至题面参考时刻提前1.5h。 |
| Q3监测预警 | 有效监测段起点2026-08-25 18:00:00，其较题面参考时刻的回放提前量16.5h；与09:00预测发布分别记录。 |
| Q3区间与检验 | 条件首次超限95%区间2026-08-26 09:00:00至2026-08-26 15:00:00；2025年相似日基准检验151次发布、67次气象超限，校准后覆盖率100.00%，平均宽度5.84h；TP/FP/FN/TN=59/23/8/61。 |
| Q3提前量 | 20次正确预警的新发超限中，提前量中位数4.00h，四分位范围0.88至4.50h。 |
| Q3空间范围 | 源区至口岸模型廊道23.71km，变断面与二维湿区融合包络9.58km²。 |
| Q4分阶段复现 | 官方DEM汇流路径长23.71km、落差3427m；初始方量103百万m³，可动输入36百万m³，按实际求解通量连接上下游。 |
| Q4口岸结果 | 源区释放后403.5s到达口岸河段，峰值流量208522m³/s，对应断面平均速度43.26m/s。 |
| Q4二维网格试验 | 30m与60m累计湿润面积分别0.6246、0.6768km²；每档输入体积均与一维累计通量一致。 |

预测预警、持续监测状态和事后过程标定分别报告。95%区间对应条件首次超限，历史覆盖率评价实际气象超限；当前历史起报驱动的回放与相似日基准检验分别列示。校准前覆盖率91.04%；校准后平均区间比完整6h窗口缩短2.70%，同时报告宽度与覆盖率。模型约定集中于第6节。

## 2. 仓库上传范围

- 计算入口、全部模型代码、配置与测试：[run.py](run.py)、[src/mathorcup](src/mathorcup)、[configs](configs)、[scripts](scripts)、[tests](tests)。
- 输入模板：[events](events)，两个仅含表头的CSV用于登记核验事件与路径。
- 成果目录仅提交80张PNG，分布于`outputs/q1`、`q2`、`q3`、`q4`、`joint`。
- `.gitignore`使生成的CSV、JSON、NPZ、GIS、SVG、GIF、PDF、图廊、缓存和重绘副本留在本地。原始附件、下载产品和字体按下一节恢复。

`configs/case_evidence.yaml`是带出处和角色的事件约束；`configs/calibrated_q4.yaml`保存本次数值实验得到的有效参数与输入指纹，属于复现模型配置。

## 3. 环境、数据与复现

推荐Python 3.12。在仓库根目录执行：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Windows使用`.venv\Scripts\activate`激活。中文绘图采用Noto Sans CJK SC，可安装系统字体，或从[Noto CJK官方目录](https://github.com/notofonts/noto-cjk/tree/main/Sans/OTF/SimplifiedChinese)取得`NotoSansCJKsc-Regular.otf`、`NotoSansCJKsc-Bold.otf`放入`assets/fonts/`；PDF生成使用这两个字体。

### 3.1 官方八类附件

将官方ZIP保留`（1）`至`（8）`的全角编号文件名放入`data/raw/`，或用`--official-dir`指定。程序解压到`data/raw/unpacked/编号/`。第③组RAR单独提取目标区域四期瓦片。

| 附件 | 内容与用途 |
|---|---|
| ① | MET/RADM/SOIL/FLUX逐小时数据：近地层结构、气象过程、能量、土壤与历史预警检验。 |
| ② | 2005至2025年县级日温：长期暖季背景。 |
| ③ | 冰川2000至2020四期高程变化：长期物源演化。 |
| ④ | 冰川高程、冰厚、历史速度和掩膜：物源分区与释放网格。 |
| ⑤ | 2020至2023年土地覆盖：有效覆盖区变化。 |
| ⑥ | 84条历史沟谷：成因类别分析与空间留出检验。 |
| ⑦ | DEM：汇流路径、实际横断面、1D与2D地形。 |
| ⑧ | 行政区人口资料：区域背景与暴露分析辅助。 |

第③组用7-Zip或Bandizip提取以下四幅至`data/raw/unpacked/3_selected/`：

```text
N28E085_2000-01-01_2005-01-01_dhdt_UTM45_UTM45_GCS.tif
N28E085_2005-01-01_2010-01-01_dhdt_UTM45_UTM45_GCS.tif
N28E085_2010-01-01_2015-01-01_dhdt_UTM45_UTM45_GCS.tif
N28E085_2015-01-01_2020-01-01_dhdt_UTM45_UTM45_GCS.tif
```

### 3.2 补充遥感与历史起报预报

Q2纳入[ICIMOD的Sentinel-1分析摘要](https://www.icimod.org/press-releases/second-major-himalayan-flood-in-seven-weeks-raises-new-questions-about-mountain-risk/)，其中2026年7月7日至8月24日的视线形变上限约0.40m，作为事后物源演化资料；公开多传感器源区解译辅助Q4定位。资料形态、发表时间与预报输入的可用时间分别登记。


Q2从[Open-Meteo Single Runs API](https://open-meteo.com/en/docs/single-runs-api)获取实际起报。当前PNG对应2026-08-25 18:00 UTC起报，按起报后6小时可用的假设检查其与09:00发布的关系。原始JSON、小时数据和provenance保存在`data/external/q2_forecast/`。Q3读取Q2已有预报；Q1、Q3、Q4不增加天气下载步骤。

离线复现时，从此前完整复现包恢复`ecmwf_source_pre_event_single_run.json`、`source_forecast_hourly.csv`、`provenance.json`到该目录。联网获取与运行：

```bash
python run.py --stage all --official-dir "/path/to/official_archives" --fetch-q2-external
```

恢复数据与预报缓存后的完整计算：

```bash
python run.py --stage all --no-fetch-q2-external
python -m pytest -q tests
python scripts/audit_provenance.py
python scripts/build_answer_pdf.py
```

支持`prepare`、`q1`、`q2`、`q3`、`q4`分阶段运行，按前后依赖顺序执行。`--force-rebuild`更新指定处理缓存，`--config`使用另一配置。输入或核心求解器变化时按指纹自动重标定；`--refit-q4`主动重复标定试验。重跑某阶段会刷新其成果，自己修改的本地图像应另存。

`--strict-route`在清理任何旧结果之前检查`route_verified`与至少两个有效核验路径点，并约束全部Q4入口；严格模式实际沿核验路径计算。它与`--allow-provisional-q4`互斥。

## 4. 四问模型与代码

### 4.1 Q1：边界层、天气过程、极端事件和五项因素

代码：[prepare_time.py](src/mathorcup/prepare_time.py)、[features.py](src/mathorcup/features.py)、[q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py)、[question_completion.py](src/mathorcup/question_completion.py)。

登记63项原始非QC字段（MET22、RADM4、SOIL12、FLUX25），分别应用可用的QC、缺失及物理范围检查。62项物理候选中本次56项满足排序条件。24h为主窗口，另比较6、72、168h及参考年份重采样。

匹配2021至2025年同期有效窗口，计算`Z=(事件值-历史中位数)/(1.4826×MAD)`，`S=0.45×min(|Z|/3,1)+0.55×(1-p尾部)`。原始前五项是VPD、10m气温、20m气温、能量余量、2m气温。作答中的五项互补机制按同物理通道及同期Spearman绝对相关0.85合并，得到VPD、10m气温、能量余量、短波辐射、土壤水分梯度。

边界层以多高度风温、矢量风切变、温度梯度及近似Richardson数描述；天气过程用逐小时观测串联；极端性明确统计同季同小时q95、q99、破纪录小时和连续小时。能量交换与土壤廓线分别出图。原始排序、机制代表项、年份稳定性和极端小时均由代码导出。

### 4.2 Q2：物源、水源和能量通道

代码：[q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[question_completion.py](src/mathorcup/question_completion.py)。

静态物源指数来自冰厚、坡度、历史速度与应变代理，驱动剪应力为`τ=ρgH sinθ`。速度梯度保持UTM北向符号，缺失保持未知。高分连通冰体提供候选源区，冰库存与实际可动体积分别记录。

主动态指数取降雨、源区7日正积温、站点正能量余量三个门控的最大值。30日热量、净短波能量和冻融转换继续单列。正通量按`W/m²×3600/334000`换算能量水当量，1kg/m²数值等于1mm；阈值来自所指定历史暖季。未来能量通道仅使用仍位于七日窗口内的灾前实测正能量贡献，作为非负累积能量的下界，不补造未来能量观测。源区高程q10/q50/q90与三档温度递减率、四组静态权重均实际重算。

历史84条沟谷采用成因类别作为目标，比较全部候选模型的空间留出ROC-AUC；时间可靠性逐年向前拟合并比较更新间隔，空间可靠性比较30、60、120m的秩相关与高分区IoU。实际预报的温度、降水和短波进入未来条件；该次零度层缺失随原始记录保存。

### 4.3 Q3：95%时间区间、预警和空间距离

代码：[q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)。

固定气象变量使用气温、短波和VPD，联合Q2水热通道；2026事件选出的因素用于补充解释。采用历史完整相似日轨迹并按最新观测作衰减调整，未来Q2通道由实际历史起报预报驱动。首次超限的频率和时间分布直接从数据轨迹计算；联合到达时间仅重采样实际计算的时间格点与传播试验结果。

预测预警以`P(未来6h条件超限)≥0.25`判定；监测状态机处理2h进入、3h低值解除、6h续期和缺测暂停。当前预测发布、已有监测段起点和各自有效期分列；决策函数没有未来事件时刻输入，题面时刻只用于事后评价。历史检验将2021至2023、2024、2025分别用于拟合、校准、检验，保存TP/FP/FN/TN、分母、提前量与覆盖过程。敏感性曲线由不同阈值和7/14/21/28天历史轨迹窗口实际重算。

空间回放以Q4河道和变断面最大水位形成廊道，融合释放网格及原始DEM上的二维湿区。研究段为源区至口岸，距离和面积来自米制几何计算。单阶段参考模型的空间范围仍保留用于实验比较，主作答采用`Q3_12`与分阶段结果。

### 4.4 Q4：分阶段级联与贴合地形的传播仿真

代码：[q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)。

1. **真实地形**：扩大计算范围以覆盖主河道的北向弯曲；Priority-Flood只辅助D8寻路，动力学保持官方原始高程。汇入口邻域结合文献河道距离与DEM集水面积变化定位。
2. **实际断面**：提取沿程DEM横剖面，以20m、60m相对高程的连通宽度拟合梯形。采用`A=b h+m h²`和静水压力矩`g(b h²/2+m h³/3)`；有效断面尺度由事件约束标定。
3. **分布式释放**：先在源区上部真实冰川像元处建立地形路径，再按UTM东、北坐标和高程作三维投影，保留逐像元位置与高程差；初速度0。该方法避免把源区坐标整体平移到更低的沟谷。事件研究提供103百万m³初始方量及0.76km²源区面积，36百万m³可动方量与67百万m³沉积分配分别约束。
4. **守恒与动量**：上游保留适用于本次强地形变化的原变断面Rusanov-Voellmy计算；下游使用静水重构MUSCL-HLL与SSP-RK2步进降低数值耗散。先计上游实际滞留，再按其累计通量分配汇入口沉积与下游输入，边界只施加相应体积与运动速度。
5. **事后标定**：上游时间采用双侧误差约束，下游联合比较传播时间、峰值流量和断面速度。以源区释放为统一时间零点，导出前锋、洪峰、5%/50%累计量时刻。文献约130s与350s分别对应两段，其约480s总经历时间另行对照。参数边界状态和全部标定偏差共同输出。
6. **真实通量耦合**：局部二维接收一维每个接受时间步的累计界面体积与速度；30m和60m网格独立计算。二维强度报告流深、速度、到达时间与单位质量动能`u²/2`。
7. **数值与地形检查**：上游15m、下游7.5m主网格，另以3.75/7.5/15m比较下游分辨率；二维30/60m独立计算。标定在15m网格进行，固定所得参数后在更细输出网格评价，不将加密结果写成重新标定。网格试验保持同一入流和物理端点，分别检查守恒、压力梯度、干湿前缘及分辨率响应。

本次下游标定参数为μ=0.003006、ξ=2711.18m/s²、有效断面尺度3.000（触及搜索边界）。全程体积相对残差1.88e-15。`Q4_14`至`Q4_25`为分阶段计算与复核图；早期编号为同一官方数据上的单阶段参考实验。

## 5. 图表数据、重绘与测试

每次`common.savefig`自动导出实际曲线、柱形、散点、区间、多边形、栅格、掩膜和坐标，保存在`outputs/redraw_data/`。首先打开`00_重绘数据索引.csv`，按PNG文件名查找数据；地图另有GeoTIFF，原始结果表在`source_tables/`。这些本地文件支持Origin、Excel或GIS自行重绘。

```bash
python scripts/replot_reconstruction.py
python scripts/check_redraw_retention.py
python scripts/export_redraw_data.py
```

第一条从缓存数组重绘分阶段结果；第二条逐图校对数值保留；第三条在隔离目录回放完整绘图并检查原结果保持一致。缓存重绘和隔离回放均需先准备完整本地输入与数值结果。`replot_q4_figures.py`、`replot_q4_1d_cached.py`对应参考单阶段图。`build_answer_pdf.py`直接读取当前表与图，README由本脚本按当前结果更新。

本次84项测试通过。测试覆盖QC和单位、未来信息隔离、预警状态、误报分母、静水平衡、真实压力梯度、干湿前缘、变断面守恒、能量范围、路径保护与累计体积耦合。解析小网格只用于验证求解器，正式个案输入与输出由来源清单单独审计。

## 6. 模型约定与尺度解释

| 环节 | 统一解释 |
|---|---|
| 关联与阈值 | Q1评分描述本次事件相对历史的异常关联；Q2 q90是条件评估阈值。历史沟谷AUC的目标是成因类别。 |
| 跨境推算链 | 站点4140m、源区代表高程约5297m，以6.5℃/km主递减率并比较5与8℃/km；实际预报格点偏移约8.38km。站点能量按站点尺度解释。 |
| 时间与覆盖率 | 95%区间是条件首次超限的时间预测区间，经验覆盖率由历史气象超限事件计算；题面10:30作为影响参考时刻。 |
| 标定资料与空间 | Q4结合事后事件研究进行过程复现；Q3天气输入与Q4事后标定资料分别记录，空间范围属于模型回放包络。 |
| 断面和网格 | 有效断面尺度3.000，触及搜索边界；一维7.5m至3.75m的峰值流量变化+2.4%，到达时间变化-1.2%，最大流深变化13.2%；二维面积变化8.4%。质量闭合与各量的分辨率响应分别报告。 |
| 传播强度 | 等效介质方程描述体积与动量演变；研究段终点设在口岸。动能按单位质量报告，结构受力可结合介质组成和构件参数进一步计算。 |

## 7. 资料出处与配置来源

- 官方数据：赛题提供的国家青藏高原科学数据中心附件①至⑧。
- [Li H等，2026，Chhochen Khola灾害链事件研究](https://jms.imde.ac.cn/en/article/doi/10.1007/s11629-026-1482-y)，DOI 10.1007/s11629-026-1482-y：分阶段方量、时间及断面重建约束。
- [Li R，arXiv:2609.04563](https://arxiv.org/abs/2609.04563)：事后源区位置资料。
- [RAMMS官方摩阻说明](https://ramms.ch/ramms-debrisflow/friction-parameters/)：Voellmy形式与事件标定方法。
- [Open-Meteo Single Runs](https://open-meteo.com/en/docs/single-runs-api)：实际历史起报响应。
- [尼泊尔DHM站点4913](https://www.dhm.gov.np/hydrology/hms-Single/4913)：口岸下游水文站地理位置辅助核对。
- [口岸附近原图地理标记](https://commons.wikimedia.org/wiki/File:Nepal_China_Border.jpg)：区域参考点；数值终点投影至DEM河槽。
- Audusse等，2004，SIAM Journal on Scientific Computing，DOI 10.1137/S1064827503431090：静水重构方法。
- Gottlieb与Shu，1998，Mathematics of Computation，DOI 10.1090/S0025-5718-98-00913-2：两阶段TVD/SSP时间步进。

输入的逐项数值、来源URL和角色见`configs/case_evidence.yaml`。完整输入、源码和配置的SHA256记录由`audit_provenance.py`生成，位于`outputs/quality/reproducibility_manifest.json`。

## 8. 审查总结逐项复核

原审查为静态判断；本次依据赛题、最新代码与实际重算进行复核。下表保留合理的模型解释要求，并明确哪些属于实现错误。

| 序号 | 审查项 | 判断 | 实际处理 |
|---|---|---|---|
| 1 | Q3预测对象 | 成立：属于预测对象与验证含义 | 主模型输出条件首次超限时间，历史检验目标也为气象条件；本次题面时刻单例命中另列。 统一图题、字段和中位时刻的离散支持；当前天气驱动回放与历史相似日检验分别登记，保留无超限的零概率输出。 |
| 2 | Q2临界条件与能量通道 | 部分成立：历史分位定义合理，能量衔接确有遗漏 | q90是明确的历史参考阈值；静态权重是声明的指数设计参数。 把实测正能量余量纳入主动态指数，补算静态权重与高程递减率敏感性；未来能量只取实测历史贡献的下界。 |
| 3 | Q4时间拟合 | 成立：旧版数值偏差与单侧时间约束均已复算 | 旧版汇入口至口岸261.71s；本次304.94s，文献约350s。 将上游时间改为双侧标定约束，统一源区释放时钟，分别导出前锋、洪峰、5%与50%累计量时间。文献参照保留为近似重建值。 |
| 4 | Q4网格敏感性 | 成立：体积守恒与分辨率响应是两项不同检查 | 本次两档较细网格峰值差2.41%，到达时间差1.19%。 下游改用静水重构MUSCL-HLL与SSP-RK2；网格试验保持同一入流、参数及两端位置，分别报告各量变化。 |
| 5 | Q4局部二维范围 | 事实正确，局部单向嵌套本身是合理实现 | 一维覆盖源区至口岸，二维覆盖口岸邻域；输入为上游实际到达局部断面的累计通量。 新增全程一维过程的地图投影，并与局部二维湿区分列，记录各计算域及对应体积收支。 |
| 6 | Q4参数测量约束 | 事实正确，标定参数不等于程序错误 | 本次有效断面尺度3.0000，边界状态见参数诊断。 保留观测、文献重建、标定参数、数值设置的角色；自动记录搜索上下界与是否触界，不用改写观测或强行消除残差。 |
| 7 | Q2源区气象代表性 | 成立：属于空间代理与资料可用性约定 | 站点、冰川代表高程与预报格点的空间尺度不同，运行时延按已有资料声明检查。 增加冰川高程q10/q50/q90与三档递减率的九组合重算；保留实际格点偏移和预报来源。 |
| 8 | Q1异常与因果 | 判断合理；当前异常排序方法可保留 | 题面要求关联性排序，单次事件与历史匹配可定义异常关联；该评分本身没有估计因果效应。 保持全部候选、原始前五、机制代表项和年份敏感性共同输出；用异常关联说明评分含义。 |
| 9 | strict-route保护 | 成立：是真实的程序控制缺口 | 旧入口在常规路径分支后仍能进入级联重建。 在清理旧图前校验，所有级联入口执行同一检查；严格模式实际使用核验路径，冲突命令参数直接拒绝。 |
| 10 | 代码副本不一致 | 在本次最新交付包中未复现 | 压缩包只有一个代码根目录；61/61个Python文件与基准GitHub版本逐一相同。 继续使用唯一源码根目录，结果清单记录输入、源码、参数、图表及版本的校验摘要。 |
| 11 | 另查出：源区释放投影 | 成立：原PCA平移改变了初始释放势能 | 旧版释放平均床高4690.3m，真实像元均值5132.0m；修正后释放平均床高5116.2m，逐像元投影平均绝对高程差29.3m。 改为源区上部真实像元起始的地形路径和UTM三维最近点投影；逐像元导出高程差，并用修正后的释放位置重新标定与求解。 |

## 9. 全部PNG的用途与源码索引

下表链接为当前PNG；“数值证据”列给出本地计算结果文件名。这些表和缓存按本次上传范围留在本地。

### 问题一（17张）

| PNG | 用途 | 生成代码 | 数值证据 |
|---|---|---|---|
| [Q1_00_data_coverage.png](outputs/q1/Q1_00_data_coverage.png) | 逐小时观测数据有效率 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_01_event_hourly_timeseries.png](outputs/q1/Q1_01_event_hourly_timeseries.png) | 吉隆泥石流灾前气象要素变化 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_02_boundary_layer_profiles.png](outputs/q1/Q1_02_boundary_layer_profiles.png) | 高寒山区近地层风温垂直结构 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_03_top15_anomaly_ranking.png](outputs/q1/Q1_03_top15_anomaly_ranking.png) | 灾前多源气象因素异常稀有性排序 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_04_multiscale_anomaly_heatmap.png](outputs/q1/Q1_04_multiscale_anomaly_heatmap.png) | 多时间尺度气象异常分布 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_05_reference_year_sensitivity.png](outputs/q1/Q1_05_reference_year_sensitivity.png) | 关键特征排序的年份敏感性 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_06_pseudo_event_controls.png](outputs/q1/Q1_06_pseudo_event_controls.png) | 灾前气象异常与历史同期对照 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_07_county_2005_2025_summer_climate.png](outputs/q1/Q1_07_county_2005_2025_summer_climate.png) | 吉隆及周边地区夏季气温变化 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_08_boundary_layer_seasonal_structure.png](outputs/q1/Q1_08_boundary_layer_seasonal_structure.png) | 近地层稳定结构的季节变化 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_09_radiation_and_turbulent_energy.png](outputs/q1/Q1_09_radiation_and_turbulent_energy.png) | 灾前地表辐射与热量通量变化 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_10_multidepth_soil_profiles.png](outputs/q1/Q1_10_multidepth_soil_profiles.png) | 灾前土壤多深度温湿廓线 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_11_all_candidate_ranking.png](outputs/q1/Q1_11_all_candidate_ranking.png) | 完整候选排序；对应Q1全部因素比较与原始前五项 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_11_part1_complete_ranking.png](outputs/q1/Q1_11_part1_complete_ranking.png) | 完整候选排序；对应Q1全部因素比较与原始前五项 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_11_part2_complete_ranking.png](outputs/q1/Q1_11_part2_complete_ranking.png) | 完整候选排序；对应Q1全部因素比较与原始前五项 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_12_all_factor_multiscale.png](outputs/q1/Q1_12_all_factor_multiscale.png) | 完整多时间尺度异常；对应Q1窗口依赖 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_13_all_factor_redundancy.png](outputs/q1/Q1_13_all_factor_redundancy.png) | 全部候选冗余结构；区分原始排序与机制代表项 | [q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_14_mechanisms_and_extreme_weather.png](outputs/q1/Q1_14_mechanisms_and_extreme_weather.png) | 五项互补机制与极端天气识别 | [question_completion.py](src/mathorcup/question_completion.py) | `Q1_direct_answer.json`；`Q1_explicit_extreme_weather.csv` |

### 问题二（19张）

| PNG | 用途 | 生成代码 | 数值证据 |
|---|---|---|---|
| [Q2_01_history_gullies_and_oof_roc.png](outputs/q2/Q2_01_history_gullies_and_oof_roc.png) | 历史沟谷成因分布与判别模型ROC曲线 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_02_spatial_holdout_auc.png](outputs/q2/Q2_02_spatial_holdout_auc.png) | 灾种判别模型空间交叉验证 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_03_elevation_drop_distributions.png](outputs/q2/Q2_03_elevation_drop_distributions.png) | 不同成因沟谷高差分布 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_04_climate_condition_timeseries.png](outputs/q2/Q2_04_climate_condition_timeseries.png) | 降雨与融水触发条件时间变化 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_05_glacier_physical_indicators.png](outputs/q2/Q2_05_glacier_physical_indicators.png) | 冰川厚度、坡度与驱动剪应力分布 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_06_dem_landcover_history_overlay.png](outputs/q2/Q2_06_dem_landcover_history_overlay.png) | 研究区地形与土地覆盖分布 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_07_glacier_historic_thinning_optional.png](outputs/q2/Q2_07_glacier_historic_thinning_optional.png) | 冰川历史高程变化空间分布 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_08_trigger_conditions_time_series.png](outputs/q2/Q2_08_trigger_conditions_time_series.png) | 降雨与融水触发条件变化 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_09_source_spatial_condition_index.png](outputs/q2/Q2_09_source_spatial_condition_index.png) | 冰川源区地形物源条件分布 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_10_spatial_resolution_reliability.png](outputs/q2/Q2_10_spatial_resolution_reliability.png) | 触发评估空间尺度可靠性 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_11_temperature_lapse_sensitivity.png](outputs/q2/Q2_11_temperature_lapse_sensitivity.png) | 源区热力条件递减率敏感性 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_12_forecast_weather_conditions.png](outputs/q2/Q2_12_forecast_weather_conditions.png) | 源区历史起报天气预报条件 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_13_glacier_velocity_and_release_candidates.png](outputs/q2/Q2_13_glacier_velocity_and_release_candidates.png) | 冰川历史运动与候选失稳区空间特征 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_14_landcover_2020_2023_comparison.png](outputs/q2/Q2_14_landcover_2020_2023_comparison.png) | 2020年与2023年土地覆盖对比 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_15_all_model_comparison.png](outputs/q2/Q2_15_all_model_comparison.png) | 全部候选模型的同一空间留出比较；非触发事件AUC | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_16_all_predictor_distributions.png](outputs/q2/Q2_16_all_predictor_distributions.png) | 全部五项判别变量分布与解释 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_17_forward_time_reliability.png](outputs/q2/Q2_17_forward_time_reliability.png) | 逐年向前阈值与不同更新间隔的监测误差 | [q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py) | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_18_energy_water_mechanism_channels.png](outputs/q2/Q2_18_energy_water_mechanism_channels.png) | 降雨、前期热量与能量通道评估 | [question_completion.py](src/mathorcup/question_completion.py) | `Q2_expanded_mechanism_answer.json`；`Q2_expanded_channel_thresholds.csv` |
| [Q2_19_source_proxy_and_weight_sensitivity.png](outputs/q2/Q2_19_source_proxy_and_weight_sensitivity.png) | 源区空间代理与静态权重敏感性 | [q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py) | `Q2_source_elevation_lapse_sensitivity.csv`；`Q2_static_weight_sensitivity.csv` |

### 问题三（13张）

| PNG | 用途 | 生成代码 | 数值证据 |
|---|---|---|---|
| [Q3_01_precursor_and_scenario_hazard.png](outputs/q3/Q3_01_precursor_and_scenario_hazard.png) | 灾前气象条件及触发风险演变 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_02_95percent_conditional_time_interval.png](outputs/q3/Q3_02_95percent_conditional_time_interval.png) | 触发条件首次超限时间与95%预测区间 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_03_warning_lead_and_backtest.png](outputs/q3/Q3_03_warning_lead_and_backtest.png) | 泥石流预警时间与气象条件回测 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_04_parameter_sensitivity.png](outputs/q3/Q3_04_parameter_sensitivity.png) | 预警模型参数敏感性分析 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_05_spatial_warning_distance.png](outputs/q3/Q3_05_spatial_warning_distance.png) | 泥石流传播廊道与空间预警距离 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_06_validation_metrics.png](outputs/q3/Q3_06_validation_metrics.png) | 气象条件时间区间与预警指标验证 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_07_Q1_feature_inheritance.png](outputs/q3/Q3_07_Q1_feature_inheritance.png) | Q3_07_Q1_feature_inheritance | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_07_feature_group_ablation.png](outputs/q3/Q3_07_feature_group_ablation.png) | Q1气象特征与Q2物理条件的贡献比较 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_08_full_forward_validation.png](outputs/q3/Q3_08_full_forward_validation.png) | 完整2025检验期，三组输入、区间与覆盖过程 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_09_warning_episodes.png](outputs/q3/Q3_09_warning_episodes.png) | 预警进入、维持、续期、解除和有效期 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_10_archived_forecast_sensitivity.png](outputs/q3/Q3_10_archived_forecast_sensitivity.png) | 实际起报天气预报进入Q2条件后的未来驱动敏感性 | [q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py) | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_11_warning_time_and_lead_distribution.png](outputs/q3/Q3_11_warning_time_and_lead_distribution.png) | 预警时间与提前量分布 | [question_completion.py](src/mathorcup/question_completion.py) | `Q3_warning_lead_answer.json`；`Q3_proxy_backtest_2025.csv` |
| [Q3_12_reconstructed_spatial_warning.png](outputs/q3/Q3_12_reconstructed_spatial_warning.png) | 源区至口岸的模型空间警戒范围 | [q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py) | `Q3_reconstructed_spatial_answer.json`；`Q3_reconstructed_warning_range.gpkg` |

### 问题四（30张）

| PNG | 用途 | 生成代码 | 数值证据 |
|---|---|---|---|
| [Q4_14_drainage_route_and_cross_sections.png](outputs/q4/Q4_14_drainage_route_and_cross_sections.png) | 源区至口岸的地形与断面 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_15_calibrated_cascade_process.png](outputs/q4/Q4_15_calibrated_cascade_process.png) | 分阶段冰岩崩与泥石流传播过程 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_16_event_constraint_comparison.png](outputs/q4/Q4_16_event_constraint_comparison.png) | 事后标定的定量对照 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_17_calibration_parameter_experiments.png](outputs/q4/Q4_17_calibration_parameter_experiments.png) | 同一输入下的标定参数试验 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_18_reconstructed_2d_depth_speed_arrival.png](outputs/q4/Q4_18_reconstructed_2d_depth_speed_arrival.png) | 口岸二维传播强度与到达时间 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_19_reconstructed_2d_six_frames.png](outputs/q4/Q4_19_reconstructed_2d_six_frames.png) | 真实通量驱动的六个传播阶段 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_20_cascade_volume_and_grid_experiments.png](outputs/q4/Q4_20_cascade_volume_and_grid_experiments.png) | 级联体积守恒与网格试验 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_21_cross_sections_and_1d_grid_experiment.png](outputs/q4/Q4_21_cross_sections_and_1d_grid_experiment.png) | 原始横断面与数值分辨率试验 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_22_time_alignment_and_residuals.png](outputs/q4/Q4_22_time_alignment_and_residuals.png) | 事件时间起点与全部标定残差 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_23_whole_route_six_stages.png](outputs/q4/Q4_23_whole_route_six_stages.png) | 源区至口岸的全程传播过程投影 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_24_resolution_response_comparison.png](outputs/q4/Q4_24_resolution_response_comparison.png) | 修订前后网格敏感性的完整比较 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_25_release_geometry_and_elevation.png](outputs/q4/Q4_25_release_geometry_and_elevation.png) | 源区三维投影与初始释放势能核查 | [q4_terrain_reconstruction.py](src/mathorcup/q4_terrain_reconstruction.py)、[q4_cascade_reconstruction.py](src/mathorcup/q4_cascade_reconstruction.py)、[q4_reconstructed_spatial.py](src/mathorcup/q4_reconstructed_spatial.py)、[q4_diagnostics.py](src/mathorcup/q4_diagnostics.py) | `Q4_cascade_reconstruction.json`；`Q4_reconstruction_geometry.json`；`Q4_reconstruction_2d_summary.json` |
| [Q4_00B_local_true_DEM_detail.png](outputs/q4/Q4_00B_local_true_DEM_detail.png) | 参考单阶段实验；下游二维计算区域地形图 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_00_actual_DEM_geographic_context.png](outputs/q4/Q4_00_actual_DEM_geographic_context.png) | 参考单阶段实验；吉隆冰崩—泥石流研究区地形图 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_01_route_geolocation_quality_gate.png](outputs/q4/Q4_01_route_geolocation_quality_gate.png) | 参考单阶段实验；冰崩—泥石流传播路径及地形纵剖面 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_02_1d_runout_timespace.png](outputs/q4/Q4_02_1d_runout_timespace.png) | 参考单阶段实验；一维泥石流传播及地形纵剖面 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_03_1d_hydrograph_and_mass_conservation.png](outputs/q4/Q4_03_1d_hydrograph_and_mass_conservation.png) | 参考单阶段实验；耦合断面流量与体积守恒检验 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_04_1d_12frames.png](outputs/q4/Q4_04_1d_12frames.png) | 参考单阶段实验；一维泥石流流深演进过程 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_05_local_2d_depth_and_arrival.png](outputs/q4/Q4_05_local_2d_depth_and_arrival.png) | 参考单阶段实验；泥石流二维最大流深与到达时间 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_06A_onset_and_wet_area_diagnostic.png](outputs/q4/Q4_06A_onset_and_wet_area_diagnostic.png) | 参考单阶段实验；耦合断面入流及淹没面积变化 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_06B_inlet_position_ablation.png](outputs/q4/Q4_06B_inlet_position_ablation.png) | 参考单阶段实验；二维入流位置对传播范围的影响 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_06C_2x2_inlet_boundary_factorial.png](outputs/q4/Q4_06C_2x2_inlet_boundary_factorial.png) | 参考单阶段实验；二维入流位置与动量组合试验 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_06_local_2d_10frames.png](outputs/q4/Q4_06_local_2d_10frames.png) | 参考单阶段实验；泥石流二维流深演进过程 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_07_parameter_sensitivity.png](outputs/q4/Q4_07_parameter_sensitivity.png) | 参考单阶段实验；泥石流动力学参数敏感性分析 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_08_volume_friction_phase_diagram.png](outputs/q4/Q4_08_volume_friction_phase_diagram.png) | 参考单阶段实验；物源体积与摩擦系数对传播的影响 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_09_grid_resolution_reliability.png](outputs/q4/Q4_09_grid_resolution_reliability.png) | 参考单阶段实验；二维模拟网格分辨率敏感性 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_10_impact_pressure_distribution.png](outputs/q4/Q4_10_impact_pressure_distribution.png) | 参考单阶段实验；二维流动峰值动压分布 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_11_six_key_frames.png](outputs/q4/Q4_11_six_key_frames.png) | 参考单阶段实验；六个关键传播阶段；共同色标，来自本次二维数组 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_12_secondary_inflow_transfer.png](outputs/q4/Q4_12_secondary_inflow_transfer.png) | 参考单阶段实验；次级预设方量情景的真实1D通量、速度及累计量 | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_13_1d_grid_resolution.png](outputs/q4/Q4_13_1d_grid_resolution.png) | 参考单阶段实验；Q4_13_1d_grid_resolution | [q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py) | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |

### 触发至到达联合计算（1张）

| PNG | 用途 | 生成代码 | 数值证据 |
|---|---|---|---|
| [Q34_trigger_and_arrival_time.png](outputs/joint/Q34_trigger_and_arrival_time.png) | 触发条件与下游到达时间分布 | [q34_joint_arrival.py](src/mathorcup/q34_joint_arrival.py) | `Q34_status.json` |
