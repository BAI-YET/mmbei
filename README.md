# MathorCup 2026 赛道B：四问模型、代码与PNG成果（v2.0）

本仓库用于妈妈杯代码图片上传保存留档，以确保能够正确回退版本。这里保存四问的计算与绘图代码、配置、测试及本次计算生成的63张PNG图像。图像中的数值来自模型输出，调整模型或输入后应重新运行相应代码。

原README与10份历史Markdown说明已按主题整合到本文：运行与数据准备、四问方法、数值验证、外部资料、图题图注、重绘及版本差异均在此查阅。历史文档中与v2.0计算冲突的数字和结论由本文替代，不再并排保留为当前结论。

## 1. 本次上传范围

| 内容 | 位置与范围 |
|---|---|
| 计算入口 | [run.py](run.py)、[configs/base.yaml](configs/base.yaml) |
| 模型源码 | [src/mathorcup](src/mathorcup)，含预处理、四问计算、绘图和重绘数据导出 |
| 辅助脚本与测试 | [scripts](scripts)、[tests](tests) |
| 人工核验输入模板 | [events](events)，2个仅含表头的CSV模板，供填写真实事件和核验路径 |
| 已生成成果 | [outputs/q1](outputs/q1) 16张、[outputs/q2](outputs/q2) 17张、[outputs/q3](outputs/q3) 11张、[outputs/q4](outputs/q4) 18张、[outputs/joint](outputs/joint) 1张，合计63张PNG |
| 统一说明 | 本README，包含全部图像的用途、源码与本地数值证据索引 |

生成后的CSV、Excel、JSON、NPZ、GIS、SVG、GIF、PDF、图廊和重绘数据副本均不提交。原始附件、下载产品、处理缓存和字体也不随仓库上传，复现时需在本地准备。`.gitignore`只允许五个成果子目录中的PNG进入版本控制，计算所需的中间表仍会在本地正常生成。两个CSV输入模板不是计算结果表。

## 2. 环境、数据与运行

### 2.1 安装与字体

推荐Python 3.12，在仓库根目录执行：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest -q tests
```

Windows可用`.venv\Scripts\activate`激活环境。依赖包含NumPy、Pandas、SciPy、Matplotlib、scikit-learn、Rasterio、GeoPandas等，详见[requirements.txt](requirements.txt)。

中文图像使用Noto Sans CJK SC。可安装系统中文字体，或从[Noto CJK官方字体目录](https://github.com/notofonts/noto-cjk/tree/main/Sans/OTF/SimplifiedChinese)取得`NotoSansCJKsc-Regular.otf`、`NotoSansCJKsc-Bold.otf`，放入`assets/fonts/`。普通绘图会尝试系统字体回退；可选PDF脚本明确需要这两个字体文件。已上传PNG已包含绘制好的中文，不依赖查看者安装字体。

### 2.2 官方附件

保留官方附件文件名中全角括号编号`（1）`至`（8）`。将7个ZIP附件放入`data/raw/`，或用`--official-dir`指定附件目录；程序自动解压到`data/raw/unpacked/1`至`8`的相应目录。第③组RAR按下节单独处理。也可以恢复此前完整压缩包中的解压数据与缓存到相同目录。

| 附件 | 在模型中的用途 |
|---|---|
| ① MET/RADM/SOIL/FLUX逐小时观测 | Q1气象与边界层；Q2水热条件；Q3历史拟合、校准与检验 |
| ② 县级历史逐日气温 | 长期气温背景，明确各年的实际覆盖 |
| ③ 冰川历史高程变化 | 四期长期变化对照，不作为2026年瞬时变形速度 |
| ④ 冰川网格 | 冰厚、坡度、历史速度、候选物源分区 |
| ⑤ 土地覆盖 | 2020—2023年有效覆盖区比较，缺测保留NoData |
| ⑥ 历史沟谷 | 84条沟谷的成因类别分析与地理留出比较 |
| ⑦ DEM | 地形、候选路径、一维与局部二维计算 |
| ⑧ 行政区/人口空间资料 | 区域背景及暴露对象辅助表达 |

第③组只需从官方RAR提取N28E085目标区域的以下四幅栅格，放入`data/raw/unpacked/3_selected/`：

```text
N28E085_2000-01-01_2005-01-01_dhdt_UTM45_UTM45_GCS.tif
N28E085_2005-01-01_2010-01-01_dhdt_UTM45_UTM45_GCS.tif
N28E085_2010-01-01_2015-01-01_dhdt_UTM45_UTM45_GCS.tif
N28E085_2015-01-01_2020-01-01_dhdt_UTM45_UTM45_GCS.tif
```

可使用7-Zip或Bandizip解压。缺少四期栅格时，对应历史高程变化结果不能完整复现；不以其他年份数据替补。

### 2.3 完整计算与分阶段运行

首次运行并尝试获取第二问的历史起报预报：

```bash
python run.py --stage all --official-dir "/path/to/official_archives" --fetch-q2-external
```

已在本地恢复外部预报缓存时可离线运行：

```bash
python run.py --stage all --no-fetch-q2-external
```

各阶段按`prepare → q1 → q2 → q3 → q4`顺序依赖；首次复现建议使用`all`，它还生成联合结果和报告索引。分阶段示例：

```bash
python run.py --stage prepare --official-dir "/path/to/official_archives"
python run.py --stage q1
python run.py --stage q2 --fetch-q2-external
python run.py --stage q3 --no-fetch-q2-external
python run.py --stage q4 --no-fetch-q2-external
python run.py --stage tests
```

`--force-rebuild`使指定的处理缓存失效后重建；输入发生变化时应使用它并重跑依赖阶段。`--config`可指定配置。重算某一题会清理该题的旧成果及对应重绘目录，运行前应将要保留的本地自绘结果另存。`--strict-route`在路径未独立核验时仅允许数值基准；默认允许带明确假设标识的候选路径情景。

Q2是唯一获取外部补充资料的阶段。Q3复用Q2已获得且符合可用时间假设的预报，不自行下载。网络下载失败时程序记录缺失状态，结果可能与本次已下载预报的PNG不同；不能把再分析或未来实测代替历史起报产品。

## 3. 四问如何由代码回答

### 3.1 Q1：天气、边界层及前五项关联因素

主要代码：[prepare_time.py](src/mathorcup/prepare_time.py)、[features.py](src/mathorcup/features.py)、[q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py)。

完整登记63项原始非QC字段：MET 22项、RADM 4项、SOIL 12项、FLUX 25项。对有对应QC列的观测使用其质量标记；没有对应QC的通量字段按缺失值与物理范围检查。仪器诊断、碳通量及重复诊断保留审计记录，不直接参与灾害物理因素排序。

62项物理候选经过样本有效性筛选，本次56项可排序。主窗24小时，同时比较6、72、168小时；历史基准使用2021—2025年实际可用的同期窗口，各变量分别记录有效年份及样本数。

稳健异常与排序为：

$$
Z_j=\frac{X_{j,\mathrm{event}}-\operatorname{median}(X_{j,\mathrm{ref}})}{1.4826\operatorname{MAD}(X_{j,\mathrm{ref}})},\qquad
S_j=0.45\min(|Z_j|/3,1)+0.55(1-p_{j,\mathrm{tail}}).
$$

经验双侧尾部稀有度由历史样本计算；重叠窗口不视为严格独立样本。完整原始排序的前五项与按物理机制和Spearman相关去冗余得到的代表项分别输出。本次原始前五项依次为：通量塔VPD、10m气温、20m气温、地表能量收支余量、2m气温。机制代表项为通量塔VPD、10m气温、地表能量收支余量、下行短波辐射、土壤水分梯度。

边界层分析包含多高度风温、矢量风切变、温度梯度及近似Bulk Richardson数；辐射、显热、潜热、Bowen比和土壤剖面用于解释过程。能量收支余量不直接换算为冰川融水。本问给出单事件条件关联与异常稀有性，缺少多场灾害标签，不能解释为因果重要性。

### 3.2 Q2：物源、触发条件、外部补充与尺度可靠性

主要代码：[q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[extra_official_sources.py](src/mathorcup/extra_official_sources.py)。

以冰厚、坡度和历史速度构造空间物源代理。冰川驱动剪应力近似为：

$$\tau_d=\rho_i gH\sin\theta,\qquad \rho_i=900\;\mathrm{kg/m^3}.$$

静态指数组合归一化剪应力、冰厚、速度和速度梯度；北向坐标递减对速度导数的符号被显式处理，缺失应变保持未知。高分位连通分区用于候选源区，像元冰厚积分只表示冰库存，不直接等于崩塌释放量或下游混合体积。

源区温度采用高差递减近似，考察5、6.5、8℃/km；融化代理使用七日正积温，降雨候选阈值来自非零历史降雨样本。0mm降雨不激活降雨门控，缺测也不记为零风险。组合条件指数属于物源和水热条件评估，缺少岩冰强度、孔压及多事件触发标签，不能当作已验证的力学失稳临界值。

历史84条沟谷的标签为成因类别。模型在同一空间留出划分中比较，并完整展示所有候选模型和判别变量；其ROC-AUC评价成因分类，不评价新灾害是否发生。30/60/120m空间尺度、温度递减率、逐年向前拟合阈值和更新间隔误差分别检查。

外部资料分两类处理：

- [ICIMOD公开遥感分析](https://www.icimod.org/press-releases/second-major-himalayan-flood-in-seven-weeks-raises-new-questions-about-mountain-risk/)在代码中登记为事后发表的Sentinel-1 LOS测量摘要，用于回顾性佐证；它不是已下载的像元位移格网，不据此绘制虚构的位移图，也不作为当时可用的灾前预警输入。
- [Open-Meteo历史单次起报接口](https://open-meteo.com/en/docs/single-runs-api)下载实际起报响应，保留起报时间、下载地址、SHA256和变量缺失状态。当前PNG使用2026-08-25 18:00 UTC起报，按6小时发布延迟假设检查其在北京时间次日09:00前可用；历史实际分发时刻未独立观测。该次数据的零度层高度全部缺失，保持缺失。

本仓库不附下载缓存。缓存路径为`data/external/q2_forecast/`，包括原始响应`ecmwf_source_pre_event_single_run.json`、小时数据`source_forecast_hourly.csv`和`provenance.json`。恢复既有完整包中的这三项可避免下载结果或接口可用性变化；只有重新下载接口且获得同一起报数据，才能相应复现。本项目登记的[事后SAR DEM](https://zenodo.org/records/22842744)只供回顾性几何检查，不用于灾前输入。

### 3.3 Q3：条件时间分布、滚动预警及独立年份检验

主要代码：[q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_feature_bridge.py](src/mathorcup/q3_feature_bridge.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)。

主模型采用预先固定的物理气象变量和Q2条件通道。2026事件筛选出的Q1因素进入解释与敏感性分析，不反向参与2025盲测变量选择。历史划分为2021—2023年拟合、2024年校准、2025年检验；每次验证同步重算Q2阈值，并阻止未来观测或未来缺测插值进入发布时输入。Q1+Q2、仅Q1、仅Q2在相同检验时段比较。

离散危险率、生存率与首次发生时间概率满足：

$$h_k=1-e^{-\lambda_k\Delta t},\quad S_k=\prod_{j\leq k}(1-h_j),\quad P(T=k)=S_{k-1}h_k.$$

条件于未来6小时内发生的时间分布用于计算2.5%—97.5%分位区间；同时报告时域内发生的模型量，不能用条件区间代替无条件概率或实测灾害保障。以半小时为时间单元，离散分位区间保留完整时间格。

预警状态机分别处理进入、维持、续期、低值解除和缺测暂停，首次报警按观测实际可用时刻计算。发布时已解除的历史报警不再计有效提前量。回测同时保留漏报率、误报率、误报比例及各自分母。题面10:30只作灾害/影响参考时刻，未被独立确认为源区触发时刻。

空间警戒范围以Q4主情景实际可达距离约束一维廊道，仅在主情景真正向二维区输送入流时融合其湿区。独立次级二维实验不并入主情景警戒面。本问的2025覆盖率针对气象条件超限代理；真实灾害触发时间的经验覆盖率仍无法估计。

### 3.4 Q4：地形约束传播、1D—2D耦合及敏感性

主要代码：[q4_flowpath.py](src/mathorcup/q4_flowpath.py)、[q4_solver.py](src/mathorcup/q4_solver.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)。

Q4共享Q2高分物源连通区，在官方DEM上构造带逆坡惩罚的候选流路，记录源区、临时口岸坐标、路径逆坡与核验状态。实际源区与河道仍需独立遥感或现场证据，不能把候选地形路径称为已验证真实运动轨迹。核验资料可按`events/verified_route.template.csv`字段另存为`events/verified_route.csv`并配合配置使用。

一维连续方程为$\partial_t A+\partial_sQ=q_e$；二维采用深度平均浅层流守恒方程。静水重构和Rusanov有限体积通量处理地形与守恒；Voellmy型阻力包含库仑项与速度平方阻力。初始流入速度用于物质输入时的动量，不假设所有运动物质在山顶静止释放。UTM北向与栅格行向的动量符号一致，实际每步一维通量及速度传给二维入口，入流速度不人为截断。

基准参数为有效移动方量5000万m³、宽度150m、历时180s、初速度65m/s、$\mu=0.04$、$\xi=2500\;\mathrm{m/s^2}$、等效密度1800kg/m³。它们是配置中的情景假设，需要方量、摩擦、宽度、初速度和释放历时敏感性检查；并非现场实测参数。

本次5000万m³主情景未抵达下游二维断面。代码采用预设8000万m³次级一维敏感性情景的真实通量驱动局部二维实验，并始终单列其来源。未到达情景的到达时间保持缺失；绘图叉号只表示未到达，不能转换为数值到达时间。

图像检查统一地理范围、时间原点与可比色标。二维过程帧按真实入流后时刻选取，避免把尚未入流的空白帧误解为数值故障。入口横向位置与方向动量作2×2对照，侧枝只能解释为边界敏感性线索；不得仅凭模型分枝认定天然分流。60m与120m的比较用于评估尺度差异，仅两档网格不能证明收敛。

动压代理使用$q=\tfrac12\rho u^2$，在流深大于0.2m处记录全时段最大值并转换为kPa。显示色标分位数上限不截断底层数值；它不是计入建筑结构、冲击放大和颗粒相互作用的真实荷载。固定沟宽、陡坡干湿前沿及粗DEM可能导致局部流深偏大，体积守恒通过也不等于真实水深、流路已经验证。

联合模块[q34_joint_arrival.py](src/mathorcup/q34_joint_arrival.py)计算$T_{\mathrm{arrival}}=T_{\mathrm{trigger}}+\tau_{\mathrm{travel}}$的参数情景分布；仅到达情景参与条件到达时间，未到达比例单独报告，不把两个边缘95%区间端点直接相加。该图文件名中的`MATHEMATICAL_DEMO`仍表示其数学情景性质。

## 4. 本次PNG对应的当前计算结论

下列数字摘自本地v2.0重算结果，说明当前提交图像的计算状态；不是重新取得的现场观测。

| 项目 | 当前结果与解释 |
|---|---|
| Q1候选与排序 | 62项候选，56项有效；原始前五项与机制代表项分开 |
| Q2降雨候选阈值 | 非零24h累计降雨历史q90：10.2mm |
| Q2热力候选阈值 | 源区七日正积温暖季q90：39.0546℃·d |
| Q2历史成因分类 | Logistic空间OOF AUC约0.63585；评价沟谷成因类别 |
| Q2外部预报 | 本次为`ACTUAL_DOWNLOADED_ARCHIVED_FORECAST`；零度层仍缺失 |
| Q3发布时间 | 北京时间2026-08-26 09:00 |
| Q3条件95%时间区间 | 同日09:00—15:00；属于触发参数情景 |
| Q3有效预警 | 09:00无有效黄色预警，有效提前量为空 |
| Q3覆盖率 | 2025气象代理检验0.9625；真实灾害时间PICP为空 |
| Q3主情景空间包络 | 约2.1563km²，以主情景可达性约束 |
| Q4主情景最远传播 | 14.52km，未达下游二维断面 |
| Q4主情景1D相对体积残差 | 约1.79×10⁻¹⁵，仅为数值收支检查 |
| Q4二维图的来源 | 8000万m³次级1D情景的通量；不是5000万m³主情景的二维结果 |

旧文档中“预报尚未下载”“物理去冗余代表项就是原始Top5”“09:00仍可引用前一日报警提前量15.5小时”“90.6%代理覆盖率”“约4.09km²主警戒范围”及早期二维面积、时刻等数值均不再作为本版结论。当前判断以相同版本代码、配置和本地重算状态共同确定。

## 5. 图像重绘、图注与辅助脚本

正常绘图通过[common.savefig](src/mathorcup/common.py)保存PNG、SVG，并由[redraw_export.py](src/mathorcup/redraw_export.py)导出面板数据；[redraw_catalog.py](src/mathorcup/redraw_catalog.py)整合导入友好的曲线、单位、栅格、掩膜、坐标、GIS和原始数值表。所有非PNG结果只在本地保留。

如果已恢复完整数值结果，可执行：

```bash
python scripts/export_redraw_data.py
python scripts/check_redraw_retention.py
```

回放脚本在`tmp/redraw_replay_outputs/`隔离计算并校对原结果不变，生成`outputs/redraw_data/00_重绘数据索引.csv`。按原PNG文件名查找面板数据，即可在其他软件重绘。当前完整包曾核查63张图、143个面板、748个逻辑绘图数据集、46组共用栅格和99个栅格图层；这些导出数据没有随本次GitHub提交上传。仅克隆仓库而未恢复原始数据和计算结果时，不能直接运行缓存重绘或完整保留检查。

| 脚本 | 用途与前置条件 |
|---|---|
| `replot_q4_1d_cached.py` | 从本地1D缓存重新绘制传播图 |
| `replot_q4_figures.py` | 从本地2D快照及边界试验缓存重新绘图 |
| `compute_inlet_factorial.py` | 单独运行局部二维入口2×2对照，需官方数据与1D结果 |
| `make_q4_figure_atlas.py` | 汇集已生成的Q4图像；派生图册不提交 |
| `finalize_outputs.py` | 完整本地结果的索引、联合结果和来源哈希刷新 |
| `export_redraw_data.py` / `check_redraw_retention.py` | 导出并检查逐图数据，需完整本地计算环境 |
| `build_answer_pdf.py` | 可选生成本地作答PDF，需完整结果及两个Noto字体 |
| `build_release.py` | 可选打包本地完整数据、结果及PDF，需报告和通过测试的质量记录 |
| `package_redraw_data.py` | 原完整交付包的增补核验工具，另需已有完整ZIP与父目录`tmp/tests_redraw_full.log`，不属于克隆后的主运行流程 |

图题陈述研究对象与指标；图注说明数据时段、站点与源区差异、时间时区、坐标、单位、模型来源和适用范围。Q1标注关联性质；Q2区分历史状态、候选条件与灾后佐证；Q3区分条件区间、气象代理覆盖率与灾害验证；Q4每张二维图应说明次级情景及其参数。地图底图来自官方DEM、冰川网格和沟谷资料，模拟面不能当作独立遥感灾后范围。

## 6. 测试与复现边界

本次完整数据环境测试为60项通过。此次上传范围的干净副本测试为51项通过、9项因缺少官方附件或数值结果而跳过；跳过项需恢复资料并完成计算后再验证。主要检查涵盖QC、派生量单位、非零降雨门控、缺测传播、未来数据隔离、预警解除及暂停、误报分母、地理北向、静水平衡、湿坡加速、1D/2D守恒、未到达绘图、主次情景隔离和重绘数值保留。

这些检查支持实现与数值一致性，不提供独立灾害预测精度。仍需补充源区和流路核验、实际触发时刻、多次灾害及未灾害记录、现场摩阻和断面资料，才能开展事件级校准与长期可靠性评价。[events/verified_events.template.csv](events/verified_events.template.csv)仅提供真实事件登记字段，不含虚构样本。

## 7. 理论与资料出处

- Caine N. (1980), *The rainfall intensity-duration control of shallow landslides and debris flows*：降雨阈值研究背景；本项目不直接移植其固定阈值。
- Hock R. (2003), *Temperature index melt modelling in mountain areas*：正积温方法的物理含义与适用范围。
- Iverson R.M. (1997), *The physics of debris flows*：泥石流动力学与连续介质近似。
- Voellmy A. (1955), *Über die Zerstörungskraft von Lawinen*；另参见[RAMMS官方摩阻参数说明](https://ramms.ch/ramms-debrisflow/friction-parameters/)。
- 官方八组附件为计算主要来源；外部预报、事后遥感摘要及事后DEM的链接与用途见第3.2节。历史文档登记的额外文献和外部图件不作为本版模型验证真值。

## 8. Markdown整合来源

以下11份原说明的有效内容已归并到本文，不另上传重复Markdown。版本与证据范围冲突处按第4节统一：

| 原说明 | 整合位置 |
|---|---|
| `README.md` | 上传范围、运行、结果解释、重绘流程 |
| `DATA3_RAR.md` | 第2.2节四期目标瓦片 |
| `COMPETITION_Q4_MODEL_NOTES.md` | 第3.4节方程、动量输入、未到达情景及局限 |
| `Q4_FIGURES_METHOD_AND_CAPTIONS_v1.2.md` | 第3.4、5节边界试验、时序、图注与重绘 |
| `Q2外部数据获取与运行说明_v1.6.md` / `v1.7.md` | 第2.3、3.2节外部数据和可用时间 |
| `模型理论与赛题逐问验收_v1.6.md` / `v1.7.md` | 第3、4、6、7节四问方法、当前结果与验证边界 |
| `v1.7_四问题目符合性及数据利用说明.md` | 第2、3节数据用途及四问代码对应 |
| `论文图题与图注清单_v1.6.md` / `v1.7.md` | 第5、9节图注原则与完整图像索引 |

## 9. 全部63张PNG图像索引

下列图片链接可直接查看PNG；“本地数值证据”是运行后生成的文件名，用于定位推导与重绘，文件本身不在此次提交中。源码路径均相对仓库根目录。

### 问题一：气象与边界层（16张）

生成代码：[q1_event_signatures.py](src/mathorcup/q1_event_signatures.py)、[task_evidence.py](src/mathorcup/task_evidence.py)。

| PNG | 用途/图注要点 | 本地数值证据 |
|---|---|---|
| [Q1_00_data_coverage.png](outputs/q1/Q1_00_data_coverage.png) | 逐小时观测数据有效率 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_01_event_hourly_timeseries.png](outputs/q1/Q1_01_event_hourly_timeseries.png) | 吉隆泥石流灾前气象要素变化 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_02_boundary_layer_profiles.png](outputs/q1/Q1_02_boundary_layer_profiles.png) | 高寒山区近地层风温垂直结构 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_03_top15_anomaly_ranking.png](outputs/q1/Q1_03_top15_anomaly_ranking.png) | 灾前多源气象因素异常稀有性排序 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_04_multiscale_anomaly_heatmap.png](outputs/q1/Q1_04_multiscale_anomaly_heatmap.png) | 多时间尺度气象异常分布 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_05_reference_year_sensitivity.png](outputs/q1/Q1_05_reference_year_sensitivity.png) | 关键特征排序的年份敏感性 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_06_pseudo_event_controls.png](outputs/q1/Q1_06_pseudo_event_controls.png) | 灾前气象异常与历史同期对照 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_07_county_2005_2025_summer_climate.png](outputs/q1/Q1_07_county_2005_2025_summer_climate.png) | 吉隆及周边地区夏季气温变化 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_08_boundary_layer_seasonal_structure.png](outputs/q1/Q1_08_boundary_layer_seasonal_structure.png) | 近地层稳定结构的季节变化 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_09_radiation_and_turbulent_energy.png](outputs/q1/Q1_09_radiation_and_turbulent_energy.png) | 灾前地表辐射与热量通量变化 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_10_multidepth_soil_profiles.png](outputs/q1/Q1_10_multidepth_soil_profiles.png) | 灾前土壤多深度温湿廓线 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_11_all_candidate_ranking.png](outputs/q1/Q1_11_all_candidate_ranking.png) | 完整候选排序；对应Q1全部因素比较与原始前五项 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_11_part1_complete_ranking.png](outputs/q1/Q1_11_part1_complete_ranking.png) | 完整候选排序；对应Q1全部因素比较与原始前五项 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_11_part2_complete_ranking.png](outputs/q1/Q1_11_part2_complete_ranking.png) | 完整候选排序；对应Q1全部因素比较与原始前五项 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_12_all_factor_multiscale.png](outputs/q1/Q1_12_all_factor_multiscale.png) | 完整多时间尺度异常；对应Q1窗口依赖 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |
| [Q1_13_all_factor_redundancy.png](outputs/q1/Q1_13_all_factor_redundancy.png) | 全部候选冗余结构；区分原始排序与机制代表项 | `Q1_all_candidate_ranking.csv`；`Q1_all_factors_effect_and_stability.csv` |


### 问题二：物源与触发条件（17张）

生成代码：[q2_spatial_mechanism.py](src/mathorcup/q2_spatial_mechanism.py)、[q2_trigger_evaluation.py](src/mathorcup/q2_trigger_evaluation.py)、[q2_external.py](src/mathorcup/q2_external.py)、[task_evidence.py](src/mathorcup/task_evidence.py)。

| PNG | 用途/图注要点 | 本地数值证据 |
|---|---|---|
| [Q2_01_history_gullies_and_oof_roc.png](outputs/q2/Q2_01_history_gullies_and_oof_roc.png) | 历史沟谷成因分布与判别模型ROC曲线 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_02_spatial_holdout_auc.png](outputs/q2/Q2_02_spatial_holdout_auc.png) | 灾种判别模型空间交叉验证 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_03_elevation_drop_distributions.png](outputs/q2/Q2_03_elevation_drop_distributions.png) | 不同成因沟谷高差分布 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_04_climate_condition_timeseries.png](outputs/q2/Q2_04_climate_condition_timeseries.png) | 降雨与融水触发条件时间变化 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_05_glacier_physical_indicators.png](outputs/q2/Q2_05_glacier_physical_indicators.png) | 冰川厚度、坡度与驱动剪应力分布 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_06_dem_landcover_history_overlay.png](outputs/q2/Q2_06_dem_landcover_history_overlay.png) | 研究区地形与土地覆盖分布 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_07_glacier_historic_thinning_optional.png](outputs/q2/Q2_07_glacier_historic_thinning_optional.png) | 冰川历史高程变化空间分布 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_08_trigger_conditions_time_series.png](outputs/q2/Q2_08_trigger_conditions_time_series.png) | 降雨与融水触发条件变化 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_09_source_spatial_condition_index.png](outputs/q2/Q2_09_source_spatial_condition_index.png) | 冰川源区地形物源条件分布 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_10_spatial_resolution_reliability.png](outputs/q2/Q2_10_spatial_resolution_reliability.png) | 触发评估空间尺度可靠性 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_11_temperature_lapse_sensitivity.png](outputs/q2/Q2_11_temperature_lapse_sensitivity.png) | 源区热力条件递减率敏感性 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_12_forecast_weather_conditions.png](outputs/q2/Q2_12_forecast_weather_conditions.png) | 源区历史起报天气预报条件 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_13_glacier_velocity_and_release_candidates.png](outputs/q2/Q2_13_glacier_velocity_and_release_candidates.png) | 冰川历史运动与候选失稳区空间特征 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_14_landcover_2020_2023_comparison.png](outputs/q2/Q2_14_landcover_2020_2023_comparison.png) | 2020年与2023年土地覆盖对比 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_15_all_model_comparison.png](outputs/q2/Q2_15_all_model_comparison.png) | 全部候选模型的同一空间留出比较；非触发事件AUC | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_16_all_predictor_distributions.png](outputs/q2/Q2_16_all_predictor_distributions.png) | 全部五项判别变量分布与解释 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |
| [Q2_17_forward_time_reliability.png](outputs/q2/Q2_17_forward_time_reliability.png) | 逐年向前阈值与不同更新间隔的监测误差 | `Q2_primary_trigger_assessment.json`；`Q2_all_models_comparison.csv` |


### 问题三：预警时间与范围（11张）

生成代码：[q3_task_aligned.py](src/mathorcup/q3_task_aligned.py)、[warning_validation.py](src/mathorcup/warning_validation.py)、[q3_spatial_footprint.py](src/mathorcup/q3_spatial_footprint.py)、[forecast_sensitivity.py](src/mathorcup/forecast_sensitivity.py)。

| PNG | 用途/图注要点 | 本地数值证据 |
|---|---|---|
| [Q3_01_precursor_and_scenario_hazard.png](outputs/q3/Q3_01_precursor_and_scenario_hazard.png) | 灾前气象条件及触发风险演变 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_02_95percent_conditional_time_interval.png](outputs/q3/Q3_02_95percent_conditional_time_interval.png) | 泥石流触发时间分布与95%预测区间 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_03_warning_lead_and_backtest.png](outputs/q3/Q3_03_warning_lead_and_backtest.png) | 泥石流预警时间与气象条件回测 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_04_parameter_sensitivity.png](outputs/q3/Q3_04_parameter_sensitivity.png) | 预警模型参数敏感性分析 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_05_spatial_warning_distance.png](outputs/q3/Q3_05_spatial_warning_distance.png) | 泥石流传播廊道与空间预警距离 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_06_validation_metrics.png](outputs/q3/Q3_06_validation_metrics.png) | 气象条件时间区间与预警指标验证 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_07_Q1_feature_inheritance.png](outputs/q3/Q3_07_Q1_feature_inheritance.png) | Q1事件因素传递及解释性敏感性 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_07_feature_group_ablation.png](outputs/q3/Q3_07_feature_group_ablation.png) | Q1气象特征与Q2物理条件的贡献比较 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_08_full_forward_validation.png](outputs/q3/Q3_08_full_forward_validation.png) | 完整2025检验期，三组输入、区间与覆盖过程 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_09_warning_episodes.png](outputs/q3/Q3_09_warning_episodes.png) | 预警进入、维持、续期、解除和有效期 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |
| [Q3_10_archived_forecast_sensitivity.png](outputs/q3/Q3_10_archived_forecast_sensitivity.png) | 实际起报天气预报进入Q2条件后的未来驱动敏感性 | `Q3_primary_competition_answer.json`；`Q3_all_models_forward_backtest.csv`；`Q3_warning_state_hourly.csv` |


### 问题四：动力学与敏感性（18张）

生成代码：[q4_solver.py](src/mathorcup/q4_solver.py)、[q4_visuals.py](src/mathorcup/q4_visuals.py)、[q4_grid_review.py](src/mathorcup/q4_grid_review.py)、[publication_assets.py](src/mathorcup/publication_assets.py)。

| PNG | 用途/图注要点 | 本地数值证据 |
|---|---|---|
| [Q4_00B_local_true_DEM_detail.png](outputs/q4/Q4_00B_local_true_DEM_detail.png) | 下游二维计算区域地形图 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_00_actual_DEM_geographic_context.png](outputs/q4/Q4_00_actual_DEM_geographic_context.png) | 吉隆冰崩—泥石流研究区地形图 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_01_route_geolocation_quality_gate.png](outputs/q4/Q4_01_route_geolocation_quality_gate.png) | 冰崩—泥石流传播路径及地形纵剖面 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_02_1d_runout_timespace.png](outputs/q4/Q4_02_1d_runout_timespace.png) | 一维泥石流传播及地形纵剖面 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_03_1d_hydrograph_and_mass_conservation.png](outputs/q4/Q4_03_1d_hydrograph_and_mass_conservation.png) | 耦合断面流量与体积守恒检验 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_04_1d_12frames.png](outputs/q4/Q4_04_1d_12frames.png) | 一维泥石流流深演进过程 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_05_local_2d_depth_and_arrival.png](outputs/q4/Q4_05_local_2d_depth_and_arrival.png) | 泥石流二维最大流深与到达时间；须注明次级情景来源 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_06A_onset_and_wet_area_diagnostic.png](outputs/q4/Q4_06A_onset_and_wet_area_diagnostic.png) | 耦合断面入流及淹没面积变化；须注明次级情景来源 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_06B_inlet_position_ablation.png](outputs/q4/Q4_06B_inlet_position_ablation.png) | 二维入流位置对传播范围的影响；须注明次级情景来源 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_06C_2x2_inlet_boundary_factorial.png](outputs/q4/Q4_06C_2x2_inlet_boundary_factorial.png) | 二维入流位置与动量组合试验；须注明次级情景来源 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_06_local_2d_10frames.png](outputs/q4/Q4_06_local_2d_10frames.png) | 泥石流二维流深演进过程；须注明次级情景来源 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_07_parameter_sensitivity.png](outputs/q4/Q4_07_parameter_sensitivity.png) | 泥石流动力学参数敏感性分析 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_08_volume_friction_phase_diagram.png](outputs/q4/Q4_08_volume_friction_phase_diagram.png) | 物源体积与摩擦系数对传播的影响 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_09_grid_resolution_reliability.png](outputs/q4/Q4_09_grid_resolution_reliability.png) | 二维模拟网格分辨率敏感性；须注明次级情景来源 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_10_impact_pressure_distribution.png](outputs/q4/Q4_10_impact_pressure_distribution.png) | 二维流动峰值动压分布；须注明次级情景来源 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_11_six_key_frames.png](outputs/q4/Q4_11_six_key_frames.png) | 六个关键传播阶段；共同色标，来自本次二维数组；须注明次级情景来源 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_12_secondary_inflow_transfer.png](outputs/q4/Q4_12_secondary_inflow_transfer.png) | 次级预设方量情景的真实1D通量、速度及累计量；须注明次级情景来源 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |
| [Q4_13_1d_grid_resolution.png](outputs/q4/Q4_13_1d_grid_resolution.png) | 一维网格分辨率比较 | `Q4_1d_settings_and_summary.json`；`Q4_2d_status.json`；`Q4_parameter_sensitivity.csv` |


### 联合时间情景（1张）

生成代码：[q34_joint_arrival.py](src/mathorcup/q34_joint_arrival.py)。

| PNG | 用途/图注要点 | 本地数值证据 |
|---|---|---|
| [Q34_joint_uncertainty_MATHEMATICAL_DEMO.png](outputs/joint/Q34_joint_uncertainty_MATHEMATICAL_DEMO.png) | 触发时间与下游到达时间联合分布 | `Q34_status.json` |

