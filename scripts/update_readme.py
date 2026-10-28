"""Build the single project README from current experiment results and figure index."""
from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs'
def j(path):return json.loads((OUT/path).read_text())

def main():
    q1=j('q1/Q1_direct_answer.json');c1=j('q1/Q1_complete_factor_contract.json');q2=j('q2/Q2_expanded_mechanism_answer.json')
    q3=j('q3/Q3_primary_competition_answer.json');m=q3['proxy_backtest'];lead=j('q3/Q3_warning_lead_answer.json')
    q4=j('q4/Q4_cascade_reconstruction.json');geo=j('q4/Q4_reconstruction_geometry.json');two=j('q4/Q4_reconstruction_2d_summary.json')
    space=j('q3/Q3_reconstructed_spatial_answer.json');deposition=j('q4/Q4_confluence_deposition_geometry.json')
    rows=list(csv.DictReader((OUT/'FIGURE_INDEX.csv').open(encoding='utf-8-sig')))
    energy=next(x for x in q2['channels'] if x['channel']=='station_positive_energy_residual_7d_mmwe')
    vpd=next(x for x in q1['extreme_summary'] if x['feature']=='vpd')
    temp=next(x for x in q1['extreme_summary'] if x['feature']=='temp_2m')
    limits=j('q2/Q2_primary_trigger_assessment.json')['thresholds']
    advance=lead['historical_new_onset_lead']
    numerical=list(csv.DictReader((OUT/'q4/Q4_reconstruction_1d_grid_experiment.csv').open()))
    by_dx={float(r['dx_m']):r for r in numerical}
    fine,base=by_dx[min(by_dx)],by_dx[sorted(by_dx)[1]]
    flow_change=100*(float(fine['peak_discharge_m3_s'])/float(base['peak_discharge_m3_s'])-1)
    time_change=100*(float(fine['arrival_s'])/float(base['arrival_s'])-1)
    area_change=100*abs(two['grid_experiments'][1]['wet_area_km2']/two['grid_experiments'][0]['wet_area_km2']-1)
    review=j('quality/review_verdicts.json');checks=j('quality/test_result.json')
    interval='至'.join(x for x in q3['conditional_95pct_interval_times_bjt'] if x is not None) or '该时域无条件超限轨迹'
    forecast_time=q3['forecast_warning_time_bjt'] or '未达到预测预警阈值'
    forecast_lead=q3['forecast_warning_lead_h_to_task_reference']
    boundary='触及搜索边界' if q4['parameter_bound_diagnostics']['at_search_bound'] else '位于搜索范围内'
    depth_change=100*abs(float(fine['peak_depth_m'])/float(base['peak_depth_m'])-1)
    document=f'''# MathorCup 2026赛道B：四问代码、方法与PNG成果（v2.2）

本仓库用于妈妈杯代码图片上传保存留档，以确保能够正确回退版本。本文整合数据准备、四问作答、模型原理、运行方法、图题图注与全部图像索引。当前成果由实际输入、分析计算和数值实验生成；参数的观测、文献、标定和实验设计角色分别记录。

## 1. 当前计算回答了什么

| 问题 | 计算结果 |
|---|---|
| Q1关联特征 | {c1['predefined']}项候选、{c1['eligible']}项有效；按机制合并相关变量后，五项代表通道为水汽VPD、近地层热力、地表能量余量、短波辐射、土壤水分。原始逐变量排序同时保留。 |
| Q1极端过程 | 事件前七日VPD有{vpd['hours_above_matched_q95']}小时达到匹配历史q95，其中{vpd['new_matched_record_hours']}小时超过匹配历史记录；2m气温有{temp['hours_above_matched_q95']}小时达到q95，最长连续{temp['longest_consecutive_q95_hours']}小时。 |
| Q2临界条件 | 非零24h降雨q90为{limits['rain_24h_nonzero_q90_mm']:.2f}mm，源区7日正积温q90为{limits['source_pdd7_warmseason_q90_Cday']:.2f}℃·d；新增30日热量、辐射能量和冻融通道。 |
| Q2本次能量通道 | 正能量余量7日当量为{energy['value_at_issue']:.2f}mm，历史q90为{energy['historical_q90']:.2f}mm，达到历史第{100*energy['historical_percentile']:.1f}百分位；现已进入主动态指数。 |
| Q3预测预警 | 6h条件超限频率{q3['predicted_condition_cross_probability']:.3f}；按0.25阈值，预测预警时间{forecast_time}，至题面参考时刻提前{forecast_lead}h。 |
| Q3监测预警 | 有效监测段起点{q3['monitoring_episode_start_bjt']}，其较题面参考时刻的回放提前量{q3['first_monitoring_warning_advance_h_to_reference']}h；与09:00预测发布分别记录。 |
| Q3区间与检验 | 条件首次超限95%区间{interval}；2025年相似日基准检验{m['n_backtest']}次发布、{m['n_proxy_crossings']}次气象超限，校准后覆盖率{100*m['proxy_interval_coverage']:.2f}%，平均宽度{m['mean_interval_width_h']:.2f}h；TP/FP/FN/TN={m['TP']}/{m['FP']}/{m['FN']}/{m['TN']}。 |
| Q3提前量 | {advance['n']}次正确预警的新发超限中，提前量中位数{advance['median_h']:.2f}h，四分位范围{advance['q25_h']:.2f}至{advance['q75_h']:.2f}h。 |
| Q3空间范围 | 源区至口岸模型廊道{space['source_to_port_corridor_km']:.2f}km，变断面与二维湿区融合包络{space['warning_envelope_km2']:.2f}km²。 |
| Q4分阶段复现 | 官方DEM汇流路径长{geo['length_m']/1000:.2f}km、落差{geo['drop_m']:.0f}m；初始方量103百万m³，可动输入36百万m³，按实际求解通量连接上下游。 |
| Q4口岸结果 | 源区释放后{q4['downstream']['arrival_s']:.1f}s到达口岸河段，峰值流量{q4['downstream']['peak_discharge_m3_s']:.0f}m³/s，对应断面平均速度{q4['downstream']['velocity_at_peak_discharge_m_s']:.2f}m/s。 |
| Q4二维网格试验 | 30m与60m累计湿润面积分别{two['grid_experiments'][0]['wet_area_km2']:.4f}、{two['grid_experiments'][1]['wet_area_km2']:.4f}km²；每档输入体积均与一维累计通量一致。 |

预测预警、持续监测状态和事后过程标定分别报告。95%区间对应条件首次超限，历史覆盖率评价实际气象超限；当前历史起报驱动的回放与相似日基准检验分别列示。校准前覆盖率{m['raw_proxy_interval_coverage']:.2%}；校准后平均区间比完整6h窗口缩短{m['mean_width_reduction_vs_full_horizon']:.2%}，同时报告宽度与覆盖率。模型约定集中于第6节。

## 2. 仓库上传范围

- 计算入口、全部模型代码、配置与测试：[run.py](run.py)、[src/mathorcup](src/mathorcup)、[configs](configs)、[scripts](scripts)、[tests](tests)。
- 输入模板：[events](events)，两个仅含表头的CSV用于登记核验事件与路径。
- 成果目录仅提交{len(rows)}张PNG，分布于`outputs/q1`、`q2`、`q3`、`q4`、`joint`。
- `.gitignore`使生成的CSV、JSON、NPZ、GIS、SVG、GIF、PDF、图廊、缓存和重绘副本留在本地。原始附件、下载产品和字体按下一节恢复。

`configs/case_evidence.yaml`是带出处和角色的事件约束；`configs/calibrated_q4.yaml`保存本次数值实验得到的有效参数与输入指纹，属于复现模型配置。

## 3. 环境、数据与复现

推荐Python 3.12。在仓库根目录执行：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Windows使用`.venv\\Scripts\\activate`激活。中文绘图采用Noto Sans CJK SC，可安装系统字体，或从[Noto CJK官方目录](https://github.com/notofonts/noto-cjk/tree/main/Sans/OTF/SimplifiedChinese)取得`NotoSansCJKsc-Regular.otf`、`NotoSansCJKsc-Bold.otf`放入`assets/fonts/`；PDF生成使用这两个字体。

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
7. **数值与地形检查**：上游15m、下游{q4['downstream']['dx_m']:g}m主网格，另以{float(fine['dx_m']):g}/{float(base['dx_m']):g}/{max(by_dx):g}m比较下游分辨率；二维30/60m独立计算。标定在15m网格进行，固定所得参数后在更细输出网格评价，不将加密结果写成重新标定。网格试验保持同一入流和物理端点，分别检查守恒、压力梯度、干湿前缘及分辨率响应。

本次下游标定参数为μ={q4['stage2_parameters']['mu']:.6f}、ξ={q4['stage2_parameters']['xi']:.2f}m/s²、有效断面尺度{q4['stage2_parameters']['section_capacity_scale']:.3f}（{boundary}）。全程体积相对残差{q4['cascade_volume_residual']:.2e}。`Q4_14`至`Q4_25`为分阶段计算与复核图；早期编号为同一官方数据上的单阶段参考实验。

## 5. 图表数据、重绘与测试

每次`common.savefig`自动导出实际曲线、柱形、散点、区间、多边形、栅格、掩膜和坐标，保存在`outputs/redraw_data/`。首先打开`00_重绘数据索引.csv`，按PNG文件名查找数据；地图另有GeoTIFF，原始结果表在`source_tables/`。这些本地文件支持Origin、Excel或GIS自行重绘。

```bash
python scripts/replot_reconstruction.py
python scripts/check_redraw_retention.py
python scripts/export_redraw_data.py
```

第一条从缓存数组重绘分阶段结果；第二条逐图校对数值保留；第三条在隔离目录回放完整绘图并检查原结果保持一致。缓存重绘和隔离回放均需先准备完整本地输入与数值结果。`replot_q4_figures.py`、`replot_q4_1d_cached.py`对应参考单阶段图。`build_answer_pdf.py`直接读取当前表与图，README由本脚本按当前结果更新。

本次{checks.get('passed','见记录')}项测试通过。测试覆盖QC和单位、未来信息隔离、预警状态、误报分母、静水平衡、真实压力梯度、干湿前缘、变断面守恒、能量范围、路径保护与累计体积耦合。解析小网格只用于验证求解器，正式个案输入与输出由来源清单单独审计。

## 6. 模型约定与尺度解释

| 环节 | 统一解释 |
|---|---|
| 关联与阈值 | Q1评分描述本次事件相对历史的异常关联；Q2 q90是条件评估阈值。历史沟谷AUC的目标是成因类别。 |
| 跨境推算链 | 站点4140m、源区代表高程约5297m，以6.5℃/km主递减率并比较5与8℃/km；实际预报格点偏移约8.38km。站点能量按站点尺度解释。 |
| 时间与覆盖率 | 95%区间是条件首次超限的时间预测区间，经验覆盖率由历史气象超限事件计算；题面10:30作为影响参考时刻。 |
| 标定资料与空间 | Q4结合事后事件研究进行过程复现；Q3天气输入与Q4事后标定资料分别记录，空间范围属于模型回放包络。 |
| 断面和网格 | 有效断面尺度{q4['stage2_parameters']['section_capacity_scale']:.3f}，{boundary}；一维{float(base['dx_m']):g}m至{float(fine['dx_m']):g}m的峰值流量变化{flow_change:+.1f}%，到达时间变化{time_change:+.1f}%，最大流深变化{depth_change:.1f}%；二维面积变化{area_change:.1f}%。质量闭合与各量的分辨率响应分别报告。 |
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
'''
    for item in review['findings']:
        document+=f"| {item['id']} | {item['finding']} | {item['judgement']} | {item['evidence']} {item['action']} |\n"
    document+='''
## 9. 全部PNG的用途与源码索引

下表链接为当前PNG；“数值证据”列给出本地计算结果文件名。这些表和缓存按本次上传范围留在本地。
'''
    sections=[('Q1','问题一'),('Q2','问题二'),('Q3','问题三'),('Q4','问题四'),('JOINT','触发至到达联合计算')]
    for group,title in sections:
        selected=[r for r in rows if r['question']==group]
        if group=='Q4':selected=sorted(selected,key=lambda r:(0 if any(r['figure'].split('/')[-1].startswith(f'Q4_{i}') for i in range(14,26)) else 1,r['figure']))
        document+=f'\n### {title}（{len(selected)}张）\n\n| PNG | 用途 | 生成代码 | 数值证据 |\n|---|---|---|---|\n'
        for r in selected:
            paths=[x if '/' in x else 'src/mathorcup/'+x for x in r['generating_code'].split(';')]
            code='、'.join(f'[{Path(x).name}]({x})' for x in paths)
            evidence='；'.join('`'+x+'`' for x in r['numerical_evidence'].split(';'))
            role='参考单阶段实验；' if group=='Q4' and not any(Path(r['figure']).name.startswith(f'Q4_{i}') for i in range(14,26)) else ''
            document+=f"| [{Path(r['figure']).name}](outputs/{r['figure']}) | {role+r['purpose']} | {code} | {evidence} |\n"
    (ROOT/'README.md').write_text(document,encoding='utf-8')
    print({'README':'README.md','figures':len(rows)})

if __name__=='__main__':main()
