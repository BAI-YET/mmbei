#!/usr/bin/env python3
"""Build the answer report solely from current model outputs and their figures."""
from pathlib import Path
import json,sys,math,hashlib,html
import numpy as np
import pandas as pd
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4,landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont as ReportFont
from reportlab.platypus import BaseDocTemplate,PageTemplate,Frame,Paragraph,Spacer,Image,Table,TableStyle,PageBreak,NextPageTemplate,KeepTogether
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs';sys.path.insert(0,str(ROOT/'src'))
from mathorcup.common import LABELS

def j(path):return json.loads((OUT/path).read_text(encoding='utf-8'))
def read(path):return pd.read_csv(OUT/path)
def pct(x):return '--' if x is None or pd.isna(x) else f'{x:.1%}'
def val(x,nd=2):return '--' if x is None or pd.isna(x) else f'{x:,.{nd}f}'
def E(t):return html.escape(str(t))

def make_font(src,dst,charset):
    """Convert only the needed CFF outlines to embedded TrueType outlines."""
    from fontTools.ttLib import TTFont
    from fontTools import subset
    from fontTools.fontBuilder import FontBuilder
    from fontTools.pens.ttGlyphPen import TTGlyphPen
    from fontTools.pens.cu2quPen import Cu2QuPen
    f=TTFont(src);op=subset.Options();op.recalc_timestamp=False
    sub=subset.Subsetter(options=op);sub.populate(unicodes={ord(c) for c in charset});sub.subset(f)
    glyphset=f.getGlyphSet();order=f.getGlyphOrder();fb=FontBuilder(f['head'].unitsPerEm,isTTF=True);fb.setupGlyphOrder(order)
    glyphs={}
    for name in order:
        pen=TTGlyphPen(None);glyphset[name].draw(Cu2QuPen(pen,1.,reverse_direction=True));glyphs[name]=pen.glyph()
    fb.setupGlyf(glyphs);fb.setupHorizontalMetrics(f['hmtx'].metrics)
    fb.setupHorizontalHeader(ascent=f['hhea'].ascent,descent=f['hhea'].descent)
    fb.setupCharacterMap(f.getBestCmap());fb.setupNameTable({'familyName':'NotoReport','styleName':'Regular','uniqueFontIdentifier':dst.stem,'fullName':dst.stem,'psName':dst.stem})
    fb.setupOS2(sTypoAscender=f['hhea'].ascent,sTypoDescender=f['hhea'].descent,usWinAscent=max(0,f['hhea'].ascent),usWinDescent=abs(f['hhea'].descent))
    fb.setupPost();fb.setupMaxp();fb.save(dst)

TMP=ROOT/'tmp/report';TMP.mkdir(parents=True,exist_ok=True)
charset=Path(__file__).read_text()+json.dumps(LABELS,ensure_ascii=False)+''.join(chr(i) for i in range(32,256))+'αβγδλμπρτξσΔΣ∂∇≤≥≠≈×√²³⁴⁵⁶⁷⁸⁹⁰−°'
for text_source in [ROOT/'README.md',OUT/'quality/review_verdicts.json']:
    if text_source.exists():charset+=text_source.read_text(encoding='utf-8')
for name,source in [('NR','NotoSansCJKsc-Regular.otf'),('NB','NotoSansCJKsc-Bold.otf')]:
    target=TMP/f'{name}.ttf';make_font(ROOT/'assets/fonts'/source,target,charset)
    pdfmetrics.registerFont(ReportFont(name,str(target)))
BLUE=colors.HexColor('#164C68');TEAL=colors.HexColor('#207F7A');INK=colors.HexColor('#253641');GRAY=colors.HexColor('#5D6C77')
styles={
 'body':ParagraphStyle('body',fontName='NR',fontSize=10,leading=16.5,textColor=INK,wordWrap='CJK',spaceAfter=8),
 'small':ParagraphStyle('small',fontName='NR',fontSize=8.5,leading=13,textColor=GRAY,wordWrap='CJK',spaceAfter=5),
 'cap':ParagraphStyle('cap',fontName='NR',fontSize=8.3,leading=12.3,textColor=GRAY,wordWrap='CJK',spaceAfter=9),
 'h1':ParagraphStyle('h1',fontName='NB',fontSize=19,leading=27,textColor=BLUE,wordWrap='CJK',spaceAfter=13),
 'h2':ParagraphStyle('h2',fontName='NB',fontSize=12.5,leading=19,textColor=TEAL,wordWrap='CJK',spaceBefore=4,spaceAfter=8),
 'cover':ParagraphStyle('cover',fontName='NB',fontSize=28,leading=40,textColor=BLUE,wordWrap='CJK',spaceAfter=20),
 'table':ParagraphStyle('table',fontName='NR',fontSize=8.8,leading=13.2,textColor=INK,wordWrap='CJK'),
 'thead':ParagraphStyle('thead',fontName='NB',fontSize=9,leading=13,textColor=colors.white,wordWrap='CJK'),
 'formula':ParagraphStyle('formula',fontName='NR',fontSize=10,leading=17,textColor=BLUE,wordWrap='CJK',backColor=colors.HexColor('#F0F5F7'),borderPadding=9,spaceBefore=3,spaceAfter=12)
}
STORY=[];PAGE_TITLES=[];curw=174*mm;first=True

def P(text,kind='body'):
    STORY.append(Paragraph(text,styles[kind]))

def page(title,kicker,wide=False):
    global first,curw
    if not first:STORY.extend([NextPageTemplate('W' if wide else 'P'),PageBreak()])
    first=False;curw=(261 if wide else 174)*mm
    PAGE_TITLES.append(title);P(kicker,'small');P(title,'h1')

def table(headers,rows,widths=None):
    data=[[Paragraph(E(x),styles['thead']) for x in headers]]+[[Paragraph(E(x),styles['table']) for x in row] for row in rows]
    widths=[x*curw for x in widths] if widths else [curw/len(headers)]*len(headers)
    t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),BLUE),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#F0F5F7'),colors.white]),('LINEBELOW',(0,0),(-1,0),.5,BLUE),('LINEBELOW',(0,-1),(-1,-1),.5,colors.HexColor('#CCD7DD'))]))
    STORY.extend([t,Spacer(1,10)])

def image(path,caption,maxh=145*mm,width=None):
    p=OUT/path
    if not p.exists():raise FileNotFoundError(f'Required computed figure missing: {p}')
    w,h=PILImage.open(p).size;scale=min((width or curw)/w,maxh/h)
    STORY.append(Image(str(p),width=w*scale,height=h*scale,hAlign='CENTER'))
    P(caption,'cap')

def foot(canvas,doc):
    w,h=canvas._pagesize;canvas.saveState();canvas.setStrokeColor(colors.HexColor('#D5DEE3'));canvas.setLineWidth(.5);canvas.line(18*mm,16*mm,w-18*mm,16*mm)
    canvas.setFont('NR',7.5);canvas.setFillColor(GRAY);canvas.drawString(18*mm,11.5*mm,'2026 MathorCup B | v2.2 数据分析与四问作答');canvas.drawRightString(w-18*mm,11.5*mm,str(doc.page));canvas.restoreState()

q1=read('q1/Q1_answer_mechanism_top5.csv');raw1=read('q1/Q1_top5_with_stability.csv')
answer1=j('q1/Q1_direct_answer.json');extreme=read('q1/Q1_explicit_extreme_weather.csv')
q2=j('q2/Q2_primary_trigger_assessment.json');limits=q2['thresholds'];extra2=j('q2/Q2_expanded_mechanism_answer.json')
remote2=j('q2/Q2_external_data_provenance.json')['satellite_summary'];models=read('q2/Q2_all_models_comparison.csv');forecast=j('q2/Q2_archived_forecast_quality.json')
q3=j('q3/Q3_primary_competition_answer.json');metrics=q3['proxy_backtest'];lead=j('q3/Q3_warning_lead_answer.json')
space=j('q3/Q3_reconstructed_spatial_answer.json')
q4=j('q4/Q4_cascade_reconstruction.json');geometry=j('q4/Q4_reconstruction_geometry.json')
two=j('q4/Q4_reconstruction_2d_summary.json');deposition=j('q4/Q4_confluence_deposition_geometry.json');grid=read('q4/Q4_reconstruction_1d_grid_experiment.csv')
observation=q4['published_constraints'];down=q4['downstream'];up=q4['upstream'];par=q4['stage2_parameters']
fine_grid,primary_grid=grid.sort_values('dx_m').iloc[:2].to_dict('records')
review=j('quality/review_verdicts.json');residual=read('q4/Q4_calibration_residuals.csv')
time_interval=' 至 '.join(x for x in q3['conditional_95pct_interval_times_bjt'] if x is not None) or '该时域无条件超限轨迹'
forecast_lead=q3['forecast_warning_lead_h_to_task_reference']
boundary='触及搜索边界' if q4['parameter_bound_diagnostics']['at_search_bound'] else '位于搜索范围内'

page('吉隆泥石流灾害：四问作答','2026 MathorCup 大数据竞赛 / 赛道B / v2.2')
STORY.append(Spacer(1,8));P('气象关联、风险预警<br/>与灾害链过程仿真','cover')
P('本文围绕四道题依次建立气象特征分析、触发条件评估、时间及空间预警、分阶段传播仿真。观测数据取自八类官方附件，第二问使用历史起报天气预报；第四问结合事件研究开展事后标定。正文中的数值与图像由本版代码计算并自动读入。')
table(['题目','核心作答'],[
 ['问题一',f'在62项候选、56项有效特征的比较中，将同机制变量合并，得到水汽、近地层热力、地表能量余量、短波辐射及土壤水分五项互补通道，并识别事件前极端小时。'],
 ['问题二',f'建立物源、水源和能量共同评估模型；七日正能量余量当量为{extra2["channels"][3]["value_at_issue"]:.2f} mm，位于历史第{extra2["channels"][3]["historical_percentile"]*100:.1f}百分位。'],
 ['问题三',f'预测预警较题面10:30提前{val(forecast_lead,1)}小时；条件时间区间{time_interval}。历史气象检验覆盖率{pct(metrics["proxy_interval_coverage"])}，平均区间宽{metrics["mean_interval_width_h"]:.2f}小时。'],
 ['问题四',f'在{geometry["length_m"]/1000:.2f} km地形河道上进行分阶段计算，传递{q4["mobile_input_derived_m3"]/1e6:.1f}百万m³可动体积，模拟口岸流量、流速、流深及二维传播。']],[.16,.84])
P('报告日期：2026年10月28日。时间统一为北京时间，空间计算采用UTM 45N米制坐标。','small')

page('数据与计算路线','数据 / 从附件到四问结果')
table(['数据','进入模型的内容'],[
 ['① 逐小时观测','MET多层风温湿、RADM辐射、SOIL多深度温湿、FLUX湍流与能量通量；依质量标记、缺测及物理范围筛查。'],
 ['② 区域气温','2005至2025年县级日温，提供长期气候背景。'],
 ['③ 冰川高程变化','2000至2020年四期目标区域高程变化，描述长期冰川演化。'],
 ['④ 冰川网格','高程、坡度、冰厚、历史速度及冰川掩膜，构造物源条件并定位释放网格。'],
 ['⑤至⑧ 空间资料','土地覆盖、历史沟谷、DEM及行政人口资料，服务于空间条件、分类比较、路径和仿真环境。'],
 ['补充资料','ECMWF历史起报响应及来源记录；第四问的事件研究方量、传播时间和断面结果作为标定约束。']],[.25,.75])
P('Q1以2021至2025年同期窗口构造稳健历史基准。Q3将2021至2023年用于拟合、2024年用于区间校准、2025年用于检验；本次发布仅使用09:00之前已经取得的观测及满足起报可用时间假设的天气预报。Q4先从官方DEM计算河道和横断面，再进行分阶段守恒求解与局部二维计算。')
P('数据来源分为官方观测、文献重建值、模型标定参数、数值设置四类。原始文件、预报响应、配置和源码均保存哈希；图像的面板数值由绘图函数同步导出。')

page('问题一：边界层结构与天气过程','Q1 / 多高度结构与灾前演变')
P('以MET的多层风温观测描述近地层垂直结构，计算温度梯度、矢量风切变及稳定度代理；将RADM与FLUX联系起来分析辐射吸收和湍流交换，再通过SOIL剖面检验地表与浅层土壤的响应。事件前七日序列用于识别天气过程，历史同季同小时数据用于评价其异常程度。')
image('q1/Q1_02_boundary_layer_profiles.png','图1  多高度气温和风场构成近地层结构的直接观测依据。',maxh=84*mm)
image('q1/Q1_01_event_hourly_timeseries.png','图2  灾前逐小时气象过程；所有时刻采用同一北京时间时钟。',maxh=94*mm)

page('问题一：五项因素与极端事件','Q1 / 合并规则、排序与量化结果')
P('对24小时主窗口计算MAD稳健异常与经验尾部稀有度。先登记完整原始排序，再把多高度气温归入同一热力机制；同机制或历史同期Spearman相关绝对值达到0.85的变量不重复入选。机制代表项沿原始评分由高到低选择，因此同时保留题目要求的排序与互补解释。')
P('Z=(事件窗口值-历史中位数)/(1.4826×MAD)；S=0.45 min(|Z|/3,1)+0.55(1-p尾部)。','formula')
table(['顺序','机制代表项','关联评分'],[[str(int(r.answer_rank)),r.label,f'{r.rarity_anomaly_score:.3f}'] for r in q1.itertuples()],[.12,.63,.25])
image('q1/Q1_14_mechanisms_and_extreme_weather.png','图3  左图为五项互补机制；右图逐小时标出同季同小时历史百分位。',maxh=86*mm)
vpd=extreme.loc[extreme.feature.eq('vpd')].iloc[0];temp=extreme.loc[extreme.feature.eq('temp_2m')].iloc[0]
P(f'事件前七日，通量塔VPD有{int(vpd.hours_above_matched_q95)}小时达到历史95%分位，其中{int(vpd.new_matched_record_hours)}小时超过匹配历史记录；2m气温有{int(temp.hours_above_matched_q95)}小时达到95%分位，最长持续{int(temp.longest_consecutive_q95_hours)}小时。原始前五项的逐高度排序另随数值结果保留。')

page('问题一：能量交换与土壤响应','Q1 / 四类观测的物理联系')
P('辐射决定地表可接收的能量，显热和潜热反映与大气的交换，土壤温湿廓线反映蓄热和水分状态。将这些通道放在共同时间轴上，可以区分天气背景、地表能量积累与土壤响应，而不是把三个高度气温重复解释为三个独立机制。')
image('q1/Q1_09_radiation_and_turbulent_energy.png','图4  净短波、净长波及湍流和土壤热通量的共同演变。',maxh=93*mm)
image('q1/Q1_10_multidepth_soil_profiles.png','图5  事件前24小时多深度土壤温度和体积含水率。',maxh=79*mm)

page('问题二：临界条件与补充通道','Q2 / 物源、水源与能量共同评估')
P('物源条件由冰厚、坡度、历史速度及速度梯度共同刻画，驱动剪应力采用τ=ρgH sinθ。动态评估包括24小时降水、源区高差修正后的7日和30日正积温、站点净短波与正能量余量当量、冻融转换。各通道沿用历史q90阈值，事件计算后不回调阈值。')
names={'source_PDD7_Cday':'源区7日正积温 / ℃·d','source_PDD30_Cday':'源区30日正积温 / ℃·d','station_net_shortwave_7d_mmwe':'净短波7日当量 / mm','station_positive_energy_residual_7d_mmwe':'正能量余量7日当量 / mm','source_freeze_thaw_transitions_7d':'冻融转换7日次数','rain24_mm':'24小时降水 / mm'}
table(['参数','发布前值','历史q90','历史百分位'],[[names[x['channel']],f"{x['value_at_issue']:.2f}",f"{x['historical_q90']:.2f}",f"{100*x['historical_percentile']:.1f}%"] for x in extra2['channels']],[.46,.18,.18,.18])
image('q2/Q2_18_energy_water_mechanism_channels.png','图6  新增能量通道与持续热量、降雨、冻融通道同时比较。',maxh=88*mm)
P('正能量余量通道超过历史q90，为本次无强降雨背景提供了可量化的热能积累线索。能量当量按正通量×3600/334000换算，单位为mm水当量；它描述站点可用能量，源区水热条件通过高差修正和空间物源条件共同解释。')

page('问题二：空间条件与尺度可靠性','Q2 / 空间物源、天气预报与检验')
image('q2/Q2_09_source_spatial_condition_index.png','图7  官方冰川网格上的坡度、驱动应力代理与静态物源条件。',maxh=80*mm)
P(f'历史84条沟谷以成因类别为目标进行空间留出比较，四地貌变量Logistic的OOF AUC为{float(models.loc[models.model.eq("L2_Logistic_4geo"),"auc_oof"].iloc[0]):.3f}。该分类检验与动态触发条件分别报告；30、60、120m空间尺度比较高分区交并比和排序相关，逐年向前拟合阈值并检验不同更新间隔的监测误差。')
image('q2/Q2_17_forward_time_reliability.png','图8  逐年向前阈值、条件超限及更新间隔误差。',maxh=89*mm)
P(f"补充遥感资料采用ICIMOD公开的Sentinel-1分析摘要，登记2026年7月7日至8月24日观测到的最高约{remote2['los_max_displacement_m']:.2f}m视线形变，作为事后物源演化线索；公开多传感器源区解译另用于第四问地理定位。",'small')
P(f'本次采用2026年8月25日18:00 UTC起报的ECMWF历史预报。预报格点相对申请源点偏移{forecast["model_grid_offset_km"]:.2f}km，并按格点与源区高程差调整温度；零度层高度缺失状态随原响应保存。跨境推算链的参数敏感性集中列于文末。','small')

page('问题二：能量衔接与源区代理','Q2 / 主指数与空间代表性检验')
P('主动态指数采用降雨、源区七日正积温和站点正能量余量三条门控的最大值。站点能量仍按站点尺度解释；它提供水热支持条件。未来七日能量只保留预测窗口内尚未过期的灾前实测正贡献，作为非负能量累积的下界。')
image('q2/Q2_19_source_proxy_and_weight_sensitivity.png','补图A  由官方冰川高程与观测序列重算九组高程-递减率组合，并比较四组声明的静态权重。',maxh=99*mm)
weights=read('q2/Q2_static_weight_sensitivity.csv')
table(['权重试验','相对参考空间秩相关','高分区IoU'],[[r.weight_experiment,f'{r.spearman_vs_reference:.3f}',f'{r.top20_IoU_vs_reference:.3f}'] for r in weights.itertuples()],[.4,.3,.3])
P('历史q90作为统一的条件参考阈值。相同模型在不同高程代理、递减率与权重下的响应一并保留，用于判断结论对空间代表性及参数设计的依赖。')

page('问题三：预测区间与预警时间','Q3 / 由实际预报轨迹给出前瞻预警')
P('将Q1预先固定的气温、短波和VPD异常通道与Q2水热条件组合，利用完整历史相似日轨迹保留变量间及时间上的相关结构；本次未来水热条件由实际历史起报预报驱动。对每条轨迹计算首次达到条件阈值的时间，直接以轨迹频数构造分布。')
P(f'条件首次超限的95%预测区间为{time_interval}。未来6小时超限频率为{q3["predicted_condition_cross_probability"]:.3f}，不超限频率为{q3["predicted_no_condition_cross_probability"]:.3f}。预设阈值{q3["forecast_warning_probability_threshold"]:.2f}，预测发布至题面10:30的提前量为{val(forecast_lead,1)}小时。')
P(f'持续监测的有效预警段起点另记为{q3["monitoring_episode_start_bjt"]}，回放至题面参考时刻的提前量{q3["first_monitoring_warning_advance_h_to_reference"]}小时。两类预警的发布、起点和有效期分别保存；决策仅使用发布时可得信息。','small')
image('q3/Q3_02_95percent_conditional_time_interval.png','图9  首次条件超限时间分布及其95%预测区间，条件于预测时域内发生超限。',maxh=77*mm)
image('q3/Q3_11_warning_time_and_lead_distribution.png','图10  历史新发超限的预警提前量及本次发布判断。',maxh=76*mm)

page('问题三：覆盖率与提前量检验','Q3 / 拟合、校准、检验分别进行')
P(f'采用2021至2023年拟合、2024年校准、2025年检验。独立检验期有{metrics["n_backtest"]}次可评价发布，其中{metrics["n_proxy_crossings"]}次出现气象条件超限；区间覆盖率为{pct(metrics["proxy_interval_coverage"])}，平均区间宽度为{metrics["mean_interval_width_h"]:.2f}小时。')
P(f'校准前覆盖率{pct(metrics["raw_proxy_interval_coverage"])}；校准后区间平均比完整6小时窗口缩短{pct(metrics["mean_width_reduction_vs_full_horizon"])}。覆盖率与区间宽度共同评价，当前天气驱动回放与历史相似日基准检验分别登记。','small')
table(['统计','次数或结果'],[['TP / FP / FN / TN',f'{metrics["TP"]} / {metrics["FP"]} / {metrics["FN"]} / {metrics["TN"]}'],['召回率 / 预警精确率',f'{pct(metrics["proxy_recall"])} / {pct(metrics["proxy_precision"])}'],['新发超限正确预警数',str(lead['historical_new_onset_lead']['n'])],['提前量中位数及四分位范围',f'{lead["historical_new_onset_lead"]["median_h"]:.2f} h；{lead["historical_new_onset_lead"]["q25_h"]:.2f}至{lead["historical_new_onset_lead"]["q75_h"]:.2f} h']],[.48,.52])
image('q3/Q3_08_full_forward_validation.png','图11  完整2025年检验期的预警概率、区间覆盖过程和三组输入比较。',maxh=124*mm)
P('覆盖率的评价对象为实际气象条件首次超限。监测状态机与预测预警分别记录：前者处理持续观测、解除和续期，后者依据未来轨迹概率发布，因此可以直接回答本次预警时间。','small')

page('问题三：空间距离与警戒廊道','Q3 / 将预测结果落实到地形空间')
P(f'空间范围以源区至口岸的河道连通关系为主线，依据Q4计算的最大水位与变断面几何确定沿程廊道，再融合原始DEM上的局部二维湿区。得到源区至口岸{space["source_to_port_corridor_km"]:.2f}km的警戒距离，模型空间包络{space["warning_envelope_km2"]:.2f}km²，其中局部二维累计湿润范围{space["nested_2d_wet_km2"]:.2f}km²。')
image('q3/Q3_12_reconstructed_spatial_warning.png','图12  释放网格、变断面廊道与口岸局部二维传播区共同组成空间警戒范围。',maxh=139*mm)
P('这一空间包络服务于本事件回放中的警戒分区。研究段终点设在口岸河段，所报距离表示应关注的连通廊道长度。坐标和每个组成部分均保存为EPSG:32645的GeoPackage，可与其他地理软件中的道路、建筑和行政区叠加。')

page('问题四：仿真环境与实际地形','Q4 / 从源区、河流汇入口到口岸',wide=True)
P(f'采用官方DEM构造覆盖主河道北向弯曲部分的计算范围，经Priority-Flood与D8汇流关系提取路径，沿原始高程进行动力学计算。路径长{geometry["length_m"]/1000:.2f}km，总落差{geometry["drop_m"]:.0f}m；口岸参考点与河槽终点相距{geometry["port_snap_distance_m"]:.0f}m。沿程断面由DEM横剖面提取并拟合梯形几何。')
image('q4/Q4_14_drainage_route_and_cross_sections.png','图13  河道连通路径、原始地形剖面与沿程断面；填洼用于寻路，水动力计算保留原始高程。',maxh=91*mm)
P('对梯形断面采用A=b h+m h²，静水压力矩为g(b h²/2+m h³/3)。计算中横断面可随地形变化，避免全程固定150m沟宽造成的局部过深。','formula')

page('问题四：源区投影与初始势能','Q4 / 让释放位置保持实际地形含义')
image('q4/Q4_25_release_geometry_and_elevation.png','补图D  源区实际像元与沿程投影的空间及高程对照。采用三维最近点，而非把源区长度平移至沟谷下段。',maxh=103*mm)
P(f'所选官方源区网格的平均高程为{geometry["selected_source_mean_elevation_m"]:.1f}m，投影后初始释放物质的平均床高为{geometry["assigned_release_mean_bed_elevation_m"]:.1f}m；逐像元平均绝对高程差为{geometry["release_projection_mean_abs_elevation_error_m"]:.1f}m。路径起点取源区上部真实像元，初始体积仍按实际像元面积分配。')
P('该检查约束初始重力势能，使传播时间和摩阻标定建立在实际释放高度之上。逐像元三维位置、投影路径距离及高程差均由代码导出。')

page('问题四：分阶段参数与体积传递','Q4 / 释放、滞留、再动员与下游传播')
P('将灾害链分为上游冰岩崩和下游泥石流两段。文献重建的初始方量、释放面积和可动方量作为事件约束；官方冰川像元的沿坡投影给出分布式释放范围，初始速度取零，由地形势能驱动运动。上游已经滞留的部分先计入沉积，再在汇入口分配剩余沉积与下游可动体积。')
table(['参数或结果','数值','取值依据'],[
 ['初始方量 / 释放面积',f'{observation["release_volume_m3"]/1e6:.1f}百万m³ / {observation["source_area_m2"]/1e6:.2f}km²','事件研究重建值'],
 ['上游数值滞留量',f'{up["stored_m3"]/1e6:.2f}百万m³','上游动力学结果'],
 ['汇入口再滞留量',f'{q4["deposit_partition_m3"]/1e6:.2f}百万m³','文献重建与上游收支共同确定'],
 ['下游可动输入',f'{q4["mobile_input_derived_m3"]/1e6:.2f}百万m³','实际上游累计通量按分配比例传递'],
 ['上游 μ / ξ',f'{q4["stage1_parameters"]["mu"]:.3f} / {q4["stage1_parameters"]["xi"]:.0f}','参数标定试验'],
 ['下游 μ / ξ',f'{par["mu"]:.5f} / {par["xi"]:.1f}','传播时间、流量、速度联合约束'],
 ['有效断面尺度',f'{par["section_capacity_scale"]:.3f}','基于DEM断面及口岸Q/u标定'],
 ['上游 / 下游 / 二维主网格',f'{up["dx_m"]:g} / {down["dx_m"]:g} / 30m','15m标定后固定参数，以更细网格评价']],[.34,.32,.34])
P('连续与动量方程采用有限体积方法。上游保留变断面Rusanov-Voellmy计算，下游采用MUSCL-HLL静水重构与两阶段SSP-RK2，以降低数值耗散。断面和床面共同进入重构，时间步满足CFL条件。阶段间保存累计体积与运动速度，下游直接读取相应结果。')

P('交汇区的堆积与回水条件','h2')
table(['相对填高','每100m河段的占据体积','占模型汇入口沉积比例'],[[f"{r['relative_fill_height_m']:.0f}m",f"{r['volume_per_100m_channel_m3']:,.0f}m³",pct(r['fraction_of_model_confluence_deposit'])] for r in deposition['results']],[.22,.44,.34])
P('上述体积由官方DEM横断面直接积分，100m仅为统一河长单位。与模型汇入口的沉积分配量相比，该结果表明交汇邻域具备明显的河槽占据潜势，为题面所述堆积阻塞和回水现象提供地形与物质条件上的解释。','small')

page('问题四：冰岩崩至泥石流的过程','Q4 / 动力学求解的流量、流速与流深')
image('q4/Q4_15_calibrated_cascade_process.png','图14  上游冰岩崩与下游泥石流的通量过程、口岸速度及沿程流深时空分布。',maxh=154*mm)
P(f'上游流动到达过渡断面的计算时刻为源区释放后{up["arrival_s"]:.1f}秒；下游通过真实累计通量接收{q4["mobile_input_derived_m3"]/1e6:.2f}百万m³，口岸河段首次湿润时刻为释放后{down["arrival_s"]:.1f}秒。口岸峰值流量为{down["peak_discharge_m3_s"]:.0f}m³/s，对应断面平均速度为{down["velocity_at_peak_discharge_m_s"]:.2f}m/s。全程体积相对残差为{q4["cascade_volume_residual"]:.2e}。')
P('各阶段均保存时空流深、最大流速、逐步流量及体积收支。源区至汇入口采用冰岩崩等效介质，下游采用可动泥石流等效介质，方量的转换由事件研究的沉积与再动员分配约束。')

page('问题四：事件约束与标定结果','Q4 / 报告计算值及其定量偏差')
P('上游以滞留方量和双侧时间误差约束，下游联合评价传播时间、峰值流量及其断面速度。以源区释放为时间零点，文献约130秒与350秒对应两段，总经历时间约480秒。数值前锋定义为局部流深首次超过0.2m；洪峰和累计量时刻另行导出。')
table(['对照量','文献重建','本次计算','相对差异'],[
 ['河道传播时间 / s',f'{observation["confluence_to_port_elapsed_s"]:.1f}',f'{q4["downstream_transit_s"]:.1f}',f'{100*(q4["downstream_transit_s"]/observation["confluence_to_port_elapsed_s"]-1):+.1f}%'],
 ['口岸峰值流量 / m³/s',f'{observation["port_peak_discharge_m3_s"]:.0f}',f'{down["peak_discharge_m3_s"]:.0f}',f'{100*(down["peak_discharge_m3_s"]/observation["port_peak_discharge_m3_s"]-1):+.1f}%'],
 ['断面平均速度 / m/s',f'{observation["port_cross_section_mean_velocity_m_s"]:.2f}',f'{down["velocity_at_peak_discharge_m_s"]:.2f}',f'{100*(down["velocity_at_peak_discharge_m_s"]/observation["port_cross_section_mean_velocity_m_s"]-1):+.1f}%']],[.37,.21,.21,.21])
image('q4/Q4_22_time_alignment_and_residuals.png','图15  统一时间起点，比较两段及全程时间，并展示五项标定量的全部相对偏差。',maxh=85*mm)
P('本次标定使可动体积与事件重建一致，并形成可以到达口岸的连续计算链。时间、流量和流速保留各自的计算偏差，便于评估有效断面、阻力与前缘判据的影响。参数和每次目标函数评价随计算结果保存；校验采用静水平衡、质量收支、能量范围与网格试验。')

page('问题四：口岸局部二维传播','Q4 / 原始DEM、实际入流与空间结果')
P('从下游一维模型提取口岸上游约2.5km断面的每步累计体积和速度，驱动局部二维浅层流方程。入口分布由横断面地形和一维计算水位确定；所有二维时刻沿用源区释放时钟。')
image('q4/Q4_18_reconstructed_2d_depth_speed_arrival.png','图16  同一地理范围上的最大流深、首次到达、最大流速及单位质量动能。动能按u²/2报告。',maxh=164*mm)
primary=two['grid_experiments'][0]
P(f'{primary["resolution_m"]:.0f}m网格计算的累计湿润面积为{primary["wet_area_km2"]:.3f}km²，二维最大等效流深为{primary["max_depth_m"]:.2f}m；二维体积相对残差为{primary["mass_relative_error"]:.2e}。空间数组、掩膜与坐标随绘图数据保存。','small')

page('问题四：六个传播阶段','Q4 / 连续演变的可视化',wide=True)
image('q4/Q4_19_reconstructed_2d_six_frames.png','图17  六个时刻采用共同地理范围和流深色标，红色方形标出口岸河段参考位置。每一帧来自实际时间步进数组。',maxh=120*mm)
P('过程图展示入流后沿谷传播、局部扩展及回落。与上游通量过程配合，可以检查何时进入局部区域、何时形成较大流深，以及网格与入口断面对传播范围的影响。','small')

page('问题四：全程传播与局部空间结果','Q4 / 源区至口岸的连续计算链',wide=True)
image('q4/Q4_23_whole_route_six_stages.png','补图B  把全程一维计算断面的流深投影到真实路线。各时刻由已保存的求解快照插值；图中着色表示沿程断面状态，口岸局部二维湿区由上一图单独给出。',maxh=119*mm)

page('问题四：修订前后的网格响应','Q4 / 固定输入与端点的分辨率比较')
image('q4/Q4_24_resolution_response_comparison.png','补图C  两版均分别固定各自标定参数，比较峰值流量与到达时间的网格响应。修订版保持同一物理起终点。',maxh=99*mm)
P('该图比较两版数值方案与各自标定结果的整体变化。当前版本内部的网格试验使用相同入流和阻力参数，并单列各网格的质量残差、流量与到达时间。数值分辨率响应与质量守恒分别核查。')

page('问题四：收支与数值分辨率','Q4 / 对仿真结果进行定量检查')
image('q4/Q4_20_cascade_volume_and_grid_experiments.png','图18  全程一维与局部二维分别核对体积收支，避免重复累加重叠计算域；同时比较二维网格。',maxh=74*mm)
table(['网格 / m','二维湿润面积 / km²','最大流深 / m','体积相对残差'],[[f'{r["resolution_m"]:.0f}',f'{r["wet_area_km2"]:.3f}',f'{r["max_depth_m"]:.2f}',f'{r["mass_relative_error"]:.2e}'] for r in two['grid_experiments']],[.16,.32,.26,.26])
image('q4/Q4_21_cross_sections_and_1d_grid_experiment.png','图19  从实际DEM提取横断面，并在固定地形和参数下改变一维数值步长。',maxh=75*mm)
P(f'一维采用{grid.dx_m.min():g}、{float(grid.dx_m.median()):g}和{grid.dx_m.max():g}m步长，分别记录到达时间、峰值流量和最大流深。二维以{two["grid_experiments"][0]["resolution_m"]:g}m与{two["grid_experiments"][1]["resolution_m"]:g}m对照。所有结果均为独立求解所得，细分一维网格增加数值分辨率，原始DEM信息保持一致。','small')

page('模型假设与不确定性','统一说明 / 解释数值的适用对象')
table(['环节','建模约定及检查方式'],[
 ['Q1关联因素','以事件窗口相对历史同期的异常评分定义关联性；多高度气温归入同一热力通道；完整原始排序保留以供复核。'],
 ['Q2跨境推算链',f'站点海拔{limits["station_elevation_m"]:.0f}m，源区代表高程{limits["source_elevation_m"]:.0f}m；主递减率{limits["lapse_rate_C_per_km"]:.1f}℃/km，另比较5和8℃/km。预报格点偏移{forecast["model_grid_offset_km"]:.2f}km，站点能量指标按站点尺度解释。'],
 ['Q3时间区间','题目所需95%时间区间按条件首次超限的预测区间构造；覆盖率以历史气象超限事件评价。本次10:30按题面灾害影响参考时刻计算预警提前量。'],
 ['Q4事后标定','方量与口岸断面约束来自事后事件研究，用于过程复现；Q3回放预报的天气输入与Q4标定资料分别记录。空间警戒范围为模型回放包络。'],
 ['地形与断面',f'官方DEM提供地形依据；寻路填洼与动力学原始高程分开。断面采用梯形近似，有效尺度{par["section_capacity_scale"]:.3f}，{boundary}；参数试验与原始断面共同展示其影响。'],
 ['传播与强度','等效介质模型描述体积和动量演变；沉积与可动体积采用分阶段约束。二维图显示流深、速度和单位质量动能，结构受力需进一步结合介质组成与构件参数。']],[.23,.77])
P(f"下游一维{primary_grid['dx_m']:g}m加密至{fine_grid['dx_m']:g}m后，到达时间变化{100*(fine_grid['arrival_s']/primary_grid['arrival_s']-1):+.1f}%，峰值流量变化{100*(fine_grid['peak_discharge_m3_s']/primary_grid['peak_discharge_m3_s']-1):+.1f}%，最大流深变化{100*(fine_grid['peak_depth_m']/primary_grid['peak_depth_m']-1):+.1f}%。二维30m与60m湿润面积差异为{100*abs(two['grid_experiments'][1]['wet_area_km2']/two['grid_experiments'][0]['wet_area_km2']-1):.1f}%。局部最大值与整体范围分别解释，质量收支与各量的分辨率响应分别评价。")
P('以上约定集中用于解释四问结果。正文按模型、计算结果和物理解释作答；来源分类、参数试验、完整图表数据及数值测试记录均随代码保留。')

for start in [0,6]:
    page(f'审查总结逐项复核（{1 if start==0 else 2}）','代码复核 / 对照赛题、源码与实际重算')
    table(['审查项','判断','依据与处理'],[[f"{r['id']}. {r['finding']}",r['judgement'],r['evidence']+' '+r['action']] for r in review['findings'][start:start+6]],[.19,.25,.56])

page('来源与复现','参考 / 数据、方法与代码对应')
refs=[
 '2026年MathorCup大数据竞赛赛道B初赛题面；国家青藏高原科学数据中心官方附件①至⑧。',
 'Li H, Li P, Hu K, et al. The 26 August 2026 Chhochen Khola ice-rock avalanche-debris flow cascade on the Nepal-China border. Journal of Mountain Science, 2026. DOI:10.1007/s11629-026-1482-y. 作为第四问分阶段方量、时间及断面约束来源。',
 'Li R. When a high-mountain slope failure cascades downstream: reconstructing the 26 August 2026 Gyirong mixed rock-ice disaster. arXiv:2609.04563. 作为事后源区位置的资料来源。',
 'ICIMOD. Second major Himalayan flood in seven weeks raises alarm over compounding mountain hazards. 2026-10-13. https://www.icimod.org/press-releases/second-major-himalayan-flood-in-seven-weeks-raises-new-questions-about-mountain-risk/. Sentinel-1研究摘要。',
 'Open-Meteo Single Runs API: https://open-meteo.com/en/docs/single-runs-api. 起报响应、变量缺测和来源哈希分别保存。',
 'RAMMS Debrisflow friction model: https://ramms.ch/ramms-debrisflow/friction-parameters/. 用于Voellmy形式与事件标定方法。',
 'Audusse E, Bouchut F, Bristeau MO, et al. A Fast and Stable Well-Balanced Scheme with Hydrostatic Reconstruction for Shallow Water Flows. SIAM J Sci Comput, 2004;25:2050-2065. DOI:10.1137/S1064827503431090.',
 'Gottlieb S, Shu CW. Total variation diminishing Runge-Kutta schemes. Mathematics of Computation, 1998;67:73-85. DOI:10.1090/S0025-5718-98-00913-2.',
 'Department of Hydrology and Meteorology, Nepal: Bhotekoshi at Rasuwagadi, station 4913, https://www.dhm.gov.np/hydrology/hms-Single/4913. 站点地理位置辅助核对。',
 'Wikimedia Commons, Nepal China Border.jpg: 原图作者的口岸附近相机地理标记，作为区域参考点；数值终点投影至DEM河槽。'
]
for i,r in enumerate(refs,1):P(f'[{i}] {E(r)}','small')
P('复现步骤','h2')
P('python -m pip install -r requirements.txt<br/>python run.py --stage all --no-fetch-q2-external<br/>python -m pytest -q tests<br/>python scripts/build_answer_pdf.py','formula')
P('原始附件与预报缓存按README中的目录恢复。模型参数文件记录标定输入指纹；数据或核心数值代码变化后重新标定。用--refit-q4可主动重跑标定实验。GitHub保留代码、统一README及PNG成果；本地计算同时导出数值表与逐图重绘数据。')

pdfpath=ROOT/'docs/MathorCup_B_v2_四问作答与原理分析.pdf';pdfpath.parent.mkdir(parents=True,exist_ok=True)
doc=BaseDocTemplate(str(pdfpath),pagesize=A4,leftMargin=18*mm,rightMargin=18*mm,topMargin=18*mm,bottomMargin=20*mm,title='MathorCup B 四问作答与计算结果 v2.2',author='数据分析与模型实验')
for name,size in [('P',A4),('W',landscape(A4))]:
    frame=Frame(18*mm,21*mm,size[0]-36*mm,size[1]-39*mm,id=name,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
    doc.addPageTemplates(PageTemplate(id=name,frames=[frame],onPage=foot,pagesize=size))
doc.build(STORY)
(ROOT/'outputs/quality').mkdir(parents=True,exist_ok=True)
(ROOT/'outputs/quality/pdf_build_metadata.json').write_text(json.dumps({'pdf':str(pdfpath),'planned_pages':len(PAGE_TITLES),'page_titles':PAGE_TITLES,'source':'current v2.2 model outputs and provenance'},ensure_ascii=False,indent=2))
print(pdfpath)
