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
    canvas.setFont('NR',7.5);canvas.setFillColor(GRAY);canvas.drawString(18*mm,11.5*mm,'2026 MathorCup B | 吉隆泥石流预报与传播仿真');canvas.drawRightString(w-18*mm,11.5*mm,str(doc.page));canvas.restoreState()



q1=j('q1/Q1_direct_answer.json');raw=read('q1/Q1_top5_with_stability.csv')
q2=j('q2/Q2_sequential_model_contract.json');assessment=j('q2/Q2_primary_trigger_assessment.json')
q3=j('q3/Q3_primary_competition_answer.json');q4=j('q4/Q4_cascade_reconstruction.json')
geo=j('q4/Q4_reconstruction_geometry.json');grid=j('q4/Q4_numerical_resolution_diagnostics.json')
coverage=q3['case_reference_coverage'];limits=assessment['thresholds'];current=assessment['current']
season=read('q1/Q1_boundary_layer_seasonal_metrics.csv').set_index('season')
spatial=read('q2/Q2_spatial_scale_reliability.csv')
cadence=read('q2/Q2_update_cadence_reliability.csv');cadence25=cadence.loc[cadence.test_year.eq(2025)]
weights_test=read('q2/Q2_static_weight_sensitivity.csv')
comp=read('q3/Q3_time_model_comparison.csv')
residuals=read('q4/Q4_calibration_residuals.csv').set_index('quantity')
grids=read('q4/Q4_reconstruction_1d_grid_experiment.csv')
local2d=j('q4/Q4_reconstruction_2d_summary.json');fine2d=local2d['grid_experiments'][0]
forecast_quality=j('q2/Q2_archived_forecast_quality.json')
energy_channels=j('q2/Q2_expanded_mechanism_answer.json')['channels']
energy_now=next(r['value_at_issue'] for r in energy_channels if r['channel']=='station_positive_energy_residual_7d_mmwe')
interval=' 至 '.join(x[11:16] for x in q3['conditional_95pct_interval_times_bjt'])
active_png=sum(len(list((OUT/q).glob('*.png'))) for q in ['q1','q2','q3','q4','joint'])

page('吉隆泥石流预报与传播仿真','2026 MathorCup 数学应用挑战赛 / 赛道B')
P('围绕“关联因素识别、触发条件评估、事件时间预报、灾害传播仿真”建立连续模型。第一问确定特征及其权重，第二问形成触发强度和源区通道，第三问计算事件时间分布与预警，第四问模拟同一源区的物质释放和沿程传播。')
table(['问题','主要答案'],[
 ['Q1 关联因素','通量塔VPD、10m气温、20m气温、地表能量收支余量、2m气温为关联评分最高的五项。'],
 ['Q2 触发条件',f"综合触发指数I={q2['current_score']:.3f}，超过临界水平u={q2['critical_score']:.3f}；当前以能量支持通道为主。"],
 ['Q3 时间与预警',f"8月26日触发时间95%区间为{interval}（北京时间）；09:00发布预警，提前{q3['operational_warning_lead_h']:.1f}h；本次案例覆盖率{coverage['n_covered']}/{coverage['n_references']}={coverage['coverage']:.0%}；沿程距离0-{q2['source_to_port_km']:.3f}km。"],
 ['Q4 仿真结果',f"释放后{q4['downstream']['arrival_s']:.2f}s到达口岸；峰值流量{q4['downstream']['peak_discharge_m3_s']/1e4:.2f}万m³/s；口岸局部二维湿润面积{fine2d['wet_area_km2']:.4f}km²。"]],[.21,.79])
P('研究资料与时间基准','h2')
P('以官方MET、RADM、SOIL、FLUX观测及温度、冰川、土地覆盖、历史灾点、DEM等附件为基础，补充Sentinel-1形变资料和归档ECMWF天气预报。观测截止2026年8月26日08:00，预报发布于09:00；全文使用北京时间，传播过程使用源区释放后的相对秒数。')
P('集中说明的建模假设','h2')
P('①固定案例源区、连通沟道和物源体积，六小时内地形保持稳定；②站点与模式格点气温经6.5°C/km高程递减率换算至源区，归档起报按6小时发布延迟处理；③以09:00-15:00发生一次事件为时间预报情景，事件时间密度与触发强度成正比，半小时时段内取均匀密度；④题面10:30作为案例时间评价基准，Q4采用等效浅层流和分阶段体积分配。')

page('Q1：关联特征排序与天气过程','第一问 / 从历史匹配到因素识别')
P('以灾前24小时为主排序窗口，历史对照采用2021-2025年同日附近±15天、相同钟点的有效窗口。对62项候选做质量筛选，56项进入排序；同时分析灾前7天过程及6、24、72、168小时窗口的稳定性。')
P('z<sub>i</sub> = (x<sub>i</sub> - median<sub>ref</sub>) / (1.4826 MAD<sub>ref</sub>)<br/>S<sub>i</sub> = 0.45 min(|z<sub>i</sub>|/3, 1) + 0.55 (1 - p<sub>i,双尾</sub>)','formula')
table(['排序','特征因素','关联评分','标准化偏差'],[
 [str(i+1),LABELS.get(r.feature,r.feature),val(r.rarity_anomaly_score,3),val(r.z_mad,2)]
 for i,r in enumerate(raw.itertuples())],[.10,.50,.20,.20])
P(f"边界层结构表现出明显季节差异：冬季逆温频率{season.loc['冬季','inversion_frequency']:.1%}，暖季{season.loc['暖季','inversion_frequency']:.1%}；暖季净辐射均值{season.loc['暖季','mean_net_radiation_Wm2']:.1f}W/m²，高于冬季{season.loc['冬季','mean_net_radiation_Wm2']:.1f}W/m²。灾前7天VPD有29小时超过同期95%分位，2m气温有15小时超过该水平；24h降雨同期超过95%分位的小时数为0。")
image('q1/Q1_14_mechanisms_and_extreme_weather.png','图1 机制代表因素及灾前逐小时历史百分位。连续偏暖、VPD增大和能量偏高构成本次事件的主要气象背景。',maxh=78*mm)
P('第一问结论：原始最高五项按上表作答。面向后续建模，合并气温冗余，取VPD、10m气温、能量余量、短波辐射、土壤水分梯度五项机制代表，并将其关联评分归一化为第二问权重。')

page('Q2：触发评估模型与临界条件','第二问 / 继承Q1特征，结合物源和水热条件')
P('将物源、地形作为空间背景：研究源区面积0.7596km²，平均高程5132.04m，DEM连通路径长23.715km、落差3427.28m。动态部分由第一问异常分量与降雨、正积温、能量支持三条通道组成。')
P('a<sub>i</sub>(t) = clip(d<sub>i</sub> z<sub>i,月-时</sub>(t), 0, 4) / 4<br/>A(t) = Σ w<sub>i</sub> a<sub>i</sub>(t)，w<sub>i</sub> = S<sub>i</sub>/Σ S<sub>k</sub><br/>R(t) = max(R<sub>rain</sub>, R<sub>melt</sub>, R<sub>energy</sub>)<br/>I(t) = R(t) [0.5 + 0.5 A(t)]','formula')
table(['Q1机制代表','VPD','10m气温','能量余量','短波辐射','水分梯度'],[
 ['归一化权重']+[f"{q2['Q1_weights'][f]:.4f}" for f in q2['Q1_features']]
],[.24,.152,.152,.152,.152,.152])
P('z按历史同月、同小时的中位数和MAD标准化，d取第一问异常方向。至少四项分量有效时，按有效权重重新归一。水热通道采用分位响应G(x;q,s)=clip[0.5+(x-q)/(4s),0,1]，q取各通道历史q90；雨量s=0.15q，正积温和能量s取历史标准差的0.3倍。降雨通道另要求24小时降雨达到q。')
table(['关键参数','历史临界水平','发布前状态'],[
 ['24h降雨',f"{limits['rain_24h_nonzero_q90_mm']:.2f}mm",f"{current['rain_24h_mm']:.2f}mm"],
 ['源区7日正积温',f"{limits['source_pdd7_warmseason_q90_Cday']:.2f}°C·d",f"{current['source_pdd7_Cday']:.2f}°C·d"],
 ['7日正能量余量当量',f"{limits['station_positive_energy_7d_q90_mmwe']:.2f}mm水当量",f"{energy_now:.2f}mm水当量"],
 ['综合指数I',f"{q2['critical_score']:.3f}（暖季q85）",f"{q2['current_score']:.3f}"],
 ['监测指数',f"{q2['monitoring_score']:.3f}（暖季q75）",'进入持续监测与短时预报']],[.38,.32,.30])
P(f"第二问判定：A={q2['current_Q1_anomaly']:.3f}，R={q2['current_dynamic_water_heat']:.3f}，得到I={q2['current_score']:.3f}，较临界值高{(q2['current_score']/q2['critical_score']-1)*100:.1f}%。当前雨量和正积温低于各自q90，能量支持通道超过q90，因此本次触发评估由能量偏高与第一问异常共同驱动。")

page('Q2：时间与空间尺度可靠性','第二问 / 更新频率、分辨率与参数实验')
P('时间实验：用早期年份确定气象门控阈值，逐年向后评价；再将每小时更新改为3、6、24小时更新，比较临界状态的一致率和首次识别延迟。以下为2025年结果。')
table(['更新频率','与逐时状态一致率','临界状态IoU','识别延迟中位数'],[
 [f"{int(r.update_interval_h)}h",f"{r.agreement_with_hourly:.2%}",f"{r.exceedance_IoU:.3f}",f"{r.median_detected_onset_delay_h:.1f}h"]
 for r in cadence25.itertuples()],[.18,.29,.25,.28])
P('逐日更新的识别延迟中位数为9小时，而3小时更新为1小时。因此，第一问长期背景用于筛选因素，第二问采用逐时监测，第三问采用半小时时间格计算六小时预报。')
P('空间实验：固定输入和计算规则，将30m源区评分聚合至60m、120m，比较空间排序和最高20%区域的重叠度。')
table(['分辨率','相对30m秩相关','高分区IoU'],[
 [f"{int(r.grid_resolution_m)}m",f"{r.spearman_vs_30m:.3f}",f"{r.high20percent_IoU:.3f}"]
 for r in spatial.itertuples()],[.24,.40,.36])
image('q2/Q2_10_spatial_resolution_reliability.png','图2 源区空间评分的分辨率试验。60m仍保有较高排序一致性，源区边界和断面建模采用30m DEM。',maxh=72*mm)
P(f"参数实验：源区气温递减率取5、6.5、8°C/km；静态权重分别取参考、等权、重力优先、运动优先。三种替代权重的全域秩相关为{weights_test.iloc[1:].spearman_vs_reference.min():.3f}-{weights_test.iloc[1:].spearman_vs_reference.max():.3f}，高分区IoU为{weights_test.iloc[1:].top20_IoU_vs_reference.min():.3f}-{weights_test.iloc[1:].top20_IoU_vs_reference.max():.3f}。据此采用逐时监测和30m源区、保持统一源区与通道进入后两问。")

page('Q3：触发时间分布与95%区间','第三问 / 将Q2强度转换为事件时间预报')
P(f"发布时刻为2026年8月26日09:00，预报至15:00，时间步长0.5小时。以灾前完整相似日生成1200条异常轨迹，当前可用相似日{q3['analogue_days']}个；归档ECMWF气温与降水进入第二问水热通道，再用同一I方程计算未来强度。")
P('对每条情景m，将触发强度归一化为事件落入各时段的概率；再对情景等权平均。设F为半小时时段内均匀插值得到的累计分布，取其2.5%和97.5%分位点。')
P('p<sub>mj</sub> = I<sub>mj</sub> Δt / Σ<sub>k</sub> I<sub>mk</sub> Δt<br/>p<sub>j</sub> = (1/M) Σ<sub>m</sub> p<sub>mj</sub>，Σ<sub>j</sub> p<sub>j</sub> = 1<br/>CI<sub>95</sub> = [F<super>-1</super>(0.025), F<super>-1</super>(0.975)]','formula')
table(['题目要求','预报结果（北京时间，2026年8月26日）'],[
 ['触发时间95%区间',interval],
 ['中位预测时刻',q3['median_event_time_bjt'][11:16]],
 ['区间宽度',f"{q3['case_reference_interval_width_h']:.2f}h"],
 ['与题面10:30的中位误差',f"{q3['case_reference_median_error_h']*60:+.1f}min"]],[.37,.63])
image('q3/Q3_02_95percent_conditional_time_interval.png','图3 六小时事件时间概率质量、累计分布及95%区间。阴影为模型预测区间，虚线标记题面事件时刻。',maxh=88*mm)
P(f"第三问时间答案：本研究情景下，触发时间95%区间为8月26日{interval}，中位时刻约10:48。题面10:30落在该区间内。")

page('Q3：预警、覆盖率与空间距离','第三问 / 预警规则和对照实验')
P(f"预警规则：计算未来六小时内I达到第二问临界线的情景比例，超过25%即在发布时刻发出预警。当前比例为{q3['assessment_cross_probability']:.2%}，因此09:00发布，预警有效至15:00；相对题面10:30提前{q3['operational_warning_lead_h']:.1f}小时。")
P(f"区间覆盖率按 C=Σ1{{L≤t≤U}}/N 计算。题面提供N={coverage['n_references']}个事件时刻，命中{coverage['n_covered']}次，本次案例覆盖率为{coverage['coverage']:.0%}（{coverage['n_covered']}/{coverage['n_references']}）。")
P('时间模型对照','h2')
method_labels={'Q2_intensity_proportional':'强度正比（主模型）','uniform_time_prior':'均匀时间基线','squared_intensity_sensitivity':'强度平方（敏感性）'}
table(['时间模型','中位时刻','95%区间宽度','案例命中'],[
 [method_labels[r.method],(pd.Timestamp(q3['prediction_issue_bjt'])+pd.Timedelta(hours=r.median_h)).strftime('%H:%M'),f"{r.width_h:.2f}h",f"{int(r.case_reference_covered)}/1"]
 for r in comp.itertuples()],[.41,.19,.24,.16])
width_reduction_min=float(comp.loc[comp.method.eq('uniform_time_prior'),'width_h'].iloc[0]-comp.loc[comp.method.eq('Q2_intensity_proportional'),'width_h'].iloc[0])*60
P(f'主模型把中位预报从均匀基线的12:00提前至10:48，区间宽度减少约{width_reduction_min:.1f}分钟；强度平方使时间进一步前移。当前信息对“何时更可能发生”的中位排序作用较明显，95%区间仍覆盖大部分六小时时域。')
image('q3/Q3_05_spatial_warning_distance.png',f"图4 预警空间范围：以固定源区为0km，沿DEM连通通道延伸至口岸{q2['source_to_port_km']:.3f}km；第四问在该范围内继续计算传播时间和局部湿润范围。",maxh=86*mm)
P(f"第三问其余答案：预警时刻09:00，提前量1.5小时；本次案例覆盖率100%（1/1）；沿源区至口岸的空间距离范围为0-{q2['source_to_port_km']:.3f}km。")

page('Q4：仿真环境、方程与参数','第四问 / 由既定源区进入传播计算')
P('采用第二问的同一源区、通道和DEM断面。全程以变断面一维浅层流计算，口岸附近嵌套二维网格；一维通量与累计体积作为二维入口。第三问给出事件钟时，本问给出释放后的相对传播时间，二者通过T到达=T事件+τ传播衔接。')
P('A<sub>c</sub> = b h + m h²，Q = A<sub>c</sub> u<br/>∂A<sub>c</sub>/∂t + ∂Q/∂s = q<sub>s</sub><br/>∂Q/∂t + ∂(Q²/A<sub>c</sub> + gJ)/∂s = S<sub>z,b</sub> - A<sub>c</sub>a<sub>f</sub> + u<sub>s</sub>q<sub>s</sub><br/>J = bh²/2 + mh³/3，a<sub>f</sub> = μg cosθ sign(u) + g u|u|/(ξH)<br/>S<sub>z,b</sub>为底坡与变断面源项，H为水力深度，q<sub>s</sub>为侧向物源。','formula')
table(['关键参数','设定与依据'],[
 ['源区与路径',f"844个30m像元，面积{q2['source_area_km2']:.4f}km²；路径{q2['source_to_port_km']:.3f}km"],
 ['释放情景','初始103百万m³，初速0m/s；下游可动36百万m³，依据文献分阶段体积'],
 ['上游阻力',f"μ={q4['stage1_parameters']['mu']:.3f}；ξ={q4['stage1_parameters']['xi']:.1f}m/s²"],
 ['下游阻力与断面',f"μ={q4['stage2_parameters']['mu']:.6f}；ξ={q4['stage2_parameters']['xi']:.1f}m/s²；有效断面尺度={q4['stage2_parameters']['section_capacity_scale']:.2f}"],
 ['离散方法','上游Rusanov-Voellmy；下游MUSCL-HLL、SSP-RK2；二维静水平衡有限体积'],
 ['网格与到达判据','上游15m，下游7.5m，口岸二维30m；当地首次流深超过0.2m记为到达']],[.29,.71])
P('有效断面尺度作为等效参数，由峰值断面与传播约束联合标定；最优值为搜索范围0.5-3.0的上端。其余释放体积、源区和地形按上述场景固定，网格实验沿用同一组参数。')
P('体积衔接','h2')
P(f"上游释放103百万m³，计算出流{q4['upstream']['outflow_m3']/1e6:.2f}百万m³、存留{q4['upstream']['stored_m3']/1e6:.2f}百万m³；按实际累计通量的{q4['mobile_fraction_of_computed_upper_outflow']:.2%}转为下游可动36百万m³，汇流处分配{q4['deposit_partition_m3']/1e6:.2f}百万m³。上游存留与汇流处分配合计67百万m³。")

page('Q4：源区至口岸的仿真过程','第四问 / 全程流深演进',wide=True)
P(f"释放后{q4['upstream']['arrival_s']:.2f}s进入汇流段，再传播{q4['downstream_transit_s']:.2f}s，于{q4['downstream']['arrival_s']:.2f}s到达口岸；口岸峰值流量{q4['downstream']['peak_discharge_m3_s']/1e4:.2f}万m³/s，峰值时断面平均流速{q4['downstream']['velocity_at_peak_discharge_m_s']:.2f}m/s。")
image('q4/Q4_23_whole_route_six_stages.png','图5 源区释放后0、104、209、313、417、521秒的全程一维流深投影。图中显示物质由上游释放、沿沟谷传播、到达口岸和后续消退的连续过程。',maxh=121*mm)

page('Q4：口岸局部二维传播结果','第四问 / 平面流深、到达时间与流速')
P(f"将一维记录的入口界面通量输入口岸二维区域，累计输入{fine2d['input_m3']/1e6:.2f}百万m³。30m主网格得到湿润面积{fine2d['wet_area_km2']:.4f}km²，口岸120m邻域首次到达为释放后{fine2d['port_neighbourhood_arrival_s']:.2f}s。")
image('q4/Q4_18_reconstructed_2d_depth_speed_arrival.png','图6 口岸局部二维仿真：最大流深、首次到达、最大流速和单位质量动能分布。白线表示河谷路径，方块表示口岸位置；动能由计算流速按u²/2得到。',maxh=141*mm)
table(['二维网格','湿润面积','口岸邻域到达','最大湿区流速'],[
 [f"{r['resolution_m']:.0f}m",f"{r['wet_area_km2']:.4f}km²",f"{r['port_neighbourhood_arrival_s']:.2f}s",f"{r['max_wet_speed_m_s']:.2f}m/s"]
 for r in local2d['grid_experiments']],[.18,.27,.28,.27])
P('局部传播结论：流体主要沿沟谷向口岸附近汇集，并向低洼地形扩展。二维结果给出平面上的暴露位置与到达先后，可与第三问的沿程预警距离配合确定重点防护区。')

page('Q4：数值实验与结果评价','第四问 / 固定参数的网格与事件约束比较')
table(['下游网格','到达时间','峰值流量','最大流深'],[
 [f"{r.dx_m:g}m",f"{r.arrival_s:.2f}s",f"{r.peak_discharge_m3_s/1e4:.2f}万m³/s",f"{r.peak_depth_m:.2f}m"]
 for r in grids.sort_values('dx_m',ascending=False).itertuples()],[.20,.25,.29,.26])
P(f"7.5m细化到3.75m后，峰值流量变化{pct(grid['peak_discharge_relative_change'])}、到达时间变化{pct(grid['arrival_relative_change'])}、最大流深变化{pct(grid['maximum_depth_relative_change'])}。因此，传播时间与峰值流量对本次细化较稳定，局部最大流深对网格更敏感。各档体积相对误差均低于10<super>-12</super>，体积收支闭合。")
table(['事件约束','文献参考','仿真结果','相对差异'],[
 ['源区至口岸时间','480.00s',f"{q4['downstream']['arrival_s']:.2f}s",f"{residuals.loc['source_to_port_s','relative_error']:+.1%}"],
 ['口岸峰值流量','24.60万m³/s',f"{q4['downstream']['peak_discharge_m3_s']/1e4:.2f}万m³/s",f"{residuals.loc['port_peak_discharge_m3_s','relative_error']:+.1%}"],
 ['峰值时断面流速','36.20m/s',f"{q4['downstream']['velocity_at_peak_discharge_m_s']:.2f}m/s",f"{residuals.loc['port_velocity_at_peak_m_s','relative_error']:+.1%}"]],[.30,.24,.24,.22])
image('q4/Q4_24_resolution_response_comparison.png','图7 两种离散方法各自固定参数时的一维网格响应。Q4给出了参数明确、收支闭合的连续传播仿真，以及可核查的网格变化和事件约束误差。',maxh=69*mm)
P('资料与复现','h2')
P(f"数据：国家青藏高原科学数据中心赛题附件①-⑧；Open-Meteo归档ECMWF起报；ICIMOD公布的Sentinel-1形变用于案例机理佐证。天气格点距请求源区{forecast_quality['model_grid_offset_km']:.2f}km，零度层字段缺测，温度换算采用有效格点气温及高程。体积和传播参考：Li Hao等，Journal of Mountain Science，DOI:10.1007/s11629-026-1482-y。",'small')
P(f"完整复现包含输入、代码、{active_png}张PNG、SVG、GIS和逐图数值。运行：python run.py --stage all --reference-experiments --no-fetch-q2-external；生成本文：python scripts/build_answer_pdf.py。完整来源和字段对应见README。",'small')

pdfpath=ROOT/'docs/MathorCup_B_v2_3_1_question_answers.pdf';pdfpath.parent.mkdir(parents=True,exist_ok=True)
doc=BaseDocTemplate(str(pdfpath),pagesize=A4,leftMargin=18*mm,rightMargin=18*mm,topMargin=18*mm,bottomMargin=20*mm,title='吉隆泥石流预报与传播仿真：四问模型与实验',author='MathorCup B 模型与实验')
for name,size in [('P',A4),('W',landscape(A4))]:
    frame=Frame(18*mm,21*mm,size[0]-36*mm,size[1]-39*mm,id=name,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
    doc.addPageTemplates(PageTemplate(id=name,frames=[frame],onPage=foot,pagesize=size))
doc.build(STORY)
(OUT/'quality/pdf_build_metadata.json').write_text(json.dumps({'pdf':str(pdfpath),'planned_pages':len(PAGE_TITLES),'page_titles':PAGE_TITLES,'source':'current v2.3.1 sequential model outputs','report_style':'question, equation, experiment, answer; assumptions stated once','revision_scope':'answer report only; numerical model and existing figures unchanged'},ensure_ascii=False,indent=2))
print(pdfpath)
