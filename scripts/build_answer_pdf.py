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
def pct(x):return '不可估计' if x is None or pd.isna(x) else f'{x:.1%}'
def val(x,nd=2):return '未到达 / 不可估计' if x is None or pd.isna(x) else f'{x:,.{nd}f}'
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
    canvas.setFont('NR',7.5);canvas.setFillColor(GRAY);canvas.drawString(18*mm,11.5*mm,'2026 MathorCup B | v2.0 代码重算与四问作答');canvas.drawRightString(w-18*mm,11.5*mm,str(doc.page));canvas.restoreState()

q1=read('q1/Q1_top5_with_stability.csv');all1=read('q1/Q1_all_factors_effect_and_stability.csv');c1=j('q1/Q1_complete_factor_contract.json')
q2=j('q2/Q2_primary_trigger_assessment.json');lim=q2['thresholds'];now=q2['current'];models=read('q2/Q2_all_models_comparison.csv');fq=j('q2/Q2_archived_forecast_quality.json')
q3=j('q3/Q3_primary_competition_answer.json');m=q3['proxy_backtest'];space=j('q3/Q3_spatial_range_results.json')
q4=j('q4/Q4_1d_settings_and_summary.json');mass=q4['mass'];two=j('q4/Q4_2d_status.json');route=j('q4/Q4_route_quality_gate.json');grid=read('q4/Q4_grid_resolution_reliability.csv');sweep=read('q4/Q4_parameter_sensitivity.csv')
secondary=j('q4/Q4_secondary_coupling_case.json') if (OUT/'q4/Q4_secondary_coupling_case.json').exists() else None

page('四问作答与原理分析','数学建模 · 气象条件识别 / 时间预报 / 灾害链仿真')
STORY.append(Spacer(1,10));P('从代码修正到<br/>可复核的作答证据','cover')
P('基于题面、八类官方附件及当前源码重算。参考成品论文的论证结构，所有数值和图片均从本次运行结果读取。', 'body')
table(['问题','当前可交付结论','证据边界'],[
 ['Q1 特征关联',f'{c1["predefined"]}项候选、{c1["eligible"]}项有效；完整排序、原始前五项和机制代表项','一次事件与历史同期的关联分析'],
 ['Q2 触发评估','物源-能量-水源条件模型；实际历史天气预报；时空尺度检验','气候候选阈值尚非现场失稳强度'],
 ['Q3 时间预报',f'09:00-15:00名义95%情景区间；气象先兆检验覆盖率{pct(m["proxy_interval_coverage"])}','真实灾害时间覆盖率不可估计；发布时无有效黄色预警'],
 ['Q4 过程仿真',f'基准传播{mass["max_runout_m"]/1000:.2f} km；另列可到达参数情景的1D-2D传播','候选路径、等效断面及释放方量待现场核验']],[.16,.46,.38])
P('本次修订的重点是：使特征计算、样本划分、预警状态、地形源项和跨问传递共同服务于题目作答。程序通过数值测试，并不自动意味着对真实灾害具有相同精度。')
P('报告日期：2026年10月24日。时钟统一为北京时间（UTC+8）。', 'small')

page('数据如何进入四问','01 / 数据与证据')
table(['附件','实际使用内容','作用与限制'],[
 ['① MET/RADM/SOIL/FLUX','2021-2026逐小时；63个原始物理及诊断字段登记','Q1/Q2/Q3；不对长缺测插值'],
 ['② 区域日温','2005-2025逐日温度','长期暖季背景，不能代表冰川源区逐小时温度'],
 ['③ 冰川高程变化','目标区域2000-2020四期变化','多年物源背景，不提供2026失稳时间'],
 ['④ 冰川地形、厚度、速度','冰体格网、历史速度、厚度、坡度和速度梯度','Q2物源指标及候选源区；速度为2014-2024背景'],
 ['⑤ 土地覆盖','2020-2023四年有效覆盖区统计','空间背景；NoData不补成某类地表'],
 ['⑥ 历史灾点','84条沟谷的成因、地形和形态属性','成因分类；日期并非经核实的灾害触发时刻'],
 ['⑦ DEM','地形重投影、路径及局地二维计算','回顾性模拟；不能宣称已验证灾前地形'],
 ['⑧ 人口普查','县域人口背景','不按栅格均匀分配为精确受灾人口']],[.23,.36,.41])
P('时间边界','h2');P('当前案例统一在2026年8月26日09:00发布，输入观测截至08:00。题面10:30用于灾害影响时间参考；该时刻不等于已经独立确认的源区破裂时刻。历史模型拟合、区间校准和检验分别使用2021-2023、2024和2025年。')
P('公开资料中，ICIMOD的Sentinel-1形变综合报告发表于2026年10月13日，因此只作事后证据说明。ECMWF实际起报文件及其原始JSON、变量有效数和哈希均保留。')

page('问题一：关联如何定义','Q1 / 原理、变量与排序')
P('作答思路','h2');P('题面要求同时讨论边界层结构、天气过程、极端事件与前五项因素。为避免仅从少量最终变量反推原因，先登记全部候选，再统一比较事件前窗口与历史同期窗口，最后解释前五项及其冗余关系。')
P('对每个因素，计算灾前24小时的均值或规定累计量；在2021-2025年同一时钟、日历日期前后15天内建立参考窗口。每个窗口至少75%小时有效，参考至少35个窗口且来自3个年份。6、72和168小时用于尺度敏感性。')
P('z = (x事件 - median(x参考)) / [1.4826 × MAD(x参考)]<br/>S = 0.45 × min(|z|/3, 1) + 0.55 × (1 - p双侧经验尾部)','formula')
P('S同时反映偏离幅度和在历史参考中的稀有性。正负方向由z保留；S不是与多次泥石流发生标签计算的相关系数，也不是因果效应。零MAD、有效年数不足及离散诊断变量不强行评分。')
table(['原始名次','因素','得分 S','有向 MAD Z'],[[int(r.rank),r.label,f'{r.rarity_anomaly_score:.3f}',f'{r.z_mad:+.3f}'] for r in q1.itertuples()],[.13,.43,.2,.24])
P('这五项是按全部因素得分直接取前五。10 m、20 m和2 m气温存在明显物理冗余，因此另列机制代表项，用于描述水汽、热力、能量及土壤等通道，不能用代表项偷偷替换原始排序。')

page('全部因素：第1至28项','Q1 / 完整比较，而非只展示最终因素')
image('q1/Q1_11_part1_complete_ranking.png','图1. 全部有效因素的前半部分。橙色标明原始前五项。每个条形均对应全因素CSV中的一个评分。',maxh=190*mm)
P('前列以大气干燥度、不同高度气温和地表能量余量为主。多高度气温重复进入前五，反映相同热力背景，不能据此将其解释为三条独立致灾机制。')

page('全部因素：第29至56项','Q1 / 弱异常、负向异常与未入选项也需解释')
image('q1/Q1_11_part2_complete_ranking.png','图2. 全部有效因素的后半部分。低分表示本窗口中的异常或稀有性较弱，不代表在所有泥石流机制中不重要。',maxh=190*mm)
P(f'预定义{c1["predefined"]}项中，本次{c1["eligible"]}项达到评分条件；另有{c1["predefined"]-c1["eligible"]}项保留在排除清单。降水等变量在灾前可能为零、参考MAD也可能为零，此时必须报告原始降水量，而不应制造连续异常分数。')

page('从数值排序转向机制解释','Q1 / 边界层、能量交换与土壤结构')
image('q1/Q1_02_boundary_layer_profiles.png','图3. 观测层间的风温结构。20 m以内的廓线描述近地层，不能单凭该高度范围识别整层边界层顶。',maxh=88*mm)
image('q1/Q1_09_radiation_and_turbulent_energy.png','图4. 灾前辐射、感热、潜热与土壤热通量。能量收支余量描述站点闭合剩余，不等于冰川实际融化能。',maxh=100*mm)
P('VPD偏高与温度正异常共同反映站点空气较干、热力状态偏暖；辐射输入和湍流交换决定能量在地表与空气间的分配。站点与高海拔冰川并非同一地表，下垫面和海拔差异必须由第二问的源区条件分析承接。')

page('排序稳定性与第一问正式作答','Q1 / 多尺度证据与结论')
image('q1/Q1_05_reference_year_sensitivity.png','图5. 去冗余机制代表项的逐年剔除敏感性。完整因素另提供以参考年份为块的600次重采样得分、区间及前五频率。',maxh=78*mm)
image('q1/Q1_10_multidepth_soil_profiles.png','图6. 土壤多深度温湿结构，用于解释浅层与深层响应差异。',maxh=85*mm)
P('可用于论文的作答','h2');P('基于四类逐小时观测，采用历史同季同小时匹配，识别到灾前24小时最突出的因素依次为通量塔VPD、10 m气温、20 m气温、地表能量收支余量和2 m气温。它们表明站点存在偏暖、偏干及能量收支变化的共同背景。空间上独立的高海拔源区是否因此发生失稳，仍需物源与水源条件共同判断。')
P('窗口由24小时改为6、72或168小时后，名次会变化，说明这些因素对事件的关联具有尺度依赖。单次事件不能支持“高温是唯一触发原因”或“排除所有源区降水”的结论。')

page('问题二：临界条件如何组成模型','Q2 / 地形能量、物源与水源')
P('将冰厚、坡度、历史运动及速度梯度作为静态物源-能量条件；将海拔修正后的正积温与累计降水作为动态水源条件。应力代理衡量重力驱动力，不将其与未知抗剪强度直接相减形成伪安全系数。')
P('τ驱动 = ρ冰 g H sinθ<br/>T源区 = T站点 - Γ (z源区 - z站点)/1000<br/>PDD7 = Σ max(T源区, 0)/24<br/>I(x,t) = S物源地形(x) × max(G融水(t), G降雨(t))','formula')
P('静态权重为驱动剪应力0.40、厚度0.25、历史速度0.20、速度梯度0.15；每项以有效冰川像元的5%-95%分位归一化。速度梯度按真实x、y坐标求导，避免把随行向南递增误认为北向递增；缺少应变的像元按可用指标重新归一化，不能赋零风险。')
table(['关键量','计算结果','含义'],[
 ['站点 / 源区代表海拔',f'{lim["station_elevation_m"]:.0f} / {lim["source_elevation_m"]:.0f} m','源区使用冰川有效像元的代表高程'],
 ['气温垂直递减率','6.5 °C/km；另检验5和8','给定情景，非现场逐日实测'],
 ['源区7日正积温候选阈值',f'{lim["source_pdd7_warmseason_q90_Cday"]:.2f} °C·d','历史暖季90%分位'],
 ['24h累计降水候选阈值',f'{lim["rain_24h_nonzero_q90_mm"]:.2f} mm','历史非零降水窗口90%分位'],
 ['发布前条件',f'降水{now["rain_24h_mm"]:.2f} mm；源区正积温{now["source_pdd7_Cday"]:.2f} °C·d','本次未达到上述两项候选阈值']],[.32,.33,.35])
P('这些值回答“达到何种气候警戒条件”，尚不能回答岩冰体在何种强度、裂隙水压力或基底状态下必然破坏。缺测时动态指标保持未知，不能按无风险处理。')

page('全部模型及全部判别变量','Q2 / 成因判别作为空间辅助证据')
image('q2/Q2_15_all_model_comparison.png','图7. 四种候选模型及多数类基线的同一空间留出比较。分类对象是已知沟谷的成因，不是某时刻是否发生泥石流。',maxh=85*mm)
image('q2/Q2_16_all_predictor_distributions.png','图8. 全部五项候选地形、形态变量在两类历史沟谷中的分布。所有比较保留原始CSV及各空间折的标准化回归系数。',maxh=84*mm)
P(f'84条沟谷含43条冰相关型和41条降雨型。四地形项模型的样本外AUC为{models.loc[models.model.eq("L2_Logistic_4geo"),"auc_oof"].iloc[0]:.3f}；将沟比降替换为距断裂距离的四形态项模型为{models.loc[models.model.eq("L2_Logistic_4geomorph"),"auc_oof"].iloc[0]:.3f}。完整比较表明模型效果依赖变量组合，不能只引用最有利的单项结果。')
P('按纬度分为四个留出块可以降低随机拆分造成的空间混合，但相邻块仍可能有关联，且样本仅覆盖已知沟谷。该结果不能代替动态触发模型的独立灾害验证。')

page('物源空间证据与多年变化','Q2 / 冰体储量不等于释放量')
image('q2/Q2_13_glacier_velocity_and_release_candidates.png','图9. 历史冰川速度、速度梯度和连续候选源区。框定的是高分冰体分区，不是已经测得的2026破裂边界。',maxh=83*mm)
image('q2/Q2_07_glacier_historic_thinning_optional.png','图10. 官方冰川高程变化资料的空间背景。四个时期的数值均另存CSV，不仅使用最后一期。',maxh=93*mm)
P(f'所选高分连续分区的冰体库存约为{route["q2_source_candidate_ice_region"]["candidate_inventory_ice_volume_m3"]/1e6:.2f}百万m³。库存必须经过释放比例、岩石参与、水量及沿程裹挟等过程才能转换为下游混合物流量；不能直接把库存赋给仿真释放方量。')

page('时间与空间尺度可靠性','Q2 / 用误差和变化解释可靠性')
image('q2/Q2_17_forward_time_reliability.png','图11. 每个检验年份仅使用更早年份拟合阈值；改变监测更新间隔后，报告相对小时基准漏掉的超限小时和发现延迟。',maxh=84*mm)
image('q2/Q2_10_spatial_resolution_reliability.png','图12. 30、60、120 m物源指标的秩相关与高分区重叠。这是尺度一致性，不是灾害边界识别精度。',maxh=90*mm)
P('时间可靠性分成两层：候选阈值随参考年份变化是否稳定，以及降低更新频率会遗漏多少小时级异常。空间可靠性也分成两层：历史沟谷类别的空间留出表现，以及冰川指标粗化后高分区是否保留。四层证据分别解释，避免把抽稀或重采样的一致性写成真实触发准确率。')

page('实际历史天气预报及第二问作答','Q2 / 外部数据已进入计算')
image('q2/Q2_12_forecast_weather_conditions.png','图13. 实际2026-08-25 18:00 UTC起报的源区附近预报。零度层高度未返回有效值，图中明确留空。',maxh=131*mm)
table(['质量核查','结果'],[
 ['原始时段 / 发布后有效小时',f'48小时 / {fq["n_future_hours"]}小时'],
 ['温度、降水、短波有效记录',f'{fq["valid_count_by_variable"]["temperature_2m"]}、{fq["valid_count_by_variable"]["precipitation"]}、{fq["valid_count_by_variable"]["shortwave_radiation"]}'],
 ['零度层有效记录','0；不据此计算零度层阈值命中数'],
 ['预报格点相对申请源点偏移',f'{fq["model_grid_offset_km"]:.2f} km；需考虑山地代表性'],
 ['未来24小时格点降水',f'{fq["precip_future_24h_mm"]:.2f} mm；不能与站点零雨直接等同']],[.47,.53])
P('正式作答：已构建静态物源地形与动态水热条件联合评估模型，并补入真实起报预报。灾前站点24小时无雨、源区正积温代理未越过暖季90%候选阈值，说明该气象阈值体系未能独立解释本次失稳，仍需裂隙、冻融、基底强度及源区运动监测。该结论比由单站无雨直接断言高温融水致灾更符合证据。')

page('问题三：时间区间是什么量','Q3 / 气象预报与灾害时间情景分开')
P('模型组成','h2');P('使用固定的2 m气温、短波辐射和2 m气象VPD，以历史同月同小时中位数与MAD标准化，组合为气象异常项；再与第二问的融水/降雨条件按0.55和0.45组合。2026案例选出的前五项只作事后解释与敏感性，不进入历史检验的特征选择。')
P('取发布前21天中完整的相似时钟轨迹，并根据发布前最后状态进行递减修正，生成未来6小时、0.5小时步长的条件集合。每条轨迹的全部输入早于发布时刻。Q2实际历史天气预报另形成一条驱动敏感性，不将单个起报当成已有多年精度验证。')
P('λj = λ0 exp[β(Sj - c)]<br/>hj = 1 - exp(-λj Δt)<br/>P(T落在第j段) = hj × Πk&lt;j(1 - hk)<br/>P(T落在第j段 | T≤6h) = Pj / ΣPj','formula')
P(f'基准强度λ0由“参考暴露期假设1次事件”设定，β=2，尚未由独立灾害样本校准。本案例6小时内发生的模型情景概率为{q3["scenario_p_event_in_next_6h"]*100:.4f}%。它不应被当成已有实证支持的真实概率。')
table(['区间或指标','本次结果','能否作为真实灾害验证'],[
 ['名义95%触发时间情景区间','2026-08-26 09:00-15:00','否；条件是给定6小时内发生和强度假设'],
 ['真实灾害触发时间覆盖率','不可估计','缺少独立、可靠触发时刻目录'],
 ['气象超限时间经验覆盖率',f'{pct(m["proxy_interval_coverage"])}，n={m["n_proxy_intervals_evaluable"]}','只针对气象先兆，不能改名为灾害覆盖率'],
 ['题面10:30的区间位置','位于上述情景区间内','仅是参考时刻落入区间，不是独立覆盖率']],[.30,.32,.38])
P('题面使用“95%置信区间”一词。这里预测的是未来发生时间，而非总体参数，所以数学上采用条件预测区间。现有数据不足以声称该区间具备对真实灾害反复成立的95%覆盖保证。')

page('本次事件的条件演变与时间分布','Q3 / 当前案例的数值输出')
image('q3/Q3_01_precursor_and_scenario_hazard.png','图14. 灾前条件历史与发布后的相似日条件集合。任何发布后的观测都不参与预测。',maxh=108*mm)
image('q3/Q3_02_95percent_conditional_time_interval.png','图15. 危险率与生存概率生成的条件时间分布。区间覆盖整个6小时时域，说明当前信息未将发生时间约束到更窄范围。',maxh=87*mm)
P('不能因10:30落在09:00-15:00内便将模型评价为精准。需要同时报告区间宽度、未发生概率、误报情况和独立检验样本。为显示良好命中而事后缩窄区间，会破坏检验意义。')

page('完整历史检验：不只展示命中案例','Q3 / 2021-2023拟合、2024校准、2025检验')
image('q3/Q3_08_full_forward_validation.png','图16. 完整2025检验期：预测超限概率、实际超限时刻、校准区间，以及三组输入的累计覆盖。',maxh=135*mm)
table(['检验要素','结果'],[
 ['计划发布 / 可评价发布',f'{m["n_scheduled_issues"]} / {m["n_backtest"]}；另2次因输入或验证时段不完整而排除'],
 ['TP / FP / FN / TN',f'{m["TP"]} / {m["FP"]} / {m["FN"]} / {m["TN"]}'],
 ['气象先兆召回率 / 精确率',f'{pct(m["proxy_recall"])} / {pct(m["proxy_precision"])}'],
 ['气象超限FPR / 报警中的误报比例',f'{pct(m["proxy_false_alarm"])} / {pct(m["proxy_false_alarm_ratio"])}'],
 ['校准前 / 校准后覆盖率',f'{pct(m["raw_proxy_interval_coverage"])} / {pct(m["proxy_interval_coverage"])}'],
 ['校准后平均区间宽度',f'{m["mean_interval_width_h"]:.2f} h；检验目标为80次气象超限']],[.48,.52])
P('区间校准只使用2024年的63次气象超限，向两侧扩展1小时后再检验2025。数据存在时间相关，7日块重采样给出的覆盖率95%诊断范围约为91.4%-100%；该结果不能外推为灾害事件的共形覆盖保证。','small')

page('预警是否仍然有效','Q3 / 把首次报警与有效提前量区分开')
image('q3/Q3_09_warning_episodes.png','图17. 完整预警状态：连续2小时达标进入、达标时续期6小时、连续3小时低于阈值解除、缺测立即转未知。决策使用前一小时观测。',maxh=105*mm)
P('09:00发布时没有仍然有效的黄色预警，因此有效提前量不赋数值。历史上首次报警与10:30的时间差不能在报警已解除后继续作为有效提前量。09:00至10:30相隔1.5小时，只表示发布与题面影响时刻的距离，不能自动称作1.5小时成功预警。')
image('q3/Q3_10_archived_forecast_sensitivity.png','图18. Q2实际起报预报作为未来水热条件输入的敏感性，与相似日基准比较。它已进入运算，但单个起报不构成其独立预报精度证据。',maxh=85*mm)

page('空间范围与第三问正式作答','Q3 / 按物理可达性划定条件范围')
image('q3/Q3_05_spatial_warning_distance.png','图19. 基准实际传播段与完整候选通道分开。独立或次级参数情景的二维结果，不自动并入基准警戒面。',maxh=128*mm)
table(['空间量','结果与定义'],[
 ['候选路径长度',f'{space["route_length_km"]:.2f} km；不是经灾后影像核验的真实中心线'],
 ['基准情景传播距离',f'{space["primary_reached_distance_km"]:.2f} km'],
 ['参数敏感性传播范围',f'{space["sensitivity_runout_range_km"][0]:.2f}-{space["sensitivity_runout_range_km"][1]:.2f} km；不对应95%概率包络'],
 ['基准等效警戒面积',f'{space["combined_warning_envelope_area_km2"]:.2f} km²；假设150 m等效宽度'],
 ['基准是否接入下游二维范围','是' if space['2d_included_only_if_primary_1d_inflow'] else '否；基准未到耦合断面，故排除次级二维面']],[.39,.61])
P('正式作答：已构建条件驱动的短时预报与时间危险率模型，给出09:00-15:00的名义95%触发时间情景区间。真实灾害时间覆盖率仍不可估计；气象先兆的独立年份校准后覆盖率为96.3%。发布时无有效黄色预警，应报告该规则对本次事件存在漏警风险。空间范围由已计算传播段及明确参数情景限定，不能当成实测灾害边界。')

page('问题四：仿真环境与明确参数','Q4 / 数值环境、方程和情景')
P('先将官方DEM投影至UTM45N，利用有向地形代价路径连接Q2连续候选源区与暂定下游点，再沿通道建立一维有效混合物流动；有通量到达二维断面时，保留体积和入流动量，驱动局地二维模型。冰崩、裹挟和融化未被分别求解，不能将其表述为多相灾害链的完整机理重建。')
P('∂h/∂t + ∇·(hu) = q入流<br/>∂(hu)/∂t + ∇·(hu⊗u + gh²I/2)<br/>= -gh∇z - [μgh + g|u|²/ξ] u/|u| + q入流 u入流<br/>动压代理：p = ρ|u|²/2','formula')
table(['参数','基准值','比较范围或来源'],[
 ['有效混合物方量',f'{q4["parameters"]["volume"]/1e6:.0f}百万m³','30 / 50 / 80百万m³；不等于冰川库存'],
 ['等效沟宽','150 m','100 / 150 / 250 m；非实测逐断面地形'],
 ['释放历时 / 初速度','180 s / 65 m/s','180-600 s / 40-90 m/s'],
 ['库仑摩擦 μ / 湍流参数 ξ','0.040 / 2500 m/s²','0.025-0.060 / 1500-3500 m/s²'],
 ['密度 / 计算时长','1800 kg/m³ / 1800 s','密度仅用于动压诊断'],
 ['空间步长 / CFL','60 m / 0.35（一维）','二维另比较60与120 m']],[.32,.29,.39])
P('q入流为按网格面积折算的厚度源项，单位m/s。采用Rusanov有限体积通量；采用Audusse静水重构处理地形源项，使非平整底床上的静水保持平衡。该实现仍属一阶有效连续介质近似，需以网格和参数检验限定解释。')

page('地形、路径及一维演进','Q4 / 先看通道假设，再看仿真结果')
image('q4/Q4_01_route_geolocation_quality_gate.png','图20. 候选源区至下游的地理路径和高程纵剖面。路径长约22.24 km，累计局地爬升约323 m，提示存在路线及DEM粗化误差。',maxh=93*mm)
image('q4/Q4_02_1d_runout_timespace.png','图21. 基准一维流深的时空变化。未到达下游是模型在当前假设下的结果，应完整呈现。',maxh=107*mm)
P(f'基准传播距离为{mass["max_runout_m"]/1000:.2f} km，1800秒内未到达耦合断面及下游终点。最大等效流深{mass["max_depth_over_time_m"]:.1f} m，相对150 m假定沟宽过大，表明固定断面和局地壅滞会显著放大流深。该值必须作为模型适用性问题讨论，而非直接用于工程设计。')

page('守恒检查与数值正确性','Q4 / 通过守恒不等于真实事件复现')
image('q4/Q4_03_1d_hydrograph_and_mass_conservation.png','图22. 一维断面流量及累计体积收支。基准下游流量为空并非作图漏项，而是没有到达。',maxh=90*mm)
table(['检查','当前结果或含义'],[
 ['一维相对体积残差',f'{mass["relative_absolute_mass_residual"]:.2e}'],
 ['负流深截断体积',f'{mass["numerical_negative_clip_m3"]:.2e} m³'],
 ['地理方向','UTM北向速度映射为负栅格行向动量；新增北向传播测试'],
 ['静水与地形源项','非平整底床静水不产生虚假流动；一维、二维均有回归测试'],
 ['陡坡浅层加速','湿润坡面的加速度误差随网格细化降低；陡坡干湿前沿的粗网格误差另列'],
 ['一维-二维通量','每个接受时间步保留断面通量；二维按入流曲线分段积分'],
 ['动压解释','只对h>0.2 m的湿润像元计算，不能直接视作结构冲击荷载']],[.38,.62])
P('这些检查解决的是符号、源项、体积和过程传递的计算正确性。候选路线是否经过真实沟谷、释放物组成是否合理、二维边界是否足够大，必须另由地理和现场证据判断。')
P('在未到达的基准情景外，程序使用预设80百万m³上界方量情景演示有物理通量的二维传播。该情景被单独存档，不改变基准情景未到达的结论。')

page('可到达情景的二维结果','Q4 / 与基准分开报告的参数敏感性')
image('q4/Q4_05_local_2d_depth_and_arrival.png','图23. 预设上界方量情景的一维通量驱动二维最大流深及首达时间。所有帧沿用同一尺度和地理范围。',maxh=130*mm)
table(['二维情景量','结果'],[
 ['入流所属方量情景',f'{two.get("coupled_1d_volume_m3",0)/1e6:.0f}百万m³；与50百万m³基准分开'],
 ['最大湿润面积',f'{two["mass"]["wet_area_h_gt_0p2_km2"]:.3f} km²（h>0.2 m）'],
 ['最大等效流深',f'{two["mass"]["peak_depth_m"]:.2f} m'],
 ['最大动压代理',f'{two["mass"]["max_dynamic_pressure_kPa"]:.1f} kPa'],
 ['二维相对体积残差',f'{two["mass"]["relative_mass_error"]:.2e}']],[.44,.56])
P('颜色表示数值模型中的最大值；到达时刻相对源区一维释放计时。它展示满足较大方量条件后局地传播的可能过程，不是实际淹没范围的遥感复原。')

page('二维传播关键阶段','Q4 / 由保存数组重绘的过程图',wide=True)
image('q4/Q4_11_six_key_frames.png','图24. 从本次10帧计算数组中选取6个关键阶段。使用共同色标，标注源区释放后的时间。完整10帧图和动画保留在代码包中。',maxh=143*mm)

page('所有参数情景都必须看见','Q4 / 六项参数的完整比较')
image('q4/Q4_07_parameter_sensitivity.png','图25. 传播距离、全时段最大等效流深和到达时间一并比较。未到达者保留并标记，不能仅对已到达结果作趋势图。',maxh=86*mm)
image('q4/Q4_08_volume_friction_phase_diagram.png','图26. 方量与摩擦的3×3联合情景。未到达的格子显示“未到达”，不伪造到达时间。',maxh=93*mm)
P(f'六项参数各3个水平的全部18行记录均保留。传播距离覆盖{space["sensitivity_runout_range_km"][0]:.2f}-{space["sensitivity_runout_range_km"][1]:.2f} km，说明结果对方量、断面及阻力敏感；它不是一个经过校准的概率分布。联合时间图仅对可到达情景计算传播延迟，同时报告未到达数量。')

page('网格与入流边界的不确定性','Q4 / 第四问正式作答')
image('q4/Q4_09_grid_resolution_reliability.png','图27. 相同入流条件、同一地理方向下的60/120 m网格比较。两个网格只能说明敏感性，不能证明已经数值收敛。',maxh=89*mm)
rows=[]
for r in grid.itertuples():rows.append([f'{r.resolution_m:.0f} m',f'{r.wet_area_km2:.3f}',f'{r.max_depth_m:.2f}',f'{r.mass_balance_rel_error:.2e}'])
table(['网格','湿润面积 km²','最大流深 m','体积残差'],rows,[.2,.27,.27,.26])
P('边界位置与注入动量另外进行2×2试验，检查旧式3×3坡面注入、横断面定位以及沿沟方向动量是否导致人工分叉。位置调整、速度方向及源项计算均写入代码，不能只靠裁图掩盖分叉。')
P('正式作答','h2');P('已建立官方DEM约束的一维-二维有限体积仿真环境，明确方量、沟宽、释放历时、初速度、Voellmy阻力、密度、网格及边界条件，并输出时空流深、流量、首达、动压和过程帧。基准情景未到达下游；可到达的预设参数情景另外形成二维传播。全部情景、守恒及尺度检验共同表明：仿真能够给出条件演进过程，但当前不能声称已经准确重建真实灾害。')

page('一维网格与可达性诊断','Q4 / 固定地形，单独检验计算步长')
image('q4/Q4_13_1d_grid_resolution.png','图28. 同一候选高程路径上的30、60和120 m一维计算步长比较。只改变数值步长，不增加地形信息。',maxh=105*mm)
onegrid=read('q4/Q4_1d_grid_resolution_reliability.csv')
table(['1D步长','传播 km','下游首达 min','最大流深 m'],[[f'{r.resolution_m:.0f} m',f'{r.runout_km:.2f}',val(r.arrival_min),f'{r.max_depth_m:.1f}'] for r in onegrid.itertuples()],[.20,.25,.30,.25])
P('三种一维步长下，基准均未到达下游；传播距离略有变化，粗网格的流深峰值更低。因此不能仅凭某一网格的峰值来描述真实壅水。二维网格比较还叠加了DEM粗化和局地边界的影响，应与本页的一维数值误差分开解释。')
P('静水重构能避免不平整底床制造虚假静水流动，但当陡坡单元的床面高差超过局部流深时，干湿前沿仍存在粗网格局限。这也是保留30/60/120 m比较、最大等效流深和未到达情景的原因。模型的数值稳定、守恒和现场适用性是不同层次的结论。')

page('代码修订与题目验收对应关系','验收 / 从源数据到作答')
table(['题目要素','主代码','可复核产物'],[
 ['Q1全因素及前五项','features.py；q1_event_signatures.py；task_evidence.py','全62项纳入/排除；56项排序；原始前五项；机制代表项；年份重采样'],
 ['Q2临界条件与补充预报','q2_trigger_evaluation.py；q2_external.py','实际起报JSON；未来条件序列；变量缺测及格点偏移；阈值及物源栅格'],
 ['Q2时空可靠性','q2_spatial_mechanism.py；task_evidence.py','所有模型及变量；空间折系数；逐年向前阈值；更新频率误差'],
 ['Q3区间与回测','q3_task_aligned.py；warning_validation.py','拟合/校准/检验分离；完整发布序列；区间宽度及覆盖分母'],
 ['Q3有效预警与距离','warning_validation.py；q3_spatial_footprint.py','进入/解除/有效期；仅使用发布前输入；可达性约束GIS'],
 ['Q4数值与情景','q4_solver.py；q4_grid_review.py','方向、静水及加速测试；完整参数情景；同入流网格比较'],
 ['四问统一交付','publication_assets.py；report.py；本PDF脚本','全部图表用途索引；运行状态；源码及数据哈希；图文作答']],[.25,.34,.41])
P('仍不能由当前数据补出的信息','h2');P('源区真实破裂时刻及多场独立灾害目录；岩冰体抗剪强度与裂隙水压力；经核实的灾前源区/河道和逐断面几何；可靠的释放组成及沿程裹挟参数。程序保留对应结果为空或情景值，避免用气象阈值、影像解译日期或人工负样本替代真实观测。')
P('成品论文可借鉴的问题展开顺序、图表组织和论述形式已经吸收；其伪标签AUC、90%区间、人口圈层及未经核实的模拟数值没有直接移植。')

page('来源、复现与交付说明','附录 / 可核查来源')
P('数据与参考','h2')
refs=[
 '赛题：2026年MathorCup大数据挑战赛赛道B初赛；用户提供的八类官方数据附件及成品论文。成品论文仅作为组织参考。',
 '国家青藏高原科学数据中心：吉隆灾害专题及附件①-⑧。各附件用途、时段和单位在数据清单与质量记录中保留。',
 'Audusse E, Bouchut F, Bristeau MO, et al. A Fast and Stable Well-Balanced Scheme with Hydrostatic Reconstruction for Shallow Water Flows. SIAM J Sci Comput. 2004;25:2050-2065. DOI:10.1137/S1064827503431090。用于地形与静水平衡的数值思想；本代码另经静水、正流深和守恒行为测试。',
 'RAMMS Debrisflow：Voellmy阻力模型、方量/入流及参数含义。https://ramms.ch/ramms-debrisflow/friction-parameters/ 。模型形式依据，非本地参数实测。',
 'Li H, Li P, Hu K, et al. The 26 August 2026 Chhochen Khola ice-rock avalanche-debris flow cascade on the Nepal-China border. Journal of Mountain Science. 2026. DOI:10.1007/s11629-026-1482-y。仅借鉴103百万m³初始冰岩体与36百万m³可动碎屑的量级区分。',
 'ICIMOD. Second major Himalayan flood in seven weeks raises alarm over compounding mountain hazards. 2026-10-13。Sentinel-1最多0.40 m视线形变的事后综合证据，不能当作灾前已发布监测。',
 'Open-Meteo Single Runs API：https://open-meteo.com/en/docs/single-runs-api 。本次实际取得2026-08-25 18:00 UTC起报。以起报加6小时作为发布可用性假设，实际历史传播延迟未独立证实。'
]
for i,r in enumerate(refs,1):P(f'[{i}] {E(r)}','small')
P('复现命令','h2');P('python -m pip install -r requirements.txt<br/>python run.py --stage all --no-fetch-q2-external<br/>python -m pytest -q tests<br/>python scripts/build_answer_pdf.py','formula')
P('本PDF直接读取当前CSV、JSON和计算图像。代码包包含实际使用的解压数据、必要缓存、已下载的历史预报以及所有图表。运行全部流程会重新生成各题输出目录；无需手工改图。')
P('图表使用说明见outputs/FIGURE_INDEX.csv。每项产物均记录所属问题、用途、生成代码及主要数值表；可据此将作答拆入正式论文的各节。')

pdfpath=ROOT/'docs/MathorCup_B_v2_四问作答与原理分析.pdf'
doc=BaseDocTemplate(str(pdfpath),pagesize=A4,leftMargin=18*mm,rightMargin=18*mm,topMargin=18*mm,bottomMargin=20*mm,title='MathorCup B 四问作答与原理分析 v2.0',author='代码重算与证据分析')
for name,size in [('P',A4),('W',landscape(A4))]:
    frame=Frame(18*mm,21*mm,size[0]-36*mm,size[1]-39*mm,id=name,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
    doc.addPageTemplates(PageTemplate(id=name,frames=[frame],onPage=foot,pagesize=size))
doc.build(STORY)
write={'pdf':str(pdfpath),'planned_pages':len(PAGE_TITLES),'page_titles':PAGE_TITLES,'source':'current outputs only'}
(ROOT/'outputs/quality/pdf_build_metadata.json').write_text(json.dumps(write,ensure_ascii=False,indent=2))
print(pdfpath)
