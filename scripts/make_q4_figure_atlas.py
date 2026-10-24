"""Generate Chinese Q4 visual-method atlas PDF (optional extra dependency: reportlab)."""
from pathlib import Path
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4, landscape,portrait
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.pdfbase import pdfmetrics
from reportlab.lib.colors import HexColor
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs'/'Q4_图像问题诊断与论文图版_v1.2.pdf'
DATA=ROOT/'outputs'/'q4'
USGS='https://www.usgs.gov/media/images/2026-nepal-debris-avalanche-and-flash-flood-map'
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
NAVY=HexColor('#18374B');SECONDARY=HexColor('#597082')

def fit(c,path,x,y,w,h):
    im=Image.open(path)
    iw,ih=im.size
    s=min(w/iw,h/ih)
    aw,ah=iw*s,ih*s
    c.drawImage(str(path),x+(w-aw)/2,y+(h-ah)/2,width=aw,height=ah,mask='auto',preserveAspectRatio=True)

def t(c,x,y,lines,size=10,leading=17,color=NAVY,max_chars=95):
    c.setFillColor(color);c.setFont('STSong-Light',size)
    ty=y
    for line in lines:
        while len(line)>max_chars:
            c.drawString(x,ty,line[:max_chars]);ty-=leading;line=line[max_chars:]
        c.drawString(x,ty,line);ty-=leading
    return ty

def page(c,pagesize,section,figfile,figure_area,statements,page_number,source='TPDC 官方数据集④、⑥、⑦'):
    W,H=pagesize;c.setPageSize(pagesize)
    c.setFillColor(HexColor('#f5f8f9'));c.rect(0,H-8,W,8,fill=1,stroke=0)
    c.setFillColor(NAVY);c.setFont('STSong-Light',18)
    c.drawString(36,H-39,'MathorCup 2026 - Q4 模型边界与可视化修订')
    c.setFillColor(SECONDARY);c.setFont('STSong-Light',10)
    c.drawString(36,H-57,section)
    fit(c,DATA/figfile,*figure_area)
    pty=figure_area[1]-13
    pty=t(c,36,pty,statements,size=9.3,leading=15,color=NAVY,max_chars=89 if W>700 else 62)
    c.setStrokeColor(HexColor('#d3dfe4'));c.line(36,55,W-36,55)
    c.setFillColor(SECONDARY);c.setFont('STSong-Light',8.1)
    c.drawString(36,40,f'资料：{source}   |   编程复现：附带的 Python 代码和数值结果')
    c.setFillColor(HexColor('#146884'));c.drawString(36,25,'独立实际灾害地图（USGS 2026-08-27，Public Domain）- 点击打开')
    c.linkURL(USGS,(36,20,400,35),relative=0,thickness=0)
    c.setFillColor(SECONDARY);c.drawRightString(W-36,27,str(page_number))
    c.showPage()

c=canvas.Canvas(str(OUT),pagesize=landscape(A4))
WL,HL=landscape(A4)
page(c,(WL,HL),'图 A - 真实地理底图（参考附件成品论文第35页地形传播图的地理表达）',
     'Q4_00_actual_DEM_geographic_context.png',(27,145,WL-54,384),[
       '地图来自官方 DEM、冰川边界与历史沟谷，而不是装饰性合成卫星图。橙色虚线为本团队的候选传播线。',
       '二维计算区域的位置用矩形标出；真实灾后影响范围可另与 USGS 公布的地图独立核对。'],1)
page(c,(WL,HL),'图 B - 原来3.3分钟空白的时间原因（绝对时间 vs 二维首次入流）',
     'Q4_06A_onset_and_wet_area_diagnostic.png',(30,180,WL-60,330),[
       '一维流体到达二维耦合边界发生在约343.2s（5.72min），所以旧3.3min无流深并非求解器错误。',
       '旧6.7min图片跨过了首次高流量入流阶段。改进后直接显示流量和湿润面积的连续演化。'],2)
page(c,(WL,HL),'图 C - 相同范围、同一色标、统一大小的十帧动态流深图',
     'Q4_06_local_2d_10frames.png',(24,155,WL-48,370),[
       '用二维首次显著入流后的关键时刻代替原来0–30min等间隔帧。所有子图的地理范围、色阶和边框尺度一致。',
       '第一帧约5.9min已出现上游入流；不同深度采用共同对数色标以避免帧间视觉误导。'],3)
WP,HP=portrait(A4)
page(c,(WP,HP),'图 D - 2×2 因素试验：位置与初始动量对分叉的作用',
     'Q4_06C_2x2_inlet_boundary_factorial.png',(20,165,WP-40,590),[
       '两因素分别为入流位置（旧侧坡/谷底低槽）与入流动量（静止/沿沟），四组模拟共享同一 DEM 和流量。',
       '谷底注入显著削弱上游侧枝；单独赋予定向动量难以消除侧坡注入时的分叉。',
       '因此原双分支应优先归因于不合理的数值入流边界，而非直接解释为真实天然分流。'],4)
page(c,(WL,HL),'图 E - 二维最大流深与首次淹没时间',
     'Q4_05_local_2d_depth_and_arrival.png',(30,165,WL-60,355),[
       '左图展示最大等效流深，右图展示首次超过0.2m的时刻，统一空间范围但分别使用各自有单位的物理色标。',
       '修正二维入流位置后的累计淹没范围约1.06km²；旧侧坡均匀入流为1.23km²。两者均为物理情景结果。'],5)
c.save();print(OUT)
