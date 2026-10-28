"""Small, transparent helpers shared across Q1--Q4."""
from __future__ import annotations
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LABELS={
 'temp_1m':'1 m 气温','temp_4m':'4 m 气温','temp_10m':'10 m 气温',
 'wind_1m':'1 m 风速','wind_4m':'4 m 风速','wind_10m':'10 m 风速',
 'rh_1m':'1 m 相对湿度','rh_4m':'4 m 相对湿度','rh_10m':'10 m 相对湿度',
 'rh_20m':'20 m 相对湿度', 'humidity_gradient':'空气湿度垂向梯度',
 'vector_shear_2_20m':'2–20 m矢量风切变', 'bulk_richardson':'总体Richardson数',
 'soil_temp_gradient':'土壤温度垂向梯度', 'soil_heat_gradient':'土壤热通量梯度',
 'soil_freeze_depth_layers':'土壤冻结层数',
 'energy_balance_residual':'地表能量收支余量',
 'vapour_pressure':'水汽分压', 'water_vapour_density':'水汽密度','rh_gradient':'湿度垂向梯度',
 'vector_wind_shear':'矢量风切变','wind_profile_exponent':'风速幂律指数',
 'wind_direction_turn':'上下层风向转角','bulk_richardson_proxy':'Richardson数近似',
 'net_shortwave':'净短波辐射','net_longwave':'净长波辐射',
 'sw_up':'向上短波辐射','lw_down':'向下长波辐射','lw_up':'向上长波辐射',
 'vwc_0p2':'0.2 m含水率','vwc_0p4':'0.4 m含水率','vwc_1p6':'1.6 m含水率',
 'soilT_0p2':'0.2 m土壤温度','soilT_0p4':'0.4 m土壤温度',
 'soilT_0p8':'0.8 m土壤温度','soilT_1p6':'1.6 m土壤温度',
 'vwc_deep_gradient':'深层土壤含水梯度',
 'soil_temperature_gradient':'土壤温度垂向梯度',
 'soil_heat':'0.1 m土壤热通量','soil_heat_0p2':'0.2 m土壤热通量',
 'temperature_inversion':'温度逆温','soil_surface_frozen':'浅层冻结指示',
 'turbulent_total_heat':'感热与潜热之和','surface_energy_residual':'地表能量闭合剩余',
 'bowen_ratio':'Bowen比','z_over_L_proxy':'Monin–Obukhov稳定度近似',
 'tau':'动量通量TAU','tstar':'湍流温度尺度TSTAR','sonic_temp':'声温',
 'flux_ux':'湍流x方向速度','flux_uy':'湍流y方向速度','flux_uz':'湍流垂直速度',
 'water_density':'水汽密度',
 'temp_2m':'2 m 气温', 'temp_20m':'20 m 气温', 'temp_gradient':'气温垂向梯度',
 'wind_2m':'2 m 风速','wind_20m':'20 m 风速','wind_gradient':'风速垂向梯度',
 'rh_2m':'2 m 相对湿度','precip_mm':'小时降水',
 'pressure':'气压','net_radiation':'净辐射', 'sw_down':'向下短波辐射',
 'albedo':'反照率','vwc_0p1':'表层含水率','vwc_0p8':'深层含水率',
 'vwc_gradient':'土壤水分梯度','soilT_0p1':'浅层土壤温度',
 'sensible_heat':'感热通量','latent_heat':'潜热通量','ustar':'摩擦速度',
 'tke':'湍流动能','vpd':'通量塔VPD','vpd_2m':'2米气象VPD','pdd_hour':'正积温贡献',
 'temp_max_24h':'24h最高气温','pdd_7d':'7日正积温','rain_24h':'24h累计降水',
 'rain_72h':'72h累计降水','rain_7d':'7日累计降水'
}
COLORS=['#205C92','#DA783F','#248F81','#9276A5','#C7AA41','#587A89','#CB5564','#6B8D4E']

# Stable paper-style titles used by the plotting functions and the figure gallery.
# Methodological assumptions belong in figure captions/README, not in the heading.
FIGURE_TITLES={
 'Q1_14_mechanisms_and_extreme_weather':'五项互补机制与极端天气识别',
 'Q2_18_energy_water_mechanism_channels':'降雨、前期热量与能量通道评估',
 'Q3_11_warning_time_and_lead_distribution':'预警时间与提前量分布',
 'Q3_12_reconstructed_spatial_warning':'源区至口岸的模型空间警戒范围',
 'Q4_14_drainage_route_and_cross_sections':'源区至口岸的地形与断面',
 'Q4_15_calibrated_cascade_process':'分阶段冰岩崩与泥石流传播过程',
 'Q4_16_event_constraint_comparison':'事后标定的定量对照',
 'Q4_17_calibration_parameter_experiments':'同一输入下的标定参数试验',
 'Q4_18_reconstructed_2d_depth_speed_arrival':'口岸二维传播强度与到达时间',
 'Q4_19_reconstructed_2d_six_frames':'真实通量驱动的六个传播阶段',
 'Q4_20_cascade_volume_and_grid_experiments':'级联体积守恒与网格试验',
 'Q4_21_cross_sections_and_1d_grid_experiment':'原始横断面与数值分辨率试验',
 'Q4_22_time_alignment_and_residuals':'事件时间起点与全部标定残差',
 'Q4_23_whole_route_six_stages':'源区至口岸的全程传播过程投影',
 'Q4_24_resolution_response_comparison':'修订前后网格敏感性的完整比较',
 'Q4_25_release_geometry_and_elevation':'源区三维投影与初始释放势能核查',
 'Q2_19_source_proxy_and_weight_sensitivity':'源区空间代理与静态权重敏感性',
 'Q34_trigger_and_arrival_time':'触发条件与下游到达时间分布',
 'Q1_00_data_coverage':'逐小时观测数据有效率',
 'Q1_01_event_hourly_timeseries':'吉隆泥石流灾前气象要素变化',
 'Q1_02_boundary_layer_profiles':'高寒山区近地层风温垂直结构',
 'Q1_03_top15_anomaly_ranking':'灾前多源气象因素异常稀有性排序',
 'Q1_04_multiscale_anomaly_heatmap':'多时间尺度气象异常分布',
 'Q1_05_reference_year_sensitivity':'关键特征排序的年份敏感性',
 'Q1_06_pseudo_event_controls':'灾前气象异常与历史同期对照',
 'Q1_07_county_2005_2025_summer_climate':'吉隆及周边地区夏季气温变化',
 'Q1_08_boundary_layer_seasonal_structure':'近地层稳定结构的季节变化',
 'Q1_09_radiation_and_turbulent_energy':'灾前地表辐射与热量通量变化',
 'Q1_10_multidepth_soil_profiles':'灾前土壤多深度温湿廓线',
 'Q2_01_history_gullies_and_oof_roc':'历史沟谷成因分布与判别模型ROC曲线',
 'Q2_02_spatial_holdout_auc':'灾种判别模型空间交叉验证',
 'Q2_03_elevation_drop_distributions':'不同成因沟谷高差分布',
 'Q2_04_climate_condition_timeseries':'降雨与融水触发条件时间变化',
 'Q2_05_glacier_physical_indicators':'冰川厚度、坡度与驱动剪应力分布',
 'Q2_06_dem_landcover_history_overlay':'研究区地形与土地覆盖分布',
 'Q2_07_glacier_historic_thinning_optional':'冰川历史高程变化空间分布',
 'Q2_08_trigger_conditions_time_series':'降雨与融水触发条件变化',
 'Q2_09_source_spatial_condition_index':'冰川源区地形物源条件分布',
 'Q2_10_spatial_resolution_reliability':'触发评估空间尺度可靠性',
 'Q2_11_temperature_lapse_sensitivity':'源区热力条件递减率敏感性',
 'Q2_12_forecast_weather_conditions':'源区历史起报天气预报条件',
 'Q2_13_glacier_velocity_and_release_candidates':'冰川历史运动与候选失稳区空间特征',
 'Q2_14_landcover_2020_2023_comparison':'2020年与2023年土地覆盖对比',
 'Q3_01_precursor_and_scenario_hazard':'灾前气象条件及触发风险演变',
 'Q3_02_95percent_conditional_time_interval':'触发条件首次超限时间与95%预测区间',
 'Q3_03_warning_lead_and_backtest':'泥石流预警时间与气象条件回测',
 'Q3_04_parameter_sensitivity':'预警模型参数敏感性分析',
 'Q3_05_spatial_warning_distance':'泥石流传播廊道与空间预警距离',
 'Q3_06_validation_metrics':'气象条件时间区间与预警指标验证',
 'Q3_07_feature_group_ablation':'Q1气象特征与Q2物理条件的贡献比较',
 'Q4_10_impact_pressure_distribution':'二维流动峰值动压分布',
 'Q4_00_actual_DEM_geographic_context':'吉隆冰崩—泥石流研究区地形图',
 'Q4_00B_local_true_DEM_detail':'下游二维计算区域地形图',
 'Q4_01_route_geolocation_quality_gate':'冰崩—泥石流传播路径及地形纵剖面',
 'Q4_02_1d_runout_timespace':'一维泥石流传播及地形纵剖面',
 'Q4_03_1d_hydrograph_and_mass_conservation':'耦合断面流量与体积守恒检验',
 'Q4_04_1d_12frames':'一维泥石流流深演进过程',
 'Q4_05_local_2d_depth_and_arrival':'泥石流二维最大流深与到达时间',
 'Q4_06A_onset_and_wet_area_diagnostic':'耦合断面入流及淹没面积变化',
 'Q4_06_local_2d_10frames':'泥石流二维流深演进过程',
 'Q4_06B_inlet_position_ablation':'二维入流位置对传播范围的影响',
 'Q4_06C_2x2_inlet_boundary_factorial':'二维入流位置与动量组合试验',
 'Q4_07_parameter_sensitivity':'泥石流动力学参数敏感性分析',
 'Q4_08_volume_friction_phase_diagram':'物源体积与摩擦系数对传播的影响',
 'Q4_09_grid_resolution_reliability':'二维模拟网格分辨率敏感性',
 'Q34_joint_uncertainty_MATHEMATICAL_DEMO':'触发时间与下游到达时间联合分布'
}

def configure_plotting():
    from matplotlib import font_manager
    fonts=Path(__file__).resolve().parents[2]/'assets'/'fonts'
    for font in fonts.glob('*.otf'):
        font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({
        'font.family':'sans-serif','font.sans-serif':['Noto Sans CJK SC','Noto Sans CJK JP','DejaVu Sans'],
        'axes.unicode_minus':False,'font.size':10,'axes.titlesize':11.5,'axes.labelsize':9.5,
        'figure.dpi':110,'savefig.dpi':170,'axes.spines.top':False,'axes.spines.right':False,
        'axes.grid':False,'grid.alpha':0.16,
        'figure.facecolor':'white','axes.facecolor':'white'
    })

def savefig(fig,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    canonical=FIGURE_TITLES.get(path.stem)
    if canonical:
        if getattr(fig,'_suptitle',None) is not None:
            fig._suptitle.set_text(canonical)
            fig._suptitle.set_fontsize(13)
            fig._suptitle.set_fontweight('medium')
        elif len(fig.axes)==1:
            fig.axes[0].set_title(canonical,fontsize=12,fontweight='medium',pad=11)

    import os
    if os.environ.get('MATHORCUP_DISABLE_REDRAW_EXPORT')!='1':
        from .redraw_export import export_figure_data
        export_figure_data(fig,path)
    fig.savefig(path,dpi=170,bbox_inches='tight',facecolor='white')
    fig.savefig(path.with_suffix('.svg'),bbox_inches='tight',facecolor='white')
    plt.close(fig)

def write_json(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    def clean(o):
        if isinstance(o,(np.integer,np.floating)):return o.item()
        if isinstance(o,pd.Timestamp):return o.isoformat()
        if isinstance(o,np.ndarray):return o.tolist()
        if isinstance(o,Path):return str(o)
        raise TypeError(str(type(o)))
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=clean),encoding='utf-8')

def robust_z(x,ref):
    a=np.array(ref,dtype=float);a=a[np.isfinite(a)]
    if len(a)<15 or not np.isfinite(x):return np.nan
    med=np.median(a);mad=np.median(np.abs(a-med));scale=1.4826*mad
    if scale<1e-12:return np.nan
    return float((x-med)/scale)

def quantile_tail(x,ref):
    a=np.array(ref,dtype=float);a=a[np.isfinite(a)]
    if not np.isfinite(x) or len(a)<15:return np.nan
    left=(1+int(np.sum(a<=x)))/(len(a)+1)
    right=(1+int(np.sum(a>=x)))/(len(a)+1)
    return float(min(1.0,2*min(left,right)))

def record_warning(message):
    print('[LIMIT]',message)

def find_first(base,pattern):
    matches=list(Path(base).rglob(pattern))
    if not matches:raise FileNotFoundError(f'Not found: {pattern} under {base}')
    return matches[0]

class Context:
    def __init__(self,root,config=None):
        self.root=Path(root).resolve();self.raw=self.root/'data'/'raw';self.processed=self.root/'data'/'processed'
        self.out=self.root/'outputs';self.cfgpath=Path(config) if config else self.root/'configs'/'base.yaml'
        self.cfg=yaml.safe_load(self.cfgpath.read_text(encoding='utf-8'))
        for d in [self.raw,self.processed,self.out,self.root/'events']:d.mkdir(parents=True,exist_ok=True)
    def folder(self,name):
        out=self.out/name;out.mkdir(parents=True,exist_ok=True);return out
