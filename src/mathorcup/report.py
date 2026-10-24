"""Current-output-only gallery and task-status summary; no frozen release numbers."""
import json,html
import pandas as pd
from .common import write_json,FIGURE_TITLES


def make_report(ctx):
    def j(p):return json.loads((ctx.out/p).read_text())
    q1=pd.read_csv(ctx.out/'q1/Q1_top5_with_stability.csv');q2=j('q2/Q2_model_comparison.json');phys=j('q2/Q2_primary_trigger_assessment.json')
    q3=j('q3/Q3_primary_competition_answer.json');space=j('q3/Q3_spatial_range_results.json');main=j('q4/Q4_1d_settings_and_summary.json');two=j('q4/Q4_2d_status.json')
    summary={'version':'v2.1','official_data_first':True,'extra_data_only_question_2':True,
      'Q1_top5_event_anomalies':q1.feature.tolist(),'Q1_normalized_feature_ranking_is_not_proven_causality':True,
      'Q1_raw_top5_separate_from_mechanism_representatives':True,
      'Q2_spatial_OOF_logistic_auc':q2['L2_Logistic_4geo_summary']['auc_oof'],
      'Q2_all_models_included':True,'Q2_trigger_thresholds':phys['thresholds'],
      'Q2_external_forecast_status':phys['external_forecast_status'],
      'Q3_issue_bjt':q3['prediction_issue_bjt'],'Q3_conditional_95_time_interval':q3['conditional_95pct_interval_times_bjt'],
      'Q3_warning_lead_hours':q3['operational_warning_lead_h'],
      'Q3_warning_time_bjt':q3['operational_warning_time_bjt'],'Q3_forecast_warning_issued':q3['forecast_warning_issued'],
      'Q3_active_monitoring_warning_at_issue':q3['active_warning_at_issue'],
      'Q3_disaster_event_PICP':q3['empirical_disaster_time_PICP'],'Q3_meteorological_proxy_PICP':q3['proxy_backtest']['proxy_interval_coverage'],
      'Q3_fit_calibration_test_years':[[2021,2022,2023],[2024],[2025]],
      'Q3_spatial_conditional_envelope_km2':space['combined_warning_envelope_area_km2'],
      'Q3_primary_reached_distance_km':space['primary_reached_distance_km'],
      'Q4_downstream_reached_by_main_1D':main['mass']['first_arrival_to_downstream_end_s'] is not None,
      'Q4_2d_origin':two.get('inflow_source'),'Q4_2d_primary_coupled':two.get('primary_1d_coupled',False),
      'Q4_1D_mass_relative_error':main['mass']['relative_absolute_mass_residual'],
      'Q4_grid_status':j('q4/Q4_grid_resolution_status.json')['status'],
      'scope':'observed mechanism analysis / multichannel trigger conditions / retrospective forecast / calibrated cascade reconstruction'}
    reconstruction=ctx.out/'q4/Q4_cascade_reconstruction.json'
    if reconstruction.exists():
        rec=j('q4/Q4_cascade_reconstruction.json');geo=j('q4/Q4_reconstruction_geometry.json')
        summary['Q4_reference_single_stage_reached']=summary['Q4_downstream_reached_by_main_1D']
        summary['Q4_downstream_reached_by_main_1D']=rec['downstream']['arrival_s'] is not None
        summary['Q4_1D_mass_relative_error']=max(rec['upstream']['relative_mass_error'],rec['downstream']['relative_mass_error'])
        summary['Q4_2d_origin']='calibrated 1D cumulative interface volume and velocity'
        summary['Q4_2d_primary_coupled']=True
        summary['Q4_grid_status']='COMPLETED_NATIVE_30_AND_60_M'
        summary['Q4_reconstructed_route_km']=geo['length_m']/1000
        summary['Q4_cascade_volume_residual']=rec['cascade_volume_residual']
        summary['Q4_calibrated_port_arrival_s']=rec['downstream']['arrival_s']
        summary['Q4_calibrated_port_peak_discharge_m3_s']=rec['downstream']['peak_discharge_m3_s']
    spatial_answer=ctx.out/'q3/Q3_reconstructed_spatial_answer.json'
    if spatial_answer.exists():
        result=j('q3/Q3_reconstructed_spatial_answer.json')
        summary['Q3_reference_single_stage_envelope_km2']=summary['Q3_spatial_conditional_envelope_km2']
        summary['Q3_spatial_conditional_envelope_km2']=result['warning_envelope_km2']
        summary['Q3_primary_reached_distance_km']=result['source_to_port_corridor_km']
    summary['Q1_answer_mechanism_top5']=j('q1/Q1_direct_answer.json')['mechanism_top5']
    summary['Q2_extra_exceeding_channels']=j('q2/Q2_expanded_mechanism_answer.json')['channels_exceeding_reference']
    write_json(ctx.out/'SUMMARY_STATUS.json',summary)
    items=[]
    for folder,title in [('q1','Q1 全因素、过程与关联'),('q2','Q2 全模型、触发条件与时空可靠性'),('q3','Q3 全检验期、有效预警与时间情景'),('q4','Q4 分阶段复现与参考参数实验'),('joint','触发与传播条件联合情景')]:
        cards=[]
        for p in sorted((ctx.out/folder).glob('*.png')):
            label=FIGURE_TITLES.get(p.stem,p.stem)
            cards.append(f'<a class="card" href="{folder}/{html.escape(p.name)}"><img src="{folder}/{html.escape(p.name)}" loading="lazy"><span>{html.escape(label)}</span></a>')
        items.append('<h2>'+title+'</h2><div class="grid">'+''.join(cards)+'</div>')
    coverage=q3['proxy_backtest']['proxy_interval_coverage']
    content=f'''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>MathorCup B v2.1</title><style>
    body{{max-width:1280px;margin:auto;padding:30px;font-family:system-ui,'Microsoft YaHei',sans-serif;color:#233642;background:#f2f6f8}}h1{{color:#164c68}}p{{line-height:1.8}}.note{{padding:16px;background:#fff;border-left:5px solid #247e7b}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(350px,1fr));gap:16px}}.card{{background:white;border:1px solid #dce4e8;border-radius:8px;padding:12px;color:#264957;text-decoration:none}}img{{width:100%;height:260px;object-fit:contain}}span{{font-size:12px;overflow-wrap:anywhere}}</style>
    <h1>四问代码重算与结果图廊 v2.1</h1><p>每张图均由本版源代码及已保存的数值表生成。完整用途与证据对应见FIGURE_INDEX.csv。</p>
    <p class="note">Q1给出五项互补机制与明确的极端天气识别；Q2补充辐射能量、前期热量和冻融通道；Q3用2021-2023拟合、2024校准、2025检验，气象条件区间覆盖率{coverage:.1%}，本次回放预警提前{q3['operational_warning_lead_h']:.1f}小时。Q4新增DEM汇流路径、实际断面、分阶段方量与真实通量耦合，并以事件研究约束标定和进行网格试验。Q4_14起为分阶段复现主模型，早期编号保留为单阶段参考实验。</p>{''.join(items)}</html>'''
    (ctx.out/'results_gallery.html').write_text(content,encoding='utf-8')
    return summary
