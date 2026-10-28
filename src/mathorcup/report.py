"""Question-ordered model answers and the existing figure gallery."""
import html,json
from .common import write_json,FIGURE_TITLES


def make_report(ctx):
    def j(path):return json.loads((ctx.out/path).read_text())
    q1=j('q1/Q1_direct_answer.json');q2=j('q2/Q2_sequential_model_contract.json')
    q3=j('q3/Q3_primary_competition_answer.json');q4=j('q4/Q4_cascade_reconstruction.json')
    geo=j('q4/Q4_reconstruction_geometry.json');two=j('q4/Q4_reconstruction_2d_summary.json')
    assert q2['source_geometry_sha256']==q3['source_geometry_sha256']==geo['geometry_contract_sha256']
    summary={'version':'v2.3.1','extra_data_only_question_2':True,
        'question_order':['Q1 features','Q2 trigger assessment and static source/terrain','Q3 conditional event time and warning','Q4 flow simulation'],
        'Q1_answer_mechanism_top5':q1['mechanism_top5'],'Q2_inherited_features':q2['Q1_features'],
        'Q2_critical_score':q2['critical_score'],'Q2_score_at_issue':q2['current_score'],
        'Q3_target':q3['forecast_target'],'Q3_conditional_interval_bjt':q3['conditional_95pct_interval_times_bjt'],
        'Q3_warning_time_bjt':q3['operational_warning_time_bjt'],'Q3_warning_lead_h':q3['operational_warning_lead_h'],
        'Q3_case_reference_PICP':q3['case_reference_coverage'],'Q3_empirical_disaster_PICP':None,
        'Q3_validated_complete_event_forecast':False,'Q3_model_based_question_outputs_available':True,
        'Q3_spatial_distance_range_km':q3['spatial_distance_range_km'],'Q3_uses_Q4_results':False,
        'Q4_uses_Q2_geometry':True,'shared_source_geometry_sha256':q2['source_geometry_sha256'],
        'Q4_reconstructed_route_km':geo['length_m']/1000,'Q4_arrival_s':q4['downstream']['arrival_s'],
        'Q4_peak_discharge_m3_s':q4['downstream']['peak_discharge_m3_s'],
        'Q4_cascade_volume_residual':q4['cascade_volume_residual'],'Q4_local_2d_grid_experiments':two['grid_experiments']}
    write_json(ctx.out/'SUMMARY_STATUS.json',summary)
    groups=[]
    for folder,title in [('q1','Q1 关联特征'),('q2','Q2 触发评估与源区'),('q3','Q3 时间预报、预警与空间距离'),('q4','Q4 传播仿真与原有验证'),('joint','Q3时间与Q4传播的组合')]:
        cards=[]
        for p in sorted((ctx.out/folder).glob('*.png')):
            label=FIGURE_TITLES.get(p.stem,p.stem)
            cards.append(f'<a class="card" href="{folder}/{html.escape(p.name)}"><img src="{folder}/{html.escape(p.name)}" loading="lazy"><span>{html.escape(label)}</span></a>')
        groups.append('<h2>'+title+'</h2><div class="grid">'+''.join(cards)+'</div>')
    interval=' 至 '.join(q3['conditional_95pct_interval_times_bjt']);coverage=q3['case_reference_coverage']
    notice=f"按Q1→Q2→Q3→Q4顺序计算。模型条件95%事件时间区间：{interval}；预警：{q3['operational_warning_time_bjt']}；题面单例覆盖{coverage['n_covered']}/{coverage['n_references']}；静态通道距离0-{q2['source_to_port_km']:.2f}km。时间密度与评估强度成正比是明确建模假设，单例覆盖不等于多事件泛化验证。"
    body='''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>MathorCup B v2.3.1</title><style>body{max-width:1280px;margin:auto;padding:28px;font-family:system-ui,sans-serif;color:#233642;background:#f2f6f8}p{line-height:1.8}.note{background:white;padding:18px;border-left:5px solid #248f81}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:15px}.card{background:white;padding:12px;text-decoration:none;color:#233642;border-radius:6px}img{width:100%;height:260px;object-fit:contain}span{font-size:12px;overflow-wrap:anywhere}</style><h1>四问顺序建模与结果 v2.3.1</h1>'''
    (ctx.out/'results_gallery.html').write_text(body+'<p class="note">'+html.escape(notice)+'</p>'+''.join(groups)+'</html>',encoding='utf-8')
    return summary
