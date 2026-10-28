"""Resolve the ten supplied review findings against this exact executable release."""
from __future__ import annotations
import json
from .common import write_json


def build_review_verdict(ctx):
    def read(path):return json.loads((ctx.out/path).read_text())
    baseline_path=ctx.folder('quality')/'review_baseline.json'
    if not baseline_path.exists():
        import yaml
        write_json(baseline_path,yaml.safe_load((ctx.root/'configs/review_baseline.yaml').read_text()))
    baseline=read('quality/review_baseline.json');q4=read('q4/Q4_cascade_reconstruction.json')
    grid=read('q4/Q4_numerical_resolution_diagnostics.json');q3=read('q3/Q3_primary_competition_answer.json')
    geometry=read('q4/Q4_reconstruction_geometry.json')
    rows=[
      (1,'Q3预测对象','成立：属于预测对象与验证含义',
       '主模型输出条件首次超限时间，历史检验目标也为气象条件；本次题面时刻单例命中另列。',
       '统一图题、字段和中位时刻的离散支持；当前天气驱动回放与历史相似日检验分别登记，保留无超限的零概率输出。'),
      (2,'Q2临界条件与能量通道','部分成立：历史分位定义合理，能量衔接确有遗漏',
       'q90是明确的历史参考阈值；静态权重是声明的指数设计参数。',
       '把实测正能量余量纳入主动态指数，补算静态权重与高程递减率敏感性；未来能量只取实测历史贡献的下界。'),
      (3,'Q4时间拟合','成立：旧版数值偏差与单侧时间约束均已复算',
       f"旧版汇入口至口岸{baseline['Q4_confluence_to_port_s']:.2f}s；本次{q4['downstream_transit_s']:.2f}s，文献约350s。",
       '将上游时间改为双侧标定约束，统一源区释放时钟，分别导出前锋、洪峰、5%与50%累计量时间。文献参照保留为近似重建值。'),
      (4,'Q4网格敏感性','成立：体积守恒与分辨率响应是两项不同检查',
       f"本次两档较细网格峰值差{grid['peak_discharge_relative_change']:.2%}，到达时间差{grid['arrival_relative_change']:.2%}。",
       '下游改用静水重构MUSCL-HLL与SSP-RK2；网格试验保持同一入流、参数及两端位置，分别报告各量变化。'),
      (5,'Q4局部二维范围','事实正确，局部单向嵌套本身是合理实现',
       '一维覆盖源区至口岸，二维覆盖口岸邻域；输入为上游实际到达局部断面的累计通量。',
       '新增全程一维过程的地图投影，并与局部二维湿区分列，记录各计算域及对应体积收支。'),
      (6,'Q4参数测量约束','事实正确，标定参数不等于程序错误',
       f"本次有效断面尺度{q4['stage2_parameters']['section_capacity_scale']:.4f}，边界状态见参数诊断。",
       '保留观测、文献重建、标定参数、数值设置的角色；自动记录搜索上下界与是否触界，不用改写观测或强行消除残差。'),
      (7,'Q2源区气象代表性','成立：属于空间代理与资料可用性约定',
       '站点、冰川代表高程与预报格点的空间尺度不同，运行时延按已有资料声明检查。',
       '增加冰川高程q10/q50/q90与三档递减率的九组合重算；保留实际格点偏移和预报来源。'),
      (8,'Q1异常与因果','判断合理；当前异常排序方法可保留',
       '题面要求关联性排序，单次事件与历史匹配可定义异常关联；该评分本身没有估计因果效应。',
       '保持全部候选、原始前五、机制代表项和年份敏感性共同输出；用异常关联说明评分含义。'),
      (9,'strict-route保护','成立：是真实的程序控制缺口',
       '旧入口在常规路径分支后仍能进入级联重建。',
       '在清理旧图前校验，所有级联入口执行同一检查；严格模式实际使用核验路径，冲突命令参数直接拒绝。'),
      (10,'代码副本不一致','在本次最新交付包中未复现',
       f"压缩包只有一个代码根目录；{baseline['remote_python_files_matched']}/{baseline['remote_python_files_checked']}个Python文件与基准GitHub版本逐一相同。",
        '继续使用唯一源码根目录，结果清单记录输入、源码、参数、图表及版本的校验摘要。'),
      (11,'另查出：源区释放投影','成立：原PCA平移改变了初始释放势能',
        f"旧版释放平均床高{baseline['Q4_release_assigned_mean_bed_m']:.1f}m，真实像元均值{baseline['Q4_selected_source_pixels_mean_elevation_m']:.1f}m；修正后释放平均床高{geometry['assigned_release_mean_bed_elevation_m']:.1f}m，逐像元投影平均绝对高程差{geometry['release_projection_mean_abs_elevation_error_m']:.1f}m。",
        '改为源区上部真实像元起始的地形路径和UTM三维最近点投影；逐像元导出高程差，并用修正后的释放位置重新标定与求解。')]
    verdict=[dict(zip(['id','finding','judgement','evidence','action'],r)) for r in rows]
    result={'release':'v2.2','review_basis':'official question, delivered source, input data and rerun experiments',
        'reference_commit':baseline['baseline_github_commit'],'findings':verdict,
        'q4_time_and_peak_diagnostics':'q4/Q4_calibration_residuals.csv',
        'q4_resolution_diagnostics':grid,'forecast_target':q3['forecast_target'],
        'historical_validation_target':'meteorological condition exceedance',
        'single_task_reference_hit':q3['single_task_reference_hit']}
    write_json(ctx.folder('quality')/'review_verdicts.json',result)
    return result
