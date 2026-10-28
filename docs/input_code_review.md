# 四问代码审查总结

> **审查范围：**仅对照官方赛题与代码、配置及代码生成的结果文件；不依据论文表述。本次为静态审查，未重跑程序或测试。

## 总体结论

代码包含四问的主要计算流程：Q1灾前异常排序、Q2触发条件评估、Q3条件超阈预报、Q4一维级联与局部二维仿真。发现的问题主要在模型是否充分回答题目、结果可靠性及代码交付一致性；不能仅凭现有计算把Q2阈值解释为物理失稳临界值，把Q3结果解释为泥石流发生概率，或把Q4结果解释为已准确复现事件。

## 按严重程度列出的问题

### 高：影响题目结论

1. **Q3预测对象与题目要求存在差距。**
   题目要求预报泥石流触发时间并给出95%区间。代码实际预测的是气象条件指数首次超过历史阈值的时间，回测目标也是气象条件是否超阈，不是泥石流是否发生。因此，超阈概率、区间覆盖率不能直接解释为泥石流发生概率或泥石流触发时间的验证结果。本次灾害时间的实际覆盖率没有可计算值。
   - 代码：[q3_task_aligned.py](代码与论文/代码/MathorCup_B_v2/src/mathorcup/q3_task_aligned.py#L126)、[warning_validation.py](代码与论文/代码/MathorCup_B_v2/src/mathorcup/warning_validation.py#L85)
   - 输出：[Q3_primary_competition_answer.json](代码与论文/代码/MathorCup_B_v2/outputs/q3/Q3_primary_competition_answer.json)

2. **Q2使用的历史分位阈值不是物理失稳临界值。**
   代码以历史q90构造降雨和正积温阈值，再组合冰川静态评分与水热门控。它能衡量条件相对历史是否偏高，但不能证明达到该值就会触发泥石流。静态权重为预设值，模型没有用真实灾害发生标签校准为成灾概率。能量余量通道虽有单独计算，却没有进入主要触发指数。
   - 代码：[q2_trigger_evaluation.py](代码与论文/代码/MathorCup_B_v2/src/mathorcup/q2_trigger_evaluation.py#L35)、[q2_spatial_mechanism.py](代码与论文/代码/MathorCup_B_v2/src/mathorcup/q2_spatial_mechanism.py#L169)
   - 输出：[Q2_primary_trigger_assessment.json](代码与论文/代码/MathorCup_B_v2/outputs/q2/Q2_primary_trigger_assessment.json)、[Q2_expanded_mechanism_answer.json](代码与论文/代码/MathorCup_B_v2/outputs/q2/Q2_expanded_mechanism_answer.json)

3. **Q4到达时间与事件约束有明显差距。**
   配置记录汇入口时间约130秒、汇入口至口岸约350秒，合计约480秒。模型输出上游到达汇入口约80.8秒、下游传播约261.7秒、源区至口岸约342.5秒。下游传播时间未达到350秒；源区至口岸模拟时间也比按配置合计的480秒早约137.5秒。比较时应统一时间起点。
   - 配置：[case_evidence.yaml](代码与论文/代码/MathorCup_B_v2/configs/case_evidence.yaml#L29)
   - 代码：[q4_cascade_reconstruction.py](代码与论文/代码/MathorCup_B_v2/src/mathorcup/q4_cascade_reconstruction.py#L166)
   - 输出：[Q4_cascade_reconstruction.json](代码与论文/代码/MathorCup_B_v2/outputs/q4/Q4_cascade_reconstruction.json)

### 中：影响模型可靠性或结论范围

4. **Q4一维结果对网格分辨率敏感。**
   沿程网格从15米细化到7.5米后，峰值流量由约20.1万升至约26.0万立方米/秒；30米网格结果约为14.3万。峰值流量尚未稳定。质量守恒误差小只能说明离散计算的体积收支闭合，不能证明流量、流深和到达时间已经收敛。
   - 输出：[Q4_reconstruction_1d_grid_experiment.csv](代码与论文/代码/MathorCup_B_v2/outputs/q4/Q4_reconstruction_1d_grid_experiment.csv)

5. **Q4二维结果仅代表口岸附近局部区域。**
   二维模型接收一维通量进行单向耦合，局部计算没有反向更新上游；接受的通量只是下游可动方量的一部分。因此，局部湿润面积不能代表整条通道的完整淹没范围。
   - 代码：[q4_reconstructed_spatial.py](代码与论文/代码/MathorCup_B_v2/src/mathorcup/q4_reconstructed_spatial.py#L13)
   - 输出：[Q4_reconstruction_2d_summary.json](代码与论文/代码/MathorCup_B_v2/outputs/q4/Q4_reconstruction_2d_summary.json)

6. **Q4部分参数缺少独立测量约束。**
   下游有效断面尺度标定值达到搜索上界3.0，且不是实测河宽；这说明模拟结果仍依赖参数设定。模型主要复现物质释放后的传播过程，没有完整模拟失稳起动、碎屑沿程侵蚀补给及冰岩崩转化为泥石流的全部机制。
   - 配置及参数：[calibrated_q4.yaml](代码与论文/代码/MathorCup_B_v2/configs/calibrated_q4.yaml)
   - 代码：[q4_cascade_reconstruction.py](代码与论文/代码/MathorCup_B_v2/src/mathorcup/q4_cascade_reconstruction.py#L211)

7. **Q2气象条件的源区代表性有限。**
   动态积温使用冰川区中位海拔；天气预报使用的模式格点距请求源区约8.38公里，并依据假定的发布时间延迟判断是否可用。它们可作为空间代理，不能直接当作源区现场实测。
   - 代码：[q2_external.py](代码与论文/代码/MathorCup_B_v2/src/mathorcup/q2_external.py#L101)
   - 输出：[Q2_primary_trigger_assessment.json](代码与论文/代码/MathorCup_B_v2/outputs/q2/Q2_primary_trigger_assessment.json)

8. **Q1排序表示灾前异常程度，不等于统计因果关系。**
   代码比较单次事件窗口与历史同期，计算异常分数并选前五项。这适合识别灾前较异常的观测因素；没有多次灾害事件和非灾害对照标签，不能据此证明这些因素导致泥石流。
   - 代码：[q1_event_signatures.py](代码与论文/代码/MathorCup_B_v2/src/mathorcup/q1_event_signatures.py#L45)

### 代码控制与交付一致性

9. **strict-route没有阻止所有级联仿真。**
   该参数会影响常规路径分支，但程序之后仍调用级联重建。因此，它不能作为“禁止使用未核验路径”的完整保护条件。
   - 代码：[run.py](代码与论文/代码/MathorCup_B_v2/run.py#L79)

10. **工作区存在不完全一致的代码副本。**
    两个 MathorCup_B_v2 目录的Python文件数量和部分文件内容不同。核心求解模块大体一致，但最终提交前仍应指定唯一代码版本，并确保清单、输入、输出和版本对应。

## 逐问判断

| 问题 | 判断 |
|---|---|
| Q1 | 异常排序流程基本完整；结论应限于灾前异常，不应扩展为因果结论。 |
| Q2 | 有静态物源指数、水热阈值和时空敏感性分析；阈值是历史气候代理，尚非物理失稳临界值。 |
| Q3 | 能预测气象条件超阈和相应预警；尚未验证泥石流发生概率或实际触发时间区间。 |
| Q4 | 已建立并运行事件传播仿真；到达时间拟合、网格稳定性、局部二维范围及参数约束仍需谨慎解释。 |

**总体结论：**代码实现了四问的主要计算，但Q2和Q3的模型含义不能扩大到真实失稳阈值或泥石流发生概率；Q4的事件时间拟合和网格稳定性是当前最明确的定量问题。此次为静态审查，未重跑代码或测试。
