# MathorCup B：吉隆泥石流预报与仿真

按 **Q1特征分析 → Q2触发评估与源区 → Q3时间预报 → Q4传播仿真** 顺序计算。仓库包含完整复现所需的输入、代码、计算结果、85张PNG、SVG、GIS和逐图重绘数据。

## 运行

使用Python 3.12，在仓库根目录执行：

```bash
python -m pip install -r requirements.txt
python run.py --stage all --reference-experiments --no-fetch-q2-external
```

该命令使用已归档的天气输入，依次完成四问，并生成已有对照实验和全部图件。结果写入`outputs/`。

## 四问结果

| 问题 | 当前结果 |
|---|---|
| Q1 关联因素 | 原始前五项：vpd、temp_10m、temp_20m、energy_balance_residual、temp_2m。 |
| Q2 触发评估 | 发布前综合指数0.383，临界水平0.320；机制代表为vpd、temp_10m、energy_balance_residual、sw_down、vwc_gradient。 |
| Q3 触发时间 | 95%时间区间：2026-08-26 09:04:00 至 2026-08-26 14:44:00；中位时刻约2026-08-26 10:48。 |
| Q3 预警与覆盖 | 2026-08-26 09:00:00发布预警，提前1.5h；本次案例区间覆盖率1/1=100%。 |
| Q3 空间距离 | 源区至口岸沿程0-23.715km。 |
| Q4 传播仿真 | 释放后403.47s到达口岸，峰值流量208522m³/s；输出全程一维和口岸局部二维过程。 |

完整模型、参数和实验见[四问作答PDF](docs/MathorCup_B_v2_3_1_question_answers.pdf)。

## 数据与模型约定

官方附件提供MET、RADM、SOIL、FLUX观测、历史温度、冰川、土地覆盖、历史灾点、DEM及人口资料；第二问补充归档ECMWF天气和Sentinel-1形变资料。所有事件钟时使用北京时间，预报发布于2026年8月26日09:00，观测截止08:00。

时间模型采用固定源区、未来六小时内发生一次事件、时间密度与Q2触发强度成正比的情景；题面10:30用于案例评价。Q4使用文献约束的103百万m³释放体积和36百万m³下游可动体积。各问继承前问的特征、指数和同一源区，参数及完整公式见作答PDF。

| 位置 | 内容 |
|---|---|
| `run.py`、`src/mathorcup/` | 运行入口和四问模型 |
| `configs/`、`events/` | 参数、事件与计算设置 |
| `tests/` | 实现及数值检查 |
| `data/raw/unpacked/1`至`8`、`3_selected/` | 实际使用的官方附件和目标冰川瓦片 |
| `data/external/q2_forecast/` | 归档天气原始JSON及来源记录 |
| `data/processed/` | 已处理时序与空间输入 |
| `assets/fonts/` | 绘图和PDF使用的中文字体 |
| `docs/` | 当前作答PDF、赛题及审查原文、历史报告 |
| `outputs/q1/`至`q4/`、`outputs/joint/` | 各问结果、图件和过程数组 |
| `outputs/redraw_data/` | 全部图件的数值输入、元数据与重绘材料 |
| `outputs/quality/` | 检查记录及输入、代码哈希清单 |

全部图件及其用途见[图件索引](outputs/FIGURE_INDEX.csv)，可直接打开[结果图廊](outputs/results_gallery.html)。

## 检查与重新打包

```bash
python -m pytest -q tests
python scripts/update_readme.py
python scripts/build_answer_pdf.py
python scripts/build_release.py
```

完整复现ZIP与代码图件ZIP生成到仓库上一级的`output/`目录。单独重算使用`--stage q1`至`--stage q4`；按照题目顺序更新前置结果。`--fetch-q2-external`更新第二问归档天气，`--refit-q4`重新标定传播参数。

数据来源：[国家青藏高原科学数据中心](https://www.tpdc.ac.cn/)及官方赛题八类附件；归档ECMWF天气由Open-Meteo提供；形变资料来自ICIMOD。Q4体积与事件传播参考：[Journal of Mountain Science, DOI:10.1007/s11629-026-1482-y](https://jms.imde.ac.cn/en/article/doi/10.1007/s11629-026-1482-y)。逐文件来源与用途记录在`outputs/quality/reproducibility_manifest.json`。
