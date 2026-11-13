# staggered-rollout

园区能源分析以站点周记录保存能耗、产出和改造周，按区域核对单位产出能耗。

## 使用

运行需要 Python 3.11 或更新版本，没有第三方运行依赖。数据在本地处理，不请求外部接口。

```bash
python3 -m analysis_app summary
python3 -m analysis_app summary --input data/observations.csv --output output/summary.json

python3 -m analysis_app rollout
python3 -m analysis_app rollout --targets consumption,intensity --covariate temperature
python3 -m analysis_app rollout --weighting baseline_consumption --confidence 0.9
```

标准输出是完整 JSON，输出文件使用临时文件原子替换；失败退出码为 2，不覆盖既有结果。空数据集保留零计数与空分组，缺失数值不进入数值分母。数据指纹由按 observation_id 排序后的解析记录计算。

## rollout 分析

站点不是同一周上线，先上线的站点规模也更大，所以不能用简单的上线前后平均数，也不能把已改造站点继续当未处理对照。`rollout` 命令按“采用批次 × 事件周”组织比较：

- 每个批次用采用前最后一个观察周作基期；每个事件周与基期做一次窗口化双重差分。
- 对照只取在基期周与比较周都观察到、且截至比较周尚未采用或始终未采用的站点；已采用站点不会在后续窗口里继续充当对照。
- 没有合格对照或没有观测的批次周保留为缺口，并给出原因码，不做插补。
- 估计目标可选能耗总量与单位产出能耗。产出为零的单位产出窗口标为 `undefined_outcome`，不换用其他分母。
- 批次效应、动态效应与自然周汇总都按显式权重计算（`site_count` 或 `baseline_consumption`），输出里列出各参与方的归一化权重。
- 协变量只用同一比较窗口内的变化做调整；窗口内无站点间差异时标为 `not_identified`，不会声称“已控制温度”。
- 区间与检验以站点为聚类；预趋势诊断按批次给出采用前斜率与跨批次联合检验，输出同时说明“未拒绝不等于识别假设成立”。

输出中的 `panel.site_inventory` 与每个窗口的 `roster` 构成从原始站点到比较组的可追踪清单，能沿站点、角色与原因解释某一周为什么没有总体值。分析只用于复盘，不作为节能承诺。

## 验证与容器

```bash
python3 -m unittest discover -s tests -v
python3 -m compileall -q analysis_app
docker build -t staggered-rollout .
docker run --rm staggered-rollout
```

容器内也可执行 `python -m unittest discover -s tests -v`。镜像构建阶段取得 Python 镜像，产品启动、测试和演示均无需下载。小规模样例的完整单元测试和 summary 可以连续执行，不需要等待定时任务。

## 代码与数据

- `analysis_app/io.py`：有类型的 CSV 读取、完整性校验及内容指纹。
- `analysis_app/schema.py`：字段类型、主体与聚合维度。
- `analysis_app/domain.py`：当前业务核算函数。
- `analysis_app/summary.py`：确定性的分组与数值口径。
- `analysis_app/rollout.py`：批次×事件周比较、权重汇总、预趋势诊断与站点台账。
- `analysis_app/stats.py`：最小二乘、聚类稳健协方差与分布函数（仅标准库）。
- `analysis_app/cli.py`：本地命令及原子输出。
- `docs/data.md`：数据契约与 rollout 字段约定。

演示数据全部为本地合成，不对应现实机构、设备或人员；不能用来推断真实业务效果。
