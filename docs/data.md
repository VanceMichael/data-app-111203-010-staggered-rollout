# 数据契约

园区能源分析以站点周记录保存能耗、产出和改造周，按区域核对单位产出能耗。

主测量 `consumption` 的单位为千瓦时，主体字段是 `site`，核算维度是 `region`。`observation_id` 是文件内唯一观测编号，主体允许重复观测。数值只接受有限值；空单元格解析为 None，不能直接当成零。额外列不进入当前核算。

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `observation_id` | str | 必填 |
| `site` | str | 必填 |
| `region` | str | 必填 |
| `week` | int | 必填 |
| `consumption` | float | 必填 |
| `output` | float | 必填 |
| `adopt_week` | int | 必填 |
| `temperature` | float | 必填 |

时间字段 day、week、tick、hour、minute、cycle、start_day、reported_day、known_at、visible_at、opened、observed 均以各字段名所述的相对时间单位表示整数；不同字段并不必然代表同一种时间。缺少的登记码与浓度等信息保持为空。示例只表达产品数据结构，不能代替实际调查资料。

## 业务字段语义

consumption 为千瓦时，output 为产出件数，temperature 为摄氏度。adopt_week=-1 表示从未改造，否则为首次采用周；site 在 week 观察，region 是当前数据所属区域。

## rollout 分析的字段约定

`python -m analysis_app rollout` 把同一份面板按“采用批次 × 事件周”重新组织，不改变原始记录。

- **采用批次（cohort）**：`adopt_week` 的取值。`adopt_week=-1` 的站点不构成批次，只作为对照候选。同一站点在所有行的 `adopt_week` 必须一致，否则拒绝。
- **事件周（event）**：`week - adopt_week`。事件 `-1` 是参照周，不报告效应；每个批次以采用前最后一个观察周为基期。
- **比较窗口**：批次 `g`、事件周 `e` 对应基期周 `b` 与比较周 `w = g + e`。对照只取在 `b` 与 `w` 都能观察到、且截至 `w` 尚未采用（`adopt_week > w`）或始终未采用（`adopt_week = -1`）的站点。已采用站点不会在后续窗口里继续被当作对照。
- **估计目标**：`consumption`（千瓦时）与 `intensity`（千瓦时/件）。`intensity` 在 `output = 0` 时没有定义，对应窗口标为 `undefined_outcome`，不换用其他分母、也不退回总量。
- **权重口径**：`--weighting site_count`（默认，参与站点等权）或 `baseline_consumption`（按站点基期平均能耗加权）。两种口径都通过 WLS 进入窗口估计，跨窗口汇总使用同一口径并明确归一化。
- **协变量**：`temperature` 只用该比较窗口基期周与比较周的站点内变化（`Δtemperature`）进入窗口模型；窗口外数据不参与。若窗口内所有站点变化相同，该项被窗口截距吸收，标记 `not_identified`，不会声称“已控制温度”。
- **区间与检验**：站点是聚类。每个窗口先单独估计，再用站点层面影响函数合成协方差；预趋势联合检验另用 wild cluster bootstrap-t。聚类数少时 p 值粒度粗，输出会标注。

输出中的 `panel.site_inventory` 是站点级台账（采用周、观察周、基期周、产出为零的周），`targets.<target>.windows[*].roster` 是每个比较窗口的逐站点角色与排除原因，`reason_codes` 给出缺口原因的对照表。任何“某周没有总体值”都能沿这两张表追到具体站点与原因。

分析只用于回顾性复盘，不作为节能承诺。
