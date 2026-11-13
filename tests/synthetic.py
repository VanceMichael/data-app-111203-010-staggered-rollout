"""合成站点面板生成器，用于 rollout 分析的边界测试。

所有场景都是本地构造的，不代表真实机构、设备或人员。
"""

import random


def build_sites(adopt_weeks, **options):
    """由采用周列表构造站点规格；站点名按顺序为 S00、S01……

    可选参数（可为常数、按站点序号的列表或 callable(index, name)）：
      baseline: 基期能耗水平
      slope: 每周自然趋势
      output: 周产出件数
      effect: 采用后每事件的效应
      pre_effect: 采用前每事件的提前反应（用于预趋势诊断）
      temperature_offset: 逐站温度偏移（叠加在共同的周温度上）
      temperature_slope: 逐站每周温度趋势（让窗口内温度变化随站点不同）
    """
    sites = []
    for index, adopt in enumerate(adopt_weeks):
        name = f'S{index:02d}'
        sites.append({
            'site': name,
            'adopt_week': adopt,
            'region': f'R{index % 3}',
            'baseline': _value(options.get('baseline', 100.0), index, name),
            'slope': _value(options.get('slope', 0.0), index, name),
            'output': _value(options.get('output', 10.0), index, name),
            'temperature_offset': _value(options.get('temperature_offset', 0.0),
                                         index, name),
            'temperature_slope': _value(options.get('temperature_slope', 0.0),
                                        index, name),
            'effect': options.get('effect'),
            'pre_effect': options.get('pre_effect'),
        })
    return sites


def _value(spec, index, name):
    if callable(spec):
        return spec(index, name)
    if isinstance(spec, (list, tuple)):
        return spec[index]
    return spec


def make_rows(sites, weeks, temperature=None, noise=0.0, seed=7):
    """展开为记录。

    temperature 为 week -> 共同温度 的映射，缺省 20.0；每个站点再叠加自己的
    temperature_offset，形成站点间差异（协变量识别需要这种差异）。
    """
    rng = random.Random(seed)
    temperature = temperature or {}
    rows = []
    for spec in sites:
        adopt = spec['adopt_week']
        for week in weeks:
            event = week - adopt if adopt is not None and adopt >= 0 else None
            level = spec['baseline'] + spec['slope'] * week
            if event is not None and event >= 0:
                level += _event_effect(spec['effect'], spec['site'], event)
            if event is not None and event < 0:
                level += _event_effect(spec['pre_effect'], spec['site'], event)
            if noise:
                level += rng.gauss(0.0, noise)
            rows.append({
                'observation_id': f"{spec['site']}-w{week:02d}",
                'site': spec['site'],
                'region': spec['region'],
                'week': week,
                'consumption': level,
                'output': spec['output'],
                'adopt_week': -1 if adopt is None else adopt,
                'temperature': (temperature.get(week, 20.0)
                                + spec['temperature_offset']
                                + spec['temperature_slope'] * week),
            })
    return rows


def _event_effect(spec, site, event):
    if spec is None:
        return 0.0
    if callable(spec):
        return spec(site, event)
    return spec
