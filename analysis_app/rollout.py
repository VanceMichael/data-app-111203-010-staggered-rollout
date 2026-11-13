"""交叠上线（staggered rollout）的批次×事件周比较分析。

估计目标：能耗总量与单位产出能耗。每个批次在每个相对事件周上，用“当时尚未
采用或始终未采用”的站点作为对照，做以采用前一周为基期的双重差分。比较窗口
之外的行不进入该比较；已采用站点在任何比较窗口中都不会被当作对照。

产出包括：批次×事件周效应、按明确权重汇总的动态效应、按自然周汇总、预趋势
诊断、协变量调整状态，以及从原始站点到比较组的可追踪清单。本模块只做回顾性
核算，不构成任何节能承诺。
"""

import math

from .io import fingerprint
from .stats import SingularDesign, chi2_sf, f_sf, invert, interval

TARGET_UNITS = {'consumption': 'kWh', 'intensity': 'kWh/件'}
TARGETS = tuple(TARGET_UNITS)
WEIGHTINGS = {
    'site_count': '参与站点等权',
    'baseline_consumption': '参与站点基期能耗加权',
}
DEFAULT_WEIGHTING = 'site_count'
REFERENCE_EVENT = -1
DISCLAIMER = ('回顾性核算：估计的是改造后相对同期未改造站点的变化，'
              '不是节能承诺，也不能外推为未来保证。')
PRETREND_NOTE = (
    '预趋势检验未拒绝零假设，只表示在当前样本与检验功效下没有发现可检测的'
    '差异化趋势；它不能证明平行趋势成立，也不能当作识别假设已经验证。'
)
COVARIATE_NOTE = (
    '协变量只用同一比较窗口（基期周与比较周）内的变化做调整；窗口外的数据'
    '不参与该窗口的调整。若同一窗口内所有站点的协变量变化相同，它会被窗口'
    '截距吸收，无法与窗口效应分开识别。'
)


class UndefinedOutcome(Exception):
    """该站点周在指定估计目标上没有可用结果值。"""

    def __init__(self, reason):
        super().__init__(reason)
        self.reason = reason


def _is_never(adopt_week):
    return adopt_week is not None and adopt_week < 0


def build_panel(records):
    """把原始记录整理为站点面板；所有顺序都显式排序，保证可复现。"""
    sites = {}
    for row in records:
        site = row['site']
        entry = sites.get(site)
        if entry is None:
            entry = sites[site] = {'site': site, 'regions': set(),
                                   'adopt_week': row['adopt_week'], 'weeks': {}}
        elif entry['adopt_week'] != row['adopt_week']:
            raise ValueError('站点采用周不一致：' + site)
        entry['regions'].add(row['region'])
        week = row['week']
        if week in entry['weeks']:
            raise ValueError('站点周重复：' + site + '/' + str(week))
        entry['weeks'][week] = row
    weeks = sorted({row['week'] for row in records})
    cohorts = sorted({row['adopt_week'] for row in records
                      if not _is_never(row['adopt_week'])})
    never = sorted(site for site, entry in sites.items()
                   if _is_never(entry['adopt_week']))
    inventory = []
    for site in sorted(sites):
        entry = sites[site]
        observed = sorted(entry['weeks'])
        baseline = _baseline_weeks(entry)
        energy = [entry['weeks'][week]['consumption'] for week in baseline]
        inventory.append({
            'site': site,
            'region': sorted(entry['regions'])[0],
            'adopt_week': entry['adopt_week'],
            'cohort': None if _is_never(entry['adopt_week']) else entry['adopt_week'],
            'weeks_observed': observed,
            'baseline_weeks': baseline,
            'output_zero_weeks': [week for week in observed
                                  if entry['weeks'][week]['output'] == 0],
            'baseline_consumption_mean': math.fsum(energy) / len(energy) if energy else None,
        })
    return {'sites': sites, 'weeks': weeks, 'cohorts': cohorts, 'never': never,
            'inventory': inventory}


def _baseline_weeks(entry):
    observed = sorted(entry['weeks'])
    if _is_never(entry['adopt_week']):
        return observed
    return [week for week in observed if week < entry['adopt_week']]


def _outcome(panel, target, site, week):
    row = panel['sites'][site]['weeks'].get(week)
    if row is None:
        raise UndefinedOutcome('missing_observation')
    if target == 'consumption':
        return row['consumption']
    output = row['output']
    if output == 0:
        raise UndefinedOutcome('output_zero')
    if output < 0:
        raise UndefinedOutcome('output_negative')
    return row['consumption'] / output


def _base_week(panel, cohort):
    earlier = [week for week in panel['weeks'] if week < cohort]
    return earlier[-1] if earlier else None


def site_weights(panel, weighting):
    """参与权重；返回 site -> 权重或 None（不可参与）。显式、可审计。"""
    if weighting == 'site_count':
        return {site: 1.0 for site in panel['sites']}
    weights = {}
    for site, entry in panel['sites'].items():
        weeks = _baseline_weeks(entry)
        if not weeks:
            weights[site] = None
            continue
        energy = math.fsum(entry['weeks'][week]['consumption'] for week in weeks)
        mean = energy / len(weeks)
        weights[site] = mean if mean > 0 else None
    return weights


def _window_roster(panel, target, cohort, base, week, weights):
    """逐站点列出该比较窗口的角色与理由，供估计和台账共用。"""
    roster, treated, controls = [], [], []
    excluded = {'treated_at_or_before_base': [],
                'treated_at_or_before_comparison_week': [],
                'missing_observation': [], 'undefined_outcome': [],
                'nonpositive_weight': []}
    for site in sorted(panel['sites']):
        entry = panel['sites'][site]
        adopt = entry['adopt_week']
        if adopt == cohort:
            role = 'treated'
        elif not _is_never(adopt) and adopt <= base:
            roster.append({'site': site, 'role': 'excluded',
                           'reason': 'treated_at_or_before_base'})
            excluded['treated_at_or_before_base'].append(site)
            continue
        elif not _is_never(adopt) and adopt <= week:
            roster.append({'site': site, 'role': 'excluded',
                           'reason': 'treated_at_or_before_comparison_week'})
            excluded['treated_at_or_before_comparison_week'].append(site)
            continue
        else:
            role = 'control'
        if base not in entry['weeks'] or week not in entry['weeks']:
            roster.append({'site': site, 'role': 'excluded',
                           'reason': 'missing_observation'})
            excluded['missing_observation'].append(site)
            continue
        if not weights.get(site):
            roster.append({'site': site, 'role': 'excluded',
                           'reason': 'nonpositive_weight'})
            excluded['nonpositive_weight'].append(site)
            continue
        try:
            _outcome(panel, target, site, base)
            _outcome(panel, target, site, week)
        except UndefinedOutcome as error:
            roster.append({'site': site, 'role': 'excluded', 'reason': error.reason})
            excluded['undefined_outcome'].append(site)
            continue
        roster.append({'site': site, 'role': role, 'reason': None})
        (treated if role == 'treated' else controls).append(site)
    if treated and controls:
        status, reason = 'estimated', None
    else:
        status = 'not_estimable'
        if not treated and excluded['undefined_outcome']:
            reason = 'undefined_outcome'
        elif not controls and excluded['nonpositive_weight']:
            reason = 'nonpositive_weight'
        elif not treated:
            reason = 'no_treated_observation'
        else:
            reason = 'no_clean_control'
    return {'status': status, 'reason': reason, 'roster': roster,
            'treated_sites': treated, 'control_sites': controls,
            'excluded': excluded}


def enumerate_windows(panel, target, weights=None):
    """枚举全部（批次，观察周）；不支持的窗口保留为显式缺口。"""
    weights = weights or site_weights(panel, 'site_count')
    windows = []
    for cohort in panel['cohorts']:
        base = _base_week(panel, cohort)
        for week in panel['weeks']:
            event = week - cohort
            if base is None:
                window = _blank_window(cohort, event, None, week,
                                       'not_estimable', 'no_pre_period')
            elif event == REFERENCE_EVENT:
                window = _blank_window(cohort, event, base, week, 'reference', None)
            else:
                window = _window_roster(panel, target, cohort, base, week, weights)
                window.update({'cohort': cohort, 'event': event,
                               'base_week': base, 'week': week})
            windows.append(window)
    return windows


def _blank_window(cohort, event, base, week, status, reason):
    return {'cohort': cohort, 'event': event, 'base_week': base, 'week': week,
            'status': status, 'reason': reason, 'roster': [],
            'treated_sites': [], 'control_sites': [],
            'excluded': {'treated_at_or_before_base': [],
                         'treated_at_or_before_comparison_week': [],
                         'missing_observation': [], 'undefined_outcome': [],
                         'nonpositive_weight': []}}


def _temperature_changes(panel, base, week, sites):
    changes = {}
    for site in sites:
        before = panel['sites'][site]['weeks'].get(base)
        after = panel['sites'][site]['weeks'].get(week)
        if before is not None and after is not None:
            changes[site] = after['temperature'] - before['temperature']
    return changes


def _window_change(panel, target, site, window):
    return (_outcome(panel, target, site, window['week'])
            - _outcome(panel, target, site, window['base_week']))


def _weighted_mean(pairs):
    total = math.fsum(weight for weight, _ in pairs)
    if total <= 0:
        return None
    return math.fsum(weight * value for weight, value in pairs) / total


def _weight_share(weights, sites):
    total = math.fsum(weights[site] for site in sites)
    if total <= 0:
        return None
    return max(weights[site] for site in sites) / total


def _window_estimate(panel, target, window, covariate, weights):
    """单个比较窗口的加权双重差分，附带站点层面的影响函数。

    模型为 Δ_i = alpha + tau * treated_i + beta * Δtemperature_i + e_i，
    用站点权重做 WLS，只用该窗口基期周与比较周的数据。tau 就是该窗口的
    处理效应；影响函数 psi_i = (X'WX)^{-1} x_i w_i e_i 的 tau 分量用于在
    站点聚类下合成窗口间的协方差（共用站点会带来跨窗口相关）。
    """
    term = (window['cohort'], window['event'])
    sites = window['treated_sites'] + window['control_sites']
    treated = set(window['treated_sites'])
    changes = _temperature_changes(panel, window['base_week'], window['week'], sites)
    values = list(changes.values())
    identify_covariate = bool(
        covariate != 'none' and values and
        max(values) - min(values) > 1e-12 * max(1.0, max(map(abs, values))))
    columns = 2 + (1 if identify_covariate else 0)
    if columns == 2:
        design = [[1.0, 1.0 if site in treated else 0.0] for site in sites]
    else:
        design = [[1.0, 1.0 if site in treated else 0.0, changes.get(site, 0.0)]
                  for site in sites]
    outcomes = [_window_change(panel, target, site, window) for site in sites]
    site_weights = [weights[site] for site in sites]

    cross = [[0.0] * columns for _ in range(columns)]
    cross_y = [0.0] * columns
    for row, value, weight in zip(design, outcomes, site_weights):
        for first in range(columns):
            cross_y[first] += weight * row[first] * value
            for second in range(columns):
                cross[first][second] += weight * row[first] * row[second]
    inverse = invert(cross)
    beta = [math.fsum(inverse[i][j] * cross_y[j] for j in range(columns))
            for i in range(columns)]
    residuals = [value - math.fsum(coef * x for coef, x in zip(beta, row))
                 for row, value in zip(design, outcomes)]
    # 组内残差只围绕组均值波动，直接用会低估方差；按 n_g/(n_g-1) 校正。
    # 无协变量时这恰好还原两样本均值差的方差 s_T²/n_T + s_C²/n_C。
    scale = {}
    for group in (True, False):
        members = [site for site in sites if (site in treated) is group]
        scale[group] = (math.sqrt(len(members) / (len(members) - 1.0))
                        if len(members) > 1 else 1.0)
    influence = {}
    for row, weight, residual, site in zip(design, site_weights, residuals, sites):
        factor = weight * residual * scale[site in treated]
        influence[site] = [math.fsum(inverse[i][j] * row[j] for j in range(columns))
                           * factor for i in range(columns)]

    treated_changes = [(weights[site], value) for site, value in
                       zip(sites, outcomes) if site in treated]
    control_changes = [(weights[site], value) for site, value in
                       zip(sites, outcomes) if site not in treated]
    return {
        'term': term,
        'estimate': beta[1],
        'tau_influence': {site: entry[1] for site, entry in influence.items()},
        'covariate_slope': beta[2] if identify_covariate else None,
        'identify_covariate': identify_covariate,
        'window_rows': len(sites),
        'n_treated': len(treated_changes),
        'n_control': len(control_changes),
        'treated_mean_change': _weighted_mean(treated_changes),
        'control_mean_change': _weighted_mean(control_changes),
        'raw_difference_in_means': (_weighted_mean(treated_changes)
                                    - _weighted_mean(control_changes)),
        'treated_weight_share': _weight_share(weights, window['treated_sites']),
        'sites': sites,
    }


def _combine_window_influence(estimates, influence_key):
    """把各窗口的站点影响函数合成聚类稳健协方差。

    每个站点是一个聚类；窗口 w 的方差用 G_w/(G_w-1) 校正，窗口 w 与 v 的
    协方差用 sqrt(f_w f_v)。共用站点的窗口会自动体现相关性。
    """
    terms = [entry['term'] for entry in estimates]
    index = {term: position for position, term in enumerate(terms)}
    all_sites = sorted({site for entry in estimates
                        for site in entry[influence_key]})
    factors = {entry['term']: (len(entry[influence_key]) /
                               (len(entry[influence_key]) - 1.0)
                               if len(entry[influence_key]) > 1 else 1.0)
               for entry in estimates}
    covariance = [[0.0] * len(terms) for _ in terms]
    for first, first_entry in enumerate(estimates):
        for second, second_entry in enumerate(estimates):
            if second < first:
                continue
            factor = math.sqrt(factors[first_entry['term']] *
                               factors[second_entry['term']])
            total = 0.0
            for site in all_sites:
                total += (first_entry[influence_key].get(site, 0.0)
                          * second_entry[influence_key].get(site, 0.0))
            covariance[first][second] = total * factor
            covariance[second][first] = total * factor
    df = max(len(all_sites) - 1, 0)
    return {'terms': terms, 'term_index': index, 'covariance': covariance,
            'n_clusters': len(all_sites), 'df': df}


def _pool(entries, covariance, df, confidence, weight_map, weight_basis):
    """按显式权重汇总系数；方差来自同一聚类稳健协方差。"""
    if not entries:
        return {'status': 'no_effects', 'estimate': None, 'weights': {}}
    total = math.fsum(weight for _, weight, _ in entries)
    if total <= 0:
        return {'status': 'no_effects', 'estimate': None, 'weights': {},
                'reason': '参与权重为零'}
    normalized = [(index, weight / total, value) for index, weight, value in entries]
    estimate = math.fsum(weight * value for _, weight, value in normalized)
    variance = 0.0
    for first_index, first_weight, _ in normalized:
        for second_index, second_weight, _ in normalized:
            variance += first_weight * second_weight * covariance[first_index][second_index]
    summary = interval(estimate, variance, df, confidence)
    summary.update({
        'status': 'estimated', 'estimate': estimate,
        # 键统一转成字符串，使库结果与 JSON 输出逐位一致。
        'weights': {str(key): weight / total for key, weight in
                    sorted(weight_map.items(), key=lambda item: str(item[0]))},
        'weight_basis': weight_basis,
    })
    return summary


def _functional_covariance(functionals, model):
    index = model['term_index']
    covariance = model['covariance']
    size = len(functionals)
    matrix = [[0.0] * size for _ in range(size)]
    for first in range(size):
        for second in range(first, size):
            value = math.fsum(
                first_weight * second_weight
                * covariance[index[first_term]][index[second_term]]
                for first_term, first_weight in functionals[first].items()
                for second_term, second_weight in functionals[second].items())
            matrix[first][second] = value
            matrix[second][first] = value
    return matrix


def _reference_df(estimates, terms_used):
    """识别变异来自处理组站点，因此参考分布用处理组聚类数。"""
    by_term = {entry['term']: entry for entry in estimates}
    counts = [by_term[term]['n_treated'] for term in terms_used if term in by_term]
    return max((min(counts) - 1) if counts else 0, 0)


def _quadratic(vector, inverse):
    return math.fsum(vector[i] * math.fsum(inverse[i][j] * vector[j]
                                           for j in range(len(vector)))
                     for i in range(len(vector)))


def joint_test(model, estimates, functionals):
    """线性函数向量=0 的联合检验。

    统计量用站点聚类协方差；参考分布用 F(q, df) 而不是卡方，因为处理组
    站点非常少，卡方近似会明显过度拒绝。df 取各批次处理组站点数的最小值
    减一：识别变异来自处理组，控制组的重复使用不提供同等的独立信息。

    这是刻意的保守处理；即使不拒绝，也只说明当前数据没有检测到差异，
    不说明平行趋势成立。
    """
    if not functionals:
        return {'status': 'not_estimable', 'p_value': None, 'restrictions': 0,
                'reason': '没有可检验的限制'}
    index = model['term_index']
    vector = [math.fsum(weight * model['coefficients'][index[term]]
                        for term, weight in functional.items())
              for functional in functionals]
    covariance = _functional_covariance(functionals, model)
    try:
        inverse = invert(covariance)
    except SingularDesign:
        return {'status': 'not_estimable', 'p_value': None,
                'restrictions': len(functionals),
                'reason': '线性函数的协方差秩亏，无法做联合检验'}
    statistic = _quadratic(vector, inverse)
    restrictions = len(functionals)
    terms_used = {term for functional in functionals for term in functional}
    df = _reference_df(estimates, terms_used)
    result = {
        'status': 'estimated', 'restrictions': restrictions,
        'estimates': vector, 'statistic': statistic,
        'chi2_p_value': chi2_sf(statistic, restrictions),
        'df_residual': df,
        'p_value_basis': 'F（按处理组聚类数取自由度的保守参考分布）',
    }
    if df > 0:
        result['f_statistic'] = statistic / restrictions
        result['f_df'] = [restrictions, df]
        result['f_p_value'] = f_sf(statistic / restrictions, restrictions, df)
        result['p_value'] = result['f_p_value']
        result['trustworthy'] = restrictions < df
        if not result['trustworthy']:
            result['reason'] = ('限制个数不少于处理组聚类数对应的自由度，'
                                '联合检验的 p 值只能当作粗略参考，不能据此'
                                '宣称预趋势成立或未成立')
    else:
        result['p_value'] = result['chi2_p_value']
        result['trustworthy'] = False
        result['reason'] = '没有足够的处理组聚类估计参考自由度'
    return result


def _restriction_functionals(pre_terms, term_index):
    """按批次构造事件时间线性斜率，作为窗口系数的线性函数。"""
    functionals = []
    cohorts = []
    details = []
    for cohort in sorted({term[0] for term in pre_terms}):
        subset = [term for term in pre_terms if term[0] == cohort]
        events = [term[1] for term in subset]
        centered = [event - math.fsum(events) / len(events) for event in events]
        denominator = math.fsum(value * value for value in centered)
        if len(subset) < 2 or denominator == 0.0:
            details.append({'cohort': cohort, 'status': 'not_estimable',
                            'reason': '采用前事件周少于两个，无法估计预趋势斜率'})
            continue
        functionals.append({term: value / denominator
                            for term, value in zip(subset, centered)})
        cohorts.append(cohort)
        details.append({'cohort': cohort, 'events': events,
                        'status': 'estimable'})
    return functionals, cohorts, details


def _pretrend(model, estimates, pre_terms):
    """预趋势诊断：批次内事件时间线性斜率 + 跨批次斜率联合检验。

    斜率把每个批次多个采用前系数压成一个自由度，检验处理组与对照组在采用
    前是否有差异化趋势。协方差由站点聚类影响函数合成；参考分布用处理组
    聚类数对应的 F 分布（见 joint_test），避免处理组站点很少时卡方近似过度
    拒绝。高维“全部采用前系数为零”检验也保留，但限制个数多时明确标注只能
    作参考。
    """
    term_index = model['term_index']
    pre_terms = [term for term in pre_terms if term in term_index]
    functionals, cohorts, details = _restriction_functionals(pre_terms, term_index)
    by_cohort = {entry['cohort']: entry for entry in details}

    per_cohort = []
    for functional, cohort in zip(functionals, cohorts):
        covariance = _functional_covariance([functional], model)[0][0]
        estimate = math.fsum(weight * model['coefficients'][term_index[term]]
                             for term, weight in functional.items())
        summary = interval(estimate, covariance, model['df'], 0.95)
        summary.update({
            'cohort': cohort, 'slope': estimate,
            'events': by_cohort[cohort]['events'],
            'status': 'estimated', 'note': PRETREND_NOTE,
            'description': '采用前事件周线性斜率是否为零'})
        per_cohort.append(summary)
    for entry in details:
        if entry['status'] == 'not_estimable':
            per_cohort.append({'cohort': entry['cohort'],
                               'status': 'not_estimable', 'p_value': None,
                               'reason': entry['reason'], 'note': PRETREND_NOTE})
    per_cohort.sort(key=lambda item: item['cohort'])

    joint = joint_test(model, estimates, functionals)
    joint.update({'cohorts': cohorts,
                  'description': '各批次采用前斜率同时为零（每个批次一个自由度）',
                  'note': PRETREND_NOTE})

    # 高维联合检验：限制个数通常远超处理组聚类数，只保留为参考，不参与结论。
    zero_test = joint_test(model, estimates, [{term: 1.0} for term in pre_terms])
    if zero_test.get('status') == 'estimated':
        zero_test['trustworthy'] = False
        zero_test['note'] = ('限制个数超过处理组聚类数提供的自由度，p 值不可解读；'
                             '这里列出只为说明“并非没有系数为零”，不参与预趋势'
                             '结论。' + PRETREND_NOTE)
    zero_test['description'] = ('全部采用前系数同时为零（高维联合检验，仅作参考）')

    detected = [entry['cohort'] for entry in per_cohort
                if entry.get('status') == 'estimated'
                and (entry.get('p_value') or 1.0) < 0.05]
    if not functionals:
        verdict = 'unavailable'
    elif joint.get('status') != 'estimated' or joint.get('trustworthy') is False:
        verdict = 'not_trustworthy'
    elif (joint.get('p_value') or 1.0) < 0.05:
        verdict = 'detected_differential_trend'
    else:
        verdict = 'no_detectable_differential_trend'
    if detected:
        joint['cohorts_with_detected_pretrend'] = detected
    return {'joint': joint, 'zero_test': zero_test, 'per_cohort': per_cohort,
            'note': PRETREND_NOTE, 'verdict': verdict,
            'interpretation': ('未拒绝零假设只表示在当前样本与检验功效下没有发现'
                               '可检测的差异化预趋势，不能证明平行趋势成立，'
                               '也不能当作识别假设已经验证。若拒绝零假设，则'
                               '预趋势诊断未通过，事后效应应谨慎解读。')}


def _concentration(panel, cohort):
    values = {}
    for site, entry in panel['sites'].items():
        if entry['adopt_week'] != cohort:
            continue
        weeks = _baseline_weeks(entry)
        if weeks:
            values[site] = math.fsum(entry['weeks'][week]['consumption']
                                     for week in weeks) / len(weeks)
    total = math.fsum(values[site] for site in sorted(values))
    largest = max(sorted(values), key=lambda site: values[site]) if values else None
    share = values[largest] / total if largest and total else None
    return {
        'largest_site': largest,
        'largest_site_share_of_baseline_consumption': share,
        'dominant_site_flag': bool(share is not None and share > 0.5),
        'note': '站点数权重下各站点等权；基期能耗权重下大站点会主导。'
                '此份额用于判断汇总是否被单个站点带动。',
    }


def _public_windows(windows, effects=None):
    effects = effects or {}
    public = []
    for window in windows:
        key = (window['cohort'], window['event'])
        entry = {
            'cohort': window['cohort'], 'event': window['event'],
            'week': window['week'], 'base_week': window['base_week'],
            'status': window['status'], 'reason': window['reason'],
            'roster': window['roster'],
            'counts': {
                'treated': len(window['treated_sites']),
                'control': len(window['control_sites']),
                'excluded_treated_at_or_before_base':
                    len(window['excluded']['treated_at_or_before_base']),
                'excluded_treated_at_or_before_comparison_week':
                    len(window['excluded']['treated_at_or_before_comparison_week']),
                'excluded_missing_observation':
                    len(window['excluded']['missing_observation']),
                'excluded_undefined_outcome':
                    len(window['excluded']['undefined_outcome']),
                'excluded_nonpositive_weight':
                    len(window['excluded']['nonpositive_weight']),
            },
        }
        if key in effects:
            entry['estimate'] = effects[key]['estimate']
        public.append(entry)
    return public


def _target_result(panel, target, windows, covariate, confidence, weighting, weights):
    result = {
        'target': target,
        'unit': TARGET_UNITS[target],
        'weighting': weighting,
        'weighting_label': WEIGHTINGS[weighting],
        'covariate': {'requested': covariate, 'status': 'not_requested',
                      'note': COVARIATE_NOTE},
        'windows': [], 'cohort_effects': [], 'dynamic_effects': [],
        'calendar_weeks': [], 'post_average': None, 'pretrend': None,
    }
    unavailable_pretrend = {
        'joint': {'status': 'not_estimable', 'restrictions': 0, 'p_value': None},
        'per_cohort': [], 'note': PRETREND_NOTE, 'verdict': 'unavailable',
        'interpretation': '没有可估计窗口时无法做预趋势检验，不能据此判断识别假设。'}
    estimated_windows = [window for window in windows
                         if window['status'] == 'estimated']
    if not estimated_windows:
        result.update({'status': 'not_estimable',
                       'reason': '没有可估计的比较窗口',
                       'pretrend': unavailable_pretrend})
        _fill_ledgers(panel, result, windows, None, None, weighting, weights)
        result['windows'] = _public_windows(windows)
        return result

    estimates = []
    try:
        for window in estimated_windows:
            estimates.append(_window_estimate(panel, target, window, covariate, weights))
    except SingularDesign:
        result.update({'status': 'not_estimable',
                       'reason': '单个比较窗口内处理组与对照组无法区分（设计秩亏）',
                       'pretrend': unavailable_pretrend})
        _fill_ledgers(panel, result, windows, None, None, weighting, weights)
        result['windows'] = _public_windows(windows)
        return result

    _describe_covariate(result['covariate'], estimates, covariate)
    model = _combine_window_influence(estimates, 'tau_influence')
    model['coefficients'] = [entry['estimate'] for entry in estimates]
    df = model['df']
    result.update({
        'status': 'estimated',
        'n_comparison_rows': sum(entry['window_rows'] for entry in estimates),
        'n_clusters': model['n_clusters'],
        'df': df,
        'inference': {
            'method': '逐窗口加权双重差分，站点层面影响函数合成聚类稳健协方差',
            'n_clusters': model['n_clusters'],
            'df': df,
            'note': ('每个比较窗口单独估计，只有窗口内的基期周与比较周进入该'
                     '窗口；共用站点的窗口通过影响函数保留相关性。'),
        },
    })

    term_index = model['term_index']
    covariance = model['covariance']
    effects = {}
    for estimate in estimates:
        term = estimate['term']
        position = term_index[term]
        summary = interval(estimate['estimate'], covariance[position][position],
                           df, confidence)
        summary['estimate'] = estimate['estimate']
        for name in ('treated_mean_change', 'control_mean_change',
                     'raw_difference_in_means', 'treated_weight_share',
                     'window_rows'):
            summary[name] = estimate[name]
        effects[term] = summary

    pre_terms = [term for term in model['terms'] if term[1] < 0]
    result['pretrend'] = _pretrend(model, estimates, pre_terms)
    _fill_ledgers(panel, result, windows, effects, model, weighting, weights)
    result['windows'] = _public_windows(windows, effects)
    return result


def _describe_covariate(covariate_result, estimates, covariate):
    values = [entry for entry in estimates if entry['identify_covariate']]
    covariate_result['windows_total'] = len(estimates)
    covariate_result['windows_with_cross_site_variation'] = len(values)
    if values:
        covariate_result['status'] = 'identified'
        covariate_result['adjusted'] = True
        covariate_result['coefficient_scope'] = (
            '每个比较窗口各自估计的协变量斜率；这里给出各窗口斜率的算术'
            '汇总，仅作描述，不进入处理效应。协变量只用该窗口基期周与比较周'
            '的变化。')
        sloped = [entry['covariate_slope'] for entry in values]
        covariate_result['coefficient'] = math.fsum(sloped) / len(sloped)
        covariate_result['coefficient_range'] = [min(sloped), max(sloped)]
        return
    covariate_result['adjusted'] = False
    if covariate == 'none':
        covariate_result['status'] = 'not_requested'
        covariate_result['detail'] = '未请求协变量调整。'
    else:
        covariate_result['status'] = 'not_identified'
        covariate_result['detail'] = (
            '所有比较窗口内站点的协变量变化相同，已被窗口截距吸收，'
            '无法与窗口效应分开识别；因此该目标未做协变量调整，报告中不应写'
            '“已控制温度”。')


def _fill_ledgers(panel, result, windows, effects, model, weighting, weights,
                  confidence=0.95):
    """批次×事件周、动态效应、自然周与事后平均；无估计时也给出完整缺口。"""
    effects = effects or {}
    model = model or {}
    by_cohort = {}
    for window in windows:
        key = (window['cohort'], window['event'])
        record = {
            'event': window['event'], 'week': window['week'],
            'base_week': window['base_week'], 'status': window['status'],
            'reason': window['reason'], 'estimate': None, 'std_error': None,
            'ci_low': None, 'ci_high': None, 'p_value': None,
            'treated_mean_change': None, 'control_mean_change': None,
            'raw_difference_in_means': None,
            'treated_sites': len(window['treated_sites']),
            'control_sites': len(window['control_sites']),
            'reference': window['event'] == REFERENCE_EVENT,
        }
        if key in effects:
            record.update({name: effects[key].get(name) for name in
                           ('estimate', 'std_error', 'ci_low', 'ci_high', 'p_value',
                            'treated_mean_change', 'control_mean_change',
                            'raw_difference_in_means')})
        by_cohort.setdefault(window['cohort'], []).append(record)
    result['cohort_effects'] = [
        {'cohort': cohort,
         'events': sorted(by_cohort[cohort], key=lambda item: item['event']),
         'concentration': _concentration(panel, cohort)}
        for cohort in sorted(by_cohort)]

    def window_weight(window):
        return math.fsum(weights[site] for site in window['treated_sites'])

    covariance = model.get('covariance')
    df = model.get('df')
    term_index = model.get('term_index', {})

    # 动态效应：同一事件周跨批次，按参与站点权重汇总。
    events = sorted({window['event'] for window in windows
                     if window['event'] != REFERENCE_EVENT})
    for event in events:
        contributions, weight_map, missing = [], {}, []
        for window in windows:
            if window['event'] != event:
                continue
            key = (window['cohort'], event)
            if effects and key in effects:
                weight = window_weight(window)
                contributions.append((term_index[key], weight,
                                      effects[key]['estimate']))
                weight_map[window['cohort']] = weight
            elif window['status'] == 'not_estimable':
                missing.append({'cohort': window['cohort'],
                                'status': window['status'],
                                'reason': window['reason']})
        if covariance is None:
            summary = _no_effects('没有可估计的比较窗口')
        else:
            summary = _pool(contributions, covariance, df, confidence,
                            weight_map,
                            WEIGHTINGS[weighting] + '；同一事件周内归一化')
        summary['event'] = event
        summary['contributing_cohorts'] = sorted(weight_map)
        summary['missing_cohorts'] = missing
        result['dynamic_effects'].append(summary)

    # 自然周：把各批次当周的窗口按同一显式权重汇总；缺口给出原因。
    for week in panel['weeks']:
        adopted = [cohort for cohort in panel['cohorts'] if cohort <= week]
        if not adopted:
            result['calendar_weeks'].append(_calendar_gap(
                week, 'no_adopted_cohort',
                '该周没有已采用站点，没有可汇总的批次效应'))
            continue
        contributions, weight_map, missing = [], {}, []
        for cohort in adopted:
            key = (cohort, week - cohort)
            window = next((item for item in windows
                           if item['cohort'] == cohort and item['week'] == week), None)
            if effects and key in effects and window is not None:
                weight = window_weight(window)
                contributions.append((term_index[key], weight,
                                      effects[key]['estimate']))
                weight_map[cohort] = weight
            else:
                missing.append({
                    'cohort': cohort,
                    'status': window['status'] if window else 'no_window',
                    'reason': window['reason'] if window else 'no_comparison_window'})
        if covariance is None:
            summary = _calendar_gap(
                week, 'no_estimable_effect',
                '没有可估计的比较窗口，无法汇总该周；逐批次原因见 missing_cohorts')
            if any(item['reason'] == 'no_clean_control' for item in missing):
                summary['reason'] += '；已采用批次在该周没有未采用对照站点'
        else:
            summary = _pool(contributions, covariance, df, confidence,
                            weight_map, WEIGHTINGS[weighting] + '；同一自然周内'
                            '归一化。各批次沿用自己的基期周，所以这是不同基期'
                            '效应的并置平均，不是单一事件研究')
            summary['week'] = week
            if not contributions:
                summary['status'] = 'no_estimable_effect'
                summary['reason'] = '该周所有已采用批次都没有合格比较窗口或对照站点'
        summary['contributing_cohorts'] = sorted(weight_map)
        summary['missing_cohorts'] = missing
        result['calendar_weeks'].append(summary)

    # 事后平均：全部事后（批次×事件周）系数按同一显式权重汇总，仅供参考。
    contributions, weight_map = [], {}
    for window in windows:
        key = (window['cohort'], window['event'])
        if window['event'] < 0 or not effects or key not in effects:
            continue
        weight = window_weight(window)
        contributions.append((term_index[key], weight, effects[key]['estimate']))
        weight_map[f'{window["cohort"]}@{window["event"]}'] = weight
    if covariance is None:
        summary = _no_effects('没有可汇总的事后事件周效应')
        summary['note'] = '没有可汇总的事后事件周效应'
    else:
        summary = _pool(contributions, covariance, df, confidence, weight_map,
                        WEIGHTINGS[weighting] + '；全部事后窗口内归一化')
        summary['note'] = ('事后事件周的回顾性平均，不是节能承诺。'
                           if summary.get('status') == 'estimated'
                           else '没有可汇总的事后事件周效应')
    result['post_average'] = summary


def _no_effects(reason):
    return {'status': 'no_effects', 'estimate': None, 'weights': {},
            'reason': reason}


def _calendar_gap(week, status, reason):
    return {'week': week, 'status': status, 'reason': reason,
            'estimate': None, 'weights': {}, 'contributing_cohorts': [],
            'missing_cohorts': []}


def diagnose_data(panel, targets=('consumption', 'intensity'), weighting=DEFAULT_WEIGHTING):
    """数据级诊断：什么情况下哪些估计目标会没有值。"""
    zero_weeks = sorted((site, week) for site, entry in panel['sites'].items()
                        for week, row in entry['weeks'].items()
                        if row['output'] == 0)
    weights = site_weights(panel, weighting)
    diagnostics = {
        'output_zero_rows': [{'site': site, 'week': week} for site, week in zero_weeks],
        'output_zero_sites': sorted({site for site, _ in zero_weeks}),
        'nonpositive_weight_sites': sorted(site for site, weight in weights.items()
                                           if not weight),
        'fraction_never_treated': (len(panel['never']) / len(panel['sites'])
                                   if panel['sites'] else None),
        'notes': [],
    }
    if not panel['cohorts']:
        diagnostics['notes'].append('所有站点均未采用（adopt_week 全为 -1），'
                                    '没有处理批次可供估计。')
    elif panel['cohorts'][0] <= panel['weeks'][0]:
        diagnostics['notes'].append('最早批次在首个观察周或更早已采用，该批次没有'
                                    '采用前数据，相关窗口保留缺口。')
    if diagnostics['output_zero_sites']:
        diagnostics['notes'].append('存在产出为零的站点周；单位产出目标在这些窗口'
                                    '标为 undefined_outcome，不换用其他分母。')
    if not panel['never']:
        diagnostics['notes'].append('所有站点最终都会采用；较晚批次必须依赖仍未'
                                    '采用的站点，最后一个批次的后期事件周可能缺口。')
    windows = {target: enumerate_windows(panel, target, weights) for target in targets}
    diagnostics['windows_total_by_target'] = {
        target: len(windows[target]) for target in targets}
    diagnostics['windows_estimated_by_target'] = {
        target: sum(1 for window in windows[target]
                    if window['status'] == 'estimated') for target in targets}
    diagnostics['windows_not_estimable_by_target'] = {
        target: sum(1 for window in windows[target]
                    if window['status'] == 'not_estimable') for target in targets}
    return diagnostics


REASON_CODES = {
    'no_pre_period': '该批次之前没有观察周，无法定义基期',
    'no_treated_observation': '该窗口内没有任何已采用站点同时拥有基期周与比较周观测',
    'undefined_outcome': '该估计目标在基期周或比较周没有有限值（如产出为零）',
    'nonpositive_weight': '该权重口径下站点没有正权重（如基期能耗为零）',
    'no_clean_control': '该窗口内没有截至比较周仍未采用、且有观测的对照站点',
    'no_adopted_cohort': '该周尚未有站点采用，没有可汇总的批次效应',
    'no_estimable_effect': '该周所有已采用批次都没有合格比较窗口',
    'no_window': '该批次在该周没有比较窗口',
    'no_comparison_window': '该批次在该周没有比较窗口',
}


def analyze(records, targets=TARGETS, covariate='temperature',
            confidence=0.95, weighting=DEFAULT_WEIGHTING):
    """对给定记录执行 rollout 分析。"""
    if not 0.0 < confidence < 1.0:
        raise ValueError('置信水平必须在 (0, 1) 内')
    if weighting not in WEIGHTINGS:
        raise ValueError('未知权重口径：' + str(weighting))
    targets = tuple(targets)
    unknown = [target for target in targets if target not in TARGET_UNITS]
    if unknown:
        raise ValueError('未知估计目标：' + ','.join(map(str, unknown)))
    if not targets:
        raise ValueError('至少需要一个估计目标')
    panel = build_panel(records)
    if not panel['sites']:
        raise ValueError('没有可用观测')
    weights = site_weights(panel, weighting)
    result = {
        'product': 'staggered-rollout',
        'schema_version': 2,
        'analysis': 'staggered-rollout-att',
        'unit_consumption': 'kWh',
        'unit_output': '件',
        'confidence': confidence,
        'covariate_requested': covariate,
        'weighting': weighting,
        'weighting_label': WEIGHTINGS[weighting],
        'reference_event': REFERENCE_EVENT,
        'disclaimer': DISCLAIMER,
        'pretrend_note': PRETREND_NOTE,
        'input_fingerprint': fingerprint(records),
        'panel': {
            'rows': len(records),
            'sites': len(panel['sites']),
            'weeks': panel['weeks'],
            'cohorts': panel['cohorts'],
            'never_treated': panel['never'],
            'site_inventory': panel['inventory'],
            'comparison_rule': ('对照只取在基期周与比较周都能观察到、且截至各自'
                                '观察时点尚未采用或始终未采用的站点；已采用站点'
                                '在任何比较窗口中都不会被当作对照。'),
            'reference_rule': '每个批次以采用前最后一个观察周为基期；事件 -1 为'
                              '参照周，不报告效应。',
        },
        'targets': {},
        'warnings': [],
    }
    for target in targets:
        windows = enumerate_windows(panel, target, weights)
        result['targets'][target] = _target_result(
            panel, target, windows, covariate, confidence, weighting, weights)
    result['diagnostics'] = diagnose_data(panel, targets, weighting)
    result['reason_codes'] = REASON_CODES
    _warnings(result, panel, weights)
    return result


def _warnings(result, panel, weights):
    warnings = result['warnings']
    if not panel['never']:
        warnings.append('没有始终未采用的站点：较早批次之后只能依赖尚未采用的站点'
                        '作为对照，最后批次的后期事件周可能没有对照而保留缺口。')
    missing_weight = sorted(site for site, weight in weights.items() if not weight)
    if missing_weight:
        warnings.append('以下站点在所选权重口径下没有正权重，已从估计中排除：'
                        + ','.join(missing_weight) + '。')
    for target, outcome in result['targets'].items():
        if outcome.get('status') == 'not_estimable':
            warnings.append(f'{target} 目标不可估计：{outcome.get("reason")}。')
        covariate = outcome.get('covariate', {})
        if covariate.get('status') == 'not_identified':
            warnings.append(f'{target} 目标的协变量在全部比较窗口内没有站点间差异，'
                            '未做协变量调整。')
        elif covariate.get('status') == 'dropped_singular':
            warnings.append(f'{target} 目标加入协变量后设计矩阵秩亏，已退回无协变量估计。')
        weighted_by_energy = result['weighting'] == 'baseline_consumption'
        for cohort in outcome.get('cohort_effects', []):
            concentration = cohort['concentration']
            if concentration.get('dominant_site_flag') and weighted_by_energy:
                warnings.append(
                    f'{target} 目标批次 {cohort["cohort"]} 的基期能耗集中于站点 '
                    f'{concentration["largest_site"]}（份额 '
                    f'{concentration["largest_site_share_of_baseline_consumption"]:.3f}），'
                    '按基期能耗加权时汇总结果可能被单个大站点主导。')
        pretrend = outcome.get('pretrend') or {}
        verdict = pretrend.get('verdict')
        if verdict == 'detected_differential_trend':
            warnings.append(f'{target} 目标的预趋势诊断未通过：批次间存在可检测的'
                            '差异化预趋势，事后效应应谨慎解读。')
        elif verdict == 'not_trustworthy':
            warnings.append(f'{target} 目标的预趋势联合检验受聚类数限制，p 值不可靠；'
                            '不能据此宣称平行趋势成立。')
        elif verdict == 'unavailable':
            warnings.append(f'{target} 目标没有可用的预趋势联合检验。')
    result['warnings'] = warnings
