"""矩阵求逆、t/F 分布与区间计算。

这里只依赖标准库，保证镜像内不需要下载数值库。rollout 分析的加权窗口
估计在 rollout.py 内组装，本模块只提供通用的数值与小样本分布工具。
"""

import math


class SingularDesign(ValueError):
    """设计矩阵秩亏，无法给出可解释的估计。"""


def invert(matrix):
    """高斯－约当消元求逆，主元过小视为秩亏。"""
    size = len(matrix)
    if size == 0:
        raise SingularDesign('空矩阵不可求逆')
    scale = max(max(abs(value) for value in row) for row in matrix)
    tolerance = 1e-10 * max(1.0, scale)
    augmented = [list(row) + [1.0 if i == j else 0.0 for j in range(size)]
                 for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= tolerance:
            raise SingularDesign('设计矩阵秩亏')
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        row = augmented[column]
        width = 2 * size
        for index in range(column, width):
            row[index] /= divisor
        for other in range(size):
            if other == column:
                continue
            target = augmented[other]
            factor = target[column]
            if factor == 0.0:
                continue
            for index in range(column, width):
                target[index] -= factor * row[index]
    return [row[size:] for row in augmented]


def _betacf(a, b, x):
    max_iterations = 300
    epsilon = 3e-14
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    result = d
    for step in range(1, max_iterations + 1):
        twice = 2 * step
        numerator = step * (b - step) * x / ((qam + twice) * (a + twice))
        d = 1.0 + numerator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        result *= d * c
        numerator = -(a + step) * (qab + step) * x / ((a + twice) * (qap + twice))
        d = 1.0 + numerator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < epsilon:
            break
    return result


def _betai(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_beta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                + a * math.log(x) + b * math.log1p(-x))
    factor = math.exp(log_beta)
    if x < (a + 1.0) / (a + b + 2.0):
        return factor * _betacf(a, b, x) / a
    return 1.0 - factor * _betacf(b, a, 1.0 - x) / b


def _gammaincc(a, x):
    """正则化上不完全伽马函数 Q(a, x)。"""
    if x <= 0.0:
        return 1.0
    if x < a + 1.0:
        term = 1.0 / a
        total = term
        ap = a
        for _ in range(500):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 3e-14:
                break
        return 1.0 - total * math.exp(-x + a * math.log(x) - math.lgamma(a))
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    result = d
    for step in range(1, 500):
        an = -step * (step - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < 3e-14:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * result


def t_two_sided_p(statistic, df):
    if statistic is None or df is None or df <= 0:
        return None
    if not math.isfinite(statistic):
        return None
    value = df / (df + statistic * statistic)
    return _betai(df / 2.0, 0.5, value)


def t_critical(df, confidence):
    if df is None or df <= 0:
        return None
    target = 1.0 - (1.0 - confidence) / 2.0
    low, high = 0.0, 1.0
    while _t_cdf(high, df) < target:
        high *= 2.0
        if high > 1e6:
            return None
    for _ in range(200):
        middle = 0.5 * (low + high)
        if _t_cdf(middle, df) < target:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def _t_cdf(value, df):
    if value == 0.0:
        return 0.5
    x = df / (df + value * value)
    tail = 0.5 * _betai(df / 2.0, 0.5, x)
    return 1.0 - tail if value > 0 else tail


def chi2_sf(statistic, df):
    if df <= 0 or statistic is None or not math.isfinite(statistic):
        return None
    return _gammaincc(df / 2.0, statistic / 2.0)


def f_sf(statistic, df1, df2):
    """F 分布上尾概率：I_{d1 f/(d1 f + d2)}(d1/2, d2/2) 的补。"""
    if df1 <= 0 or df2 <= 0 or statistic is None or not math.isfinite(statistic) \
            or statistic < 0:
        return None
    x = df1 * statistic / (df1 * statistic + df2)
    return 1.0 - _betai(df1 / 2.0, df2 / 2.0, x)


def interval(estimate, variance, df, confidence):
    """由方差与聚类自由度给出 t 区间。"""
    critical = t_critical(df, confidence)
    if estimate is None or critical is None:
        return {'std_error': None, 'ci_low': None, 'ci_high': None,
                't_statistic': None, 'p_value': None}
    if variance is None:
        return {'std_error': None, 'ci_low': None, 'ci_high': None,
                't_statistic': None, 'p_value': None}
    if variance < 0.0:
        if variance > -1e-8:
            variance = 0.0
        else:
            return {'std_error': None, 'ci_low': None, 'ci_high': None,
                    't_statistic': None, 'p_value': None}
    standard_error = math.sqrt(variance)
    if standard_error == 0.0:
        return {'std_error': 0.0, 'ci_low': estimate, 'ci_high': estimate,
                't_statistic': None, 'p_value': None}
    statistic = estimate / standard_error
    return {'std_error': standard_error,
            'ci_low': estimate - critical * standard_error,
            'ci_high': estimate + critical * standard_error,
            't_statistic': statistic,
            'p_value': t_two_sided_p(statistic, df)}
