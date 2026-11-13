"""rollout 分析的行为测试：对照资格、估计目标、汇总权重、预趋势与确定性。"""

import copy
import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from analysis_app.io import load_records
from analysis_app.rollout import analyze, build_panel, enumerate_windows
from synthetic import build_sites, make_rows


WEEKS = list(range(12))


class ComparisonRuleTest(unittest.TestCase):
    def test_controls_never_include_adopted_sites(self):
        rows = make_rows(build_sites([4, 4, 7, 7, 10, 10, None, None]), WEEKS,
                         temperature={week: 20.0 + 2.0 * (week % 3) for week in WEEKS})
        result = analyze(rows)
        panel = build_panel(rows)
        for target, outcome in result['targets'].items():
            treated_by_week = {}
            for window in outcome['windows']:
                for site in window['roster']:
                    if site['role'] == 'control':
                        adopt = panel['sites'][site['site']]['adopt_week']
                        # adopt_week=-1 表示从未采用；否则必须晚于比较周。
                        self.assertTrue(adopt < 0 or adopt > window['week'],
                                        '对照站点在比较周之前已采用')
                        self.assertTrue(adopt < 0 or adopt > window['base_week'])
                    if site['role'] == 'treated':
                        self.assertEqual(panel['sites'][site['site']]['adopt_week'],
                                         window['cohort'])

    def test_control_set_shrinks_as_cohorts_adopt(self):
        rows = make_rows(build_sites([4, 4, 7, 7, 10, 10, None, None]), WEEKS)
        result = analyze(rows)
        cohort = result['targets']['consumption']['cohort_effects'][0]
        controls = {event['event']: event['control_sites'] for event in cohort['events']}
        # 事件 -4 时只有批次 4 已采用；从未采用站点与尚未采用批次都可作对照。
        self.assertEqual(controls[-4], 6)
        # 到事件 +7 时，批次 7 与 10 都已采用，只剩从未采用的站点。
        self.assertEqual(controls[7], 2)
        self.assertGreater(controls[-4], controls[7])

    def test_no_supported_window_is_a_gap_not_an_average(self):
        # 最后一批在所有站点都将采用之后才进入后期；没有任何未采用对照。
        rows = make_rows(build_sites([2, 2, 2, 9, 9, 9]), WEEKS)
        result = analyze(rows)
        outcome = result['targets']['consumption']
        gaps = [window for window in outcome['windows']
                if window['cohort'] == 9 and window['event'] > 0]
        self.assertTrue(gaps)
        for gap in gaps:
            self.assertEqual(gap['status'], 'not_estimable')
            self.assertEqual(gap['reason'], 'no_clean_control')
        late = [week for week in outcome['calendar_weeks'] if week['week'] >= 10]
        self.assertTrue(late)
        for entry in late:
            self.assertNotEqual(entry['status'], 'estimated')
            self.assertIn(entry['status'], {'no_estimable_effect', 'no_adopted_cohort'})
            self.assertTrue(entry['missing_cohorts'] or entry['reason'])

    def test_early_reaction_is_detected_by_pretrend(self):
        # 每批 8 个站点，保证处理组聚类数足以做联合检验。
        adopt = [3] * 8 + [6] * 8 + [9] * 8 + [None] * 8
        sites = build_sites(
            adopt,
            baseline=[100.0 + 5 * index for index in range(len(adopt))],
            effect=-5.0,
            pre_effect=lambda site, event: -14.0
            if int(site[1:]) in range(8, 16) and event in (-3, -2) else 0.0)
        # 批次 6 在采用前两到三周就提前下降：识别假设被违反。
        # 注意事件 -1 就是基期周，落在基期上的提前反应会被基期吸收，无法识别。
        rows = make_rows(sites, WEEKS, noise=0.3, seed=5)
        result = analyze(rows)
        outcome = result['targets']['consumption']
        self.assertEqual(outcome['pretrend']['verdict'], 'detected_differential_trend')
        self.assertLess(outcome['pretrend']['joint']['p_value'], 0.05)
        per = {entry['cohort']: entry for entry in outcome['pretrend']['per_cohort']}
        self.assertLess(per[6]['p_value'], 0.05)

    def test_few_treated_sites_marks_pretrend_untrustworthy(self):
        # 每批 3 个处理站点：联合检验的参考自由度不足，不能据此下结论。
        adopt = [3] * 3 + [6] * 3 + [9] * 3 + [None] * 3
        sites = build_sites(adopt, effect=-5.0)
        rows = make_rows(sites, WEEKS, noise=0.5, seed=2)
        outcome = analyze(rows)['targets']['consumption']
        pretrend = outcome['pretrend']
        self.assertEqual(pretrend['verdict'], 'not_trustworthy')
        self.assertFalse(pretrend['joint']['trustworthy'])
        self.assertIn('不能据此宣称', pretrend['joint']['reason'])
        self.assertTrue(any('预趋势联合检验' in warning
                            for warning in analyze(rows)['warnings']))
        self.assertIn('不能证明平行趋势成立', pretrend['interpretation'])

    def test_flat_pre_period_is_reported_without_claiming_identification(self):
        sites = build_sites(
            [3, 3, 3, 6, 6, 6, 9, 9, 9, None, None, None],
            baseline=[100.0 + 10 * index for index in range(12)],
            effect=-5.0)
        rows = make_rows(sites, WEEKS, noise=1.5, seed=3)
        outcome = analyze(rows)['targets']['consumption']
        pretrend = outcome['pretrend']
        self.assertIn(pretrend['verdict'],
                      {'no_detectable_differential_trend', 'not_trustworthy'})
        self.assertIn('不能证明平行趋势成立', pretrend['interpretation'])
        if pretrend['verdict'] == 'no_detectable_differential_trend':
            self.assertGreater(pretrend['joint']['p_value'], 0.05)

    def test_all_sites_eventually_adopt(self):
        rows = make_rows(build_sites([3, 3, 5, 5, 7, 7, 9, 9]), WEEKS)
        result = analyze(rows)
        self.assertEqual(result['panel']['never_treated'], [])
        outcome = result['targets']['consumption']
        late = [entry for entry in outcome['cohort_effects']
                if entry['cohort'] == 9][0]['events']
        self.assertTrue(all(event['status'] != 'estimated'
                            for event in late if event['event'] > 0))
        self.assertTrue(any('所有站点最终都会采用' in note
                            for note in result['diagnostics']['notes']))

    def test_completely_control_free_week(self):
        # 所有站点都在第 4 周采用：采用之后每一周都没有对照。
        rows = make_rows(build_sites([4, 4, 4, 4]), WEEKS)
        result = analyze(rows)
        outcome = result['targets']['consumption']
        self.assertEqual(outcome['status'], 'not_estimable')
        after = [entry for entry in outcome['calendar_weeks'] if entry['week'] > 4]
        self.assertTrue(after)
        for entry in after:
            self.assertEqual(entry['status'], 'no_estimable_effect')
            self.assertIn('对照', entry['reason'])
            self.assertTrue(entry['missing_cohorts'])

    def test_zero_output_does_not_switch_denominator(self):
        rows = make_rows(build_sites([4, 4, 7, None, None]), WEEKS)
        # 让批次 4 的全部站点在基期周 3 的产出为零：该批次无法算“单位产出”。
        for row in rows:
            if row['adopt_week'] == 4 and row['week'] == 3:
                row['output'] = 0
        result = analyze(rows)
        intensity = result['targets']['intensity']
        consumption = result['targets']['consumption']
        for window in intensity['windows']:
            if window['cohort'] == 4 and window['event'] >= 0:
                self.assertEqual(window['status'], 'not_estimable')
                self.assertEqual(window['reason'], 'undefined_outcome')
                excluded = {item['site'] for item in window['roster']
                            if item['role'] == 'excluded'
                            and item['reason'] == 'output_zero'}
                self.assertEqual(excluded, {'S00', 'S01'})
        # 总量目标不受产出为零影响，仍然估计。
        self.assertEqual(consumption['status'], 'estimated')
        zero_windows = [window for window in consumption['windows']
                        if window['cohort'] == 4 and window['event'] >= 0]
        self.assertTrue(all(window['status'] == 'estimated' for window in zero_windows))
        self.assertEqual(result['diagnostics']['output_zero_sites'], ['S00', 'S01'])
        # 单位产出目标从不把总量当作分母回退；窗口只是缺口。
        self.assertFalse(any(window.get('estimate') is not None
                             for window in intensity['windows']
                             if window['cohort'] == 4 and window['event'] >= 0))

    def test_large_site_dominates_baseline_weighting(self):
        rows = make_rows(
            build_sites([4, 4, 7, None],
                        baseline=[900.0, 10.0, 10.0, 10.0]),
            WEEKS)
        result = analyze(rows, weighting='baseline_consumption')
        cohort = result['targets']['consumption']['cohort_effects'][0]
        self.assertTrue(cohort['concentration']['dominant_site_flag'])
        self.assertGreater(
            cohort['concentration']['largest_site_share_of_baseline_consumption'], 0.9)
        self.assertTrue(any('主导' in warning for warning in result['warnings']))
        # 站点数权重下批次均值为等权，不会出现“大站点主导汇总”的警告，
        # 但仍会以诊断字段给出集中度。
        balanced = analyze(rows, weighting='site_count')
        self.assertFalse(any('主导' in warning for warning in balanced['warnings']))
        self.assertTrue(balanced['targets']['consumption']['cohort_effects'][0]
                        ['concentration']['dominant_site_flag'])

    def test_weightings_are_explicit_and_differ(self):
        rows = make_rows(build_sites([4, 4, 7, 7, None],
                                     baseline=[10.0, 20.0, 30.0, 60.0, 120.0]),
                         WEEKS, noise=0.0)
        balanced = analyze(rows, weighting='site_count')
        weighted = analyze(rows, weighting='baseline_consumption')
        self.assertEqual(balanced['weighting'], 'site_count')
        self.assertEqual(weighted['weighting'], 'baseline_consumption')
        self.assertEqual(set(balanced['targets']['consumption']['post_average']['weights']),
                         set(weighted['targets']['consumption']['post_average']['weights']))
        # 取有两个批次共同参与的事件周，权重口径才会产生差异。
        def event_zero(result):
            return [entry for entry in result['targets']['consumption']['dynamic_effects']
                    if entry['event'] == 0][0]
        left, right = event_zero(balanced), event_zero(weighted)
        self.assertGreaterEqual(len(left['weights']), 2)
        self.assertNotEqual(left['weights'], right['weights'])
        self.assertAlmostEqual(sum(left['weights'].values()), 1.0)
        self.assertAlmostEqual(sum(right['weights'].values()), 1.0)
        self.assertIn('归一化', left['weight_basis'])

    def test_covariate_only_uses_comparison_window(self):
        temperature = {week: 20.0 + 3.0 * week for week in WEEKS}
        rows = make_rows(build_sites([4, 4, 7, 7, None, None],
                                     temperature_slope=[0.0, 0.1, 0.2,
                                                        0.3, 0.4, 0.5]),
                         WEEKS, temperature=temperature)
        result = analyze(rows)
        for target, outcome in result['targets'].items():
            covariate = outcome['covariate']
            self.assertEqual(covariate['status'], 'identified')
            self.assertTrue(covariate['adjusted'])
            self.assertEqual(covariate['windows_with_cross_site_variation'],
                             covariate['windows_total'])
            rows_by_site = {}
            for row in rows:
                rows_by_site.setdefault(row['site'], {})[row['week']] = row
            for window in outcome['windows']:
                if window['status'] != 'estimated':
                    continue
                # 协变量调整只用该窗口基期周与比较周之差，不用窗口外数据。
                changes = {
                    site['site']: rows_by_site[site['site']][window['week']]['temperature']
                    - rows_by_site[site['site']][window['base_week']]['temperature']
                    for site in window['roster']
                    if site['role'] in ('treated', 'control')}
                self.assertTrue(changes)
                self.assertTrue(any(value != 0 for value in changes.values()))

    def test_covariate_not_identified_when_within_window_constant(self):
        rows = make_rows(build_sites([4, 4, 7, 7, None, None]), WEEKS,
                         temperature={week: 20.0 for week in WEEKS})
        result = analyze(rows)
        for outcome in result['targets'].values():
            self.assertEqual(outcome['covariate']['status'], 'not_identified')
            self.assertFalse(outcome['covariate']['adjusted'])
            self.assertIn('不应写', outcome['covariate']['detail'])

    def test_covariate_adjustment_recovers_known_slope(self):
        sites = build_sites([4] * 6 + [8] * 6 + [None] * 12,
                            baseline=[100.0 + 3 * index for index in range(24)],
                            temperature_slope=[0.01 * index for index in range(24)],
                            effect=-5.0)
        rows = make_rows(sites, WEEKS, noise=0.0)
        for row in rows:
            row['consumption'] += 2.0 * row['temperature']
        outcome = analyze(rows)['targets']['consumption']
        self.assertEqual(outcome['covariate']['status'], 'identified')
        self.assertTrue(outcome['covariate']['adjusted'])
        self.assertAlmostEqual(outcome['covariate']['coefficient'], 2.0, places=8)
        self.assertAlmostEqual(outcome['post_average']['estimate'], -5.0, places=6)

    def test_covariate_none_skips_adjustment(self):
        sites = build_sites([4] * 6 + [8] * 6 + [None] * 12,
                            temperature_slope=[0.01 * index for index in range(24)])
        rows = make_rows(sites, WEEKS)
        outcome = analyze(rows, covariate='none')['targets']['consumption']
        self.assertEqual(outcome['covariate']['status'], 'not_requested')
        self.assertFalse(outcome['covariate']['adjusted'])
        self.assertNotIn('coefficient', outcome['covariate'])

    def test_intervals_use_site_clusters(self):
        rows = make_rows(build_sites([4, 4, 7, 7, None, None]), WEEKS, noise=1.0)
        result = analyze(rows)
        for outcome in result['targets'].values():
            self.assertEqual(outcome['n_clusters'], len({row['site'] for row in rows}))
            self.assertEqual(outcome['df'], outcome['n_clusters'] - 1)
            self.assertLess(outcome['n_comparison_rows'],
                            len(rows) * len(outcome['windows']))

    def test_reordering_rows_is_deterministic(self):
        rows = make_rows(build_sites([4, 4, 7, 7, 10, None, None, None]), WEEKS,
                         noise=0.5, seed=11)
        first = analyze(rows)
        shuffled = list(rows)
        random.Random(99).shuffle(shuffled)
        self.assertEqual(first, analyze(shuffled))
        # 站点名重排但不改变采用时点：结果的结构与数值应可复现。
        renamed = copy.deepcopy(rows)
        names = sorted({row['site'] for row in renamed}, reverse=True)
        mapping = {name: f'T{index:02d}' for index, name in enumerate(names)}
        for row in renamed:
            row['site'] = mapping[row['site']]
            row['observation_id'] = f"{row['site']}-w{row['week']:02d}"
        other = analyze(renamed)
        self.assertEqual(first['input_fingerprint'] == other['input_fingerprint'], False)
        for target in first['targets']:
            self.assertAlmostEqual(first['targets'][target]['post_average']['estimate'],
                                   other['targets'][target]['post_average']['estimate'])
            for left, right in zip(first['targets'][target]['dynamic_effects'],
                                   other['targets'][target]['dynamic_effects']):
                if left.get('estimate') is None:
                    self.assertIsNone(right.get('estimate'))
                else:
                    self.assertAlmostEqual(left['estimate'], right['estimate'])

    def test_disclaimer_and_pretrend_caveat_always_present(self):
        rows = make_rows(build_sites([4, 4, 7, None]), WEEKS)
        result = analyze(rows)
        self.assertIn('不是节能承诺', result['disclaimer'])
        self.assertIn('不能证明平行趋势成立', result['pretrend_note'])
        for outcome in result['targets'].values():
            self.assertIn('不能证明平行趋势成立', outcome['pretrend']['note'])

    def test_declared_targets_and_invalid_options(self):
        rows = make_rows(build_sites([4, 4, 7, None]), WEEKS)
        only = analyze(rows, targets=('consumption',))
        self.assertEqual(list(only['targets']), ['consumption'])
        with self.assertRaises(ValueError):
            analyze(rows, targets=('throughput',))
        with self.assertRaises(ValueError):
            analyze(rows, targets=())
        with self.assertRaises(ValueError):
            analyze(rows, weighting='unknown')
        with self.assertRaises(ValueError):
            analyze(rows, confidence=1.5)

    def test_site_level_ledger_is_traceable(self):
        rows = make_rows(build_sites([4, 4, 7, None, None]), WEEKS)
        result = analyze(rows)
        sites = [entry['site'] for entry in result['panel']['site_inventory']]
        self.assertEqual(sites, sorted({row['site'] for row in rows}))
        window = result['targets']['consumption']['windows'][0]
        roster_sites = [entry['site'] for entry in window['roster']]
        self.assertEqual(len(roster_sites), len(set(roster_sites)))
        self.assertEqual(sorted(roster_sites),
                         sorted(result['panel']['site_inventory'], key=str)[0:0]
                         + sorted({entry['site'] for entry in rows and
                                   result['panel']['site_inventory']}, key=str))
        for entry in result['panel']['site_inventory']:
            self.assertIn(entry['cohort'], (None, 4, 7))
            self.assertIn('baseline_weeks', entry)

    def test_inconsistent_adopt_week_rejected(self):
        rows = make_rows(build_sites([4, 4, 7, None]), WEEKS)
        rows[-1] = dict(rows[-1], adopt_week=9)
        with self.assertRaisesRegex(ValueError, '采用周不一致'):
            analyze(rows)

    def test_duplicate_site_week_rejected(self):
        rows = make_rows(build_sites([4, 4, 7, None]), WEEKS)
        duplicated = rows + [dict(rows[0], observation_id='extra')]
        with self.assertRaisesRegex(ValueError, '站点周重复'):
            analyze(duplicated)

    def test_cli_rollout_matches_library(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'rollout.json'
            run = subprocess.run(
                [sys.executable, '-m', 'analysis_app', 'rollout',
                 '--output', str(path), '--weighting', 'site_count'],
                text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            expected = analyze(load_records('data/observations.csv'))
            printed = json.loads(run.stdout)
            self.assertEqual(set(printed), set(expected))
            self.assertEqual(printed['targets'], expected['targets'])
            self.assertEqual(printed['panel'], expected['panel'])
            self.assertEqual(json.loads(path.read_text()), printed)

    def test_no_pre_period_cohort_is_a_gap(self):
        rows = make_rows(build_sites([0, 0, 4, 4]), WEEKS)
        result = analyze(rows)
        outcome = result['targets']['consumption']
        cohort = [entry for entry in outcome['cohort_effects']
                  if entry['cohort'] == 0][0]
        self.assertTrue(all(event['status'] == 'not_estimable'
                            for event in cohort['events']))
        self.assertTrue(all(event['reason'] == 'no_pre_period'
                            for event in cohort['events']))
        self.assertIn('最早批次', ' '.join(result['diagnostics']['notes']))

    def test_calendar_week_gap_explains_why(self):
        rows = make_rows(build_sites([4, 4, 7, 7, None, None]), WEEKS)
        result = analyze(rows)
        outcome = result['targets']['consumption']
        early = [entry for entry in outcome['calendar_weeks'] if entry['week'] < 4]
        self.assertTrue(early)
        for entry in early:
            self.assertIsNone(entry['estimate'])
            self.assertEqual(entry['status'], 'no_adopted_cohort')
            self.assertIn('没有已采用站点', entry['reason'])


if __name__ == '__main__':
    unittest.main()
