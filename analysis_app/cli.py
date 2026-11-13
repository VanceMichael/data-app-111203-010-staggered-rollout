import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

from .io import load_records
from .rollout import TARGETS, analyze
from .summary import summarize


def save_json(path, result):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix='.analysis-', dir=destination.parent)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, sort_keys=True, allow_nan=False, indent=2)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _emit(result, output):
    if output:
        save_json(output, result)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, allow_nan=False))


def main():
    parser = argparse.ArgumentParser(description='本地观测数据分析')
    commands = parser.add_subparsers(dest='command', required=True)
    summary = commands.add_parser('summary', help='按业务维度核对观测口径')
    summary.add_argument('--input', type=Path, default=Path('data/observations.csv'))
    summary.add_argument('--output', type=Path)
    rollout = commands.add_parser('rollout', help='交叠上线的批次×事件周比较分析')
    rollout.add_argument('--input', type=Path, default=Path('data/observations.csv'))
    rollout.add_argument('--output', type=Path)
    rollout.add_argument('--targets', default=','.join(TARGETS),
                         help='逗号分隔的估计目标：consumption,intensity；默认两者都算')
    rollout.add_argument('--covariate', choices=['temperature', 'none'],
                         default='temperature',
                         help='比较窗口内的协变量调整；none 表示不调整')
    rollout.add_argument('--weighting', choices=['site_count', 'baseline_consumption'],
                         default='site_count',
                         help='跨窗口汇总的权重口径，默认按参与站点等权')
    rollout.add_argument('--confidence', type=float, default=0.95,
                         help='区间与检验的置信水平，默认 0.95')
    args = parser.parse_args()
    try:
        if args.command == 'summary':
            result = summarize(load_records(args.input))
        else:
            targets = tuple(item.strip() for item in args.targets.split(',') if item.strip())
            result = analyze(load_records(args.input), targets=targets,
                             covariate=args.covariate, confidence=args.confidence,
                             weighting=args.weighting)
        _emit(result, args.output)
    except (ValueError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
