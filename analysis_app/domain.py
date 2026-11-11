import math
import statistics
from collections import Counter


def describe(records):
    output = math.fsum(row['output'] for row in records)
    energy = math.fsum(row['consumption'] for row in records)
    return {'consumption': energy, 'output': output, 'intensity': energy / output if output else None}
