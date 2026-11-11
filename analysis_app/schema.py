FIELDS = {'observation_id': 'str', 'site': 'str', 'region': 'str', 'week': 'int', 'consumption': 'float', 'output': 'float', 'adopt_week': 'int', 'temperature': 'float'}
REQUIRED = ['observation_id', 'site', 'region', 'consumption']
ENTITY = 'site'
GROUP = 'region'
MEASURE = 'consumption'
PRODUCT = 'staggered-rollout'
UNIT = '千瓦时'
