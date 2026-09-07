# Timeout expressed in seconds
CAMPAIGN_TIME_LIMITS = [60, 120, 240, 480]


POOL_SIZE = 1

# constraint mode and number of workers
CONFIGURATIONS = [
    ("hard", 1),
    ("soft", 1),
    ("hard", 2),
    ("soft", 2),
    ("hard", 4),
    ("soft", 4),
]
