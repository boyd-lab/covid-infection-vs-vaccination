"""Global configuration."""
import os
import numpy as np
import random

# Struct-like: namedtuple (https://stackoverflow.com/a/45426493/130164) or simplenamespace (https://dbader.org/blog/records-structs-and-data-transfer-objects-in-python)
from types import SimpleNamespace

reference_timepoint_wuhan_comparison = "day 21 / weeks 2&3"
reference_timepoint_variant_comparison = "day 28 / week 4"

# Set global random seed for reproducibility
random_seed = 0
np.random.seed(random_seed)
random.seed(random_seed)

# Configure file paths
paths = SimpleNamespace()
paths.data_dir = "data"
paths.generated_data_dir = "data/generated"
paths.output_dir = "out"

# convert all to absolute paths - consider the above relative to where this script lives, not where it's called from
for key, relative_path in paths.__dict__.items():
    # this would be relative to where config.py is imported from!
    # full_path = os.path.abspath(relative_path)

    # apply path relative to where config.py lives, not where it's imported from:
    dirname = os.path.dirname(__file__)
    # go one more level up
    dirname = os.path.dirname(dirname)
    paths.__dict__[key] = os.path.abspath(os.path.join(dirname, relative_path))
