import os
import sys
import re
import numpy as np
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import embedding_plotter as eplt

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

DO_TSNE = True

DO_INDIVIDUAL_EVENTS = True
DO_COMBINED_EVENTS = True

TSNE_LR_LIST = 'all'
TSNE_PP_LIST = 'all'

PLOT_DIR = os.path.join(config.plot_dir, 'altsplice')
EMBED_DIR = os.path.join(config.embed_dir, 'altsplice')
EVENT_LIST = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'concatenated']

TSNE_PP_PLOT_SET  = 'all'
TSNE_LR_PLOT_SET  = 'all'


