

import numpy as np
from scipy.io import loadmat  # this is the SciPy module that loads mat-files
import matplotlib.pyplot as plt
from datetime import datetime, date, time
import pandas as pd

mat = loadmat("/home/aaron/git/Umea/semanticPredictors/NewGOA_code/NewGOA_Demo/YeastGOA_R.mat") # load mat-file



for key, value in mat.items():
    print(key,":",value)