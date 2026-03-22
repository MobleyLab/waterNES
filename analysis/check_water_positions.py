import pandas as pd
import os, sys
import pickle
import pathlib
import numpy as np
np.set_printoptions(suppress=True)
from typing import Optional, Tuple, Dict, List

class Water_Positions_Analysis:
    def __init__(self, df_xvg_file, current_stage):
        self.df_xvg_file = df_xvg_file
        self.current_stage = current_stage
        print(f"checking water positions for stage {current_stage}\n")
        self.df_xvg_file = self.engine()



    def engine(self):
        return self.df_xvg_file
