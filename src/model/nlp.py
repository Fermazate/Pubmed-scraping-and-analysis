import sys
import os
import pandas as pd
from nltk.corpus import stopwords

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..')))

from utils.data_processing import df_resample

df = pd.read_csv('data/raw/pubmed_data.csv')

def nlp():
    df = df_resample(df, save=False, **{'frac': 0.5})
    