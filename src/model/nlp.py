import pandas as pd
from nltk.corpus import stopwords
from utils.data_processing import df_resample

df = pd.read_csv('../data/interim/pubmed_data.csv')

def nlp():
    df = df_resample(df, save=False, **{'random_state': 20})
    pass