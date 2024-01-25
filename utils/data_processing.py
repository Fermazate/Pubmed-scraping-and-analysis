import numpy as np
import pandas as pd
import random

random.seed(10)


def reduce_data(id_list, perc:float=0.75):
    """
    This function randomly reduces the number of samples in a dataframe by a percentage.
    """
    df_reduced = random.sample(id_list, int(perc*len(id_list)))
    return df_reduced

def monts_standardizing(df:pd.DataFrame):
    """
    This function standardizes the data in a dataframe.
    """
    df['Month'].replace('Jan', '01', inplace=True)
    df['Month'].replace('Feb', '02', inplace=True)
    df['Month'].replace('Mar', '03', inplace=True)
    df['Month'].replace('Apr', '04', inplace=True)
    df['Month'].replace('May', '05', inplace=True)
    df['Month'].replace('Jun', '06', inplace=True)
    df['Month'].replace('Jul', '07', inplace=True)
    df['Month'].replace('Aug', '08', inplace=True)
    df['Month'].replace('Sep', '09', inplace=True)
    df['Month'].replace('Oct', '10', inplace=True)
    df['Month'].replace('Nov', '11', inplace=True)
    df['Month'].replace('Dec', '12', inplace=True)
    df['Month'].replace('No Data', np.nan, inplace=True)
    
    return df

def year_numeric(df:pd.DataFrame):
    """
    This function converts the year column in a dataframe to numeric.
    """
    df['Year'].replace('No Data', np.nan, inplace=True)
    df['Year'] = round(pd.to_numeric(df['Year'], errors='coerce'), 0)
    return df

def df_resample(df:pd.DataFrame, save:bool=True, save_path:str='../data/raw/pubmed_data_resampled.csv', **kwargs):
    df = df.sample(frac=1, random_state=20 **kwargs)
    if save == True:
        df.to_csv(save_path, index=False)
    else:
        pass
    return df