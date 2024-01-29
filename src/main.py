import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))

from utils.text_processing import cutter
from utils.data_processing import reduce_data
from utils.data_processing import monts_standardizing
from utils.data_processing import year_numeric
from scraping.scraper import search
from processing.build_df import df_pubmed

def input():
    #Introducing query
    mesh_string = input('introduce mesh terms separated by "/": ')
    
    #Processing introduced query
    mesh_terms = cutter(mesh_string)
    return mesh_terms

def main():
    mesh_terms = input()
    #Searching in PubMed
    studies = search(mesh_terms)
    studiesIdList = studies['IdList']
    
    #Reducing the number of samples
    studiesIdList = reduce_data(studiesIdList, 0.75)
    
    # Making a dataframe
    df = df_pubmed(studiesIdList)
    
    #Standardizing months values
    df = monts_standardizing(df)
    df = year_numeric(df)
    df.to_csv('data/raw/pubmed_data.csv')
    print("CSV file created")

if __name__ == "__main__":
    main()