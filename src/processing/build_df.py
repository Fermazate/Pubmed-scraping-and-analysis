#Importing modules.
import pandas as pd

#Importing custom modules
from scraping.scraper import fetch_details

def df_pubmed(studies_id_list, db:str='pubmed'):
    title_list = []
    abstract_list = []
    journal_list = []
    language_list = []
    pubdate_year_list = []
    pubdate_month_list = []
    publicationtype_list = []
    articledate_list = []
    daterevised_list = []
    
    chunk_size = 10000
    
    for chunk_i in range (0, len(studies_id_list), chunk_size):
        chunk = studies_id_list[chunk_i:chunk_i + chunk_size]
        papers = fetch_details(chunk, db)
        for i, paper in enumerate (papers['PubmedArticle']):
            title_list.append(paper['MedlineCitation']['Article']['ArticleTitle'])
            try:
                abstract_list.append(paper['MedlineCitation']['Article']['Abstract']['AbstractText'])
            except:
                pass
            journal_list.append(paper['MedlineCitation']['Article']['Journal']['Title'])
            language_list.append(paper['MedlineCitation']['Article']['Language'][0])
            publicationtype_list.append(paper['MedlineCitation']['Article']['PublicationTypeList'])
            articledate_list.append(paper['MedlineCitation']['Article']['ArticleDate'])
            try:
                pubdate_year_list.append(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'])
            except:
                pubdate_year_list.append('No Data')
            try:
                pubdate_month_list.append(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month'])
            except:
                pubdate_month_list.append('No Data')
            daterevised_list.append(paper['MedlineCitation']['DateRevised'])
    df = pd.DataFrame(list(zip(
                                title_list, abstract_list, journal_list, language_list, pubdate_year_list, pubdate_month_list
                                )),
                                columns = ['title', 'abstract', 'journal', 'language', 'year', 'month']
                                )
    
    return df