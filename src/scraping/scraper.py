from http.client import IncompleteRead
from Bio import Entrez
import numpy as np
import time

def search(query:str, db:str='pubmed'):
    Entrez.email = "fernandomazate@gmail.com"
    Entrez.tool = "PubMedScraper"
    handle = Entrez.esearch(db=db,
                            sort = 'relevance',
                            retmax = '100000',
                            retmode = 'xml',
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list, db:str='pubmed', max_attempts=3, delay=10):
    attempt = 0
    while attempt<max_attempts:
        try:
            ids = ','.join(id_list)
            Entrez.email = "fernandomazate@gmail.com"
            handle = Entrez.efetch(db=db,
                                    id=ids,
                                    retmode = 'xml')
            results = Entrez.read(handle)
            return results
        except IncompleteRead as e:
            print(f'Attemt {attempt+1} failed with error: {str(e)}. Retrying in {delay} seconds...')
            time.sleep(delay)
            attempt += 1
    raise Exception('Maximum number of attempts reached')