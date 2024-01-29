from Bio import Entrez
import numpy as np

def search(query:str, db:str='pubmed'):
    Entrez.email = "fernandomazate@gmail.com"
    handle = Entrez.esearch(db=db,
                            sort = 'relevance',
                            retmax = '250000',
                            retmode = 'xml',
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list, db:str='pubmed'):
    ids = ','.join(id_list)
    Entrez.email = "fernandomazate@gmail.com"
    handle = Entrez.efetch(db=db,
                            id=ids,
                            retmode = 'xml')
    results = Entrez.read(handle)
    return results