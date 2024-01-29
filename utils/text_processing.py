import re
import spacy
import string
import pandas as pd
import numpy as np
from collections import Counter
from nltk.corpus import stopwords
from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer

nlp = spacy.load("en_core_web_sm")

def cutter(string:str):
    """
    This function cuts a string in multiple parts.
    """
    cutted_string = string.split("/")
    cutted_string = [cutted.strip() for cutted in cutted_string]
    return cutted_string

def count_words(df:pd.DataFrame, column:str, chunk_size=10000):
    """
    This function counts the words in a dataframe.
    """
    word_counts = Counter()
    
    num_chunks = max(1, len(df) // chunk_size)
    
    for text in df[column].dropna():
        doc = nlp(text)
        words = [token.text.lower() for token in doc if token.text.isalpha()]
        word_counts.update(words)
        
    freq_df = pd.DataFrame(word_counts.items(), columns=['Word', 'Frequency'])
    freq_df = freq_df.sort_values(by='Frequency', ascending=False).reset_index(drop=True)
    
    return freq_df

def process_text(text, nlp_model):
    #Split the text into chunks
    max_lenght = 1000000
    chunks = [text[i:i+max_lenght] for i in range(0,len(text), max_lenght)]
    
    words = []
    for chunk in chunks:
        doc = nlp_model(chunk)
        words.extend([token.text.lower() for token in doc if token.is_alpha])
    return words

def stopwords_removal(text, words:list = [], synonyms:list = []):
    #Combine stopwords with other excluded words
    remove_list = set(stopwords.words('english') + words)
    
    #Pre-precessing synonyms
    synms_dict = {}
    
    for synms_list in synonyms:
        first_word = synms_list[0].lower()
        for synm in synms_list[1:]:
            synms_dict[synm.lower()] = first_word
    
    #Processing each word
    
    non_stopwords = []
    
    nopunc = [char for char in text if char not in string.punctuation]
    nopunc = ''.join(nopunc)
    
    for word in nopunc.split():
        word_lower = word.lower()
        if word_lower in synms_dict:
            word_lower = synms_dict[word_lower]
        if word_lower not in remove_list:
            non_stopwords.append(word_lower)
    return non_stopwords

def counting(df, column:str, words:list=[], synonyms:list = []):
    vect = CountVectorizer(analyzer= (lambda text: stopwords_removal(text, words=words, synonyms=synonyms))).fit(df[column])
    
    bow = vect.transform(df[column])
    
    feature_names = vect.get_feature_names_out()
    
    df_bow = pd.DataFrame(bow.todense(), columns=feature_names)
    
    return df_bow, bow, feature_names

def make_weighted_table(bow, feature_names):
    tfid_transformer = TfidfTransformer().fit(bow)
    tfidf = tfid_transformer.transform(bow)
    
    df_tfidf = pd.DataFrame(tfidf.todense(), columns=feature_names)
    
    word_cnts = np.asarray(bow.sum(axis=0)).ravel().tolist()
    df_cnts = pd.DataFrame({'word':feature_names,'count':word_cnts})
    df_cnts = df_cnts.sort_values(by='count', ascending = False)
    
    weights = np.asarray(df_tfidf.mean(axis = 0)).ravel().tolist()
    df_weights = pd.DataFrame({'word':feature_names,'weight':weights})
    df_weights = df_weights.sort_values(by='weight', ascending = False)
    
    df_weights = df_weights.merge(df_cnts, on='word', how = 'left')
    
    return df_weights