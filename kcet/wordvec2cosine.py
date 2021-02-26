import pandas as pd
import numpy as np
from itertools import islice
from collections import defaultdict
from scipy.spatial.distance import cosine

class Wordvec2Cosine:


    def __init__(self, embeddings, words) -> None:
        embedding = np.load(embeddings, mmap_mode=None, allow_pickle=False, fix_imports=True, encoding='ASCII')
        word_list = []
        with open(words) as f:
            for line in f:
                word = line[2:-3]
                word_list.append(word)
        self._df = pd.DataFrame(data=embedding,index = word_list)

    def get_embeddings(self) -> pd.DataFrame:
        return self._df

    @staticmethod
    def _take(n, iterable):
        """
        Return the first n items of the iterable as a listsource kcet
        """
        return list(islice(iterable, n))

    def n_most_similar_words(self, target_word, n):
        """
        Returns a list with the top n words most similar to the target word
        """
        cosine_similarities = defaultdict()
        for word in self._df.index:
            cosine_similarity = 1 - cosine(self._df.loc[target_word],self._df.loc[word])
            cosine_similarities[word] = cosine_similarity
        sorted_cosin_similarities = {k: v for k, v in sorted(cosine_similarities.items(), key=lambda item: item[1], reverse=True)}
        n_items = Wordvec2Cosine._take(n, sorted_cosin_similarities.items())
        return n_items 

    def n_most_similar_words_df(self, target_word, n):
        n_items = self.n_most_similar_words(target_word=target_word, n=n)
        return pd.DataFrame(n_items, columns=["word","similarity"])

    def n_least_similar_words(self, target_word, n): 
        cosine_similarities = defaultdict()
        for word in self._df.index:
            cosine_similarity = 1 - cosine(self._df.loc[target_word],self._df.loc[word])
            cosine_similarities[word] = cosine_similarity
        sorted_cosin_similarities = {k: v for k, v in sorted(cosine_similarities.items(), key=lambda item: item[1], reverse=False)}
        n_items = Wordvec2Cosine._take(n, sorted_cosin_similarities.items())
        return n_items   

    def n_least_similar_words_df(self, target_word, n):
        n_items = self.n_least_similar_words(target_word=target_word, n=n)
        return pd.DataFrame(n_items, columns=["word","similarity"])

    def n_close_to_zero_similar_words(self, target_word, n, e):
        cosine_similarities = defaultdict()
        for word in self._df.index:
            cosine_similarity = 1 - cosine(self._df.loc[target_word],self._df.loc[word])
            if np.abs(cosine_similarity) < e:
                cosine_similarities[word] = cosine_similarity

        sorted_cosin_similarities = {k: v for k, v in sorted(cosine_similarities.items(), key=lambda item: item[1], reverse=False)}
        if n < len(sorted_cosin_similarities):
            n_items = Wordvec2Cosine._take(n, sorted_cosin_similarities.items())
        else:
            n_items = Wordvec2Cosine._take(len(sorted_cosin_similarities), sorted_cosin_similarities.items())
        return n_items

    def n_close_to_zero_similar_words_df(self, target_word, n, e):
        n_items = self.n_close_to_zero_similar_words(target_word=target_word, n=n, e=e)
        return pd.DataFrame(n_items, columns=["word","similarity"])