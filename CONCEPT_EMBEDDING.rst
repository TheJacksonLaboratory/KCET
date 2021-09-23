######################
Word/concept embedding
######################

The word embedding method based on the word2vec algorithm is performed on the preprocessed corpus to embed words to vectors. 
We used the `EMBeddInG GENerator (embiggen)<https://pypi.org/project/embiggen/>`_, a Python 3 software library developed by our group for word and 
graph embedding. 



.. code-block:: python
  
  import silence_tensorflow.auto # Import needed to avoid TensorFlow warnings and general useless infos.
  import pandas as pd
  import numpy as np
  from embiggen import CorpusTransformer
  from embiggen import Word2VecSequence
  from tensorflow.distribute import MirroredStrategy
  from embiggen import SkipGram
  from tensorflow.keras.optimizers import Nadam
  from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
  
Parameters
^^^^^^^^^^
  
In the current project, the skip-gram model was used for word2vec with the parameters window size = 5, 
minimum count (minimum word frequency) = 5, batch size = 128, negative samples = 20 and dimension = 100. 
Word embedding on the total corpus resulted in embeddings of 293,274 words each with dimension 100.

.. code-block:: python

  batch_size = 128
  embedding_size = 100
  window_size=5
  negative_samples=20
  delta=0.0001
  patience=20
  min_counts=5
  
PubMed abstract data
^^^^^^^^^^^^^^^^^^^^

We used our `marea <https://github.com/TheJacksonLaboratory/marea>`_ package to filter PubMed articles for relevance, apply PubTator Central concept recognition to the titles and abstracts of relevant articles, to remove stopwords, and to 
format the abstract data in preparation for word embedding. See the documentation of marea for further details. In this project, we used marea to extract
and process pubmed data up to and including three different target dates: 2010 (shown here), 2014, and 2020. The target date should match with the data for
training/validation data data used for creating the random forest classifier in later steps.

The output of marea is shown here as ``pubmed_cr.tsv`` -- adjust the path as required.
  
  
.. code-block:: python

  year = 2010
  pubmed_data = pd.read_csv("pubmed_cr.tsv",header = None,sep='\t')
  pubmed_year = pubmed_data.loc[pubmed_data.loc[:,1] <= year]
  pubmed = pubmed_year.loc[:,2].tolist()
  
  
Word/concept embedding
^^^^^^^^^^^^^^^^^^^^^^

The following code uses embiggen to perform embedding.

.. code-block:: python

  transformer = CorpusTransformer(tokenizer_method = "space",
        apply_stemming=False, 
        verbose=False,
        remove_stop_words=False,
        remove_punctuation = False,
        min_word_length=0,
        to_lower_case=False,
        min_count = min_counts, 
        min_sequence_length=window_size*2+1)
  transformer.fit(pubmed)
  encoded_pubmed = transformer.transform(pubmed)
  np.save(f"encoded_pubmed.npy",encoded_pubmed)
  print(f"The transformer will use {transformer.vocabulary_size} different terms.")
  word2vec_sequence = Word2VecSequence(
    encoded_pubmed,
    batch_size=batch_size,
    window_size=window_size,
    support_mirror_strategy = True
  )
  strategy = MirroredStrategy()
  with strategy.scope():
    model = SkipGram(
        vocabulary_size=transformer.vocabulary_size,
        embedding_size=embedding_size,
        optimizer=Nadam(0.01),
        window_size=window_size,
        negative_samples=negative_samples,
    )
  history = model.fit(
    word2vec_sequence,
    steps_per_epoch=word2vec_sequence.steps_per_epoch,
    epochs=1000,
    callbacks=[
        EarlyStopping(
            "loss",
            min_delta=delta,
            patience=patience,
            restore_best_weights=True
        ),
        ReduceLROnPlateau(monitor="loss", patience=patience//4)
    ]
  )


Saving the embeddings
^^^^^^^^^^^^^^^^^^^^^

Finally, we save the results of embedding to files.

.. code-block:: python

  np.save(f"embedding.npy", model.embedding)
  with open("words.txt", "w") as f_write:
    for i in range(model.embedding.shape[0]):
        f_write.write("{}\n".format(transformer.reverse_transform([[i]])))
