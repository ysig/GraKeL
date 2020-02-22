"""
==============================================================================
Retrieval of most similar document using the Weisfeiler-Lehman subtree kernel.
==============================================================================
Script makes use of :class:`grakel.WeisfeilerLehman`, :class:`grakel.VertexHistogram`
"""
from __future__ import print_function
print(__doc__)

import numpy as np
import time

from nltk import word_tokenize
from nltk.corpus import sentence_polarity

from grakel.kernels import WeisfeilerLehman, VertexHistogram
from grakel import Graph

sents = sentence_polarity.sents()
sents = [sent for sent in sents if len(sent) > 1]
n_sents = 3000
sents = sents[:n_sents]
print("Loaded %d sentences\n" % n_sents)

print("Creating word co-occurrence networks\n")
word_networks = list()
for sent in sents:

    node_labels = dict()
    tokens_to_ids = dict()
    for token in sent:
        if token not in tokens_to_ids:
            tokens_to_ids[token] = len(tokens_to_ids)
            node_labels[tokens_to_ids[token]] = token
     
    edges = list()
    for i in range(len(sent)-1):
        edges.append((tokens_to_ids[sent[i]], tokens_to_ids[sent[i+1]]))

    word_networks.append(Graph(edges, node_labels=node_labels))

query_sent_id = 54
query_sent = [word_networks[query_sent_id]]

# Initialize Weisfeiler-Lehman subtree kernel
gk = WeisfeilerLehman(niter=2, normalize=True, base_kernel=VertexHistogram)

print("Computing similarities\n")
t0 = time.time()
gk.fit(query_sent)
K = gk.transform(word_networks)
print("done in %0.3fs\n" % (time.time() - t0))

print("Query sentence")
print("--------------")
print(" ".join(sents[query_sent_id]))
print()
print("Most similar sentence")
print("---------------------")
print(" ".join(sents[np.argsort(K[:,0])[-2]]))