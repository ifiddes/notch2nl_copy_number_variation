"""
Calculates kmer weights for a cohort of individuals.
Expects that there are jellyfish Counts files in the output directory structure.
"""
import sys, os
import cPickle as pickle
from sonLib.bioio import fastaRead
from collections import defaultdict

#hard coded allele fractions seen in 281 TCGA individuals
avg_c_frac = 1.0 * (2 * 250 + 1 * 18) / 281
avg_d_frac = 1.0 * (2 * 130 + 1 * 121) / 281
avg_frac_dict = {"A":2.0, "B":2.0, "C":avg_c_frac, "D":avg_d_frac, "2":2.0}

G = pickle.load(open("data/graphs/OriginalWithOffsets.pickle"))

count_files = [os.path.join("output", x, x + ".Counts.fa") for x in os.listdir("output")]
counts = defaultdict(int)
for i, f in enumerate(count_files):
    print "working on {}. {} to go.".format(f, len(count_files) - i)
    for count, kmer in fastaRead(f):
        counts[kmer] += int(count)

with open("raw_population_counts.pickle", "w") as outf:
    pickle.dump(counts, outf)

weights = {}
for kmer in G.G.nodes():
    input_sequences = [x.split("_")[0][-1] for x in G.G.node[kmer]['label'].split(", ")]
    population_counts = counts[kmer]
    weight = (1.0 * population_counts) / sum(avg_frac_dict[x] for x in input_sequences)
    weights[kmer] = weight

with open("original_graph_weights_ONLY.pickle", "w") as outf:
    pickle.dump(weights, outf)

for kmer in G.G.nodes():
    G.G.node[kmer]['weight'] = weights[kmer]

with open("OriginalWithOffsets_WITH_WEIGHTS.pickle", "w") as outf:
    pickle.dump(G, outf)