from os import path, mkdir, listdir
from subprocess import call
from shutil import rmtree, copyfile
from Bio import SeqIO
from argparse import ArgumentParser
from timeit import default_timer
from itertools import combinations, product, chain
from scipy.spatial.distance import pdist, cdist
from scipy.cluster.hierarchy import linkage, fcluster
import pandas as pd
from gc import collect
# from numpy import average, std
# import matplotlib.pyplot as plt
import random
import string
from warnings import warn
from numpy import where
from math import ceil


def read_fasta_sequences(seqdir):
    list_of_files = listdir(seqdir)
    sequences = {}
    counter = 1
    for filename in list_of_files:
        fileid = 'fid' + str(counter).zfill(7) + '_'
        with open(seqdir + '/' + filename, 'rU') as fasta_file:
            # sequences.update(SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta')))
            records = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
            for key in list(records.keys()):
                newkey = fileid + key
                sequences[newkey] = records[key]
                sequences[newkey].id = newkey
        counter += 1
    for key in list(sequences.keys()):
        if len(sequences[key]) > 10000:
            del sequences[key]
            warn('WARNING: Sequence ' + key + ' is too long (>10000). It will be ignored in the analysis.')
    return sequences


def extract_features(sequences, chunksize, storepath):
    keys = list(sequences.keys())
    AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    cutoffs = list(range(0, len(keys), chunksize))
    maxlength = 40
    if len(keys) not in cutoffs:
        cutoffs.append(len(keys))
    # with pd.get_store(storepath) as store:
    with pd.HDFStore(storepath) as store:
        for i, s in enumerate(cutoffs[:-1]):
            e = cutoffs[i + 1]
            features = pd.DataFrame(0, index=keys[s:e], columns=AAs)
            lengths = pd.Series(0, index=keys[s:e])
            for key in keys[s:e]:
                seq = sequences[key].seq.upper()
                length = len(seq)
                for item in seq:
                    if item in AAs:
                        features.at[key, item] += 1
                    elif item == 'B':
                        features.at[key, 'N'] += 0.5
                        features.at[key, 'D'] += 0.5
                    elif item == 'Z':
                        features.at[key, 'Q'] += 0.5
                        features.at[key, 'E'] += 0.5
                    elif item == 'J':
                        features.at[key, 'I'] += 0.5
                        features.at[key, 'L'] += 0.5
                    else:
                        length -= 1
                lengths.at[key] = length
            features = features.divide(lengths, axis='index')
            store.append('features', features, min_itemsize={'index': maxlength})
            # store.append('lengths', lengths, min_itemsize={'index': maxlength})


def calculate_distances(storepath, lenfeatures, chunksize, thresh, outpath):
    eps = 1e-5
    cutoffs = list(range(0, lenfeatures, chunksize))
    if lenfeatures not in cutoffs:
        cutoffs.append(lenfeatures)
    # with pd.get_store(storepath) as store:
    with pd.HDFStore(storepath) as store:
        for i, s1 in enumerate(cutoffs[:-1]):
            e1 = cutoffs[i + 1]
            for j, s2 in enumerate(cutoffs[i:-1]):
                e2 = cutoffs[i + j + 1]
                if s1 == s2:
                    chunk = store.select('features', start=s1, stop=e1)
                    distances = pdist(chunk, 'euclidean')
                    tuples = combinations(chunk.index, 2)
                else:
                    chunk1 = store.select('features', start=s1, stop=e1)
                    chunk2 = store.select('features', start=s2, stop=e2)
                    distances = cdist(chunk1, chunk2, 'euclidean')
                    distances = list(chain.from_iterable(distances))
                    tuples = product(chunk1.index, chunk2.index)
                tuplesdf = pd.DataFrame.from_records(tuples)
                distancesdf = pd.Series(distances) + eps
                tupledist = pd.concat([tuplesdf, distancesdf], axis=1)
                tupledist.columns = ['g1', 'g2', 'dist']
                tuplesmalldist = tupledist.loc[tupledist['dist'] < thresh]
                tuplesmalldist.to_csv(outpath, sep='\t', header=False, index=False, mode='a')
                collect()  # Added to avoid memory leak and consequently freezing of system


def run_mcl(inpath, outpath):
    try:
        call(['mcl', inpath, '--abc', '-o', outpath], shell=False)
    except:
        raise SystemExit('Unable to run MCL. Check if:\n 1- ' + inpath + ' exists.\n 2- MCL is installed.')


def separate_sequences(sequences, clusters, outpath):
    counter = 0
    with open(clusters, 'r') as fromfile:
        write = SeqIO.write
        for line in fromfile:
            cells = line.split()
            if len(cells) > 1:
                selectedsequences = [sequences[x] for x in cells]
                address = outpath + '/seq' + str(counter) + '.txt'
                write(selectedsequences, address, 'fasta')
                counter += 1


def run_phmmer(evalue, seqs, phmmers):
    list_of_files = listdir(seqs)
    for filename in list_of_files:
        seqdb = seqs + '/' + filename
        outpath = phmmers + '/' + filename
        try:
            call(['phmmer', '-o', '/dev/null', '-E', str(evalue), '--tblout', outpath, seqdb, seqdb], shell=False)
        except:
            raise SystemExit('Unable to run phmmer. Check if:\n 1- ' + seqdb + ' exists.\n 2- HMMER is installed.')


def parse_tblouts(inpath, outpath, storepath, chunksize): #rewrite without using pandas
    maxlength = 40
    with open(outpath, 'w') as tofile:
            list_of_files = listdir(inpath)
            write = tofile.write
            for filename in list_of_files:
                edgelist = list()
                address = inpath + '/' + filename
                with open(address, 'r') as fromfile:
                    for line in fromfile:
                        if not line.startswith('#'):
                            cells = line.split()
                            seq1 = min(cells[0], cells[2])
                            seq2 = max(cells[0], cells[2])
                            score = float(cells[5])
                            edgelist.append([seq1, seq2, score])
                edgelist.sort()
                i = 0
                while i < len(edgelist):
                    seq1 = edgelist[i][0]
                    seq2 = edgelist[i][1]
                    if seq1 == seq2:
                        i += 1
                        continue
                    score = float(edgelist[i][2])
                    j = i + 1
                    while j < len(edgelist) and edgelist[i][0:2] == edgelist[j][0:2]:
                        score += float(edgelist[j][2])
                        j += 1
                    steps = j - i
                    score /= steps
                    score = 1 / (score + 1)
                    i += steps
                    write(seq1 + '\t' + seq2 + '\t' + str(score) + '\n')

    with pd.HDFStore(storepath) as store:
        with open(outpath, 'r') as fromfile:
            counter = 1
            gene_dict = {'g1':[], 'g2':[], 'score':[]}
            for line in fromfile:
                cells = line.split()
                gene_dict['g1'].append(cells[0])
                gene_dict['g2'].append(cells[1])
                gene_dict['score'].append(float(cells[2]))
                counter += 1
                if not(counter % chunksize):
                    distances = pd.DataFrame(gene_dict)
                    store.append('distances', distances, min_itemsize={'g1': maxlength, 'g2': maxlength}, index=False,
                                data_columns=['g1', 'g2'])
                    gene_dict = {'g1': [], 'g2': [], 'score': []}
            if counter % chunksize:
                distances = pd.DataFrame(gene_dict)
                store.append('distances', distances, min_itemsize={'g1': maxlength, 'g2': maxlength}, index=False,
                             data_columns=['g1', 'g2'])
###################edit
                    # index = pd.MultiIndex.from_tuples([(seq1, seq2)], names=['g1', 'g2'])
                    # distances = pd.Series(score, index=index)
                    # store.append('distances', distances, min_itemsize={'g1': maxlength, 'g2': maxlength})
                    #https://stackoverflow.com/questions/22934996/pandas-pytables-append-performance-and-increase-in-file-size
                    #https://stackoverflow.com/questions/22777284/improve-query-performance-from-a-large-hdfstore-table-with-pandas


def run_hierarchical(storepath, mcl, outpath):
    with pd.HDFStore(storepath) as store:
        with open(mcl, 'r') as fromfile:
            with open(outpath, 'w') as tofile:
                for line in fromfile:
                    cells = line.split()
                    fileids = {item[0:10] for item in cells}
                    num_genomes = len(fileids)
                    n = len(cells)
                    if num_genomes == 1:
                        for item in cells:
                            tofile.write(item[11:] + '\t\n')
                    elif n > num_genomes:
                        comb = int(n * (n - 1) / 2)
                        distarray = [0] * comb
                        cellsproduct = list(combinations(cells, 2))
                        i = 0
                        for item in cellsproduct:
                            minim = min(item[0], item[1])
                            maxim = max(item[0], item[1])
                            selectres = store.select('distances', where=("g1==minim", "g2==maxim"), columns=['score'])
                            if selectres.empty:
                                distarray[i] = 1
                            else:
                                distarray[i] = float(selectres['score'])
                            i += 1
                        linkageres = linkage(distarray)
                        clusters = fcluster(linkageres, ceil(n/num_genomes), criterion='maxclust')
                        for i in range(1, max(clusters) + 1):
                            current_cluster = where(clusters == i)[0]
                            # fileids_c = {cells[item][0:10] for item in current_cluster}
                            # n_c = len(current_cluster)
                            # num_genomes_c = len(fileids_c)
                            # if num_genomes_c == 1:
                            #     for item in current_cluster:
                            #         tofile.write(cells[item][11:] + '\t\n')
                            # else:
                            for item in current_cluster:
                                tofile.write(cells[item][11:] + '\t')
                            tofile.write('\n')
                    else:
                        for item in cells:
                            tofile.write(item[11:] + '\t')
                        tofile.write('\n')

def add_singletons(sequences, clust_no_sings, outpath):
    genes = set()
    update = genes.update
    with open(clust_no_sings, 'r') as fromfile:
        for line in fromfile:
            cells = line.split()
            update(cells)
    copyfile(clust_no_sings, outpath)
    with open(outpath, 'a') as tofile:
        write = tofile.write
        keys = sequences.keys()
        keys = {item[11:] for item in keys}
        for key in keys:
            if key not in genes:
                write(key + '\n')


def pairwise_orthology_call(seq, outpath, chunksize, thresh, evalue, temp):
    sequences = read_fasta_sequences(seq)
    storepath = temp + '/porc.h5'
    extract_features(sequences, chunksize, storepath)
    distances = temp + '/distances.txt'
    lenfeatures = len(list(sequences.keys()))
    calculate_distances(storepath, lenfeatures, chunksize, thresh, distances)
    print('Calculated distances')
    initclusters = temp + '/' + 'initial_clusters.txt'
    run_mcl(distances, initclusters)
    seqs = temp + '/sequences'
    mkdir(seqs)
    separate_sequences(sequences, initclusters, seqs)
    phmmers = temp + '/phmmers'
    mkdir(phmmers)
    run_phmmer(evalue, seqs, phmmers)
    edges = temp + '/edges.txt'
    parse_tblouts(phmmers, edges, storepath, chunksize)
    mcl = temp + '/mcl.txt'
    run_mcl(edges, mcl)
    clust_no_sings = temp + '/no-singletons.txt'
    run_hierarchical(storepath, mcl, clust_no_sings)
    add_singletons(sequences, clust_no_sings, outpath)


def main():
    start = default_timer()
    parser = ArgumentParser(description='Clusters orthologous proteins. Needs Python 2.7 or higher, MCL, and HMMER 3')
    parser.add_argument('p', help='Directory containing all proteomes in fasta format. Each proteome should be in a separate file.')
    parser.add_argument('o', help='Path to the output file')
    parser.add_argument('-e', '--eva', help='E-value threshold for phmmer', default=1e-10)
    parser.add_argument('-s', '--thr', help='Similarity threshold', default=0.05)
    parser.add_argument('-c', '--chnk', help='Chunk size', default=5000)
    args = parser.parse_args()
    seq = args.p
    outpath = args.o
    evalue = float(args.eva)
    thresh = float(args.thr)
    chunksize = int(args.chnk)
    directory, filename = path.split(outpath)
    if directory == '':
        directory = '.'
    temp = directory + '/tmp' + ''.join(random.choice(string.digits) for _ in range(6))
    mkdir(temp)
    pairwise_orthology_call(seq, outpath, chunksize, thresh, evalue, temp)
    #rmtree(temp)
    stop = default_timer()
    print(stop - start)


if __name__ == '__main__':
    main()###http://pandas-docs.github.io/pandas-docs-travis/io.html#io-hdf5

# Submit to cluster with memory requirements: qsub -l h_vmem=32g -l h_rt=480:00:00 run.sh