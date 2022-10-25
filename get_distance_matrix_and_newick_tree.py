import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import math
import shlex, subprocess
from collections import Counter
from scipy.spatial.distance import hamming
import time
import argparse, sys
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA, TruncatedSVD
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform, hamming
import scipy
from scipy.sparse import csr_matrix, dok_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from joblib import Parallel, delayed, parallel_backend
from scipy.cluster import hierarchy


def read_file(filepath):
    f = open(filepath, 'r')
    lines = f.readlines()
    f.close()
    return lines


def write_file(filename, out_str):
    f = open(filename, 'w+')
    f.write(out_str)
    f.close()


def update_weighted_hamming_distance(assembly_1_idx, assembly_2_idx):
    col_1 = np.array(metablocks_pa_df[metablocks_pa_df.columns[assembly_1_idx+2]].astype(bool))
    col_2 = np.array(metablocks_pa_df[metablocks_pa_df.columns[assembly_2_idx+2]].astype(bool))
    pos = np.logical_xor(col_1, col_2)
    val = np.sum(metablocks_lengths[pos])

    col_1 = np.array(retained_blocks_pa_df[retained_blocks_pa_df.columns[assembly_1_idx+2]].astype(bool))
    col_2 = np.array(retained_blocks_pa_df[retained_blocks_pa_df.columns[assembly_2_idx+2]].astype(bool))
    pos = np.logical_xor(col_1, col_2)
    val += round(0.5*np.sum(retained_blocks_lengths[pos]))

    hamming_distance_matrix[assembly_1_idx][assembly_2_idx] += val
    hamming_distance_matrix[assembly_2_idx][assembly_1_idx] += val


def get_pricipal_components_count(model, min_explained_variance=0.8):
    cumsum = 0
    for idx, val in enumerate(model.explained_variance_ratio_):
        cumsum += val
        if(cumsum>=min_explained_variance):
            print(idx, cumsum)
            return idx
            break


def generate_hamming_distance_dendrogram(hamming_distance_mat, labels, save_plot=True, reference='', outdir='',
                                         show_plot=False, figsize=(12,60), leaf_font_size=9, logtransform=False, return_linkage=False):
    if(logtransform):
        hamming_distance_mat = np.log(hamming_distance_mat+np.ones(hamming_distance_mat.shape))
    Z = linkage(hamming_distance_mat)
    fig, ax = plt.subplots(figsize=figsize)
    dendrogram(Z, labels=labels, orientation='left', ax=ax, leaf_font_size=leaf_font_size)
    colored_labels = ax.get_ymajorticklabels()
    for label in colored_labels:
        text = label.get_text()
        if(':R' in text):
            color = 'red'
        elif(':S' in text):
            color = 'blue'
        else:
            color = 'black'
        label.set_color(color)
    plt.tight_layout()
    if(show_plot):
        plt.show()
    if(save_plot):
        outfile = '{}hamming_distance_{}'.format(outdir, len(labels))
        if(len(reference)>0):
            outfile += '_{}'.format(reference)
        if(logtransform):
            outfile += '_logtransformed'
        outfile += '.png'
        fig.savefig(outfile, dpi=fig.dpi, bbox_inches="tight")
    if(return_linkage):
        return Z


def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "{}:{:.2f}{}".format(leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if(len(newick) > 0):
            newick = "):{:.2f}{}".format(parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",{}".format(newick), node.dist, leaf_names)
        newick = "({}".format(newick)
        return newick


def str2bool(v):
    return v.lower() in ("y", "yes", "true", "t", "1")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Distance matrix computation for PRAWNS')
    parser.add_argument('-i', '--input', required=True, help="Input csv file")
    parser.add_argument('-d', '--prawns_dir', required=True, help="PRAWNS results directory (e.g. PRAWNS_results/)")
    parser.add_argument('-n', '--ncores', type=int, nargs='?', default=8, help="Number of cores to be used (default: 8)")
    parser.add_argument('-a', '--annotations', type=str2bool, nargs='?', default=False, const=True, help="Use given labels for NEWICK tree; " + 
                                        "if True, the last column in input csv should be contain the annotations for each sample (default: False)")

    args = parser.parse_args()

    ncores = args.ncores
    results_dir = args.prawns_dir

    input_pd = pd.read_csv(args.input, header=None)
    input_assemblies = input_pd[input_pd.columns[0]].values
    fasta_filepaths = input_pd[input_pd.columns[1]].values

    assembly_count = len(input_assemblies)

    metablocks_pa_df = pd.read_csv('{}/metablock_presence_absence.csv'.format(results_dir), header=None)
    retained_blocks_pa_df = pd.read_csv('{}/retained_block_presence_absence.csv'.format(results_dir), header=None)

    metablocks_lengths = np.array(metablocks_pa_df[metablocks_pa_df.columns[1]])
    retained_blocks_lengths = np.array(retained_blocks_pa_df[retained_blocks_pa_df.columns[1]])


    hamming_distance_matrix = np.ones((assembly_count, assembly_count))

    start = time.time()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(update_weighted_hamming_distance)(  assembly_1_idx, assembly_2_idx )
            for assembly_1_idx in range(1,assembly_count) for assembly_2_idx in range(assembly_1_idx))
    end = time.time() # timeit.timeit()
    print('TIME taken to generate distance matrix:', end-start)

    np.save('{}weighted_hamming_distance_matrix.npy'.format(results_dir), hamming_distance_matrix)

    log_hamming_distance_matrix = np.log(hamming_distance_matrix)
    np.save('{}log_weighted_hamming_distance_matrix.npy'.format(results_dir), log_hamming_distance_matrix)


    start = time.time()

    model_pca = PCA(n_components=log_hamming_distance_matrix.shape[1])
    X_dim_reduce = model_pca.fit_transform(log_hamming_distance_matrix)

    n_principal_components = get_pricipal_components_count(model_pca, 0.9)
    print(n_principal_components)

    previous_slope = np.inf
    previous_sum_of_sq = np.inf
    k = 0
    max_iter = math.floor(assembly_count/2)
    while(True and k<max_iter):
        k += 1
        km = KMeans(n_clusters=k)
        current_sum_sq_list = []
        for i in range(3):
            km = km.fit(X_dim_reduce[:,:n_principal_components])
            current_sum_sq_list.append(km.inertia_)
        current_sum_of_sq = np.mean(current_sum_sq_list)
        print('\t', k, current_sum_of_sq)
        if(k==1):
            previous_sum_of_sq = current_sum_of_sq
            continue
        else:
            current_slope = previous_sum_of_sq - current_sum_of_sq
            print(k, previous_slope, current_slope, previous_slope/current_slope)
            previous_sum_of_sq = current_sum_of_sq
            if(current_slope > 1.4*previous_slope): # filter anomalies in elbow point selection
                print(k, previous_slope, current_slope, previous_slope/current_slope, current_slope/previous_slope)
                continue
            if(previous_slope <= 1.4*current_slope):
                n_clusters = max(1, k-2)
                print()
                print(n_clusters)
                break
            previous_slope = current_slope

    kmeans_pca = KMeans(n_clusters=n_clusters)
    kmeans_pca = kmeans_pca.fit(X_dim_reduce[:,:n_principal_components])
    np.save('{}kmeans_labels.npy'.format(results_dir), kmeans_pca.labels_)
    print(Counter(kmeans_pca.labels_))
    end = time.time() # timeit.timeit()
    print('TIME taken to get genome groupings:', end-start)


    if(args.annotations):
        labels = input_pd[input_pd.columns[-1]].values
    else:
        labels = input_assemblies

    fig_ht = max(5, math.ceil(assembly_count/6.8))
    fig_wt = min( max(7,fig_ht), 12)
    tree_labels = list(labels + '    ' + list(map(str, np.array(kmeans_pca.labels_))))
    Z = generate_hamming_distance_dendrogram(   log_hamming_distance_matrix, tree_labels, save_plot=True, show_plot=False,
                                                figsize=(fig_wt,fig_ht), outdir=results_dir, logtransform=False, return_linkage=True)

    tree = hierarchy.to_tree(Z,False)
    f = getNewick(tree, "", tree.dist, tree_labels)
    write_file('{}phylogenetic_tree.newick'.format(results_dir), f)