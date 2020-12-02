import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
import glob
import os, time
from collections import Counter
import math
import shlex, subprocess
from io import StringIO
from queue import Queue
from multiprocessing import Pool
import random
from scipy.spatial.distance import hamming
import pickle
import timeit
import time
from datetime import datetime
import argparse, sys
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA, TruncatedSVD
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform, hamming
import psutil
import scipy
from scipy.sparse import csr_matrix, dok_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from joblib import Parallel, delayed, parallel_backend


def read_file(filepath):
    f = open(filepath, 'r')
    lines = f.readlines()
    f.close()
    return lines


def write_file(filename, out_str):
    f = open(filename, 'w+')
    f.write(out_str)
    f.close()


def run_cpp_binaries(binary_file, *args):
    cmd = binary_file
    for arg in args:
        cmd += ' {}'.format(arg)
    print('run_cpp_binaries: ', cmd)
    try:
        p = subprocess.Popen(shlex.split(cmd))#, stdout=open(mum_results_file, 'w'))
        # print('waiting...', cmd)
        p.wait()
        print(cmd, ' output obtained')
    except Exception as e:
        print(cmd, e)


def create_dir(directory_name):
    cmd = "mkdir {}".format(directory_name)
    p = subprocess.Popen(shlex.split(cmd))
    p.wait()
    print(directory_name, " created")


def remove_file(filename):
    cmd = "rm -rf {}".format(filename)
    try:
        p = subprocess.Popen(shlex.split(cmd), shell=True)
        p.wait()
    except Exception as e:
        print("remove_file error", cmd, e)


def remove_file2(filename):
    cmd = "rm -rf {}".format(filename)
    try:
        p = subprocess.Popen(shlex.split(cmd))
        p.wait()
    except Exception as e:
        print("remove_file error", cmd, e)


def cpp_dir_input_setup(assemblies, fasta_file_list, ncores, block_pair_binned_partitions, min_presence_count,
                        use_oriented_links, oriented_links_file_list):
    outdir = "PRAWNS_results"
    if(os.path.exists(outdir)):
        outdir += "_{}".format(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    outdir += '/'
    create_dir(outdir)
    create_dir("{}binned_kmers/".format(outdir))
    create_dir("{}retained_binned_kmers_{}/".format(outdir, min_presence_count))
    create_dir("{}retained_binned_kmers_{}/assemblywise/".format(outdir, min_presence_count))
    create_dir("{}contig_lengths/".format(outdir))
    create_dir("{}kmer_pairs_{}/".format(outdir, min_presence_count))
    create_dir("{}grouped_pairs/".format(outdir))
    create_dir("{}collinear_blocks/".format(outdir))
    create_dir("{}oll/".format(outdir))
    create_dir("{}assemblywise_blocks/".format(outdir))
    create_dir("{}neighbour_pairs/".format(outdir))
    create_dir("{}k_neighbours/".format(outdir))
    create_dir("{}metablocks/".format(outdir))
    create_dir("{}structural_variants/".format(outdir))
    create_dir("{}paired_variants/".format(outdir))


    assembly_count = len(assemblies)
    step = math.ceil(assembly_count/(float)(ncores))
    assembly_partition_pos_arr = np.arange(step, assembly_count, step)
    core_idx = 0
    partition_start = 0

    for partition_end in assembly_partition_pos_arr:
        opstr = '\n'.join(fasta_file_list[partition_start:partition_end])+'\n'
        write_file("{}input_binned_assemblies_{}.txt".format(outdir, core_idx), opstr)
        core_idx += 1
        partition_start = partition_end
    if(partition_start<assembly_count):
        opstr = '\n'.join(fasta_file_list[partition_start:assembly_count])+'\n'
        write_file("{}input_binned_assemblies_{}.txt".format(outdir, core_idx), opstr)

    opstr = '\n'.join(fasta_file_list) + '\n'
    write_file("{}all_assembly_filepaths.txt".format(outdir), opstr)


    block_pair_step = math.ceil(assembly_count/(float)(block_pair_binned_partitions))
    block_pair_partition_pos_list = [0]
    core_idx = 0
    partition_start = 0
    partition_end = block_pair_step
    
    while(partition_end < assembly_count):
        if(use_oriented_links):
            opstr = '\n'.join(oriented_links_file_list[partition_start:partition_end])+'\n'
            write_file("{}ol_files_{}.txt".format(outdir, core_idx), opstr)
            core_idx += 1
        block_pair_partition_pos_list.append(partition_end)
        partition_start = partition_end
        partition_end += block_pair_step
    # if(partition_end!=assembly_count):
    if(use_oriented_links):
        opstr = '\n'.join(oriented_links_file_list[partition_start:assembly_count])+'\n'
        write_file("{}ol_files_{}.txt".format(outdir, core_idx), opstr)
    block_pair_partition_pos_list.append(assembly_count)

    # opstr = '\n'.join(assemblies)+'\n'
    # write_file("{}input_assemblies.txt".format(outdir), opstr)
    return outdir, np.concatenate([[0],assembly_partition_pos_arr]), np.array(block_pair_partition_pos_list)


def make_oriented_links_file_groups(oriented_links_file_list, partitions, assembly_count, outdir):
    step = math.ceil(assembly_count/(float)(partitions))
    # ol_partition_pos_arr = np.arange(step, assembly_count, step)
    ol_partition_pos_list = []
    group_idx = 0
    partition_start = 0

    while(partition_start < assembly_count):
        opstr = '\n'.join(oriented_links_file_list[partition_start:min((partition_start+step),assembly_count)])+'\n'
        write_file("{}ol_files_{}.txt".format(outdir, group_idx), opstr)
        group_idx += 1
        ol_partition_pos_list.append(partition_start)
        partition_start += step
    # if((partition_start+1) < assembly_count):
    #     opstr = '\n'.join(oriented_links_file_list[partition_start:assembly_count])+'\n'
    #     write_file("{}ol_files_{}.txt".format(outdir, group_idx), opstr)
    ol_partition_pos_list.append(assembly_count)

    return np.array(ol_partition_pos_list)


def presence_vector_files(hash_string, hash_string_idx, in_dir, outdir='.'):
    cmd = "ls {}{}* > {}{}_group_files".format(in_dir, hash_string, outdir, hash_string_idx)
    print("presence_vector_files : ", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("presence_vector_files created:", hash_string)
    except Exception as e:
        print("presence_vector_files :", cmd, e)


def combine_filemaps(in_dir, outdir='.'):
    cmd = "echo {}filemap_* | xargs cat > {}combined_map.csv 2>/dev/null".format(in_dir, outdir)
    print("combine_filemaps : ", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("combine_filemaps combined")
    except Exception as e:
        print("combine_filemaps error :", cmd, e)


def get_intermediate_filemaps(filenames_str, outfilename):
    cmd = "cat {} > {} 2>/dev/null".format(filenames_str, outfilename)
    rm_cmd = "rm -f {}".format(filenames_str)
#     print("get_intermediate_filemaps : ", cmd)
    print("get_intermediate_filemaps : ", outfilename)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("get_intermediate_filemaps combined")
        p = subprocess.Popen(rm_cmd, shell=True)
        p.wait()
    except Exception as e:
        print("get_intermediate_filemaps error :", cmd, e)


def combine_intermediate_filemaps(in_dir, outdir='.'):
    cmd = "echo {}intermediate_combined_* | xargs cat > {}final_combined_map.csv 2>/dev/null".format(in_dir, outdir)
    rm_cmd = "echo {}intermediate_combined_* | xargs rm -rf".format(in_dir)
#     print("combine_filemaps : ", cmd)
    print("combine_filemaps : ", outdir)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("combine_intermediate_filemaps combined")
        p = subprocess.Popen(rm_cmd, shell=True)
        p.wait()
    except Exception as e:
        print("combine_intermediate_filemaps error :", cmd, e)


def vertical_combine_filemaps(filenames_str, outfilename, merge_separator=','):
    cmd = "echo {}| xargs paste -d{} > {}".format(filenames_str, merge_separator, outfilename)
    print("vertical_combine_filemaps : ", cmd)
    # rm_cmd = "rm -f {}".format(filenames_str)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("vertical_combine_filemaps combined")
        # p = subprocess.Popen(rm_cmd, shell=True)
        # p.wait()
    except Exception as e:
        print("vertical_combine_filemaps error :", cmd, e)


def sv_cov_combine_header(input_file, cov_file, outfilename, merge_separator=','):
    cmd = "cut -d, -f1 {} | paste -d, - {} | sed '1i Sample_name,SV_genome_coverage(%),SV_count' > {}".format(input_file, cov_file, outfilename)
    rm_cmd = "rm -rf {}".format(cov_file)
    print("sv_cov_combine_header : ", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("sv_cov_combine_header combined")
        p = subprocess.Popen(rm_cmd, shell=True)
        p.wait()
    except Exception as e:
        print("sv_cov_combine_header error :", cmd, e)


def get_hamming_distance(df, assembly_count, start_idx=1, end_idx=-1):
    nrows = df.shape[0]
    hamming_dist_mat = np.zeros((assembly_count, assembly_count))
    
    for i in range(start_idx-1, end_idx):
        for j in range(start_idx, assembly_count):
            pos = np.logical_xor(df[df.columns[i]], df[df.columns[j]])
            val = np.sum(df.blocks_count[pos])
            hamming_dist_mat[i][j] = val
            hamming_dist_mat[j][i] = val
    print("get_hamming_distance", start_idx, end_idx)
    # print(hamming_dist_mat)
    return hamming_dist_mat


def get_weighted_hamming_distance(df, assembly_count, start_idx=1, end_idx=-1):
    nrows = df.shape[0]
    hamming_dist_mat = np.zeros((assembly_count, assembly_count))
    
    for i in range(start_idx-1, end_idx):
        for j in range(start_idx, assembly_count):
            pos = np.logical_xor(df[df.columns[i]], df[df.columns[j]])
            val = np.sum(df.total_blocks_length[pos])
            hamming_dist_mat[i][j] = val
            hamming_dist_mat[j][i] = val
    print("get_weighted_hamming_distance", start_idx, end_idx)
    # print(hamming_dist_mat)
    return hamming_dist_mat


def generate_column_vector(filename, col_idx, outfilename):
    cmd = "cut -d, -f{} {} > {}".format(col_idx, filename, outfilename)
    print("generate_column_vector started:", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("generate_column_vector ended:", col_idx)
    except Exception as e:
        print("generate_column_vector error:", cmd, e)


def generate_start_indices(filename, col_idx, outfilename):
    cmd = "cut -d, -f{} {} | awk '{{total += $0; $0 = total - $0}}1' > {}".format(col_idx, filename, outfilename)
    print("generate_start_indices started:", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("generate_start_indices ended:", col_idx)
    except Exception as e:
        print("generate_start_indices error:", cmd, e)


def get_total_block_count(filename, col_idx):
    cmd = "cut -d, -f{} {} | awk '{{ sum += $1 }} END {{ print sum }}'".format(col_idx, filename)
    print("get_total_block_count: ", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, error = p.communicate()
        print(out, error)
        return out.decode("utf-8").strip()
    except Exception as e:
        print("get_total_block_count ERROR: ", cmd, e)


def retain_presence_vectors(filename, col_idx_str, outfilename):
    cmd = "cut -d, -f{} --complement {} > {}".format(col_idx_str, filename, outfilename)
    print("retain_presence_vectors started:", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("retain_presence_vectors ended:", col_idx_str)
    except Exception as e:
        print("retain_presence_vectors error:", cmd, e)


def get_pricipal_components_count(model, min_explained_variance=0.8):
    cumsum = 0
    for idx, val in enumerate(model.explained_variance_ratio_):
        cumsum += val
        if(cumsum>=min_explained_variance):
            print(idx, cumsum)
            return idx
            break


def get_kmeans_elbow_point(dim_reduced_data, linearity_limiting_threshold=1.1):
    previous_slope = np.inf
    previous_sum_of_sq = np.inf
    k = 0
    while(True):
        k += 1
        print(k)
        km = KMeans(n_clusters=k)
        km = km.fit(X)
        current_sum_of_sq = km.inertia_
        if(k==1):
            previous_sum_of_sq = current_sum_of_sq
            continue
        else:
            current_slope = previous_sum_of_sq - current_sum_of_sq
            if(previous_slope <= linearity_limiting_threshold*current_slope):
                return k-2
            previous_slope = current_slope
            previous_sum_of_sq = current_sum_of_sq


def kmeans_wrapper(X, dim_reduce_model, min_explained_variance=0.8, linearity_limiting_threshold=1.1):
    
    dim_reduce_model.fit(X)
    
    n_principal_components = get_pricipal_components_count(dim_reduce_model, min_explained_variance)
    X_dim_reduce = dim_reduce_model.transform(X)
    n_clusters = get_kmeans_elbow_point(X_dim_reduce[:,:n_principal_components], linearity_limiting_threshold)
    
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans = kmeans.fit(X_dim_reduce[:,:n_principal_components])
    
    return kmeans, n_clusters


def generate_hamming_distance_dendrogram(hamming_distance_mat, labels, save_plot=True, reference='', outdir='',
                                         show_plot=False, figsize=(12,60), leaf_font_size=9, logtransform=False):
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


def str2bool(v):
    return v.lower() in ("y", "yes", "true", "t", "1")


def load_distance_edge_weighted_sparse_connectivity_submatrix(knn_subgraph_dir, partition_idx, current_partition_block_start,
                                                            total_blocks_count, blocks_per_partition):
    knn_subgraph_file = knn_subgraph_dir + str(partition_idx)
    lines = read_file(knn_subgraph_file)
    current_block_count = blocks_per_partition
    if(current_partition_block_start >= total_blocks_count):
        current_block_count = 0
    elif(current_partition_block_start+blocks_per_partition > total_blocks_count):
        current_block_count = total_blocks_count - current_partition_block_start
    
    print("load_distance_edge_weighted_sparse_connectivity_submatrix: ",
          knn_subgraph_file, len(lines), current_block_count)
    
    connectivity_submatrix = dok_matrix((current_block_count, total_blocks_count), dtype=np.uint16)
    if(len(lines)==0):
        return partition_idx, connectivity_submatrix
    
#     block_idx_offset = int(lines[0].split(',')[0])
    print("\t\t",lines[0])
    print("\t\t", lines[-1])
    for line in lines:
        lsplit = line.split(',')
        block_1_idx = int(lsplit[0])-current_partition_block_start #block_idx_offset
        loc = 1
        while(loc < len(lsplit)):
            neighbour_block = int(lsplit[loc])
            edge_weight = int(lsplit[loc+1])
            connectivity_submatrix[block_1_idx, neighbour_block] = edge_weight
            loc += 2
            if(len(connectivity_submatrix.keys())%1000 == 0):
                print('\t\t', knn_subgraph_file, len(connectivity_submatrix.keys()))
    print('\t', knn_subgraph_file, len(connectivity_submatrix.keys()))
    return partition_idx, connectivity_submatrix


def run_bfs_in_parallel(start_candidates):
    for node in start_candidates:
        if(node in remaining_nodes):
            current_queue = [node]
            current_component_nodes = set()
            while(len(current_queue)>0):
                current_component_nodes.add(current_queue[0])
                edges = mst_birectional[current_queue[0], :]
                for e in edges.indices:
                    if(e not in current_component_nodes):
                        current_queue.append(e)
                current_queue.pop(0)
            
            if(len(current_component_nodes)>0):
                min_node = min(current_component_nodes)
                if(min_node not in active_nodes):
                    active_nodes.add(min_node)
#                     lock.acquire()
                    remaining_nodes.difference_update(current_component_nodes)
                    all_connected_components.append(list(current_component_nodes))
#                     active_nodes.remove(min_node)
#                     lock.release()


def save_component_groups(comp_start_idx, comp_end_idx, filename):
    outstring = "{}\n".format(comp_start_idx)
    for idx in range(comp_start_idx, comp_end_idx):
        outstring += ",".join(np.sort(all_connected_components[idx]).astype(str)) + "\n"
        outstring += ",".join(all_connected_components_presence_mat[idx].astype(np.uint8).astype(str)) + "\n"
    write_file(filename, outstring)


if __name__ == "__main__":

    # parser = argparse.ArgumentParser(description='PRAWNS: Pan-genome RepresentAtion of Whole geNomeS tool')
    parser = argparse.ArgumentParser(description='PRAWNS: Pan-genome representation of whole genomes tool')
    parser.add_argument('-i', '--input', required=True, help="Input csv file")
    parser.add_argument('-n', '--ncores', type=int, nargs='?', default=8, help="Number of cores to be used (default: 8)")
    parser.add_argument('-K', '--kmer_len', type=int, nargs='?', default=25, help="Length of kmers (default: 25)")
    parser.add_argument('-p', '--min_perc', type=float, nargs='?', default=5.0, help="Minimum %% of genomes a variant would be present in (default: 5.0)")
    parser.add_argument('-l', '--use_oriented_links', type=str2bool, nargs='?', default=False, const=True, help="Use MetaCarvel oriented links; " + 
                                        "if True, 3rd column in input csv should be path to oriented links of corresponding assembly (default: False)")
    parser.add_argument('-b', '--min_group_blocks', type=int, nargs='?', default=3, help="Minimum number of exact matching regions (blocks) that " +
                                        "can be grouped into metablocks across the genomes (default: 3)")
    parser.add_argument('-M', '--max_metablock_mismatch', type=int, nargs='?', default=25, help="Maximum number of mismatches permitted to allow " +
                                        "merger and extension of metablocks across the genomes (default: 25)")
    parser.add_argument('-s', '--min_block_size', type=int, nargs='?', default=50, help="Smallest size of a block that is to be retained as a " +
                                        "structural variant (default: 50)")
    # parser.add_argument('-R', '--max_pairing_range', type=int, nargs='?', default=100, help="Maximum number of bases between the structural variants " +
    #                                     "from a genome for paired analysis (default: 100)")
    parser.add_argument('-S', '--max_intervariant_separation', type=int, nargs='?', default=50, help="Maximum number of bases between the structural variants " +
                                        "from a genome for paired analysis (default: 50)")
    parser.add_argument('-m', '--mem', nargs='?', default="36000MB", help="Upper limit for RAM memory usage.  " +
                                        "Can be in mb/MB/gb/GB/tb/TB (case insensitive), default unit is MB. (default: 36000MB)")
    parser.add_argument('-g', '--genome_len', nargs='?', default="4M", help="Average genome length.  " +
                                        "Can be in k/K/m/M/g/G (case insensitive), default unit is M, i.e. 1x10^6 nt. (default: 4M)")

    args = parser.parse_args()

# def main(args):

    ncores = args.ncores #8 # 5

    input_pd = pd.read_csv(args.input, header=None)
    input_assemblies = input_pd[input_pd.columns[0]].values
    fasta_filepaths = input_pd[input_pd.columns[1]].values


    kmer_len = args.kmer_len
    prefix_len = 5 ## Fixed kmer prefix length
    min_presence_perc = args.min_perc
    prefix_count = pow(4, prefix_len)
    min_presence_count = math.ceil(len(input_assemblies)*min_presence_perc/100.0)
    use_oriented_links = args.use_oriented_links #True
    assembly_count = len(input_assemblies)
    # feature_partitions = max(ncores, math.ceil(assembly_count/10))

    genome_len = args.genome_len
    if(genome_len.isdigit()):
        genome_len_mb = int(genome_len)
    else:
        val = genome_len[:-1]
        units = genome_len[-1]
        assert(val.isdigit())
        val = int(val)
        assert(units in set(['k', 'K', 'm', 'M', 'g', 'G']) )
        # assert(units[0].isalpha() and units[1].isalpha())
        units = units.upper()
        if(units=='M'):
            genome_len_mb = val
        elif(units=='G'):
            genome_len_mb = val*1000
        else: # K
            genome_len_mb = val/1000

    mem = args.mem
    if(mem.isdigit()):
        available_memory = int(mem)
    else:
        assert(len(mem)>2)
        val = mem[:-2]
        units = mem[-2:]
        assert(val.isdigit())
        val = int(val)
        assert(units[0] in set(['m', 'M', 'g', 'G', 't', 'T']) and units[1] in set(['b', 'B']))
        # assert(units[0].isalpha() and units[1].isalpha())
        units = units.upper()
        if(units[0]=='M'):
            available_memory = val
        elif(units[0]=='G'):
            available_memory = val*1000
        else:
            available_memory = val*1000000

    # max_assemblies_per_partition = available_memory/ncores - 2500
    # if(max_assemblies_per_partition < 10):
    #     max_assemblies_per_partition = 10
    # else:
    #     max_assemblies_per_partition = math.floor(math.sqrt(max_assemblies_per_partition/genome_len_mb))
    #     max_assemblies_per_partition = max(10, max_assemblies_per_partition)

    # max_assemblies_per_partition = (available_memory/ncores - 2250)
    max_assemblies_per_partition = (available_memory/ncores - 1024)
    # max_assemblies_per_partition = (available_memory/ncores - 512)
    if(max_assemblies_per_partition < 2):
        max_assemblies_per_partition = 2
    else:
        # max_assemblies_per_partition = math.floor(math.sqrt(max_assemblies_per_partition/genome_len_mb))
        # max_assemblies_per_partition = math.floor(math.sqrt(max_assemblies_per_partition/(4*genome_len_mb)))
        # max_assemblies_per_partition = math.floor(max_assemblies_per_partition/16*genome_len_mb)
        max_assemblies_per_partition = math.floor(max_assemblies_per_partition/(20*genome_len_mb))
        max_assemblies_per_partition = max(2, max_assemblies_per_partition)

    feature_partitions = max(ncores, math.ceil(assembly_count/max_assemblies_per_partition))

    assemblies_per_partition = math.ceil(assembly_count/feature_partitions)
    min_component_blocks = args.min_group_blocks
    # max_inter_block_pair_separation = args.max_pairing_range
    max_inter_block_pair_separation = args.max_intervariant_separation
    min_block_size = args.min_block_size
    max_metablock_mismatch = args.max_metablock_mismatch
    k_neighbours = 4
    max_neighbour_separation = 5
    # block_pair_partitions = max(feature_partitions, math.ceil(assembly_count/4.0))
    # block_pair_partitions = max(feature_partitions, math.ceil(assembly_count/ncores))
    block_pair_partitions = feature_partitions
    seq2write_batch = 1000

    if(use_oriented_links):
        oriented_links_paths = input_pd[input_pd.columns[2]].values
        uol=1
    else:
        uol=0
        oriented_links_paths = []

    print(ncores, input_pd.shape, kmer_len, min_presence_perc, min_presence_count, use_oriented_links, assembly_count, feature_partitions,
            block_pair_partitions, available_memory, max_inter_block_pair_separation, k_neighbours, max_neighbour_separation)
            

    results_dir, binned_assembly_arr, block_pair_partition_pos_arr = cpp_dir_input_setup( input_assemblies, fasta_filepaths, ncores,
                                                                                            feature_partitions, min_presence_count, use_oriented_links,
                                                                                            oriented_links_paths)
    print(binned_assembly_arr)
    print(block_pair_partition_pos_arr)

    out_str = "#genomes: {}\n".format(len(input_assemblies))
    out_str += "kmer_len: {}\n".format(kmer_len)
    out_str += "min_perc: {}\n".format(min_presence_perc)
    out_str += "#cores: {}\n".format(ncores)
    out_str += "min_presence_count: {}\n".format(min_presence_count)
    out_str += "use_oriented_links: {}\n".format(use_oriented_links)
    out_str += "available_memory: {}\n".format(available_memory)
    out_str += "min_group_blocks: {}\n".format(min_component_blocks)
    out_str += "max_pairing_range: {}\n".format(max_inter_block_pair_separation)
    out_str += "max_metablock_mismatch: {}\n".format(max_metablock_mismatch)
    out_str += "min_block_size: {}\n".format(min_block_size)
    out_str += "k_neighbours: {}\n".format(k_neighbours)
    out_str += "max_neighbour_separation: {}\n".format(max_neighbour_separation)
    out_str += "genome_len: {}\n".format(genome_len)


    pool = Pool(processes=ncores)
    start = time.time() # timeit.timeit()
    # for core_idx in range(min(ncores, len(binned_assembly_arr))):
    for core_idx in range(len(binned_assembly_arr)):
        pool.apply_async(run_cpp_binaries, args=("./kmer_positional_binning.o", "{}input_binned_assemblies_{}.txt".format(results_dir, core_idx),
                                                kmer_len, "{}binned_kmers/".format(results_dir), binned_assembly_arr[core_idx],
                                                "{}contig_lengths/".format(results_dir)))
    pool.close()
    pool.join()
    end = time.time() # timeit.timeit()
    print('TIME taken to bin kmers from all files:', end-start)


    # start = time.time()
    # Parallel(n_jobs=ncores, prefer="threads")(
    #     delayed(run_cpp_binaries)(  "./kmer_pair_positional_binning.o", "{}input_binned_assemblies_{}.txt".format(results_dir, core_idx),
    #                                 kmer_len, "{}kmer_pairs_{}/".format(results_dir, min_presence_count), binned_assembly_arr[core_idx],
    #                                 "{}contig_lengths/".format(results_dir), prefix_count, feature_partitions)
    #         for core_idx in range(len(binned_assembly_arr)))
    # end = time.time() # timeit.timeit()
    # print('TIME taken to generate binned kmer pairs:', end-start)


    pool = Pool(processes=ncores)
    start = time.time() # timeit.timeit()
    for prefix_idx in range(prefix_count):
        pool.apply_async(run_cpp_binaries, args=("./kmer_filtering.o", assembly_count, prefix_idx, "{}binned_kmers/".format(results_dir),
                                                min_presence_count, "{}retained_binned_kmers_{}/".format(results_dir, min_presence_count)))
    pool.close()
    pool.join()
    end = time.time() # timeit.timeit()
    print('TIME taken to get kmers present in at least {} perc of input files:'.format(min_presence_perc), end-start)


    pool = Pool(processes=ncores)
    start = time.time() # timeit.timeit()
    for assembly_idx in range(assembly_count):
        pool.apply_async(run_cpp_binaries, args=("./kmer_pair_generator.o", assembly_idx, prefix_count,
                                                "{}binned_kmers/".format(results_dir),
                                                "{}retained_binned_kmers_{}/".format(results_dir, min_presence_count),
                                                "{}retained_binned_kmers_{}/assemblywise/".format(results_dir, min_presence_count),
                                                "{}kmer_pairs_{}/".format(results_dir, min_presence_count), feature_partitions))
    pool.close()
    pool.join()
    end = time.time() # timeit.timeit()
    print('TIME taken to generate kmer pairs files:', end-start)


    pool = Pool(processes=ncores)
    start = time.time() # timeit.timeit()
    for partition_idx in range(feature_partitions): # higher partition indices tend to have higher number of kmer pairs mapped
        pool.apply_async(run_cpp_binaries, args=("./kmer_pair_grouping.o", feature_partitions - 1 - partition_idx, 
                                                "{}kmer_pairs_{}/".format(results_dir, min_presence_count), assembly_count,
                                                "{}grouped_pairs/".format(results_dir), min_presence_count, kmer_len))
    pool.close()
    pool.join()
    end = time.time() # timeit.timeit()
    print('TIME taken to generate grouped kmer pairs per partition:', end-start)


    # start = time.time()
    # offset = max(math.ceil(min(feature_partitions/ncores, available_memory/10*ncores)), ncores)
    # offset_range = math.ceil(feature_partitions/offset)
    # Parallel(n_jobs=ncores, prefer="threads")(
    #     delayed(run_cpp_binaries)(  "./kmer_pair_grouper.o", feature_partitions, offset_idx, offset, 
    #                                 "{}kmer_pairs_{}/".format(results_dir, min_presence_count), assembly_count,
    #                                 "{}grouped_pairs/".format(results_dir), min_presence_count, kmer_len)
    #         for offset_idx in range(offset_range))
    # end = time.time() # timeit.timeit()
    # print('TIME taken to generate grouped kmer pairs per partition:', end-start)


    hash_strings = np.unique([filename.split('_')[-1] for filename in glob.glob("{}grouped_pairs/group_files_*".format(results_dir))])

    pool = Pool(processes=ncores)
    start = time.time() # timeit.timeit()
    for hash_idx, hash_string in enumerate(hash_strings):
        pool.apply_async(run_cpp_binaries, args=("./kmer_block_constructor.o", hash_string,
                                                "{}grouped_pairs/group_files_{}".format(results_dir, hash_string),
                                                assembly_count, kmer_len, "{}collinear_blocks/".format(results_dir)))
    pool.close()
    pool.join()
    end = time.time() # timeit.timeit()
    print('TIME taken to merge kmer pairs to form blocks:', end-start)


    pool = Pool(processes=ncores)
    start = time.time() # timeit.timeit()
    part_size = min(1250, len(hash_strings)/ncores) ## TO DO: Choice for optimal number of files to concat simultaneously
    hash_partitions = np.arange(0, len(hash_strings), part_size)
    hash_partitions = [math.floor(t) for t in hash_partitions]
    hash_partitions.append(len(hash_strings))
    hash_partitions = np.unique(hash_partitions)
    print(hash_partitions)
    for idx in range(1,len(hash_partitions)):
        instr = ' '.join(['{}collinear_blocks/filemap_{}'.format(results_dir, hash_str)
                              for hash_str in hash_strings[hash_partitions[idx-1]:hash_partitions[idx]]])
        outfile = '{}collinear_blocks/intermediate_combined_{}'.format(results_dir, idx)
        pool.apply_async(get_intermediate_filemaps, args=(instr, outfile))
    pool.close()
    pool.join()
    combine_intermediate_filemaps('{}collinear_blocks/'.format(results_dir), '{}collinear_blocks/'.format(results_dir))
    end = time.time() # timeit.timeit()
    print('TIME taken to merge filemaps:', end-start)


    filemap_combined_filename = '{}collinear_blocks/final_combined_map.csv'.format(results_dir)
    pool = Pool(processes=ncores)
    for assembly_idx in range(assembly_count):
        pool.apply_async(generate_column_vector, args=(filemap_combined_filename, assembly_idx+1, "{}col_{}.txt".format(results_dir, assembly_idx)))
    pool.close()
    pool.join()

    paths_filename = '{}merged_filemap_paths.txt'.format(results_dir)
    start_block_indices_filename = '{}merged_filemap_start_block_indices.txt'.format(results_dir)
    if ncores>1:
        w = 2
    else:
        w = 1
    pool = Pool(processes=w)
    pool.apply_async(generate_column_vector, args=(filemap_combined_filename, assembly_count+1, paths_filename))
    pool.apply_async(generate_start_indices, args=(filemap_combined_filename, assembly_count+2, start_block_indices_filename))
    pool.close()
    pool.join()
    
    total_blocks_count = int(get_total_block_count(filemap_combined_filename, assembly_count+2))
    # block_pair_partitions = max(block_pair_partitions, math.ceil(total_blocks_count*(2 + np.log2(ncores))/available_memory)) #/5000))
    length_based_adjustment = 1
    if(genome_len_mb > 4):
        length_based_adjustment = genome_len_mb/4
    block_pair_partitions = max(block_pair_partitions, math.ceil(total_blocks_count*(1.5 + np.log2(ncores))*length_based_adjustment/available_memory)) #/5000))
    # block_pair_partitions = max(block_pair_partitions, math.ceil((total_blocks_count*4*ncores)/(available_memory*max_assemblies_per_partition))) #/5000))
    # block_pair_partitions = max(block_pair_partitions, math.ceil((total_blocks_count*8*ncores)/(available_memory*100))) #/5000))
    # block_pair_partitions = max(block_pair_partitions, math.ceil((total_blocks_count*ncores)/(available_memory*4))) #/5000))
    # block_pair_partitions = max(block_pair_partitions, math.ceil((total_blocks_count*ncores*assemblies_per_partition)/(available_memory*64))) #/5000))
    block_pair_partitions = math.ceil(block_pair_partitions/ncores)*ncores
    print("Total individual blocks count:", total_blocks_count)
    print("Modified block_pair_partitions count:", block_pair_partitions)
    print("Available Memory:", available_memory)

    out_str += "\n\nTotal individual block count: {}\n".format(total_blocks_count)

    # max_pairs_before_dump_limit = math.floor((available_memory - total_blocks_count*ncores*1000)/(1000))
    # print("max_pairs_before_dump_limit Memory:", max_pairs_before_dump_limit)
    # max_pairs_before_dump_limit = math.floor((36000000000 - total_blocks_count*ncores*1000)/(1000))
    # Another approximation: if memory is given in MB ==> max pairs = (mem*1000*1000 - total_blocks_count*ncores*1000)/(ncores*100)

    # max_pairs_before_dump_limit = math.floor((available_memory-1500)*8000/ncores) - math.ceil(total_blocks_count*10*(1 + np.log10(ncores)))
    # print("max_pairs_before_dump_limit Memory:", max_pairs_before_dump_limit)


    pool = Pool(processes=ncores)
    start = time.time() # timeit.timeit()
    for idx in range(1, len(block_pair_partition_pos_arr)):
        if(use_oriented_links):
            pool.apply_async(run_cpp_binaries, args=("./kmer_block_and_neighbours_locator.o", "{}col_".format(results_dir), block_pair_partition_pos_arr[idx-1],
                                                block_pair_partition_pos_arr[idx]-1, paths_filename, start_block_indices_filename, feature_partitions,
                                                block_pair_partitions, total_blocks_count, "{}kmer_pairs_{}/".format(results_dir, min_presence_count),
                                                "{}assemblywise_blocks/".format(results_dir), kmer_len, "{}contig_lengths/".format(results_dir), k_neighbours-1,
                                                max_neighbour_separation, "{}neighbour_pairs/".format(results_dir), min_block_size, uol,
                                                "{}ol_files_{}.txt".format(results_dir, idx-1), "{}oll/".format(results_dir)))
        else:
            pool.apply_async(run_cpp_binaries, args=("./kmer_block_and_neighbours_locator.o", "{}col_".format(results_dir), block_pair_partition_pos_arr[idx-1],
                                                block_pair_partition_pos_arr[idx]-1, paths_filename, start_block_indices_filename, feature_partitions,
                                                block_pair_partitions, total_blocks_count, "{}kmer_pairs_{}/".format(results_dir, min_presence_count),
                                                "{}assemblywise_blocks/".format(results_dir), kmer_len, "{}contig_lengths/".format(results_dir), k_neighbours-1,
                                                max_neighbour_separation, "{}neighbour_pairs/".format(results_dir), min_block_size, uol))
    pool.close()
    pool.join()
    end = time.time() # timeit.timeit()
    print('TIME taken to locate neighbour pairs per partition:', end-start)

    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(run_cpp_binaries)("./remove_files.o", paths_filename, core_idx, ncores)
        for core_idx in range(ncores))
        
    end = time.time() # timeit.timeit()
    print('TIME taken to remove collinear block files:', end-start)


    pool = Pool(processes=ncores)
    start = time.time() # timeit.timeit()
    for partition_idx in range(block_pair_partitions):
        pool.apply_async(run_cpp_binaries, args=("./kmer_k_neigbour_subgraph_generator.o", k_neighbours, assembly_count, partition_idx,
                                                min_presence_count, "{}neighbour_pairs/".format(results_dir), "{}k_neighbours/".format(results_dir)))
    pool.close()
    pool.join()
    end = time.time() # timeit.timeit()
    print('TIME taken to k nearest neighbour subgraphs:', end-start)


    blocks_per_partition = math.ceil(total_blocks_count/block_pair_partitions)
    part_arr = np.arange(0, total_blocks_count, blocks_per_partition)
    pool = Pool(processes=ncores)
    start = time.time() # timeit.timeit()
    subgraphs_results = []
    for partition_idx, current_partition_block_start in enumerate(part_arr):
        subgraphs_results.append(pool.apply_async(load_distance_edge_weighted_sparse_connectivity_submatrix,
                                                  args=("{}k_neighbours/".format(results_dir), partition_idx,
                                                        current_partition_block_start, total_blocks_count,
                                                        blocks_per_partition)))
    pool.close()
    pool.join()
    end = time.time() # timeit.timeit()
    print('TIME taken to load k-nearest connectivity submatrices:', end-start)

    stack_order = np.argsort([mat.get()[0] for mat in subgraphs_results])
    print(stack_order)

    connectivity_matrix = subgraphs_results[stack_order[0]].get()[1]
    for idx in stack_order[1:]:
        connectivity_matrix = scipy.sparse.vstack([connectivity_matrix, subgraphs_results[idx].get()[1]])
        print(idx, subgraphs_results[stack_order[idx]].get()[1].shape, connectivity_matrix.shape)

    connectivity_mat_csr = connectivity_matrix.tocsr()
    # Asymetry removal eliminates the cases where block 1 is a neighbor for block 2 but not vice versa
    # Particularly required in cases such as:
    # Assemblies having: -------[REGION_1]--[REGION_2]-------
    # Assemblies having: -------[REGION_1]-[REGION_3]-[REGION_2]-------
    # Blocks at the ends of region 3 should not be considered to be neighbors of blocks from region 1 and 2
    asymetry_remover = (connectivity_mat_csr!=0).astype(np.uint16)
    connectivity_mat_csr_2 = connectivity_mat_csr.multiply(asymetry_remover.T)

    print("Connectivity Matrix:", connectivity_mat_csr.nnz, connectivity_mat_csr_2.nnz)

    mst = minimum_spanning_tree(connectivity_mat_csr_2)
    mst_birectional = mst + mst.T
    mst_birectional = mst_birectional.astype(np.uint16)

    print("MST:", mst.nnz, mst_birectional.nnz)

    start = time.time()
    singleton_check = mst_birectional.sum(axis=0)
    all_start_candidates = np.where((singleton_check!=0).tolist()[0])[0]
    batch_size = math.ceil(len(all_start_candidates)/ncores)
    print(len(all_start_candidates), batch_size)
    print(np.sum((singleton_check==0).tolist()[0]))
    start_candidates_arr = np.arange(0, len(all_start_candidates), batch_size)
    start_candidates_arr = np.concatenate([start_candidates_arr, [len(all_start_candidates)]])
    print(start_candidates_arr)

    active_nodes = set()
    remaining_nodes = set(all_start_candidates)
    all_connected_components = []

    # TODO: Can have collisions; resolve the issue of same component being traversed over separate cores
    Parallel(n_jobs=ncores, prefer="threads", require='sharedmem')(
        delayed(run_bfs_in_parallel)(
            all_start_candidates[start_candidates_arr[core_idx]:start_candidates_arr[core_idx+1]])
        for core_idx in range(ncores))
        
    end = time.time() # timeit.timeit()
    print('TIME taken to fetch all connected components:', end-start)

    out_str += "\n\nTotal connected component count: {}\n".format(len(all_connected_components))

    print(len(all_connected_components), len(active_nodes))
    with open('{}all_connected_components.pkl'.format(results_dir), 'wb') as fp:
            pickle.dump(all_connected_components, fp)

    component_lengths = np.array([len(comp) for comp in all_connected_components])
    print(len(component_lengths))
    print(np.mean(component_lengths), np.median(component_lengths), np.std(component_lengths))
    print(np.min(component_lengths), np.max(component_lengths))
    print(Counter(component_lengths))
    print()
    print(np.sum(component_lengths>=10), np.sum(component_lengths>=100))
    print()

    individual_blocks_expanded = []
    lines = read_file('{}collinear_blocks/final_combined_map.csv'.format(results_dir))
    for line in lines:
        lsplit = line.split(',')
        current_presence = np.array([int(val) for val in lsplit[:-3]])
        occurrence = int(lsplit[-2])
        while(occurrence>0):
            individual_blocks_expanded.append(current_presence)
            # if(len(individual_blocks_expanded)%10000 == 0):
            #     print(len(individual_blocks_expanded))
            occurrence -= 1
    individual_blocks_expanded = np.array(individual_blocks_expanded, dtype=np.int8)
    print(individual_blocks_expanded.shape)

    start = time.time()
    all_connected_components_presence_mat = []
    min_component_blocks_presence_ratio = 0.75
    for comp in all_connected_components:
        comp_len = len(comp)
        counts = np.sum(individual_blocks_expanded[np.array(comp)], axis=0)
        valid_check_pos = np.logical_and(counts>0, counts<comp_len)
        if(np.sum(valid_check_pos)>0):
            quantile = np.percentile(counts[valid_check_pos], 90)
            threshold = round(0.85*comp_len)
            if(quantile < threshold and quantile > min_component_blocks_presence_ratio*comp_len):
                threshold = round(quantile)
        else:
            threshold = comp_len
        all_connected_components_presence_mat.append(counts>=threshold)
    all_connected_components_presence_mat = np.array(all_connected_components_presence_mat)
    print(all_connected_components_presence_mat.shape)
    np.save('{}all_connected_components_presence_mat.npy'.format(results_dir), all_connected_components_presence_mat)
    # with open('{}all_connected_components_presence_mat.pkl'.format(results_dir), 'wb') as fp:
    #         pickle.dump(all_connected_components_presence_mat, fp)
    end = time.time() # timeit.timeit()
    print('TIME taken to assign assemblywise presence of all connected components:', end-start)

    # trimmed for components containing at least a certain number of blocks per component (min_group_blocks)
    print("min_component_blocks: ", min_component_blocks)
    all_connected_components = [component for component in all_connected_components if len(component)>=min_component_blocks]
    print(len(all_connected_components))
    all_connected_components_presence_mat = all_connected_components_presence_mat[component_lengths>=min_component_blocks]
    print(all_connected_components_presence_mat.shape)

    out_str += "\n\nConnected components retained count: {}\n".format(len(all_connected_components))
    write_file('{}params.txt'.format(results_dir), out_str)


    group_count = max(ncores, math.ceil(len(all_connected_components)*2/available_memory))
    print(group_count)
    grouped_components_count = math.ceil(len(all_connected_components)/group_count)
    max_component_group_size = min(2000, math.floor(available_memory*1.5/ncores))
    while(grouped_components_count > max_component_group_size):
        group_count += ncores
        grouped_components_count = math.ceil(len(all_connected_components)/group_count)
    print(group_count, grouped_components_count)

    comp_idx_arr = np.arange(0, all_connected_components_presence_mat.shape[0], grouped_components_count)
    comp_idx_arr = np.concatenate([comp_idx_arr, [all_connected_components_presence_mat.shape[0]]])
    print(comp_idx_arr)

    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(save_component_groups)( comp_idx_arr[idx-1], comp_idx_arr[idx],
                                        '{}comp_{}_trimmed_new_group_{}.txt'.format(results_dir, min_component_blocks, idx-1))
    for idx in range(1, len(comp_idx_arr)))


    start = time.time()
    if(use_oriented_links):
        Parallel(n_jobs=ncores, prefer="threads")(
            delayed(run_cpp_binaries)("./kmer_metablock_constructor.o",
                                    '{}comp_{}_trimmed_new_group_{}.txt'.format(results_dir, min_component_blocks, idx-1),
                                    "{}assemblywise_blocks/".format(results_dir), max_metablock_mismatch, "{}contig_lengths/".format(results_dir),
                                    assembly_count, math.ceil(assembly_count/feature_partitions), "{}metablocks/".format(results_dir), uol,
                                    "{}oll/".format(results_dir))
            for idx in range(1, len(comp_idx_arr)))
    else:
        Parallel(n_jobs=ncores, prefer="threads")(
            delayed(run_cpp_binaries)("./kmer_metablock_constructor.o",
                                    '{}comp_{}_trimmed_new_group_{}.txt'.format(results_dir, min_component_blocks, idx-1),
                                    "{}assemblywise_blocks/".format(results_dir), max_metablock_mismatch, "{}contig_lengths/".format(results_dir),
                                    assembly_count, math.ceil(assembly_count/feature_partitions), "{}metablocks/".format(results_dir), uol)
            for idx in range(1, len(comp_idx_arr)))
    end = time.time() # timeit.timeit()
    print('TIME taken to generate metablocks:', end-start)


    start = time.time()
    if(use_oriented_links):
        Parallel(n_jobs=ncores, prefer="threads")(
            delayed(run_cpp_binaries)("./kmer_feature_filtering_and_pairing.o", len(comp_idx_arr)-1, "{}metablocks/".format(results_dir),
                                    min_block_size, "{}assemblywise_blocks/".format(results_dir), idx-1, block_pair_partition_pos_arr[idx-1],
                                    block_pair_partition_pos_arr[idx]-1, max_inter_block_pair_separation, "{}contig_lengths/".format(results_dir),
                                    block_pair_partitions, total_blocks_count, len(all_connected_components), "{}structural_variants/".format(results_dir),
                                    "{}paired_variants/".format(results_dir), uol, "{}oll/".format(results_dir))
            for idx in range(1, len(block_pair_partition_pos_arr)))
    else:
        Parallel(n_jobs=ncores, prefer="threads")(
            delayed(run_cpp_binaries)("./kmer_feature_filtering_and_pairing.o", len(comp_idx_arr)-1, "{}metablocks/".format(results_dir),
                                    min_block_size, "{}assemblywise_blocks/".format(results_dir), idx-1, block_pair_partition_pos_arr[idx-1],
                                    block_pair_partition_pos_arr[idx]-1, max_inter_block_pair_separation, "{}contig_lengths/".format(results_dir),
                                    block_pair_partitions, total_blocks_count, len(all_connected_components), "{}structural_variants/".format(results_dir),
                                    "{}paired_variants/".format(results_dir), uol)
            for idx in range(1, len(block_pair_partition_pos_arr)))
    end = time.time() # timeit.timeit()
    print('TIME taken to locate structural variants:', end-start)


    # start = time.time()
    # Parallel(n_jobs=ncores, prefer="threads")(
    #     delayed(run_cpp_binaries)(  "./kmer_feature_aggregation_and_filtering.o", partition_idx, assembly_count, assemblies_per_partition,
    #                                 "{}structural_variants/".format(results_dir), "{}paired_variants/".format(results_dir), min_presence_count)
    #         for partition_idx in range(block_pair_partitions))
    # end = time.time() # timeit.timeit()
    # print('TIME taken to aggregate partitioned variants matrices:', end-start)


    start = time.time()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(run_cpp_binaries)(  "./kmer_feature_aggregation_and_pair_filtering.o", partition_idx, assembly_count, assemblies_per_partition,
                                    "{}structural_variants/".format(results_dir), "{}paired_variants/".format(results_dir), min_presence_count)
            for partition_idx in range(block_pair_partitions))
    end = time.time() # timeit.timeit()
    print('TIME taken to aggregate partitioned variants matrices:', end-start)


    # start = time.time()
    # Parallel(n_jobs=ncores, prefer="threads")(
    #     delayed(run_cpp_binaries)(  "./kmer_feature_aggregation_and_pair_filtering.o", partition_idx, assembly_count, assemblies_per_partition,
    #                                 "{}structural_variants/".format(results_dir), "{}paired_variants/".format(results_dir), min_presence_count,
    #                                 max_inter_block_pair_separation)
    #         for partition_idx in range(block_pair_partitions))
    # end = time.time() # timeit.timeit()
    # print('TIME taken to aggregate partitioned variants matrices:', end-start)


    start = time.time()
    if(use_oriented_links):
        Parallel(n_jobs=ncores, prefer="threads")(
            delayed(run_cpp_binaries)(  "./kmer_distant_feature_pair_updating.o", block_pair_partition_pos_arr[idx-1], block_pair_partition_pos_arr[idx]-1,
                                        block_pair_partitions, max_inter_block_pair_separation, "{}structural_variants/".format(results_dir),
                                        "{}paired_variants/".format(results_dir), "{}contig_lengths/".format(results_dir), uol, "{}oll/".format(results_dir))
            for idx in range(1, len(block_pair_partition_pos_arr)))
    else:
        Parallel(n_jobs=ncores, prefer="threads")(
            delayed(run_cpp_binaries)(  "./kmer_distant_feature_pair_updating.o", block_pair_partition_pos_arr[idx-1], block_pair_partition_pos_arr[idx]-1,
                                        block_pair_partitions, max_inter_block_pair_separation, "{}structural_variants/".format(results_dir),
                                        "{}paired_variants/".format(results_dir), "{}contig_lengths/".format(results_dir), uol)
            for idx in range(1, len(block_pair_partition_pos_arr)))
    end = time.time() # timeit.timeit()
    print('TIME taken to update partitioned paired variants:', end-start)


    start = time.time()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(vertical_combine_filemaps)( ' '.join([ "{}paired_variants/{}_{}_{}_updated_separation".format(results_dir,
                                                            block_pair_partition_pos_arr[idx-1], block_pair_partition_pos_arr[idx]-1, partition_idx)
                                                    for idx in range(1, len(block_pair_partition_pos_arr)) ]),
                                            '{}paired_variants/intrapair_separation_{}'.format(results_dir, partition_idx))
        for partition_idx in range(block_pair_partitions))

    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(vertical_combine_filemaps)( ' '.join([ "{}paired_variants/{}_{}_{}_updated_presence".format(results_dir,
                                                            block_pair_partition_pos_arr[idx-1], block_pair_partition_pos_arr[idx]-1, partition_idx)
                                                    for idx in range(1, len(block_pair_partition_pos_arr)) ]),
                                            '{}paired_variants/pair_presence_absence_{}'.format(results_dir, partition_idx))
        for partition_idx in range(block_pair_partitions))
    end = time.time() # timeit.timeit()
    print('TIME taken to merge partitioned paired variant matrices:', end-start)


    start = time.time() # timeit.timeit()
    part_size = math.ceil(block_pair_partitions/ncores)
    concat_groups = list(np.arange(0, block_pair_partitions, part_size))
    if(concat_groups[-1] != block_pair_partitions):
        concat_groups.append(block_pair_partitions)
    print(concat_groups)
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(get_intermediate_filemaps)( ' '.join(['{}structural_variants/filtered_block_idx_{}'.format(results_dir, partition_idx)
                                                    for partition_idx in range(concat_groups[group_idx-1], concat_groups[group_idx])]),
                                            '{}structural_variants/merged_block_idx_{}'.format(results_dir, group_idx))
            for group_idx in range(1, len(concat_groups)))
    instr = ' '.join(['{}structural_variants/merged_block_idx_{}'.format(results_dir, group_idx)
                              for group_idx in range(1, len(concat_groups))])
    get_intermediate_filemaps(instr, '{}structural_variants/merged_block_indices'.format(results_dir))
    end = time.time() # timeit.timeit()
    print('TIME taken to merge filtered block indices:', end-start)


    start = time.time() # timeit.timeit()

    instr = ' '.join(['{}contig_lengths/cov_{}'.format(results_dir, idx-1) for idx in range(1, len(block_pair_partition_pos_arr))])
    get_intermediate_filemaps(instr, '{}sv_coverage.txt'.format(results_dir))
    sv_cov_combine_header(args.input, '{}sv_coverage.txt'.format(results_dir), '{}sv_coverages.txt'.format(results_dir))


    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(get_intermediate_filemaps)( ' '.join(['{}structural_variants/metablock_coords_{}'.format(results_dir, partition_idx)
                                                    for partition_idx in range(concat_groups[group_idx-1], concat_groups[group_idx])]),
                                            '{}structural_variants/merged_metablock_coords_{}'.format(results_dir, group_idx))
            for group_idx in range(1, len(concat_groups)))

    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(get_intermediate_filemaps)( ' '.join(['{}structural_variants/metablock_presence_absence_{}'.format(results_dir, partition_idx)
                                                    for partition_idx in range(concat_groups[group_idx-1], concat_groups[group_idx])]),
                                            '{}structural_variants/merged_metablock_presence_absence_{}'.format(results_dir, group_idx))
            for group_idx in range(1, len(concat_groups)))

    instr_list = [' '.join(['{}structural_variants/merged_metablock_coords_{}'.format(results_dir, group_idx)
                                  for group_idx in range(1, len(concat_groups))])]
    outstr_list = ['{}metablock_coords.csv'.format(results_dir)]

    instr_list.append(' '.join(['{}structural_variants/merged_metablock_presence_absence_{}'.format(results_dir, group_idx)
                              for group_idx in range(1, len(concat_groups))]))
    outstr_list.append('{}metablock_presence_absence.csv'.format(results_dir))

    Parallel(n_jobs=w, prefer="threads")(
        delayed(get_intermediate_filemaps)(instr_list[idx], outstr_list[idx]) for idx in range(w))

    end = time.time() # timeit.timeit()
    print('TIME taken to merge metablock coords and presence:', end-start)


    start = time.time() # timeit.timeit()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(get_intermediate_filemaps)( ' '.join(['{}structural_variants/block_coords_{}'.format(results_dir, partition_idx)
                                                    for partition_idx in range(concat_groups[group_idx-1], concat_groups[group_idx])]),
                                            '{}structural_variants/merged_block_coords_{}'.format(results_dir, group_idx))
            for group_idx in range(1, len(concat_groups)))

    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(get_intermediate_filemaps)( ' '.join(['{}structural_variants/block_presence_absence_{}'.format(results_dir, partition_idx)
                                                    for partition_idx in range(concat_groups[group_idx-1], concat_groups[group_idx])]),
                                            '{}structural_variants/merged_block_presence_absence_{}'.format(results_dir, group_idx))
            for group_idx in range(1, len(concat_groups)))

    instr_list = [' '.join(['{}structural_variants/merged_block_coords_{}'.format(results_dir, group_idx)
                                  for group_idx in range(1, len(concat_groups))])]
    outstr_list = ['{}retained_block_coords.csv'.format(results_dir)]

    instr_list.append(' '.join(['{}structural_variants/merged_block_presence_absence_{}'.format(results_dir, group_idx)
                              for group_idx in range(1, len(concat_groups))]))
    outstr_list.append('{}retained_block_presence_absence.csv'.format(results_dir))

    Parallel(n_jobs=w, prefer="threads")(
        delayed(get_intermediate_filemaps)(instr_list[idx], outstr_list[idx]) for idx in range(w))

    end = time.time() # timeit.timeit()
    print('TIME taken to merge filtered block coords and presence:', end-start)


    start = time.time() # timeit.timeit()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(get_intermediate_filemaps)( ' '.join(['{}paired_variants/intrapair_separation_{}'.format(results_dir, partition_idx)
                                                    for partition_idx in range(concat_groups[group_idx-1], concat_groups[group_idx])]),
                                            '{}paired_variants/merged_intrapair_separation_{}'.format(results_dir, group_idx))
            for group_idx in range(1, len(concat_groups)))

    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(get_intermediate_filemaps)( ' '.join(['{}paired_variants/pair_presence_absence_{}'.format(results_dir, partition_idx)
                                                    for partition_idx in range(concat_groups[group_idx-1], concat_groups[group_idx])]),
                                            '{}paired_variants/merged_pair_presence_absence_{}'.format(results_dir, group_idx))
            for group_idx in range(1, len(concat_groups)))

    instr_list = [' '.join(['{}paired_variants/merged_intrapair_separation_{}'.format(results_dir, group_idx)
                                  for group_idx in range(1, len(concat_groups))])]
    outstr_list = ['{}intrapair_separation.csv'.format(results_dir)]

    instr_list.append(' '.join(['{}paired_variants/merged_pair_presence_absence_{}'.format(results_dir, group_idx)
                              for group_idx in range(1, len(concat_groups))]))
    outstr_list.append('{}pair_presence_absence.csv'.format(results_dir))

    Parallel(n_jobs=w, prefer="threads")(
        delayed(get_intermediate_filemaps)(instr_list[idx], outstr_list[idx]) for idx in range(w))
    
    end = time.time() # timeit.timeit()
    print('TIME taken to merge filtered block coords and presence:', end-start)


    start = time.time()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(run_cpp_binaries)(  "./kmer_feature_fasta_generator.o", "{}all_assembly_filepaths.txt".format(results_dir), seq2write_batch,
                                    "{}metablocks/".format(results_dir), "{}structural_variants/merged_block_indices".format(results_dir),
                                    "{}assemblywise_blocks/".format(results_dir), "{}structural_variants/".format(results_dir), core_idx, ncores)
            for core_idx in range(ncores))

    instr_list = [  ' '.join(['{}structural_variants/metablocks_{}.fasta'.format(results_dir, core_idx) for core_idx in range(ncores)]),
                    ' '.join(['{}structural_variants/blocks_{}.fasta'.format(results_dir, core_idx) for core_idx in range(ncores)])]
    outstr_list = ['metablocks.fasta', 'retained_blocks.fasta']

    Parallel(n_jobs=w, prefer="threads")(
        delayed(get_intermediate_filemaps)(instr_list[idx], '{}{}'.format(results_dir, outstr_list[idx]))
            for idx in range(w))
    end = time.time() # timeit.timeit()
    print('TIME taken to generate variant fasta sequences:', end-start)


    start = time.time()
    clean_up_list  = [  'assemblywise_blocks/', 'binned_kmers/', 'collinear_blocks/', 'grouped_pairs/',
                        'kmer_pairs_{}/'.format(min_presence_count), 'k_neighbours/', 'merged_filemap_paths.txt',
                        'merged_filemap_start_block_indices.txt', 'metablocks/', 'neighbour_pairs/', 'oll', 'paired_variants/',
                        'retained_binned_kmers_{}/'.format(min_presence_count), 'structural_variants/']

    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(remove_file2)('{}{}'.format(results_dir, clean_up)) for clean_up in clean_up_list)
    end = time.time() # timeit.timeit()
    print('TIME taken for clean-up:', end-start)



# if __name__ == "__main__":

#     # parser = argparse.ArgumentParser(description='PRAWNS: Pan-genome RepresentAtion of Whole geNomeS tool')
#     parser = argparse.ArgumentParser(description='PRAWNS: Pan-genome representation of whole genomes tool')
#     parser.add_argument('-i', '--input', required=True, help="Input csv file")
#     parser.add_argument('-n', '--ncores', type=int, nargs='?', default=8, help="Number of cores to be used (default: 8)")
#     parser.add_argument('-K', '--kmer_len', type=int, nargs='?', default=25, help="Length of kmers (default: 25)")
#     parser.add_argument('-p', '--min_perc', type=float, nargs='?', default=5.0, help="Minimum %% of genomes a variant would be present in (default: 5.0)")
#     parser.add_argument('-l', '--use_oriented_links', type=str2bool, nargs='?', default=False, const=True, help="Use MetaCarvel oriented links; " + 
#                                         "if True, 3rd column in input csv should be path to oriented links of corresponding assembly (default: False)")
#     parser.add_argument('-b', '--min_group_blocks', type=int, nargs='?', default=3, help="Minimum number of exact matching regions (blocks) that " +
#                                         "can be grouped into metablocks across the genomes (default: 3)")
#     parser.add_argument('-M', '--max_metablock_mismatch', type=int, nargs='?', default=25, help="Maximum number of mismatches permitted to allow " +
#                                         "merger and extension of metablocks across the genomes (default: 25)")
#     parser.add_argument('-s', '--min_block_size', type=int, nargs='?', default=50, help="Smallest size of a block that is to be retained as a " +
#                                         "structural variant (default: 50)")
#     parser.add_argument('-R', '--max_pairing_range', type=int, nargs='?', default=100, help="Maximum number of bases between the structural variants " +
#                                         "from a genome for paired analysis (default: 100)")
#     parser.add_argument('-m', '--mem', nargs='?', default="36000MB", help="Upper limit for RAM memory usage.  " +
#                                         "Can be in mb/MB/gb/GB/tb/TB (case insensitive), default unit is MB. (default: 36000MB)")
#     parser.add_argument('-g', '--genome_len', nargs='?', default="4M", help="Average genome length.  " +
#                                         "Can be in k/K/m/M/g/G (case insensitive), default unit is M, i.e. 1x10^6 nt. (default: 4M)")

#     args = parser.parse_args()

#     main(args)