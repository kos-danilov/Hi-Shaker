import os.path as op
import matplotlib.pyplot as plt
import numpy as np
import pandas
import h5py
import cooler
from scipy.sparse import coo_matrix
from scipy.sparse import bmat
import random
import os
import sys
from Bio import SeqIO

def read_hic(file):
    # reading data
    h5 = h5py.File(file, 'r')

    # save resolutions/x path for further adding of attributes
    reso_path = []

    def save_reso_path(path):
        if path.count("/") == 1:
            reso_path.append(path)

    h5.visit(save_reso_path)

    # save only paths with datasets (on depth 3)
    data_paths = []

    def save_paths(path):
        if path.count("/") == 3:
            data_paths.append(path)

    h5.visit(save_paths)

    # create corresponding arrays in memory
    #bins_chrom = h5.get(data_paths[0])[:]
    #bins_end = h5.get(data_paths[1])[:]
    #bins_start = h5.get(data_paths[2])[:]
    #bins_weight = h5.get(data_paths[3])[:]

    #chroms_length = h5.get(data_paths[4])[:]
    chroms_name = h5.get(data_paths[5])[:]

    #indexes_bin1_offset = h5.get(data_paths[6])[:]
    indexes_chrom_offset = h5.get(data_paths[7])[:]

    pixels_bin1_id = h5.get(data_paths[8])[:]
    pixels_bin2_id = h5.get(data_paths[9])[:]
    pixels_count = h5.get(data_paths[10])[:]

    #attr = h5[reso_path[0]].attrs
    
    return chroms_name, indexes_chrom_offset, pixels_bin1_id, pixels_bin2_id, pixels_count

# move contigs with coordinates changes
def move_contig(contigs, contig, before_contig):
    if contig > before_contig:
        contigs.insert(before_contig, contigs.pop(contig))
        
        contigs[before_contig]["curr_start"] = contigs[before_contig + 1]["curr_start"]
        contigs[before_contig]["curr_end"] = contigs[before_contig]["curr_start"] + contigs[before_contig]["length"]
        
        #changed_contigs = contigs[before_contig: contig]

        for i in range(before_contig + 1, contig + 1):
            contigs[i]["curr_start"] = contigs[i]["curr_start"] + contigs[before_contig]["length"]
            contigs[i]["curr_end"] = contigs[i]["curr_end"] + contigs[before_contig]["length"]

    else:
        contigs.insert(before_contig, contigs.pop(contig))
        
        
        if before_contig == len(contigs) - 1:
            # if we are moving contig to the end
            contigs[before_contig]["curr_end"] = contigs[before_contig - 1]["curr_end"]
            contigs[before_contig]["curr_start"] = contigs[before_contig]["curr_end"] - contigs[before_contig]["length"]
        else:
            # other case
            contigs[before_contig]["curr_end"] = contigs[before_contig + 1]["curr_start"]
            contigs[before_contig]["curr_start"] = contigs[before_contig]["curr_end"] - contigs[before_contig]["length"]
        
        if contig == 0:
            # for moving firts contigs
            contigs[0]["curr_start"] = 0
            contigs[0]["curr_end"] = contigs[0]["length"]
            
            for i in range(1, before_contig):
                contigs[i]["curr_start"] = contigs[i - 1]["curr_end"]
                contigs[i]["curr_end"] = contigs[i]["curr_start"] + contigs[i]["length"]
        else:
            # for any other contig 
            for i in range(contig, before_contig):
                contigs[i]["curr_start"] = contigs[i - 1]["curr_end"]
                contigs[i]["curr_end"] = contigs[i]["curr_start"] + contigs[i]["length"]

# change indexes to get inverse state
def inverse(arr, l):
    inv = np.array([], dtype = int)

    for i in arr:
        inv = np.append(inv, [l - i - 1])

    return inv

# change status of contig
def reverse_contig(actual_contigs, contig):
    if actual_contigs[contig]['orientation'] == 'forward':
        actual_contigs[contig]['orientation'] = 'reverse'
    else:
        actual_contigs[contig]['orientation'] = 'forward'

# find contigs within provided bin coordinates
def get_contig_subset(contigs, i1, i2, j1, j2):
    ith = []
    jth = []
    ith_ori = []
    jth_ori = []
    
    for i in range(len(contigs)):
        if contigs[i]["curr_end"] > i1 and contigs[i]["curr_start"] <= i2:
            ith.append(contigs[i]["contig"])
            ith_ori.append(i)
        if contigs[i]["curr_end"] > j1 and contigs[i]["curr_start"] <= j2:
            jth.append(contigs[i]["contig"])
            jth_ori.append(i)

    # return 1) contigs within coordinates with original names (indexes) 2) current locations of these contigs
    return ith, jth, ith_ori, jth_ori

# funciton for get T/F array of contigs. T - if no additional modification is nedeed from original Hi-C map
def get_changes_matrix(contig_subset):
    return np.reshape([x <= y for x in contig_subset[0] for y in contig_subset[1]], (len(contig_subset[0]), len(contig_subset[1])))

# return involved slices of sparse matrix by i and j interaction according to contig_subset with ZERO coordinates
def get_coo_matrix(actual_contigs, changes_matrix, contig_subset, i, j, pixels_bin1_id, pixels_bin2_id, pixels_count):
    if not changes_matrix[i][j]:
        subset = np.logical_and(np.logical_and(pixels_bin2_id >= actual_contigs[contig_subset[3][i]]["ori_start"], pixels_bin2_id < actual_contigs[contig_subset[3][i]]["ori_end"]), 
               np.logical_and(pixels_bin1_id >= actual_contigs[contig_subset[2][j]]["ori_start"], pixels_bin1_id < actual_contigs[contig_subset[2][j]]["ori_end"]))
        
        if actual_contigs[contig_subset[3][i]]["orientation"] == "reverse" and actual_contigs[contig_subset[2][j]]["orientation"] == "reverse":
            ss_bin1 = inverse(pixels_bin1_id[subset] - actual_contigs[contig_subset[2][j]]["ori_start"], actual_contigs[contig_subset[2][j]]["length"])
            ss_bin2 = inverse(pixels_bin2_id[subset] - actual_contigs[contig_subset[3][i]]["ori_start"], actual_contigs[contig_subset[3][i]]["length"])
            ss_count = pixels_count[subset]
            
            return ss_count, ss_bin2, ss_bin1
        
        elif actual_contigs[contig_subset[3][i]]["orientation"] == "reverse" and actual_contigs[contig_subset[2][j]]["orientation"] == "forward":
            ss_bin1 = pixels_bin1_id[subset] - actual_contigs[contig_subset[2][j]]["ori_start"]
            ss_bin2 = inverse(pixels_bin2_id[subset] - actual_contigs[contig_subset[3][i]]["ori_start"], actual_contigs[contig_subset[3][i]]["length"])
            ss_count = pixels_count[subset]
            
            return ss_count, ss_bin2, ss_bin1
        
        elif actual_contigs[contig_subset[3][i]]["orientation"] == "forward" and actual_contigs[contig_subset[2][j]]["orientation"] == "reverse":
            ss_bin1 = inverse(pixels_bin1_id[subset] - actual_contigs[contig_subset[2][j]]["ori_start"], actual_contigs[contig_subset[2][j]]["length"])
            ss_bin2 = pixels_bin2_id[subset] - actual_contigs[contig_subset[3][i]]["ori_start"]
            ss_count = pixels_count[subset]
            
            return ss_count, ss_bin2, ss_bin1
        
        else:
            ss_bin1 = pixels_bin1_id[subset] - actual_contigs[contig_subset[2][j]]["ori_start"]
            ss_bin2 = pixels_bin2_id[subset] - actual_contigs[contig_subset[3][i]]["ori_start"]
            ss_count = pixels_count[subset]

            return ss_count, ss_bin2, ss_bin1
        
    else:
        subset = np.logical_and(np.logical_and(pixels_bin2_id >= actual_contigs[contig_subset[3][j]]["ori_start"], pixels_bin2_id < actual_contigs[contig_subset[3][j]]["ori_end"]), 
               np.logical_and(pixels_bin1_id >= actual_contigs[contig_subset[2][i]]["ori_start"], pixels_bin1_id < actual_contigs[contig_subset[2][i]]["ori_end"]))
        
        if actual_contigs[contig_subset[3][i]]["orientation"] == "reverse" and actual_contigs[contig_subset[2][j]]["orientation"] == "reverse":
            ss_bin1 = inverse(pixels_bin1_id[subset] - actual_contigs[contig_subset[2][i]]["ori_start"], actual_contigs[contig_subset[2][i]]["length"])
            ss_bin2 = inverse(pixels_bin2_id[subset] - actual_contigs[contig_subset[3][j]]["ori_start"], actual_contigs[contig_subset[3][j]]["length"])
            ss_count = pixels_count[subset]
            
            if contig_subset[0][i] == contig_subset[1][j]:
                return np.append(ss_count, ss_count), np.append(ss_bin1, ss_bin2), np.append(ss_bin2, ss_bin1)
            else:
                return ss_count, ss_bin1, ss_bin2
        
        elif actual_contigs[contig_subset[3][i]]["orientation"] == "reverse" and actual_contigs[contig_subset[2][j]]["orientation"] == "forward":
            ss_bin1 = inverse(pixels_bin1_id[subset] - actual_contigs[contig_subset[2][i]]["ori_start"], actual_contigs[contig_subset[2][i]]["length"])
            ss_bin2 = pixels_bin2_id[subset] - actual_contigs[contig_subset[3][j]]["ori_start"]

            ss_count = pixels_count[subset]
            
            if contig_subset[0][i] == contig_subset[1][j]:
                return np.append(ss_count, ss_count), np.append(ss_bin1, ss_bin2), np.append(ss_bin2, ss_bin1)
            else:
                return ss_count, ss_bin1, ss_bin2
        
        elif actual_contigs[contig_subset[3][i]]["orientation"] == "forward" and actual_contigs[contig_subset[2][j]]["orientation"] == "reverse":
            ss_bin1 = pixels_bin1_id[subset] - actual_contigs[contig_subset[2][i]]["ori_start"]
            ss_bin2 = inverse(pixels_bin2_id[subset] - actual_contigs[contig_subset[3][j]]["ori_start"], actual_contigs[contig_subset[3][j]]["length"])
            ss_count = pixels_count[subset]
            
            if contig_subset[0][i] == contig_subset[1][j]:
                return np.append(ss_count, ss_count), np.append(ss_bin1, ss_bin2), np.append(ss_bin2, ss_bin1)
            else:
                return ss_count, ss_bin1, ss_bin2
        
        else:
            ss_bin1 = pixels_bin1_id[subset] - actual_contigs[contig_subset[2][i]]["ori_start"]
            ss_bin2 = pixels_bin2_id[subset] - actual_contigs[contig_subset[3][j]]["ori_start"]
            ss_count = pixels_count[subset]

            if contig_subset[0][i] == contig_subset[1][j]:
                return np.append(ss_count, ss_count), np.append(ss_bin1, ss_bin2), np.append(ss_bin2, ss_bin1)
            else:
                return ss_count, ss_bin1, ss_bin2

        #ss_bin1 = pixels_bin1_id[subset] - actual_contigs[contig_subset[2][i]]["ori_start"]
        #ss_bin2 = pixels_bin2_id[subset] - actual_contigs[contig_subset[3][j]]["ori_start"]
        #ss_count = pixels_count[subset]
        
        #if contig_subset[0][i] == contig_subset[1][j]:
        #    return np.append(ss_count, ss_count), np.append(ss_bin1, ss_bin2), np.append(ss_bin2, ss_bin1)
        #else:
        #    return ss_count, ss_bin1, ss_bin2

# return dense matrix
def get_dense_matrix(actual_contigs, i1, i2, j1, j2, pixels_bin1_id, pixels_bin2_id, pixels_count):
    contig_subset = get_contig_subset(actual_contigs, i1, i2, j1, j2)
    changes_matrix = get_changes_matrix(contig_subset)
    
    d = np.array([0], dtype = "int32")
    r = np.array([0], dtype = "int32")
    c = np.array([0], dtype = "int32")
    
    shift_y = 0

    for i in range(len(contig_subset[0])):
        shift_x = 0
    
        for j in range(len(contig_subset[1])):
            try:
                curr = get_coo_matrix(actual_contigs, changes_matrix, contig_subset, i, j, pixels_bin1_id, pixels_bin2_id, pixels_count)
                d = np.append(d, curr[0])
                r = np.append(r, curr[1] + shift_y)
                c = np.append(c, curr[2] + shift_x)

                ###!!!
                shift_x += actual_contigs[contig_subset[3][j]]["length"]
            except:
                pass

        ###!!!
        shift_y += actual_contigs[contig_subset[2][i]]["length"]
    
    return coo_matrix((d, (r, c))).toarray()

def save_fasta(original_fasta, edited_fasta, actual_contigs, chroms_name):
    if not os.path.exists(original_fasta):
        return "No original fasta found"
    else:
        my_seq = SeqIO.to_dict(SeqIO.parse(original_fasta, "fasta"))
        
        my_file = open(edited_fasta + ".fasta", "a")

        for i in actual_contigs:
            if i['orientation'] == "reverse":
                reverse = my_seq[str(chroms_name[i['contig']])].reverse_complement()
                reverse.id = my_seq[str(chroms_name[i['contig']])].id
                reverse.description = my_seq[str(chroms_name[i['contig']])].description

                SeqIO.write(reverse, my_file, "fasta")
            else:
                SeqIO.write(my_seq[str(chroms_name[i['contig']])], my_file, "fasta")

        my_file.close()

def save_png(name, i1, i2, j1, j2, actual_contigs, pixels_bin1_id, pixels_bin2_id, pixels_count):
    fig = plt.figure(figsize=(25, 25))
    ax = fig.add_subplot(111)
    im = ax.matshow(np.log10(get_dense_matrix(actual_contigs, i1, i2, j1, j2, pixels_bin1_id, pixels_bin2_id, pixels_count)), cmap='YlOrRd')
    fig.colorbar(im)
    fig.savefig(name + '.png')

def show_matrix(i1, i2, j1, j2, actual_contigs, pixels_bin1_id, pixels_bin2_id, pixels_count):
    fig = plt.figure(figsize=(25, 25))
    ax = fig.add_subplot(111)
    im = ax.matshow(np.log10(get_dense_matrix(actual_contigs, i1, i2, j1, j2, pixels_bin1_id, pixels_bin2_id, pixels_count)), cmap='YlOrRd')
    fig.colorbar(im)
