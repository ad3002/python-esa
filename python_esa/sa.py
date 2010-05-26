#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 11.04.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

import os
import pickle

from trseeker.tools.seqfile import gen_block_sequences 
from trseeker.tools.sequence_tools import clear_sequence

def fasta_to_sa_input(fasta_file, sa_input_file, index_file_name=None, file_name=None, start_id=None):
    
    if not index_file_name:
        index_file_name = fasta_file + ".pos_index"
    
    if os.path.isfile(sa_input_file):
        os.unlink(sa_input_file)
    if os.path.isfile(index_file_name):
        os.unlink(index_file_name)
    
    with open(fasta_file, "rb") as fasta_fh:
        with open(index_file_name, "a") as index_fw:
            with open(sa_input_file, "a") as fa_fw:
                if start_id:
                    i = start_id
                else:
                    i = 0
                for head, sequence, head_start, next_head in gen_block_sequences(">", fasta_fh):
                    if i>0:
                        fa_fw.write("$")
                    i += 1
                    sequence = clear_sequence(sequence)
                    fa_fw.write(sequence)
                    
                    if file_name:
                        head_info = "%s:%s" % (file_name, head.strip())
                    else:
                        head_info = head.strip()
                    index_fw.write("%s\t%s\t%s\n" % (i, head_info, len(sequence) ))
    return i

def clean_suffix_array_data(sa_file):
    '''
    Смысл - убрать из lcp суффиксу содержащие
    терминальный символ и неизвестные нуклеотиды
    @param sa_file:
    '''
    
    temp_sa_file = sa_file + ".tmp"
    with open(sa_file, "r") as fh:
        with open(temp_sa_file, "a") as fw:
            for line in fh:
                if line.startswith("type"):
                    continue
                data = line.strip().split("\t")
                
                suffix = data[7]
                
                if "n" in suffix or "N" in suffix or "$" in suffix:
                    continue
                fw.write(line)
    
    os.remove(sa_file)
    os.rename(temp_sa_file, sa_file)
    
def translate_docid_to_trfid(contr_sa_doc_index, doc_to_trf_file, trf_to_doc_file):
    
    doc_to_trf = {}
    trf_to_doc = {}
    
    with open(contr_sa_doc_index) as fh:
        
        for line in fh:
            data = line.strip().split("\t")
            
            docid = int(data[0]) - 1
            
            id = data[1].split(">")[-1]
            try:
                trfid = int(id)
            except:
                trfid = docid + 1
            
            doc_to_trf[docid] = trfid
            trf_to_doc[trfid] = docid
    
    with open(doc_to_trf_file, "w") as dffh:
        with open(trf_to_doc_file, "w") as tdfh:
            pickle.dump(doc_to_trf, dffh)
            pickle.dump(trf_to_doc, tdfh)
            
def compute_doc_suffix_tf(sa_file, trfid_suffix_tf_file, doc_tf_filedoc_to_trf_file):
    
    cutoff = 0
    length_cutoff = 0
    
    with open(doc_tf_filedoc_to_trf_file) as fh:
        docid_to_trfid = pickle.load(fh)
    
    result = {}
    all_ids = set()
    with open(sa_file, "rb") as fh:
        with open(trfid_suffix_tf_file, "a") as fw:
            for line in fh:
                data = line.strip().split("\t")
                
                suffix = data[7]
                length = int(data[6])
                tf = int(data[4])
                df = int(data[5])
                ids = data[9].split()
                
                for id in ids:
                    all_ids.add(id)
                    result.setdefault(suffix, {}).setdefault(id, 0)
                    result[suffix][id] += 1
                    
            for id in all_ids:
                for suffix in result:
                    if id in result[suffix]:
                        tf = result[suffix][id]
                        if tf > cutoff and len(suffix)>length_cutoff:
                            fw.write("%s\t%s\t%s\t%s\n" % (docid_to_trfid[int(id)],
                                                           suffix, 
                                                           tf, 
                                                           df) )
    
def reduce_sa(sa_file, cont_sa_reduced,  params=None):
    
    
    tf_cutoff = 2
    df_cutoff = 2
    length_cutoff = 15
    
    if params:
        tf_cutoff = params["tf_cutoff"]
        df_cutoff = params["df_cutoff"]
        length_cutoff = params["length_cutoff"]
        
    
    with open(sa_file, "r") as fh:
        with open(cont_sa_reduced, "a") as fw:
            for line in fh:
                if line.startswith("type"):
                    continue
                data = line.strip().split("\t")
                
                suffix = data[7]
                tf = int(data[4])
                df = int(data[5])
                length = int(data[6])
                
                if tf >= tf_cutoff and df >= df_cutoff and length >= length_cutoff:
                    fw.write(line)
    
def split_sa(sa_file, cont_sa_reduced, other_sa_file,  params=None):
    
    
    tf_cutoff = 2
    df_cutoff = 1
    length_cutoff = 15
    
    if params:
        tf_cutoff = params["tf_cutoff"]
        df_cutoff = params["df_cutoff"]
        length_cutoff = params["length_cutoff"]
        up_length_cutoff = params["up_length_cutoff"]
        
    
    with open(sa_file, "r") as fh:
        with open(cont_sa_reduced, "a") as fw:
            with open(other_sa_file, "a") as fw2:
                for line in fh:
                    if line.startswith("type"):
                        continue
                    data = line.strip().split("\t")
                    
                    suffix = data[7]
                    tf = int(data[4])
                    df = int(data[5])
                    length = int(data[6])
                    
                    if tf >= tf_cutoff and df >= df_cutoff and length >= length_cutoff and length <= up_length_cutoff:
                        fw.write(line)
                    else:
                        fw2.write(line)
    
def compute_doc_suffix_stat(sa_file, sa_counts_file, verbose=True):
    
    counts = {}
    
    with open(sa_file) as fh:
        for line in fh:
            if line.startswith("type"):
                continue
            
            
            data = line.strip().split("\t")
            
            suffix = data[7]
    
            if verbose:
                print suffix
    
            items = data[-1].split()
            elements_a = set(items)
            
            df = int(data[5])
            tf = int(data[4])
            suffix_len = int(data[6])

            counts[suffix] = {}
            
            for item in elements_a:
                
                n = items.count(item)
                counts[suffix][item] = (n, suffix_len*n, tf, df, suffix_len)
                
    with open(sa_counts_file,"w") as fw:
        for suffix in counts:
            string = "%s\t" % (suffix)
            for docid in counts[suffix]:
                string += "%s,%s\t" % (docid, str(counts[suffix][docid]))
            string += "\n"
            fw.write(string)
                
    
def compute_documnet_distances(sa_file, output_file, index_file, method="weighted_min_tf"):
    
    
    said_to_doclength = {}
    
    with open(index_file) as fh:
        for line in fh:
            data = line.strip().split("\t")
            id = int(data[0]) - 1
            l = int(data[2])
            
            said_to_doclength[id] = l
    
    
    
    
    dataset = {}
    
    n = 0
    with open(sa_file) as fh:
        for line in fh:
            n += 1
    k = 0
    with open(sa_file) as fh:
        with open(output_file, "a") as fw:
            for line in fh:
                if line.startswith("type"):
                    continue
                k += 1
                print k,"/", n
                data = line.strip().split("\t")
                items = data[-1].split()
                elements_a = set(items)
#                df = float(data[5])
                
                suffix_len = int(data[6])
                
                counts = {}
                for item in elements_a:
                    if method == "raw_min_tf":
                        counts[item] = items.count(item)
                    elif method == "weighted_min_tf":
                        counts[item] = suffix_len * items.count(item) / float(said_to_doclength[int(item)])
                    
                
                elements_b = set(items)
                
#                idf = (math.log( 1012.0/(df), 2 ))
#                
                for a in elements_a:
                    elements_b.remove(a)
                    for b in elements_b:
                        dataset.setdefault( (a,b), 0)
                        dataset[ (a,b) ] += min( counts[a], counts[b] )
                        
            for key, value in dataset.items():
                fw.write( "%s\t%s\t%s\n" % (key[0], key[1], value) )
                    
def filter_sa_dataset(sa_file, output_file, min_tf, min_df):
    dataset = []
    with open(sa_file) as fh:
        for line in fh:
            if line.startswith("type"):
                continue
            data = line.strip().split("\t")
            tf = int(data[4])
            td = int(data[5])
            if tf<min_tf and td<min_df: 
                continue
            dataset.append( (td, line) )
        dataset.sort(reverse=True)
        with open(output_file, "a") as fw:
            for item in dataset:
                fw.write(item[1])
    


def get_sa_intersection_meta(intersection_file, distance_file):
    
    intersection_file = r"M:\home\ad3002\work\comp_rat_mouse\result.sa_ma_mi"
    left_id_to_doc = r"M:\home\ad3002\work\comp_rat_mouse\doc_to_trf_file.pickle.left"
    right_id_to_doc = r"M:\home\ad3002\work\comp_rat_mouse\doc_to_trf_file.pickle.right"
    
    
    distance_file = r"M:\home\ad3002\work\comp_rat_mouse\result_trfid.sa"
    
    result = []
    
    with open(intersection_file, "r") as fh:
        for line in fh:
            if line.startswith("type"):
                continue
            data = line.strip().split("\t")
            suffix = data[0]
            left_ids = data[1].split()
            right_ids = data[2].split()
            
            result.append( (suffix, left_ids, right_ids) )
            
    # Нам нужно открыть файл с переводами id -> doc id
    # и перевести файлы
    
    left_doc_to_trf = pickle.load(open(left_id_to_doc))
    right_doc_to_trf = pickle.load(open(right_id_to_doc))
    
    with open(distance_file, "a") as fw:
        for item in result:
            
            string = "%s\t%s\t%s\n" % (item[0],  " ".join( [ str(left_doc_to_trf[int(id)]) for id in item[1] ]), " ".join([str(right_doc_to_trf[int(id)]) for id in item[2]]) )
            
            fw.write(string)
 
def get_sa_intersection(sa_a_file, sa_b_file, output_sa):
    
    suf = {}
    result = {}
        
    i = 0
    with open(sa_a_file, "r") as fh:
        for line in fh:
            i += 1
            if line.startswith("type"):
                continue
            data = line.strip().split("\t")
            suf[data[7]] = data[-1]
            print "Get suffix data %s" % i

            
    
            
    j = 0
    k = 0
    with open(sa_b_file, "r") as fh:
        for line in fh:
            j += 1
            if line.startswith("type"):
                continue
            data = line.strip().split("\t")
            if data[7] in suf:
                result[data[7]] = ["",""]
                result[data[7]][0] = suf[data[7]]
                result[data[7]][1] = data[-1]
                k += 1  
            print "Get intersection total %s got %s" % (j, k)

            
    
    
    n = 0
    with open(output_sa, "a") as fw:
        for suffix in result:
            n += 1
            print "Write result %s/%s" % (n, k)

            string = "%s\t%s\t%s\n" % (suffix, result[suffix][0], result[suffix][1])
            
            fw.write(string)
       
    
if __name__ == '__main__':
#    get_sa_intersection(r"M:\home\ad3002\work\comp_rat_mouse\out.sa.reduced_l15.left", 
#                        r"M:\home\ad3002\work\comp_rat_mouse\out.sa.reduced_l15.right", 
#                        r"M:\home\ad3002\work\comp_rat_mouse\result.sa_ma_mi")
#    get_sa_intersection_meta("","")
    pass