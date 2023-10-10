#!/usr/bin/python

import sys
import os
import pandas as pd
from Bio import SeqIO
class CIGAR(object):
    def __init__(self, name, ref_name, start, CIGAR, read): 
        self.name = name
        self.ref_name = ref_name
        self.start = start
        self.CIGAR = CIGAR
        self.read = read
        
    def Ref(self): 
         with open('All_refs.fasta', 'r') as refs:
            line = refs.readline()
            ref_all = []
            while line:
                if self.ref_name in line: 
                    line = refs.readline()
                    while not line.startswith(">"):
                        ref = list(line)[:-1]
                        ref_all = ref_all + ref
                        line = refs.readline()
                line = refs.readline()  
         return ref_all
    
    def ReadySeq(self):
        '''For input the self object, for output the seq is ready for comparing with reference'''
        cigar = list(self.CIGAR.strip()) 

        for i in range(len(cigar)):
            try:
                cigar[i] = int(cigar[i])
            except ValueError:
                cigar[i] = cigar[i]

        cigar_new = []
        for i in range(len(cigar)): 
            if type(cigar[i]) == str:
                cigar_new.append(cigar[i])                                    
            try:
                if type(cigar[i]) == type(cigar[i+1]) == type(cigar[i+2]) :
                    cigar_new.append(int(str(cigar[i])+str(cigar[i+1])+str(cigar[i+2])))
                if type(cigar[i]) == type(cigar[i+1]) and type(cigar[i]) != type(cigar[i+2]) and type(cigar[i-1]) ==str:    
                    cigar_new.append(int(str(cigar[i])+str(cigar[i+1])))
                if type(cigar[i+1]) == str and type(cigar[i-1]) == str: 
                    cigar_new.append(cigar[i])
            except IndexError: 
                continue
                
        
        numbers = cigar_new[::2] 
        clipping_type = [] 
        
        for i in range(1, len(cigar_new), 2): 
            clipping_type.append(cigar_new[i]);

        cigar_dict = [] 
        for i in range(len(numbers)):
            cigar_dict.append((numbers[i],clipping_type[i]))
            
        
        nucl_number = 0 
        seq = self.read
        for e in cigar_dict:
            if e[1] == 'S':
                seq = seq[:nucl_number] + 'N' * e[0] + seq[nucl_number + e[0]:]
            if e[1] == 'D':
                seq = seq[:nucl_number] + ' ' * e[0] + seq[nucl_number:]
            if e[1] == 'H':
                continue 
            if e[1] == 'I':
                seq = seq[:nucl_number] + 'N' * e[0] + seq[nucl_number + e[0]:]
            nucl_number += e[0]
        seq = seq.replace('N', '') 
        new_seq = (self.start-1) * ' ' + seq 

        return new_seq
        
def RefCover(refs, reads):

    ref_cover_list = []
    for i in range(len(refs) + 1):
        cover_dict = {} 
        cover_dict['A'] = 0
        cover_dict['C'] = 0
        cover_dict['G'] = 0
        cover_dict['T']= 0
        ref_cover_list.append(cover_dict)

    for read in reads:
        for nucl_number, nucl in enumerate(read):
                #print(nucl, nucl_number)
                if nucl == 'A':
                    ref_cover_list[nucl_number]['A']  += 1
                if nucl == 'T':
                    ref_cover_list[nucl_number]['T']  += 1
                if nucl == 'G':
                    ref_cover_list[nucl_number]['G']  += 1
                if nucl == 'C':
                    ref_cover_list[nucl_number]['C']  += 1

    
    cover_list = [] 
    for i in range(len(refs) + 1):
        cover_dict = {} 
        cover_dict['A'] = 0
        cover_dict['C'] = 0
        cover_dict['G'] = 0
        cover_dict['T'] = 0
        cover_list.append(cover_dict)


    cover_abs = []
    for nucl_number, nucl in enumerate(ref_cover_list):
        cover = 0
        for key, value in nucl.items():
            cover += value
        for key, value in nucl.items():    
            if key == 'A':
                if cover == 0:
                    cover_list[nucl_number][key] =  0
                else:
                    cover_list[nucl_number][key] =  value * 100 / cover
            if key == 'T':
                if cover == 0:
                    cover_list[nucl_number][key] =  0
                else:
                    cover_list[nucl_number][key] =  value * 100 / cover
            if key == 'G':
                if cover == 0:
                    cover_list[nucl_number][key] =  0
                else:
                    cover_list[nucl_number][key] =  value * 100 / cover
            if key == 'C':
                if cover == 0:
                    cover_list[nucl_number][key] =  0
                else:
                    cover_list[nucl_number][key] =  value * 100 / cover
        cover_abs.append(cover)     
    return cover_list, cover_abs
    
def AllChanges(refs, reads, Names):

    changing = []
    for nucl_number, nucl in enumerate(refs): 
        for seq_number, seq in enumerate(reads): 
            #print(seq[nucl_number])
            try:
                if seq[nucl_number] != ' ': 
                    if nucl != seq[nucl_number]: 
                        changing.append([Names[seq_number], 
                                         '{0:s}'.format(seq[nucl_number], nucl_number +1), len(seq.replace(' ', ''))])                 
            except IndexError:
                continue
    return changing
    
def ReadWithChanges(AllChanges):

    read_with_change = []
    for e in AllChanges: 
        read_with_change.append(e[0]) 

    change = [] 
    for read in read_with_change:
        change.append(('{0:s}'.format(read), read_with_change.count(read))) 
    read_with_change = list(set(change)) 
    return read_with_change
    
def Cleaning(ReadWithChanges, AllChanges, RefCover):
    candidate_for_remove = [] 
    for read_tuple in ReadWithChanges: 
        for i in range(len(AllChanges)): 
            if read_tuple[0] == AllChanges[i][0]: 
                if read_tuple[1]*100/AllChanges[i][-1] >= 10: 
                    candidate_for_remove.append(read_tuple)


    count =[]
    for read_number, read in enumerate(candidate_for_remove):
        k = 0
        for i, e in enumerate(AllChanges):
            if read[0] == e[0]: 
                for cover_i, cover in enumerate(RefCover[0]): 
                    for key, value in cover.items(): 
                        if key == e[1]: 
                            if value * RefCover[1][cover_i] / max(RefCover[1]) < 40: 
                              
                                k+=1
        if k !=0:
            count.append((read,k)) 
    reads_for_remove = []
    for i in count:
        if i[1] / i[0][1] > 0.4:
            reads_for_remove.append(i[0][0])
        
    return reads_for_remove
    



curpath = os.path.abspath(os.path.curdir)    
files_sam_list = [x for x in os.listdir(curpath) if x.startswith("cart_")] 

for sam_file in files_sam_list:

    data = pd.read_csv(sam_file, sep = '\t',header=None, usecols=[0,3,9,5,1,2], index_col = False)
    data.columns = ['name','Flag','Ref','start', 'CIGAR','read']
    data = data.loc[data['Flag'] != 2048] 
    data = data.loc[data['Flag'] != 2064]
    data = data.reset_index(drop=True)
    refs_in_sam = set(data["Ref"]) 
    
    for i, ref in enumerate(refs_in_sam): 
        reads = [] 
        Cigar = []
        Names = []
        for i in range(len(data['Ref'])):
            if data['Ref'][i].encode() == ref.encode():
                reads.append(CIGAR(data.iloc[i][0], data.iloc[i][2], data.iloc[i][3], data.iloc[i][4], data.iloc[i][5]).ReadySeq())
                refs = CIGAR(data.iloc[i][0], data.iloc[i][2], data.iloc[i][3], data.iloc[i][4], data.iloc[i][5]).Ref() 
                Cigar.append(data.iloc[i][4])
                Names.append(data.iloc[i][0])
        
        bad_reads = Cleaning(ReadWithChanges(AllChanges(refs, reads, Names)), AllChanges(refs, reads, Names), RefCover(refs, reads))
        good_reads = []
        for read in Names:
            if read not in bad_reads:
                good_reads.append(read + '\t')
        
        with open("Clear_{0:s}".format(sam_file), 'a') as output:
            with open(sam_file, "r") as file:
                for lines in file:
                    for read in good_reads:
                        if read in lines:
                            if ref in lines:
                                #print(read, lines)
                                output.write(lines)
                                
for sam_file in files_sam_list:
    data = pd.read_csv(sam_file, sep = '\t',header=None, usecols=[0,3,9,5,1,2], index_col = False)
    data.columns = ['name','Flag','Ref','start', 'CIGAR','read']
    data = data.loc[data['Flag'] != 2048] 
    data = data.loc[data['Flag'] != 2064]
    data = data.reset_index(drop=True)
    refs_in_sam = set(data["Ref"]) 
    for ref in refs_in_sam:
        reads = [] 
        refs = []
        for i in range(len(data['Ref'])):
            if data['Ref'][i].encode() == ref.encode():
                reads.append(CIGAR(data.iloc[i][0], data.iloc[i][2], data.iloc[i][3], data.iloc[i][4], data.iloc[i][5]).ReadySeq())
                refs = CIGAR(data.iloc[i][0], data.iloc[i][2], data.iloc[i][3], data.iloc[i][4], data.iloc[i][5]).Ref()
        with open("Seq_by_{0:s}.fasta".format(sam_file), 'a') as output:
            fragment = ""
            for nucleotide in RefCover(refs, reads)[0]:
                if nucleotide['A'] >= 61:
                    fragment += "A"
                if nucleotide['C'] >= 61:
                    fragment += "C"
                if nucleotide['G'] >= 61:
                    fragment += "G"
                if nucleotide['T'] >= 61:
                    fragment += "T"    
                if (nucleotide['A'] >= 40) & (nucleotide['A'] <= 60) & (nucleotide['C'] >= 40) & (nucleotide['C'] <= 60): # A C
                    fragment += "M"
                if (nucleotide['A'] >= 40) & (nucleotide['A'] <= 60) & (nucleotide['G'] >= 40) & (nucleotide['G'] <= 60): # A G
                    fragment += "R"
                if (nucleotide['A'] >= 40) & (nucleotide['A'] <= 60) & (nucleotide['T'] >= 40) & (nucleotide['T'] <= 60): # A T
                    fragment += "W"
                if (nucleotide['G'] >= 40) & (nucleotide['G'] <= 60) & (nucleotide['C'] >= 40) & (nucleotide['C'] <= 60): # c G
                    fragment += "S"
                if (nucleotide['T'] >= 40) & (nucleotide['T'] <= 60) & (nucleotide['C'] >= 40) & (nucleotide['C'] <= 60): # c t
                    fragment += "Y"
                if (nucleotide['T'] >= 40) & (nucleotide['T'] <= 60) & (nucleotide['G'] >= 40) & (nucleotide['G'] <= 60): # T G
                    fragment += "K"
                #if (nucleotide['C'] == 0) & (nucleotide['A'] == 0) & (nucleotide['T'] == 0) & (nucleotide['G'] == 0):
                #    fragment += "N"
            #print(ref, fragment)
            output.write("> " + ref + "\n")
            output.write(fragment + '\n')
            
files_fasta_list = [x for x in os.listdir(curpath) if x.startswith("Seq_")]
for spicies in files_fasta_list: 
    with open(spicies, 'r') as sp:
        for line in sp: 
            with open("All_refs.fasta", 'r') as ref: 
                for i in ref: 
                    if i.startswith('>'): 
                        with open('{0:s}.fasta'.format(i[2:]), 'a') as UCE:
                            #print(i[2:])
                            if i[2:] in line: 
                                UCE.write('> ' + i[2:-1] + "_sp_"  + spicies[53:56] + " |"  +i[2:-1] + "\n") 
                                UCE.write(sp.readline()) 

