#!/bin/python
import os 
curpath = os.path.abspath(os.path.curdir) # зафиксируем папку 
files_fasta_list = [x for x in os.listdir(curpath) if x.startswith("mafft_")] # список видов
for file in files_fasta_list:
    with open(file, 'r') as UCE:
        with open("New_" +file[6:], 'a') as out:
            for line in UCE:
                if line.startswith(">"):
                    out.write(">" + line[-7:])
                    line = UCE.readline()
                out.write(line)
