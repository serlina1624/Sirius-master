#!/usr/bin/python
# 1. Обрезка адаптеров
# Вход: файлы по видам, fastq
# Выход: fastq - файлы по видам, с обрезанными адаптерами
import os
curpath = os.path.abspath(os.path.curdir) # зафиксируем папку    
files_fastq_list = [x for x in os.listdir(curpath) if x.startswith("IonXpress_")] # endswith - кончается на 
with open ("Barcodes.fasta", 'r') as file:
    for line in file: 
        if line.startswith(">"):
            line = line[1:-1]
            for fastq in files_fastq_list:
                if line in fastq:
                    adapter = file.readline()[:-1]
                    os.system('cutadapt -a {0:s} -o Non_adapters_{1:s} {1:s}'.format(adapter, fastq))

