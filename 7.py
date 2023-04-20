#!/usr/bin/python
# 7. CAP3. Сборка псевдореференсев
# Вход: фаста файлы, где в каждом лежат контиги с UCE, файлы разбиты по кластерам
# Выход: фаста файлы с псевдореверенсами, на которые будет делаться выравнивание

UCE_pull = [i for i in range(1, 1500)]
for uce_num in UCE_pull:
    try:
        with open('../uce-{0:s}_Cluster_0.fasta.cap.contigs'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
            with open('../uce-{0:s}_Cluster_0_copy.fasta.cap.contigs'.format(str(uce_num)), 'w') as contig1:
                for line in contig:
                    if not line.startswith('>'):
                        contig1.write(line.replace('\n',''))
                    else:
                        contig1.write('\n' + line)
                contig1.write('\n')
    except: FileNotFoundError
     
    try:
        with open('../uce-{0:s}_Cluster_1.fasta.cap.contigs'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
            with open('../uce-{0:s}_Cluster_1_copy.fasta.cap.contigs'.format(str(uce_num)), 'w') as contig1:
                for line in contig:
                    if not line.startswith('>'):
                        contig1.write(line.replace('\n',''))
                    else:
                        contig1.write('\n' + line)
                contig1.write('\n')
    except: FileNotFoundError
    
    try:
        with open('../uce-{0:s}_Cluster_2.fasta.cap.contigs'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
            with open('../uce-{0:s}_Cluster_2_copy.fasta.cap.contigs'.format(str(uce_num)), 'w') as contig1:
                for line in contig:
                    if not line.startswith('>'):
                        contig1.write(line.replace('\n',''))
                    else:
                        contig1.write('\n' + line)
                contig1.write('\n')
    except: FileNotFoundError
    
    
    try:
        with open('../uce-{0:s}_Cluster_0.fasta.cap.singlets'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
            with open('../uce-{0:s}_Cluster_0_copy.fasta.cap.singlets'.format(str(uce_num)), 'w') as contig1:
                for line in contig:
                    if not line.startswith('>'):
                        contig1.write(line.replace('\n',''))
                    else:
                        contig1.write('\n' + line)
                contig1.write('\n')
    except: FileNotFoundError
     
    try:
        with open('../uce-{0:s}_Cluster_1.fasta.cap.singlets'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
            with open('../uce-{0:s}_Cluster_1_copy.fasta.cap.singlets'.format(str(uce_num)), 'w') as contig1:
                for line in contig:
                    if not line.startswith('>'):
                        contig1.write(line.replace('\n',''))
                    else:
                        contig1.write('\n' + line)
                contig1.write('\n')
    except: FileNotFoundError
    
    try:
        with open('../uce-{0:s}_Cluster_2.fasta.cap.singlets'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
            with open('../uce-{0:s}_Cluster_2_copy.fasta.cap.singlets'.format(str(uce_num)), 'w') as contig1:
                for line in contig:
                    if not line.startswith('>'):
                        contig1.write(line.replace('\n',''))
                    else:
                        contig1.write('\n' + line)
                contig1.write('\n')
    except: FileNotFoundError
with open('All_refs.fasta', 'a') as ref:
    for uce_num in UCE_pull:
        try:
            with open('../uce-{0:s}_Cluster_0_copy.fasta.cap.contigs'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
                line = contig.readline()
                line = contig.readline()
                if line.startswith('>'): # проверка на пустой файл
                    while line: 
                        if line.startswith('>'): # проверка на названия контига
                            ref.writelines('> UCE{0:s}_Claster_0_'.format(str(uce_num)) + line[1:]) # записываем название контига и его uce
                            line = contig.readline()
                            ref.writelines(line)
                            line = contig.readline()
                else:
                    with open('../uce-{0:s}_Cluster_0_copy.fasta.cap.singlets'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
                        line = contig.readline()
                        if line.startswith('>'): # проверка на пустой файл
                            ref.writelines('> UCE{0:s}_Claster_0 \n'.format(str(uce_num))) # записываем название контига и его uce
                            while line:
                                line = contig.readline()
                                ref.writelines(line)
                      
        except: FileNotFoundError
        
        try:
            with open('../uce-{0:s}_Cluster_1_copy.fasta.cap.contigs'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
                line = contig.readline()
                line = contig.readline()
                if line.startswith('>'): # проверка на пустой файл
                    while line: 
                        if line.startswith('>'): # проверка на названия контига
                            ref.writelines('> UCE{0:s}_Claster_1_'.format(str(uce_num)) + line[1:]) # записываем название контига и его uce
                            line = contig.readline()
                            ref.writelines(line)
                            line = contig.readline()
                else:
                    with open('../uce-{0:s}_Cluster_1_copy.fasta.cap.singlets'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
                        line = contig.readline()
                        if line.startswith('>'): # проверка на пустой файл
                            ref.writelines('> UCE{0:s}_Claster_1 \n'.format(str(uce_num))) # записываем название контига и его uce
                            while line:
                                line = contig.readline()
                                ref.writelines(line)
                      
        except: FileNotFoundError
        
        try:
            with open('../uce-{0:s}_Cluster_2_copy.fasta.cap.contigs'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
                line = contig.readline()
                line = contig.readline()
                if line.startswith('>'): # проверка на пустой файл
                    while line: 
                        if line.startswith('>'): # проверка на названия контига
                            ref.writelines('> UCE{0:s}_Claster_2_'.format(str(uce_num)) + line[1:]) # записываем название контига и его uce
                            line = contig.readline()
                            ref.writelines(line)
                            line = contig.readline()
                else:
                    with open('../uce-{0:s}_Cluster_2.fasta_copy.cap.singlets'.format(str(uce_num)), 'r') as contig: # переписать для всех файлов
                        line = contig.readline()
                        if line.startswith('>'): # проверка на пустой файл
                            ref.writelines('> UCE{0:s}_Claster_2 \n'.format(str(uce_num))) # записываем название контига и его uce
                            while line:
                                line = contig.readline()
                                ref.writelines(line)
                      
        except: FileNotFoundError
