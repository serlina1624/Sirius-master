#!/usr/bin/python
# 7. CAP3. Сборка псевдореференсев
# Вход: фаста файлы, где в каждом лежат контиги с UCE, файлы разбиты по кластерам
# Выход: фаста файлы с псевдореверенсами, на которые будет делаться выравнивание
import os
UCE_pull = [i for i in range(1, 1500)]
clst_pull = [0, 1, 2, 3, 4, 5, 6]
with open('output/All_refs.fasta', 'a') as ref:
    for uce_num in UCE_pull:
        for clst_num in clst_pull:
            try:
                isempty = os.stat('uce-{0:s}_Cluster_{1:s}.fasta.cap.contigs'.format(str(uce_num), str(clst_num))).st_size == 0
                if isempty == False:
                    with open('uce-{0:s}_Cluster_{1:s}.fasta.cap.contigs'.format(str(uce_num), str(clst_num)), 'r') as contig: 
                        for line in contig:
                            if line.startswith('>'): # проверка на названия контига
                                ref.writelines('>UCE{0:s}_Claster_{1:s}_'.format(str(uce_num), str(clst_num)) + line[1:]) # записываем название контига и его uce
                            else:
                                ref.writelines(line)


                else:
                    with open('uce-{0:s}_Cluster_{1:s}.fasta.cap.singlets'.format(str(uce_num), str(clst_num)), 'r') as contig:
                        contig_number = 0
                        for line in contig:
                            if line.startswith('>'): # проверка на названия контига
                                contig_number += 1
                                ref.writelines('>UCE{0:s}_Claster_{1:s}_Contig{2:s}'.format(str(uce_num), str(clst_num), str(contig_number)) + "\n") # записываем название контига и его uce
                            else:
                                ref.writelines(line)

            except: FileNotFoundError
