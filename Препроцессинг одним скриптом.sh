#!/bin/bash
#!/usr/bin/python
#1. Объединение видов в общий пул
# Вход: Фаста файлы видов
# Выход: Единый пул, фаста файл
cat IonXpress_*.fasta > All_spicies.fasta

#2. Работа с trinity (без доступа к кластеру работаем с готовым файлом по пяти видам)
# Вход: единый пул, фаста
# Выход: Фаста файл с контигами
mkdir 2_Trinity_de_novo
#cd 2_Trinity_de_novo
# ln -s ../All_spicies
#Trinity --seqType fa --single All_spicies.fasta
# cd ../

# 3. Аннотация
# Вход: фаста-файл с контигами, список UCE
# Выход: outfmt - таблица со списком UCE, и контигов, содержащих их 
mkdir 3_Blastn_annotation
cd 3_Blastn_annotation
cp  ../2_Trinity_de_novo/Trinity.fasta Trinity.fasta
cp ../UCE500_final.fasta UCE500_final.fasta
makeblastdb -in 'Trinity.fasta' -dbtype nucl -out blast_index
blastn -query UCE500_final.fasta -db blast_index -outfmt 6 -out UCE_in_Trinity.outfmt6
cd ../

# 4. Поиск контигов, содержащих UCE
# Вход: outfmt - таблица, фаста файл с контигами(Trinity.fasta)
# Выход: фаста файлы, по числу UCE. С контигами, содержащими эти UCE.
mkdir 4_Fasta_with_contigs_with_UCE
cd 4_Fasta_with_contigs_with_UCE
cp  ../3_Blastn_annotation/UCE_in_Trinity.outfmt6 UCE_in_Trinity.outfmt6
cp  ../2_Trinity_de_novo/Trinity.fasta Trinity.fasta
mkdir output
cd output
script=`cat <<_EOF_
refs = {}
with open('../UCE_in_Trinity.outfmt6', 'r') as file:
    for line in file:
        refs[line[line.find('\t') + 1:line.find('\t',line.find('\t')+1)]] = line[:line.find('\t')]

with open('../Trinity.fasta', 'r') as trin:
    with open('../Trinity_copy.fasta', 'w') as trin1:
        for line in trin:
            if not line.startswith('>'):
                trin1.write(line.replace('\n',''))
            else:
                trin1.write('\n' + line.replace('\n', '\t'))
                
for value, key in refs.items():
    with open('{0:s}.fasta'.format(key), 'a') as out:
        with open('../Trinity_copy.fasta', 'r') as trin:
            for i, line in enumerate(trin):
                if value in line:
                    try:
                        out.write('>' + str(i) + ' contig ' + key + '\n' + line[line.find('\t') +1:])
                    except IndexError:
                        continue
_EOF_`
echo "$script" | python
cd ../../

# 5. Кластеризация в CD-HIT-EST
# Вход: фаста файлы uce-*.fasta
# Выход: файлы двух форматов, фаста и clstr, на каждый UCE.
cp -r 4_Fasta_with_contigs_with_UCE/output 5_CD_HIT_clasterization
cd 5_CD_HIT_clasterization
for i in uce-*.fasta;
do
	cd-hit-est -i $i -c 0.8 -l 200 -r 1 -d 30 -o cluster_$i
done
cd ../

# 6. Фаста файлы с кластерами. Подготовка для cap3.
# Вход: фаста файлы с контигами, содержащими UCE, файлы clstr с разбиением на кластеры
# Выход: фаста файлы с контигами, разложенными по кластерам
cp -r  5_CD_HIT_clasterization/ 6_preparing_CAP3
cd 6_preparing_CAP3
mkdir output
cd output
script=`cat <<_EOF_
UCE_pull = [i for i in range(1, 1500)]
for uce_num in UCE_pull:

    try:
        seq_num = {}
        with open('../cluster_uce-{0:s}.fasta.clstr'.format(str(uce_num)), 'r') as clstr:
            for line in clstr:
                if line.startswith('>'):
                    cluster_num = line[1:].replace('\n','').replace(' ', '_')
                if not line.startswith('>'):
                    seq_num[(line[line.find('>') + 1: line.find('.')])] = cluster_num
    except: FileNotFoundError
    
    try: 
        with open('../uce-{0:s}.fasta'.format(str(uce_num)), 'r') as uce:
            line1 = uce.readline()
            for key, value in seq_num.items():
                with open('{0:s}'.format('uce-{0:s}_'.format(str(uce_num)) + value + '.fasta'), 'a') as out:
                    if line1.startswith('>'):
                        if line1[1:line1.find(' ')] in seq_num.keys():
                            out.write(line1)
                            line1 = uce.readline()
                            out.write(line1)
                            line1 = uce.readline()     
    except: FileNotFoundError
_EOF_`
echo "$script" | python
cd ../../

# 7. CAP3. Сборка псевдореференсев
# Вход: фаста файлы, где в каждом лежат контиги с UCE, файлы разбиты по кластерам
# Выход: фаста файлы с псевдореверенсами, на которые будет делаться выравнивание
cp -r  6_preparing_CAP3/output 7_CAP3_assembl
cd 7_CAP3_assembl
mkdir output
cd output
for i in ../*.fasta;
do
cap3 $i -p 75 -f 1500 -e 500 -z 1 -o 100 -s 500
done
script=`cat <<_EOF_
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
_EOF_`
echo "$script" | python
cd ../../

# 8. Выравнивание в BWA
# Вход: файл с референсами, риды видов
# Выход: SAM файлы по числу видов
mkdir 8_BWA_alignment
cp  7_CAP3_assembl/output/All_refs.fasta 8_BWA_alignment/All_refs.fasta
mkdir 8_BWA_alignment/output
bwa index 8_BWA_alignment/All_refs.fasta
for i in Ion*.fasta;
do
bwa mem 8_BWA_alignment/All_refs.fasta $i > 8_BWA_alignment/output/$i.sam
done

# 9. Индексция и фильтрация
# Вход: Sam - файлы
# Выход: отфильтрованные sam файлы для дальнейшей обработки, bam файлы для визуализации
cp -r 8_BWA_alignment/output/ 9_SAMtools_filter
cd 9_SAMtools_filter
for i in *.sam;
do
samtools view -F 4 $i -o cart_$i # для парсинга и чистки ридов
samtools sort $i -o $i.bam
samtools index $i.bam
done





























 

