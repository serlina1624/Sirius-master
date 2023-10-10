#!/bin/bash
# 3. Аннотация
# Вход: фаста-файл с контигами, список UCE
# Выход: outfmt - таблица со списком UCE, и контигов, содержащих их 
mkdir 3_Blastn_annotation
cd 3_Blastn_annotation
cp  ../Trinity.fasta Trinity.fasta
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
cp  ../Trinity.fasta Trinity.fasta
mkdir output
cd output
python ../../4.py
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
python ../../6.py
cd ../../

# 7. CAP3. Сборка псевдореференсев
# Вход: фаста файлы, где в каждом лежат контиги с UCE, файлы разбиты по кластерам
# Выход: фаста файлы с псевдореверенсами, на которые будет делаться выравнивание
cp -r  6_preparing_CAP3/output 7_CAP3_assembl
cp 7.py 7_CAP3_assembl/7.py
cd 7_CAP3_assembl
mkdir output
cd output
for i in ../*.fasta;
do
cap3 $i -p 75 -f 1500 -e 500 -z 1 -o 100 -s 500
done
cd ../
python 7.py
cd ../

# 8. Выравнивание в BWA
# Вход: файл с референсами, риды видов
# Выход: SAM файлы по числу видов
mkdir 8_BWA_alignment
cp  7_CAP3_assembl/output/All_refs.fasta 8_BWA_alignment/All_refs.fasta
mkdir 8_BWA_alignment/output
bwa index 8_BWA_alignment/All_refs.fasta
for i in Non_adapters*.fastq;
do
bwa mem 8_BWA_alignment/All_refs.fasta $i > 8_BWA_alignment/output/$i.sam

done

# 9. Картирование и фильтрация
# Вход: картированные Sam - файлы
# Выход: отфильтрованные sam файлы для дальнейшей обработки, #bam файлы для визуализации
cp -r 8_BWA_alignment/output/ 9_filter
cp  8_BWA_alignment/All_refs.fasta 9_filter/All_refs.fasta
cd 9_filter
for i in *.sam;
do
samtools view -F 4 $i -o cart_$i # для парсинга и чистки ридов
done
python ../9.py






























 

