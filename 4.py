#!/usr/bin/python
# 4. Поиск контигов, содержащих UCE
# Вход: outfmt - таблица, фаста файл с контигами(Trinity.fasta)
# Выход: фаста файлы, по числу UCE. С контигами, содержащими эти UCE.

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

