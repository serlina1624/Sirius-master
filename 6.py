#!/usr/bin/python
# 6. Фаста файлы с кластерами. Подготовка для cap3.
# Вход: фаста файлы с контигами, содержащими UCE, файлы clstr с разбиением на кластеры
# Выход: фаста файлы с контигами, разложенными по кластерам

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
