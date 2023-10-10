#!/usr/bin/python
# 9. Индексция и фильтрация
# Вход: Sam - файлы
# Выход: отфильтрованные sam файлы для дальнейшей обработки, bam файлы для визуализации

import sys
import os
import pandas as pd
from Bio import SeqIO
class CIGAR(object):
    '''Класс для парсинга SAM файлов'''
    def __init__(self, name, ref_name, start, CIGAR, read): # задаем объект, где как атрибуты будут имя рида, флаг и тд
        self.name = name
        self.ref_name = ref_name
        self.start = start
        self.CIGAR = CIGAR
        self.read = read
        
    def Ref(self): # Подгрузим соответствующий референс в виде списка букв
         with open('All_refs.fasta', 'r') as refs:
            line = refs.readline()
            ref_all = []
            while line:
                if self.ref_name in line: # Если название референса из атрибутов совпадает с названием строки в фаста файле с референсами, следующую строку записываем как референс
                    line = refs.readline()
                    while not line.startswith(">"):
                        ref = list(line)[:-1]
                        ref_all = ref_all + ref
                        line = refs.readline()
                line = refs.readline()  
         return ref_all
    
    def ReadySeq(self):
        '''For input the self object, for output the seq is ready for comparing with reference'''
        cigar = list(self.CIGAR.strip()) #сделали список
        # для начала нужно сигар разбить на блоки и переделать цифры в инт. По идее олжно получиться int(61) str(S) int(37) str(M)...
        # Каждый символ пытаюсь перевести в инт, обработка исключений на буквы
        for i in range(len(cigar)): #пробуем каждый символ списка сделать интом
            try:
                cigar[i] = int(cigar[i])
            except ValueError:
                cigar[i] = cigar[i]

        # Потом объединяю цифры в одно число.
        cigar_new = []
        for i in range(len(cigar)): # если символ - строка - записываем его отдельно
            if type(cigar[i]) == str:
                cigar_new.append(cigar[i])                                    
            try:
                if type(cigar[i]) == type(cigar[i+1]) == type(cigar[i+2]) : # если три символа подряд это инт, объединяем их
                    cigar_new.append(int(str(cigar[i])+str(cigar[i+1])+str(cigar[i+2])))
                if type(cigar[i]) == type(cigar[i+1]) and type(cigar[i]) != type(cigar[i+2]) and type(cigar[i-1]) ==str:  # если два символа подряд одного типа (такое может быть только с инт), объединяем их  
                    cigar_new.append(int(str(cigar[i])+str(cigar[i+1])))
                if type(cigar[i+1]) == str and type(cigar[i-1]) == str: # если один инт - тоже записываем его
                    cigar_new.append(cigar[i])
            except IndexError: # это нужно на случай окончания строки, list index out of range
                continue
                
        # удобнее сделать словарь
        numbers = cigar_new[::2] # список с числом нуклеотидов
        clipping_type = [] # список с типом совпадений
        
        for i in range(1, len(cigar_new), 2): 
            clipping_type.append(cigar_new[i]);

        cigar_dict = [] # делаем словарь
        for i in range(len(numbers)):
            cigar_dict.append((numbers[i],clipping_type[i]))
            
        # Теперь прогоняем сам seq по строке cigar, для каждой буквы свое действие
        nucl_number = 0 # сразу столбец из data
        seq = self.read
        for e in cigar_dict:
            if e[1] == 'S':
                seq = seq[:nucl_number] + 'N' * e[0] + seq[nucl_number + e[0]:]
            if e[1] == 'D':
                seq = seq[:nucl_number] + ' ' * e[0] + seq[nucl_number:]
            if e[1] == 'H':
                continue # это пока, может быть hard clipping не вырезается
            if e[1] == 'I':
                seq = seq[:nucl_number] + 'N' * e[0] + seq[nucl_number + e[0]:]
            nucl_number += e[0]
        seq = seq.replace('N', '') 
        new_seq = (self.start-1) * ' ' + seq #Тут может быть ошибка, мб не надо -1

        return new_seq
        
def RefCover(refs, reads):
    '''Делает словарь покрытий по конкретному референсу и сборке.
    На вход: Референс по буквам и Список ридов, 
    на выходе список словарей с покрытием для каждого нуклеотида в референсе'''
    ref_cover_list = [] # большой словарь для всех нуклеотидов, по всей длине референса
    for i in range(len(refs) + 1):
        cover_dict = {} # для каждого нуклеотида
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

    # посчитаем покрытие для каждого нуклеотида
    cover_list = [] # большой словарь для всех нуклеотидов, по всей длине референса
    for i in range(len(refs) + 1):
        cover_dict = {} # для каждого нуклеотида
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
    '''Смотрим на все существующие замены в этой сборке для этой хромосомы.
        На вход берем референс и риды сборки. На выходе получаем список списков, для каждого рида свой список.
        Внутри два элемента, номер рида и номер нуклеотида с типом замены'''
    changing = []
    for nucl_number, nucl in enumerate(refs): # для всех нуклеотидов референса
        for seq_number, seq in enumerate(reads): # для каждого рида сборки
            #print(seq[nucl_number])
            try:
                if seq[nucl_number] != ' ': # если нуклеотид рида существует в этой позиции
                    if nucl != seq[nucl_number]: # если нуклеотид референса не равен нуклеотиду рида
                        changing.append([Names[seq_number], 
                                         '{0:s}'.format(seq[nucl_number], nucl_number +1), len(seq.replace(' ', ''))])                 
            except IndexError:
                continue
    return changing
    
def ReadWithChanges(AllChanges):
    '''Считаем количество ридов с заменами и считаем число замен в риде.
        На вход подаем список списков из функции AllChanges.
        На выход список кортежей с номером рида и числом замен в этом риде'''
    read_with_change = [] # Риды с заменами
    for e in AllChanges: # Для каждого элемента из общего файла существующих замен:
        read_with_change.append(e[0]) # Добавляем сюда номер рида

    change = [] 
    for read in read_with_change:
        change.append(('{0:s}'.format(read), read_with_change.count(read))) # Считаем сколько раз этот рид встретился в наборе всех замен
    read_with_change = list(set(change)) # Берем сет от этого набора:
    return read_with_change
    
def Cleaning(ReadWithChanges, AllChanges, RefCover):
    candidate_for_remove = [] # список кандидатов на удаление
    for read_tuple in ReadWithChanges: # для всех ридов с заменами:
        for i in range(len(AllChanges)): # пройдемся по списку замен
            if read_tuple[0] == AllChanges[i][0]: # если рид с заменой в списке всех изменений, берем его длину
                if read_tuple[1]*100/AllChanges[i][-1] >= 10: # Более скольки процентов ридов
                    candidate_for_remove.append(read_tuple)

        # те риды, где есть инсерция и делеция, тоже попадают в список
       # if 'I' in data.loc[data['Название рида'] == read_tuple[1]]['CIGAR'].iloc[0] or 'D' in data.loc[data['Название рида'] == read_tuple[0]]['CIGAR'].iloc[0]:
        #    candidate_for_remove.append(read_tuple)
    #candidate_for_remove = list(set(candidate_for_remove)) # сет, чтобы один рид не попал дважды
    #print(candidate_for_remove)

    count =[]
    for read_number, read in enumerate(candidate_for_remove): # в каждом риде с большим числом замен
        k = 0 # Счетчик редких замен
        for i, e in enumerate(AllChanges): # иду по списку всех возможных замен
            if read[0] == e[0]: # если совпадает название рида с ошибкой и название кандидата на удаление
                for cover_i, cover in enumerate(RefCover[0]): # смотрим словарь покрытия
                    for key, value in cover.items(): # key - тип нуклеотида с заменой и его номер, value - процент покрытия для этого типа замены
                        if key == e[1]: # когда нуклеотид совпадает с нуклеотидом из списка замен
                            if value * RefCover[1][cover_i] / max(RefCover[1]) < 40: # для скольки процентов ридов характнрна такая замена
                                #print('Процент замен в риде, относительно его длины:', read[1],'\n', 'Информация о риде и замене:',e, '\n','Процент такой замены в целом:',value)
                                k+=1
        if k !=0:
            count.append((read,k)) # сколько в риде всего замен, и сколько из них "Редкие"
    reads_for_remove = []
    for i in count:
        if i[1] / i[0][1] > 0.4: #если относительно всех замен рида более 40% из них редкие,
            reads_for_remove.append(i[0][0])
        
    return reads_for_remove
    



curpath = os.path.abspath(os.path.curdir) # зафиксируем папку    
files_sam_list = [x for x in os.listdir(curpath) if x.startswith("cart_")] # endswith - кончается на (можно так сделать для начала)

for sam_file in files_sam_list:
    #with open("Clear_{0:s}".format(sam_file), 'a') as output:
     #   with open("HEADER_{0:s}".format(sam_file), "r") as file:
      #      for line in file:
       #         if line.startswith("@"):
        #            output.write(line)
    print("Идет обработка файла: ", sam_file)
    data = pd.read_csv(sam_file, sep = '\t',header=None, usecols=[0,3,9,5,1,2], index_col = False)
    data.columns = ['Название рида','Flag','Ref','Начало рида', 'CIGAR','Рид']
    data = data.loc[data['Flag'] != 2048] # чистка по флагам
    data = data.loc[data['Flag'] != 2064]
    data = data.reset_index(drop=True)
    refs_in_sam = set(data["Ref"]) # в sam файле смотрим сколько всего референсов.
    for ref in refs_in_sam: # открываем один из референсов: (весь алгоритм был расчитан на такую сборку)
        reads = [] # в списке будут лежать риды для этого конкретного референса
        Cigar = []
        Names = []
        for i in range(len(data['Ref'])):
            if data['Ref'][i].encode() == ref.encode():
                reads.append(CIGAR(data.iloc[i][0], data.iloc[i][2], data.iloc[i][3], data.iloc[i][4], data.iloc[i][5]).ReadySeq())
                refs = CIGAR(data.iloc[i][0], data.iloc[i][2], data.iloc[i][3], data.iloc[i][4], data.iloc[i][5]).Ref() # сам референс по буквам, перенос строки не берем
                Cigar.append(data.iloc[i][4])
                Names.append(data.iloc[i][0])
        print('Очистка сборки на референс: ', ref)
        bad_reads = Cleaning(ReadWithChanges(AllChanges(refs, reads, Names)), AllChanges(refs, reads, Names), RefCover(refs, reads))
        good_reads = []
        for read in Names:
            if read not in bad_reads:
                good_reads.append(read + '\t')
        print("Число плохих ридов",len(bad_reads), ", Число хороших ридов", len(good_reads), " Для референса ", ref, "В виде: ", sam_file)
        with open("Clear_{0:s}".format(sam_file), 'a') as output:
            with open(sam_file, "r") as file:
                for lines in file:
                    for read in good_reads:
                        if read in lines:
                            if ref in lines:
                                #print(read, lines)
                                output.write(lines)
                                
curpath = os.path.abspath(os.path.curdir) # зафиксируем папку    
files_sam_list = [x for x in os.listdir(curpath) if x.startswith("Clear_")] # endswith - кончается на (можно так сделать для начала)

for sam_file in files_sam_list:
    data = pd.read_csv(sam_file, sep = '\t',header=None, usecols=[0,3,9,5,1,2], index_col = False)
    data.columns = ['Название рида','Flag','Ref','Начало рида', 'CIGAR','Рид']
    data = data.loc[data['Flag'] != 2048] # чистка по флагам
    data = data.loc[data['Flag'] != 2064]
    data = data.reset_index(drop=True)
    refs_in_sam = set(data["Ref"]) # в sam файле смотрим сколько всего референсов.
    for ref in refs_in_sam: # открываем один из референсов: (весь алгоритм был расчитан на такую сборку)
        reads = [] # в списке будут лежать риды для этого конкретного референса
        refs = []
        for i in range(len(data['Ref'])):
            if data['Ref'][i].encode() == ref.encode():
                reads.append(CIGAR(data.iloc[i][0], data.iloc[i][2], data.iloc[i][3], data.iloc[i][4], data.iloc[i][5]).ReadySeq())
                refs = CIGAR(data.iloc[i][0], data.iloc[i][2], data.iloc[i][3], data.iloc[i][4], data.iloc[i][5]).Ref()[:-1] # сам референс по буквам, перенос строки не берем
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
            output.write(">" + ref + "\n")
            output.write(fragment + '\n')
            
            
import os 
curpath = os.path.abspath(os.path.curdir) # зафиксируем папку 
files_fasta_list = [x for x in os.listdir(curpath) if x.startswith("Seq_")] # список видов
with open("All_refs.fasta", 'r') as refs:
    for ref_line in refs:
        if ref_line.startswith(">"):
            with open('{0:s}.fasta'.format(ref_line[1:-1]), 'a') as UCE:
                for spicies in files_fasta_list:
                    with open(spicies, 'r') as sp:
                        for line_sp in sp:
                            if line_sp.startswith(">"):
                                if line_sp == ref_line:
                                    UCE.write(ref_line[:-1] + "_Sp_" + spicies[41:44] + "\n")
                                    UCE.write(sp.readline())


