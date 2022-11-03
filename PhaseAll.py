#!/usr/bin/env python
# coding: utf-8

import pysam
import collections
from operator import itemgetter
import re
import os
import sys

default_sam = sys.argv[1]

# сделать фасту из файла из сиквенсов разной длины (после find_subseq)
def make_fasta(haplo_in):
    haplo_out = []
    seq = ''.join(haplo_in)
    seq_length = len(seq)
    counter = (seq_length // 61) + 1
    start = 0
    end = 61
    i = 0
    test = []
    while i < counter:
        subline = seq[start:end] + '\n'
        start += 61
        end += 61
        i += 1
        haplo_out.extend(subline)
    if i == counter:
        end = seq_length - start + 1
        start += 61
        subline = seq[start:end] + '\n'
        haplo_out.extend(subline)
    return ''.join(haplo_out)

def letter_occurence(my_list):
    string = ''.join(my_list)
    total = len(string)
    my_list_of_nucl = [['A', string.count('A') / total], ['G', string.count('G') / total],                        ['C', string.count('C') / total], ['T', string.count('T') / total]]
    result = sorted(my_list_of_nucl, key=itemgetter(1))
    return result
    
my_list = [['A', 'A', 'A', 'A', 'G', 'G'], ['A', 'A', 'A', 'G', 'G', 'G']]
letter_occurence(my_list[1])

def find_allele(query_name1_stack, query_name2_stack, query_name1, query_name2, haplo1, haplo2, top1_nucl, top2_nucl):
    # Первое вхождение полиморфной позиции
    if len(query_name1_stack) == 0 and len(query_name2_stack) == 0:
        # в лист записываются листы с названиями ридов        
        query_name1_stack.append(query_name1)
        query_name2_stack.append(query_name2)
        haplo1.append(top1_nucl[0])
        haplo2.append(top2_nucl[0])
        return haplo1, haplo2, query_name1_stack, query_name2_stack
    # Второе+ вхождение полиморфной позиции
    else:
        # если убрать этот фильтр, можно получить запись всех ридов, относящихся к
        # аллелям, и собрать их потом в ассемблере
        query_name1_stack = [query_name1_stack[-1]]
        query_name2_stack = [query_name2_stack[-1]]
        # проверка на пробел между аллелями
        past = set(query_name1_stack[-1] + query_name2_stack[-1])
        present = set(query_name1 + query_name2)
        if len(past.intersection(present)) == 0:
            query_name1_stack.append(query_name1)
            query_name2_stack.append(query_name2)
            haplo1.append('X')
            haplo2.append('X')
            haplo1.append(top1_nucl[0])
            haplo2.append(top2_nucl[0])
            return haplo1, haplo2, query_name1_stack, query_name2_stack
        else:
            flag = 0
            for name in query_name1:
                if name in query_name1_stack[-1]:
                    query_name1_stack.append(query_name1)
                    query_name2_stack.append(query_name2)
                    haplo1.append(top1_nucl[0])
                    haplo2.append(top2_nucl[0])
                    flag = 1
                    return haplo1, haplo2, query_name1_stack, query_name2_stack

                elif name in query_name2_stack[-1]:
                    query_name2_stack.append(query_name1)
                    query_name1_stack.append(query_name2)
                    haplo2.append(top1_nucl[0])
                    haplo1.append(top2_nucl[0])
                    flag = 1
                    return haplo1, haplo2, query_name1_stack, query_name2_stack
            if flag == 0:
                for name in query_name2:
                    if name in query_name2_stack[-1]:   
                        query_name2_stack.append(query_name2)
                        query_name1_stack.append(query_name1)
                        haplo2.append(top2_nucl[0])
                        haplo1.append(top1_nucl[0])
                        return haplo1, haplo2, query_name1_stack, query_name2_stack
                    elif name in query_name1_stack[-1]:
                        query_name1_stack.append(query_name2)
                        query_name2_stack.append(query_name1)
                        haplo1.append(top2_nucl[0])
                        haplo2.append(top1_nucl[0])
                        return haplo1, haplo2, query_name1_stack, query_name2_stack

# ищем вставки
def find_insertions(column_feature, pileupcolumn, column_names):
    column_names1 = []
    column_names2 = []
    insertion1_count = 0
    insertion1_1_count = 0
    insertion2_count = 0
    insertion1 = 0 # вставка 
    # если вставки на всех ридах, и они оданаковые, или разные = записываем в разные insertion1/2
    if ''.join(column_feature).count('+') > int(pileupcolumn.nsegments) * 0.9:
        for i in range(0, len(column_feature)):
            if i == 0:
                try:
                    insertion1 = re.search('[0-9](.*?)$', column_feature[i]).group(1)
                    insertion1_count += 1
                    column_names1.append(column_names[i])
                except AttributeError:
                    continue
            else:
                try:
                    insertion2 = re.search('[0-9](.*?)$', column_feature[i]).group(1)
                    insertion2_count += 1
                except AttributeError:
                    continue
                if insertion2 != insertion1:
                    column_names2.append(column_names[i])
                else:
                    column_names1.append(column_names[i])
        
    # если есть вставка, но не у всех ридов:
    else:
        insertion2 = 'P' # костыль
        # нужно выявить вставку. Если присутствует ошибочная минорная вставка с мутацией, её
        # необходимо выявить и удалить
        features = pileupcolumn.get_query_sequences(add_indels=True)
        for i in range(0, len(features)):
            try:
                current_insertion = re.search('[0-9](.*?)$', features[i]).group(1).upper()
            except AttributeError:
                continue
            if insertion1 == 0:
                insertion1 = current_insertion
            else:
                if insertion1 == current_insertion:
                    insertion1_count += 1
                    continue
                else:
                    insertion1_1 = current_insertion
                    insertion1_1_count += 1
        for i in range(0, len(column_feature)):
            # следующая строка выполняет ту же функцию что try except, не дает возникнуть
            # ошибке при чтении 'T'
            if re.search('\+', column_feature[i]):
                if insertion1 == re.search('[0-9](.*?)$', column_feature[i]).group(1):
                    insertion1_count += 1
                    column_names1.append(column_names[i])
            else:
                insertion2 = '-' * len(insertion1)
                insertion2_count += 1
                column_names2.append(column_names[i])
    return insertion1, column_names1, insertion2, column_names2, insertion1_count, insertion2_count


# ищем делеции        
def find_deletions(column_feature, column_names):
    column_names1 = []
    column_names2 = []
    deletion_count = 0
    non_deletion_count = 0
    for i in range(0, len(column_feature)):
        if re.search('\*', column_feature[i]):
            deletion1 = '-'
            deletion_count += 1
            column_names1.append(column_names[i])
        # есть баг, если это последний нуклеотид в риде, т.е. рид 'прилегает' к инделю, то он 
        # тоже учтется в подсчете
        else: 
            column_names2.append(column_names[i])
            non_deletion_count += 1
    return deletion1, column_names1, column_names2, deletion_count, non_deletion_count

samfile = pysam.AlignmentFile(default_sam, "r" )

query_name1 = []
query_name2 = []
query_name1_stack = []
query_name2_stack = []

name_nucl = []
name_nucl_pileupcolumn = []
name_nucl_stack = []

haplo1 = []
haplo2 = []

column_feature = []
column_names = []

deletion_flag = 0
insertion_flag = 0

for pileupcolumn in samfile.pileup(flag_require = False):
    column_feature = []
    column_names = []
    column_feature.extend(pileupcolumn.get_query_sequences(add_indels=True))
    column_names.extend(pileupcolumn.get_query_names())
    indel_count = ''.join(column_feature).count('+')
    if indel_count / len(column_feature) >= 0.2: # or indel_count >= 3: 
        insertion_flag = 1
        # это костыль, проверка повторяется дальше, пока лень переделывать. + надо
        # переделать весь блок индлей, чтобы просто считывалось с column_feature *(=индель)
        insertion1, column_names1, insertion2, column_names2, insertion1_count,         insertion2_count = find_insertions(column_feature, pileupcolumn, column_names)
        
    # ищем делеции
    indel_count = ''.join(column_feature).count('*')
    if indel_count / len(column_feature) >= 0.2: # or indel_count >= 3:
        deletion_flag = 1
        deletion1, column_names1, column_names2, deletion_count, non_deletion_count =         find_deletions(column_feature, column_names)
    nucleotide = []
    name_nucl_pileupcolumn = []
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            # добавляем для каждой позиции [имя рида], [нуклеотид]
            name_nucl.append(pileupread.alignment.query_name)
            name_nucl.append(pileupread.alignment.query_sequence[pileupread.query_position])
            # создаем лист из нуклеотидов в каждой позиции
            nucleotide.append(pileupread.alignment.query_sequence[pileupread.query_position])
            # создаем лист листов со всеми позициями
            name_nucl_pileupcolumn.append(name_nucl)
            name_nucl = []
    # проверка на глубину прочтения > 0, а также общие для двух аллелей делеции
    if len(nucleotide) == 0:
        deletion_flag = 0
        haplo1.append('-')
        haplo2.append('-')
        continue
    # проверка на делеции. Определяем будут ли следующие позиции инделями, и с 
    # следующей итерации начинаем их записывать столько раз, сколько длина делеции. Записываем
    # в один аллель делецию, в другую нуклеотид. 
    if deletion_flag == 1:
        deletion_flag = 0
        top1_nucl = [deletion1, deletion_count]
        top2_nucl = letter_occurence(nucleotide)[-1]
        top2_nucl[1] = non_deletion_count
        # проверка на ошибочные нуклеотиды в позиции где есть делеции
        if top2_nucl[1] / pileupcolumn.nsegments < 0.2 and top2_nucl[1] < 3:
            # если нуклеотидов меньше 20% или 3 штук, то в аллели мы записываем делецию
            haplo1.append(top1_nucl[0])
            haplo2.append(top1_nucl[0])
            continue
        else:
            query_name1 = column_names1
            query_name2 = column_names2
            # intersection(query_name1_stack, query_name2_stack, query_name1, query_name2)
            haplo1, haplo2, query_name1_stack, query_name2_stack = find_allele(query_name1_stack,             query_name2_stack, query_name1, query_name2, haplo1, haplo2, top1_nucl, top2_nucl)
            continue

    else:
        # два наиболее часто встречающихся нуклеотида
        top1_nucl = letter_occurence(nucleotide)[-1]
        top2_nucl = letter_occurence(nucleotide)[-2]
        # если в позиции только один вариант
        if top1_nucl[1] == 1:
            haplo1.append(top1_nucl[0])
            haplo2.append(top1_nucl[0])
        elif top2_nucl[1] < 0.2: # and nucleotide.count(top2_nucl[0]) < 3:
            haplo1.append(top1_nucl[0])
            haplo2.append(top1_nucl[0])
        else:
            # создаем лист в котором хранятся листы в виде 'колонкок' из [имя рида], [нуклеотид] 
            # для всех полиморфных позиций
            name_nucl_stack.append(name_nucl_pileupcolumn)
            # записываем названия ридов, содержащие разные нуклеотиды
            query_name1 = []
            query_name2 = []
            for i in name_nucl_stack[-1]:
                # i = [имя рида], [нуклеотид]
                if i[1] == top1_nucl[0]:
                    query_name1.append(i[0])
                elif i[1] == top2_nucl[0]:
                    query_name2.append(i[0])
                else:
                    continue
            # intersection(query_name1_stack, query_name2_stack, query_name1, query_name2)
            haplo1, haplo2, query_name1_stack, query_name2_stack = find_allele(query_name1_stack,             query_name2_stack, query_name1, query_name2, haplo1, haplo2, top1_nucl, top2_nucl)

        # проверка, являются ли следующие позиции инделями
        # индели показываются в позиции после настоящей, поэтому мы обрабатываем их 
        # после того как закончили обрабатывать нуклеотиды. Вставки обрабатываем также 
        # как и snp, записывая в переменные top1(2)_nucl [вставку], [частоту], а в переменные
        # query_name1(2) названия ридов, ассоциированных с вставками
        if insertion_flag == 1:
            insertion_flag = 0
            if len(column_names1) == pileupcolumn.nsegments:
                haplo1.append(insertion1)
                haplo2.append(insertion1)
            else:
                top1_nucl = [insertion1, insertion1_count]
                top2_nucl = [insertion2, insertion2_count]
                if min(top1_nucl[1], top2_nucl[1]) / pileupcolumn.nsegments < 0.2: 
                    if top1_nucl[1] > top2_nucl[1]:
                        haplo1.append(top1_nucl[0])
                        haplo2.append(top1_nucl[0])
                    else:
                        haplo1.append(top2_nucl[0])
                        haplo2.append(top2_nucl[0])
                else:
                    query_name1 = column_names1
                    query_name2 = column_names2
                    haplo1, haplo2, query_name1_stack, query_name2_stack = find_allele(query_name1_stack,                     query_name2_stack, query_name1, query_name2, haplo1, haplo2, top1_nucl, top2_nucl)

output_file = default_sam + '_alleles.fas'

with open(output_file, 'w') as output:
    haplo_out1 = make_fasta(haplo1)
    haplo_out2 = make_fasta(haplo2)
    output.write('>allele_1\n' + haplo_out1 + '\n' + '>allele_2\n' + haplo_out2)

samfile.close()
