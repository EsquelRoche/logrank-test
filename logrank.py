__author__ = 'bioinf'
#-*- coding: utf-8 -*-

####################################################################
import requests
import plotly.tools as tls
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
####################################################################
## Параметры

gene_list = open("genes.txt","r").readlines()
hpv_neg_list = open("hpv_neg.txt","r").readlines()

gene_list = [gene_list[i].replace("\r\n","") for i in range(len(gene_list))]
hpv_neg_list = [string.split("\r\n")[0] for string in hpv_neg_list]

case_set_id = "hnsc_tcga_all"
genetic_profile_id = "hnsc_tcga_rna_seq_v2_mrna_median_Zscores"

expr_cutoff = [0.8,1.5,0.1] # Порог экспрессии

up_cutoff = 0
down_cutoff = -1.0

proxy ={"http":"http://nakberov:hemongoo@proxy.ksu.ru:8080"}

####################################################################
## Транспонирование матрицы
def transpose_list(lst):
    x = len(lst)
    y = len(lst[0])

    return_list = []
    for i in range(y):
        return_list.append([])
        for j in range(x):
            return_list[i].append(lst[j][i])
    return return_list

####################################################################
## Получение данных об экспрессии
def get_information_about_expression(gene_list,hpv_neg_list,case_set_id,genetic_profile_id,proxy):

    dic = {}
    gene_list.sort()
    gene_string = '+'.join(gene_list)

    link = "http://www.cbioportal.org/webservice.do?cmd=getProfileData&case_set_id="+case_set_id+"&genetic_profile_id="+genetic_profile_id+"&gene_list="+gene_string
    data = requests.get(link,proxies=proxy).text

    data = data.split("\n")
    data = data[2:len(data)-1]

    for i in range(len(data)):
        data[i] = data[i].split("\t")[2::]

    data = transpose_list(data)

    for elem in data:
        if elem[0] in hpv_neg_list:
            dic[elem[0]] = {}
            for i in range(1,len(elem)):
                if elem[i] != "NaN":
                    dic[elem[0]][gene_list[i-1]] = elem[i]

    return dic
####################################################################
## Получение информации о пациентах
def get_information_about_patient(case_set_id,hpv_neg_list,proxy):
    dic = {}
    data = requests.get("http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id="+case_set_id,proxies=proxy).text

    data = data.split("\n")


    for i in range(0,len(data)):
        data[i] = data[i].split("\t")
        if data[i][0] in hpv_neg_list:

            if "" not in [data[i][47],data[i][48]] and "NA" not in [data[i][47],data[i][48]]:
                if "LIVING" == data[i][48]:
                    dic[data[i][0]] = {"MONTH":float(data[i][47]),"STATUS":0}
                elif "DECEASED" == data[i][48]:
                    dic[data[i][0]] = {"MONTH":float(data[i][47]),"STATUS":1}
    return dic
####################################################################
## Вычисдение p-value
def get_logrank(alt_time, not_alt_time, alt_status, not_alt_status):
    summary = logrank_test(alt_time, not_alt_time, alt_status, not_alt_status,alpha=0.95)
    return summary.p_value

####################################################################
## Комбинирование
def combine(lst):
    new_lst = []
    for elem in lst:
        new_lst.append([elem])

    return_list = new_lst.copy()

    for i in range(len(lst)-1):
        temp = []
        for j in range(len(new_lst)):
            for k in range(lst.index(new_lst[j][-1])+1,len(lst)):
                if lst[k] not in new_lst[j]:
                    if k == len(lst)-1:
                        return_list.append(new_lst[j]+[lst[k]])
                    else:
                        return_list.append(new_lst[j]+[lst[k]])
                        temp.append(new_lst[j]+[lst[k]])
        new_lst = temp
    return_list = [elem for elem in return_list if len(elem)>1]
    return return_list

####################################################################
def vefify(lst,cutoff,sign):
    count = 0
    if sign == "+":
        for elem in lst:
            if float(elem) >= cutoff:
                count += 1
    elif sign == "-":
        for elem in lst:
            if float(elem) <= cutoff:
                count += 1
    return count
####################################################################
def frange(start, stop, step):
    i = start
    while i < stop:
         yield i
         i += step

####################################################################
## Получение информации из БД
gene_string = "+".join([string.split("\n")[0] for string in gene_list])

patient_list = get_information_about_patient(case_set_id,proxy)
expr_list = get_information_about_expression(gene_string,hpv_neg_list,case_set_id,genetic_profile_id,proxy)
case_list = list(expr_list.keys())

variant_genes = combine(gene_list)
####################################################################
out = open("out.txt","w")


for numb in frange(expr_cutoff[0],expr_cutoff[1],expr_cutoff[2]):
    for gene_lst in variant_genes:
        for count_up in range(0,len(gene_lst)):
            for count_down in range(0,len(gene_lst)):
                alt_time = []
                alt_status = []
                not_alt_time = []
                not_alt_status = []
                i = 0
                for case in patient_list:
                    i += 1
                    temp = []
                    for gene in gene_lst:
                        temp.append(expr_list[case][gene])
                    if vefify(temp,up_cutoff,"+") >= count_up or vefify(temp,down_cutoff,"-") >= count_down:
                        #print(i,vefify(temp,up_cutoff),count_up,vefify(temp,down_cutoff),count_down,temp)
                        if vefify(temp,numb,"+") >= 1:
                            alt_time.append(patient_list[case][0])
                            alt_status.append(patient_list[case][1])
                        else:
                            not_alt_time.append(patient_list[case][0])
                            not_alt_status.append(patient_list[case][1])
                if len(alt_time) > 0 and len(not_alt_time) > 0:
                    logrank = get_logrank(alt_time, not_alt_time, alt_status, not_alt_status)
                    #print(temp,str(count_up),str(count_down),str(numb),str(logrank),str(len(alt_time)+len(not_alt_time)))
                    out.write(" ".join(gene_lst)+"\t"+str(count_up)+"\t"+str(count_down)+"\t"+str(numb)+"\t"+str(logrank)+"\t"+str(len(alt_time)+len(not_alt_time))+"\n")


out.close()
'''













