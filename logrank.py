#-*- coding: utf-8 -*-

####################################################################
import requests
from lifelines.statistics import logrank_test
####################################################################
## Параметры

gene_list = open("genes.txt","r").readlines()
hpv_neg_list = open("hpv_neg.txt","r").readlines()

gene_list = [gene_list[i].replace("\r\n","") for i in range(len(gene_list))]
hpv_neg_list = [string.split("\r\n")[0] for string in hpv_neg_list]

case_set_id = "hnsc_tcga_all"
genetic_profile_id = "hnsc_tcga_rna_seq_v2_mrna_median_Zscores"

down = 1
up = 3

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
def get_information_about_expression(gene_list,case_set_id,genetic_profile_id,proxy):

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
        dic[elem[0]] = {}
        for i in range(1,len(elem)):
            if elem[i] != "NaN":
                dic[elem[0]][gene_list[i-1]] = elem[i]
        if len(dic[elem[0]]) == 0:
            del dic[elem[0]]

    return dic
####################################################################
def get_information_about_mutation(gene_list,case_set_id,proxy):

    dic = {}

    gene_list.sort()
    gene_string = '+'.join(gene_list)

    link = "http://www.cbioportal.org/webservice.do?cmd=getProfileData&case_set_id="+case_set_id+"&genetic_profile_id=hnsc_tcga_mutations&gene_list="+gene_string

    data = requests.get(link,proxies=proxy).text
    data = data.split("\n")
    data = data[2:len(data)-1]

    for i in range(len(data)):
        data[i] = data[i].split("\t")[2::]

    data = transpose_list(data)


    for string in data:
        dic[string[0]] = {}
        for i in range(1,len(string)):
           dic[string[0]][gene_list[i-1]] = string[i]

    return dic

def get_information_about_amplification(gene_list,case_set_id,proxy):

    dic = {}

    gene_list.sort()
    gene_string = '+'.join(gene_list)

    link = "http://www.cbioportal.org/webservice.do?cmd=getProfileData&case_set_id="+case_set_id+"&genetic_profile_id=hnsc_tcga_gistic&gene_list="+gene_string

    data = requests.get(link,proxies=proxy).text
    data = data.split("\n")
    data = data[2:len(data)-1]

    for i in range(len(data)):
        data[i] = data[i].split("\t")[2::]

    data = transpose_list(data)


    for string in data:
        dic[string[0]] = {}
        for i in range(1,len(string)):
           dic[string[0]][gene_list[i-1]] = string[i]

    return dic

####################################################################
## Получение информации о пациентах
def get_information_about_patient(case_set_id,proxy):
    dic = {}
    data = requests.get("http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id="+case_set_id,proxies=proxy).text
    data = data.split("\n")

    for i in range(1,len(data)-1):
        data[i] = data[i].split("\t")
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

    return_list = new_lst[:]

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
def verify(lst,cutoff):
    count = 0
    for elem in lst:
        if abs(float(elem)) >= abs(cutoff):
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
expr_list = get_information_about_expression(gene_list,case_set_id,genetic_profile_id,proxy)

mut_list = get_information_about_mutation(gene_list,case_set_id,proxy)
amp_list = get_information_about_amplification(gene_list,case_set_id,proxy)



expr_list_new = {}
for case in expr_list:
    up_count = 0
    down_count = 0
    values = expr_list[case].values()
    for elem in values:
        if elem != "NaN":
            if float(elem)  <= -1.0:
                down_count += 1
            elif float(elem) >= 0:
                up_count += 1
    if up_count >= up and down_count >= down:
        expr_list_new[case] = expr_list[case]


variant_list = combine(gene_list)
expr_list = expr_list_new

write_file = open("out.txt","w")

for lst in variant_list:
    alt_time = []
    alt_status = []
    not_alt_time = []
    not_alt_status = []
    for key in patient_list:
        if key in expr_list and key in mut_list and key in amp_list and key in hpv_neg_list:
            expr_status = False
            mut_status = False
            amp_status = False

            for gene in lst:
                if abs(float(expr_list[key][gene])) >= 2:
                    expr_status = True
                    break

            for gene in lst:
                if mut_list[key][gene] != "NaN":
                    mut_status = True
                    break

            for gene in lst:
                if amp_list[key][gene] != "NaN":
                    if abs(int(amp_list[key][gene])) == 2 or int(amp_list[key][gene]) == -1:
                        amp_status = True
                        break

            if expr_status == True:
                alt_status.append(patient_list[key]["STATUS"])
                alt_time.append(patient_list[key]["MONTH"])
            elif mut_status == True:
                alt_status.append(patient_list[key]["STATUS"])
                alt_time.append(patient_list[key]["MONTH"])
            elif amp_status == True:
                alt_status.append(patient_list[key]["STATUS"])
                alt_time.append(patient_list[key]["MONTH"])
            else:
                not_alt_status.append(patient_list[key]["STATUS"])
                not_alt_time.append(patient_list[key]["MONTH"])

    write_file.write('%d\t%s\t%f\n' % (len(lst),", ".join(lst),get_logrank(alt_time, not_alt_time, alt_status, not_alt_status)))
write_file.close()







