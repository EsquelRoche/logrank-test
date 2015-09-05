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
for i in range(len(gene_list)):
    gene_list[i] = gene_list[i].replace("\n","")
    
hpv_neg_list = [string.split("\n")[0] for string in hpv_neg_list]

case_set_id = "hnsc_tcga_all"
genetic_profile_id = "hnsc_tcga_rna_seq_v2_mrna_median_Zscores"
expr_cutoff = [0.8,1.5,0.1]

up_cutoff = 0
down_cutoff = -1.0

proxy ={"http":"http://nakberov:hemongoo@proxy.ksu.ru:8080",
        "https":"https://nakberov:hemongoo@proxy.ksu.ru:8080"}

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
def get_information_about_expression(gene_string,hpv_neg_list,case_set_id,genetic_profile_id,proxy):
    link = "http://www.cbioportal.org/webservice.do?cmd=getProfileData&case_set_id="+case_set_id+"&genetic_profile_id="+genetic_profile_id+"&gene_list="+gene_string
    get_data = requests.get(link,proxies=proxy)
    data = get_data.text  
    data = data.split("\n")
    data = data[2:len(data)-1]

    for i in range(len(data)):
        data[i] = data[i].split("\t")

    dic = {}
    for i in range(2,len(data[0])):
        if data[0][i] in hpv_neg_list:
            dic[data[0][i]] = {}
            for j in range(1,len(data)):
                
                if data[j][i] != "NaN":
                    dic[data[0][i]][data[j][1]] = float(data[j][i])
    return dic

####################################################################
## Получение информации о пациентах    
def get_information_about_patient(case_set_id,proxy):
    link = "http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id="+case_set_id
    get_data = requests.get(link,proxies=proxy)
    data = get_data.text

    data = data.split("\n")

    dic = {}
    
    for i in range(len(data)-1):
        data[i] = data[i].split("\t")
        if data[i][0] in hpv_neg_list:
            if "" not in [data[i][30],data[i][31]] and "NA" not in [data[i][30],data[i][31]]:
                if "LIVING" == data[i][31]:
                    dic[data[i][0]] = [float(data[i][30]),0]
                elif "DECEASED" == data[i][31]:
                    dic[data[i][0]] = [float(data[i][30]),1]
    return dic

####################################################################
## Вычисдение p-value
def get_logrank(alt_time, not_alt_time, alt_status, not_alt_status):
    summary = logrank_test(alt_time, not_alt_time, alt_status, not_alt_status,alpha=0.95)
    return summary.p_value

####################################################################
## Построение графика
def get_plot(alt_time, not_alt_time, alt_status, not_alt_status,name):
    kmf = KaplanMeierFitter()

    ax = plt.subplot(111)

    kmf.fit(alt_time, event_observed=alt_status,label=['Cases with Alteration(s)'])
    kmf.plot(ax=ax,flat=True,ci_show=False)
    kmf.fit(not_alt_time, event_observed=not_alt_status,label=['Cases without Alteration(s)'])
    kmf.plot(ax=ax,flat=True,ci_show=False)
    
    ax.set_xlim(0,250)
    plt.savefig(name+'.png')

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

   












