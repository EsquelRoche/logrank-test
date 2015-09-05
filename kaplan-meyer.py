__author__ = 'bioinf'


from lifelines.statistics import logrank_test

results = logrank_test(T[dem], T[~dem], C[dem], C[~dem], alpha=.99 )

results.print_summary()


T = data["duration"]
C = data["observed"]


dic = get_information_about_patient(case_set_id,hpv_neg_list,proxy)

alt_time = []
not_alt_time = []
alt_status = []
not_alt_status = []

for key in dic:
    if dic[key]["STATUS"] == 0:
        alt_time.append(dic[key]["MONTH"])
        alt_status.append(dic[key]["STATUS"])
    else:
        not_alt_time.append(dic[key]["MONTH"])
        not_alt_status.append(dic[key]["STATUS"])

print(get_logrank(alt_time, not_alt_time, alt_status, not_alt_status))