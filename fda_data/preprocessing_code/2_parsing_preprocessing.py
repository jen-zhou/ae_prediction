# PREPROCESSING
import pickle
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import glob
from tqdm import tqdm
from datetime import date, datetime

# maps noisy term i.e. 11-beta-hydroxylase deficiency to a list of MedDRA codes and PTs for all levels of hierarchical classification
se_dic = pickle.load(open('/PHShome/jz1082/AE/fda_study/fda_data/preprocessing_code/AE_dic.pk', 'rb'))
# maps drug names to list [DrugBank ID, some internal identifier?]
drug_dic = pickle.load(open('/PHShome/jz1082/AE/fda_study/fda_data/preprocessing_code/drug_mapping.pk', 'rb'))
# dataframe with cols PT, PT_name, HLT, HLT_name, HLGFT, HLGFT_name, SOC, SOC_name, SOC_abbr
meddra_pd_all =pickle.load(open('/PHShome/jz1082/AE/fda_study/fda_data/preprocessing_code/AE_mapping.pk', 'rb'))

def days_to_date(days):
    stand_date = date(2000, 1, 1)

    if int(days)<0:
        days = 1
    
    dt = datetime.fromordinal(int(days))
    return dt.strftime('%Y-%m-%d')

# initial setup
def date_normalize(formate, dat): 
    stand_date = date(2000, 1, 1)
    if formate=='102':  # the date is formed as yyyymmdd
        current_date = date(int(dat[:4]), int(dat[4:6]), int(dat[6:8])) 
    elif formate=='610':  # formed as yyyymm  
        current_date = date(int(dat[:4]), int(dat[4:6]), 1)
    elif formate=='602':  #formed as yyyy      
        current_date = date(int(dat[:4]), 1, 1)
    delta = current_date - stand_date
    return delta.days

n_reports = []
miss_count = {}

for yr in range(2015, 2025):
    qtr_list=[1,2,3,4]
    for qtr in qtr_list:
        qtr_name = str(yr)+'q'+ str(qtr)
        print('I am parsing:',qtr_name)
        
        lab_storage = '/PHShome/jz1082/AE/fda_study/fda_data/raw_data/'

        files = lab_storage + qtr_name + '/**/**'
        xml_files = glob.glob(files +"/*.xml", recursive=True)
        unique_files = list(set(xml_files))  # only keep the unique values, remove duplicated files.
        xml_files = unique_files
        xml_files.sort()
        print('find {} files'.format(len(xml_files)))
        print(xml_files)

        root = None
        for xml_file in xml_files:
            print(xml_file)
            data = ET.parse(xml_file).getroot()
            if root is None:
                root = data
            else:
                root.extend(data)
                print('finished merge',xml_file)
        nmb_reports = len(root)
        print(nmb_reports)

        count = 0
        patient_ID = 0
        dic = {}
        
        miss_admin = miss_patient = miss_reaction = miss_drug =0
        for report in tqdm(root.findall('safetyreport')):
            """Administrative Information"""
#             report.find('').text
            try:  # Mandatory Information: report_id
                try:
                    version = report.find('safetyreportversion').text
                except:
                    version = '1'
                    
                report_id = report.find('safetyreportid').text
                
                try:
                    case_id = report.find('companynumb').text
                except:
                    case_id = '0'  # unknown case id
                    
                try:
                    country = report.find('primarysource')[0].text
                except:
                    country = 'unknown'          

                    
                if country =='COUNTRY NOT SPECIFIED':
                    country = 'unknown'
                    
                    
                try:
                    qualify = report.find('primarysource')[1].text
                except:
                    qualify = '6'  # the qualify is unknown
                    
#                 qualify = report.find('primarysource')[1].text
                    
                if qualify not in {'1', '2', '3', '4', '5', '6','7'}:
                    qualify = '0'
                                      
                    
                try:
                    serious = report.find('serious').text
                except:
                    serious = '-1'
                
                try:
                    s_1 = report.find('seriousnessdeath').text
                except:
                    s_1 = '0'
                try:
                    s_2 = report.find('seriousnesslifethreatening').text
                except:
                    s_2 = '0'
                try:
                    s_3 = report.find('seriousnesshospitalization').text
                except:
                    s_3 = '0'
                try:
                    s_4 = report.find('seriousnessdisabling').text
                except:
                    s_4 = '0'
                try:
                    s_5 = report.find('seriousnesscongenitalanomali').text
                except:
                    s_5 = '0'
                try:
                    s_6 = report.find('seriousnessother').text
                except:
                    s_6 = '0'
                serious_subtype = [s_1, s_2, s_3, s_4, s_5, s_6]
            except:
                miss_admin +=1
                continue

            try:  # Optional information
                # receivedate: Date when the report was the FIRST received
                receivedateformat, receivedate = report.find('receivedateformat').text, report.find('receivedate').text
                receivedate = date_normalize(receivedateformat, receivedate)
            except:
                receivedate = '0'
            
            try:
                # receiptdate: Date of most RECENT report received
                receiptdateformat, receiptdate = report.find('receiptdateformat').text, report.find('receiptdate').text
                receiptdate = date_normalize(receiptdateformat, receiptdate)
            except:
                 receiptdate =  '0'

            for patient in report.findall('patient'):
                """Demographic Information"""                
                try:
                    age = patient.find('patientonsetage').text
                except:
                    age = -1 # unknown age
                try:
                    ageunit = patient.find('patientonsetageunit').text
                except:
                    ageunit = '801' 
                # normalize age
                try:
                    age = int(age)  
                    if age!= -1:
                        if ageunit == '800':  # Decade 
                            age = '-1'
                        elif ageunit == '801':  # Year
                            age = age
                        elif ageunit == '802':  # Month
                            age = int(age/12)
                        elif ageunit == '803':  # Week
                            age = int(age/52)
                        elif ageunit == '804':  # Day
                            age = int(age/365)
                        elif ageunit == '805':  # Hour
                            age = int(age/(24*365))
                except:
                    age = -1
                    
                      
                try:
                    gender = patient.find('patientsex').text
                except:
                    gender = '0'
                try:
                    weight = patient.find('patientweight').text
                except:
                    weight = '0'

                reaction_list = []
                for side_ in patient.findall('reaction'):
                    try:  # outcome: 1-6, 6 levels in total
                        try: 
                            PT_code = side_[0].text
                        except:
                            PT_code = '0'
                        try:
                            outcome = side_[2].text
                        except:
                            outcome = '6'
                        try:
                            PT = side_[1].text
                        except:
                            PT = 'none'
                        reaction = [PT_code, PT, outcome]
                    except:
                        continue
                    reaction_list.append(reaction) 
                if reaction_list.__len__() == 0:  # Mandatory condition: at least has one reaction
                    miss_reaction += 1
                    continue

                drug_list = []
                for drug_ in patient.findall('drug'):
                    try:
                        try:
                            char =  drug_.find('drugcharacterization').text  # drugcharacterization: 1(suspect)/2(concomitant)/3(interacting)
                        except:
                            char = '0'
                        try:
                            product =  drug_.find('medicinalproduct').text  # drug brand
                        except:
                            product = 'none'
                        """Dosage are generally fixed according to the indication"""
                        try: 
                            dorse, unit=  drug_.find('drugstructuredosagenumb').text, drug_.find('drugstructuredosageunit').text
                            drugseparatedosagenumb, drugintervaldosageunitnumb, drugintervaldosagedefinition = \
                                drug_.find('drugseparatedosagenumb').text, drug_.find('drugintervaldosageunitnumb').text, \
                                drug_.find('drugintervaldosagedefinition').text
                            form = drug_.find('drugdosageform').text  # tablet or capsule or sth 
                        except:
                            dorse, unit, drugseparatedosagenumb,drugintervaldosageunitnumb, drugintervaldosagedefinition, form =\
                            '0', '0', '0','0','0', '0'
                        try:
                            route = drug_.find('drugadministrationroute').text
                            if route == '048':
                                route = '1'  # oral 
                            elif route == '061':
                                route = '2'  # Topical
                        except:
                            route = '0'  # no information of route
                        
                        try:
                            indication = drug_.find('drugindication').text  # indication (disease): super important
                        except:
                            indication = 'none'

                        try:
                            start_format, start_date = drug_.find('drugstartdateformat').text, drug_.find('drugstartdate').text
                            start_date = date_normalize(start_format, start_date)
                        except:
                            start_date = '0'
                        try:
                            end_format, end_date = drug_.find('drugenddateformat').text, drug_.find('drugenddate').text
                            end_date = date_normalize(end_format, end_date)
                        except:
                            try:
                                end_date = receiptdate
                            except:
                                end_date = '0'
                            
                        try:
                            action = drug_.find('actiondrug').text
                        except:
                            action = '5'
                        try:
                            additional = drug_.find('drugadditional').text
                        except:
                            additional = '3'
                        try:
                            readm = drug_.find('drugrecurreadministration').text
                        except:
                            readm = '3'
                        try:
                            substance = drug_.find('activesubstance')[0].text
                        except:
                            substance = 'none'
                    except:  # Mandatory condition: if none of the above information is provided, ignore this report
                        continue
                    drug = [char, product, dorse, unit, drugseparatedosagenumb, drugintervaldosageunitnumb,
                            drugintervaldosagedefinition, form, route, indication, start_date, end_date, action,
                            readm, additional, substance]
                    drug_list.append(drug)
                if drug_list.__len__() ==0:
                    miss_drug += 1
                    continue

                """for patient_ID"""
                dic[count] = [version, report_id, case_id, country, qualify, serious, 
                              s_1, s_2, s_3, s_4, s_5, s_6, 
                              receivedate, receiptdate,  
                              age, gender, weight, reaction_list, drug_list]
                count += 1
                break

        pickle.dump(dic, open('../parsed_data/'+ qtr_name+'.pk', 'wb'))
        
        n_reports.append(len(dic))
        print(qtr_name+' file saved. with', len(dic), 'reports')
        miss_count[qtr_name] = [nmb_reports, miss_admin, miss_patient, miss_reaction, miss_drug]

print ('All data saved')

ho_combos = {}

nmb = 0 
n_reports = []
for yr in range(2015, 2025):
    qtr_list = [1, 2, 3, 4]
    for qtr in qtr_list:
        qtr_name = str(yr) + 'q' + str(qtr)
        file_path = './parsed_data/' + qtr_name + '.pk'
        print('I am loading {} from {}'.format(qtr_name, file_path))
        dic = pickle.load(open(file_path, 'rb'))
        n_reports.append(len(dic))
        print('loaded', len(dic))
        values = np.array(list(dic.values()))
        print(values.shape)   

        admin_demo = values[:, :17]
        se = values[:, 17] 
        drugs = values[:, 18]              

        for i in tqdm(range(values.shape[0]), desc=f"Processing {qtr_name} reports"): 
            new_se = []
            new_drug = [] 
            new_indication = [] 

            se_report = se[i] 
            drugs_report = drugs[i]

            for j in range(len(se_report)):
                se_reaction = se_report[j][1].lower() 
                se_key = None
                if '\\' in se_reaction:
                    se_reaction = se_reaction.split('\\')[0]

                if se_reaction in se_dic:
                    se_key = se_reaction
                else:
                    tokens = se_reaction.split(' ')
                    for token in tokens:
                        if token in se_dic:
                            se_key = token
                            break

                try: 
                    new_se.append(se_dic[se_key][0]) 
                except:  
                    continue

            for k in range(len(drugs_report)):
                drug = drugs_report[k][-1].lower()   
                indication = drugs_report[k][9].lower()  

                if '\\' in drug:
                    drug = drug.split('\\')[0] 
                if drug in drug_dic:
                    drug_key = drug
                elif ' ' in drug:
                    key_sets = drug.split(' ')
                    if key_sets[0] in drug_dic:
                        drug_key = key_sets[0]
                    elif len(key_sets) > 1 and key_sets[1] in drug_dic:
                        drug_key = key_sets[1]
                    else:
                        continue
                else:
                    continue

                try:
                    new_drug.append(drug_dic[drug_key][0])
                except:
                    continue

                indication_key = None
                if '\\' in indication:
                    indication = indication.split('\\')[0]

                match = meddra_pd_all[meddra_pd_all['PT_name'].str.lower() == indication.lower()]
                if not match.empty:
                    indication_key = match.iloc[0]['PT']
                else:
                    for token in indication.split(' '):
                        match = meddra_pd_all[meddra_pd_all['PT_name'].str.lower() == token.lower()]
                        if not match.empty:
                            indication_key = match.iloc[0]['PT']
                            break

                if indication_key:
                    new_indication.append(indication_key)

            new_se.sort()
            new_drug.sort()
            new_indication.sort()

            xx = list(admin_demo[i])
            xx.append(new_se)
            xx.append(new_drug)
            xx.append(new_indication)
            ho_combos[nmb] = xx
            nmb += 1
            
    print('#- all reports', sum(n_reports))
    print('The No. of high order combos:', len(ho_combos))
    pickle.dump(ho_combos, open(f'../parsed_data/reports_v4_{yr}.pk', 'wb'))

print('ho_combos saved', ho_combos.get(0))

# # Years you want to combine
# years = list(range(2015, 2025))  # Change this range to include all the years you have

# combined_combos = {}
# nmb = 0

# for yr in years:
#     try:
#         file_path = f'/PHShome/jz1082/ae_prediction/fda_study/fda_data/parsed_data/reports_v4_{yr}.pk'
#         print(f"Loading {file_path}")
#         yearly_data = pickle.load(open(file_path, 'rb'))

#         for _, value in yearly_data.items():
#             combined_combos[nmb] = value
#             nmb += 1

#     except FileNotFoundError:
#         print(f"File for year {yr} not found, skipping...")

# print(f"Total entries combined: {len(combined_combos)}")
# output_path = '/PHShome/jz1082/ae_prediction/fda_study/fda_data/parsed_data/reports_v4.pk'
# pickle.dump(combined_combos, open(output_path, 'wb'))
# print(f"Saved combined data to {output_path}")

# # import time
# print("Starting to load reports_v4.pk...")
# # start = time.time()

for year in range(2015,2025):

    with open(f'/PHShome/jz1082/ae_prediction/fda_study/fda_data/parsed_data/reports_v4_{year}.pk', 'rb') as f:
        reports_v4 = pickle.load(f)

    print(f"reports_v4 has {len(reports_v4)} entries")

    reports_pd = pd.DataFrame(reports_v4.values(), 
                            columns=['version','report_id','case_id','country','qualify','serious',
                                    's1','s2','s3','s4','s5','s6','receivedate','receiptdate',
                                    'age','gender','weight','SE','drugs','indications'])

    tqdm.pandas()
    reports_pd['date'] = reports_pd.progress_apply(
        lambda row: str(pd.Period('2000-01-01') + int(row['receivedate'])), axis=1
    )

    print('#-of all reports',reports_pd.shape)
    reports_pd = reports_pd.drop_duplicates(['report_id', 'case_id','receivedate'], keep="last")
    print('After remove duplicate reports',reports_pd.shape)

    # Convert with error handling
    reports_pd['receiptdate'] = pd.to_datetime(reports_pd['receiptdate'], errors='coerce')
    reports_pd['receivedate'] = pd.to_datetime(reports_pd['receivedate'], errors='coerce')

    # Now subtract safely (NaT will propagate)
    reports_pd['lastingdays'] = (reports_pd['receiptdate'] - reports_pd['receivedate']).dt.days

    pickle.dump(reports_pd, open(f'/PHShome/jz1082/ae_prediction/fda_study/fda_data/parsed_data/patient_safety_{year}.pk', 'wb'))


