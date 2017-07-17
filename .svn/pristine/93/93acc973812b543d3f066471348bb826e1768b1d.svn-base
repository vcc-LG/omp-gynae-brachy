"""
Code to query OP database, retrieve hex XML, process + convert to BrachyPlan object
"""
import pyodbc
from lxml import etree
import pyTG43.pyTG43 as pyTG43

def connect_to_database():
    cnxn = pyodbc.connect("Driver={SQL Server};Server=10.69.174.20,53092;Database=Nucletron;Trusted_connection=yes")
    cur = cnxn.cursor()
    return cur


def get_plan_from_database(cur, patient_ID, study_ID, plan_number):
    cur.execute("SELECT channels FROM [Nucletron].[dbo].[BrachyTreatmentPlan.ApplicationSetups] applsetuptable INNER JOIN [Nucletron].[dbo].[PatientAdmin.SeriesItems] seriesitemstable ON applsetuptable.parentOID = seriesitemstable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Series] seriestable ON seriesitemstable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.BaseSeries] basetable ON basetable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Studies] studytable ON seriestable.parentOID = studytable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Patients] patienttable ON studytable.parentOID = patienttable.OID WHERE patienttable.patientID = ? AND studytable.studyID = ? AND basetable.modality = 'RTPLAN' AND basetable.number = ?",(patient_ID,study_ID,plan_number))
    cursor_content = [row for row in cur]
    plan_content =  cursor_content[0][0].decode()
    return plan_content

def process_plan_xml(plan_content):
    tree = etree.fromstring(plan_content.encode())
    plan_data_dict = {}
    plan_data_dict['sources'] = {}
    plan_data_dict['sources']['coordinates'] = []
    plan_data_dict['sources']['dwell_times'] = []
    plan_data_dict['sources']['cumulative_dwell_times'] = []
    for channel_number in range(len(tree.xpath('.//Channel'))):
        channel_total_time = float(tree.xpath('.//TotalTime')[channel_number].text)
        for control_point_number in range(len(tree.xpath('.//Channel')[channel_number].xpath('.//BrachyControlPoint'))):
            dwell_relative_time = float(tree.xpath('.//Channel')[channel_number].xpath('.//BrachyControlPoint')[control_point_number].getchildren()[5].text)
            plan_data_dict['sources']['cumulative_dwell_times'].append((dwell_relative_time/100) * channel_total_time)
            plan_data_dict['sources']['coordinates'].append([0.1*float(i) for i in tree.xpath('.//Channel')[channel_number].xpath('.//BrachyControlPoint')[control_point_number].getchildren()[4:5][0].text.split(':')])
    plan_data_dict['sources']['dwell_times'] = [x - y for x,y in zip(plan_data_dict['sources']['cumulative_dwell_times'][1::2],plan_data_dict['sources']['cumulative_dwell_times'][0::2])]
    plan_data_dict['sources']['coordinates'] = [i for i in plan_data_dict['sources']['coordinates'][0::2]]
    return plan_data_dict

def calculate_doses_from_plan(plan_data_dictionary, points_of_interest):
    my_source_train = []
    for i in range(len(plan_data_dictionary['sources']['coordinates'])):
        my_source_train.append(pyTG43.SourcePosition(x=plan_data_dictionary['sources']['coordinates'][i][0],
                                              y=plan_data_dictionary['sources']['coordinates'][i][1],
                                              z=plan_data_dictionary['sources']['coordinates'][i][2],
                                              dwell_time=plan_data_dictionary['sources']['dwell_times'][i],
                                              apparent_activity=10,
                                              Sk=40820,
                                              dose_rate_constant=1.109,
                                              L=0.36,
                                              t_half=73.83))
    for point in points_of_interest:
        point = [0.1* i for i in point]
        # import ipdb; ipdb.set_trace()
        my_dose = pyTG43.calculate_dose(my_source_train, point)
        print(my_dose)



if __name__ == "__main__":
    cur = connect_to_database()
    patient_ID = 'zzHDRcommissioning'
    study_ID = 'HDRCommissioning'
    plan_number = '18'
    plan_content = get_plan_from_database(cur, patient_ID, study_ID, plan_number)
    plan_data_dictionary = process_plan_xml(plan_content)

    points_of_interest = [[-9.8,-18,-20.71],[-9.8,-23,-20.71],[-9.8,-28,-30.71],[-9.8,-28,-10.71]]
    calculate_doses_from_plan(plan_data_dictionary, points_of_interest)
