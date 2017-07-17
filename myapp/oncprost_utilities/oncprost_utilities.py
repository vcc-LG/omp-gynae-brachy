import pyodbc
from lxml import etree

def connect_to_database():
    cnxn = pyodbc.connect("Driver={SQL Server};Server=10.69.173.148,53092;Database=Nucletron;Trusted_connection=yes")
    cur = cnxn.cursor()
    return cur

def get_patient_IDs_and_names(patient_name):
    """Return patient names and IDs from search query"""
    cur = connect_to_database()
    cur.execute("""SELECT pt.patientID, pt.name FROM [Nucletron].[dbo].[PatientAdmin.Patients] pt WHERE pt.name LIKE '%s' OR pt.name LIKE '%s'""" % ('%' + patient_name + '%','%' + patient_name + '%'))
    cursor_content = [row for row in cur]
    patient_IDs = [x[0] for x in cursor_content]
    patient_names = [x[1] for x in cursor_content]
    return patient_IDs, patient_names

def get_patient_studies(patient_ID):
    cur = connect_to_database()
    cur.execute("""SELECT st.studyID FROM [Nucletron].[dbo].[PatientAdmin.Studies] st INNER JOIN [Nucletron].[dbo].[PatientAdmin.Patients] pt ON st.parentOID = pt.OID WHERE pt.patientID = '%s'""" % patient_ID)
    cursor_content = [row for row in cur]
    studies = [i[0] for i in cursor_content]
    return studies


def get_patient_name(patient_ID):
    """Return patient name from ID"""
    cur = connect_to_database()
    cur.execute("""SELECT pt.name FROM [Nucletron].[dbo].[PatientAdmin.Patients] pt WHERE pt.patientID = '%s'""" % patient_ID)
    cursor_content = [row for row in cur]
    patient_name = cursor_content[0][0]
    return patient_name

def get_plans_from_study(patient_ID, study_ID):
    """Returns list of plans for patient specified ID, case"""
    cur = connect_to_database()
    cur.execute("""SELECT basetable.modality, basetable.number, basetable.description FROM [Nucletron].[dbo].[PatientAdmin.BaseSeries] basetable
                INNER JOIN [Nucletron].[dbo].[PatientAdmin.Series] seriestable ON basetable.parentOID = seriestable.OID
                INNER JOIN [Nucletron].[dbo].[PatientAdmin.Studies] studytable ON seriestable.parentOID = studytable.OID
                INNER JOIN [Nucletron].[dbo].[PatientAdmin.Patients] patienttable ON studytable.parentOID = patienttable.OID
                WHERE basetable.modality = 'RTPLAN' AND patienttable.patientID = '%s' AND studytable.studyID = '%s'
                ORDER BY basetable.number DESC """ % (patient_ID, study_ID))
    cursor_content = [row for row in cur]
    plans = []
    for row in cursor_content:
        plans.append(' '.join([str(i) for i in list(row)]))
    plans = [plan.replace('/','-') for plan in plans]
    return plans, cursor_content

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
    plan_data_dict['sources']['dwell_positions'] = []
    plan_data_dict['sources']['channel_number'] = []
    plan_data_dict['sources']['cumulative_dwell_times'] = []
    plan_data_dict['channels'] = []
    for channel_number in range(len(tree.xpath('.//Channel'))):
        channel_total_time = float(tree.xpath('.//TotalTime')[channel_number].text)
        temp_dict = {}
        temp_dict['channel_number'] = float(tree.xpath('.//Channel')[channel_number].xpath('.//Number')[0].text)
        temp_dict['channel_time_total'] = float(tree.xpath('.//Channel')[channel_number].xpath('.//TotalTime')[0].text)
        temp_dict['reference_length'] = float(tree.xpath('.//Channel')[channel_number].xpath('.//Length')[0].text)
        temp_dict['step_size'] = float(tree.xpath('.//Channel')[channel_number].xpath('.//SourceApplicatorStepSize')[0].text)
        plan_data_dict['channels'].append(temp_dict)
        for control_point_number in range(len(tree.xpath('.//Channel')[channel_number].xpath('.//BrachyControlPoint'))):
            dwell_relative_time = float(tree.xpath('.//Channel')[channel_number].xpath('.//BrachyControlPoint')[control_point_number].getchildren()[5].text)
            plan_data_dict['sources']['cumulative_dwell_times'].append((dwell_relative_time/100) * channel_total_time)
            plan_data_dict['sources']['channel_number'].append(channel_number)
            plan_data_dict['sources']['dwell_positions'].append(int(float(tree.xpath('.//Channel')[channel_number].xpath('.//BrachyControlPoint')[control_point_number].xpath('.//RelativePosition')[0].text)/5)+1)
            plan_data_dict['sources']['coordinates'].append([0.1*float(i) for i in tree.xpath('.//Channel')[channel_number].xpath('.//BrachyControlPoint')[control_point_number].getchildren()[4:5][0].text.split(':')])
    plan_data_dict['sources']['dwell_times'] = [x - y for x,y in zip(plan_data_dict['sources']['cumulative_dwell_times'][1::2],plan_data_dict['sources']['cumulative_dwell_times'][0::2])]
    plan_data_dict['sources']['coordinates'] = [i for i in plan_data_dict['sources']['coordinates'][0::2]]
    plan_data_dict['sources']['channel_number'] = [i+1 for i in plan_data_dict['sources']['channel_number'][0::2]]
    plan_data_dict['sources']['dwell_positions'] =  [i for i in plan_data_dict['sources']['dwell_positions'][0::2]]
    return plan_data_dict


def fetch_plan(patient_ID, study_ID, plan_number):
    """
    Fetch xml data from OP database
    """
    cur = connect_to_database()
    cur.execute("SELECT channels FROM [Nucletron].[dbo].[BrachyTreatmentPlan.ApplicationSetups] applsetuptable INNER JOIN [Nucletron].[dbo].[PatientAdmin.SeriesItems] seriesitemstable ON applsetuptable.parentOID = seriesitemstable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Series] seriestable ON seriesitemstable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.BaseSeries] basetable ON basetable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Studies] studytable ON seriestable.parentOID = studytable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Patients] patienttable ON studytable.parentOID = patienttable.OID WHERE patienttable.patientID = ? AND studytable.studyID = ? AND basetable.modality = 'RTPLAN' AND basetable.number = ?",(patient_ID,study_ID,plan_number))
    cursor_content = [row for row in cur]
    plan_content =  cursor_content[0][0].decode()
    plan_data_dictionary = process_plan_xml(plan_content)
    plan_data_dictionary['patient_ID'] = patient_ID
    plan_data_dictionary['study_ID'] = study_ID
    plan_data_dictionary['plan_number'] = plan_number
    cur.execute("SELECT settable.prescribedDose FROM [Nucletron].[dbo].[BrachyTreatmentPlan.SettingsSwifts] settable INNER JOIN [Nucletron].[dbo].[PatientAdmin.SeriesItems] seriesitemstable ON settable.parentOID = seriesitemstable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Series] seriestable ON seriesitemstable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.BaseSeries] basetable ON basetable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Studies] studytable ON seriestable.parentOID = studytable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Patients] patienttable ON studytable.parentOID = patienttable.OID WHERE patienttable.patientID = ? AND studytable.studyID = ? AND basetable.modality = 'RTPLAN' AND basetable.number = ?",(patient_ID,study_ID,plan_number))
    cursor_content = [row for row in cur]
    plan_data_dictionary['prescription_Gy'] = cursor_content[0][0]/100

    cur.execute("SELECT sourcetable.referenceAirKermaRate FROM [Nucletron].[dbo].[BrachyTreatmentPlan.Sources] sourcetable INNER JOIN [Nucletron].[dbo].[PatientAdmin.SeriesItems] seriesitemstable ON sourcetable.parentOID = seriesitemstable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Series] seriestable ON seriesitemstable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.BaseSeries] basetable ON basetable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Studies] studytable ON seriestable.parentOID = studytable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Patients] patienttable ON studytable.parentOID = patienttable.OID WHERE patienttable.patientID = ? AND studytable.studyID = ? AND basetable.modality = 'RTPLAN' AND basetable.number = ?",(patient_ID,study_ID,plan_number))
    cursor_content = [row for row in cur]
    air_kerma_rate_constant = 4.082
    plan_data_dictionary['apparent_activity'] = (cursor_content[0][0]/air_kerma_rate_constant)/1000

    cur.execute("SELECT dosetable.description, dosetable.pointCoordinates, targetPrescriptionDose FROM [Nucletron].[dbo].[TreatmentPlan.DoseReferences] dosetable INNER JOIN [Nucletron].[dbo].[PatientAdmin.SeriesItems] seriesitemstable ON dosetable.parentOID = seriesitemstable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Series] seriestable ON seriesitemstable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.BaseSeries] basetable ON basetable.parentOID = seriestable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Studies] studytable ON seriestable.parentOID = studytable.OID INNER JOIN [Nucletron].[dbo].[PatientAdmin.Patients] patienttable ON studytable.parentOID = patienttable.OID WHERE patienttable.patientID = ? AND studytable.studyID = ? AND basetable.modality = 'RTPLAN' AND basetable.number = ?",(patient_ID,study_ID,plan_number))
    cursor_content = [row for row in cur]
    my_POIs = [list(POI) for POI in cursor_content]
    # import ipdb; ipdb.set_trace()
    # calculate_doses_from_plan(plan_data_dictionary, points_of_interest)
    return plan_data_dictionary, my_POIs



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
