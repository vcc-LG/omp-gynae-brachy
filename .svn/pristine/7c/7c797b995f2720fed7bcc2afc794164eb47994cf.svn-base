from django.test import TestCase
import myapp.pyTG43.pyTG43 as pyTG43
# import os.path, sys
# sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import csv
import json


class pyTG43TestCase(TestCase):
    """
    This test simply looks at the accuracy of the pyTG43 calculation
    with reference to a QA data set provided by ESTRO. I have set the
    tolerance to be < 0.2 % different to the QA data.
    """
    def test_pyTG43(self):
        qa_file = r'myapp\\pyTG43\\source_data\\qa_data.csv'
        with open(qa_file, 'r') as f:
            reader = csv.reader(f)
            qa_list = list(reader)
        z_list = range(7,-8,-1)
        x_list = range(0,8)
        y_list = range(0,8)
        qa_list_norm_point = float(qa_list[10][5])
        poi_coord = [0,1,0]
        json_dict = {}
        json_dict['sources']={}
        json_dict['sources']['coordinates'] = [[0,0,0]]
        json_dict['sources']['dwell_times'] = [[10]]
        json_dict['POIs']=[poi_coord]
        file_out_path = r'myapp\\pyTG43\\source_data\\test_data.json'
        with open(file_out_path, 'w') as outfile:
            json.dump(json_dict, outfile)
        pytg43_norm_dose = pyTG43.open_json_and_calc(file_out_path)
        for i in range(1,len(qa_list)):
            for j in range(1,len(qa_list[0])):
                qa_point_in = float(qa_list[i][j])
                qa_ratio = qa_point_in/qa_list_norm_point

                z_coord = float(qa_list[i][0])
                y_coord =  float(qa_list[0][j])
                x_coord = 0
                poi_coord = [x_coord,y_coord,z_coord]
                if y_coord != 0 and z_coord != 0:
                    json_dict = {}
                    json_dict['sources']={}
                    json_dict['sources']['coordinates'] = [[0,0,0]]
                    json_dict['sources']['dwell_times'] = [[10]]
                    json_dict['POIs']=[poi_coord]
                    file_out_path = r'myapp\\tests\\testdata.json'
                    with open(file_out_path, 'w') as outfile:
                        json.dump(json_dict, outfile)
                    pytg43_dose_out = pyTG43.open_json_and_calc(file_out_path)
                    pytg43_ratio = pytg43_dose_out/pytg43_norm_dose
                    self.assertTrue(100*abs(1-(pytg43_ratio/qa_ratio)) < 0.2)
