from django.test import TestCase
import hashlib
import csv
# def md5(fname):
#     hash_md5 = hashlib.md5()
#     with open(fname, "rb") as f:
#         for chunk in iter(lambda: f.read(4096), b""):
#             hash_md5.update(chunk)
#     return hash_md5.hexdigest()

def md5(fname):
    x = csv.reader(fname)
    return hashlib.md5(str(list(x)).encode('utf-8')).hexdigest()

class SourceDataTestCase(TestCase):
    """
    Checks md5 digests of source data files
    """

    def test_radial_dose(self):
        self.assertEqual(md5(r'myapp\\pyTG43\\source_data\\v2r_ESTRO_radialDose.csv'), '9a48b234c77c8d1336acd9f41448ea91')

    def test_anisotropy_function(self):
        self.assertEqual(md5(r'myapp\\pyTG43\\source_data\\v2r_ESTRO_anisotropyFunction.csv'), 'abec4a4e5e767978cc04f19131160630')

    def test_qa_data(self):
        self.assertEqual(md5(r'myapp\\pyTG43\\source_data\\qa_data.csv'), '8fcf5e4631157628bbc4188f5a751f79')

    def test_qa_data_json(self):
        self.assertEqual(md5(r'myapp\\pyTG43\\source_data\\test_data.json'), '1bfdb8da4dc35e9cb43b9d2cb9f530af')
