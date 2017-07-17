from django.test import TestCase
import myapp.omputilities.omp_connect as omp
from myapp.views import fetch_plan

class OMPConnectionTestCase(TestCase):
    """
    Test ability to connect to OMP database
    """
    def test_connection_to_omp(self):
        self.assertTrue(omp.connect_to_db())

class TestGetPatientCases(TestCase):
    """
    Test ability to retrieve list of patient cases
    """
    def test_get_patient_cases(self):
        self.assertTrue(omp.get_patient_cases('V117108T'))

class TestGetPatientPlans(TestCase):
    """
    Test ability to retrieve list of patient plans
    """
    def test_get_patient_plans(self):
        self.assertTrue(omp.get_plans_from_case('V117108T', 'Insertion 1'))

class TestGetPatientRTPlan(TestCase):
    """
    Test ability to fetch BLOB RTPlan file and open it
    """
    def test_get_patient_rtplan(self):
        self.assertTrue(fetch_plan('V117108T', 'Insertion 1','Pre Ins 1'))

class TestGetPatientIDs(TestCase):
    """
    Test ability to retrieve list of patient plans
    """
    def test_get_patient_IDs(self):
        self.assertTrue(omp.get_patient_IDs_and_names('white'))
