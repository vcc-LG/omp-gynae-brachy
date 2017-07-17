from django.test import TestCase

class TestViews(TestCase):
    def test_index(self):
        """
        Is the index page loaded?
        """
        response = self.client.get('/', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/index.html')

    def test_patientid_not_found(self):
        """
        Is the error page loaded with a nonexistent patient ID?
        """
        response = self.client.get('/view_patient/AAA', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/error_cases_not_found.html')

    def test_list_cases(self):
        """
        Is the view_patient page loaded?
        """
        response = self.client.get('/view_patient/V117108T', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/view_patient.html')

    def test_list_ids(self):
        """
        Is the view_ids page loaded?
        """
        response = self.client.get('/view_ids/%25WHITE%25/', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/view_ids.html')

    def test_patientname_not_found(self):
        """
        Is the error page loaded with a nonexistent patient name?
        """
        response = self.client.get('/view_ids/%25ASDASDASDASD%25', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/error_patient_not_found.html')

    def test_list_plans(self):
        """
        Is the view_case page loaded?
        """
        response = self.client.get('/view_case/V117108T/Insertion%201/', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/view_case.html')

    def test_view_plan(self):
        """
        Is the view_plan page loaded?
        """
        response = self.client.get('/view_plan/V117108T/Insertion%201/Pre%20Ins%201/', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/view_plan.html')


    def test_view_dvh_dump(self):
        """
        Is the dvh_dump page loaded?
        """
        response = self.client.get('/dvh_dump', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/dvh_dump.html')

    def test_dvh_dump_error(self):
        """
        Is the dvh_dump page loaded?
        """
        response = self.client.get('/dvh_dump', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/dvh_dump.html')
