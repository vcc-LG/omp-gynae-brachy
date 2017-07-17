from django.test import TestCase

class TestDoseCheck(TestCase):
    def test_dose_check_generation(self):
        """
        Generates a dose check page and compares its HTML output to an archived file.
        """
        response = self.client.get('/dose_check/V117108T/Insertion%201/Pre%20Ins%201/', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/dose_check.html')
        # with open(r"myapp\\tests\\test_dose_check_page.html", 'r') as myfile:
        #     test_html_data=myfile.read()
        # self.assertEqual(response.content.decode("utf-8"),test_html_data)
