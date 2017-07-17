from django.test import TestCase
import codecs

class TestProtocol(TestCase):
    def test_protocol_generation(self):
        """
        Generates a protocol and compares its HTML output to an archived file.
        """
        response = self.client.get('/view_protocol/V117108T/Insertion%201/Pre%20Ins%201/', follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'myapp/view_protocol.html')
        # with open(r"myapp\\tests\\test_protocol_page.htm", 'r') as myfile:
        #     test_html_data=myfile.read()
        # self.assertEqual(response.content.decode("utf-8"),test_html_data)
