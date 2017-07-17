from django import forms
from django.utils.translation import ugettext_lazy as _
from .models import Patient, DVHDump, PatientName

class PatientForm(forms.ModelForm):
    class Meta:
        model = Patient
        fields = ('patient_ID',)

class PatientNameForm(forms.ModelForm):
    class Meta:
        model = PatientName
        fields = ('patient_name',)

class DVHDumpForm(forms.ModelForm):
    class Meta:
        model = DVHDump
        fields = ('dump',)
        labels = {
            'dump': _(''),
        }
