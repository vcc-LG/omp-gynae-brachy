from django.db import models

class Patient(models.Model):
    patient_ID = models.CharField(max_length=200)
    def __str__(self):
        return self.patient_ID

class PatientName(models.Model):
    patient_name = models.CharField(max_length=200)
    def __str__(self):
        return self.patient_name

class DVHDump(models.Model):
    dump = models.TextField()
