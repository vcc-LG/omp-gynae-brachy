from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^view_patient/(?P<patient_ID>.*)/$', views.view_patient, name='view_patient'),
    url(r'^view_case/(?P<patient_ID>.*)/(?P<case_label>.*)/$', views.view_case, name='view_case'),
    url(r'^view_ids/(?P<patient_name>.*)/$', views.view_ids, name='view_ids'),
    url(r'^view_plan/(?P<patient_ID>.*)/(?P<case_label>.*)/(?P<plan_name>.*)/$', views.view_plan, name='view_plan'),
    url(r'^view_protocol/(?P<patient_ID>.*)/(?P<case_label>.*)/(?P<plan_name>.*)/$', views.view_protocol, name='view_protocol'),
    url(r'^dose_check/(?P<patient_ID>.*)/(?P<case_label>.*)/(?P<plan_name>.*)/$', views.dose_check, name='dose_check'),
    url(r'^dvh_dump', views.dvh_dump, name='dvh_dump'),
]
