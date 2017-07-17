# Oncentra Prostate protocol generator and dose check web app

## Introduction

This application produces radiographer protocols and independent dose check calculations for HDR prostate brachytherapy treatments using data pulled from the Oncentra Prostate database.

## Usage

The Oncentra Prostate laptop must be switched on for this app to work.

Create Special Points in the Dose Evaluation module for use in the dose check. Try to ensure that the points are more than a few mm away from any dwell to avoid significant discrepancies in the dose calculation. 

Activate the virtual environment and start the server:
```
venv\Scripts\activate
(venv) python manage.py runserver
```

And open a browser to: [localhost:8000](localhost:8000)

## Tests

Tests cover the pyTG43 calculation itself as well as the views, protocol generation, and dose check display:

```
(venv) python manage.py test
```

And view test coverage with:

```
(venv) coverage run --source='.' manage.py test
(venv) coverage report
```
