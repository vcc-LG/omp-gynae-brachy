"""Perform TG43 calculations"""
from __future__ import division
import sys
import os
import getopt
import numpy as np
from scipy.spatial.distance import pdist
from scipy import interpolate
from lmfit import minimize, Parameters, Parameter, report_fit
import json
from pprint import pprint

np.seterr(divide='raise')    #this suppresses a RunTimeWarning to a more
                             #meaningful exception

import csv
import numpy as np


def read_file(full_path):
    """
    Read in CSV files
    """
    in_file = open(full_path, "r")
    reader = csv.reader(in_file)
    input_data = []
    for row in reader:
        input_data.append(row)
    in_file.close()
    return input_data


class BrachyPlan(object):
    def __init__(self, ds):
        self.ds = ds
        self.applicator = self.ds.ApplicationSetupSequence[0][0x300b, 0x100f].value
        self.points = self.get_poi()
        self.channel_numbers = self.get_channel_numbers()
        self.prescription = float(
            ds.FractionGroupSequence[0].ReferencedBrachyApplicationSetupSequence[0].BrachyApplicationSetupDose)
        self.treatment_model = ds.TreatmentMachineSequence[0].TreatmentMachineName
        self.ref_air_kerma_rate = float(ds.SourceSequence[0].ReferenceAirKermaRate)
        self.channels, self.total_channel_times = self.get_channel_dwell_times()
        self.total_number_dwells = sum([len(i) for i in self.channels])
        self.half_life = float(ds.SourceSequence[0].SourceIsotopeHalfLife)
        self.patient_ID = ds.PatientID
        self.plan_name = ds.RTPlanLabel
        self.applicator_setup_type = ds.ApplicationSetupSequence[0].ApplicationSetupType
        self.reference_lengths = self.get_reference_lengths()
        self.step_sizes = self.get_step_sizes()
        self.air_kerma_rate_constant = 4.082
        self.apparent_activity = self.ref_air_kerma_rate / self.air_kerma_rate_constant * 0.001
        self.channel_names = self.get_channel_names()

    def get_channel_names(self):
        channel_names = []
        for i in range(len(self.ds[0x300f, 0x1000][0].StructureSetROISequence)):
            channel_names.append(self.ds[0x300f, 0x1000][0].StructureSetROISequence[i][0x3006,0x0026].value)
        return channel_names

    def get_reference_lengths(self):
        reference_lengths = []
        for c in range(len(self.ds.ApplicationSetupSequence[0].ChannelSequence)):
            reference_lengths.append(round( float(self.ds.ApplicationSetupSequence[0].ChannelSequence[c].ChannelLength)))
        return reference_lengths

    def get_step_sizes(self):
        step_sizes = []
        for c in range(len(self.ds.ApplicationSetupSequence[0].ChannelSequence)):
            step_sizes.append(round( float(self.ds.ApplicationSetupSequence[0].ChannelSequence[c].SourceApplicatorStepSize),1))
        return step_sizes

    def get_channel_numbers(self):
        return [int(x.SourceApplicatorID) for x in self.ds.ApplicationSetupSequence[0].ChannelSequence]

    def get_poi(self):
        points = []
        for p in self.ds.DoseReferenceSequence:
            points.append(self.Point(p))
        return points

    def get_channel_dwell_times(self):
        channel_dwells = []
        total_channel_times = []
        for c in range(len(self.ds.ApplicationSetupSequence[0].ChannelSequence)):
            total_channel_times.append(round(float(self.ds.ApplicationSetupSequence[0].ChannelSequence[c].ChannelTotalTime),2))
            dwell_weights = []
            dwells_list = []
            dwell_positions = []
            dwell_times = []
            time_factor = self.ds.ApplicationSetupSequence[0].ChannelSequence[c].ChannelTotalTime  / self.ds.ApplicationSetupSequence[0].ChannelSequence[c].FinalCumulativeTimeWeight
            number_of_dwells = int(self.ds.ApplicationSetupSequence[0].ChannelSequence[c].NumberOfControlPoints / 2)
            for i in range(0, len(self.ds.ApplicationSetupSequence[0].ChannelSequence[c].BrachyControlPointSequence), 2):
                d1 = float(
                    self.ds.ApplicationSetupSequence[0].ChannelSequence[c].BrachyControlPointSequence[
                        i].CumulativeTimeWeight)
                d2 = float(
                    self.ds.ApplicationSetupSequence[0].ChannelSequence[c].BrachyControlPointSequence[
                        i + 1].CumulativeTimeWeight)
                dwell_weights.append((d2 - d1))
                dwell_times.append(time_factor*(d2-d1))
                dwells_list.append(self.ds.ApplicationSetupSequence[0].ChannelSequence[c].BrachyControlPointSequence[
                        i])
                dwell_positions.append(int(1+self.ds.ApplicationSetupSequence[0].ChannelSequence[c].BrachyControlPointSequence[i].ControlPointRelativePosition/self.ds.ApplicationSetupSequence[0].ChannelSequence[c].SourceApplicatorStepSize))
            # dwell_times = [(total_channel_time / number_of_dwells) * x for x in dwell_weights]
            dwells = []
            for i in range(len(dwells_list)):
                dwells.append(self.Dwell(dwells_list[i],dwell_times[i],dwell_weights[i], dwell_positions[i]))
            channel_dwells.append(dwells)
        return channel_dwells, total_channel_times

    class Point(object):
        def __init__(self, ds_sequence):
            self.name = ds_sequence.DoseReferenceDescription
            self.coords = [float(x) for x in ds_sequence.DoseReferencePointCoordinates]
            self.dose = float(ds_sequence.TargetPrescriptionDose)

    class Dwell(object):
        def __init__(self, control_sequence, d_time, d_weight, d_pos):
            self.coords = [float(x) for x in control_sequence.ControlPoint3DPosition]
            self.time_weight = d_weight
            self.dwell_time = d_time
            self.dwell_position = d_pos


class PointComparison(object):
    def __init__(self, point_name, omp_dose, pytg43_dose):
        self.point_name = point_name
        self.omp_dose = omp_dose
        self.pytg43_dose = pytg43_dose
        self.abs_difference = omp_dose-pytg43_dose
        self.percentage_difference = 100*((self.omp_dose/self.pytg43_dose)-1)

class SourcePosition:
    """
    Class to hold source description data
    """

    def __init__(
            self,
            x,
            y,
            z,
            apparent_activity,
            dwell_time,
            Sk,
            dose_rate_constant,
            L,
            t_half):
        self.x = x
        self.y = y
        self.z = z
        self.apparent_activity = apparent_activity
        self.dwellTime = dwell_time
        self.coords = [x, y, z]
        self.Sk = Sk  # air kerma strength cGy.cm2/hr
        # dose rate constant cGy/(h.U)
        self.dose_rate_constant = dose_rate_constant
        self.L = L
        self.t_half = t_half

def make_source_trains(source_class):
    source_train = []
    for channel in source_class.channels:
        for source in channel:
            source_train.append(
                SourcePosition(
                    x=source.coords[0] / 10,  # lateral
                    y=source.coords[2] / 10,  # sup-inf
                    z=source.coords[1] / 10,  # ant-post
                    apparent_activity=10,
                    dwell_time=source.dwell_time,
                    Sk=source_class.ref_air_kerma_rate,
                    dose_rate_constant=1.108,
                    L=0.35,
                    t_half=source_class.half_life))
    return source_train


def make_radial_dose(radial_dose_raw):
    """
    Create radial dose function from raw input data
    """
    r_cm = []
    gL = []
    for i in range(1, len(radial_dose_raw)):
        r_cm.append(float(radial_dose_raw[i][0]))
        gL.append(float(radial_dose_raw[i][1]))
    return RadialDoseClass(r_cm, gL)


class RadialDoseClass:
    """
    Class to hold radial dose function
    """

    def __init__(self, r_cm, gL):
        self.r_cm = r_cm
        self.gL = gL


def make_anisotropy_function(anisotropy_function_raw):
    """
    Create anisotropy function from raw input data
    """
    A = [[row for row in anisotropy_function_raw[i][0]]
         for i in range(2, len(anisotropy_function_raw))]
    theta = [float("".join(A[i])) for i in range(len(A))]
    B = anisotropy_function_raw[1][1:]
    r_cm = [float(i) for i in B]
    C = [[row for row in anisotropy_function_raw[i][1:]]
         for i in range(2, len(anisotropy_function_raw))]
    F = np.zeros([len(C), len(C[0])])
    for i in range(len(C)):
        for j in range(len(C[i])):
            if i != -1:
                try:
                    F[i][j] = float(C[i][j])
                except ValueError:
                    F[i][j] = None
            elif i == -1:
                F[i][j] = None
    return AnisotropyFunctionClass(r_cm, theta, F)


class AnisotropyFunctionClass:
    """
    Class to hold anisotropy function
    """

    def __init__(self, r_cm, theta, F):
        self.r_cm = r_cm
        self.theta = theta
        self.F = F


def find_nearest(array, value):
    """
    Find the index of the closest value in an array
    """
    idx = (np.abs(array - value)).argmin()
    return array[idx]

radialDose = make_radial_dose(
    read_file(r'myapp\\pyTG43\\source_data\\v2r_ESTRO_radialDose.csv'))
anisotropyFunc = make_anisotropy_function(
    read_file(r'myapp\\pyTG43\\source_data\\v2r_ESTRO_anisotropyFunction.csv'))

class PointPosition:
    """
    Class to hold special point location
    """

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.coords = [x, y, z]


class SourcePosition:
    """
    Class to hold source description data
    """

    def __init__(
            self,
            x,
            y,
            z,
            apparent_activity,
            dwell_time,
            Sk,
            dose_rate_constant,
            L,
            t_half):
        self.x = x
        self.y = y
        self.z = z
        self.apparent_activity = apparent_activity
        self.dwellTime = dwell_time
        self.coords = [x, y, z]
        self.Sk = Sk  # air kerma strength cGy.cm2/hr
        # dose rate constant cGy/(h.U)
        self.dose_rate_constant = dose_rate_constant
        self.L = L
        self.t_half = t_half


def get_geometry_function(my_source, my_point):
    """
    Calculate the geometry function
    x = out of plane, y = vertical, z = source axis
    """
    rRef = 1  # cm
    thetaRef = np.pi / 2  # 90 degrees
    betaRef = 2 * np.arctan((my_source.L / 2) / rRef)
    GlRef = betaRef / (my_source.L * rRef * np.sin(thetaRef))

    Gl = []
    R2 = pdist([[my_source.x, my_source.y, my_source.z],
                [my_point.x, my_point.y, my_point.z - (my_source.L / 2)]])
    R1 = pdist([[my_source.x, my_source.y, my_source.z],
                [my_point.x, my_point.y, my_point.z + (my_source.L / 2)]])
    R = pdist([[my_source.x, my_source.y, my_source.z],
               [my_point.x, my_point.y, my_point.z]])

    theta1 = np.arccos((my_point.z - my_source.z + (my_source.L / 2)) / R1)
    theta2 = np.arccos((my_point.z - my_source.z - (my_source.L / 2)) / R2)

    theta = np.arccos((my_point.z - my_source.z) / R)

    if theta == 0 or theta == np.pi:
        Gl = 1 / (R ** 2 - (my_source.L ** 2 / 4))
    else:
        beta = np.abs(theta2 - theta1)
        Gl = beta / (my_source.L * R * np.sin(theta))

    Glout = Gl / GlRef
    return Glout


def nan_helper(y):
    """
    Return the indexes of NaNs in list
    """
    return np.isnan(y), lambda z: z.nonzero()[0]


def log_interp(xdata, ydata, xnew):
    """
    Perform log linear interpolation
    """
    xdata_np = np.asarray(xdata)
    ydata_np = np.asarray(ydata)
    try:
        logx = np.log(xdata)
    except FloatingPointError:
        logx = np.log(xdata_np.clip(min=0.0000000001))
    try:
        logy = np.log(ydata)
    except FloatingPointError:
        logy = np.log(ydata_np.clip(min=0.0000000001))
    # logy = np.log(ydata)
    return np.exp(np.interp(np.log(xnew), logx, logy))


def linear_interp_2d(xdata, ydata, zdata, xnew, ynew):
    """
    Perform linear 2D interpolation
    """
    f = interpolate.interp2d(xdata, ydata, zdata, kind='linear')
    znew = f(xnew, ynew)
    return znew

def fcn2min(params, x, data):
    """
    For exponential curve fitting
    """
    pw = params['pw'].value
    adj1 = params['adj1'].value
    adj2 = params['adj2'].value
    model = adj1 * np.power(x + adj2, pw)
    return model - data

def get_radial_dose(radial_dose_in, dwell_in, point_in):
    """
    Calculate the radial dose function value
    """
    R = pdist([[dwell_in.x, dwell_in.y, dwell_in.z],
               [point_in.x, point_in.y, point_in.z]])
    if R in radial_dose_in.r_cm:
        return radial_dose_in.gL[radial_dose_in.r_cm.index(R)]
    elif R < min(radial_dose_in.r_cm):
        nrVal = find_nearest(np.array(radial_dose_in.r_cm), R)
        return radial_dose_in.gL[radial_dose_in.r_cm.index(nrVal)]
    elif R > max(radial_dose_in.r_cm):
        pw = 2
        params = Parameters()
        params.add('pw', value=pw, vary=False)
        params.add('adj1', value=1)
        params.add('adj2', value=1)
        xf = radial_dose_in.r_cm[-2:]
        yf = radial_dose_in.gL[-2:]
        result = minimize(fcn2min, params, args=(xf, yf))
        report_fit(result.params)
        adj1 = result.params['adj1']
        adj2 = result.params['adj2']
        next_y = adj1 * np.power(R + adj2, pw)
        return next_y
    else:
        return log_interp(radial_dose_in.r_cm, radial_dose_in.gL, R)


def get_anisotropy_function(anisotropy_function, my_source, my_point):
    """
    Calculate the anisotropy function value
    """
    R = pdist([[my_source.x, my_source.y, my_source.z],
               [my_point.x, my_point.y, my_point.z]])
    theta = np.degrees(np.arccos((my_point.z - my_source.z) / R))
    # print "R = %.3f, theta = %.2f"%(R,theta)
    if R in anisotropy_function.r_cm and theta in anisotropy_function.theta:
        return anisotropy_function.F[anisotropy_function.theta.index(theta)][anisotropy_function.r_cm.index(R)]
    elif R > max(anisotropy_function.r_cm) or R < min(anisotropy_function.r_cm) or theta > max(anisotropy_function.theta) or theta < min(anisotropy_function.theta):
        nrValR = find_nearest(np.array(anisotropy_function.r_cm), R)
        nrValtheta = find_nearest(np.array(anisotropy_function.theta), theta)
        return anisotropy_function.F[
            anisotropy_function.theta.index(nrValtheta)][anisotropy_function.r_cm.index(nrValR)]
    else:
        return linear_interp_2d(
            anisotropy_function.r_cm,
            anisotropy_function.theta,
            anisotropy_function.F,
            R,
            theta)


class DosePointClass:
    """
    Class to hold relevant values for dose point
    """

    def __init__(
            self,
            my_source,
            my_point,
            radial_dose_value,
            anisotropy_function_value,
            geometry_function_value,
            dose_rate_out,
            dose_total_out):
        self.my_source = my_source
        self.my_point = my_point
        self.radial_dose_value = radial_dose_value
        self.anisotropy_function_value = anisotropy_function_value
        self.geometry_function_value = geometry_function_value
        self.dose_rate_out = dose_rate_out
        self.dose_total_out = dose_total_out

    def print_values(self):
        """
        Method to print out helpful dose point descriptors
        """
        print("Source @ %s" % self.my_source.coords)
        print("Point @ %s" % self.my_point.coords)
        print("Dwell time = %.2f s" % self.my_source.dwellTime)
        print("Air kerma strength Sk = %.2f U" % self.my_source.Sk)
        print("Apparent activity Aapp = %.2f Ci" % self.my_source.Aapp)
        print("Radial dose g(r) = %.3f" % self.radial_dose_value)
        print(
            "Anisotropy function F(r,theta) = %.3f" %
            self.anisotropy_function_value)
        print(
            "Geometry function G(r,theta) = %.3f" %
            self.geometry_function_value)
        print("Dose rate = %.3f Gy/h" % self.dose_rate_out)
        print("Total dose = %.3f Gy" % self.dose_total_out)

    def print_dose(self):
        """
        Method to print out total dose
        """
        print("%.3f" % self.dose_total_out)


def calculate_my_dose(my_source, my_point, anisotropy_function, radial_dose_function):
    """
    Calculate the total dose at a point
    """
    radial_dose_val = get_radial_dose(
        radial_dose_function, my_source, my_point)
    anisotropy_func_val = get_anisotropy_function(
        anisotropy_function, my_source, my_point)
    geometry_func_val = get_geometry_function(my_source, my_point)
    dose_rate_out = my_source.Sk * my_source.dose_rate_constant * geometry_func_val * \
                    anisotropy_func_val * radial_dose_val * (1 / 100)
    dose_total_out = dose_rate_out * (my_source.dwellTime / (60 * 60))
    #import ipdb; ipdb.set_trace()
    return DosePointClass(
        my_source,
        my_point,
        radial_dose_val,
        anisotropy_func_val,
        geometry_func_val,
        dose_rate_out,
        dose_total_out)


def make_special_points(special_points_raw):
    """
    Feed raw data into special points class
    """
    x_points = []
    y_points = []
    z_points = []
    for i in range(1, len(special_points_raw)):
        x_points.append(float(special_points_raw[i][0]))
        y_points.append(float(special_points_raw[i][1]))
        z_points.append(float(special_points_raw[i][2]))
    return SpecialPointsClass(x_points, y_points, z_points)


class SpecialPointsClass:
    """
    Create Class of special points
    """

    def __init__(self, x_points, y_points, z_points):
        self.xPoints = x_points
        self.yPoints = y_points
        self.zPoints = z_points
        self.numSpecialPoints = len(x_points)

    def print_special_points(self):
        for i in range(len(self.xPoints)):
            print("x = %.1f, y = %.1f, z = %.1f" %
                  (self.xPoints, self.yPoints, self.zPoints))


"""x = out of plane, y = vertical, z = source axis"""
"""Distances in cm"""


def make_source_trains(source_class):
    source_train = []
    for channel in source_class.channels:
        for source in channel:
            source_train.append(
                SourcePosition(
                    x=source.coords[0],  # lateral
                    y=source.coords[1],  # sup-inf
                    z=source.coords[2],  # ant-post
                    apparent_activity=10,
                    dwell_time=source.dwell_time,
                    Sk=source_class.ref_air_kerma_rate,
                    dose_rate_constant=1.108,
                    L=0.35,
                    t_half=source_class.half_life))
    return source_train


def calculate_dose(source_train_in, poi_in):
    dose = 0
    # import ipdb; ipdb.set_trace()
    my_point = PointPosition(
        poi_in[0],  # lateral
        poi_in[1],  # sup-inf
        poi_in[2])  # ant-post
    for dwell in source_train_in:
        my_dose = calculate_my_dose(
            dwell,
            my_point,
            anisotropyFunc,
            radialDose)
        dose += my_dose.dose_total_out
    # import ipdb; ipdb.set_trace()
    return dose.tolist()[0]

def open_json_and_calc(file_in):
    with open(file_in) as data_file:
        data = json.load(data_file)
    # pprint(data['sources'])
    my_source_train = []
    points_of_interest = []
    for i in range(len(data['sources']['coordinates'])):
        my_source_train.append(SourcePosition(x=data['sources']['coordinates'][i][0],
                                              y=data['sources']['coordinates'][i][1],
                                              z=data['sources']['coordinates'][i][2],
                                              dwell_time=data['sources']['dwell_times'][0][i],
                                              apparent_activity=10,
                                              Sk=40820,
                                              dose_rate_constant=1.109,
                                              L=0.36,
                                              t_half=73.83))
    for point in data['POIs']:
        my_dose = calculate_dose(my_source_train, point)
        # print("Dose at point {0} = {1:.4f} Gy".format(point, my_dose))
    return my_dose

def main(argv):
    input_file = ''
    try:
        opts, args = getopt.getopt(argv, "hi:", ["input_file="])
    except getopt.GetoptError:
        print('pyTG43.py -i <input_file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('pyTG43.py -i <input_file>')
            sys.exit()
        elif opt in ("-i", "--input_file"):
            input_file = arg
    if os.path.isfile(input_file):
        open_json_and_calc(input_file)
    else:
        if not os.path.isfile(input_file):
            print("Error: Invalid file path")

if __name__ == "__main__":
    main(sys.argv[1:])
