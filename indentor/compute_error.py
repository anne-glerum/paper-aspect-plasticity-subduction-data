#!/bin/python

from io import StringIO
import numpy as np
import sys

# This program compares ASPECT's point value velocity and pressure output 
# for the indentor benchmark with an analytical
# solution at the same points and reports the percentage error. 

# read in the file name of the predicted velocities to compare with GPS
data_dir = raw_input("Enter directory name of predicted velocities file: ")
data_file = data_dir+'/point_values.txt'
print 'Comparing ', data_file, ' with analytical solution.'

# set the names of the columns we wish to use from the point_values.txt file
required_names = ['evaluation_point_x', 'evaluation_point_y', 'velocity_x', 'velocity_y', 'pressure']

# read the point_values.txt file and only select the requested columns
# default dtype is float
aspect_data = np.genfromtxt(data_file, names = True, delimiter = ' ', usecols=(required_names), dtype = float)

# compute the magnitude of the velocity
aspect_vel_all = np.sqrt(aspect_data['velocity_x']*aspect_data['velocity_x']+aspect_data['velocity_y']*aspect_data['velocity_y'])
aspect_vel = np.zeros((3))
aspect_vel [0] = aspect_vel_all[0]
aspect_vel [1] = aspect_vel_all[1]
aspect_vel [2] = aspect_vel_all[7]

# get the pressure
aspect_p = np.zeros((3))
aspect_p[0] = aspect_data['pressure'][5]
aspect_p[1] = aspect_data['pressure'][2]
aspect_p[2] = aspect_data['pressure'][6]

# the number of points in the file
n_points = float(aspect_vel.shape[0])

# the magnitude of the velocity of the analytical solution
analytical_vel = np.zeros((3))
analytical_vel[0] = 0.74246212024
analytical_vel[1] = 0.74246212024
analytical_vel[2] = 1.05

# the magnitude of the pressure of the analytical solution
analytical_p = np.zeros((3))
analytical_p[0] = 1.
analytical_p[1] = 4.141592654
analytical_p[2] = 1.

# the number of points in the file
n_points_analytical = float(analytical_vel.shape[0])

# make sure both files have the same number of points
if (not(n_points == n_points_analytical)) :
  sys.exit('ASPECT file and analytical solution have a different number of data points: stopping execution.')

# compute the absolute magnitude difference
diff_vel = np.absolute(aspect_vel - analytical_vel)/analytical_vel
diff_vel_mean = np.mean(diff_vel)
diff_p = np.absolute(aspect_p - analytical_p)/analytical_p
diff_p_mean = np.mean(diff_p)

# print mean velocity error
print 'mean vel error:', format(diff_vel_mean*100.0, '0.3f'),  '%'
print 'vel error K:', diff_vel[0]*100.0
print 'vel error L:', diff_vel[1]*100.0
print 'vel error ID:', diff_vel[2]*100.0
# print mean pressure error
print 'mean pressure error:', format(diff_p_mean*100.0, '0.3f'),  '%'
print 'P error H:', diff_p[0]*100.0
print 'P error I:', diff_p[1]*100.0
print 'P error J:', diff_p[2]*100.0

