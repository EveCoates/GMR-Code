# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 14:14:00 2022

@author: Samco
"""

### IMPORTING LIBRARIES & FUNCTIONS ###
import matplotlib
matplotlib.use('TkAgg')
from Plotting_Funcs import import_all, data_info, replace_outliers, get_ydata, get_xdata, pandas_df, analysis

### DATA INDEX ###
# 0: Resonance Separation
# 1: Amplitude
# 2: FWHM
# 3: Offset
# 4: Gamma
# 5: r2
# 6: Asymmetry 

### USER INFO INPUT ###
filename = "Sensor_1_Run_1_Results_ALLDATA.csv"
data_idx = 2
num_ROIS = 6
num_subROIS = 7
wait_time = 5
num_points = 6
num_pairs = int(num_points / 2)


# raw data import
all_data = import_all(filename, num_ROIS, num_subROIS)
ROIS = all_data[1][data_idx]
data_len = all_data[2]

# storing final data and info
info_array = data_info(data_idx)
plot_title = info_array[0]
yaxis_title = info_array[1]
data_sets = info_array[2]

y_data = get_ydata(ROIS, data_idx, data_len, num_pairs, num_ROIS)
x_data = get_xdata(data_len, wait_time)

# cleaning anomalous data points
threshold = 300
replace_outliers(y_data, data_len, data_sets, threshold)

# printing basic info
print(filename)
print(yaxis_title, "data")

### CHOICE ###
# 1) select new points
# 2) use previous points
choice = 1
analysis = analysis(choice, num_points, wait_time, x_data, y_data, data_sets, plot_title, yaxis_title)
input_points = analysis[0]
avr = analysis[1]

# displays the final data in a pandas library dataframe
pandas_df(input_points, avr, num_pairs)









