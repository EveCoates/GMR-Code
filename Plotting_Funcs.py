# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 11:01:18 2022

@author: Samco
"""

##############################################################################

# importing relevant libraries

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import statistics
import pandas as pd

##############################################################################

def import_all(filename, num_ROIS, num_subROIS):
    
    # reading and appending relevant data
    all_data = []
    file = pd.read_csv(filename)
    all_data.append(file[[" res"]])
    all_data.append([file[" amp"]])
    all_data.append(file[[" FWHM"]])
    all_data.append([file[" off"]])
    all_data.append(file[[" gamma"]])
    all_data.append([file[" r2"]])
    all_data.append(file[[" assym"]])
    
    # averaging every 7th point
    all_data_avr = []
    ROIS = []
    for i in range(7):
        data_avr = np.mean(np.array(all_data[i]).reshape(-1, num_subROIS), axis=1)
        all_data_avr.append([data_avr])
        split = int(len(data_avr) / num_ROIS)
        
        # splitting data by ROI
        array = np.array(all_data_avr[i])
        data = array[0]
        ROI_split = [data[i:i + split] for i in range(0, len(data), split)]
        ROIS.append(ROI_split)
        
    # centering data around 0 by the first value
    for i in range(7):
        for j in range(num_ROIS):
            init = ROIS[i][j][0]
            for k in range(split):
                ROIS[i][j][k] = ROIS[i][j][k] - init                   

    return all_data, ROIS, split

##############################################################################

# for resonance position, pairs up the data into resonance separation

def make_pairs(ROI_data, length, num_pairs, num_ROIS):
    pairs = []
    if num_ROIS == 6:
        for i in range(num_pairs):
            nest = [] 
            for j in range(length):
                nest.append(abs(ROI_data[i][j]) + abs(ROI_data[i + 1][j]))
                if len(nest) == length:
                    pairs.append(nest)  
    elif num_ROIS == 3:
        for i in range(num_ROIS):
            nest = []
            for j in range(length):
                nest.append(ROI_data[i][j])
                if len(nest) == length:
                    pairs.append(nest)
    return pairs
    
##############################################################################

def get_ydata(ROIS, idx, length, num_pairs, num_ROIS): 
    if idx == 0:
        y_data = make_pairs(ROIS, length, num_pairs, num_ROIS)
    else:
        y_data = ROIS
    return y_data

##############################################################################

title_list = ["Resonance Separation over time", "Amplitude over time", "FWHM over time", "Offset over time", "Gamma over time", "r2 over time", "Asymmetry over time"]
yaxis_list = ["Resonance Separation", "Amplitude", "FWHM", "Offset", "Gamma", "r2", "Asymmetry"]
data_dict = {0: title_list, 1: yaxis_list, 2: [3, 6, 6, 6, 6, 6, 6]}

def data_info(idx):
    data_info = []
    for i in range(3):
        data_info.append(data_dict[i][idx])
        
    return data_info

##############################################################################

def replace_outliers(y_data, length, data_sets, threshold):
    
    # replacing blank image peaks with data points before them
    for i in range(3):
        for j in range(length - 1): 
            if y_data[i][j] > threshold:
                y_data[i][j] = y_data[i][j+1]
     
    # removing every 11th data point due to quirk of analysis code
    for i in range(data_sets):
        j = 10
        n = 2
        while j < length:
            if n % 11 == 0:
                y_data[i][j] = y_data[i][j-1] 
                j = j + 1
                n = n + 1
            else:
                y_data[i][j] = y_data[i][j-1]
                j = j + 11
                n = n + 1
    
##############################################################################

# creates a list of x values for plotting

def get_xdata(length, wait_time): 
    x_data = []
    for i in range(length):
        x_data.append(wait_time * i)
    return x_data

##############################################################################

# asks the user for inputs, rounding to the nearest data point

def get_points(num_points, wait_time):
    inputs = plt.ginput(num_points)
    plt.close()
    rounded_points = []
    for i in range(num_points):
        rounded = wait_time * round(inputs[i][0] / wait_time)
        rounded_points.append(rounded)

    rounded_points_array = [rounded_points[i:i + 2] for i in range(0, num_points, 2)]

    return rounded_points_array
    
##############################################################################   
    
# plots the data and labels the figure

def plotter(x_data, y_data, data_range, plot_title, yaxis_title):
    
    for i in range(data_range):
        tag = "plot" + "_" + str(i) 
        plt.plot(x_data, y_data[i], label=tag)
        
        plt.xlabel("time (seconds)")
        plt.ylabel(yaxis_title)
        plt.title(plot_title)
        
        plt.xlim(0, x_data[len(x_data) - 1])
        plt.ylim()
        
        plt.legend()
        plt.grid()    
   
##############################################################################

# recieves pairs of points, finding the average value of the data between them
# as well as the uncertainty (standard deviation)

def find_avr(points, x_data, y_data, data_sets):

    # iterates over each pair of x1 and x2
    avr_array = []
    stdv_array = []
    for i in range(len(points)):
        x1 = x_data.index(points[i][0])
        x2 = x_data.index(points[i][1])
        
        # iterates for each line for given x1 and x2
        for j in range(data_sets):
            pair_avr = np.mean((y_data[j][x1:x2]), axis=0)
            avr_array.append(np.mean(pair_avr))
            stdv_array.append(statistics.stdev((y_data[j][x1:x2])))
           
    # collects and averages the results for each individual line
    averages = []    
    standard_dev = []    
    avr = np.mean(np.array(avr_array).reshape(-1, 3), axis=1)
    avr_stdv = np.mean(np.array(stdv_array).reshape(-1, 3), axis=1)
    
    averages.append(avr)
    standard_dev.append(avr_stdv)
    
    return averages[0], standard_dev[0]
    
##############################################################################

# uses the pandas library to tabulate the final data in a dataframe

def pandas_df(points, avr, num_pairs):
    
    x1_points_list = []
    x2_points_list = []
    avr_list = []
    stdv_list = []
    
    #sorting x1 and x2 data
    for i in range(num_pairs):
        x1_points_list.append(points[i][0])
        x2_points_list.append(points[i][1])
        
    # sorting average and stdv data
    for i in range(2):
        if i == 0:
            for j in range(num_pairs):
                avr_list.append(avr[i][j])
        elif i == 1:
            for j in range(num_pairs):
                stdv_list.append(avr[i][j])
        
    data_dict = {'X1': x1_points_list, 'X2': x2_points_list, 'Average': avr_list, 'STDV': stdv_list}
    print(pd.DataFrame(data=data_dict))
    
##############################################################################

def write_to_memory(points):
    
    # breaks the nested lists into a single list 
    memory_points = []
    for i in range(len(points)):
        memory_points.append(points[i][0])
        memory_points.append(points[i][1])
    
    # writes the single list to memory.txt 
    memory = open('memory.txt', 'w')
    for i in range(len(memory_points)):
        memory.write(str(memory_points[i]))
        memory.write(",")

##############################################################################

def read_from_memory():
    
    #opens memory.txt and returns the first line as a list
    memory_file = pd.read_csv('memory.txt', header=None)
    points_list = memory_file.iloc[0].tolist()
    points_list.pop()
    points_pairs = [points_list[x:x+2] for x in range(0, len(points_list), 2)]
    
    return points_pairs
    
##############################################################################
 
def analysis(choice, num_points, wait_time, x_data, y_data, data_sets, plot_title, yaxis_title):
    
    if choice == 1:
        plotter(x_data, y_data, data_sets, plot_title, yaxis_title)
        inputs = get_points(num_points, wait_time)
        write_to_memory(inputs)
        avr = find_avr(inputs, x_data, y_data, data_sets)
        return inputs, avr
    
    elif choice == 2:
        inputs = read_from_memory()
        avr = find_avr(inputs, x_data, y_data, data_sets)
        return inputs, avr
        
##############################################################################