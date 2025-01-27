""" import libraries """

import csv
import os
import numpy as np 
import pickle
import math
import matplotlib.pyplot as plt
    
""" configs """

use_intensity = 'MeanIntensity'
use_method = 'nowarp'

fps = 30 # transformation constant of X min between frames

raw_plots = True # create raw intensity value plots?
raw_plot_scatter = True # create raw short vs raw long intensity scatter?
raw_hist_norm = True # create raw short and long normalized by raw histone intensity plots?
raw_plot_scatter_hist = True # create raw short vs raw long intensity scatter where raw histone intensity determines marker size?

raw_ratio = True # create raw short over raw long long intensity plots?

ind_max_norm_hist = True # create individually max-normalized intensity plots (raw values divided by highest within that reporter in that lineage)?
ind_max_norm = True # create individually max-normalized reporter plots (raw values divided by highest within that reporter in that lineage)?

glo_max_norm_hist = True # create globally max-normalized intensity plots (raw values divided by highest within that reporter in that stack)?
glo_max_norm = True # create globally max-normalized reporter plots (raw values divided by highest within that reporter in that stack)?
glo_max_scatter = True # create globally-normalized short vs long intensity scatter?
glo_max_scatter_strict = True # create globally-normalized short vs long intensity scatter with axes locked at 0 to 1?

max_norm = ind_max_norm_hist or ind_max_norm or glo_max_norm_hist or glo_max_norm or glo_max_scatter or glo_max_scatter_strict
minmax_norm = False

marker_size_scaling = 2 # polinomial power of marker size scale
print_plot_start = True # print a message when a plot has started being made [troubleshooting boolean]


""" functions """

def lineage(frame, idi): # returns a string of frame_IDs for a given cell using the edges from the lineage graph
    found = False
    i = 0
    while not found:
        if frame < start_frame:
            mother = ''
            found = True
        elif i > len(track):
            print('LINEAGE_ERROR')
            mother = ''
            found = True
        elif track[i][1][0] == frame and track[i][1][1] == idi:
            mother = str(track[i][1][0]) + '_' + str(track[i][1][1]) + '<' + lineage(track[i][0][0] , track[i][0][1])
            found = True
        else:
            i += 1
        
    return mother

def transform(string): # transforms the lineage() string output to a dataframe
    out = []

    while len(string) > 0:
        sub = string[0:string.index('<')]
        out.append([int(sub[0:sub.find('_')]), int(sub[sub.find('_')+1:len(sub)])])
        string = string.replace(sub + '<', '')
    
    return np.sort(np.array(out), axis = 0)

def find (frame, idi): # returns the short and long camera intensities for a given frame_ID; [-1;-1] if not found
    out = [-1, -1]
    
    for i in range(len(short)):
        if frame == short[i][0] and idi == short[i][1]:
            out = [short[i][2], long[i][2]]
            break
    
    return out

def express(track, interpol = False): # reconstructs the intensity levels of a lineage based of transform() as [frame, ID, short, long, histone]; interpol is the interpolation boolean
    out = []
    
    for t in track:
        temp = find(t[0], t[1])
        temph = -1
        for i in range(len(histone)):
            if t[0] == histone[i][0] and t[1] == histone[i][1]:
                temph = histone[i][2]
                break
        if interpol:
            out.append([t[0], t[1], temp[0], temp[1], temph])
        elif temp[0] > 0:
            out.append([t[0], t[1], temp[0], temp[1], temph])
    
    return np.array(out)

def terminus(track, frame): # returns a list of cells at target frame
    out = []
    
    for t in track:
        if t[1][0] == frame:
            out.append(t[1])
    
    return out

""" start analysis """

main_directory = os.getcwd()
# find all stack data folders: must be named stackXYZ with XYZ being the unique name (not character limited) and NOT contain a file type (dot symbol)
stacks = [f for f in os.listdir(main_directory) if ('stack' in f)and(not '.' in f)]

all_data_raw = [] # master list containing raw import data for all stacks
all_tracks = [] # master list to contain all tracks

""" process individual stacks """

for stack in stacks:
    stack_directory = os.path.join(main_directory, stack)
    
    # grab the lineage data from the graph - for now I manually copied edges from the .mat to a .csv file since I couldn't make .mat reading to work
    graph = [f for f in os.listdir(stack_directory) if ('graph' in f)and(f.endswith('.csv'))]
    track = [] # contains [mother, daughter] cell [frame, ID] data
    with open(os.path.join(stack_directory, graph[0]), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            s1 = row[0].replace("'", '')
            s1_1 = int(s1[0:s1.find('_')])
            s1_2 = int(s1[s1.find('_')+1 : len(s1)])
            s2 = row[1].replace("'", '')
            s2_1 = int(s2[0:s2.find('_')])
            s2_2 = int(s2[s2.find('_')+1 : len(s2)])
            track.append([[s1_1, s1_2], [s2_1, s2_2]])
    
    """ extract intensity measurements """
    
    histone = []
    with open(os.path.join(stack_directory, 'extract_histone.csv'), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        header = True
        for row in reader:
            if header:
                frame = row.index('Frame')
                ID = row.index('ID')
                intense = row.index(use_intensity + '_' + use_method)
                header = False
            else:
                histone.append([int(row[frame]), int(row[ID]), float(row[intense])])
    histone = np.array(histone)
    
    short = []
    with open(os.path.join(stack_directory, 'extract_short.csv'), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        header = True
        for row in reader:
            if header:
                frame = row.index('Frame')
                ID = row.index('ID')
                intense = row.index(use_intensity + '_' + use_method)
                header = False
            else:
                short.append([int(row[frame]), int(row[ID]), float(row[intense])])
    short = np.array(short)
    
    long = []
    with open(os.path.join(stack_directory, 'extract_long.csv'), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        header = True
        for row in reader:
            if header:
                frame = row.index('Frame')
                ID = row.index('ID')
                intense = row.index(use_intensity + '_' + use_method)
                header = False
            else:
                long.append([int(row[frame]), int(row[ID]), float(row[intense])])
    long = np.array(long)
    
    """ format and save stack data """
    
    # set tracking limits to frames with reporter data
    start_frame = max(min(histone[:,0].astype(int)), min(short[:,0].astype(int)), min(long[:,0].astype(int)))
    end_frame = min(max(histone[:,0].astype(int)), max(short[:,0].astype(int)), max(long[:,0].astype(int)))
    
    # get the list of cells at the last frame
    cells = terminus(track, end_frame)
    
    # save raw imports
    all_data_raw.append([stack, short, long, histone])
    
    # compute and save reconstructed terminal lineage raw intensities as individual timeseries
    stack_data = []
    for cell in cells:
        stack_data.append([cell, express(transform(lineage(cell[0], cell[1])))])
    all_tracks.append([stack, len(stack_data), stack_data])
    
    """ make plots """
    
    # calculate the dimensions of a minimal rectangular grid that can contain all plots
    N = math.ceil(math.sqrt(len(cells)))
    M = (len(cells) // N) + 1
    
    # plot all intensities' raw measurements
    if raw_plots:
        if print_plot_start: 
            print(f'{raw_plots = }')
        for i in range(len(cells)):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,4], c='Black')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,2], c='Blue')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,3], c='Red')
        plt.show()
    
    # plot last frame short-vs-long raw intensity scatter
    if raw_plot_scatter:
        if print_plot_start: 
            print(f'{raw_plot_scatter = }')
        for i in range(len(cells)):
            plt.scatter(stack_data[i][1][len(stack_data[i][1])-1, 2], stack_data[i][1][len(stack_data[i][1])-1, 3], s=(stack_data[i][1][len(stack_data[i][1])-1, 4]**marker_size_scaling)/(10**(marker_size_scaling-1)), c='Black')
        plt.show()
    
    # plot histone-normalized [within each frame] raw reporter intensities
    if raw_hist_norm:
        if print_plot_start: 
            print(f'{raw_hist_norm = }')
        for i in range(len(cells)):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,2]/stack_data[i][1][:,4], c='Blue')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,3]/stack_data[i][1][:,4], c='Red')
        plt.show()
    
    # plot last frame short-vs-long raw intensity scatter with histone raw intensity determining size
    if raw_plot_scatter_hist:
        if print_plot_start: 
            print(f'{raw_plot_scatter_hist = }')
        for i in range(len(cells)):
            plt.scatter(stack_data[i][1][len(stack_data[i][1])-1, 2], stack_data[i][1][len(stack_data[i][1])-1, 3], s=stack_data[i][1][len(stack_data[i][1])-1, 4]/10, c='Black')
        plt.show()
    
    # plot short/long intensity raw intensity ratios
    if raw_ratio:
        if print_plot_start: 
            print(f'{raw_ratio = }')
        for i in range(len(cells)):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,2]/stack_data[i][1][:,3], c='Orange')
        plt.show()
    
    # plot individually max-normalized intensities
    if ind_max_norm_hist:
        if print_plot_start: 
            print(f'{ind_max_norm_hist = }')
        for i in range(len(cells)):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,4]/max(stack_data[i][1][:,4]), c='Black')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,2]/max(stack_data[i][1][:,2]), c='Blue')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,3]/max(stack_data[i][1][:,3]), c='Red')
        plt.show()
    
    # plot individually max-normalized reporters
    if ind_max_norm:
        if print_plot_start: 
            print(f'{ind_max_norm = }')
        for i in range(len(cells)):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,2]/max(stack_data[i][1][:,2]), c='Blue')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,3]/max(stack_data[i][1][:,3]), c='Red')
        plt.show()
    
    # compute global max values if needed
    if max_norm or minmax_norm:
        black_max = -1
        blue_max = -1
        red_max = -1
        for i in range(len(cells)):
            blue_max = max(blue_max, max(stack_data[i][1][:,2]))
            red_max = max(red_max, max(stack_data[i][1][:,3]))
            black_max = max(black_max, max(stack_data[i][1][:,4]))
    
    # plot globally max-normalized intensities
    if glo_max_norm_hist:
        if print_plot_start: 
            print(f'{glo_max_norm_hist = }')
        for i in range(len(cells)):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,4]/black_max, c='Black')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,2]/blue_max, c='Blue')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,3]/red_max, c='Red')
        plt.show()
    
    # plot globally max-normalized reporters
    if glo_max_norm:
        if print_plot_start: 
            print(f'{glo_max_norm = }')
        for i in range(len(cells)):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,2]/blue_max, c='Blue')
            plt.plot(stack_data[i][1][:,0], stack_data[i][1][:,3]/red_max, c='Red')
        plt.show()
    
    # plot the terminal frame globally max-normallized short vs long reporters
    if glo_max_scatter:
        if print_plot_start: 
            print(f'{glo_max_scatter = }')
        for i in range(len(cells)):
            plt.scatter(stack_data[i][1][len(stack_data[i][1])-1, 2]/blue_max, stack_data[i][1][len(stack_data[i][1])-1, 3]/red_max, s=10*stack_data[i][1][len(stack_data[i][1])-1, 4]/black_max, c='Black')
        plt.show()
    
    # plot the strict (axis from 0 to 1) terminal frame globally max-normallized short vs long reporters
    if glo_max_scatter_strict:
        if print_plot_start: 
            print(f'{glo_max_scatter_strict = }')
        for i in range(len(cells)):
            plt.scatter(stack_data[i][1][len(stack_data[i][1])-1, 2]/blue_max, stack_data[i][1][len(stack_data[i][1])-1, 3]/red_max, s=10*stack_data[i][1][len(stack_data[i][1])-1, 4]/black_max, c='Black')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.show()
     
""" pickle data """

pickle_file = open(os.path.join(main_directory, 'raw_data.pkl'), 'wb')
pickle.dump(all_data_raw, pickle_file)
pickle_file.close()

pickle_file = open(os.path.join(main_directory, 'track_data.pkl'), 'wb')
pickle.dump(all_tracks, pickle_file)
pickle_file.close()