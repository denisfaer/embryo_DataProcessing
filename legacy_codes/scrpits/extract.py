""" import libraries """

import csv
import os
import numpy as np 
import pickle
    
""" configs """

use_intensity = 'MeanIntensity'
use_method = 'nowarp'

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
     
""" pickle data """

pickle_file = open(os.path.join(main_directory, 'raw_data.pkl'), 'wb')
pickle.dump(all_data_raw, pickle_file)
pickle_file.close()

pickle_file = open(os.path.join(main_directory, 'track_data.pkl'), 'wb')
pickle.dump(all_tracks, pickle_file)
pickle_file.close()