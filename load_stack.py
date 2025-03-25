""" import libraries """

import csv
import os
import numpy as np

""" supporting functions """

def frame_IDs(given_track, frame): # returns a list of cells in a track at a chosen frame
    out = []
    
    for t in given_track:
        if t[1][0] == frame:
            out.append(t[1])
    
    return out


def lineage(given_track, frame, ID): # recursively reconstructs a cell's lineage returning a string of "frame_ID" separated by "<"
    found = False
    mother = ['',[]]
    i = 0
    
    while not found and i < len(given_track):
        if (given_track[i][1][0] == frame) and (given_track[i][1][1] == ID):
            lin = lineage(given_track, given_track[i][0][0], given_track[i][0][1])
            mother[0] = str(given_track[i][1][0]) + '_' + str(given_track[i][1][1]) + '<' + lin[0]
            mother[1] = lin[1]
            if given_track[i][1] in div:
                mother[1].append(given_track[i][1])
            found = True
        else:
            i = i + 1
      
    return mother


def find(frame, idi): # returns the short and long camera intensities for a given frame_ID; [-1;-1] if not found
    out = [-1, -1]
    
    for i in range(len(short)):
        if frame == short[i][0] and idi == short[i][1]:
            out = [short[i][2], long[i][2]]
            break
    
    return out


def express(transformed, interpol = False): # reconstructs the intensity levels of a lineage based of transform() as [frame, ID, short, long, histone]; interpol is the interpolation boolean
    out = []
    
    for t in transformed:
        temp = find(t[0], t[1])
        temph1 = -1
        temph2 = [-1, -1, -1]
        for i in range(len(histone)):
            if t[0] == histone[i][0] and t[1] == histone[i][1]:
                temph1 = histone[i][2]
                temph2 = histone[i][3]
                break
        if temp[0] > 0:
            out.append([t[0], t[1], temp[0], temp[1], temph1, temph2])
    
    return np.array(out)


def transform(string): # transforms the lineage() string output to a dataframe
    out = []

    while len(string) > 0:
        sub = string[0:string.index('<')]
        out.append([int(sub[0:sub.find('_')]), int(sub[(1 + sub.find('_')):len(sub)])])
        string = string.replace(sub + '<', '')
    
    return np.flip(np.array(out), axis=0)

""" load stack """

stack_category = "Nanog_Gata6"
stack_name = "230212_stack6"
metric_name = "MeanIntensity_nowarp"

main_directory = os.getcwd()
stack_directory = os.path.join(main_directory, "datasets", stack_category, stack_name)
extraction_directory = os.path.join(stack_directory, "extraction")

# load embryo tree graph

div_check = []
div = [] # contains all division frames
graph = [f for f in os.listdir(extraction_directory) if ('graph' in f) and (f.endswith('.csv'))] # finds the graph csv
track = [] # contains [[mother], [daughter]] cell [frame, ID] data

with open(os.path.join(extraction_directory, graph[0]), newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        s1 = row[0].replace("'", '')
        s1_1 = int(s1[0:s1.find('_')])
        s1_2 = int(s1[s1.find('_') + 1 : len(s1)])
        
        s2 = row[1].replace("'", '')
        s2_1 = int(s2[0:s2.find('_')])
        s2_2 = int(s2[s2.find('_') + 1 : len(s2)])
        
        tempM = [s1_1, s1_2]
        tempD = [s2_1, s2_2]
        
        if s2_1 < s1_1:
            temp = tempM
            tempM = tempD
            tempD = temp
            
        track.append([tempM, tempD])
        
        if tempM in div_check:
            div.append(tempM)
        else:
            div_check.append(tempM)

# extract intensity measurements

histone = []
with open(os.path.join(extraction_directory, 'extract_histone.csv'), newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    header = True
    for row in reader:
        if header:
            frame = row.index('Frame')
            ID = row.index('ID')
            intense = row.index(metric_name)
            centX = row.index("Centroid_1")
            centY = row.index("Centroid_2")
            centZ = row.index("Centroid_3")
            header = False
        else:
            histone.append([int(row[frame]), int(row[ID]), float(row[intense]), [float(row[centX]), float(row[centY]), float(row[centZ])]])
histone = np.array(histone)

short = []
with open(os.path.join(extraction_directory, 'extract_short.csv'), newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    header = True
    for row in reader:
        if header:
            frame = row.index('Frame')
            ID = row.index('ID')
            intense = row.index(metric_name)
            header = False
        else:
            short.append([int(row[frame]), int(row[ID]), float(row[intense])])
short = np.array(short)

long = []
with open(os.path.join(extraction_directory, 'extract_long.csv'), newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    header = True
    for row in reader:
        if header:
            frame = row.index('Frame')
            ID = row.index('ID')
            intense = row.index(metric_name)
            header = False
        else:
            long.append([int(row[frame]), int(row[ID]), float(row[intense])])
long = np.array(long)

# set tracking limits to frames with reporter data
start_frame = max(min(histone[:,0].astype(int)), min(short[:,0].astype(int)), min(long[:,0].astype(int)))
end_frame = min(max(histone[:,0].astype(int)), max(short[:,0].astype(int)), max(long[:,0].astype(int)))

allcells = frame_IDs(track, end_frame) # get the list of cells at the last frame

# compute and save reconstructed terminal lineage raw intensities as individual timeseries
stack_data = []
tick = 1
ticks = len(allcells)
for cell in allcells:
    print('Processing ' + stack_name + ' lineage:', tick, '/', ticks)
    lin = lineage(track, cell[0], cell[1])
    stack_data.append([cell, express(transform(lin[0])), lin[1]])
    tick += 1

# compute stack intensity max values
black_max = -1
blue_max = -1
red_max = -1
for i in range(len(stack_data)):
    blue_max = max(blue_max, max(stack_data[i][1][:,2]))
    red_max = max(red_max, max(stack_data[i][1][:,3]))
    black_max = max(black_max, max(stack_data[i][1][:,4]))

# compute stack intensity min values
black_min = 999_999_999
blue_min = 999_999_999
red_min = 999_999_999
for i in range(len(stack_data)):
    blue_min = min(blue_min, min(stack_data[i][1][:,2]))
    red_min = min(red_min, min(stack_data[i][1][:,3]))
    black_min = min(black_min, min(stack_data[i][1][:,4]))
    
# output
all_ranges = [[black_min, black_max], [blue_min, blue_max], [red_min, red_max]]
out_meta = [div, [start_frame, end_frame], all_ranges]
# stack_data