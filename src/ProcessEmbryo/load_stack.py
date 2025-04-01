""" import libraries """

import csv
import os
import numpy as np

""" supporting functions """

def frame_IDs(given_track, frame):
    """
    Returns a list of cells in a given track at the given frame:
        given_track: target track
        frame: target frame
    """
    out = []
    
    for t in given_track:
        if t[1][0] == frame:
            out.append(t[1])
    
    return out


def lineage_reconstruct(given_track, div, frame, ID):
    """
    Recursively reconstructs a cell's [frame, ID] lineage as a string spaced with '<':
        given_track: target track
        div: division history
        frame: frame
        ID: id
    """
    
    found = False
    mother = ['', []]
    i = 0
    
    if frame == given_track[0][0][0]:
        source_found = False
        j = 0
        
        while not source_found and given_track[j][0][0] == frame:
            if (given_track[j][0][0] == frame) and (given_track[j][0][1] == ID):
                mother = [str(given_track[j][0][0]) + '_' + str(given_track[j][0][1]) + '<', []]
                source_found = True
            j += 1
    
    else:
        while not found and i < len(given_track):
            if (given_track[i][1][0] == frame) and (given_track[i][1][1] == ID):
                lin = lineage_reconstruct(given_track, div, given_track[i][0][0], given_track[i][0][1])
                mother[0] = str(given_track[i][1][0]) + '_' + str(given_track[i][1][1]) + '<' + lin[0]
                mother[1] = lin[1]
                if given_track[i][1] in div:
                    mother[1].append(given_track[i][1])
                found = True
            else:
                i = i + 1
      
    return mother


def find_cells(frame, idi, source_stack):
    """
    Finds the short and long camera intensities for a given cell [frame, idi]:
        frame: target frame
        idi: target ID
        source_stack: raw stack data
    
    Returns [-1;-1] if the cell isn't found
    """
    out = [-1, -1]
    
    for i in range(len(source_stack[2])):
        if frame == source_stack[2][i][0] and idi == source_stack[2][i][1]:
            out = [source_stack[2][i][2], source_stack[1][i][2]]
            break
    
    return out


def lineage_transform(string):
    """
    Transforms a lineage string to an array of [frame, ID]:
        string: lineage_reconstruct() output string of cells in a lineage
    """
    out = []

    while len(string) > 0:
        sub = string[0:string.index('<')]
        out.append([int(sub[0:sub.find('_')]), int(sub[(1 + sub.find('_')):len(sub)])])
        string = string.replace(sub + '<', '')
    
    return np.flip(np.array(out), axis=0)


def lineage_express(transformed, source_stack):
    """
    Recostructs the lineage_transform() output as a list of [frame, ID, short, long, histone]:
        transformed: lineage_transform() list containing the sequence of tracked cells' [frame, ID]
        source_stack: raw stack data
    """
    out = []
    
    for t in transformed:
        temp = find_cells(t[0], t[1], source_stack)
        temph = -1
        temp_centroid = [-1, -1, -1]
        for i in range(len(source_stack[0])):
            if t[0] == source_stack[0][i][0] and t[1] == source_stack[0][i][1]:
                temph = source_stack[0][i][2]
                temp_centroid[0] = source_stack[0][i][3]
                temp_centroid[1] = source_stack[0][i][4]
                temp_centroid[2] = source_stack[0][i][5]
                break
        if temp[0] > 0:
            out.append([t[0], t[1], temp[0], temp[1], temph, temp_centroid[0], temp_centroid[1], temp_centroid[2]])
    
    return np.array(out)


def prune_stack(source_stack, allow = 0):
    """
    Returns a stack with shorter lineages removed:
        source_stack: raw stack data
        allow: tolerance on shorter lineages (frames)
    """
    out = []
    max_len = 0
    
    for cell in source_stack:
        if len(cell[1]) > max_len:
            max_len = len(cell[1])
    
    for cell in source_stack:
        if max_len - len(cell[1]) <= allow:
            out.append(cell)
    
    return out


""" main function """

def process_stack(stack_name, data_directory = "None", stack_category = "None", metric_name = "None", prune = [True, 0]):
    """
    Processes an extracted embryo stack:
        stack_name: name of the target stack
        data_directory: custom data folder name
        stack_category: custom category folder name
        metric_name: custom CSV column name
        prune: [bool, int] defining pruning and prune allowance
        
    Returns:
        stack_data: formatted stack lineages
        raw_data: raw stack data tracks and inesities as [track, raw_stack]
        meta: stack meta-data containing division frames, timing and min-max value ranges as [divisions, [start_frame, end_frame], all_ranges]
    """
    
    main_directory = os.getcwd()
    
    if data_directory == "None":
        data_directory = "datasets"
    
    if stack_category == "None":
        stack_category = "Nanog_Gata6"
    
    if metric_name == "None":
        metric_name = "MeanIntensity_nowarp"
    
    stack_directory = os.path.join(main_directory, data_directory, stack_category, stack_name)
    extraction_directory = os.path.join(stack_directory, "extraction")
    
    # load embryo tree graph
    
    div_check = []
    divisions = [] # contains all division frames
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
                divisions.append(tempM)
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
                histone.append([int(row[frame]), int(row[ID]), float(row[intense]), float(row[centX]), float(row[centY]), float(row[centZ])])
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
    
    # save raw data
    raw_stack = [histone, long, short]
    
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
        lin = lineage_reconstruct(track, divisions, cell[0], cell[1])
        stack_data.append([cell, lineage_express(lineage_transform(lin[0]), raw_stack), lin[1]])
        tick += 1
    
    # prune short lineages if enabled
    if prune[0]:
        stack_data = prune_stack(stack_data, prune[1])
    
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
        
    # format output
    raw_data = [track, raw_stack]
    all_ranges = [[black_min, black_max], [blue_min, blue_max], [red_min, red_max]]
    meta = [divisions, [start_frame, end_frame], all_ranges]
    
    return stack_data, raw_data, meta

