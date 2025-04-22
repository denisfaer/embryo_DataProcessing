""" import libraries """

import csv
import os
import numpy as np
from pathlib import Path

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


def lineage_reconstruct(given_track, div, frame, ID, echo = True, full_echo = False):
    """
    Recursively reconstructs a cell's [frame, ID] lineage as a string spaced with '<':
        given_track: target track
        div: division history
        frame: tagert frame
        ID: target id
        echo: report found root cells that don't trace to track start
        full_echo: report all root cells found
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
    
    if not found:
        if full_echo:
            print('REPORT: Lineage root cell identified as [' + str(frame) + '_' + str(ID) + ']')
        elif echo and frame != given_track[0][0][0]:
            print('WARNING: Lineage root cell identified as [' + str(frame) + '_' + str(ID) + ']')
    
    return mother


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


def find_cells(frame, ID, ext_stack):
    """
    Finds intensities for a given cell [frame, idi] in all channels:
        frame: target frame
        ID: target ID
        ext_stack: extracted stack data
    
    Returns [NaN] if the cell isn't found
    """
    out = []
    
    for chan in ext_stack:
        for i in range(len(chan)):
            if frame == chan[i][0] and ID == chan[i][1]:
                out.append(chan[i][2])
                break
    
    if len(out) == 0:
        out = [np.nan]
        print('ERROR: find_cell([' + str(frame) + '_' + str(ID) + ']) - cell not found')
    
    return out


def lineage_express(transformed, ext_stack):
    """
    Recostructs the lineage_transform() output as a list of [frame, ID, histone, centroid1, centroid2, centroid3, channel1, ...]:
        transformed: lineage_transform() list containing the sequence of tracked cells' [frame, ID]
        ext_stack: extracted stack data
    """
    out = []
    
    for t in transformed:
        temp = find_cells(t[0], t[1], ext_stack)
        temp_centroid = [np.nan, np.nan, np.nan]
        for i in range(len(ext_stack[0])):
            if t[0] == ext_stack[0][i][0] and t[1] == ext_stack[0][i][1]:
                temp_centroid[0] = ext_stack[0][i][3]
                temp_centroid[1] = ext_stack[0][i][4]
                temp_centroid[2] = ext_stack[0][i][5]
                break
        if temp[0] > 0:
            out_tmp = [t[0], t[1]]
            out_tmp.append(temp_centroid[0])
            out_tmp.append(temp_centroid[1])
            out_tmp.append(temp_centroid[2])
            for chan in temp:
                out_tmp.append(chan)
            out.append(out_tmp)
    
    return out


def prune_stack(recon_stack, allow = 0, echo = True):
    """
    Returns a stack with shorter lineages removed:
        recon_stack: reconstructed stack data
        allow: tolerance for shorter lineages (in frames, from longest lineage)
        echo: report lineage pruning
    """
    out = []
    max_len = 0
    
    for cell in recon_stack:
        if len(cell[1]) > max_len:
            max_len = len(cell[1])
    
    if echo:
        print()
        print('Prune incomplete lineages:')
    
    for i in range(len(recon_stack)):
        if max_len - len(recon_stack[i][1]) <= allow:
            out.append(recon_stack[i])
        elif echo:
            print('REPORT: lineage ' + str(i + 1) + ' pruned; terminal cell [' + str(recon_stack[i][0][0]) + '_' + str(recon_stack[i][0][1]) + ']')
    
    for cell in recon_stack:
        if max_len - len(cell[1]) <= allow:
            out.append(cell)
    
    return out


def process_stack(stack_dir, channels = 3, histone_metric = "MeanIntensity_nowarp", channel_metric = "MeanIntensity_rigid", histone_centroid = True, prune = [True, 0], glob_echo = True):
    """
    Processes an extracted embryo stack:
        stack_dir: string path to to directory with the target stack data
        channels: number of channels to expect; min = 1
        histone_metric: CSV column name to use for histone intensity measurement
        channel_metric: CSV column name to use for other channels' intensity measurement
        histone_centroid: extract nuclear centroid coordinates
        prune: [bool, int] defining pruning and prune allowance
        glob_echo: report progress at each key step
        
    Returns [stack_data, raw_data, meta]:
        stack_data: formatted stack lineages
        raw_data: raw extracted stack data tracks as [track, extracted_stack]
        meta: stack meta-data containing [divisions, [start_frame, end_frame], [stack_mins, stack_maxs]]
    """
    
    # grab stack name and stack directory
    last_dash = stack_dir.rfind('/')
    if last_dash == len(stack_dir) - 1:
        stack_dir = stack_dir[:-1]
        last_dash = stack_dir.rfind('/')
    stack_name = stack_dir[stack_dir.rfind('/') + 1 : len(stack_dir)]
    
    if glob_echo:
        print()
        print('Processing ' + stack_name)
    
    stack_directory = Path(stack_dir)
    extraction_directory = os.path.join(stack_directory, "extraction")
    
    """ load embryo tree graph """
    if glob_echo:
        print()
        print('Loading the lineage tree graph:')
    
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
    
    if glob_echo:
        print('Lineage tree graph loaded successfully')
    
    """ extract intensity measurements """
    if glob_echo:
        print()
        print('Loading channel extractions:')
    
    extracted_stack = []
    
    # extract histone data
    histone = []
    with open(os.path.join(extraction_directory, 'extract_histone.csv'), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        header = True
        for row in reader:
            if header:
                frame = row.index('Frame')
                ID = row.index('ID')
                intense = row.index(histone_metric)
                if histone_centroid:
                    centX = row.index("Centroid_1")
                    centY = row.index("Centroid_2")
                    centZ = row.index("Centroid_3")
                header = False
            else:
                centroid = [np.nan, np.nan, np.nan] # default non-specified centroid
                if histone_centroid:
                    centroid = [float(row[centX]), float(row[centY]), float(row[centZ])]
                histone.append([int(row[frame]), int(row[ID]), float(row[intense]), centroid[0], centroid[1], centroid[2]])
    histone = np.array(histone)
    extracted_stack.append(histone)
    if glob_echo:
        print('Histone loaded')
    
    # extract other channels
    for i in range(1, channels):
        tmp = []
        with open(os.path.join(extraction_directory, 'extract_signal_' + str(i) + '.csv'), newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            header = True
            for row in reader:
                if header:
                    frame = row.index('Frame')
                    ID = row.index('ID')
                    intense = row.index(channel_metric)
                    header = False
                else:
                    tmp.append([int(row[frame]), int(row[ID]), float(row[intense])])
        tmp = np.array(tmp)
        extracted_stack.append(tmp)
        if glob_echo:
            print('Channel_' + str(i) + ' loaded')
    
    """ reconstruct cell lineages """
    if glob_echo:
        print()
        print('Reconstructing cell lineages:')
    
    # set tracking limits to frames with reporter data
    start_frame = -1
    for chan in extracted_stack:
        start_frame = max(start_frame, min(chan[:,0].astype(int)))
    
    end_frame = 999_999_999
    for chan in extracted_stack:
        end_frame = min(end_frame, max(chan[:,0].astype(int)))
    
    # get the list of cells at the last frame
    allcells = frame_IDs(track, end_frame)
    
    # compute and save reconstructed terminal lineage raw intensities as individual timeseries
    stack_data = []
    tick = 1
    ticks = len(allcells)
    for cell in allcells:
        print('Processing ' + stack_name + ' lineage', tick, '/', ticks)
        lin = lineage_reconstruct(track, divisions, cell[0], cell[1], echo = glob_echo)
        stack_data.append([cell, lineage_express(lineage_transform(lin[0]), extracted_stack), lin[1]])
        tick += 1
    
    # prune short lineages (if enabled)
    if prune[0]:
        stack_data = prune_stack(stack_data, prune[1], echo = glob_echo)
        
    """ process data for outputting """
    if glob_echo:
        print()
        print('Computing metadata:')
    
    # compute stack intensity max values
    stack_maxs = []
    for i in range(channels + 1):
        tmp_max = -1
        for lineage in stack_data:
            for cell in lineage[1]:
                if (len(cell) - 1) <= (i + 5):
                    tmp_max = max(tmp_max, max(cell[i + 5]))
        stack_maxs.append(tmp_max)
    
    if glob_echo:
        print('Maximum intensities:')
        print(stack_maxs)
    
    # compute stack intensity min values
    stack_mins = []
    for i in range(channels + 1):
        tmp_min = 999_999_999
        for lineage in stack_data:
            for cell in lineage[1]:
                if (len(cell) - 1) <= (i + 5):
                    tmp_min = min(tmp_min, min(cell[i + 5]))
        stack_mins.append(tmp_min)
    
    if glob_echo:
        print('Minimum intensities:')
        print(stack_mins)
        
    # format output
    stack_raw = [track, extracted_stack]
    stack_meta = [divisions, [start_frame, end_frame], [stack_mins, stack_maxs]]
    
    if glob_echo:
        print()
        print('Finished processing' + stack_name)
    
    return [stack_data, stack_raw, stack_meta]