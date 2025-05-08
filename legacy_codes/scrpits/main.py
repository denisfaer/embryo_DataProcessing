""" import libraries """

import csv
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import pickle
    
""" configs """

use_intensity = 'MeanIntensity'
use_method = 'nowarp'
prune_termini = False # remove terminal cells that don't trace to the starting frame
trouble = True # print a message before each step [troubleshhoting boolean]
force_reprocessing = False # force data reprocessing even if pickled data is available
scatter_alpha = 0.5 # transparancy level of scatter plot markers
cbins = 50 # force a custom number of huistogram bins
atan_thr1 = 0.65 # arctan(short/long) threshold 1
atan_thr2 = 1.1 # arctan(short/long) threshold 2
cut = 1/2 # show only this last part of the video
plot = True # make plots
pheno_t = True # plot data separated into 
labels_max = 4 # allow title labels on a grid of no more than X by X
delin_plots = True # make deliniation plots

save_panels = False # save panels


""" functions """


def exists(frame, idi): # returns a boolean of if a mother cell with this frame-ID exists
    booly = False
    
    for i in range(len(track)):
        booly = booly or ((track[i][0][0] == frame) and (track[i][0][1] == idi))
    
    return booly


def lineage(frame, idi): # returns a string of frame_IDs for a given cell using the edges from the lineage graph
    found = False
    mother = ''
    i = 0
    
    while not found and i < len(track):
        if (track[i][1][0] == frame) and (track[i][1][1] == idi):
            mother = str(track[i][1][0]) + '_' + str(track[i][1][1]) + '<' + lineage(track[i][0][0], track[i][0][1])
            found = True
        else:
            i = i + 1
      
    return mother


def transform(string): # transforms the lineage() string output to a dataframe
    out = []

    while len(string) > 0:
        sub = string[0:string.index('<')]
        out.append([int(sub[0:sub.find('_')]), int(sub[(1+sub.find('_')):len(sub)])])
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
        if temp[0] > 0:
            out.append([t[0], t[1], temp[0], temp[1], temph])
    
    return np.array(out)


def terminus(track, frame): # returns a list of cells at target frame
    out = []
    
    for t in track:
        if t[1][0] == frame:
            out.append(t[1])
    
    return out


def prune(stack, allow):  # removes cells that don't trace to the starting frame within allowed deviation
    out = []
    
    for i in range(len(stack)):
        if stack[i][1][0, 0] - start_frame <= allow:
            out.append(stack[i])
    
    return out


def minimal_rectangle(L): # calculates the dimensions of a minimal rectangular grid that can contain L elements
    N = math.ceil(math.sqrt(L))
    M = (L // N) + 1
    
    if N**2 == L:
        M = M - 1
    
    return N, M


def deliniation_plot(label): # makes a large text panel deliniating what the following panels show
    x = 0
    y = 0
    plt.figure(figsize = (100, 100))
    plt.scatter(x, y, marker = 'x', alpha = 0)
    plt.annotate(label, xy=(x, y), ha='center', fontsize = 1000)
    plt.axis('off')
    plt.show()


def make_plots(lineages): # make a set of plots for a set of tracked cells
    
    # format key variables
    stack = lineages
    cells = stack[1]
    N, M = minimal_rectangle(cells)
    titles = (N <= labels_max)
    
    deliniation_plot(stack[0])
    
    """ globally minmax-normalized plots """
    
    deliniation_plot(stack[0] + ' globally minmax-normalized plots')
    
    # plot globally minmax-normalized intensities
    for i in range(cells):
        plt.subplot(N, M, i + 1)
        plt.axis('off')
        if titles:
            plt.title(stack[2][i][0], fontsize = 10)
        plt.ylim(0, 1)
        plt.xlim(start_frame, end_frame)
        plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,4] - black_min)/(black_max - black_min), c='Black')
        plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min), c='Blue')
        plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,3] - red_min)/(red_max - red_min), c='Red')
    plt.show()

    # plot globally minmax-normalized reporters
    for i in range(cells):
        plt.subplot(N, M, i + 1)
        plt.axis('off')
        plt.ylim(0, 1)
        plt.xlim(start_frame, end_frame)
        plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min), c='Blue')
        plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,3] - red_min)/(red_max - red_min), c='Red')
    plt.show()

    # plot the terminal frame globally minmax-normallized short vs long reporters
    for i in range(cells):
        plt.scatter((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min), (stack[2][i][1][len(stack[2][i][1])-1, 3] - red_min)/(red_max - red_min), s=50*(stack[2][i][1][len(stack[2][i][1])-1, 4] - black_min)/(black_max - black_min), c='Black', alpha=scatter_alpha)
    plt.title('Globally minmax-normallized short vs long intensities')
    plt.show()
    
    # plot the terminal frame globally minmax-normallized short vs long reporters with strict axes
    for i in range(cells):
        plt.scatter((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min), (stack[2][i][1][len(stack[2][i][1])-1, 3] - red_min)/(red_max - red_min), s=50*(stack[2][i][1][len(stack[2][i][1])-1, 4] - black_min)/(black_max - black_min), c='Black', alpha=scatter_alpha)
    plt.title('Globally minmax-normallized short vs long intensities')
    plt.ylim(0, 1)
    plt.xlim(0, 1)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.show()
    
    """ raw data angle plots (from x = short, y = long, calculated as arctan) """
    
    deliniation_plot(stack[0] + ' raw angle plots')
    
    # plot the short-long angle of reporter intensities
    for i in range(cells):
        plt.subplot(N, M, i + 1)
        plt.ylim(0, math.pi/2)
        plt.axis('off')
        if titles:
            plt.title(stack[2][i][0], fontsize = 10)
        temp = []
        for j in range(len((stack[2][i][1]))):
            temp.append(math.pi/2 - math.atan(stack[2][i][1][j, 2] / stack[2][i][1][j, 3]))
        plt.plot(stack[2][i][1][:,0], temp, c='Orange')
    plt.show()
    
    # plot the short-long angle of globally minmax-normallized intensity ratio histogram
    temp = []
    for i in range(cells):
        temp.append(math.pi/2 - math.atan(stack[2][i][1][len(stack[2][i][1]) - 1, 2] / stack[2][i][1][len(stack[2][i][1]) - 1, 3]))
    plt.hist(temp, bins = cbins)
    plt.title('Globally minmax-normallized short/long ratio distribution (arctan)')
    plt.show()
    
    # plot the short-long angle of reporter intensities with threshold lines
    for i in range(cells):
        plt.subplot(N, M, i + 1)
        plt.ylim(0, math.pi/2)
        plt.axis('off')
        temp = []
        for j in range(len((stack[2][i][1]))):
            temp.append(math.pi/2 - math.atan(stack[2][i][1][j, 2] / stack[2][i][1][j, 3]))
        plt.plot(stack[2][i][1][:,0], temp, c='Orange')
        plt.axhline(y = atan_thr1, color = 'Black', linestyle = '-', linewidth = 0.5) 
        plt.axhline(y = atan_thr2, color = 'Black', linestyle = '-', linewidth = 0.5) 
    plt.show()
        
    # plot the short-long angle of globally minmax-normallized intensity ratios with threshold lines (cut marker-dash version)
    for i in range(cells):
        plt.subplot(N, M, i + 1)
        plt.ylim(0, math.pi/2)
        plt.axis('off')
        temp = []
        for j in range(round(len((stack[2][i][1])) * (1 - cut)), len((stack[2][i][1]))):
            temp.append(math.pi/2 - math.atan(stack[2][i][1][j, 2] / stack[2][i][1][j, 3]))
        plt.plot(stack[2][i][1][round(len((stack[2][i][1])) * (1 - cut)):,0], temp, marker='o', markersize=1, linestyle='-', linewidth = 0.2, c='Orange')
        plt.axhline(y = atan_thr1, color = 'Black', linestyle = '-', linewidth = 0.5) 
        plt.axhline(y = atan_thr2, color = 'Black', linestyle = '-', linewidth = 0.5) 
    plt.show()
    
    """ globally minmax-normallized angle plots (from x = short, y = long, calculated as arctan) """
    
    deliniation_plot(str(stack[0] + ' globally minmax-normallized angle plots'))
    
    # plot the short-long angle of globally minmax-normallized intensity ratios
    for i in range(cells):
        plt.subplot(N, M, i + 1)
        plt.ylim(0, math.pi/2)
        plt.axis('off')
        if titles:
            plt.title(stack[2][i][0], fontsize = 10)
        temp = []
        for j in range(len((stack[2][i][1]))):
            temp.append(math.pi/2 - math.atan(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
        plt.plot(stack[2][i][1][:,0], temp, c='Orange')
    plt.show()
    
    # plot the short-long angle of globally minmax-normallized intensity ratio histogram
    temp = []
    for i in range(cells):
        temp.append(math.pi/2 - math.atan(((stack[2][i][1][len(stack[2][i][1]) - 1, 2] - blue_min)/(blue_max - blue_min)) / ((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min))))
    plt.hist(temp, bins = cbins)
    plt.title('Globally minmax-normallized short/long ratio distribution (arctan)')
    plt.show()
    
    # plot the short-long angle of globally minmax-normallized intensity ratios with threshold lines
    for i in range(cells):
        plt.subplot(N, M, i + 1)
        plt.ylim(0, math.pi/2)
        plt.axis('off')
        temp = []
        for j in range(len((stack[2][i][1]))):
            temp.append(math.pi/2 - math.atan(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
        plt.plot(stack[2][i][1][:,0], temp, c='Orange')
        plt.axhline(y = atan_thr1, color = 'Black', linestyle = '-', linewidth = 0.5) 
        plt.axhline(y = atan_thr2, color = 'Black', linestyle = '-', linewidth = 0.5) 
    plt.show()
        
    # plot the short-long angle of globally minmax-normallized intensity ratios with threshold lines (cut marker-dash version)
    for i in range(cells):
        plt.subplot(N, M, i + 1)
        plt.ylim(0, math.pi/2)
        plt.axis('off')
        temp = []
        for j in range(round(len((stack[2][i][1])) * (1 - cut)), len((stack[2][i][1]))):
            temp.append(math.pi/2 - math.atan(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
        plt.plot(stack[2][i][1][round(len((stack[2][i][1])) * (1 - cut)):,0], temp, marker='o', markersize=1, linestyle='-', linewidth = 0.2, c='Orange')
        plt.axhline(y = atan_thr1, color = 'Black', linestyle = '-', linewidth = 0.5) 
        plt.axhline(y = atan_thr2, color = 'Black', linestyle = '-', linewidth = 0.5) 
    plt.show()


def colsel(p): # select colour based on called phenotype
    if p == 'EPI':
        col = 'Red'
    elif p == 'PE':
        col = 'Blue'
    elif p == 'ICM':
        col = 'Purple'
    else:
        col = 'Black'
        print('Color selection failed, defaulting to: ' + col)
    
    return col

def scatter_size(val): # scatter plot marker size
    out = 50 * (stack[2][i][1][len(stack[2][i][1])-1, 4] - black_min)/(black_max - black_min)
    return out

""" start analysis """

main_directory = os.getcwd()

# find all stack data folders: must be named stackXYZ with XYZ being the unique name (not character limited) and NOT contain a file type (dot symbol)
given_stacks = [f for f in os.listdir(main_directory) if ('stack' in f) and (not '.' in f)]
pickles = [f for f in os.listdir(main_directory) if '.pkl' in f]

all_data = [['stack_name', 'short_intensity', 'long_intensity', 'histone_intensity']] # master list containing raw import data for all stacks
all_tracks = [['stack_name', 'num_cells', 'data']] # master list containing all tracks' raw intensities
all_phenotypes = [['stack_name', 'data']] # master list containing all phenotype-sorted tracks' raw intensities


""" process individual stacks """

### if len(pickles) < 3 or force_reprocessing:

for given_stack in given_stacks:
    stack_directory = os.path.join(main_directory, given_stack)
    
    # grab the lineage data from the graph - for now I manually copied edges from the .mat to a .csv file since I couldn't make .mat reading to work
    graph = [f for f in os.listdir(stack_directory) if ('graph' in f)and(f.endswith('.csv'))]
    track = [] # contains [mother, daughter] cell [frame, ID] data
    with open(os.path.join(stack_directory, graph[0]), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            s1 = row[0].replace("'", '')
            s1_1 = int(s1[0:s1.find('_')])
            s1_2 = int(s1[s1.find('_') + 1 : len(s1)])
            
            s2 = row[1].replace("'", '')
            s2_1 = int(s2[0:s2.find('_')])
            s2_2 = int(s2[s2.find('_') + 1 : len(s2)])
            
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
    
    """ format and globally min-max normallize stack data """
    
    # set tracking limits to frames with reporter data
    start_frame = max(min(histone[:,0].astype(int)), min(short[:,0].astype(int)), min(long[:,0].astype(int)))
    end_frame = min(max(histone[:,0].astype(int)), max(short[:,0].astype(int)), max(long[:,0].astype(int)))
    
    # get the list of cells at the last frame
    allcells = terminus(track, end_frame)
    
    # compute and save reconstructed terminal lineage raw intensities as individual timeseries
    stack_data = []
    for cell in allcells:
        if trouble:
            print(lineage(cell[0], cell[1]))
            print()
        stack_data.append([cell, express(transform(lineage(cell[0], cell[1])))])
    
    if trouble:
        temp = []
        for i in range(len(stack_data)):
            temp.append(len(stack_data[i][1]))
        
        plt.hist(temp)
        plt.title('Track length in ' + given_stack)
        plt.show()
    
    if prune_termini:
        stack_data = prune(stack_data)
    
    # save data
    all_data.append([given_stack, short, long, histone])
    all_tracks.append([given_stack, len(stack_data), stack_data])
    
    stack_data = [given_stack, len(stack_data), stack_data]
    
    # compute stack intensity max values
    black_max = -1
    blue_max = -1
    red_max = -1
    for i in range(stack_data[1]):
        blue_max = max(blue_max, max(stack_data[2][i][1][:,2]))
        red_max = max(red_max, max(stack_data[2][i][1][:,3]))
        black_max = max(black_max, max(stack_data[2][i][1][:,4]))
        
    # compute stack intensity min values
    black_min = 999_999_999
    blue_min = 999_999_999
    red_min = 999_999_999
    for i in range(stack_data[1]):
        blue_min = min(blue_min, min(stack_data[2][i][1][:,2]))
        red_min = min(red_min, min(stack_data[2][i][1][:,3]))
        black_min = min(black_min, min(stack_data[2][i][1][:,4]))
    
    """ get and format called phenotypes """
    
    if pheno_t:
        
        # get the list of cells with called phenotypes
        cells_called = []
        with open(os.path.join(stack_directory, 'phenotypes.csv'), newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for row in reader:
                s1 = row[0].replace("'", '')
                s1_1 = int(s1[0:s1.find('_')])
                s1_2 = int(s1[s1.find('_') + 1 : len(s1)])
                cells_called.append([s1_1, s1_2, row[1]])
        cells_called = np.array(cells_called)
        
        # make a list of called phenotypes
        pheno = []
        for i in range(len(cells_called)):
            if cells_called[i][2] not in pheno:
                pheno.append(cells_called[i][2])
        
        # sort called cells by phenotype
        pheno_data = []
        
        for p in pheno:
            temp_cells = []
        
            for i in range(len(cells_called)):
                if cells_called[i][2] == p:
                    temp_cells.append([int(cells_called[i][0]), int(cells_called[i][1])])
            
            temp_data = []
            for cell in temp_cells:
                temp_data.append([cell, express(transform(lineage(cell[0], cell[1])))])
                
            pheno_data.append([p, [given_stack + '-' + p, len(temp_data), temp_data]])
        
        all_phenotypes.append([given_stack, pheno_data])
    
    """ plot stack data """
    if plot:
        
        # make all-cell plots
        make_plots(stack_data)
        
        # make phenotype-specific plots
        if pheno_t:
            
            # make plots for each phenotype
            for p in range(len(pheno)):
                make_plots(pheno_data[p][1])
            
            # plot the terminal frame globally minmax-normallized short vs long reporters
            for p in range(len(pheno)):
                stack = pheno_data[p][1]
                for i in range(stack[1]):
                    plt.scatter((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min), (stack[2][i][1][len(stack[2][i][1])-1, 3] - red_min)/(red_max - red_min), s=50*(stack[2][i][1][len(stack[2][i][1])-1, 4] - black_min)/(black_max - black_min), c=colsel(pheno_data[p][0]), alpha=scatter_alpha)
                
            plt.title(given_stack + ' short vs long intensities of called cells')
            plt.xlabel('Minmax-normalized Gata6')
            plt.ylabel('Minmax-normalized Nanog')
            plt.show()
     
""" pickle data """

pickle_file = open(os.path.join(main_directory, 'raw_data.pkl'), 'wb')
pickle.dump(all_data, pickle_file)
pickle_file.close()

pickle_file = open(os.path.join(main_directory, 'track_data.pkl'), 'wb')
pickle.dump(all_tracks, pickle_file)
pickle_file.close()

pickle_file = open(os.path.join(main_directory, 'phenotype_data.pkl'), 'wb')
pickle.dump(all_phenotypes, pickle_file)
pickle_file.close()

""" add plotting-only option """

if plot and False:
    N, M = minimal_rectangle(len(all_phenotypes))
    
    for i in range(len(all_phenotypes)):
        plt.subplot(N, M, i + 1)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        
        for p in range(len(pheno)):
            stack = all_phenotypes[i + 1][1][p][1]
            for i in range(stack[1]):
                plt.scatter((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min), (stack[2][i][1][len(stack[2][i][1])-1, 3] - red_min)/(red_max - red_min), s=50*(stack[2][i][1][len(stack[2][i][1])-1, 4] - black_min)/(black_max - black_min), c=colsel(pheno_data[p][0]), alpha=scatter_alpha)
            
        plt.title(given_stack)
    
    plt.show()
