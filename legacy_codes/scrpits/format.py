""" import libraries """

import pickle
import os
import math
import matplotlib.pyplot as plt

""" configs """

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

""" load extracted data """

pickle_file = open(os.path.join(os.getcwd(), 'track_data.pkl'), 'rb')
stacks = pickle.load(pickle_file)
pickle_file.close()

""" process individual stacks """

for stack in stacks:
    # get the number of cells at the terminus
    cells = stack[1]
    
    # calculate the dimensions of a minimal rectangular grid that can contain all plots
    N = math.ceil(math.sqrt(cells))
    M = (cells // N) + 1
    
    # plot all intensities' raw measurements
    if raw_plots:
        if print_plot_start: 
            print(f'{raw_plots = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,4], c='Black')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2], c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3], c='Red')
        plt.show()
    
    # plot last frame short-vs-long raw intensity scatter
    if raw_plot_scatter:
        if print_plot_start: 
            print(f'{raw_plot_scatter = }')
        for i in range(cells):
            plt.scatter(stack[2][i][1][len(stack[2][i][1])-1, 2], stack[2][i][1][len(stack[2][i][1])-1, 3], s=(stack[2][i][1][len(stack[2][i][1])-1, 4]**marker_size_scaling)/(10**(marker_size_scaling-1)), c='Black')
        plt.show()
    
    # plot histone-normalized [within each frame] raw reporter intensities
    if raw_hist_norm:
        if print_plot_start: 
            print(f'{raw_hist_norm = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/stack[2][i][1][:,4], c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3]/stack[2][i][1][:,4], c='Red')
        plt.show()
    
    # plot last frame short-vs-long raw intensity scatter with histone raw intensity determining size
    if raw_plot_scatter_hist:
        if print_plot_start: 
            print(f'{raw_plot_scatter_hist = }')
        for i in range(cells):
            plt.scatter(stack[2][i][1][len(stack[2][i][1])-1, 2], stack[2][i][1][len(stack[2][i][1])-1, 3], s=stack[2][i][1][len(stack[2][i][1])-1, 4]/10, c='Black')
        plt.show()
    
    # plot short/long intensity raw intensity ratios
    if raw_ratio:
        if print_plot_start: 
            print(f'{raw_ratio = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/stack[2][i][1][:,3], c='Orange')
        plt.show()
    
    # plot individually max-normalized intensities
    if ind_max_norm_hist:
        if print_plot_start: 
            print(f'{ind_max_norm_hist = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,4]/max(stack[2][i][1][:,4]), c='Black')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/max(stack[2][i][1][:,2]), c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3]/max(stack[2][i][1][:,3]), c='Red')
        plt.show()
    
    # plot individually max-normalized reporters
    if ind_max_norm:
        if print_plot_start: 
            print(f'{ind_max_norm = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/max(stack[2][i][1][:,2]), c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3]/max(stack[2][i][1][:,3]), c='Red')
        plt.show()
    
    # compute global max values if needed
    if max_norm or minmax_norm:
        black_max = -1
        blue_max = -1
        red_max = -1
        for i in range(cells):
            blue_max = max(blue_max, max(stack[2][i][1][:,2]))
            red_max = max(red_max, max(stack[2][i][1][:,3]))
            black_max = max(black_max, max(stack[2][i][1][:,4]))
    
    # plot globally max-normalized intensities
    if glo_max_norm_hist:
        if print_plot_start: 
            print(f'{glo_max_norm_hist = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,4]/black_max, c='Black')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/blue_max, c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3]/red_max, c='Red')
        plt.show()
    
    # plot globally max-normalized reporters
    if glo_max_norm:
        if print_plot_start: 
            print(f'{glo_max_norm = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/blue_max, c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3]/red_max, c='Red')
        plt.show()
    
    # plot the terminal frame globally max-normallized short vs long reporters
    if glo_max_scatter:
        if print_plot_start: 
            print(f'{glo_max_scatter = }')
        for i in range(cells):
            plt.scatter(stack[2][i][1][len(stack[2][i][1])-1, 2]/blue_max, stack[2][i][1][len(stack[2][i][1])-1, 3]/red_max, s=10*stack[2][i][1][len(stack[2][i][1])-1, 4]/black_max, c='Black')
        plt.show()
    
    # plot the strict (axis from 0 to 1) terminal frame globally max-normallized short vs long reporters
    if glo_max_scatter_strict:
        if print_plot_start: 
            print(f'{glo_max_scatter_strict = }')
        for i in range(cells):
            plt.scatter(stack[2][i][1][len(stack[2][i][1])-1, 2]/blue_max, stack[2][i][1][len(stack[2][i][1])-1, 3]/red_max, s=10*stack[2][i][1][len(stack[2][i][1])-1, 4]/black_max, c='Black')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.show()