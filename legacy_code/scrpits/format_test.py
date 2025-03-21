""" import libraries """

import pickle
import os
import math
import matplotlib.pyplot as plt
import statistics

""" configs """

# global boolean configs

all_true = False # create all plots regardless of individual settings
print_plot_start = True # print a message when a plot has started being made [troubleshooting boolean]
same_scale_timeseries = True # force timeseries with a variable scale (e.g., absolute or ratio) to set the same y limits (observed extremes)

# specific boolean configs

raw_plots = False # create raw intensity value timeseries plots?
raw_plot_scatter = False # create raw short vs raw long intensity scatter?
raw_hist_norm = False # create raw short and long (normalized by raw histone) intensity timeseries plots?
raw_plot_scatter_hist = False # create raw short vs raw long intensity scatter where raw histone intensity determines marker size?

raw_ratio = False # create raw short over raw long intensity timeseries plots?
raw_ratio_distrib = False # create a histogram of terminal cell short/long raw intensities?

ind_max_norm_hist = False # create individually max-normalized intensity timeseries plots (raw values divided by highest within that reporter in that lineage)?
ind_max_norm = False # create individually max-normalized reporter timeseries plots (raw values divided by highest within that reporter in that lineage)?

glo_max_norm_hist = False # create globally max-normalized intensity timeseries plots (raw values divided by highest within that reporter in that stack)?
glo_max_norm = False # create globally max-normalized reporter timeseries plots (raw values divided by highest within that reporter in that stack)?
glo_max_scatter = False # create globally max-normalized short vs long intensity scatter?
glo_max_scatter_strict = False # create globally max-normalized short vs long intensity scatter with axes locked at 0 to 1?

glo_minmax_norm_histone = True # create globally minmax-normalized intensity timeseries plots (raw values divided by highest within that reporter in that stack)?
glo_minmax_norm = True # create globally minmax-normalized reporter timeseries plots (raw values divided by highest within that reporter in that stack)?
glo_minmax_scatter = True # create globally minmax-normalized short vs long intensity scatter?

minmax_ratio = True # create globally minmax-normallized short vs long intensity timeseries plots?
glo_minmax_scatter_strict = True
minmax_ratio_rev = False
minmax_scatter = False # create globally minmax-normalized short vs long intensity scatter?
minmax_distrib = True # create a histogram of terminal cell short/long globally minmax-normallized intensities?
minmax_distrib_arctan = True
minmax_ratio_log = True
minmax_ratio_arctan = True
minmax_distrib_log = True # create a histogram of log10 terminal cell short/long globally minmax-normallized intensities?
minmax_distrib_log_cbin = True
minmax_distrib_arctan_cbin = True

glo_minmax_ratio_scatter = False
glo_minmax_arctan_scatter = True
glo_minmax_ratio_scatter_avghis = False
glo_minmax_ratio_scatter_avgshort = False
glo_minmax_ratio_scatter_avglong = False

# specific-dependant boolean computation

max_norm = ind_max_norm_hist or ind_max_norm or glo_max_norm_hist or glo_max_norm or glo_max_scatter or glo_max_scatter_strict
minmax_norm = glo_minmax_norm_histone or glo_minmax_norm or glo_minmax_scatter or minmax_ratio or minmax_scatter

# numerical constant configs

scatter_alpha = 0.5 # transparancy level of scatter plot markers
cbins = 50 # force a custom number of huistogram bins

atan_thr1 = 0.65 # arctan(short/long) threshold 1
atan_thr2 = 1.1 # arctan(short/long) threshold 2

ex = 2 # marker size polinomial exaggeration constant

# temporary
minmax_ratio_arctan_wthr = True #
minmax_ratio_arctan_wthr_cut = True
cut = 1/2 # show only this last part of the video

""" load extracted data """

pickle_file = open(os.path.join(os.getcwd(), 'track_data.pkl'), 'rb')
stacks = pickle.load(pickle_file)
pickle_file.close()

""" process individual stacks """

for stack in stacks:
    
    # make a delineation 'plot'
    plt.plot(0, 0)
    plt.axis('off')
    plt.title(stack[0])
    plt.show()
    
    # get the number of cells at the terminus
    cells = stack[1]
    
    # calculate the dimensions of a minimal rectangular grid that can contain all plots
    N = math.ceil(math.sqrt(cells))
    M = (cells // N) + 1
    
    # plot all intensities' raw measurements
    if same_scale_timeseries:
        ymin = 999_999_999
        ymax = -ymin
        for i in range(cells):
            ymin = min(ymin, min(stack[2][i][1][:,2]), min(stack[2][i][1][:,3]), min(stack[2][i][1][:,4]))
            ymax = max(ymax, max(stack[2][i][1][:,2]), max(stack[2][i][1][:,3]), max(stack[2][i][1][:,4]))
    
    if raw_plots or all_true:
        if print_plot_start: 
            print(f'{raw_plots = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            if same_scale_timeseries:
                plt.ylim(ymin, ymax)
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,4], c='Black')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2], c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3], c='Red')
        plt.show()
    
    # plot last frame short-vs-long raw intensity scatter
    if raw_plot_scatter or all_true:
        if print_plot_start: 
            print(f'{raw_plot_scatter = }')
        for i in range(cells):
            plt.scatter(stack[2][i][1][len(stack[2][i][1])-1, 2], stack[2][i][1][len(stack[2][i][1])-1, 3], c='Black', alpha=scatter_alpha)
        plt.title('Raw short vs long intensities')
        plt.show()
    
    # plot histone-normalized [within each frame] raw reporter intensities
    if raw_hist_norm or all_true:
        if print_plot_start: 
            print(f'{raw_hist_norm = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/stack[2][i][1][:,4], c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3]/stack[2][i][1][:,4], c='Red')
        plt.show()
    
    # plot last frame short-vs-long raw intensity scatter with histone raw intensity determining size
    if raw_plot_scatter_hist or all_true:
        if print_plot_start: 
            print(f'{raw_plot_scatter_hist = }')
        for i in range(cells):
            plt.scatter(stack[2][i][1][len(stack[2][i][1])-1, 2], stack[2][i][1][len(stack[2][i][1])-1, 3], s=(stack[2][i][1][len(stack[2][i][1])-1, 4] ** ex)/(15 ** ex), c='Black', alpha=scatter_alpha)
        plt.title('Raw short vs long intensities (markers sized by histone intensity)')
        plt.show()
    
    # plot short/long raw intensity ratios
    if raw_ratio or all_true:
        if print_plot_start: 
            print(f'{raw_ratio = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/stack[2][i][1][:,3], c='Orange')
        plt.show()
    
    # plot short/long raw intensity ratio histogram
    if raw_ratio_distrib or all_true:
        if print_plot_start: 
            print(f'{raw_ratio_distrib = }')
        temp = []
        for i in range(cells):
            temp.append((stack[2][i][1][len(stack[2][i][1])-1, 2]/stack[2][i][1][len(stack[2][i][1])-1, 3]))
        plt.hist(temp)
        plt.title('Raw short/long ratio distribution')
        plt.show()
    
    # plot individually max-normalized intensities
    if ind_max_norm_hist or all_true:
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
    if ind_max_norm or all_true:
        if print_plot_start: 
            print(f'{ind_max_norm = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/max(stack[2][i][1][:,2]), c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3]/max(stack[2][i][1][:,3]), c='Red')
        plt.show()
    
    # compute global max values
    if max_norm or minmax_norm or all_true:
        black_max = -1
        blue_max = -1
        red_max = -1
        for i in range(cells):
            blue_max = max(blue_max, max(stack[2][i][1][:,2]))
            red_max = max(red_max, max(stack[2][i][1][:,3]))
            black_max = max(black_max, max(stack[2][i][1][:,4]))
    
    # plot globally max-normalized intensities
    if glo_max_norm_hist or all_true:
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
    if glo_max_norm or all_true:
        if print_plot_start: 
            print(f'{glo_max_norm = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,2]/blue_max, c='Blue')
            plt.plot(stack[2][i][1][:,0], stack[2][i][1][:,3]/red_max, c='Red')
        plt.show()
    
    # plot the terminal frame globally max-normallized short vs long reporters
    if glo_max_scatter or all_true:
        if print_plot_start: 
            print(f'{glo_max_scatter = }')
        for i in range(cells):
            plt.scatter(stack[2][i][1][len(stack[2][i][1])-1, 2]/blue_max, stack[2][i][1][len(stack[2][i][1])-1, 3]/red_max, s=(stack[2][i][1][len(stack[2][i][1])-1, 4]**ex)/black_max, c='Black', alpha=scatter_alpha)
        plt.title('Globally max-normallized short vs long intensities')
        plt.show()
    
    # plot the strict (axis from 0 to 1) terminal frame globally max-normallized short vs long reporters
    if glo_max_scatter_strict or all_true:
        if print_plot_start: 
            print(f'{glo_max_scatter_strict = }')
        for i in range(cells):
            plt.scatter(stack[2][i][1][len(stack[2][i][1])-1, 2]/blue_max, stack[2][i][1][len(stack[2][i][1])-1, 3]/red_max,  s=(1/20)*(stack[2][i][1][len(stack[2][i][1])-1, 4]**ex)/black_max, c='Black', alpha=scatter_alpha)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.title('Globally max-normallized short vs long intensities')
        plt.show()
    
    # compute global min values
    if max_norm or minmax_norm or all_true:
        black_min = 999_999_999
        blue_min = 999_999_999
        red_min = 999_999_999
        for i in range(cells):
            blue_min = min(blue_min, min(stack[2][i][1][:,2]))
            red_min = min(red_min, min(stack[2][i][1][:,3]))
            black_min = min(black_min, min(stack[2][i][1][:,4]))
    
    # plot globally minmax-normalized intensities
    if glo_minmax_norm_histone or all_true:
        if print_plot_start: 
            print(f'{glo_minmax_norm_histone = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.ylim(0, 1)
            plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,4] - black_min)/(black_max - black_min), c='Black')
            plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min), c='Blue')
            plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,3] - red_min)/(red_max - red_min), c='Red')
        plt.show()
    
    # plot globally minmax-normalized reporters
    if glo_minmax_norm or all_true:
        if print_plot_start: 
            print(f'{glo_minmax_norm = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.ylim(0, 1)
            plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min), c='Blue')
            plt.plot(stack[2][i][1][:,0], (stack[2][i][1][:,3] - red_min)/(red_max - red_min), c='Red')
        plt.show()
    
    # plot the terminal frame globally minmax-normallized short vs long reporters
    if glo_minmax_scatter or all_true:
        if print_plot_start: 
            print(f'{glo_minmax_scatter = }')
        for i in range(cells):
            plt.scatter((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min), (stack[2][i][1][len(stack[2][i][1])-1, 3] - red_min)/(red_max - red_min), s=50*(stack[2][i][1][len(stack[2][i][1])-1, 4] - black_min)/(black_max - black_min), c='Black', alpha=scatter_alpha)
        plt.title('Globally minmax-normallized short vs long intensities')
        plt.show()
    
    # plot the terminal frame globally minmax-normallized short vs long reporters
    if glo_minmax_scatter_strict or all_true:
        if print_plot_start: 
            print(f'{glo_minmax_scatter_strict = }')
        for i in range(cells):
            plt.scatter((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min), (stack[2][i][1][len(stack[2][i][1])-1, 3] - red_min)/(red_max - red_min), s=50*(stack[2][i][1][len(stack[2][i][1])-1, 4] - black_min)/(black_max - black_min), c='Black', alpha=scatter_alpha)
        plt.title('Globally minmax-normallized short vs long intensities')
        plt.ylim(0, 1)
        plt.xlim(0, 1)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.show()
    
    # plot the short/long globally minmax-normallized intensity ratios 
    if minmax_ratio or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], ((stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min))/((stack[2][i][1][:,3] - red_min)/(red_max - red_min) + 0.000001), c='Orange')
        plt.show() 
    
    if same_scale_timeseries:
        ymin = 999_999_999
        ymax = -ymin
        for i in range(cells):
            ymin = min(ymin, min(((stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min))/((stack[2][i][1][:,3] - red_min)/(red_max - red_min) + 0.000001)))
            ymax = max(ymax, max(((stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min))/((stack[2][i][1][:,3] - red_min)/(red_max - red_min) + 0.000001)))
    
    if minmax_ratio or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            if same_scale_timeseries:
                plt.ylim(ymin, ymax)
            plt.plot(stack[2][i][1][:,0], ((stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min))/((stack[2][i][1][:,3] - red_min)/(red_max - red_min) + 0.000001), c='Orange')
        plt.show() 
    
    # plot the long/short globally minmax-normallized intensity ratios 
    if minmax_ratio_rev or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio_rev = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            plt.plot(stack[2][i][1][:,0], ((stack[2][i][1][:,3] - red_min)/(red_max - red_min)) / ((stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min)), c='Orange')
        plt.show()
    
    if same_scale_timeseries:
        ymin = 999_999_999
        ymax = -ymin
        for i in range(cells):
            ymin = min(ymin, min(((stack[2][i][1][:,3] - red_min)/(red_max - red_min))/((stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min) + 0.000001)))
            ymax = max(ymax, max(((stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min))/((stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min) + 0.000001)))
    
    if minmax_ratio_rev or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio_rev = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            if same_scale_timeseries:
                plt.ylim(ymin, ymax)
            plt.plot(stack[2][i][1][:,0], ((stack[2][i][1][:,3] - red_min)/(red_max - red_min)) / ((stack[2][i][1][:,2] - blue_min)/(blue_max - blue_min)), c='Orange')
        plt.show()
    
    # plot the short/long globally minmax-normallized intensity ratio histogram
    if minmax_distrib or all_true:
        if print_plot_start: 
            print(f'{minmax_distrib = }')
        temp = []
        for i in range(cells):
            temp.append(((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min))/((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min)))
        plt.hist(temp)
        plt.title('Globally minmax-normallized short/long ratio distribution')
        plt.show()
    
    # plot the log10() short/long globally minmax-normallized intensity ratios
    if minmax_ratio_log or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio_log = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.axis('off')
            temp = []
            for j in range(len((stack[2][i][1]))):
                temp.append(math.log10(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
            plt.plot(stack[2][i][1][:,0], temp, c='Orange')
        plt.show()
    
    if same_scale_timeseries:
        ymin = 999_999_999
        ymax = -ymin
        for i in range(cells):
            for j in range(len((stack[2][i][1]))):
                ymin = min(ymin, math.log10(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
                ymax = max(ymax, math.log10(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
    
    if minmax_ratio_log or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio_log = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            if same_scale_timeseries:
                plt.ylim(ymin, ymax)
            plt.axis('off')
            temp = []
            for j in range(len((stack[2][i][1]))):
                temp.append(math.log10(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
            plt.plot(stack[2][i][1][:,0], temp, c='Orange')
        plt.show()
    
    
    # plot the log10() of short/long globally minmax-normallized intensity ratio histogram
    if minmax_distrib_log or all_true:
        if print_plot_start: 
            print(f'{minmax_distrib_log = }')
        temp = []
        for i in range(cells):
            temp.append(math.log10(((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min)) / ((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min))))
        plt.hist(temp)
        plt.title('Globally minmax-normallized short/long ratio distribution (log10)')
        plt.show()
    
    # plot the log10() of short/long globally minmax-normallized intensity ratio histogram
    if minmax_distrib_log_cbin or all_true:
        if print_plot_start: 
            print(f'{minmax_distrib_log_cbin = }')
        temp = []
        for i in range(cells):
            temp.append(math.log10(((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min)) / ((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min))))
        plt.hist(temp, bins = cbins)
        plt.title('Globally minmax-normallized short/long ratio distribution (log10)')
        plt.show()
    
    # plot the arctan() short/long globally minmax-normallized intensity ratios
    if minmax_ratio_arctan or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio_arctan = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.ylim(0, math.pi/2)
            plt.axis('off')
            temp = []
            for j in range(len((stack[2][i][1]))):
                temp.append(math.atan(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
            plt.plot(stack[2][i][1][:,0], temp, c='Orange')
        plt.show()
    
    # plot the arctan() of short/long globally minmax-normallized intensity ratio histogram
    if minmax_distrib_arctan or all_true:
        if print_plot_start: 
            print(f'{minmax_distrib_arctan = }')
        temp = []
        for i in range(cells):
            temp.append(math.atan(((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min)) / ((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min))))
        plt.hist(temp)
        plt.title('Globally minmax-normallized short/long ratio distribution (arctan)')
        plt.show()
    
    # plot the arctan() of short/long globally minmax-normallized intensity ratio histogram
    if minmax_distrib_arctan_cbin or all_true:
        if print_plot_start: 
            print(f'{minmax_distrib_arctan_cbin = }')
        temp = []
        for i in range(cells):
            temp.append(math.atan(((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min)) / ((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min))))
        plt.hist(temp, bins = cbins)
        plt.title('Globally minmax-normallized short/long ratio distribution (arctan)')
        plt.show()
    
    # plot a scatter of log10() short/long globally minmax-normallized intensity ratio against terminal frame histone intensity
    if glo_minmax_ratio_scatter or all_true:
        if print_plot_start: 
            print(f'{glo_minmax_ratio_scatter = }')
        for i in range(cells):
            plt.scatter(math.log10(((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min)) / ((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min))), (stack[2][i][1][len(stack[2][i][1])-1, 4] - black_min)/(black_max - black_min), c='Black', alpha=scatter_alpha)
        plt.title('log10(short/long) vs terminal histone intensity')
        plt.show()
    
    # plot a scatter of log10() short/long globally minmax-normallized intensity ratio against average histone intensity
    if glo_minmax_ratio_scatter_avghis or all_true:
        if print_plot_start: 
            print(f'{glo_minmax_ratio_scatter_avghis = }')
        for i in range(cells):
            plt.scatter(math.log10(((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min)) / ((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min))), statistics.mean((stack[2][i][1][:, 4] - black_min)/(black_max - black_min)), c='Black', alpha=scatter_alpha)
        plt.title('log10(short/long) vs average histone intensity')
        plt.show()
    
    # plot a scatter of log10() short/long globally minmax-normallized intensity ratio against average short intensity
    if glo_minmax_ratio_scatter_avgshort or all_true:
        if print_plot_start: 
            print(f'{glo_minmax_ratio_scatter_avgshort = }')
        for i in range(cells):
            plt.scatter(math.log10(((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min)) / ((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min))), statistics.mean((stack[2][i][1][:, 2] - blue_min)/(blue_max - blue_min)), c='Blue', alpha=scatter_alpha)
        plt.title('log10(short/long) vs average short intensity')
        plt.show()
    
    # plot a scatter of log10() short/long globally minmax-normallized intensity ratio against average long intensity
    if glo_minmax_ratio_scatter_avglong or all_true:
        if print_plot_start: 
            print(f'{glo_minmax_ratio_scatter_avglong = }')
        for i in range(cells):
            plt.scatter(math.log10(((stack[2][i][1][len(stack[2][i][1])-1, 2] - blue_min)/(blue_max - blue_min)) / ((stack[2][i][1][len(stack[2][i][1])-1, 3] + 0.000001 - red_min)/(red_max - red_min))), statistics.mean((stack[2][i][1][:, 3] - red_min)/(red_max - red_min)), c='Red', alpha=scatter_alpha)
        plt.title('log10(short/long) vs average long intensity')
        plt.show()
    
    # plot the arctan() short/long globally minmax-normallized intensity ratios with threshold lines
    if minmax_ratio_arctan_wthr or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio_arctan_wthr = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.ylim(0, math.pi/2)
            plt.axis('off')
            temp = []
            for j in range(len((stack[2][i][1]))):
                temp.append(math.atan(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
            plt.plot(stack[2][i][1][:,0], temp, c='Orange')
            plt.axhline(y = atan_thr1, color = 'Black', linestyle = '-', linewidth = 0.5) 
            plt.axhline(y = atan_thr2, color = 'Black', linestyle = '-', linewidth = 0.5) 
        plt.show()
    
    # plot the arctan() short/long globally minmax-normallized intensity ratios with threshold lines
    if minmax_ratio_arctan_wthr or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio_arctan_wthr = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.ylim(0, math.pi/2)
            plt.axis('off')
            temp = []
            for j in range(round(len((stack[2][i][1])) * (1 - cut)), len((stack[2][i][1]))):
                temp.append(math.atan(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
            plt.plot(stack[2][i][1][round(len((stack[2][i][1]))*1/2):,0], temp, c='Orange')
            plt.axhline(y = atan_thr1, color = 'Black', linestyle = '-', linewidth = 0.5) 
            plt.axhline(y = atan_thr2, color = 'Black', linestyle = '-', linewidth = 0.5) 
        plt.show()
        
    # plot the arctan() short/long globally minmax-normallized intensity ratios with threshold lines (markers only)
    if minmax_ratio_arctan_wthr or all_true:
        if print_plot_start: 
            print(f'{minmax_ratio_arctan_wthr = }')
        for i in range(cells):
            plt.subplot(N, M, i + 1)
            plt.ylim(0, math.pi/2)
            plt.axis('off')
            temp = []
            for j in range(round(len((stack[2][i][1])) * (1 - cut)), len((stack[2][i][1]))):
                temp.append(math.atan(((stack[2][i][1][j, 2] + 0.000001 - blue_min)/(blue_max - blue_min))/((stack[2][i][1][j, 3] + 0.000001 - red_min)/(red_max - red_min))))
            plt.plot(stack[2][i][1][round(len((stack[2][i][1]))*1/2):,0], temp, marker='o', markersize=1, linestyle='-', linewidth = 0.2, c='Orange')
            plt.axhline(y = atan_thr1, color = 'Black', linestyle = '-', linewidth = 0.5) 
            plt.axhline(y = atan_thr2, color = 'Black', linestyle = '-', linewidth = 0.5) 
        plt.show()
        