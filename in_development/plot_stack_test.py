""" import libraries """

import math
import os
import numpy as np
import matplotlib.pyplot as plt
import ProcessEmbryo

colors = ["black", "red", "blue", "green"]
stack_path = os.path.join(os.getcwd(), '230212_stack6')
proc_stack = ProcessEmbryo.get_stack(os.path.join(stack_path, 'processed.stack'))

    
def minimal_rectangle(N):
    """
    Calculates the dimensions of a minimal X by Y rectangular grid that can contain N elements
    """
    X = int(round(math.sqrt(N), 0))
    Y = N // X
    
    if N % X > 0:
        Y += 1
    
    return [X, Y]

def list_transform(in_list):
    """
    Conversts the joint [in_list] into 2 single-variable lists [list_x];[list_y]
    """
    list_x = []
    list_y = []
    
    for obj in in_list:
        list_x.append(obj[0])
        list_y.append(obj[1])
    
    return [list_x, list_y]


stack_name = proc_stack[2][3][0]
stack_data = proc_stack[0]
channels = len(proc_stack[2][2][0])
channel_mins = proc_stack[2][2][0]
channel_maxs = proc_stack[2][2][1]

plot_channels = []
for c in range(channels):
    channel_tmp = []
    for lineage in stack_data:
        lineage_tmp = []
        for cell in lineage[1]:
            if (c + 5) <= (len(cell) - 1):
                temp = [cell[0], cell[c + 5]]
                lineage_tmp.append(temp)
        channel_tmp.append(lineage_tmp)
    plot_channels.append(channel_tmp)


test = plot_channels[0]

test_traces = []
for trace in test:
    temp_x = []
    temp_y = []
    for cell in trace:
        temp_x.append(cell[0])
        temp_y.append((cell[1] - channel_mins[0]) / (channel_maxs[0] - channel_mins[0]))
    test_traces.append([temp_x, temp_y])
    
for trace in test_traces:
    plt.plot(trace[0], trace[1], color = "k", linewidth = 0.1)
    
plt.ylabel("Normalized Histone")
plt.xlabel("Frame")
plt.show()

nanog = plot_channels[1]

nanog_traces = []
for trace in nanog:
    temp_x = []
    temp_y = []
    for cell in trace:
        temp_x.append(cell[0])
        temp_y.append((cell[1] - channel_mins[1]) / (channel_maxs[1] - channel_mins[1]))
    nanog_traces.append([temp_x, temp_y])
    
for trace in nanog_traces:
    plt.plot(trace[0], trace[1], color = "r", linewidth = 0.1)
    
plt.ylabel("Normalized Nanog")
plt.xlabel("Frame")
plt.show()

gata6 = plot_channels[2]

gata6_traces = []
for trace in gata6:
    temp_x = []
    temp_y = []
    for cell in trace:
        temp_x.append(cell[0])
        temp_y.append((cell[1] - channel_mins[2]) / (channel_maxs[2] - channel_mins[2]))
    gata6_traces.append([temp_x, temp_y])
    
for trace in gata6_traces:
    plt.plot(trace[0], trace[1], color = "b", linewidth = 0.1)
    
plt.ylabel("Normalized Gata6")
plt.xlabel("Frame")
plt.show()