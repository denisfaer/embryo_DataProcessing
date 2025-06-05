import math
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import ProcessEmbryo

stack_path = os.path.join(os.getcwd(), '230212_stack6')
stack = ProcessEmbryo.process_stack(stack_path)

ProcessEmbryo.save_stack(stack)

out = ProcessEmbryo.get_stack(os.path.join(stack_path, 'processed.stack'))

IDs = []

for cell in out[0]:
    IDs.append(cell[0][1])
    
unique_IDs = list(set(IDs))

print(IDs)
print()
print(len(IDs))
print()
print(len(IDs) == len(unique_IDs))
