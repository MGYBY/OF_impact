# import readforce function from fluidfoam package
from fluidfoam.readpostpro import readforce
import matplotlib.pyplot as plt
import numpy as np
import csv
import os

# sol = './postProcessing/forces/0/force_0.dat'
sol = './'

# force = readforce(sol, time_name = '0', name='force_0')
force = readforce(sol, time_name = '0', name='force')

plt.figure()

plt.plot(force[:, 0], force[:, 1])

# Setting axis labels
plt.xlabel('t (s)')
plt.ylabel('force')

# add grid
plt.grid()
plt.savefig("force.png")
plt.show()

file_name = "force_formatted"
if (not os.path.exists('./%s' % file_name)):
    with open(file_name,'w') as file:
        print("Writing filtered data to file ... ...")
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(zip(np.transpose(force[:, 0]),np.transpose(force[:, 1])))