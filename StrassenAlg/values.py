import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import csv     
matplotlib.rc('xtick', labelsize=8)     
matplotlib.rc('ytick', labelsize=8)


# standard = [0.0004, 0.0006, 0.002, 0.0239, 0.1504, 1.5296, 13.8573, 99.9051, 717.519, 15928.7, 240078]
# strassen = [0.0007, 0.0078, 0.0265, 0.2088, 1.4183, 9.2479, 64.5065, 347.411, 2563.22, 20558.1, 140073]
# m_1 = [0.0003, 0.0011, 0.0023, 0.0096, 0.0563, 0.3545, 2.6341, 31.283, 310.546, 3008.77, 89253.2]

# Read in data from file
m = []
standard = []
strassen = []
hybrid1 = []
hybrid2 = []
with open('runtime_data.csv') as csvfile: #, newline=''
    reader = csv.reader(csvfile, delimiter=';')
    first_row = True
    for row in reader:
        if first_row:
            first_row = False
            continue
        m.append(float(row[0]))
        standard.append(float(row[1]))
        strassen.append(float(row[2]))
        hybrid1.append(float(row[3]))
        hybrid2.append(float(row[4]))

# Scale values
standard = [x/1000 for x in standard]
strassen = [x/1000 for x in strassen]
hybrid1 = [x/1000 for x in hybrid1]
hybrid2 = [x/1000 for x in hybrid2]

# Plot
plt.semilogy(m, standard, '-o', label='Standard')
plt.semilogy(m, strassen, '-o', label='Strassen')
plt.suptitle('Runtimes of different\nmatrix-matrix-multiplication algorithms', size = 16)
plt.xlabel('m --> n=2^m', size=10)
plt.ylabel('average runtime [s]', size = 10)
plt.grid()
plt.legend()

fig = plt.gcf()
fig.set_size_inches(8,6)
scriptdir = os.path.dirname(__file__)
plt.savefig(scriptdir+'/compare_std_vs_strassen.png', dpi=100)
plt.clf()
# Plot
plt.semilogy(m, standard, '-o', label='Standard')
plt.semilogy(m, strassen, '-o', label='Strassen')
plt.semilogy(m, hybrid1, '-o', label='Hybrid (min_size=3)')
plt.semilogy(m, hybrid2, '-o', label='Hybrid (min_size=6)')
plt.suptitle('Runtimes of different\nmatrix-matrix-multiplication algorithms', size = 16)
plt.xlabel('m --> n=2^m', size=10)
plt.ylabel('average runtime [s]', size = 10)
plt.grid()
plt.legend()

fig = plt.gcf()
fig.set_size_inches(8,6)
scriptdir = os.path.dirname(__file__)
plt.savefig(scriptdir+'/compare_all.png', dpi=100)
