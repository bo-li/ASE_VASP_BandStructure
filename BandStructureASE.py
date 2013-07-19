from ase.structure import bulk
from ase.dft.kpoints import *
import numpy as np
import os
from ase.calculators.vasp import Vasp
from ase.lattice.spacegroup import crystal
import matplotlib.pyplot as plt


#Define the coordinates
a = 5.459
si = crystal('Si', [(0,0,0)], spacegroup=227, cellpar=[a, a, a, 90, 90, 90])
# Define KPOINTS for special points
si_s = bulk('Si', 'diamond', a=5.459)
points = ibz_points['fcc']
G = points['Gamma']
X = points['X']
W = points['W']
K = points['K']
L = points['L']
kpoints, x, X = get_bandpath([W, L, G, X, W, K], si_s.cell,10)
# Point names for plot axis
point_names = ['W','L','G','W','K']

# FROM HERE ON NO REAL INPUT IS NEEDED, UNLESS ONE WISHES TO
# CHANGE SOME VASP INCAR PARAMETERS

# Define the first calculation
os.system('rm KPOINTS')
calc_single = Vasp(system = "Generic System Name",
               istart = 0,iniwav = 1,icharg = 0,
               prec="Accurate", lreal = False, algo = "Normal", encut = 500.00,
               nelm = 200, ediff = 1e-6, gga = "PS",kpts=(4,4,4),
               ediffg = 1e-3, nsw = 0, ibrion = 1, isif = 3, isym = 1,
               ismear = -5)

si.set_calculator(calc_single)
energy = si.get_potential_energy()
os.system('rm KPOINTS')
print("Energy: ",energy)
print("Moving to band structure calculation")
#Dirty printing of the KPOINTS file
f_handle = file('KPOINTS', 'a')
f_handle_b = file('tmp', 'w')
f_handle.write('KPOINTS FILE\n')
f_handle.write(str(len(kpoints)))
f_handle.write('\n')
f_handle.write('Reciprocal \n')
np.savetxt(f_handle_b, kpoints)
f_handle_b.close()
f_handle.close()
os.system('awk \'{printf \"%s %10.8f %10.8f %10.8f %5.2f \\n\","  ",$1,$2,$3,1}\' tmp >> KPOINTS')
os.system('rm tmp')
# Define the band structure calculation
calc_band = Vasp(system = "Generic System Name X",
               istart = 1,iniwav = 1,icharg = 11,
               prec="Accurate", lreal = False, algo = "Normal", encut = 500.00,
               nelm = 200, ediff = 1e-6, gga = "PS",
               ediffg = 1e-3, nsw = 0, ibrion = 1, isif = 3, isym = 1,
               ismear = 0, sigma = 0.05)

si.set_calculator(calc_band)

#si.get_potential_energy()
bands = si.get_band_plot()

# Realign to Fermi level
for line in open('OUTCAR', 'r'):
 if line.rfind('E-fermi') > -1:
  Ef=float(line.split()[2])
print("Fermi level: ",Ef)
for a in range (len(bands)):
 for y in range (len(bands[0])):
  bands[a,y] = bands[a,y] - Ef

# Set up the plotting action
print("The range range (0-",len(bands[0]),")")
Range_L = int(raw_input("Which is the lower bound of the band range to plot? "))
Range_U = int(raw_input("Which is the upper bound of the band range to plot? "))

e_min = bands[:,Range_L].min() - 1
e_max = bands[:,Range_U].max() + 1
for n in range(Range_L,Range_U):
 plt.plot(x, bands[:,n],'b')
for p in X:
 plt.plot([p,p],[e_min,e_max], 'k-')
plt.plot([0,X[-1]],[0,0],'k-')
plt.xticks(X, ['$%s$' % n for n in point_names])
plt.axis(xmin=0, xmax=X[-1], ymin=e_min, ymax=e_max)
plt.show()
