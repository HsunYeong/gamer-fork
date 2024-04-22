import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys

file_dir = '../Record__Correlation/'
files = glob.glob(file_dir + 'correlation_function_t=*.txt')

filename = np.array(files)

r = []
C = []
t = []
for f in filename:
   Corr = np.loadtxt(f, skiprows=1, dtype=float)
   if not r:
      r.append(Corr[:,0])
   else:
      if not np.array_equal(r[0], Corr[:,0]):
         print('radius bin not matched!! filename=\"%s\"'%f)
         sys.exit(1)

   C.append(Corr[:,1])
   match = re.search(r'correlation_function_t=(\d+\.\d+e[+-]\d+)', f)
   if match:
      time = float(match.group(1))
      t.append(time)
   else:
      print('time pattern not matched!! filename=\"%s\"'%f)
      sys.exit(1)


t = np.array(t)
C = np.array(C)

color = plt.cm.turbo(np.linspace(0.1, 0.9, len(r[0])))

plt.figure()
for i in range(len(r[0])):
   plt.plot(t, C[:,i], c=color[i], label = 'r = %1.3e'%(r[0][i]))
plt.xlabel('$t$ (Myr)')
plt.ylabel(r'$C(t)$')
plt.legend(bbox_to_anchor=(1.03,0.03), loc='lower left')
plt.savefig('fig_correlation.png', dpi = 150, bbox_inches="tight")
plt.close()


