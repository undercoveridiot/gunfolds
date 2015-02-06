import seaborn as sb
import zickle as zkl
import graphkit as gk
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append('./tools/')
densities = [0.1]#0.2, 0.25, 0.3, 0.35]



d01 = zkl.load('neptune_nodes_8_samples_1000_noise_0.1_OCE_b_svar.zkl')
d10 = zkl.load('saturn_nodes_8_samples_1000_noise_1.0_OCE_b_svar.zkl')

d = d01
density = np.sort(d.keys())
OE = [[gk.oerror(x) for x in d[dd]] for dd in density]
COE = [[gk.cerror(x) for x in d[dd]] for dd in density]


shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

plt.figure(figsize=[10,2])

# g = sb.boxplot(alltimes_old,names=map(lambda x: str(int(x*100))+"%",
#                                      densities),
#               widths=wds, color="Reds", fliersize=fliersz, linewidth=lwd,
#               **{'positions':np.arange(len(densities))-shift,
#                  'label':'naive approach'})

g = sb.boxplot(OE,names=map(lambda x: str(int(x*100))+"%",
                                      densities),
               widths=wds, color="Blues",fliersize=fliersz,
               linewidth=lwd,
               **{'positions':np.arange(len(densities))+shift,
                  'label':'MSL'})

# plt.plot(np.arange(len(densities))-shift,
#         map(np.median,alltimes_old), 'ro-', lw=0.5, mec='k')
# plt.plot(np.arange(len(densities))+shift,
#         map(np.median,alltimes_new), 'bo-', lw=0.5, mec='k')
g.figure.get_axes()[0].set_yscale('log')
plt.xlabel('density (% of 64 total possible edges)')
plt.ylabel('Omission error')
plt.title('100 6 node graphs per density\n$G_2 \\rightarrow G_1$',
          multialignment='center')
plt.subplots_adjust(right=0.99, left=0.2)
plt.legend(loc=0)
plt.show()
