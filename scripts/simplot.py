import seaborn as sb
import zickle as zkl
import graphkit as gk
import bfutils as bfu
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append('./tools/')


l = ['hooke_nodes_10_samples_1000_density_0.15_noise_0.1_OCE_b_svar.zkl',
     'hooke_nodes_10_samples_1000_density_0.2_noise_0.1_OCE_b_svar.zkl',
     'hooke_nodes_10_samples_1000_density_0.25_noise_0.1_OCE_b_svar.zkl']

d10_01 = {0.15: zkl.load(l[0]),
          0.2: zkl.load(l[1]),
          0.25: zkl.load(l[2])}

d01 = zkl.load('neptune_nodes_8_samples_1000_noise_0.1_OCE_b_svar.zkl')
d10 = zkl.load('leibnitz_nodes_8_samples_1000_noise_1.0_OCE_b_svar_flat.zkl')


def estOE(d):
    gt= d['gt']['graph']
    gt=bfu.undersample(gt,1)
    e = gk.OCE(d['estimate'],gt)
    N = np.double(len(gk.edgelist(gt))) +\
        np.double(len(gk.bedgelist(gt)))    
    return (e['directed'][0]+e['bidirected'][0])/N

def estCOE(d):
    gt= d['gt']['graph']
    gt=bfu.undersample(gt,1)    
    e = gk.OCE(d['estimate'],gt)
    n = len(gt)
    N = np.double(n**2+(n-1)**2/2.0\
                  -len(gk.edgelist(gt))
                  -len(gk.bedgelist(gt)))                  
    return (e['directed'][1]+e['bidirected'][1])/N

d = d10
density = np.sort(d.keys())
OE = [[gk.oerror(x) for x in d[dd]] for dd in density]
COE = [[gk.cerror(x) for x in d[dd]] for dd in density]

eOE = [[estOE(x) for x in d[dd]] for dd in density]
eCOE = [[estCOE(x) for x in d[dd]] for dd in density]


shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

plt.figure(figsize=[10,6])

plt.subplot(121)
g = sb.boxplot(eOE,names=map(lambda x: str(int(x*100))+"%",
                                     density),
              widths=wds, color="Reds", fliersize=fliersz, linewidth=lwd,
              **{'positions':np.arange(len(density))-shift,
                 'label':'$G_2$-space SVAR estimation error'})

g = sb.boxplot(OE,names=map(lambda x: str(int(x*100))+"%",
                                      density),
               widths=wds, color="Blues",fliersize=fliersz,
               linewidth=lwd,
               **{'positions':np.arange(len(density))+shift,
                  'label':'$G_1$-space error'})

# plt.plot(np.arange(len(densities))-shift,
#         map(np.median,alltimes_old), 'ro-', lw=0.5, mec='k')
# plt.plot(np.arange(len(densities))+shift,
#         map(np.median,alltimes_new), 'bo-', lw=0.5, mec='k')
#g.figure.get_axes()[0].set_yscale('log')
plt.ylim([-0.1,1.1])
plt.xlabel('density (% of 100 total possible edges)')
plt.ylabel('Edge omission error')
plt.title('100 10-node graphs per density',
          multialignment='center')
plt.legend(loc=0)

plt.subplot(122)
g = sb.boxplot(eCOE,names=map(lambda x: str(int(x*100))+"%",
                                     density),
              widths=wds, color="Reds", fliersize=fliersz, linewidth=lwd,
              **{'positions':np.arange(len(density))-shift,
                 'label':'$G_2$-space SVAR estimation error'})

g = sb.boxplot(COE,names=map(lambda x: str(int(x*100))+"%",
                                      density),
               widths=wds, color="Blues",fliersize=fliersz,
               linewidth=lwd,
               **{'positions':np.arange(len(density))+shift,
                  'label':'$G_1$-space error'})

# plt.plot(np.arange(len(densities))-shift,
#         map(np.median,alltimes_old), 'ro-', lw=0.5, mec='k')
# plt.plot(np.arange(len(densities))+shift,
#         map(np.median,alltimes_new), 'bo-', lw=0.5, mec='k')
#g.figure.get_axes()[0].set_yscale('log')
plt.ylim([-0.1,1.1])
plt.xlabel('density (% of 100 total possible edges)')
plt.ylabel('Edge comission error')
plt.title('100 10-node graphs per density',
          multialignment='center')
plt.legend(loc=0)

plt.subplots_adjust(right=0.99, left=0.2)

plt.show()
