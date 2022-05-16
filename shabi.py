import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from functools import partial
import pandas as pd
import mpl_scatter_density
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde
from moviepy.video.io.bindings import mplfig_to_npimage

def search_string_in_file(file_name, string_to_search):
    line_number = 0
    list_of_results = []
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            line_number += 1
            if string_to_search in line:
                list_of_results.append((line_number))
    return list_of_results
    

N=search_string_in_file('testall1.txt','Y-Y0')
heatmap=0
def kanade(a,b):
    with open('testall1.txt') as f:
        f=f. readlines()[a: a+b+2]
        
    line0=f[0]
    arrays0 = list(map(float, line0.split()[1:]))
    arrays = [list(map(float, line.split()[2:])) for line in f[2:]]
    bb=np.array(arrays)
    SB0=np.array([(arrays0[2]+bb[i][2])/1000*200000000+arrays0[0]+bb[i][0] for i in range(0,b)])
    SB1=np.array([(arrays0[3]+bb[i][3])/1000*200000000+arrays0[1]+bb[i][1] for i in range(0,b)])
    SB2=np.array([bb[i][4] for i in range(0,50)])
    heatmap, xedges, yedges = np.histogram2d(SB0, SB1, bins=512,range=[[-4e6,4e6],[-4e6,4e6]])
    heatmap =np.flipud(gaussian_filter(heatmap, sigma=2))
    return heatmap

for p in range(0,1):
    heatmap=heatmap+kanade(N[p],199000)


plt.imshow(heatmap,cmap=mpl.cm.nipy_spectral,norm = mpl.colors.Normalize(vmin=0, vmax=100))
plt.colorbar()
np.savetxt('0.1000_149.75.txt',heatmap)
