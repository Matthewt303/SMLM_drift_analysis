# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 10:53:02 2023

@author: Administrator
"""
#%% read all CSV.
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt


#%%
df = pd.read_csv(r'C:\Users\Administrator\Desktop\commerical_Elyra\92_to_130_sigma\11.csv')
df.drop(columns=['sigma [nm]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty [nm]'], axis=1, inplace=True)
#%% DBscan
"""
import numpy as np
from sklearn.cluster import dbscan
 # the following is a simulation of two clusters (beads) with 150 points each (xy coordinates) 
ds = np.concatenate([np.random.randn(150, 2), np.random.randn(150, 2)+ [20, 8]])
 # dbscan has two parameters: eps is the maximum distance between two neighbours (larger values give larger clusters) (you could try eps=10 to eps=50)
 # min_samples is the minimum number of point per cluster, this can be used to avoid clustering spurious detection (I would you 80% of the number of frames 
 # i.e. 250 frames -> min_samples=200) 
 # every point gets a label; points that are not in a cluster get the label "-1". 
_, labels = dbscan(ds, eps=2, min_samples=100)
 #plot the results
import matplotlib.pyplot as plt
plt.scatter(*ds.T, c=labels)
"""

from sklearn.cluster import dbscan
_, labels = dbscan(df[['x [nm]','y [nm]']], eps=100, min_samples=200)  
                                                                       
df['labels']=labels

#%% Groupby labels
x_groupby_label=[grp['x [nm]'] for l, grp in df.groupby('labels')]
y_groupby_label=[grp['y [nm]'] for l, grp in df.groupby('labels')]
#x_groupby_label=[(l,grp['x [nm]')] for l, grp in df.groupby('labels')]

#%% X-axis Drift
x_drift=[np.subtract(list(x_groupby_label[i]),float(list(x_groupby_label[i])[0])) for i in range(len(x_groupby_label))]
#stuff_sd=[np.std(stuff[j][1:]) for j in range(1,len(stuff))]   #two ways to writerange
x_sd=[np.std(x_drift[j][1:]) for j in range(len(x_drift))] 
x_mean=[np.mean(x_sd[1:])]

x_length=[len(item) for item in x_groupby_label]

#%% Y-axis Drift 
y_drift=[np.subtract(list(y_groupby_label[i]),float(list(y_groupby_label[i])[0])) for i in range(len(y_groupby_label))]
y_sd=[np.std(y_drift[j][1:]) for j in range(len(y_drift))] 
y_mean=[np.mean(y_sd[1:])]
"""
for i in range(len(x_length)):
    col_i = int(round(i/10.0,0))
    plt.scatter(df['x [nm]'][df['labels']==i],df['y [nm]'][df['labels']==i],alpha=0.1, marker='o')
plt.xlim(0,35000)
plt.ylim(0,35000)
plt.show()
i=-1
plt.scatter(df['x [nm]'][df['labels']==i],df['y [nm]'][df['labels']==i],alpha=0.1)
plt.xlim(0,35000)
plt.ylim(0,35000)
"""
#%% xy coordinate
x_hstack=np.hstack(x_drift[1:]).tolist()
y_hstack=np.hstack(y_drift[1:]).tolist()
xy_coordinate=np.dstack((x_hstack,y_hstack))
xy_coordinate= np.reshape(xy_coordinate,(-1,2))

#%% scatterplot
concat_x=np.concatenate(x_drift[1:])
concat_y=np.concatenate(y_drift[1:])

frame_id=[list(range(1,len(item)+1)) for item in x_drift[1:]]
frame_ids=np.concatenate(frame_id)
frame_ids=(frame_ids-1)/2
#%% Plt.Scatter
from matplotlib.pyplot import MultipleLocator #scale of plot setting
import matplotlib.ticker as mticker
import matplotlib.ticker as ticker

sc=plt.scatter(concat_x,concat_y,c=frame_ids, cmap=plt.cm.RdYlBu, marker='+',s=30, alpha=0.6) #s=30,alpha=0.6
cb1=plt.colorbar(sc)

"""Setting of the colour bar scale"""
tick_locator = ticker.MaxNLocator(nbins=5)  # The numbers of scales in colorbar
cb1.locator = tick_locator
cb1.set_ticks([0,20,40, 60,80, 100, 120])
cb1.update_ticks()

cb1.ax.tick_params(width=1,length=3, labelsize=24) 
cb1.set_label('Time (min)',fontdict={'family':'Times New Roman', 'size': 24})
for l in cb1.ax.yaxis.get_ticklabels():  # special way to set the font
        l.set_family('Times New Roman')

"""Setting of the plot scale"""
x_major_locator=MultipleLocator(100)
y_major_locator=MultipleLocator(100)
bwith = 1 #width of scales    5pt = 0.07in
ax=plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.spines['bottom'].set_linewidth(bwith) #width of scales
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)

ax.tick_params(which='major',width=1,length=3)


plt.xlim(-250,250) 
plt.ylim(-250,250)
plt.xlabel(u'x (nm)',fontdict={'family':'Times New Roman', 'size': 22})
plt.ylabel(u'y (nm)',fontdict={'family':'Times New Roman', 'size': 22})
plt.xticks(fontproperties = 'Times New Roman', size = 22)
plt.yticks(fontproperties = 'Times New Roman', size = 22)


ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
#%% absolute distance
result= pow((pow(np.float_(x_mean),2)+pow(np.float_(y_mean),2)),0.5)  

#%%
#from matplotlib.lines import Line2D
import io
from PIL import Image

# Save the image in memory in PNG format
png1 = io.BytesIO()
plt.savefig(png1, format="png", dpi=500, pad_inches = .1, bbox_inches = 'tight')

# Load this image into PIL
png2 = Image.open(png1)

# Save as TIFF
png2.save("ap_qt_ad.tiff")
png1.close()


