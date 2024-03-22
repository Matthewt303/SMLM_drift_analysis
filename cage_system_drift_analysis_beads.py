
"""
Created on Tue Jan 24 15:33:30 2023

@author: Administrator
"""
#%% read all CSV.
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.cluster import dbscan
from matplotlib.pyplot import MultipleLocator #scale of plot setting
import matplotlib.ticker as mticker
import matplotlib.ticker as ticker
import io
from PIL import Image


df = pd.read_csv(r'C:\Users\Administrator\Desktop\raw_data\after_sigma_filter_less_than_130_13.csv')
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
 
#error_x_sd=np.power(sum([np.power((x_sd[i]),2) for i in range(len(x_sd))]),0.5)

x_length=[len(item) for item in x_groupby_label]

#%%color
#col_liss=[['Red']*10,['Orange']*10,['Yellow','Green']*10,['Blue']*10,['Cyan']*10,['Purple']*10,['Black']*10]

#col_list=['Red','Red','Red','Red','Red','Red','Red','green','green','green','green','snow']
 
"""
for i in range(len(x_length)):
    col_i = int(round(i/10.0,0))
    plt.scatter(df['x [nm]'][df['labels']==i],df['y [nm]'][df['labels']==i],alpha=0.1, marker='o',color=col_list[col_i])
plt.xlim(0,35000)
plt.ylim(0,35000)
"""

#%%
"""
i=-1
plt.scatter(df['x [nm]'][df['labels']==i],df['y [nm]'][df['labels']==i],alpha=0.1)
plt.xlim(0,35000)
plt.ylim(0,35000)
"""

#%% Y-axis Drift 
y_drift=[np.subtract(list(y_groupby_label[i]),float(list(y_groupby_label[i])[0])) for i in range(len(y_groupby_label))]
y_sd=[np.std(y_drift[j][1:]) for j in range(len(y_drift))] 
y_mean=[np.mean(y_sd[1:])]

#%% xy coordinate
x_hstack=np.hstack(x_drift[1:]).tolist()
y_hstack=np.hstack(y_drift[1:]).tolist()
xy_coordinate=np.dstack((x_hstack,y_hstack))
xy_coordinate= np.reshape(xy_coordinate,(-1,2))

#%% scatterplot
"""
labels_normal=np.delete(labels, np.where(labels == -1))
#labels_normal=np.reshape(labels_normal,(len(labels_normal),1))
#xy_coordinate=np.concatenate((xy_coordinate,labels_normal),axis=1)
"""
concat_x=np.concatenate(x_drift[1:])
concat_y=np.concatenate(y_drift[1:])

frame_id=[list(range(1,len(item)+1)) for item in x_drift[1:]]
frame_ids=np.concatenate(frame_id)
frame_ids=(frame_ids-1)/2
#%% Plt.Scatter

sc=plt.scatter(concat_x,concat_y,c=frame_ids, cmap=plt.cm.RdYlBu, marker='+',s=30, alpha=0.6)
cb1=plt.colorbar(sc)

"""Setting of the colour bar scale"""
tick_locator = ticker.MaxNLocator(nbins=5)  # colorbar上的刻度值个数 
cb1.locator = tick_locator
cb1.set_ticks([0,20,40,60,80,100,120])
cb1.update_ticks()

cb1.ax.tick_params(width=1,length=3, labelsize=24) 
cb1.set_label('Time (min)',fontdict={'family':'Times New Roman', 'size': 24})
for l in cb1.ax.yaxis.get_ticklabels():  # special way to set the font
        l.set_family('Times New Roman')

"""Setting of the plot scale"""
x_major_locator=MultipleLocator(20)
y_major_locator=MultipleLocator(20)
bwith = 1 #width of scales    5pt = 0.07in
ax=plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.spines['bottom'].set_linewidth(bwith) #width of scales
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)

ax.tick_params(which='major',width=1,length=3)


plt.xlim(-100,100) 
plt.ylim(-100,100)
plt.xlabel(u'x (nm)',fontdict={'family':'Times New Roman', 'size': 24})
plt.ylabel(u'y (nm)',fontdict={'family':'Times New Roman', 'size': 24})
plt.xticks(fontproperties = 'Times New Roman', size = 24)
plt.yticks(fontproperties = 'Times New Roman', size = 24)


ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

#%%
#from matplotlib.lines import Line2D

# Save the image in memory in PNG format
png1 = io.BytesIO()
plt.savefig(png1, format="png", dpi=500, pad_inches = .1, bbox_inches = 'tight')

# Load this image into PIL
png2 = Image.open(png1)

# Save as TIFF
png2.save("ap_qt_ad.tiff")
png1.close()

#%% test_coloring
"""
fig=plt.colorbar(figsize=(8,3))
ax =fig.add_axes([0.0,0.0,0.6,0.8])
bounds=[0,24,48,72,96,120]


for i in range(0,12):
    start=(i*20)+1
    end=((i+1)*20)+1
    bool_test=np.array([item >= start & item <end for item in frame_ids])
    plt.scatter(concat_x[bool_test],concat_y[bool_test],color=col_list[i], marker='+',s=25, alpha=0.1)
"""
    #plt.scatter(np.concatenate(x_drift[1:]),np.concatenate(y_drift[1:]),marker='+',s=25,alpha=0.4)
    #plt.scatter(x_drift[1],y_drift[1], c=i,cmap='coolwarm')
    #color=col_list[i]
#%% absolute distance
#result= pow((pow(np.float_(x_mean),2)+pow(np.float_(y_mean),2)),0.5)  
#result=np.power(np.power(x_mean,2)+np.power(y_mean,2),0.5) 

#%%
"""
item=2
import matplotlib.pyplot as plt
plt.scatter(x[item].index,x_drift[item])
"""
 
#%% Export
"""
import xlsxwriter

workbook = xlsxwriter.Workbook('text.xlsx')
worksheet = workbook.add_worksheet()
row=0

for col, data in enumerate(xy_coordinate):
    worksheet.write_column(row, col, data)

workbook.close()
"""

#%%

df = pd.read_csv(r'C:\Users\Administrator\Desktop\raw_data\after_sigma_filter_less_than_127_4_12.csv')
df.drop(columns=['sigma [nm]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty [nm]'], axis=1, inplace=True)
from sklearn.cluster import dbscan
_, labels = dbscan(df[['x [nm]','y [nm]']], eps=30, min_samples=20) 
df['labels']=labels
x_groupby_label=[grp['x [nm]'] for l, grp in df.groupby('labels')]
y_groupby_label=[grp['y [nm]'] for l, grp in df.groupby('labels')]
x_drift=[np.subtract(list(x_groupby_label[i]),float(list(x_groupby_label[i])[0])) for i in range(len(x_groupby_label))]
#stuff_sd=[np.std(stuff[j][1:]) for j in range(1,len(stuff))]   #two ways to writerange
x_sd=[np.std(x_drift[j][1:]) for j in range(len(x_drift))]
x_mean=[np.mean(x_sd[1:])]
y_drift=[np.subtract(list(y_groupby_label[i]),float(list(y_groupby_label[i])[0])) for i in range(len(y_groupby_label))]
y_sd=[np.std(y_drift[j][1:]) for j in range(len(y_drift))] 
y_mean=[np.mean(y_sd[1:])]

x_hstack=np.hstack(x_drift[1:]).tolist()
y_hstack=np.hstack(y_drift[1:]).tolist()
xy_coordinate=np.dstack((x_hstack,y_hstack))
xy_coordinate= np.reshape(xy_coordinate,(-1,2))
