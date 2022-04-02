#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 14:58:40 2018

@author: fujita
"""

from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import mpl_toolkits.axes_grid1
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
import itertools
import matplotlib.cm as cm
import scipy.io
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform

def mkClstAveTimeSign(save_dir,TimeSignDF,ClstMol):###
    ClstAveTimeSeriesDF = pd.DataFrame(index=TimeSignDF.columns,columns=list(ClstMol.keys()))
    for i in range(len(ClstMol.keys())):
        ForClusterMean=np.empty((0,len(TimeSignDF.columns)),int)

        target_dataname = ClstMol['cluster_' + str(i+1)]
        for j in range(len(target_dataname)):
          if i ==0:
              ForClusterMean = np.append(ForClusterMean,np.array([TimeSignDF.loc[target_dataname[j]]]),axis=0)
          else:
              ForClusterMean = np.append(ForClusterMean,np.array([TimeSignDF.loc[target_dataname[j]]]),axis=0)
        ClusterMean=[]
        for ii in range(0,len(ForClusterMean[0,:])):
            ClusterMean.append(np.mean(ForClusterMean[:,ii]))
        ClstAveTimeSeriesDF['cluster_' + str(i+1)] = ClusterMean
    ClstAveTimeSeriesDF.to_excel(save_dir + 'ClstAveTimeSeriesDF.xlsx')
    return(ClstAveTimeSeriesDF)
    
def mkClstStdTimeSign(save_dir,TimeSignDF,ClstMol):
    ClstStdTimeSeriesDF = pd.DataFrame(index=TimeSignDF.columns,columns=list(ClstMol.keys()))
    for i in range(len(ClstMol.keys())):
        ForClusterStd=np.empty((0,len(TimeSignDF.columns)),int)

        target_dataname = ClstMol['cluster_' + str(i+1)]
        for j in range(len(target_dataname)):
          if i ==0:
              ForClusterStd = np.append(ForClusterStd,np.array([TimeSignDF.loc[target_dataname[j]]]),axis=0)
          else:
              ForClusterStd = np.append(ForClusterStd,np.array([TimeSignDF.loc[target_dataname[j]]]),axis=0)
        ClusterStd=[]
        for ii in range(0,len(ForClusterStd[0,:])):
            ClusterStd.append(np.nanstd(ForClusterStd[:,ii]))
        ClstStdTimeSeriesDF['cluster_' + str(i+1)] = ClusterStd
    ClstStdTimeSeriesDF.to_excel(save_dir + 'ClstStdTimeSeriesDF.xlsx')
    return(ClstStdTimeSeriesDF)

def ChangeDictToDF_MolColor(ClstMolDict,ColorDict):
    import pandas as pd
    from itertools import chain
    DF = pd.DataFrame(index=list(chain.from_iterable(list(ClstMolDict.values()))),columns=['ClsterNo.','ClsterColor'])
    for i in range(len(ClstMolDict.values())):
        DF.loc[list(ClstMolDict.values())[i],'ClsterNo.']=list(ClstMolDict.keys())[i]
    for j in range(len(ClstMolDict.keys())):
        try:
            DF.loc[DF['ClsterNo.']==list(ClstMolDict.keys())[j],'ClsterColor'] =ColorDict[list(ClstMolDict.keys())[j]]
        except:
            DF.loc[DF['ClsterNo.']==list(ClstMolDict.keys())[j],'ClsterColor']='nan'
            
    return(DF)
    

def ChangeDictToDF(Dict):###
    import pandas as pd
    IdxNum= len(list(Dict.values())[0])
    for i in range(1,len(Dict.values())):
        
        if (len(list(Dict.values())[i]) > len(list(Dict.values())[i-1])) and (len(list(Dict.values())[i]) > IdxNum):
            IdxNum = len(list(Dict.values())[i])
            

    DF = pd.DataFrame(index = list(range(0,int(IdxNum))), columns=list(Dict.keys()))
    for i in range(len(Dict.keys())):
        DF[list(Dict.keys())[i]]=pd.Series(list(Dict[list(Dict.keys())[i]]))
    return(DF)

def draw_heatmapclustering(a,MolColorSwitch,ColorSwitch,Optindict,save_dir,cmap='PuOr_r'):###
    from matplotlib import pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage, dendrogram
    import mpl_toolkits.axes_grid1
    from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
    import itertools
    import matplotlib.cm as cm
    import scipy.io
    import re
    
    metric = Optindict['metric']
    method = 'ward'#'centroid'#'ward'#'average'
    threshold_ratio=0.7
    plt.rcParams['font.family'] = 'Arial Unicode MS'
    MolColor = Optindict['MolColor']
    colorlist = ['purple','blue','magenta','deepskyblue','red','orange','gold','cyan','mediumspringgreen','greenyellow','lime']

    vmin = min(a.values.flatten());vmax = max(a.values.flatten())
    v = max(np.abs(vmin), np.abs(vmax))
    set_link_color_palette(colorlist)
    ih = 70
    wl = 120
    fig, ax = plt.subplots( figsize=(ih,wl))#

    main_axes = plt.gca()

    divider = make_axes_locatable(main_axes)
    
    ydendro_axes = divider.append_axes("left",size=5 , pad=0.25)

    plt.sca(ydendro_axes)
    threshold,numcluster,distance_sort =  Optindict['Threshold'] ,Optindict['numcluster'],'descending'#
    

    plt.rcParams['font.size'] = 3
    #
    ylinkage = linkage(pdist(a.astype(float), metric=metric), method=method,
                       metric=metric)
    ydendro = dendrogram(ylinkage, orientation='left',leaf_font_size = 0.01,no_labels=1,
                         color_threshold=threshold,above_threshold_color='black',#
                         count_sort='False'
                         )
    ydendro_axes.invert_yaxis()

    plt.rcParams['lines.linewidth'] = 10
    plt.rcParams['axes.linewidth'] = 0.4

    plt.axvline(x=threshold, c='k',lw=0.5,ls='dashed')  

    plt.gca().tick_params(axis='x', labelsize=100)
    a = a.loc[[a.index[i] for i in ydendro['leaves']]]
    plt.sca(ax)
    a=a.astype(float)
    if (cmap == 'Reds') or (cmap == 'afmhot_r'):
        aa = ax.imshow(a,  aspect='equal', interpolation='nearest', cmap=cmap, vmin=0, vmax=v)
    else:
        aa = ax.imshow(a,  aspect='auto', interpolation='nearest', cmap=cmap, vmin=-v, vmax=v)
    ys, xs = np.meshgrid(range(a.shape[0]),range(a.shape[1]),indexing='ij')

    cluster_assign = fcluster(ylinkage,5,
                                 criterion='distance')
    clustered = fcluster(ylinkage, numcluster,
                         criterion='maxclust'
                         )
    cluster_along_id =clustered[ydendro['leaves']]
  

    for ii in range(1,numcluster+1):  
       cluster_pos = np.where(cluster_along_id==ii)[0][0]
       plt.axhline(y=cluster_pos-0.5, lw=15,color='k')  


    char = plt.colorbar(aa,pad=0.2,shrink=0.01)#

    char.ax.tick_params(labelsize=175)#
    
    plt.gca().yaxis.tick_right()#

    xl = plt.xticks(range(a.shape[1]), a.columns, rotation=270, size=130,color='k')#
    plt.xlabel('Time(min.)',size=10)
    ax.set_xticklabels(a.columns,minor=False)

    plt.tick_params(axis='both',direction='in',width=0.05)
   
    bottom, top = ax.get_ylim() 
    plt.yticks(range(a.shape[0]), list(a.index), size=100)#
    plt.ylim(top,bottom)

    b=pd.concat([a,MolColor.reindex(index=a.index)],axis=1)
    b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')

    if MolColorSwitch == 1:
        try:
            MolColor = Optindict['MolColor']
        except:
            pass
        if ('Glucose' in list(a.index)) or ('Glucose (mg/dL)' in list(a.index)):
            if '_' in list(a.index)[0] :
                NewColorList=[]
                for ij in list(a.index):
                    r = re.compile("(.*)(_)(.*)"); d = r.search(ij); NewColorList += [MolColor['MolColor'][d.group(1)] ]              
                [t.set_color(i) for (i,t) in zip(NewColorList,ax.yaxis.get_ticklabels())]
            else:
                b=pd.concat([a,MolColor.reindex(index=a.index)],axis=1)
                b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.yaxis.get_ticklabels())]
        elif len(a.index) < len(a.columns) or 'Glucose' not in list(a.index):
            if '_' in list(a.columns)[0] :
                NewColorList=[]
                for ij in list(a.columns):
                    r = re.compile("(.*)(_)(.*)"); d = r.search(ij); NewColorList += [MolColor['MolColor'][d.group(1)]]               
                [t.set_color(i) for (i,t) in zip(NewColorList,ax.xaxis.get_ticklabels())]  
            else:
                b=pd.concat([a,MolColor.reindex(index=a.index)],axis=1)
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.xaxis.get_ticklabels())]  
    else:
        ColorSwitch=''
    
    plt.gca().xaxis.set_ticks_position('none')
    plt.gca().yaxis.set_ticks_position('none')
    plt.gca().invert_yaxis()
    tempcolor = ydendro['color_list']
    
    plt.savefig(save_dir + 'dendrogram'+ColorSwitch+Optindict['title']+'.pdf',bbox_inches="tight")    
   ## plt.savefig(save_dir + 'dendrogram'+ColorSwitch+Optindict['title']+'.png',bbox_inches="tight")    

    ClstMol={}
    Colordict={}
    ii=0
    for i in range(1,numcluster+1): 
        ClstMol.update({'cluster_' + str(i):[]})
        Colordict.update({'cluster_' + str(i):[]})
        try:
            while cluster_along_id[ii]==cluster_along_id[ii+1]:
                ClstMol['cluster_' + str(i)].append(a.index[ii])
                ii+=1
            ClstMol['cluster_' + str(i)].append(a.index[ii])
            ii+=1
        except:
            ClstMol['cluster_' + str(i)].append(a.index[ii])

    rmnum=[]
    for i in range(len(tempcolor)-1):
        if tempcolor[i]=='black':
            rmnum.append(i)
    for i in rmnum:
        tempcolor[i]=tempcolor[i-1]
    i=0
    tempcolor2=[]
    for ii in range(numcluster):
        while i<len(tempcolor)-2 and tempcolor[i]==tempcolor[i+1]:
                i+=1
        try:
            tempcolor2.append(tempcolor[i])
        except:
            pass
        i+=1
    count=0
    for i in range(numcluster):
        if len(ClstMol['cluster_' + str(1+i)])==1:
            Colordict['cluster_' + str(1+i)].append('black')
        else:
            try:
                Colordict['cluster_' + str(1+i)].append(tempcolor2[count])
                count+=1
            except:
                pass
    
    ClstMolDF = ChangeDictToDF(ClstMol)
    ClstMolDF.to_excel(save_dir + 'ClstMol.xlsx')
    ClstColorDF = ChangeDictToDF(Colordict)
    ClstColorDF.to_excel(save_dir + 'ClstColor.xlsx')
    
    ColorDF =pd.DataFrame()
    
    ClstNoColorDF = ChangeDictToDF_MolColor(ClstMol,Colordict)
    tempDF = MolColor.drop('ClstColor',axis=1)
    NewDF=pd.concat([tempDF,ClstNoColorDF.reindex(index=tempDF.index)],axis=1)
    NewDF.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
    ClstNoColorDF.to_excel(save_dir+'ClstNoColorDF'+Optindict['title']+'.xlsx')
    ColorDF = pd.concat([NewDF,b[['MolColor','Metabo']].reindex(index=NewDF.index)],axis=1)
      
    return(ClstMol,ClstMolDF,Colordict,b,ClstColorDF,ColorDF,ClstNoColorDF)



def draw_heatmap(a,MolColorSwitch,ColorSwitch,Optindict,save_dir,cmap='bwr'):

    metric = 'euclidean'
    method = 'ward'#'average'
    threshold_ratio=0.7
    plt.rcParams['font.family'] = 'Arial Unicode MS'

    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,index_col=0).iloc[:83,:]
    max_TimeSign = max(a.max(axis=1))
    min_TimeSign = min(a.min(axis=1))
    colorlist = ['magenta','deepskyblue','purple','orange','lime','green','red','cyan','blue','gold','mediumspringgreen','black','greenyellow']

    vmin = min(a.values.flatten());vmax = max(a.values.flatten())
    v = max(np.abs(vmin), np.abs(vmax))
    set_link_color_palette(colorlist)
    fig, ax = plt.subplots(figsize=(90,60)) 

    threshold,numcluster,distance_sort = 5.15,13,'descending'
    plt.rcParams['font.size'] = 5

    plt.rcParams['lines.linewidth'] = 5
    plt.rcParams['axes.linewidth'] = 0.1


    a=a.astype(float)
    aa = ax.imshow(a.T,  aspect='equal', interpolation='nearest', cmap=cmap, vmin=-v, vmax=v)

    char = plt.colorbar(aa,pad=0.02)
    char.ax.tick_params(labelsize=120)

    plt.xticks(range(a.shape[0]), a.index, rotation=270, size=40,color='k')

    plt.tick_params(axis='both',direction='in',width=0.05)
    aindex=list(a.columns)
    plt.yticks(range(a.shape[1]),a.columns , size=40)
    b=pd.concat([a,MolColor],axis=1)
    b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')

    
    if MolColorSwitch == 1:    
        if len(a.index) > len(a.columns) or 'Glucose' in list(a.columns):
            if '_' in list(a.columns)[0] :
                NewColorList=[]
                for ij in list(a.columns):
                    r = re.compile("(.*)(_)(.*)"); d = r.search(ij); NewColorList += [MolColorDict[d.group(1)] ]              
                [t.set_color(i) for (i,t) in zip(NewColorList,ax.yaxis.get_ticklabels())]
            else:
                b=pd.concat([a,MolColor],axis=1)
                b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.yaxis.get_ticklabels())]
        elif len(a.index) < len(a.columns):
            if '_' in list(a.columns)[0] :
                NewColorList=[]
                for ij in list(a.columns):
                    r = re.compile("(.*)(_)(.*)"); d = r.search(ij); NewColorList += [MolColorDict[d.group(1)]]               
                [t.set_color(i) for (i,t) in zip(NewColorList,ax.xaxis.get_ticklabels())]  
            else:
                b=pd.concat([a,MolColor],axis=1)
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.xaxis.get_ticklabels())] 
    elif MolColorSwitch == 'both':
                b=pd.concat([a,MolColor],axis=1)
                b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.xaxis.get_ticklabels())]
                ColorListrevse = list(b[ColorSwitch]);
                [t.set_color(i) for (i,t) in zip(ColorListrevse,ax.yaxis.get_ticklabels())]
                xtic = plt.xticks();ytic = plt.yticks()
                turn = sorted(set(list(b['MolColor'])), key=list(b['MolColor']).index)
                pos=0               
    else:
        ColorSwitch=''
    plt.gca().xaxis.set_ticks_position('none')
    plt.gca().yaxis.set_ticks_position('none')
   
    plt.savefig(save_dir + 'heatmap'+ColorSwitch+Optindict['title']+'.pdf')
    ##plt.savefig(save_dir + 'heatmap'+ColorSwitch+Optindict['title']+'.png')
    
    ClstMol={}
    Colordict={}
    ii=0
    for i in range(numcluster):#
        ClstMol.update({'cluster_' + str(numcluster-i):[]})
        Colordict.update({'cluster_' + str(numcluster-i):[]})

        try:
            while cluster_along_id[ii]==cluster_along_id[ii+1]:
                ClstMol['cluster_' + str(numcluster-i)].append(a.index[ii])
                ii+=1

            ClstMol['cluster_' + str(numcluster-i)].append(a.index[ii])

            ii+=1

        except:
            ClstMol['cluster_' + str(numcluster-i)].append(a.index[ii])
    rmnum=[]

    return(ClstMol,Colordict)

