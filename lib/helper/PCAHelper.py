#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 22:39:38 2019

@author: fujita
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import scipy
import re
from matplotlib.colors import LinearSegmentedColormap
import os


def DelLabem(num_row1,MolLabel):
 
    r = re.compile("(.*)(_)(.*)") 
    MolLabelNew = []    

    for i in range(num_row1):

        d = r.search(MolLabel[i])
        MolLabelNew.append(d.group(1)) 

    return(MolLabelNew)
    
def drawEllipse(XX,pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from scipy.stats import chi2
    from scipy.sparse.linalg import svds

    a = plt.subplot(111, aspect='equal')
    dimention=2
    pca_score2,pca_coef,cov_ratio,pca_components_,W, v = PCA(XX,dimention)
    ScoreThetaMolNamedict = {}
    MolNameList=[]
    EucDistDict = {}
    RadianDict = {}

    for j in range(len(ClstMolDF.columns)):
        pca_score[j,cov_ratioIdx]#
        MolLabel[j]
        Con = list(ClstMolDF.columns)
        #Con.reverse()
        TargetMolList = ClstMolDF[Con[j]]
        TargetMolList = TargetMolList.dropna(how='any')

        PcaScoreAve = np.empty((0,2), int)
        for jj in range(len(TargetMolList)):
            labelloc = np.where(np.array(MolLabel)==TargetMolList[jj])
            PcaScoreForAve = np.append(np.array([pca_score[labelloc[0][0],cov_ratioIdx]]),np.array([pca_score[labelloc[0][0],cov_ratioIdx+1]]), axis=0)
            PcaScoreAve = np.append(PcaScoreAve,[PcaScoreForAve],axis=0)
            ScoreThetaMolNamedict.update({TargetMolList[jj]:PcaScoreForAve})
        #ScoreThetaMolNamedict.update({ScoreThetaMolNameList})
        x,y = np.mean(PcaScoreAve,axis=0) 
        
        PcaScoreAve = colCentering(PcaScoreAve)
        PcaScoreAve = np.array(PcaScoreAve)
        
        U, s, V = np.linalg.svd(PcaScoreAve, full_matrices=True)      

        degree = math.degrees(math.atan2(V.T[1][0],V.T[0][0]))#result1,rad
        Theta = degree
        plt.scatter(x,y, c=list(ClstColorDF[ClstMolDF.columns[j]])[0],s =1,marker='+')

        for k in range(len(TargetMolList)):
            EucDest = np.sqrt(ScoreThetaMolNamedict[TargetMolList[k]][0]**2 + ScoreThetaMolNamedict[TargetMolList[k]][1]**2)
            Radian = np.degrees(math.atan2(ScoreThetaMolNamedict[TargetMolList[k]][1],ScoreThetaMolNamedict[TargetMolList[k]][0]))
            if Radian < 0:
                Radian = Radian + 360
            MolNameList.append(list(ScoreThetaMolNamedict.keys())[k])
            EucDistDict.update({TargetMolList[k]:EucDest})
            RadianDict.update({TargetMolList[k]:Radian})
        xbar = np.linspace(0,1,len(MolNameList))
        x_tick = np.arange(1,len(MolNameList)+1)
        
        if len(s)>1:
            #print([s[0]-np.abs(x),s[1]-np.abs(y)])
            e = patches.Ellipse(xy=(x,y), width=s[0]*2, height=s[1]*2, angle=Theta,fc='none',ec=list(ClstColorDF[ClstMolDF.columns[j]])[0],linewidth=0.5)
            a.add_artist(e)
        else:
            pass
        
def PCA(X,dimention):
    from sklearn.decomposition import PCA
    #pca = PCA(n_components=dimention)
    #pca.fit(X)
    pca = PCA()
    X[np.isnan(X)] = 0
    pca_score = pca.fit_transform(X)
    pca.components_, pca.mean_, pca.get_covariance() 
    cv = pca.get_covariance()
    W, v = np.linalg.eig(cv)

    cov_ratio = pca.explained_variance_ratio_
    #print(pca.components_[0,:])    
    pca_coef= pca.components_.T* np.sqrt(pca.explained_variance_)#FactorLoading
    Xd = pca.transform(X)


    return(pca_score,pca_coef,cov_ratio,pca.components_,W, v)


def mkEachMolColor(pca_score,TimeSignDF,OptionSwitch):
    MolColorDF = OptionSwitch['MolColor']
    
    if len(pca_score[:,0]) > 0:
        #if re.compile("(.*)(_)(.*)").search(jj).group(3) in list(PosDF.columns) :
         #   ii=re.compile("(.*)(_)(.*)").search(jj).group(3)            
        alist = []
        for i in list(TimeSignDF.columns):
            #print(i.get_text())
            nowlabel = re.compile("(.*)(_)(.*)").search(i).group(1)
            alist.append(MolColorDF['MolColor'][nowlabel])
        return(alist)

def mkEachMolClsterColor(pca_score,TimeSignDF,ClstColorDF,ClstMolDF,OptionSwitch):#
    MolColorDF = OptionSwitch['MolColor']
    
    if len(pca_score[:,0]) > 0:
        #if re.compile("(.*)(_)(.*)").search(jj).group(3) in list(PosDF.columns) :
         #   ii=re.compile("(.*)(_)(.*)").search(jj).group(3)            
        alist = []
        for i in list(TimeSignDF.columns):
            #print(i.get_text())
            #nowlabel = re.compile("(.*)(_)(.*)").search(i).group(1)
            #
            try:
                Cluster = ClstMolDF[ClstMolDF.isin([i])].dropna(axis=1,how='all').columns[0]
            except:
                pass
            alist.append(ClstColorDF[Cluster][0])
        return(alist)
        



def PCAScatter(pca_score,label,XIdx,YIdx,ClstMol,MolLabel,BiplotSwitch,OptionSwitch):
    plt.axhline(color="k")
    plt.axvline(color="k")
    Color=BiplotSwitch['DotColor']

    MolColorDF =  OptionSwitch['MolColor']

    hori_align = 'right'
    vert_align = 'top'
    max_xscore = max(pca_score[:,0])*1.2
    max_yscore = max(pca_score[:,1])*1.2
    min_xscore = min(pca_score[:,0])
    min_yscore = min(pca_score[:,1])
    
    fig = plt.figure()
    #colors = [plt.cm.hsv(i/len(Xd)) for i in range(len(Xd))]
    axes = fig.add_subplot(111,aspect='equal')
    for i in range(len(pca_score)):
        if Color=='Black':#
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c='k', label=label[i],s =2)
        elif Color=='EachMol':#
            
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=BiplotSwitch['EachMolColor'][i], label=label[i],s =2)
        elif Color=='EachMolCluster':#           
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=BiplotSwitch['EachMolColor'][i], label=label[i],s =2)            
        elif Color in ['BC','BWoC','CWoB']:          
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=BiplotSwitch['EachMolColor'][i], label=label[i],s =2)            

        elif  Color=='ClstColor_direct':#
                MolColorDF =  OptionSwitch['MolColor']
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=MolColorDF.loc[MolLabel[i]]['ClsterColor'], label=label[i],s =3)
        elif Color=='ClstColor':
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=BiplotSwitch['EachMolColor'].loc[MolLabel[i]]['ClsterColor'], label=label[i],s =3)
        
        else:
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=MolColorDF.loc[MolLabel[i]][Color], label=label[i],s =3)
        if BiplotSwitch['Label']=='Annotate':
            axes.annotate(label[i], xy=(pca_score[i,XIdx],pca_score[i,YIdx]), size=3,horizontalalignment=hori_align, verticalalignment=vert_align)#
    axes.set_xlabel('PC'+str(XIdx+1),fontsize=10)
    axes.set_ylabel('PC'+str(YIdx+1),fontsize=10)
    axes.set_xlim([-max_xscore,max_xscore])
    axes.set_ylim([-max_yscore,max_yscore])
    #axes.set_yticks([-6,-4.5,-3,-1.5,0.0,1.5,3.0,4.5,6])
    #axes.set_yticklabels(['-6.0','-4.5','-3.0','-1.5','0.0','1.5','3.0','4.5','6.0'],size=15)

def plotPCAcoef(pca_coef,LabelProp,pca_score,i,j,MolLabel,BiplotSwitch,OptionSwitch, linecolor='r',fontsize=1):

    Color=BiplotSwitch['DotColor']
    #ScaleFactor= max(pca_score[:,cov_ratioIdx])
    MolColorDF =  OptionSwitch['MolColor']
    ScattColorSwitch=2
    if ScattColorSwitch==0:
        Color='Color'
    else:
        Color='ClsterColor'
    scale_ratio = (np.max(np.sqrt(np.sum(np.power(pca_score[:,:2],2),axis=1, keepdims=True)))/
                 np.max(np.sqrt(np.sum(np.power(pca_coef[:,:2],2),axis=1, keepdims=True))))
    for ii in range(pca_coef.shape[0]):
        #print([0,pca_coef[ii,0]*max_norm],[0,pca_coef[ii,1]*max_norm])
        if ScattColorSwitch==2:
            plt.plot([0,pca_coef[ii,i]*scale_ratio*0.5],[0,pca_coef[ii,j]*scale_ratio*0.5],'brown',linewidth=0.5)
        else:
            plt.plot([0,pca_coef[ii,i]*scale_ratio]*0.7,[0,pca_coef[ii,j]*scale_ratio]*0.7,MolColorDF.loc[MolLabel[ii]][Color],linewidth=0.5)
        if BiplotSwitch['Label']=='Annotate':
            plt.annotate(LabelProp[ii],(pca_coef[ii,i]*scale_ratio*0.5,pca_coef[ii,j]*scale_ratio*0.5), size=3,weight = 'bold')#,     
            #horizontalalignment=hori_align, verticalalignment=vert_align)
        #print(pca_coef.shape[1])
    max_xscore = max(np.abs(pca_coef[:,i]*scale_ratio))*1.1
    max_yscore = max(np.abs(pca_coef[:,j]*scale_ratio))

    


def colCentering(data_2d):
  row_mean = np.mean(data_2d, axis=0,keepdims=True)
  return data_2d - row_mean

def plotPCACovRatio(cov_ratio):

  num_var = len(cov_ratio)
  x_tick = np.arange(1,num_var+1)
  plt.bar(x_tick, cov_ratio)
  plt.plot(x_tick, np.cumsum(cov_ratio),'-o', mfc='none',mec='b',mew=8,linewidth=3)
  plt.xticks(x_tick, fontsize=40)#20
  plt.yticks(fontsize=40)#20

  plt.axis([1-0.4, num_var+0.5, 0,1])
  plt.xlabel("Number of PC", fontsize=40)#10
  plt.ylabel("Explained variance ratio", fontsize=40)#10
  plt.rcParams['axes.linewidth'] = 1.5#
  
def heatPCA(data_2d,label_name): 
  import mpl_toolkits.axes_grid1
  num_x = data_2d.shape[1]
  num_y = data_2d.shape[0]

  c_max = max([data_2d.max(), abs(data_2d.min())])
  colors = np.array([plt.cm.hsv(i/len(data_2d)) for i in range(len(data_2d))])
  ax = plt.imshow(data_2d,aspect='auto',interpolation='nearest',
             cmap='PuOr_r',#color_map, 

             vmin=-c_max,vmax=c_max)
  plt.gca().xaxis.tick_top()

  #cax = divider.append_axes('right', '5%', pad='3%')
  char = plt.colorbar(shrink=1)
  char.ax.tick_params(labelsize=30)
  plt.xticks(range(num_x),np.arange(1,num_x+1).astype('<U'),fontsize=20)
  #plt.ylabel('Time (min.)', fontsize=30)
  plt.yticks(range(num_y),label_name, fontsize=10)
  #plt.title('Number of principal components')
  bx = plt.gca()
  label = bx.set_xlabel('Number of PC', fontsize = 20)
  #bx.xaxis.set_label_coords(0.55,1.15)
  plt.rcParams['axes.linewidth'] = 1.5#
  plt.gca().get_xaxis().set_ticks_position('none')
  plt.gca().get_yaxis().set_ticks_position('none')