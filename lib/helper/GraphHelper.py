#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 11:32:53 2018

@author: fujita
"""

from numpy.random import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.stats import pearsonr, spearmanr
import pandas as pd
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D


def Piledbargraph(cov_ratio,numlevel):

  num_var = len(cov_ratio['1'])
  x_tick = np.arange(1,num_var+1)
  cov_ratiocum=[]
  count=0
  for i in range(numlevel):
      y2=cov_ratio[str(i+1)]
      if count>0:
          plt.bar(x_tick, y2, bottom=y1, color=cm.GnBu_r(count/numlevel))#, width=0.5)
          y1=[ y1[i] +y2[i] for i in range(len(y1))]
      else:
          plt.bar(x_tick, y2, color=cm.GnBu_r(count/numlevel))# width=0.5)
          y1=y2
      for x, y in zip(x_tick, y1):
        plt.text(x, y, count+1, ha='center', va='center',size=1)
      cov_ratiocum+=y2
      count+=1
  plt.xticks(x_tick, fontsize=40)
  plt.yticks(fontsize=40)

  plt.axis([1-0.4, num_var+0.5, 0,1])
  plt.xlabel("Number of PC", fontsize=40)
  plt.ylabel("Explained variance ratio", fontsize=40)
  plt.rcParams['axes.linewidth'] = 1.5   
    

    
def plot3D(DF,OptionDict,save_dir):
    try:
        if OptionDict['plot3DColor'] == 'MolColor':
            ColorList = list(DF['MolColor']);
        if OptionDict['plot3DColor'] == 'ClstColor':
            ColorList = list(MolColorDF['ClstColor'])
        if OptionDict['plot3DColor'] == 'PC12':
            ColorList = OptionDict['LoadingColor']
    except:
        pass
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')#ax = Axes3D(fig,projection='3d')    
    
    ColLabel = list(DF.columns)
    IdxLabel = list(DF.index)

    X = list(DF[ColLabel[1]])
    Y = list(DF[ColLabel[0]])
    Z = list(DF[ColLabel[2]])
    if OptionDict['plot3Dseparation'] == 'Discrete':
        stcount=0;edcount=13;inum=len(X)//13;count=0
        for i in range(inum):
            if count == 83:
                count=0
            ax.plot(X[stcount:edcount], Y[stcount:edcount], Z[stcount:edcount],lw=1,c=cm.jet(count/83))#c=Y,cmap='hsv')

            count+=1
            stcount=edcount;edcount+=13
    elif OptionDict['plot3Dseparation'] == 'Continuous':
        stcount=0;edcount=13;inum=len(X)//13;count=0
        Loading = OptionDict['Loading'] 
        Lmax = np.max(Loading); Lmin=np.min(Loading)
        for i in range(inum):
            if count == 8300000:
                count=0
            ax.plot(X[stcount:edcount], Y[stcount:edcount], Z[stcount:edcount],lw=1,c=cm.jet(((Loading[count]-Lmin)/(Lmax-Lmin))))
            count+=1
            stcount=edcount;edcount+=13
    else:
        stcount=0;edcount=1;inum=len(X)//13;count=0
        Loading = OptionDict['Loading'] 
        Lmax = np.max(Loading); Lmin=np.min(Loading)
        for i in range(len(X)):
            if count == 8300000:
                count=0
            ax.scatter(X[stcount:edcount], Y[stcount:edcount], Z[stcount:edcount],lw=1,c=cm.jet(((Loading[count]-Lmin)/(Lmax-Lmin))))#c=Y,cmap='hsv')
            count+=1
            stcount=edcount;edcount+=1       
    try:    
        if OptionDict['plot3DAnnotation'] == 0:
            for ij in range(len(IdxLabel)):
                ax.text(X[ij], Y[ij], Z[ij],IdxLabel[ij],size=5,zorder=100,color = ColorList[ij] )#, transfoerm=ax.transAxes)     
    except:
        pass

    ax.set_xlabel(ColLabel[1],size=10)
    ax.set_ylabel(ColLabel[0],size=10)
    ax.set_zlabel(ColLabel[2],size=10)
    
    ax.set_xlim(np.min(X), np.max(X))
    ax.set_ylim(np.min(Y), np.max(Y))
    ax.set_zlim(np.min(Z), np.max(Z))
    
    ax.tick_params(labelsize=10)
    plt.savefig(save_dir+'3D.pdf')
    for angle in range(0, 360):
        ax.view_init(30, angle)
        plt.savefig(save_dir + "figs/{0}_{1:03d}.jpg".format('aaall', angle))


    
def colCentering(data_2d):
  row_mean = np.mean(data_2d, axis=0,keepdims=True)
  return data_2d - row_mean

def drawEllipse(ax1,XX,groups, MolLabel,ColLabel1,ColLabel2):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from scipy.stats import chi2
    from scipy.sparse.linalg import svds


    a = ax1
    dimention=2
    ScoreThetaMolNamedict = {}
    MolNameList=[]
    EucDistDict = {}
    RadianDict = {}

    for status, group in groups:
        TargetMolList = list(group.index)
        Color = group['color'][0]

        PcaScoreAve = np.empty((0,2), int)
        for jj in range(len(TargetMolList)):
            PcaScoreForAve = np.append(np.array([XX.loc[TargetMolList[jj],ColLabel1]]),np.array([XX.loc[TargetMolList[jj],ColLabel2]]), axis=0)#2次元の行列(x,y)を縦に足してく
            PcaScoreAve = np.append(PcaScoreAve,[PcaScoreForAve],axis=0)
            ScoreThetaMolNamedict.update({TargetMolList[jj]:PcaScoreForAve})
        x,y = np.mean(PcaScoreAve,axis=0)

        PcaScoreAve = colCentering(PcaScoreAve)
        PcaScoreAve = np.array(PcaScoreAve)
        

        U, S, V = np.linalg.svd(PcaScoreAve, full_matrices=True)
       
        data = np.cov(PcaScoreAve)
        data = np.dot(PcaScoreAve.T,PcaScoreAve)
        W, V_pca = np.linalg.eigh(data)
        index = W.argsort()[::-1]
        V_pca = V_pca[:, index]

        degree = math.degrees(math.atan2(V.T[1][0],V.T[0][0]))

        Theta = degree
        a.plot(x,y, c=Color,markersize =10,marker='+')
        a.set_ylim([-0.1, 1.1])
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
        
        if len(S)>1:
            e = patches.Ellipse(xy=(x,y), width=S[0], height=S[1], angle=Theta,fc='none',ec=Color,linewidth=1)
            a.add_artist(e)
        else:
            pass


####     
def mkScatterWHist(list1,list2,save_dir,ColorList,Optiondict):
    xs, ys = (6,6)
    fig = plt.figure(figsize=(xs,ys))
    
    try:
        if Optiondict['calcR']=='pearson':
            r, p = pearsonr(list1,list2)
        elif Optiondict['calcR']=='spearman':
            r, p = spearmanr(list1,list2)#
    except:
            Optiondict['calcR']='pearson';r, p = pearsonr(list1,list2)
        
    RankTimeVarDF = pd.DataFrame(data=None,columns=[Optiondict['xlabel'],Optiondict['ylabel']])
    RankTimeVarDF[Optiondict['xlabel']] = list1; RankTimeVarDF[Optiondict['ylabel']]= list2

    if (xs==12) and (ys==6):
        ax1 = fig.add_axes((0, 0.25, 0.75, 0.75))
        ax2 = fig.add_axes((0, 1, 0.75, 0.25), sharex=ax1)
        ax3 = fig.add_axes((0.75, 0.25, 0.17, 0.75), sharey=ax1)
    else:
        ax1 = fig.add_axes((0, 0.25, 0.75, 0.75))
        ax2 = fig.add_axes((0, 1, 0.75, 0.25), sharex=ax1)
        ax3 = fig.add_axes((0.75, 0.25, 0.25, 0.75), sharey=ax1)      
    try:
        markersize= Optiondict['markersize'] 
    except:
        markersize= 100 
        
    ax1.tick_params(labelsize=20,direction='out')
    ax3.tick_params(labelleft=False,labelsize=20,axis='y',color='white',left=False,direction='out')#;ax3.tick_params()
    ax3.tick_params(labelleft=False,labelsize=20,axis='x')#;ax3.tick_params()

    try:
        if Optiondict['errorbar']==1:
            x_errlist=Optiondict['x_err']
            y_errlist=Optiondict['y_err']
            for x, y, x_err,y_err,c in zip(list1,list2, x_errlist, y_errlist,ColorList):
                ax1.errorbar(x, y, xerr = x_err, yerr = y_err,  fmt=c, elinewidth=0.5,capsize=1, ecolor='black')
            ax1.scatter(list1,list2,c=ColorList,s=markersize, facecolor='white')# = sns.jointplot(x=DFCol[i]+ str(num[jj]) , y='CV of log10(Param'+DF2Col[ii]+')',color=ColorList,space=0, data=SNSDF)

    except:
        ax1.scatter(list1,list2,c=ColorList,s=markersize,edgecolors='black',linewidths=0.1)# = sns.jointplot(x=DFCol[i]+ str(num[jj]) , y='CV of log10(Param'+DF2Col[ii]+')',color=ColorList,space=0, data=SNSDF) , alpha=0.4)
     
    ax2.hist(list1, bins=20,ec='black')
    ax3.hist(list2, bins=20,orientation="horizontal",ec='black')
    ax2.tick_params(labelleft="true",labelbottom=False,direction='out',labelsize=5,axis='x',color='white');
    ax2.tick_params(labelleft="true",labelsize=20,axis='y');
    plot_axis = plt.axis()


    axis=['top','bottom','left','right']
    line_width=[1,1,1,1]
    
    for a,w in zip(axis, line_width):  # change axis width
        ax1.spines[a].set_linewidth(w)
        ax2.spines[a].set_linewidth(w)        
        ax3.spines[a].set_linewidth(w)
        
    ax1.set_xlabel(Optiondict['xlabel'],fontsize=20)
    ax1.set_ylabel(Optiondict['ylabel'],fontsize=20)  
#### TPSI_TVRI
    #ax1.set_ylim([-0.1,0.8])
    #ax1.set_xlim([75,170])
    #ax1.set_ylim([-0.7,0.7])

    #ax1.set_xticks([ 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]) #[0, 0.2, 0.4, 0.6]
    #ax1.set_xticklabels(['0.25', '0.50', '0.75', '1.00', '1.25', ''])#, rotation=30, fontsize='small') #'0.0', '0.2', '0.4', '0.6'
    #ax1.set_yticks([0, 0.2, 0.4, 0.6])
    #ax1.set_yticklabels(['0.0', '0.2', '0.4', '0.6'])#, rotation=30, fontsize='small')

    xmin, xmax, ymin, ymax = ax1.axis() 
    x = np.arange(xmin, 1, 0.1)

    
    ymin=min(list2); ymax=max(list2)
    xmin=min(list1); xmax=max(list1)
    try:
        if Optiondict['calcR'] in ['pearson','spearman']:
            Roundp = round(p,100)
            Roundr = round(r,3)
            b = '%.2e'%Roundp
            TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
            p='{:e}'.format(p)
            r='{:e}'.format(r)
            TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
            TenRoundr = str(r)[0:4] + '$\it{×10^{-'+str(r)[str(r).find('e')+2:]+'}}$'
            ax1.annotate('$\it{r}$ = ' + TenRoundr + ', $\it{p}$ =' +TenRoundp,fontsize=20, xy=(-730,-10.5))
        else:
            r=0;p=0   
    except:
        r=0;p=0  
    try:
        if Optiondict['mkScatterWHist_drawEllipse']==1:
            NewDF = Optiondict['NewDF'];groups = Optiondict['groups']; MolLabel = Optiondict['MolLabel']
            drawEllipse(ax1,NewDF,groups, MolLabel,'TimeVar','AveCorr')
    except:
        pass
    if Optiondict['Annotate'] == 1:
        AnDict=dict({'Glucose':'Glc','Insulin':'Ins','C-peptide':'CRP','GIP(Active)':'GIP','Pyruvate':'Pyr','Total bile acid':'TBA',
                   'Citrate':'Cit','Cortisol':'Cor','Free fatty acid':'FFA','Total ketone body':'Ketone','Glutamic acid':'Glu',
                   'Citrulline':'Citr','Methionine':'Met','Isoleucine':'Ile','Leucine':'Leu','Tyrosine':'Tyr','4-Methyl-2-oxopentanoate':'4M2O','Glu + threo-beta-methylaspartate':'Glu+TBM','Growth hormone':'GH'})
        for i in range(len(list1)):
            if Optiondict['Label'][i] in ['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
                   'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
                   'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate','Growth hormone']:
                ax1.annotate(AnDict[Optiondict['Label'][i]],fontsize=0.1, xy=(list1[i],list2[i]))
            else:
                try:
                    ax1.annotate(Optiondict['Label'][i],fontsize=Optiondict['LabelSize'], xy=(list1[i],list2[i]))     
                except:
                    ax1.annotate(Optiondict['Label'][i],fontsize=7, xy=(list1[i],list2[i]))
        RankTimeVarDF.index=Optiondict['Label']
    try:        
        if Optiondict['y=x'] == 1:#
                x = np.linspace(xmin,xmax)
                y = x
                ax1.plot(x,y,"r-",linewidth=1)
    except:
        pass

    title = Optiondict['title']+Optiondict['calcR']

    RankTimeVarDF.to_excel(save_dir + 'DF_'+title+'.xlsx')

    plt.savefig(save_dir +'Scatter_'+title+'.pdf',format='pdf',bbox_inches="tight")
    ##plt.savefig(save_dir +'Scatter_'+title+'.png',format='png',bbox_inches="tight")

    return(r,p)

    
def mkSortedBar(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir):
    
    
    x = np.linspace( 1, len(List1), len(List1) )
    fig = plt.figure(figsize=(16,16))

    plt.bar( x, List1,width=1, color=Color,tick_label=xticks,linewidth=1,ec='black')   
    
    plt.title(Title,fontsize=Titlesize)
    plt.xlabel(xlabel,fontsize=xsize)
    plt.ylabel(ylabel,fontsize=size)
    #plt.xticks(x_tick)
    plt.xticks(rotation=270,fontsize=xsize)
    plt.yticks(fontsize=size)
    plt.savefig(save_dir+Title+'_Sort.pdf',bbox_inches="tight")
    plt.close()

def mkBarAveStd(List1, StdList, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize,save_dir):
    
    
    x = np.linspace( 1, len(List1), len(List1) )
    fig = plt.figure(figsize=(20,12))
#error_kw=dict(lw=5, capsize=5, capthick=3)
    plt.bar( x, List1,width=1,yerr = StdList, color=Color,tick_label=xticks,linewidth=0.5,ec='black',error_kw=dict(lw=1, capsize=1, capthick=1))
    
    
    plt.title(Title,fontsize=Titlesize)
    plt.xlabel(xlabel,fontsize=xsize)
    plt.ylabel(ylabel,fontsize=size)
    #plt.xticks(x_tick)
    plt.xticks(rotation=270,fontsize=xsize) 
    plt.yticks(fontsize=size)    
    xmin, xmax, ymin, ymax = plt.axis() 
    #plt.plot([xmin, xmax],[1, 1], "black", linestyle='dashed') # normal way
    plt.savefig(save_dir+Title+'_.pdf',bbox_inches="tight")
    plt.close()
    
def mkBar(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize,save_dir):
    
    
    x = np.linspace( 1, len(List1), len(List1) )
    
    fig = plt.figure(figsize=(5,3))
    ax1 = fig.add_axes((0, 1, 1, 1))

    ax1.bar( [0+0.1*i for i in range(len(xticks))], List1,width=0.05, color=Color,tick_label=xticks,linewidth=0.5,ec='black')
    
    ax1.set_ylim([-1.0,1.0])
    
    xmin, xmax, ymin, ymax = ax1.axis() 
    
    p = ax1.hlines([0], xmin, xmax, "black", linestyles='solid',linewidth=1)     # hlines
    
    axis=['top','bottom','left','right']
    line_width=[1,1,1,1]
    
    for a,w in zip(axis, line_width):  # change axis width
        ax1.spines[a].set_linewidth(w)

    ax1.set_xlim(xmin, xmax)
    ax1.set_title(Title,fontsize=Titlesize)
    ax1.set_xlabel(xlabel,fontsize=xsize)
    ax1.set_ylabel(ylabel,fontsize=size)
    #plt.xticks(x_tick)
    ax1.set_xticklabels(labels=xticks,rotation=270,fontsize=xsize) 
    #ax1.set_xticklabels(fontsize=size)    
    xmin, xmax, ymin, ymax = ax1.axis() 
    #plt.plot([xmin, xmax],[1, 1], "black", linestyle='dashed') # normal way
    plt.savefig(save_dir+Title+'_.pdf',bbox_inches="tight")
    plt.close()


