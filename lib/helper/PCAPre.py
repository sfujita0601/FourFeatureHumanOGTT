#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 20:54:49 2019
PCA
@author: fujita
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import matplotlib as mpl
from matplotlib.colors import Normalize
from . import GraphHelper as GH
import scipy
import os
import re
from . import PCAHelper as PCAHel
import matplotlib.cm as cm
 
   

def LoadingPCsimulation(timepoint,pca_coef,pca_score,currPC,othrPC,cov_ratio,ylim,MolLabel,label,BiplotSwitch,ClstColorDF,save_dir):

    Name = ''
    plt.rcParams['axes.linewidth'] = 1.0

    if BiplotSwitch['LoadingPCsimulation'] == 'PC1_2_represen':
        fig, host2 = plt.subplots(2,2,figsize=(4, 3))
    else:    
        fig, host2 = plt.subplots(2,5,figsize=(18, 18))
    currx=0;curry=0
    totalfignum=9
    maxscore=np.max(pca_score[:,currPC]);minscore=np.min(pca_score[:,currPC])
    maxscore2=np.max(pca_score[:,othrPC]);minscore2=np.min(pca_score[:,othrPC])
    if currPC == 0:
        maxscore=10;minscore=-10
        maxscore2=8;minscore2=-8
    else:
        maxscore=8;minscore=-8
        maxscore2=10;minscore2=-10
    NewDF= pd.DataFrame(data=None,index=[],columns=['Time','PC'+str(currPC+1),'Loading'])#'Loading'])'PC'+str(currPC+2)
    Col= list(NewDF.columns)
    Time=[];PC1=[];Loading=[] ;ColorList=[]    ;Loading2=[]   
    count=0
    try:
        if BiplotSwitch['LoadingPCsimulation'] == 'PC1_2':
            tryList1 = [-7.5,-2.5,2.5,7.5]
            tryList2=[4.5]
        elif BiplotSwitch['LoadingPCsimulation'] == 'PC1_2_represen': #Choose one molecule
            Name = 'Glucose'#'Growth hormone','Methionine','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
                   #'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
                   #'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate',
               #'Glu + threo-beta-methylaspartate','Growth hormone']### 'SM-C IGF-1_B', 'SM-C IGF-1_C', 'Tyrosine_C', 'Tyrosine_B'
            IndList=[list(label).index(Name)]
            
            tryList1 = [pca_score[IndList[jj],currPC] for jj in range(len(IndList))]
            tryList2=[pca_score[IndList[jj],othrPC] for jj in range(len(IndList))]
        else:
            tryList1 = np.linspace(minscore,maxscore,9)
            tryList2=np.linspace(minscore2,maxscore2,83)
    except:
        pass
    
    for ii in tryList1:
        for i in tryList2:
            count+=1
            if BiplotSwitch['LoadingPCsimulation'] == 'ClstColor':
                plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['ColorDF'].loc[MolLabel[ii]]['ClsterColor'],linewidth=0.5, marker = 'o',markersize=3)               
            elif BiplotSwitch['LoadingPCsimulation'] == 'EachMolCluster':
                plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['EachMolColor'][ii],linewidth=0.5, marker = 'o',markersize=3)
            elif BiplotSwitch['LoadingPCsimulation'] == 'plot3D':
                plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['EachMolColor'][ii],linewidth=0.5, marker = 'o',markersize=3)
                Time += list(timepoint);PC1 += [ii]*len(timepoint); Loading += [i]*len(timepoint)#
                Loading2 += list(ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC]);ColorList+=list(cm.hsv(count/83))*len(timepoint)

            elif (BiplotSwitch['LoadingPCsimulation'] == 'PC1_2') or (BiplotSwitch['LoadingPCsimulation'] == 'PC1_2_represen'):

                if (BiplotSwitch['LoadingPCsimulation'] == 'PC1_2_represen') and (count==1):
                    flag=1
                elif (BiplotSwitch['LoadingPCsimulation'] == 'PC1_2'):
                    flag=1
                else:
                    flag=0
                if flag==1:
                    #PC1+PC2
                    host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC], color='black',linewidth=0.5, marker = 'o',markersize=3)
                    #PC1
                    host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC], color='red', linestyle = "--",linewidth=0.5, marker='+',markersize=3)
                    #PC2
                    host2[curry,currx].plot(timepoint, i*pca_coef[:,othrPC], color='blue', linestyle = "--",linewidth=0.5,marker='+',markersize=3)
                    host2[curry,currx].tick_params(bottom=False, left=False,  right=False,     top=False)
                    host2[curry,currx].set_ylim([-5.5,5.5])
                    host2[curry,currx].tick_params(axis='both', colors='white')
            else: #¥
                if i==0:
                    host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC], color='black',linewidth=0.5, marker = 'o',markersize=3)
                else:
                    
                    host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC], color=cm.jet(count/83),linewidth=0.5, marker = 'o',markersize=3)
        try:
            if (ii==minscore):
                NewPanel = pd.Panel({ii:NewDF})
            else:
                NewPanel[ii] = NewDF                 
        except:
            pass
        
        if (BiplotSwitch['LoadingPCsimulation'] != 'plot3D') and (BiplotSwitch['LoadingPCsimulation'] != 'PC1_2') and (BiplotSwitch['LoadingPCsimulation'] != 'PC1_2_represen'):
    
            host2[curry,currx].set_xticks([0,60,120,180,240])
            host2[curry,currx].set_xticklabels(['0','60','120','180','240'],size=15)
            host2[curry,currx].set_xlabel('Time (min.)',size=20)
            host2[curry,currx].set_title('PC'+str(currPC+1)+'='+str(round(ii,3)),size=15)
            host2[curry,currx].set_ylabel('Loading',size=20)
            host2[curry,currx].tick_params(axis='y', labelsize=15)
            host2[curry,currx].set_ylim([-5.5,5.5])
            host2[curry,currx].set_ylim([-7,7])
        count=0#
        currx +=1;
        if currx==5:
            currx=0
            curry+=1

    
    NewDF[Col[0]] = Time; NewDF[Col[1]] = PC1; NewDF[Col[2]] = Loading2
   
    if BiplotSwitch['LoadingPCsimulation'] == 'plot3D':
        BiplotSwitch['plot3DAnnotation']=1
        BiplotSwitch['plot3DColor']= 'PC12'
        BiplotSwitch['LoadingColor']=ColorList
        BiplotSwitch['plot3Dseparation'] ='Discrete'#'Discrete'　：'Continuous'
        BiplotSwitch['Loading'] = Loading2 
        OptionDict=BiplotSwitch
        GH.plot3D(NewDF,OptionDict,save_dir+str(ii)+'_')
    if (BiplotSwitch['LoadingPCsimulation'] != 'plot3D') and (BiplotSwitch['LoadingPCsimulation'] != 'PC1_2') and (BiplotSwitch['LoadingPCsimulation'] != 'PC1_2_represen'):

        host2[1,4].tick_params(labelbottom="off",bottom="off") 
        host2[1,4].tick_params(labelleft="off",left="off")
        host2[1,4].tick_params(bottom=False, left=False,  right=False,     top=False)
        host2[1,4].spines["right"].set_color("none")
        host2[1,4].spines["left"].set_color("none")
        host2[1,4].spines["top"].set_color("none")
        host2[1,4].spines["bottom"].set_color("none") 
        t=np.linspace(0,1,pca_score.shape[0])
        x = np.linspace(0,0,pca_score.shape[0])
        yy=np.linspace(0,0,pca_score.shape[0])
        ax = plt.scatter(x,yy,c=t,cmap=cm.jet,marker='.',lw=0)

        CF = plt.colorbar(ax)
        CF.ax.tick_params(labelsize=10)
        CF.ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False, color='white')
        for j, lab in enumerate(np.linspace(minscore2,maxscore2,9)): 
            CF.ax.text(2.5, (j )/9+0.05, round(lab,2), ha='center', va='center', size=20) 
  

    plt.subplots_adjust(wspace=0.6, hspace=0.6)
    plt.savefig(save_dir + 'LoadingPlotWCluster'+str(currPC)+'_'+str(othrPC)+Name+'.pdf') 
    plt.savefig(save_dir + 'LoadingPlotWCluster'+str(currPC)+'_'+str(othrPC)+Name+'.png') 
    currplot =0 
    plt.close()
                
def FactorLoadingPlotWCluster(timepoint,pca_coef,pca_score,currPC,cov_ratio,ylim,MolLabel,label,BiplotSwitch,ClstColorDF):
    import matplotlib.colors as clr
    import matplotlib as mpl
    MolColorDF = BiplotSwitch['MolColor']

    ClstColorList=list(ClstColorDF.loc[0])
    cmap = mpl.colors.ListedColormap(ClstColorList)

    for ii in range(pca_score.shape[0]):
        if BiplotSwitch['DotColor'] == 'ClstColor':
            plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['ColorDF'].loc[MolLabel[ii]]['ClsterColor'],linewidth=0.5, marker = 'o',markersize=3)

        elif BiplotSwitch['DotColor'] == 'EachMol':
            Color=BiplotSwitch['DotColor']

            plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['EachMolColor'][ii],linewidth=0.5, marker = 'o',markersize=3)
            plotdata = pca_coef*pca_score[ii]
            if any(plotdata > 5):
                aa= np.where(plotdata > 5)[0]
                for jj in aa:
                
                    plt.text(timepoint[jj], plotdata[jj], label[ii],size=1)
            elif any(plotdata < -5):
                aa= np.where(plotdata < -5)[0]
                for jj in aa:
                
                    plt.text(timepoint[jj], plotdata[jj], label[ii],size=1)
            else:
                    plt.text(timepoint[jj], plotdata[jj], label[ii],size=1)
                
        elif BiplotSwitch['DotColor'] == 'EachMolCluster':
            plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['EachMolColor'][ii],linewidth=0.5, marker = 'o',markersize=3)
        elif BiplotSwitch['DotColor'] == 'MolColor':
            if label[ii] in ['Citrate', 'Growth hormone']:
                plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=MolColorDF['MolColor'].iloc[ii],linewidth=0.5, marker = 'o',markersize=3)
        else:
            plt.plot(timepoint, pca_coef*pca_score[ii]  ,color='k',linewidth=0.5, marker = 'o',markersize=3)
    plt.xticks([0,60,120,180,240],['0','60','120','180','240'],size=15)

    t=np.linspace(0,1,pca_score.shape[0])
    x = np.linspace(0,0,pca_score.shape[0])
    yy=np.linspace(0,0,pca_score.shape[0])
    ax = plt.scatter(x,yy,c=t,cmap=cmap,marker='.',lw=0)
    CF = plt.colorbar(ax)

    CF.ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False, color='white')
    for j, lab in enumerate(['1','2','3','4','5','6','7','8','9','10','11','12','13']): 
        CF.ax.text(1.9, ( j )/13 + 0.03 , lab, ha='center', va='center', size=20) 
    plt.title('PC'+str(currPC+1)+' ('+str(cov_ratio*100)[0:4]+'%)', size=20)

    plt.ylim([-ylim,ylim])
    plt.xlabel('Time(min)',size=20)
def mkZscore(CorrFastPropDF,label,labelProp):
    meanNP = np.array(CorrFastPropDF) -np.nanmean(np.array(CorrFastPropDF),axis=0,keepdims=True)
    stdNP = np.nanstd(np.array(CorrFastPropDF).astype(float),axis=0)
    signDF = pd.DataFrame(index=label,columns=labelProp)
    for i in range(len(labelProp)):
        signNP = meanNP[:,i]/stdNP[i]
        signDF[labelProp[i]] = signNP
    return(signDF)
    
def PCAprep(XDF,ClstAveTimeSeriesDF,ClstColorDF,MolLabel,Label,ClstMolDF,AnalSwitch,LoadingOpt,BiplotSwitch,OptionSwitch,CompareBolus,ClstNoColorDF,save_dir):
    
    NormlzDF = mkZscore(XDF,MolLabel,Label)
    X = np.array(NormlzDF)
    XDF=XDF.T
    pca_score,pca_coef,cov_ratio,pca_components_,W,v = PCAHel.PCA(X,2)
    try:    
        if CompareBolus['Target'] == 'Betw':#
            pca_score[:,0] = pca_score[:,0]*-1
            pca_coef[:,0] = pca_coef[:,0]*-1#
            pca_components_[:,0] = pca_components_[:,0]*-1
    except:
        pass
    label=MolLabel

    if BiplotSwitch['DotColor'] == 'EachMol':
        BiplotSwitch['EachMolColor'] = PCAHel.mkEachMolColor(pca_score,XDF,OptionSwitch)
    if BiplotSwitch['DotColor'] == 'EachMolCluster':
        BiplotSwitch['EachMolColor'] = PCAHel.mkEachMolClsterColor(pca_score,XDF,ClstColorDF,ClstMolDF,OptionSwitch)
    if BiplotSwitch['DotColor'] == 'ClstColor':
        BiplotSwitch['EachMolColor'] = ClstNoColorDF   
    if AnalSwitch['LengVector'] == 1:#
        CompareBolus['lengthPC12'] = PCAHel.calclenbetwScore(pca_score,BiplotSwitch,CompareBolus,MolLabel,save_dir)   
        PCAHel.calclenScore(pca_score,BiplotSwitch,CompareBolus,MolLabel,save_dir)
    else:
        CompareBolus['lengthPC12']=20
    if AnalSwitch['PlotPCACovRatio'] == 1:
        fig = plt.figure(figsize=(10,8))
        PCAHel.plotPCACovRatio(cov_ratio)
        plt.xticks(size='40');plt.yticks(size='40')
        plt.savefig(save_dir + 'PCExplained.pdf',bbox_inches="tight")
        ##plt.savefig(save_dir + 'PCExplained.png',bbox_inches="tight")
    if AnalSwitch['Biplot'] == 1:
        cov_ratioIdx = 0
        for i in range(min(np.where(np.cumsum(cov_ratio)>0.85)[0])):
            for j in range(i+1,min(np.where(np.cumsum(cov_ratio)>0.85)[0])+1):
                fig = plt.figure()
                fig.set_size_inches(pca_score.shape[1], pca_score.shape[0]/5.0)
                PCAHel.PCAScatter(pca_score,label,i,j,ClstMolDF,label,BiplotSwitch,OptionSwitch)
                fig.tight_layout() 
                plt.xlabel('PC'+str(i+1)+' ('+str(cov_ratio[i]*100)[0:4]+'%)')
                plt.ylabel('PC'+str(j+1)+' ('+str(cov_ratio[j]*100)[0:4]+'%)')
                plt.xticks([-10,-7.5,-5.0,-2.5,0,2.5,5,7.5,10],['-10','-7.5','-5.0','-2.5','0','2.5','5.0','7.5','10'],size=15)
                plt.yticks([-4,-2,0,2,4],['-4.0','-2.0','0','2.0','4.0'],size=15)

                plt.xticks(size='10');plt.yticks(size='10')
                plt.ylim([-5,5])#
                plt.xlim([-11,11])#
                axis=['top','bottom','left','right']
                line_width=[1,1,1,1]
                
                for a,w in zip(axis, line_width):  # change axis width
                    plt.gca().spines[a].set_linewidth(w)
                plt.savefig(save_dir + 'PCAScatter'+ str(i+1) + 'vs' + str(j+1) + '.pdf')
                if AnalSwitch['DrawEllipse']  == 1:#
                    XX = ClstAveTimeSeriesDF
                    PCAHel.drawEllipse(XX,pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF)
                    plt.ylim([-5.5,5.5])#
                    plt.xlabel('PC'+str(i+1)+' ('+str(cov_ratio[i]*100)[0:4]+'%)')
                    plt.ylabel('PC'+str(j+1)+' ('+str(cov_ratio[j]*100)[0:4]+'%)')
                    plt.xticks(size='10');plt.yticks(size='10')
                    if AnalSwitch['WHist']==1:
                        Optiondict={'calcR':'',#'spearman'
                        'xlabel':'',
                            'ylabel':'',
                            'Annotate':0,
                            'Label':[],
                            'title':'',
                            'y=x' :0 
                            }
                        GH.mkScatterWHist(list(pca_score[:,0]),list(pca_score[:,1]),save_dir,['k']*len(pca_score[:,0]),Optiondict)
                        plt.ylim([-5.5,5.5])
                    axis=['top','bottom','left','right']
                    line_width=[1,1,1,1]
                    for a,w in zip(axis, line_width): 
                        plt.gca().spines[a].set_linewidth(w)
                    plt.savefig(save_dir + 'PCAScatterWEllipse'+ '_PC'+str(i+1)+ 'vs' + 'PC'+str(j+1) + '.pdf',bbox_inches="tight")

        cov_ratioIdx += 1
        plt.close()     

    
    if len(BiplotSwitch['LoadingPCsimulation']) >0: #
        currplot =0  
        max_score=max(np.abs(pca_score[:,0]))*max(np.abs(pca_components_.T[:,0]))*1.2

        LoadingPCsimulation(list(XDF.index),pca_components_.T[:,0:2],pca_score,0,1,cov_ratio,max_score,MolLabel,label,BiplotSwitch,ClstColorDF,save_dir)
        LoadingPCsimulation(list(XDF.index),pca_components_.T[:,0:2],pca_score,1,0,cov_ratio,max_score,MolLabel,label,BiplotSwitch,ClstColorDF,save_dir)
     
    if AnalSwitch['LoadingPlotWCluster'] ==1:
        BiplotSwitch['MolColor'] = OptionSwitch['MolColor']
        pca_components_#
        max_score=3.5
        subplot_row = 3
        subplot_col = 3
        fig, host = plt.subplots(subplot_row,subplot_col,figsize=(18, 18))
        totalfignum = pca_score.shape[1]
        currplot =0  
        
        if OptionSwitch['Target'] == 'EachCluster':
            for ii in range(len(ClstMolDF.columns)
            ):
                locIdx=[XDF.columns.get_loc(i) for i in list(ClstMolDF['cluster_'+str(ii+1)].dropna())]
                fig, host = plt.subplots(subplot_row,subplot_col,figsize=(18, 18))

                currplot =0   
                for i in range(subplot_row):
                    for j in range(subplot_col):                
                        plt.subplot(subplot_row, subplot_col, 1+j+i*subplot_col)
                        MolLabel=list(ClstMolDF['cluster_'+str(ii+1)].dropna())
                        FactorLoadingPlotWCluster(list(XDF.index),pca_components_[:,currplot],pca_score[locIdx,currplot],currplot,cov_ratio[currplot],max_score,MolLabel,label,BiplotSwitch,ClstColorDF)
                        plt.yticks(size='20')
                        [host[i][j].spines[ii].set_linewidth(0.1) for ii in host[i][j].spines.keys()]
                        currplot +=1
                        if currplot >= totalfignum:
                            break
                    if currplot >= totalfignum:
                        break
                plt.subplots_adjust(wspace=0.6, hspace=0.6)
                plt.savefig(save_dir + 'LoadingPlotWCluster'+str(ii)+'.pdf') 
                ##plt.savefig(save_dir + 'LoadingPlotWCluster'+str(ii)+'.png') 
                currplot =0 
                plt.close()               
        else:
            for i in range(subplot_row):
                for j in range(subplot_col):                
                    plt.subplot(subplot_row, subplot_col, 1+j+i*subplot_col)
                    FactorLoadingPlotWCluster(list(XDF.index),pca_components_.T[:,currplot],pca_score[:,currplot],currplot,cov_ratio[currplot],max_score,MolLabel,label,BiplotSwitch,ClstColorDF)
                    plt.yticks(size='20')
                    fig.tight_layout() 
                    currplot +=1
                    if currplot >= totalfignum:
                        break
                if currplot >= totalfignum:
                    break

                host[i][j].spines['bottom'].set_linewidth(0.1); host[i][j].spines['left'].set_linewidth(0.1); host[i][j].spines['right'].set_linewidth(0.1); host[i][j].spines['top'].set_linewidth(0.1); 

            plt.subplots_adjust(wspace=0.4, hspace=0.4)
            
            plt.savefig(save_dir + 'LoadingPlotWCluster.pdf') 
            ##plt.savefig(save_dir + 'LoadingPlotWCluster.png') 
            plt.close()
            
    if AnalSwitch['ScoreHeatMap'] ==1:
        fig, ax = plt.subplots()
        fig.set_size_inches(pca_score.shape[1], pca_score.shape[0]/5.0)
        ScoreDF = pd.DataFrame(pca_score)
        ScoreDF.to_excel(save_dir+'Score.xlsx')
        PCAHel.heatPCA(pca_score,label)   
        MolColorDF =OptionSwitch['MolColor']
        if len(pca_score[:,0]) > len(MolColorDF['MolColor']):#        
            splabel = ax.yaxis.get_ticklabels()
            for i in splabel:
                nowlabel = re.compile("(.*)(_)(.*)").search(i.get_text()).group(1)
                i.set_color(MolColorDF['MolColor'][nowlabel])
        else:
            [t.set_color(i) for (i,t) in
             zip(list(MolColorDF['MolColor']),ax.yaxis.get_ticklabels())]
        plt.savefig(save_dir + 'PCAScoreHeatmap.pdf',bbox_inches="tight")
         ##plt.savefig(save_dir + 'PCAScoreHeatmap.png',bbox_inches="tight")

        plt.close()
        
    if AnalSwitch['FactorLoadingHeatMap']==1:
        ## FactorLoading Heatmap
        fig = plt.figure(figsize=(5,10))
        pca_coef[0,-1] = 0
        if LoadingOpt == 0:
            fig.set_size_inches(pca_coef.shape[1]*1.2, pca_coef.shape[0]/2)
            PCAHel.heatPCA(pca_coef, Label)

        fig.tight_layout() 
        plt.xticks(size='40');plt.yticks(size='40')
        plt.ylabel('Time(min.)',size=40)
        plt.savefig(save_dir + 'PCAFactorLoadingHeatmap.pdf',bbox_inches="tight")
        ##plt.savefig(save_dir + 'PCAFactorLoadingHeatmap.png',bbox_inches="tight")
        plt.close()

        