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
from matplotlib.colors import LinearSegmentedColormap
import re
from . import PCAHelper as PCAHel
import matplotlib.cm as cm
 
   

def LoadingPCsimulation(timepoint,pca_coef,pca_score,currPC,othrPC,cov_ratio,ylim,MolLabel,label,BiplotSwitch,ClstColorDF,save_dir):#あるPCのスコア入れたらある数降って、別のPCの分だけ別に描画
    ClstColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191219/RawDelta_NormalizationByStd/ClstColor.xlsx',header=0,index_col=0)

    ClstColorList=list(ClstColorDF.loc[0])
    ClstColorList.reverse()
    cmap = mpl.colors.ListedColormap(ClstColorList)
    Name = ''
    plt.rcParams['axes.linewidth'] = 1.0

    if BiplotSwitch['DotColor'] == 'PC1_2_represen':
        fig, host2 = plt.subplots(2,2,figsize=(4, 3))#2,2,figsize=(4, 3)
    else:    
        fig, host2 = plt.subplots(2,5,figsize=(18, 18))#2,2,figsize=(4, 3)
    currx=0;curry=0
    totalfignum=9
    maxscore=np.max(pca_score[:,currPC]);minscore=np.min(pca_score[:,currPC])
    maxscore2=np.max(pca_score[:,othrPC]);minscore2=np.min(pca_score[:,othrPC])
    if currPC == 0:#PA1固定するなら
        maxscore=10;minscore=-10
        maxscore2=8;minscore2=-8
        #maxscore2=10;minscore2=-10
    else:
        maxscore=8;minscore=-8
        maxscore2=10;minscore2=-10
    NewDF= pd.DataFrame(data=None,index=[],columns=['Time','PC'+str(currPC+1),'Loading'])#'Loading'])'PC'+str(currPC+2)
    Col= list(NewDF.columns)
    Time=[];PC1=[];Loading=[] ;ColorList=[]    ;Loading2=[]   
    count=0# np.linspace(-10,10,9):#
    try:
        if BiplotSwitch['DotColor'] == 'PC1_2':
            tryList1 = [-7.5,-2.5,2.5,7.5]
            tryList2=[4.5]#ここは手で変える
        elif BiplotSwitch['DotColor'] == 'PC1_2_represen':###　　ここも手で変える
            Name = 'Growth hormone'#'Methionine'#,'Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
                   #'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
                   #'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate',
               #'Glu + threo-beta-methylaspartate','Growth hormone']### 'SM-C IGF-1_B', 'SM-C IGF-1_C', 'Tyrosine_C', 'Tyrosine_B'
            IndList=[list(label).index(Name)]#+[list(label).index('Glucose_C')]
            
            #[list(label).index('Growth hormone_B')]+[list(label).index('Growth hormone_C')]+
            #[list(label).index('Growth hormone')]+[list(label).index('Free fatty acid')]+[list(label).index('Total ketone body')]+[list(label).index('Isoleucine')]#[list(label).index('Glucose')]+[list(label).index('Insulin')]+[list(label).index('GIP(Active)')]+[list(label).index('C-peptide')]#+[list(label).index('Growth hormone')]+[list(label).index('Free fatty acid')]+[list(label).index('Total ketone body')]+[list(label).index('Isoleucine')]+[list(label).index('Citrulline')]#GIP, Glc, Ins, CRP, GH, Ketone, FFA, Ile, Citのインデックス入手         
            tryList1 = [pca_score[IndList[jj],currPC] for jj in range(len(IndList))]
            tryList2=[pca_score[IndList[jj],othrPC] for jj in range(len(IndList))]#ここは手で変える
        else:
            tryList1 = np.linspace(minscore,maxscore,9)
            tryList2=np.linspace(minscore2,maxscore2,83)
    except:
        pass
    
    for ii in tryList1:#np.linspace(minscore,maxscore,9):#PC1を固定で、PC2の成分を取りうる範囲内で増減させる。
        for i in tryList2:#np.linspace(minscore2,maxscore2,83):#値を振る
            count+=1
            if BiplotSwitch['DotColor'] == 'ClstColor':
                plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['ColorDF'].loc[MolLabel[ii]]['ClsterColor'],linewidth=0.5, marker = 'o',markersize=3)               
            elif BiplotSwitch['DotColor'] == 'EachMolCluster':
                plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['EachMolColor'][ii],linewidth=0.5, marker = 'o',markersize=3)
            elif BiplotSwitch['DotColor'] == 'plot3D':
                #plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['EachMolColor'][ii],linewidth=0.5, marker = 'o',markersize=3)
                #この時点であるPC1, いろんなPC2のLpadingを3Dにplotする
                Time += list(timepoint);PC1 += [ii]*len(timepoint); Loading += [i]*len(timepoint)#
                Loading2 += list(ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC]);ColorList+=list(cm.hsv(count/83))*len(timepoint)

            elif (BiplotSwitch['DotColor'] == 'PC1_2') or (BiplotSwitch['DotColor'] == 'PC1_2_represen'):

                if (BiplotSwitch['DotColor'] == 'PC1_2_represen') and (count==1):
                    flag=1
                elif (BiplotSwitch['DotColor'] == 'PC1_2'):
                    flag=1
                else:
                    flag=0
                if flag==1:
                    #足し合わせは黒で
                    host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC], color='black',linewidth=0.5, marker = 'o',markersize=3)
                    #PC1だけは赤点線
                    host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC], color='red', linestyle = "--",linewidth=0.5, marker='+',markersize=3)
                    #PC1だけは青点線
                    host2[curry,currx].plot(timepoint, i*pca_coef[:,othrPC], color='blue', linestyle = "--",linewidth=0.5,marker='+',markersize=3)
                    host2[curry,currx].tick_params(bottom=False, left=False,  right=False,     top=False)
                    host2[curry,currx].set_ylim([-5.5,5.5])
    
                    host2[curry,currx].tick_params(axis='both', colors='white')#y : 15

            else: #PC片方を離散的に変えた時の、もう片方のPCを何通りかにフルloading波形
                if i==0:
                    host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC], color='black',linewidth=0.5, marker = 'o',markersize=3)
                #elif i==ii:
                 #   host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC], color='white',linewidth=0.5, marker = 'o',markersize=3)
                    
                else:
                    #host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC], color='black',linewidth=0.5, marker = 'o',markersize=3)
                    
                    host2[curry,currx].plot(timepoint, ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC], color=cm.jet(count/83),linewidth=0.5, marker = 'o',markersize=3)

        #NewDF['Time'] = Time; NewDF['PC1'] = PC1; NewDF['Loading'] = Loading
        try:
            if (ii==minscore):
                NewPanel = pd.Panel({ii:NewDF})
            else:
                NewPanel[ii] = NewDF#PC片方変化させた時のもう片方のPCとLading                   
        except:
            pass
        
        if (BiplotSwitch['DotColor'] != 'plot3D') and (BiplotSwitch['DotColor'] != 'PC1_2') and (BiplotSwitch['DotColor'] != 'PC1_2_represen'):#forで回して全部描画する
    
            host2[curry,currx].set_xticks([0,60,120,180,240])
            host2[curry,currx].set_xticklabels(['0','60','120','180','240'],size=15)
            host2[curry,currx].set_xlabel('Time (min.)',size=20)
            host2[curry,currx].set_title('PC'+str(currPC+1)+'='+str(round(ii,3)),size=15)
            host2[curry,currx].set_ylabel('Loading',size=20)
            host2[curry,currx].tick_params(axis='y', labelsize=15)#y : 15
            
            host2[curry,currx].set_ylim([-5.5,5.5])
            host2[curry,currx].set_ylim([-7,7])

            #host2[curry,currx].tick_params(bottom=False, left=False,  right=False,     top=False)

        count=0#
        currx +=1;
        if currx==5:
            currx=0
            curry+=1

    #PC1,Loading = np.meshgrid(np.arange(minscore,maxscore,0.5), np.arange(minscore2,maxscore2,0.5))
    #Loading2 = list(ii * pca_coef[:,currPC]+i*pca_coef[:,othrPC])
    
    NewDF[Col[0]] = Time; NewDF[Col[1]] = PC1; NewDF[Col[2]] = Loading2
    #とりあえず隠列はtimepointの数ずつ区切られている。
    
    if BiplotSwitch['DotColor'] == 'plot3D':#forで回して全部描画する

        BiplotSwitch['plot3DAnnotation']=1
        BiplotSwitch['plot3DColor']= 'PC12'
        BiplotSwitch['LoadingColor']=ColorList
        BiplotSwitch['plot3Dseparation'] ='Discrete'#離散的に'Discrete'　：'Continuous':#連続的に
        BiplotSwitch['Loading'] = Loading2 
        OptionDict=BiplotSwitch
            #for ii in np.linspace(minscore,maxscore,9):
         #   GH.plot3D(NewPanel[ii],OptionDict,save_dir+str(ii)+'_')
        GH.plot3D(NewDF,OptionDict,save_dir+str(ii)+'_')
    #plt.ylabel('Valiation Index',size=20)
    if (BiplotSwitch['DotColor'] != 'plot3D') and (BiplotSwitch['DotColor'] != 'PC1_2') and (BiplotSwitch['DotColor'] != 'PC1_2_represen'):#forで回して全部描画する

        host2[1,4].tick_params(labelbottom="off",bottom="off") # x軸の削除
        host2[1,4].tick_params(labelleft="off",left="off") # y軸の削除
        host2[1,4].tick_params(bottom=False, left=False,  right=False,     top=False)
        host2[1,4].spines["right"].set_color("none")  # 右消し
        host2[1,4].spines["left"].set_color("none")   # 左消し
        host2[1,4].spines["top"].set_color("none")    # 上消し
        host2[1,4].spines["bottom"].set_color("none") # 下消し
        #host2[1,4].set_xticklabels([]).box("off") #枠線の削除
        t=np.linspace(0,1,pca_score.shape[0])
        x = np.linspace(0,0,pca_score.shape[0])
        yy=np.linspace(0,0,pca_score.shape[0])
        ax = plt.scatter(x,yy,c=t,cmap=cm.jet,marker='.',lw=0)#vmin=min(pca_score),vmax=max(pca_score))
        #cax = fig.add_axes([0.9,0.8,2.5,0.08])        #[左端、下端、幅、高さ]
    
        #divider = make_axes_locatable(ax)
        #ax_cb = divider.new_horizontal(size="2%", pad=0.05)
        CF = plt.colorbar(ax)
        CF.ax.tick_params(labelsize=10)
        #cbar = plt.colorbar(CF)
        CF.ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False, color='white')#plt.tick_params(color='white')
        for j, lab in enumerate(np.linspace(minscore2,maxscore2,9)): 
            CF.ax.text(2.5, (j )/9+0.05, round(lab,2), ha='center', va='center', size=20) 
        #plt.title('PC'+str(currPC+1)+' ('+str(cov_ratio*100)[0:4]+'%)', size=20)
                  #を説明')

    

    plt.subplots_adjust(wspace=0.6, hspace=0.6)
    plt.savefig(save_dir + 'LoadingPlotWCluster'+str(currPC)+'_'+str(othrPC)+Name+'.pdf') 
    plt.savefig(save_dir + 'LoadingPlotWCluster'+str(currPC)+'_'+str(othrPC)+Name+'.png') 
    currplot =0 
    plt.close()
                
def FactorLoadingPlotWCluster(timepoint,pca_coef,pca_score,currPC,cov_ratio,ylim,MolLabel,label,BiplotSwitch,ClstColorDF):#FactorLoadingの時間波形をplot,クラスタごとでで色分け
    import matplotlib.colors as clr
    import matplotlib as mpl
    MolColorDF = BiplotSwitch['MolColor']
    #ClstColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/20180317/ClstColor.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)

    ClstColorList=list(ClstColorDF.loc[0])
    #ClstColorList.reverse()
    cmap = mpl.colors.ListedColormap(ClstColorList)

    for ii in range(pca_score.shape[0]):#まず主成分数を変えた時のFactorLoadingのスコアによる波形の違い,84種
        if BiplotSwitch['DotColor'] == 'ClstColor':
            #PC1_2
            #plt.plot(timepoint,pca_score[ii,0] * pca_coef.T[0,:] + pca_score[ii,1] * pca_coef.T[1,:],color=BiplotSwitch['ColorDF'].loc[MolLabel[ii]]['ClsterColor'],linewidth=0.5, marker = 'o',markersize=3)
            #allPC
            #plt.plot(timepoint, pca_coef[ii,:]  ,color=BiplotSwitch['ColorDF'].loc[MolLabel[ii]]['ClsterColor'],linewidth=0.5, marker = 'o',markersize=3)
                #EachPC
            plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['ColorDF'].loc[MolLabel[ii]]['ClsterColor'],linewidth=0.5, marker = 'o',markersize=3)
                    #axis=['top','bottom','left','right']
            #line_width=[0, 0.5, 0.5, 0]
            
            #for a,w in zip(axis, line_width):  # change axis width
             #   ax.spines[a].set_linewidth(w)
        elif BiplotSwitch['DotColor'] == 'EachMol':
            Color=BiplotSwitch['DotColor']

            plt.plot(timepoint, pca_coef*pca_score[ii]  ,color=BiplotSwitch['EachMolColor'][ii],linewidth=0.5, marker = 'o',markersize=3)
            plotdata = pca_coef*pca_score[ii]
            if any(plotdata > 5):#BiplotSwitch['dotname'][ii]
                aa= np.where(plotdata > 5)[0]
                for jj in aa:
                
                    plt.text(timepoint[jj], plotdata[jj], label[ii],size=1)
            elif any(plotdata < -5):#BiplotSwitch['dotname'][ii]
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
    #plt.ylabel('Valiation Index',size=20)

    t=np.linspace(0,1,pca_score.shape[0])
    x = np.linspace(0,0,pca_score.shape[0])
    yy=np.linspace(0,0,pca_score.shape[0])
    ax = plt.scatter(x,yy,c=t,cmap=cmap,marker='.',lw=0)#vmin=min(pca_score),vmax=max(pca_score))
    #cax = fig.add_axes([0.9,0.8,2.5,0.08])        #[左端、下端、幅、高さ]
    
    CF = plt.colorbar(ax)
    #cbar = plt.colorbar(CF)
    CF.ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False, color='white')#plt.tick_params(color='white')
    for j, lab in enumerate(['1','2','3','4','5','6','7','8','9','10','11','12','13']): 
        CF.ax.text(1.9, ( j )/13 + 0.03 , lab, ha='center', va='center', size=20) 
    plt.title('PC'+str(currPC+1)+' ('+str(cov_ratio*100)[0:4]+'%)', size=20)
              #を説明')
    plt.ylim([-ylim,ylim])
    plt.xlabel('Time(min)',size=20)
def mkZscore(CorrFastPropDF,label,labelProp):
    meanNP = np.array(CorrFastPropDF) -np.nanmean(np.array(CorrFastPropDF),axis=0,keepdims=True)
    stdNP = np.nanstd(np.array(CorrFastPropDF).astype(float),axis=0)
    signDF = pd.DataFrame(index=label,columns=labelProp)

    for i in range(len(labelProp)):
        #if i == 0:
         #   signNP = np.array([0]*len(label))
        #else:    
        signNP = meanNP[:,i]/stdNP[i]

        signDF[labelProp[i]] = signNP
    return(signDF)
    
def PCAprep(XDF,ClstAveTimeSeriesDF,ClstColorDF,MolLabel,Label,ClstMolDF,AnalSwitch,LoadingOpt,BiplotSwitch,OptionSwitch,CompareBolus,ClstNoColorDF,save_dir):
    
    NormlzDF = mkZscore(XDF,MolLabel,Label)#正規化する#列方向たてに正規化する=相関行列が返ってくる
    X = np.array(NormlzDF)
    #X=XDF
    XDF=XDF.T
    pca_score,pca_coef,cov_ratio,pca_components_,W,v = PCAHel.PCA(X,2)#この時点でscore,loading,label,LabelPropの順番は対応している


    try:    
        if CompareBolus['Target'] == 'Betw':#VIときは調整する
            pca_score[:,0] = pca_score[:,0]*-1#PC1の各分子のスコアは列方向にある
            pca_coef[:,0] = pca_coef[:,0]*-1#PC1の各分子のローディングは列方向にある
            pca_components_[:,0] = pca_components_[:,0]*-1#Factor Loadingではなく、Loading（固有ベクトル）
    except:
        pass
    label=MolLabel
    #pca_score=pca_score*-1

    if BiplotSwitch['DotColor'] == 'EachMol':
        BiplotSwitch['EachMolColor'] = PCAHel.mkEachMolColor(pca_score,XDF)#:分子x○○PCAにおける、各分子でのいろわえk
    if BiplotSwitch['DotColor'] == 'EachMolCluster':
        BiplotSwitch['EachMolColor'] = PCAHel.mkEachMolClsterColor(pca_score,XDF,ClstColorDF,ClstMolDF)#:分子x○○PCAにおける、各分子のClusterでの色で色分け
    if BiplotSwitch['DotColor'] == 'ClstColor':
        BiplotSwitch['EachMolColor'] = ClstNoColorDF
    if BiplotSwitch['DotColor'] == 'BC':#PCA時に色分けBCにする
        BiplotSwitch['EachMolColor'] = PCAHel.mkEachBCColor(pca_score,XDF)#:分子x○○PCAにおける、各分子のClusterでの色で色分け
    if BiplotSwitch['DotColor'] == 'BWoC':#PCA時に色分けBCにする
        BiplotSwitch['EachMolColor'] = PCAHel.mkEachBColorWOC(pca_score,XDF)#:分子x被験者PCAにおける、Bのみ色付き
    if BiplotSwitch['DotColor'] == 'CWoB':#PCA時に色分けBCにする
        BiplotSwitch['EachMolColor'] = PCAHel.mkEachCColorWOB(pca_score,XDF)#:分子x被験者PCAにおける、Bのみ色付き      
    if AnalSwitch['LengVector'] == 1:#PC1,2score間のベクトル長を算出する。
        CompareBolus['lengthPC12'] = PCAHel.calclenbetwScore(pca_score,BiplotSwitch,CompareBolus,MolLabel,save_dir)   
        #さらにB,Cそれぞれの原点からの距離も算出、棒グラフと3Dplot?
        PCAHel.calclenScore(pca_score,BiplotSwitch,CompareBolus,MolLabel,save_dir)
        
    else:#しないならデフォルトのマーカーサイズ->PCAHel.PlotAngle
        CompareBolus['lengthPC12']=20
        
    if AnalSwitch['PlotPCACovRatio'] == 1:#各PCの寄与率
        fig = plt.figure(figsize=(10,8))
        PCAHel.plotPCACovRatio(cov_ratio)
        plt.xticks(size='40');plt.yticks(size='40')
        plt.savefig(save_dir + 'PCExplained.pdf',bbox_inches="tight")
        plt.savefig(save_dir + 'PCExplained.png',bbox_inches="tight")
    if AnalSwitch['Biplot'] == 1:#PC平面プロっト
        cov_ratioIdx = 0
        #ClstMolDF=[]
        #while np.cumsum(x_covratio)[cov_ratioIdx] < 0.8:#ほんとはiを固定して、i+1,2...とのPC平面でプロットしたい
        for i in range(min(np.where(np.cumsum(cov_ratio)>0.85)[0])):#最後のの相手0,1,2,3
            for j in range(i+1,min(np.where(np.cumsum(cov_ratio)>0.85)[0])+1):#1,2,3,..
                fig = plt.figure()
                #ax = plt.subplot(111, aspect='equal')#fig.axes()
                #pca_coef = np.transpose(pca.components_)
                fig.set_size_inches(pca_score.shape[1], pca_score.shape[0]/5.0)
                PCAHel.PCAScatter(pca_score,label,i,j,ClstMolDF,label,BiplotSwitch)
                fig.tight_layout() 
                ##PCAHel.plotPCAcoef(pca_coef,Label,pca_score,i,j,MolLabel,BiplotSwitch, linecolor='r',fontsize=10)
                plt.xlabel('PC'+str(i+1)+' ('+str(cov_ratio[i]*100)[0:4]+'%)')
                plt.ylabel('PC'+str(j+1)+' ('+str(cov_ratio[j]*100)[0:4]+'%)')
                plt.xticks([-10,-7.5,-5.0,-2.5,0,2.5,5,7.5,10],['-10','-7.5','-5.0','-2.5','0','2.5','5.0','7.5','10'],size=15)
                #plt.yticks([-6,-4,-2,0,2,4,6],['-6.0','-4.0','-2.0','0','2.0','4.0','6.0'],size=15)
                plt.yticks([-4,-2,0,2,4],['-4.0','-2.0','0','2.0','4.0'],size=15)

                plt.xticks(size='10');plt.yticks(size='10')
                plt.ylim([-5,5])#昔は-7,6
                plt.xlim([-11,11])#昔は-7,6
                axis=['top','bottom','left','right']
                line_width=[1,1,1,1]
                
                for a,w in zip(axis, line_width):  # change axis width
                    plt.gca().spines[a].set_linewidth(w)
                plt.savefig(save_dir + 'PCAScatter'+ str(i+1) + 'vs' + str(j+1) + '.pdf')
                if AnalSwitch['DrawEllipse']  == 1:#楕円フィッティング
                    XX = ClstAveTimeSeriesDF
                    PCAHel.drawEllipse(XX,pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF)
                    plt.ylim([-5.5,5.5])#昔は-7,6
                    plt.xlabel('PC'+str(i+1)+' ('+str(cov_ratio[i]*100)[0:4]+'%)')
                    plt.ylabel('PC'+str(j+1)+' ('+str(cov_ratio[j]*100)[0:4]+'%)')
                    plt.xticks(size='10');plt.yticks(size='10')
                    if AnalSwitch['WHist']==1:#ヒストグラムつきにする
                        Optiondict={'calcR':'',#'spearman'
                        'xlabel':'',
                            'ylabel':'',
                            'Annotate':0,
                            'Label':[],
                            'title':'',
                            'y=x' :0 #y=xの線を足す
                            }
                        GH.mkScatterWHist(list(pca_score[:,0]),list(pca_score[:,1]),save_dir,['k']*len(pca_score[:,0]),Optiondict)#2つのリストの散布図+ヒストグラム
                        plt.ylim([-5.5,5.5])
                    axis=['top','bottom','left','right']
                    line_width=[1,1,1,1]
                    
                    for a,w in zip(axis, line_width):  # change axis width
                        plt.gca().spines[a].set_linewidth(w)

                    plt.savefig(save_dir + 'PCAScatterWEllipse'+ '_PC'+str(i+1)+ 'vs' + 'PC'+str(j+1) + '.pdf',bbox_inches="tight")
                try:
                    if AnalSwitch['VectorDiagram']  == 1:#Bolus->Continuouのベクトル図
                        PCAHel.VectorDiagram(pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF)                
                        plt.savefig(save_dir + 'PCAScatterWVectorDiagram'+ '_PC'+str(i+1)+ 'vs' + 'PC'+str(j+1) + '.pdf',bbox_inches="tight")
                        
                        PCAHel.PlotAngle(pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF,BiplotSwitch,CompareBolus)  ## Bolus->Rampの変化角θと仰角Φのscatter
                        plt.savefig(save_dir + 'PCAScatterWAngle'+ '_PC'+str(i+1)+ 'vs' + 'PC'+str(j+1) + '.pdf',bbox_inches="tight")
                except:
                    pass
                if AnalSwitch['calcinner_outer_product'] ==1:##各平面でのBC間の内積外積を算出
                    inner_product, outer_product,AngleList = PCAHel.calcinner_outer_product(pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF)
                    Title = 'InnerProduct'; xlabel = ''; ylabel = 'InnerProduct';  size=20; Titlesize=20; xsize=10;
                    MolColor = OptionSwitch['MolColor']
                    MolColor['IP'] = inner_product
                    MolColor['OP'] = outer_product
                    MolColor['Angle'] = AngleList
                    List1 = list(MolColor['IP'].sort_values()); xticks = list(MolColor['IP'].sort_values().index); Color = list(MolColor.sort_values('IP')['MolColor'])

                    GH.mkSortedBarWHist(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
                    
                    List2 = list(MolColor['OP'].sort_values()); xticks = list(MolColor['OP'].sort_values().index); Color = list(MolColor.sort_values('OP')['MolColor'])
                    Title = 'OuterProduct'; xlabel = ''; ylabel = 'OuterProduct';
                    
                    GH.mkSortedBarWHist(List2, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
                    GH.mkScatterWHist(List1,List2,save_dir,Color,dict({'calcR':'pearson','xlabel':'inner_product','ylabel':'outer_product','Annotate':0,'title':''}))#2つのリストの散布図+ヒストグラム

                    ### Angleの棒グラフ
                    MolColor['lengthPC12'] = CompareBolus['lengthPC12']
                    List3 = list(MolColor['Angle'].sort_values()); xticks = list(MolColor['Angle'].sort_values().index); Color = list(MolColor.sort_values('Angle')['MolColor'])
                    Title = 'Angle'; xlabel = ''; ylabel = 'Angle';
                    ListBC = list(MolColor.sort_values('Angle')['lengthPC12'])
                    
                    GH.mkSortedAngleBarWHist(List3, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
                    ### AngleとBCNormの散布図
                    Optiondict={'xlabel':'BCNorm','ylabel':'Angle','Annotate':1,'Label':xticks,'calcR':''}
                    Optiondict['title'] = 'BC vs Angle';Optiondict['y=x'] = 0#y=xの線を足す
                    GH.mkScatterYAngleWHist(ListBC,List3,save_dir,Color,Optiondict)#2つのリスト(2つ目は角度）の散布図+ヒストグラム                   
                    ###  sin(Anglue)の棒グラフ
                    MolColor['sin'] = np.sin(list(MolColor['Angle']));
                    List4 = list(MolColor['sin'].sort_values()); xticks = list(MolColor['sin'].sort_values().index); Color = list(MolColor.sort_values('sin')['MolColor'])
                    ListBC = list(MolColor.sort_values('sin')['lengthPC12'])
                    
                    Title = 'sin'; xlabel = ''; ylabel = 'sin';                    
                    GH.mkSortedBarWHist(List4, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
                    ###  sin(Anglue)とBCNormの散布図
                    #Optiondict['markersize'] = [i*100 for i in CompareBolus['lengthPC12'] ]
                    Optiondict={'xlabel':'BCNorm','ylabel':'sin(θ)','Annotate':1,'Label':xticks,'calcR':''}
                    Optiondict['title'] = 'BC vs sin(θ)';Optiondict['y=x'] = 0#y=xの線を足す

                    GH.mkScatterWHist(ListBC,List4,save_dir,Color,Optiondict)#2つのリストの散布図+ヒストグラム

                    MolColor.to_excel(save_dir+  'innner_outer_product.xlsx')
        cov_ratioIdx += 1
        plt.close()     
    #if AnalSwitch['PC1-2Table'] == 1:#PC1-2Table作成
     #   pass
      #  PCAHel.PC1_2()
    if np.sum(cov_ratio[0:2])<0.8:
        PCAHel.plot3DPCA(label,pca_score,pca_coef,Label,BiplotSwitch,save_dir)
    
    if AnalSwitch['LoadingPCsimulation'] ==1: #PC1,2固定して片方動かした時の時系列描画
        currplot =0  
        max_score=max(np.abs(pca_score[:,0]))*max(np.abs(pca_components_[:,0]))*1.2

        #for i in [0,1]:#PC1,2どちらも入れないと意味ない
        LoadingPCsimulation(list(XDF.index),pca_components_[:,0:2],pca_score,0,1,cov_ratio,max_score,MolLabel,label,BiplotSwitch,ClstColorDF,save_dir)#あるPCのスコア入れたらある数降って、別のPCの分だけ別に描画
            #pca_score[locIdx,currplot],currplot,cov_ratio[currplot],max_score,MolLabel,label,BiplotSwitch
        LoadingPCsimulation(list(XDF.index),pca_components_[:,0:2],pca_score,1,0,cov_ratio,max_score,MolLabel,label,BiplotSwitch,ClstColorDF,save_dir)#あるPCのスコア入れたらある数降って、別のPCの分だけ別に描画
     
    if AnalSwitch['LoadingPlotWCluster'] ==1:#Loadingをプロットする色々によって色分け
        BiplotSwitch['MolColor'] = OptionSwitch['MolColor']
        #subplotの諸設定
        pca_components_#Factor Loadingではなく、Loading（固有ベクトル）
        max_score=3.5#max(np.abs(pca_score[:,0]))*max(np.abs(pca_components_[:,0]))*1.2
        subplot_row = 3
        subplot_col = 3
        fig, host = plt.subplots(subplot_row,subplot_col,figsize=(18, 18))
        #fig.set_size_inches(subplot_row*3, subplot_col*5)
        #fig.set_size_inches(10,5)
        totalfignum = pca_score.shape[1]
        currplot =0  
        
        if OptionSwitch['Target'] == 'EachCluster':
            #Clusterに沿った分子のindexを取得
            for ii in range(len(ClstMolDF.columns)
            ):#Clusterの数だけ、
                locIdx=[XDF.columns.get_loc(i) for i in list(ClstMolDF['cluster_'+str(ii+1)].dropna())]
                fig, host = plt.subplots(subplot_row,subplot_col,figsize=(18, 18))

                currplot =0   
                for i in range(subplot_row):
                    for j in range(subplot_col):                
                        plt.subplot(subplot_row, subplot_col, 1+j+i*subplot_col)
                        MolLabel=list(ClstMolDF['cluster_'+str(ii+1)].dropna())

                        FactorLoadingPlotWCluster(list(XDF.index),pca_components_[:,currplot],pca_score[locIdx,currplot],currplot,cov_ratio[currplot],max_score,MolLabel,label,BiplotSwitch,ClstColorDF)
                        #plt.xticks([0,60,120,240], ['0','60','120','240'],size=20)#,rotation='20')#'vertical')
                        plt.yticks(size='20')
                        #plt.ylim([- max_score, max_score])
                        #fig.tight_layout()
                        [host[i][j].spines[ii].set_linewidth(0.1) for ii in host[i][j].spines.keys()]
                        currplot +=1
                        if currplot >= totalfignum:
                            break
                    if currplot >= totalfignum:
                        break
                plt.subplots_adjust(wspace=0.6, hspace=0.6)
                plt.savefig(save_dir + 'LoadingPlotWCluster'+str(ii)+'.pdf') 
                plt.savefig(save_dir + 'LoadingPlotWCluster'+str(ii)+'.png') 
                currplot =0 
                plt.close()               
        else:
            for i in range(subplot_row):
                for j in range(subplot_col):                
                    plt.subplot(subplot_row, subplot_col, 1+j+i*subplot_col)

                    FactorLoadingPlotWCluster(list(XDF.index),pca_components_[:,currplot],pca_score[:,currplot],currplot,cov_ratio[currplot],max_score,MolLabel,label,BiplotSwitch,ClstColorDF)
                    #plt.xticks([0,60,120,240], ['0','60','120','240'],size=20)#,rotation='20')#'vertical')
                    plt.yticks(size='20')
                    #plt.ylim([- max_score, max_score])
                    fig.tight_layout() 
                    currplot +=1
                    if currplot >= totalfignum:
                        break
                if currplot >= totalfignum:
                    break
                #[host[i][j].spines[ii].set_linewidth(0.05) for ii in host[i][j].spines.keys()]
                host[i][j].spines['bottom'].set_linewidth(0.1); host[i][j].spines['left'].set_linewidth(0.1); host[i][j].spines['right'].set_linewidth(0.1); host[i][j].spines['top'].set_linewidth(0.1); 

            plt.subplots_adjust(wspace=0.4, hspace=0.4)
            
            plt.savefig(save_dir + 'LoadingPlotWCluster.pdf') 
            plt.savefig(save_dir + 'LoadingPlotWCluster.png') 
            plt.close()
            
    if AnalSwitch['ScoreHeatMap'] ==1:#Scoreのヒートマップ
        fig, ax = plt.subplots()
        #pca_score = np.log10(pca_score)
        fig.set_size_inches(pca_score.shape[1], pca_score.shape[0]/5.0)
        ScoreDF = pd.DataFrame(pca_score)
        ScoreDF.to_excel(save_dir+'Score.xlsx')
        PCAHel.heatPCA(pca_score,label)#npでぶち込む    
        MolColorDF =OptionSwitch['MolColor']
        if len(pca_score[:,0]) > len(MolColorDF['MolColor']):#雑だが、被験者x分子とかの時
            #if re.compile("(.*)(_)(.*)").search(jj).group(3) in list(PosDF.columns) :
             #   ii=re.compile("(.*)(_)(.*)").search(jj).group(3)            
            splabel = ax.yaxis.get_ticklabels()
            for i in splabel:
                #print(i.get_text())
                nowlabel = re.compile("(.*)(_)(.*)").search(i.get_text()).group(1)
                i.set_color(MolColorDF['MolColor'][nowlabel])
        else:
            [t.set_color(i) for (i,t) in
             zip(list(MolColorDF['MolColor']),ax.yaxis.get_ticklabels())]
        #fig.tight_layout()
        plt.savefig(save_dir + 'PCAScoreHeatmap.pdf',bbox_inches="tight")
        plt.savefig(save_dir + 'PCAScoreHeatmap.png',bbox_inches="tight")

        plt.close()
        
    if AnalSwitch['FactorLoadingHeatMap']==1:#FactorLoadingのヒートマップ
        ## Loading Heatmap
        fig = plt.figure(figsize=(5,10))
        pca_coef[0,-1] = 0
        #pca_coef = np.log10(pca_coef)
        if LoadingOpt == 0:#データ（縦）×主成分（横）
            fig.set_size_inches(pca_coef.shape[1]*1.2, pca_coef.shape[0]/2)
            PCAHel.heatPCA(pca_coef, Label)
        elif LoadingOpt == 1:#主成分（縦）×データ（横）
            fig.set_size_inches(pca_coef.shape[0]/5, pca_coef.shape[1])
            heatPCAhorz(pca_coef.T, Label)
        fig.tight_layout() 
        plt.xticks(size='40');plt.yticks(size='40')
        plt.ylabel('Time(min.)',size=40)
        plt.savefig(save_dir + 'PCAFactorLoadingHeatmap.pdf',bbox_inches="tight")
        plt.savefig(save_dir + 'PCAFactorLoadingHeatmap.png',bbox_inches="tight")

        plt.close()
        
        
    if AnalSwitch['ScoreVariance'] == 1:#'ScoreVariance':各分子、被験者間のスコアの分散
            
        ScoreDF= pd.DataFrame(data=pca_score.T,index= np.arange(1,len(pca_score[0])+1),columns = list(XDF.columns))
        ScoreDF.to_excel(save_dir+'ScoreVarDF.xlsx')
        if CompareBolus['TargetDetail'] == 'EachSubject':#BC各被験者なら
            PCAHel.calcScoreVarBC(ScoreDF,label,save_dir)#DFでぶち込む 
            
        else:
            PCAHel.calcScoreVar(ScoreDF,label,save_dir)#DFでぶち込む 
        