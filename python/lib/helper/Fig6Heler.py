# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 11:00:04 2017
すべての分子の空腹値とプロパティ間の相関、スキャッターを表示
今後基本的にBolusデータのみ扱う
@author: fujita

"""
    
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import pandas as pd
from Helper import general as ge
from DataLoader import mkBolusIdx as BI
from pylab import *
import matplotlib.font_manager
from scipy.stats import pearsonr
import itertools
from DataLoader import AdjustData as aD
import matplotlib.cm as cm
import re
import MolTimeSeriesCV as MTSCV
import scipy.stats as sp
import Helper.StatCal as SC
import Helper.AnalMolTime as AmT
import MatrHelper as MatHel
import collections
import Helper.GraphHelper as GH
import Helper.MolCorrSubjmeanStdHelper as MSMSH
import Helper.ISHelper as ISH
from scipy.stats import zscore
import DataLoader.ChangeNameEng as CNE
import mkGraph 
import Helper.LabelHeler as LH


def PlotHist(SubjRmean,Col,c,ColDict,NormlSwitch):#ヒストグラムを描画
    #red = ColDict['red'];
    print(Col)
    blue = ColDict['blue']; green = ColDict['green']
    if NormlSwitch == 0:#とりあえず簡単に片方の代謝系で色分け__ヒストグラム描画            
            fig2, ax2 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = abs(SubjRmean[Col[0]:Col[(c['red'])]]).values.flatten(); tempnp4 = abs(SubjRmean[Col[green]:Col[blue]]).values.flatten(); tempnp2 = abs(SubjRmean[Col[(c['red'])]:Col[magenta]]).values.flatten(); tempnp3 = abs(SubjRmean[Col[magenta]:Col[green]]).values.flatten();  #tempnp5 = abs(SubjRmean[Col[blue]:
            ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
            plt.savefig(save_dir + 'AbsHistAve_Singgle.pdf')        
        #Stdもやる
            fig3, ax3 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = abs(SubjRstd[Col[0]:Col[(c['red'])]]).values.flatten();  tempnp4 = abs(SubjRstd[Col[green]:Col[blue]]).values.flatten(); tempnp2 = abs(SubjRstd[Col[(c['red'])]:Col[magenta]]).values.flatten();  tempnp3 = abs(SubjRstd[Col[magenta]:Col[green]]).values.flatten();   tempnp5 = abs(SubjRstd[Col[blue]:Col[purple]]).values.flatten(); tempnp6 = abs(SubjRstd[Col[purple]:Col[black]]).values.flatten()
            ax3 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
            plt.savefig(save_dir + 'AbsHistStd_Single.pdf')  
        #絶対値ではない           
            fig2, ax2 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = SubjRmean[Col[0]:Col[(c['red'])]].values.flatten(); tempnp4 = SubjRmean[Col[green]:Col[blue]].values.flatten(); tempnp2 = SubjRmean[Col[(c['red'])]:Col[magenta]].values.flatten(); tempnp3 = SubjRmean[Col[magenta]:Col[green]].values.flatten();  tempnp5 = SubjRmean[Col[blue]:Col[purple]].values.flatten(); tempnp6 = SubjRmean[Col[purple]:Col[black]].values.flatten()
            ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
            plt.savefig(save_dir + 'HistAve_Single.pdf')        
        #Stdもやる
            fig3, ax3 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = SubjRstd[Col[0]:Col[(c['red'])]].values.flatten();  tempnp4 = SubjRstd[Col[green]:Col[blue]].values.flatten(); tempnp2 = SubjRstd[Col[(c['red'])]:Col[magenta]].values.flatten();  tempnp3 = SubjRstd[Col[magenta]:Col[green]].values.flatten();  tempnp5 = SubjRstd[Col[blue]:Col[purple]].values.flatten();  tempnp6 = SubjRstd[Col[purple]:Col[black]].values.flatten()
            ax3 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
            plt.savefig(save_dir + 'HistStd_Single.pdf')  
            
            #両方がその代謝系で、それ以外はグレー
            fig4, ax4 = plt.subplots(1,1)
            ConList1 = [listr1, listm1, listg1, listb1, listpl1, listbl1, listgr1];         ConList2 = [listr2, listm2, listg2, listb2, listpl2, listbl2, listgr2];
    
            ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
    
            ax4 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black','gray'])        
            plt.savefig(save_dir + 'HistAveSameMol.pdf')  
    
            fig5, ax5 = plt.subplots(1,1)
            ax4 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])],ConList2[6][~np.isnan(ConList2[6])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black','gray'])        
    
            #ax5 = plt.hist([listr2, listm2, listg2, listb2, listbl2, listgr2], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'])        
            plt.savefig(save_dir + 'HistStdSameMol.pdf')  
    
            fig6, ax6 = plt.subplots(1,1)
            ax6 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
            plt.savefig(save_dir + 'HistAveSameMolWOGray.pdf') 
    
            fig7, ax7 = plt.subplots(1,1)
            ax7 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
            plt.savefig(save_dir + 'HistStdSameMolWOGray.pdf')  ;print(ConList2[4][~np.isnan(ConList2[4])])
                                 
            fig8, ax8 = plt.subplots(1,1)
            ConList1 = [listrabs1, listmabs1, listgabs1, listbabs1, listplabs1, listblabs1, listgrabs1];         ConList2 = [listrabs2, listmabs2, listgabs2, listbabs2,listplabs2, listblabs2, listgrabs2];
    
            ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
    
            ax8 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black','gray'])        
            plt.savefig(save_dir + 'AbsHistAveSameMol.pdf')  
    
            fig9, ax9 = plt.subplots(1,1)
            ax9 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])],ConList2[6][~np.isnan(ConList2[6])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black','gray'])        
    
            #ax5 = plt.hist([listr2, listm2, listg2, listb2, listbl2, listgr2], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'])        
            plt.savefig(save_dir + 'AbsHistStdSameMol.pdf')
            
            fig10, ax10 = plt.subplots(1,1)
            ax10 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
            plt.savefig(save_dir + 'AbsHistAveSameMolWOGray.pdf') 
    
            fig11, ax11 = plt.subplots(1,1)
            ax11 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
            plt.savefig(save_dir + 'AbsHistStdSameMolWOGray.pdf') 

def EachSubjAnal(SubjRPanel,SubjectName):
    ColorSwitchDict  = {1:'ClstColor',2:'MolColor',3:'TimeVarColor'}#1でクラスタ、2で分子,3で時間ばらつき
    SwitchDict = {'MolColorScatter':1 }#0でやらない
    LabelSum =   pd.read_excel("//Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng.xlsx",header=0,encoding = "ISO-8859-1")
    OptionDict = {'Edge':''}#'Subject'}    
    NormlSwitch=0 #0で　:#とりあえず簡単に片方の代謝系で色分け__ヒストグラム描画　1で代謝ごとの累積分布、2で分布の割合
    for i in range(20):#SubjRPanel[SubjectName[i]]
        
        SubjRmeanrev = CNE.ChnageNametoEnglish(SubjRPanel[SubjectName[i]],2);# SubjRstdrev = CNE.ChnageNametoEnglish(SubjRstdrev,2)
        SubjRmeanrev = MatHel.adjustMatrlower(SubjRmeanrev)
        PlotScatter(save_dir+SubjectName[i]+'_',SubjRmeanrev,SubjRmeanrev,ColorSwitch,NormlSwitch,SwitchDict)
        #GH.mkHist(SubjRmeanrev,save_dir)
        #CombDF = MSMSH.mkExcelupper(SubjRmeanrev,SubjRmeanrev,0.8)#ある値以上の分子名組み合わせをエクセルにまとめる           
        #CombDF.to_excel(save_dir + 'CombDf' + SubjectName[i]+ '.xlsx')
        #mkGraph.mkNetWorkGraph(SubjRmeanrev,LabelSum,0.8,ColorSwitchDict[2],OptionDict,SubjectName,save_dir+SubjectName[i]+'_')
        fig = plt.figure()
        fig.set_size_inches(np.array(SubjRmeanrev).shape[1]/5, np.array(SubjRmeanrev).shape[0]/6)
        #ASTH.heatmap(np.array(SubjRmeanrev),list(SubjRmeanrev.columns))
        fig.tight_layout() 
        #plt.savefig(save_dir+'HeatMap_'+SubjectName[i]+'.pdf')        

def perform4plot(List,nump):
    #List = np.insert(List, 0, 0)   
    #if len(np.where(List==1)[0]) > 1:  #なぜか最後に1が多ければ除く 
     #   List = np.delete(List,np.where(List==1)[0][-1])
    delIdx = check(List)
    #stIdx, endIdx, delIdx   =  ncalc4plot(List)
    
    #if len(stIdx) > 0:
     #   List = np.insert(List, stIdx, 0)
    #if len(endIdx) > 0:
     #   List = np.insert(List, endIdx, 0)
    #CumList = np.add.accumulate(List) / List.sum()    
    if len(delIdx) > 0:
        nump=np.delete(nump,delIdx)
        List =np.delete(List,delIdx)
    #CumList = np.insert(CumList, 0, 0)    
    
    return(List, nump)

def ncalc4plot(List):
    if len(List) < 21:
        first = List[0]
        end = List[-1]
        stIdx=[]
        endIdx=[]
        delIdx=[]
        count=0
        if first == 0:#始めが0なら
                stIdx.append(0)
        elif end == 0:#終わりが0なら
                endIdx.append(20)        
        for num in List[1:]:
            #if (num == 0) and (last==0) and (first==0):#始めから0が続く限り
             #   stIdx.append(count)
            #elif (num == 0) and (last!=0) and (first==0):#続かなくなったらflag
             #   flag4end=1
            #if (num == 0) and (last==0) and (flag4end==1):#途中から終わりまで0が続く限り
             #   endIdx.append(count)
            if num==0:#途中の点
                delIdx.append(count)
            last = num
            count += 1
    else:
        stIdx=[];endIdx=[];delIdx=[]
    return stIdx, endIdx, delIdx      

def check(lst):
    first = lst[0]
    second = lst[1]
    Idx=[]#重複してる箇所
    count=1#リスト中の3こ以上連続するとこは両端のみ残すようにする
    for num in lst[2:]:
        if (num == first) and (first == second):#3つ同じなら
            Idx.append(count)
        count += 1    
        first = second
        second = num
    return Idx

def CumsumCountList(List):#このリストが0.1の0.05刻みの中に何個あるかのリストを自作する、
    CumSumList=[]
    for i in [0.05*x for x in range(21)]:  
       if  np.where( (List < i) & (List > i-0.05) )[0].any() == 1:#この刻み幅中に含まれるなら
           CumSumList.append(len( np.where( (List < i) & (List > i-0.05 ) )[0] ))
       else:
           CumSumList.append(0)
    return(np.cumsum(CumSumList) / sum(CumSumList)) 
def mksevengrayList(SubjRmean,Col,CnumList):#異なる代謝物グループ同士のリストの作成
    list1=[];list2=[];list3=[];list4=[];list5=[];list6=[];list7=[];list8=[]
    for i in range(0,len(SubjRmean.columns)):#列を1~84まで
        if i < (CnumList[1]):
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[1]-1)]]);
            list1 += List;
            stloc = Col[CnumList[1]]
        elif i < (CnumList[2]):
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[2]-1)]]);
            list2 += List;            
            stloc = Col[CnumList[2]]
        elif i < (CnumList[3]):  
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[3]-1)]]);
            list3 += List;              
            stloc = Col[CnumList[3]]
        elif i < (CnumList[4]):
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[4]-1)]]);
            list4 += List;               
            stloc = Col[CnumList[4]]
        elif i < (CnumList[5]):            
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[5]-1)]]);
            list5 += List; 
            stloc = Col[CnumList[5]]
        elif i < (CnumList[6]):            
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[6]-1)]]);
            list6 += List; 
            stloc = Col[CnumList[6]-1]
        elif i < (CnumList[7]):            
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[7]-1)]]);
            list7 += List; 
            stloc = Col[CnumList[7]-1]            
        List=list(SubjRmean.loc[Col[i],stloc:]);
        list8+=List
    return([list1,list2,list3,list4,list5,list6,list7,list8])
def  mksevenabsList(SubjRmean,Col,CnumList):

    tempnp1 = abs(SubjRmean[Col[CnumList[0]]:Col[CnumList[1]-1]]).values.flatten(); 
    tempnp4 = abs(SubjRmean[Col[CnumList[3]]:Col[CnumList[4]-1]]).values.flatten(); 
    tempnp2 = abs(SubjRmean[Col[CnumList[1]]:Col[CnumList[2]-1]]).values.flatten(); 
    tempnp3 = abs(SubjRmean[Col[CnumList[2]]:Col[CnumList[3]-1]]).values.flatten();  
    tempnp5 = abs(SubjRmean[Col[CnumList[4]]:Col[CnumList[5]-1]]).values.flatten(); 
    tempnp6 = abs(SubjRmean[Col[CnumList[5]]:Col[CnumList[6]-1]]).values.flatten()  
    tempnp7 = abs(SubjRmean[Col[CnumList[5]]:Col[CnumList[7]-1]]).values.flatten()  

    return(tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6,tempnp7 )

def mksevenList(SubjRmean,Col,CnumList):
    tempnp1 = SubjRmean[Col[CnumList[0]]:Col[CnumList[1]-1]].values.flatten(); 
    tempnp4 = SubjRmean[Col[CnumList[3]]:Col[CnumList[4]-1]].values.flatten(); 
    tempnp2 = SubjRmean[Col[CnumList[1]]:Col[CnumList[2]-1]].values.flatten(); 
    tempnp3 = SubjRmean[Col[CnumList[2]]:Col[CnumList[3]-1]].values.flatten();  
    tempnp5 = SubjRmean[Col[CnumList[4]]:Col[CnumList[5]-1]].values.flatten(); 
    tempnp6 = SubjRmean[Col[CnumList[5]]:Col[CnumList[6]-1]].values.flatten()  
    tempnp7 = SubjRmean[Col[CnumList[5]]:Col[CnumList[7]-1]].values.flatten()  

    return(tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6,tempnp7)
    
def mksixgrayList(SubjRmean,Col,CnumList):#異なる代謝物グループ同士のリストの作成
    list1=[];list2=[];list3=[];list4=[];list5=[];list6=[];list7=[]
    listr1=[];listr2=[];listm1=[];listm2=[];listg1=[];listg2=[];listb1=[];listb2=[];listbl1=[];listpl1=[];listpl2=[];listbl2=[];listgr1=[];listgr2=[]
    listrabs1=[];listrabs2=[];listmabs1=[];listmabs2=[];listgabs1=[];listgabs2=[];listbabs1=[];listbabs2=[];listblabs1=[];listplabs1=[];listplabs2=[];listblabs2=[];listgrabs1=[];listgrabs2=[]
    listposr1=[];listposr2=[];listposm1=[];listposm2=[];listposg1=[];listposg2=[];listposb1=[];listposb2=[];listposbl1=[];listpospl1=[];listpospl2=[];listposbl2=[];listposgr1=[];listposgr2=[]
    for i in range(0,len(SubjRmean.columns)):#列を1~84まで
        if i < (CnumList[1]):
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[1]-1)]]);
            list1 += List;
            stloc = Col[CnumList[1]]
        elif i < (CnumList[2]):
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[2]-1)]]);
            list2 += List;            
            stloc = Col[CnumList[2]]
        elif i < (CnumList[3]):  
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[3]-1)]]);
            list3 += List;              
            stloc = Col[CnumList[3]]
        elif i < (CnumList[4]):
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[4]-1)]]);
            list4 += List;               
            stloc = Col[CnumList[4]]
        elif i < (CnumList[5]):            
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[5]-1)]]);
            list5 += List; 
            stloc = Col[CnumList[5]]
        elif i < (CnumList[6]):            
            List=list(SubjRmean.loc[Col[i],Col[i]:Col[(CnumList[6]-1)]]);
            list6 += List; 
            stloc = Col[CnumList[6]-1]
        List=list(SubjRmean.loc[Col[i],stloc:]);
        list7+=List
    return([list1,list2,list3,list4,list5,list6,list7])
    
def mksixabsList(SubjRmean,Col,CnumList):
    tempnp1 = abs(SubjRmean[Col[CnumList[0]]:Col[CnumList[1]-1]]).values.flatten(); 
    tempnp4 = abs(SubjRmean[Col[CnumList[3]]:Col[CnumList[4]-1]]).values.flatten(); 
    tempnp2 = abs(SubjRmean[Col[CnumList[1]]:Col[CnumList[2]-1]]).values.flatten(); 
    tempnp3 = abs(SubjRmean[Col[CnumList[2]]:Col[CnumList[3]-1]]).values.flatten();  
    tempnp5 = abs(SubjRmean[Col[CnumList[4]]:Col[CnumList[5]-1]]).values.flatten(); 
    tempnp6 = abs(SubjRmean[Col[CnumList[5]]:Col[CnumList[6]-1]]).values.flatten()  
    return(tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6)

def mksixList(SubjRmean,Col,CnumList):
    tempnp1 = SubjRmean[Col[CnumList[0]]:Col[CnumList[1]-1]].values.flatten(); 
    tempnp4 = SubjRmean[Col[CnumList[3]]:Col[CnumList[4]-1]].values.flatten(); 
    tempnp2 = SubjRmean[Col[CnumList[1]]:Col[CnumList[2]-1]].values.flatten(); 
    tempnp3 = SubjRmean[Col[CnumList[2]]:Col[CnumList[3]-1]].values.flatten();  
    tempnp5 = SubjRmean[Col[CnumList[4]]:Col[CnumList[5]-1]].values.flatten(); 
    tempnp6 = SubjRmean[Col[CnumList[5]]:Col[CnumList[6]-1]].values.flatten()  
    return(tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6)

def PlotMolList(Col,c,MolColor):#リスト内に対象分子のリスト入れてく
    alist = [[0] for i in range(len(c)+1)]
    for i in range(len(alist)):#MolColorの数を足しあげる
        if i==0:#最初はそのまま
            alist[i] = 0
        else:
            alist[i] = c[MolColor[i-1]] + alist[i-1]
    return(alist)
def plothistEachMetabo(SubjRmean,Col,c,MolColor,save_dir):#各代謝グループで色分けヒストグラム
        #リスト内に対象分子のリスト入れてく
            CnumList = PlotMolList(Col,c,MolColor)
            fig2, ax2 = plt.subplots(1,1)
            if len(CnumList) == 7:
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6 = mksixabsList(SubjRmean,Col,CnumList)
                ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=MolColor)        
                plt.savefig(save_dir + 'AbsHistAve.pdf')        
            #Stdもやる
                """
                fig3, ax3 = plt.subplots(1,1)
                tempDF = pd.DataFrame(data=None)
                tempnp1 = abs(SubjRstd[Col[0]:Col[(c['red'])]]).values.flatten();  tempnp4 = abs(SubjRstd[Col[green]:Col[blue]]).values.flatten(); tempnp2 = abs(SubjRstd[Col[(c['red'])]:Col[magenta]]).values.flatten();  tempnp3 = abs(SubjRstd[Col[magenta]:Col[green]]).values.flatten();   tempnp5 = abs(SubjRstd[Col[blue]:Col[purple]]).values.flatten(); tempnp6 = abs(SubjRstd[Col[purple]:Col[black]]).values.flatten()
                ax3 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
                plt.savefig(save_dir + 'AbsHistStd.pdf')  
                """
            #絶対値ではない           
                fig2, ax2 = plt.subplots(1,1)
                tempDF = pd.DataFrame(data=None)
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6 = mksixList(SubjRmean,Col,CnumList)
                #tempnp1 = SubjRmean[Col[0]:Col[(c['red'])]].values.flatten(); tempnp4 = SubjRmean[Col[green]:Col[blue]].values.flatten(); tempnp2 = SubjRmean[Col[(c['red'])]:Col[magenta]].values.flatten(); tempnp3 = SubjRmean[Col[magenta]:Col[green]].values.flatten();  tempnp5 = SubjRmean[Col[blue]:Col[purple]].values.flatten(); tempnp6 = SubjRmean[Col[purple]:Col[black]].values.flatten()
                ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=MolColor)        
                plt.savefig(save_dir + 'HistAve.pdf') 
                plt.close()
                """
            #Stdもやる
                fig3, ax3 = plt.subplots(1,1)
                tempDF = pd.DataFrame(data=None)
                tempnp1 = SubjRstd[Col[0]:Col[(c['red'])]].values.flatten();  tempnp4 = SubjRstd[Col[green]:Col[blue]].values.flatten(); tempnp2 = SubjRstd[Col[(c['red'])]:Col[magenta]].values.flatten();  tempnp3 = SubjRstd[Col[magenta]:Col[green]].values.flatten();  tempnp5 = SubjRstd[Col[blue]:Col[purple]].values.flatten();  tempnp6 = SubjRstd[Col[purple]:Col[black]].values.flatten()
                ax3 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
                plt.savefig(save_dir + 'HistStd.pdf') 
                plt.close()
                """
                """
             #正だけ＿その他付き
                fig4, ax4 = plt.subplots(1,1)
                ConListPos1 = [listposr1, listposm1, listposg1, listposb1, listpospl1, listposbl1, listposgr1];         ConList2 = [listr2, listm2, listg2, listb2, listpl2, listbl2, listgr2];
        
                ConListPos1 = [np.array(ConListPos1[x]) for x in range(len(ConListPos1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
        
                ax4 = plt.hist([ConListPos1[0][~np.isnan(ConListPos1[0])],ConListPos1[1][~np.isnan(ConListPos1[1])],ConListPos1[2][~np.isnan(ConListPos1[2])],ConListPos1[3][~np.isnan(ConListPos1[3])],ConListPos1[4][~np.isnan(ConListPos1[4])],ConListPos1[5][~np.isnan(ConListPos1[5])],ConListPos1[6][~np.isnan(ConListPos1[6])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black','gray'])        
                plt.savefig(save_dir + 'HistAve_Pos.pdf')  
                """
         
                #両方がその代謝系で、それ以外はグレー
                fig4, ax4 = plt.subplots(1,1)
                ConList1 = mksixgrayList(SubjRmean,Col,CnumList)#異なる代謝物グループ同士のリストの作成
                #ConList1 = [listr1, listm1, listg1, listb1, listpl1, listbl1, listgr1];         #ConList2 = [listr2, listm2, listg2, listb2, listpl2, listbl2, listgr2];
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         #ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
                MolColorGray = MolColor+['gray']
                ax4 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=MolColorGray)        
                plt.savefig(save_dir + 'HistAveSameMol.pdf')  
    
                #fig5, ax5 = plt.subplots(1,1)
                #ax4 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])],ConList2[6][~np.isnan(ConList2[6])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black','gray'])        
        
                #ax5 = plt.hist([listr2, listm2, listg2, listb2, listbl2, listgr2], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'])        
                #plt.savefig(save_dir + 'HistStdSameMol.pdf')  
                fig6, ax6 = plt.subplots(1,1)
                ax6 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])]], stacked=True, bins = 20,color=MolColor)        
                plt.savefig(save_dir + 'HistAveSameMolWOGray.pdf') 
        
                #fig7, ax7 = plt.subplots(1,1)
                #ax7 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
                #plt.savefig(save_dir + 'HistStdSameMolWOGray.pdf')  ;print(ConList2[4][~np.isnan(ConList2[4])])
            
                fig8, ax8 = plt.subplots(1,1)
                #[[abs(ablist[j][i]) for i in range(len(alist))] for j in range(len(ablist))]
                ConList1 = [[abs(ConList1[j][i]) for i in range(len(ConList1[j]))] for j in range(len(ConList1))]
                
                #ConList1 = [listrabs1, listmabs1, listgabs1, listbabs1, listplabs1, listblabs1, listgrabs1];         ConList2 = [listrabs2, listmabs2, listgabs2, listbabs2,listplabs2, listblabs2, listgrabs2];
        
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         #ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
        
                ax8 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=MolColorGray)        
                plt.savefig(save_dir + 'AbsHistAveSameMol.pdf') 
    
                plt.close()
                
                fig10, ax10 = plt.subplots(1,1)
                ax10 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])]], stacked=True, bins = 20,color=MolColor)        
    
                xmin, xmax, ymin, ymax = plt.axis()
                #plt.plot( [Upper5,Upper5],[ymin, ymax], linestyle="dashed" ,color='black', linewidth=1)
                plt.ylim([ymin,ymax])
                plt.savefig(save_dir + 'AbsHistAveSameMolWOGray.pdf')             
                """
                ########3#正だけ_同じ代謝物
                fig4, ax4 = plt.subplots(1,1)
                ConListPos1 = [listposr1, listposm1, listposg1, listposb1, listpospl1, listposbl1, listposgr1];         ConList2 = [listr2, listm2, listg2, listb2, listpl2, listbl2, listgr2];
        
                ConListPos1 = [np.array(ConListPos1[x]) for x in range(len(ConListPos1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
        
                ax4 = plt.hist([ConListPos1[0][~np.isnan(ConListPos1[0])],ConListPos1[1][~np.isnan(ConListPos1[1])],ConListPos1[2][~np.isnan(ConListPos1[2])],ConListPos1[3][~np.isnan(ConListPos1[3])],ConListPos1[4][~np.isnan(ConListPos1[4])],ConListPos1[5][~np.isnan(ConListPos1[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
                plt.savefig(save_dir + 'HistAveSaneMol_Pos.pdf') 
                """
                #################5%の比較#################5%の比較#################5%の比較#################5%の比較#################5%の比較
                Abs = list(np.abs(SubjRmean.values.flatten())[~np.isnan(SubjRmean.values.flatten())])
                print( '5%パーセンタイル:'+str( np.percentile(Abs,95) ))
                Upper5 = SC.select_sorted_num(Abs, 0.05)
                print( '上位5%に当たる値:' + str(Upper5) )
                xmin, xmax, ymin, ymax = plt.axis()
                #plt.plot( [Upper5,Upper5],[ymin, ymax], linestyle="dashed" ,color='black', linewidth=1)
                plt.ylim([ymin,ymax])
                plt.yticks([0, 100, 200, 300, 400, 500])
                #
                ax8 = plt.gca()
                ax8.spines['bottom'].set_linewidth(1.5); ax8.spines['left'].set_linewidth(1.5); ax8.spines['right'].set_linewidth(1.5); ax8.spines['top'].set_linewidth(1.5); 
                #plt.savefig(save_dir + 'AbsHistAveSameMol.pdf')  
        
                #fig9, ax9 = plt.subplots(1,1)
                #ax9 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])],ConList2[6][~np.isnan(ConList2[6])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black','gray'])        
                #ax5 = plt.hist([listr2, listm2, listg2, listb2, listbl2, listgr2], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'])        
                #plt.savefig(save_dir + 'AbsHistStdSameMol.pdf')
                #plt.rcParams['axes.linewidth'] = 1.5# 軸の線幅edge linewidth。囲みの太さ
                """
                ConList1 = [listrabs1, listmabs1, listgabs1, listbabs1, listplabs1, listblabs1, listgrabs1];         ConList2 = [listrabs2, listmabs2, listgabs2, listbabs2,listplabs2, listblabs2, listgrabs2];
        
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
        
                
                fig10, ax10 = plt.subplots(1,1)
                ax10 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
    
                xmin, xmax, ymin, ymax = plt.axis()
                #plt.plot( [Upper5,Upper5],[ymin, ymax], linestyle="dashed" ,color='black', linewidth=1)
                plt.ylim([ymin,ymax])
                plt.savefig(save_dir + 'AbsHistAveSameMolWOGray.pdf') 
        
                #fig11, ax11 = plt.subplots(1,1)
                #ax11 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black'])        
                #plt.savefig(save_dir + 'AbsHistStdSameMolWOGray.pdf') 
                """
            elif len(CnumList) == 8:
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6,tempnp7 = mksevenabsList(SubjRmean,Col,CnumList)
                ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)]),list(tempnp7[~np.isnan(tempnp7)])], stacked=True, bins = 20,color=MolColor)        
                plt.savefig(save_dir + 'AbsHistAve.pdf')   
            #絶対値ではない           
                fig2, ax2 = plt.subplots(1,1)
                tempDF = pd.DataFrame(data=None)
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6,tempnp7 = mksevenList(SubjRmean,Col,CnumList)
                ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)]),list(tempnp7[~np.isnan(tempnp7)])], stacked=True, bins = 20,color=MolColor)        
                plt.savefig(save_dir + 'HistAve.pdf') 
                plt.close()                
                #両方がその代謝系で、それ以外はグレー
                fig4, ax4 = plt.subplots(1,1)
                ConList1 = mksevengrayList(SubjRmean,Col,CnumList)#異なる代謝物グループ同士のリストの作成
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         #ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
                MolColorGray = MolColor+['gray']
                ax4 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])],ConList1[7][~np.isnan(ConList1[7])]], stacked=True, bins = 20,color=MolColorGray)        
                plt.savefig(save_dir + 'HistAveSameMol.pdf')  

                fig6, ax6 = plt.subplots(1,1)
                ax6 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=MolColor)        
                plt.savefig(save_dir + 'HistAveSameMolWOGray.pdf') 
            
                fig8, ax8 = plt.subplots(1,1)
                ConList1 = [[abs(ConList1[j][i]) for i in range(len(ConList1[j]))] for j in range(len(ConList1))]
                
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         #ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
        
                ax8 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])],ConList1[7][~np.isnan(ConList1[7])]], stacked=True, bins = 20,color=MolColorGray)        
                plt.savefig(save_dir + 'AbsHistAveSameMol.pdf') 
    
                plt.close()
                
                fig10, ax10 = plt.subplots(1,1)
                ax10 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=MolColor)        
    
                xmin, xmax, ymin, ymax = plt.axis()
                plt.ylim([ymin,ymax])
                plt.savefig(save_dir + 'AbsHistAveSameMolWOGray.pdf')   
def plotculmitivehistEachMetabo(SubjRmean,Col,c,MolColor,save_dir):#各代謝グループで色分け累積ヒストグラム
            msize=3; lwidth=1.5
        #リスト内に対象分子のリスト入れてく
            CnumList = PlotMolList(Col,c,MolColor)
            if len(CnumList) == 7:
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6 = mksixabsList(SubjRmean,Col,CnumList)
                ConList1 = mksixgrayList(SubjRmean,Col,CnumList)#異なる代謝物グループ同士のリストの作成

                ConList1 = [[abs(ConList1[j][i]) for i in range(len(ConList1[j]))] for j in range(len(ConList1))]
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];
                fig7, ax7 = plt.subplots(1,1)
                x = np.sort(ConList1[0][~np.isnan(ConList1[0])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0)
                ### temp_20201104
                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)
                
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[0])
                x = np.sort(ConList1[1][~np.isnan(ConList1[1])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0);
                ### temp_20201104
                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)                

                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[1])
    
                x=np.insert(x,0,0);y=np.insert(y,0,0)
                
                x = np.sort(ConList1[2][~np.isnan(ConList1[2])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0);
                ### temp_20201104
                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[2])
    
                x = np.sort(ConList1[3][~np.isnan(ConList1[3])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0);
                ### temp_20201104
                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)                
                #y=np.insert(y,0,0)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[3])
    
                x = np.sort(ConList1[4][~np.isnan(ConList1[4])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0);
                                ### temp_20201104
                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)
                #y=np.insert(y,0,0)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[4])
                
                x = np.sort(ConList1[5][~np.isnan(ConList1[5])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0);
                                ### temp_20201104
                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)
                #y=np.insert(y,0,0)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[5])
                
                x = np.sort(ConList1[6][~np.isnan(ConList1[6])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0);
                                ### temp_20201104
                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)
                #y=np.insert(y,0,0)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c='gray')
            
            
            elif len(CnumList) == 8:
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6,tempnp7 = mksevenabsList(SubjRmean,Col,CnumList)            
                ConList1 = mksevengrayList(SubjRmean,Col,CnumList)#異なる代謝物グループ同士のリストの作成
                ConList1 = [[abs(ConList1[j][i]) for i in range(len(ConList1[j]))] for j in range(len(ConList1))]
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];

                fig7, ax7 = plt.subplots(1,1)
                x = np.sort(ConList1[0][~np.isnan(ConList1[0])]); y=np.arange(1,len(x)+1) / len(x)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[0])
                x = np.sort(ConList1[1][~np.isnan(ConList1[1])]); y=np.arange(1,len(x)+1) / len(x)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[1])
    
                x = np.sort(ConList1[2][~np.isnan(ConList1[2])]); y=np.arange(1,len(x)+1) / len(x)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[2])
    
                x = np.sort(ConList1[3][~np.isnan(ConList1[3])]); y=np.arange(1,len(x)+1) / len(x)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[3])
    
                x = np.sort(ConList1[4][~np.isnan(ConList1[4])]); y=np.arange(1,len(x)+1) / len(x)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[4])
                x = np.sort(ConList1[5][~np.isnan(ConList1[5])]); y=np.arange(1,len(x)+1) / len(x)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[5])

                x = np.sort(ConList1[6][~np.isnan(ConList1[6])]); y=np.arange(1,len(x)+1) / len(x)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[6])
                x = np.sort(ConList1[7][~np.isnan(ConList1[7])]); y=np.arange(1,len(x)+1) / len(x)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c='gray')            

            
def PlotScatter(save_dir,SubjRmean,SubjRstd,ColorSwitch,NormlSwitch,SwitchDict):
    Ave = list(SubjRmean.values.flatten()[~np.isnan(SubjRmean.values.flatten().astype(float))]);        AbsAve = list(abs(SubjRmean.values.flatten()[~np.isnan(SubjRmean.values.flatten().astype(float))]))
    AbsAvemin = min(AbsAve); AbsAvemax = max(AbsAve)
    Std = list(SubjRstd.values.flatten()[~np.isnan(SubjRstd.values.flatten().astype(float))])
    Stdmin = min(Std); Stdmax = max(Std)
    plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
    plt.rcParams["font.size"] = 20
#    Ave.remove(1);Std.remove(0);#Std.remove(1)

    r, p=pearsonr(Ave,Std)
    p='{:e}'.format(p)

    TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
    if ColorSwitch == 0:

        fig, ax = plt.subplots(1,1)
        ax.scatter(Ave,Std)
        SaveName = []
        SaveName='_WNmlColor'
        fig.savefig(save_dir +'AveStdScatter' + SaveName +'.pdf',format='pdf')
    
    elif ColorSwitch == 1:#全部青色正負考慮
        MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180306/MolColorDF.xlsx',header=0,encoding = "ISO-8859-1")
        fig, ax = plt.subplots(1,1)
        ColorList =[]# MolColorDF['MolColor']
        Col = SubjRmean.columns
        for i in range(0,len(SubjRmean.columns)):
            list1=list(SubjRmean.loc[Col[i],Col[i]:]);list2 =list(SubjRstd.loc[Col[i],Col[i]:])
            list1.remove(1);list2.remove(0);#Std.remove(1)
            
            ax.scatter(list1,list2,color='Blue')#ColorList[i])        
        SaveName='_WBlueColor'
        fig.savefig(save_dir +'AveStdScatter' + SaveName +'.pdf',format='pdf')
        
    elif ColorSwitch == 2:#全部青色絶対値
        MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180306/MolColorDF.xlsx',header=0,encoding = "ISO-8859-1")
        fig, ax = plt.subplots(1,1)
        ColorList = []#MolColorDF['MolColor']
        Col = SubjRmean.columns
        for i in range(0,len(SubjRmean.columns)):
            list1=list(abs(SubjRmean.loc[Col[i],Col[i]:]));list2 =list(SubjRstd.loc[Col[i],Col[i]:])
            list1.remove(1);list2.remove(0);#Std.remove(1)
            
            ax.scatter(list1,list2,color='blue')#ColorList[i]) 
        SaveName='_WBlueColor'
        fig.savefig(save_dir +'abs(Ave)StdScatter' + SaveName +'.pdf',format='pdf')
            
    elif ColorSwitch == 3:#同じ代謝の散布図に色をつける
        MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/MolColor.xlsx',header=0) # 日本語
        #MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolColor_Eng.xlsx',header=0,encoding = "ISO-8859-1",index_col=0) # 日本語
        MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,index_col=0)
        
        fig1, ax1 = plt.subplots(1,1)#abs
        fig2, ax2 = plt.subplots(1,1)#normal

        Col = list(SubjRmean.columns)
        #Col = MolColorDF.index
        try:
            MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0)
            if SwitchDict['Check'] == 'Amino':#AAのラベルなど変える
                MolColorDF,a = LH.AACheck(MolColor,MolColor,'EAA')#'protein','ketogenic','EAA','SemiEAA'
            elif SwitchDict['Check'] == 'Glc':#糖代謝系のラベルなど変える
                MolColorDF,a = LH.TCACheck(MolColor,MolColor,'TCA')#'TCA'
        except:
            pass
        ColorList = MolColorDF['MolColor']
        ColorList = MolColorDF.loc[Col,'MolColor']
        c = collections.Counter(ColorList);listr1=[];listr2=[];listm1=[];listm2=[];listg1=[];listg2=[];listb1=[];listb2=[];listbl1=[];listpl1=[];listpl2=[];listbl2=[];listgr1=[];listgr2=[]
        listrabs1=[];listrabs2=[];listmabs1=[];listmabs2=[];listgabs1=[];listgabs2=[];listbabs1=[];listbabs2=[];listblabs1=[];listplabs1=[];listplabs2=[];listblabs2=[];listgrabs1=[];listgrabs2=[]
        listposr1=[];listposr2=[];listposm1=[];listposm2=[];listposg1=[];listposg2=[];listposb1=[];listposb2=[];listposbl1=[];listpospl1=[];listpospl2=[];listposbl2=[];listposgr1=[];listposgr2=[]
        if SwitchDict['MolColorScatter']  == 1:#代謝ごとの散布図をだす
            for i in range(0,len(SubjRmean.columns)):#列を1~84まで
                if i < (c['red']):
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[(c['red']-1)]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(c['red']-1)]])
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[(c['red']-1)]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(c['red']-1)]])
                    redNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[(c['red']-1)]]);
                    listpos1=list(redNP[redNP>0]);listpos2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(c['red']-1)]])
   
                    #list1.remove(1);list2.remove(0);#Std.remove(1)
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    
                    ax2.scatter(list1,list2,color=ColorList[i])
    
                    listrabs1 += listabs1; listrabs2 += listabs2;  listr1 += list1; listr2 += list2
                    listposr1 += listpos1; listposr2 += listpos2
                    #GH.ScatterWHeatmap(np.array(listrabs1)[~np.isnan(listrabs1)],np.array(listrabs2)[~np.isnan(listrabs2)],AbsAvemin,AbsAvemax,Stdmin,Stdmax,save_dir+'_red',bin=30)#散布図をヒートマップで
                    #GH.ScatterWColorbar(np.array(listrabs1)[~np.isnan(listrabs1)],np.array(listrabs2)[~np.isnan(listrabs2)],0,save_dir+'red')
                    #カーネル密度推定して色分け
                    plt.close()
                    
                    stloc = Col[(c['red'])]
                    magenta = c['red']+c['magenta']
                elif i < magenta:
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[magenta]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(magenta)]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[magenta]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(magenta)]])
                    magentaNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[magenta]]);
                    listpos1=list(magentaNP[magentaNP>0]);listpos2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(magenta)]])
    
                   #list1.remove(1);list2.remove(0);#Std.remove(1)
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i])
    
                    listm1 += list1; listm2 += list2; listmabs1 += listabs1; listmabs2 += listabs2;
                    listposm1 += listpos1; listposm2 += listpos2;
                    #GH.ScatterWHeatmap(np.array(listmabs1)[~np.isnan(listmabs1)],np.array(listmabs2)[~np.isnan(listmabs2)],AbsAvemin,AbsAvemax,Stdmin,Stdmax,save_dir+'_magenta',bin=30)#散布図をヒートマップで
                    #GH.ScatterWColorbar(np.array(listmabs1)[~np.isnan(listmabs1)],np.array(listmabs2)[~np.isnan(listmabs2)],0,save_dir+'magenta')
    
                    plt.close()
                    stloc = Col[magenta]
                    green = magenta + c['green']
                    #print(Col[i],ColorList[i])
                elif i < green:
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[green]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[green]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
                    greenNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[green]])
                    listpos1=list(greenNP[greenNP>0]);#list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
    
                    #list1.remove(1);list2.remove(0);#Std.remove(1)
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i]) 
                    listg1 += list1; listg2 += list2;  listgabs1 += listabs1; listgabs2 += listabs2
                    listposg1 += listpos1; listg2 += list2; 
                    #GH.ScatterWHeatmap(np.array(listgabs1)[~np.isnan(listgabs1)],np.array(listgabs2)[~np.isnan(listgabs2)],AbsAvemin,AbsAvemax,Stdmin,Stdmax,save_dir+'_green',bin=30)#散布図をヒートマップで
                    #GH.ScatterWColorbar(np.array(listgabs1)[~np.isnan(listgabs1)],np.array(listgabs2)[~np.isnan(listgabs2)],0,save_dir+'green')
                    
                    plt.close()
                    
                    stloc = Col[green]
                    blue = green + c['blue']
                    #print(Col[i],ColorList[i])
                elif i < (blue):
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[blue]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[blue]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[blue]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[blue]])
                    blueNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[blue]])
                    listpos1=list(blueNP[blueNP>0]);#list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
    
                    #list1.remove(1);list2.remove(0);#Std.remove(1)
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i]) 
                    listb1 += list1; listb2 += list2;  listbabs1 += listabs1; listbabs2 += listabs2
                    listposb1 += listpos1; listb2 += list2; 
                    #GH.ScatterWHeatmap(np.array(listbabs1)[~np.isnan(listbabs1)],np.array(listbabs2)[~np.isnan(listbabs2)],AbsAvemin,AbsAvemax,Stdmin,Stdmax,save_dir+'_blue',bin=30)#散布図をヒートマップで
                    #GH.ScatterWColorbar(np.array(listbabs1)[~np.isnan(listbabs1)],np.array(listbabs2)[~np.isnan(listbabs2)],0,save_dir+'blue')
                    
                    plt.close()                
                    stloc = Col[blue]
                    purple = blue + c['purple']
                    #print(Col[i],ColorList[i])
                elif i < (purple):
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[purple]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[purple]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[purple]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[purple]])
                    purpleNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[purple]])
                    listpos1=list(purpleNP[purpleNP>0]);#list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
    
                    #list1.remove(1);list2.remove(0);#Std.remove(1)
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i]) 
                    listpl1 += list1; listpl2 += list2;  listplabs1 += listabs1; listplabs2 += listabs2
                    listpospl1 += listpos1; listpl2 += list2;
                    #GH.ScatterWHeatmap(np.array(listbabs1)[~np.isnan(listbabs1)],np.array(listbabs2)[~np.isnan(listbabs2)],AbsAvemin,AbsAvemax,Stdmin,Stdmax,save_dir+'_blue',bin=30)#散布図をヒートマップで
                    #GH.ScatterWColorbar(np.array(listplabs1)[~np.isnan(listplabs1)],np.array(listplabs2)[~np.isnan(listplabs2)],0,save_dir+'purple')
                    
                    plt.close()                
                    stloc = Col[purple]
                    black = purple + c['black']-1
                    #print(Col[i],ColorList[i])
                elif i < (black):
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[black]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[black]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[black]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[black]])
                    blackNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[black]])
                    listpos1=list(blackNP[blackNP>0]);#list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
    
                    #list1.remove(1);list2.remove(0);#Std.remove(1)
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i])  
                    listbl1 += list1; listbl2 += list2;  listblabs1 += listabs1; listblabs2 += listabs2
                    listposbl1 += listpos1; listbl2 += list2; 
                    #GH.ScatterWHeatmap(np.array(listblabs1)[~np.isnan(listblabs1)],np.array(listblabs2)[~np.isnan(listblabs2)],AbsAvemin,AbsAvemax,Stdmin,Stdmax,save_dir+'_black',bin=30)#散布図をヒートマップで
                    #GH.ScatterWColorbar(np.array(listblabs1)[~np.isnan(listblabs1)],np.array(listblabs2)[~np.isnan(listblabs2)],0,save_dir+'black')
     
                    plt.close()
                    stloc = Col[black]
                list1=list(SubjRmean.loc[Col[i],stloc:]);list2 =list(SubjRstd.loc[Col[i],stloc:])
                listabs1=list(abs(SubjRmean.loc[Col[i],stloc:]));listabs2 =list(SubjRstd.loc[Col[i],stloc:])
                greyNP = np.array(SubjRmean.loc[Col[i],stloc:])
                listpos1=list(greyNP[greyNP>0]);#list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
    
                listgr1 += list1;listgr2 += list2;    listgrabs1 += listabs1;listgrabs2 += listabs2
                listposgr1 += listpos1;listgr2 += list2;
                #GH.ScatterWColorbar(np.array(listgrabs1)[~np.isnan(listgrabs1)],np.array(listgrabs2)[~np.isnan(listgrabs2)],0,save_dir+'gray')
    
                #listgr1.remove(1);#listgr2.remove(1);#Std.remove(1)
                
                ax1.scatter(listabs1,listabs2,color='gray',alpha=0.5)
                ax2.scatter(list1,list2,color='gray',alpha=0.5)  
                
            Std = list(SubjRstd.values.flatten()[~np.isnan(SubjRstd.values.flatten())])
            absr, absp=pearsonr(AbsAve,Std)
            p='{:e}'.format(absp)
            SaveName='_WMolColor_Abs'
                
            ax1.set_title('R=' + str(round(absr,3)) + ', p=' + p,fontsize=10);ax2.set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp,fontsize=10)
            ax1.set_xlabel('Abs(Ave)',fontsize=10);ax2.set_xlabel('Ave',fontsize=10)              
            ax1.set_ylabel('Std',fontsize=10);ax2.set_ylabel('Std',fontsize=10)#,fontproperties=prop)     fig.tight_layout()
            plot_axis = plt.axis()
            fig1.savefig(save_dir +'AbsAveStdScatter' + SaveName +'.pdf',format='pdf')
            plt.close()
            fig2.savefig(save_dir +'AveStdScatter' + SaveName +'.pdf',format='pdf')
            plt.close()       
            
        if NormlSwitch == 0:#とりあえず簡単に片方の代謝系で色分け__ヒストグラム描画
            plothistEachMetabo(SubjRmean,Col,c,list(c),save_dir)#各代謝グループで色分けヒストグラム
        elif NormlSwitch == 1:#各代謝ごとに総数で正規化して累積分布
            al=1
            hype='bar'
            FileNameTag = 'Nmlz'
            T=True;F=False
            
            """
            fig2, ax2 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = abs(SubjRmean[Col[0]:Col[(c['red'])]]).values.flatten(); tempnp2 = abs(SubjRmean[Col[(c['red'])]:Col[magenta]]).values.flatten(); tempnp3 = abs(SubjRmean[Col[magenta]:Col[green]]).values.flatten();   tempnp4 = abs(SubjRmean[Col[green]:Col[blue]]).values.flatten();tempnp5 = abs(SubjRmean[Col[blue]:Col[purple]]).values.flatten();tempnp6 = abs(SubjRmean[Col[purple]:Col[black]]).values.flatten()
            n, bins, patches = plt.hist(list(tempnp1[~np.isnan(tempnp1)]),  normed=False,cumulative=False,bins = 20,color='red',alpha=al,histtype=hype)  
            # 第2軸用値の算出
            """
        #Stdもやる
        #しばらくいらない20190514
            """   
            fig3, ax3 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = SubjRstd[Col[0]:Col[(c['red'])]].values.flatten();  tempnp4 = SubjRstd[Col[green]:Col[blue]].values.flatten(); tempnp2 = SubjRstd[Col[(c['red'])]:Col[magenta]].values.flatten();  tempnp3 = SubjRstd[Col[magenta]:Col[green]].values.flatten();   tempnp5 = SubjRstd[Col[blue]:Col[black]].values.flatten()
            ax3 = plt.hist(list(tempnp1[~np.isnan(tempnp1)]),  normed=True,cumulative=True,bins = 20,color='red',alpha=al,histtype=hype)        
            ax3 = plt.hist(list(tempnp2[~np.isnan(tempnp2)]),  normed=True,cumulative=True,bins = 20,color='magenta',alpha=al,histtype=hype)        
            ax3 = plt.hist(list(tempnp3[~np.isnan(tempnp3)]),  normed=True,cumulative=True,bins = 20,color='green',alpha=al,histtype=hype)        
            ax3 = plt.hist(list(tempnp4[~np.isnan(tempnp4)]),  normed=True,cumulative=True,bins = 20,color='blue',alpha=al,histtype=hype)        
            ax3 = plt.hist(list(tempnp5[~np.isnan(tempnp5)]),  normed=True,cumulative=True,bins = 20,color='purple',alpha=al,histtype=hype) 
            ax3 = plt.hist(list(tempnp6[~np.isnan(tempnp6)]),  normed=True,cumulative=True,bins = 20,color='black',alpha=al,histtype=hype) 
            #plt.savefig(save_dir + 'HistStd'+FileNameTag+'.pdf') 
            plt.close();#ax3.axis('off')
            
            #両方がその代謝系で、それ以外はグレー
            fig4, ax4 = plt.subplots(1,1)
            ConList1 = [listr1, listm1, listg1, listb1, listpl1, listbl1, listgr1];         ConList2 = [listr2, listm2, listg2, listb2, listpl2, listbl2, listgr2];
    
            ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
            """
            """    
            ax4 = plt.hist(ConList1[0][~np.isnan(ConList1[0])],  normed=True,cumulative=True,bins = 20,color='red',alpha=al,histtype=hype)        
            ax4 = plt.hist(ConList1[1][~np.isnan(ConList1[1])],  normed=True,cumulative=True,bins = 20,color='magenta',alpha=al,histtype=hype)        
            ax4 = plt.hist(ConList1[2][~np.isnan(ConList1[2])],  normed=True,cumulative=True,bins = 20,color='green',alpha=al,histtype=hype)        
            ax4 = plt.hist(ConList1[3][~np.isnan(ConList1[3])],  normed=True,cumulative=True,bins = 20,color='blue',alpha=al,histtype=hype)        
            ax4 = plt.hist(ConList1[4][~np.isnan(ConList1[4])],  normed=True,cumulative=True,bins = 20,color='purple',alpha=al,histtype=hype)        
            ax4 = plt.hist(ConList1[5][~np.isnan(ConList1[5])],  normed=True,cumulative=True,bins = 20,color='black',alpha=al,histtype=hype)        
            ax4 = plt.hist(ConList1[6][~np.isnan(ConList1[6])],  normed=True,cumulative=True,bins = 20,color='gray',alpha=al,histtype=hype)        

            plt.savefig(save_dir + 'HistAveSameMol'+FileNameTag+'.pdf')
            plt.close();#ax4.axis('off')
    
            fig5, ax5 = plt.subplots(1,1)
            ax5 = plt.hist(ConList2[0][~np.isnan(ConList2[0])],  normed=True,cumulative=True,bins = 20,color='red',alpha=0.3,histtype=hype)        
            ax5 = plt.hist(ConList2[1][~np.isnan(ConList2[1])],  normed=True,cumulative=True,bins = 20,color='magenta',alpha=0.3,histtype=hype)        
            ax5 = plt.hist(ConList2[2][~np.isnan(ConList2[2])],  normed=True,cumulative=True,bins = 20,color='green',alpha=0.3,histtype=hype)        
            ax5 = plt.hist(ConList2[3][~np.isnan(ConList2[3])],  normed=True,cumulative=True,bins = 20,color='blue',alpha=0.3,histtype=hype)        
            ax5 = plt.hist(ConList2[4][~np.isnan(ConList2[4])],  normed=True,cumulative=True,bins = 20,color='purple',alpha=0.3,histtype=hype) 
            ax5 = plt.hist(ConList2[5][~np.isnan(ConList2[5])],  normed=True,cumulative=True,bins = 20,color='black',alpha=0.3,histtype=hype) 
            ax5 = plt.hist(ConList2[6][~np.isnan(ConList2[6])],  normed=True,cumulative=True,bins = 20,color='gray',alpha=al,histtype=hype)        
    
            #ax5 = plt.hist([listr2, listm2, listg2, listb2, listbl2, listgr2], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'])        
            plt.savefig(save_dir + 'HistStdSameMol'+FileNameTag+'.pdf')
            plt.close();#ax5.axis('off')
    
            fig6, ax6 = plt.subplots(1,1)
            ax6 = plt.hist(ConList1[0][~np.isnan(ConList1[0])],  normed=True,cumulative=True,bins = 20,color='red',alpha=al,histtype=hype)        
            ax6 = plt.hist(ConList1[1][~np.isnan(ConList1[1])],  normed=True,cumulative=True,bins = 20,color='magenta',alpha=al,histtype=hype)        
            ax6 = plt.hist(ConList1[2][~np.isnan(ConList1[2])],  normed=True,cumulative=True,bins = 20,color='green',alpha=al,histtype=hype)        
            ax6 = plt.hist(ConList1[3][~np.isnan(ConList1[3])],  normed=True,cumulative=True,bins = 20,color='blue',alpha=al,histtype=hype)  
            ax6 = plt.hist(ConList1[4][~np.isnan(ConList1[4])],  normed=True,cumulative=True,bins = 20,color='purple',alpha=al,histtype=hype)  
            ax6 = plt.hist(ConList1[5][~np.isnan(ConList1[5])],  normed=True,cumulative=True,bins = 20,color='black',alpha=al,histtype=hype)  

            plt.savefig(save_dir + 'HistAveSameMolWOGray'+FileNameTag+'.pdf')
            plt.close();#ax6.axis('off')
            """   
            """#しばらくいらない「20190514
            hype='step'
            CList=['red','magenta','green','blue','purple','black','gray']

            fig7, ax7 = plt.subplots(1,1)
            x = np.sort(ConList1[0][~np.isnan(ConList1[0])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[0])
            x = np.sort(ConList1[1][~np.isnan(ConList1[1])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[1])

            x = np.sort(ConList1[2][~np.isnan(ConList1[2])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[2])

            x = np.sort(ConList1[3][~np.isnan(ConList1[3])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[3])

            x = np.sort(ConList1[4][~np.isnan(ConList1[4])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[4])
            x = np.sort(ConList1[5][~np.isnan(ConList1[5])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[5])

            #下はSTD
            #ax7 = plt.hist(,  normed=True,cumulative=True,bins = len(ConList2[0][~np.isnan(ConList2[0])]),color='red',alpha=al,histtype=hype)        
            #ax7 = plt.hist(ConList2[1][~np.isnan(ConList2[1])],  normed=True,cumulative=True,bins = 20,color='magenta',alpha=al,histtype=hype)        
            #ax7 = plt.hist(ConList2[2][~np.isnan(ConList2[2])],  normed=True,cumulative=True,bins = 20,color='green',alpha=al,histtype=hype)        
            #ax7 = plt.hist(ConList2[3][~np.isnan(ConList2[3])],  normed=True,cumulative=True,bins = 20,color='blue',alpha=al,histtype=hype)  
            ##ax7 = plt.hist(ConList2[4][~np.isnan(ConList2[4])],  normed=True,cumulative=True,bins = 20,color='purple',alpha=al,histtype=hype)  
            #ax7 = plt.hist(ConList2[5][~np.isnan(ConList2[5])],  normed=True,cumulative=True,bins = 20,color='black',alpha=al,histtype=hype)  

            plt.savefig(save_dir + 'HistStdSameMolWOGray'+FileNameTag+'.pdf')  
            plt.close();#ax7.axis('off')
            
            fig8, ax7 = plt.subplots(1,1)
            fig9, ax8 = plt.subplots(1,1)
            ax8.set_ylim([-0.05,1.05]);
            ax7.set_yticks(np.linspace(0,1.0,0.2))
            """
            ########################################累積分布線分；マーカーのサイズ
            plotculmitivehistEachMetabo(SubjRmean,Col,c,list(c),save_dir)#各代謝グループで色分け累積ヒストグラム
            CList=['red','magenta','green','blue','purple','black','gray']

            msize=3; lwidth=1.5
            ConList1 = [listrabs1, listmabs1, listgabs1, listbabs1,listplabs1, listblabs1, listgrabs1];         ConList2 = [listrabs2, listmabs2, listgabs2, listbabs2, listplabs2, listblabs2, listgrabs2];

            ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
            print(ConList1)
            #とりあえず削除20190514
            """
            fig7, ax7 = plt.subplots(1,1)
            x = np.sort(ConList1[0][~np.isnan(ConList1[0])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[0])
            x = np.sort(ConList1[1][~np.isnan(ConList1[1])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[1])

            x = np.sort(ConList1[2][~np.isnan(ConList1[2])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[2])

            x = np.sort(ConList1[3][~np.isnan(ConList1[3])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[3])

            x = np.sort(ConList1[4][~np.isnan(ConList1[4])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[4])
            x = np.sort(ConList1[5][~np.isnan(ConList1[5])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[5])
            x = np.sort(ConList1[6][~np.isnan(ConList1[6])]); y=np.arange(1,len(x)+1) / len(x)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[6])
            """
            """

            n, binsAll, patches =ax7.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])],ConList2[6][~np.isnan(ConList2[6])]], stacked=True, bins = 20,color=['red','magenta','green','blue','purple','black','gray'])
            #糖代謝の各binに含まれる数を足しあげる。
            n, bins, patches =ax7.hist(ConList1[0][~np.isnan(ConList1[0])],  normed=F,cumulative=F,bins = 20,color='red',alpha=al,histtype=hype)        
            #nがはじめ続くところと終わり流づくところはplot, 中のところは Idxとってきてplotしなくする

            yr = np.insert( np.add.accumulate(n) / n.sum() , -1, 1) ;xr = np.arange(0,1.05,0.05)#binsAll#np.convolve(bins, np.ones(2) / 2, mode="same")[1:]

            #yr = CumsumCountList(ConList1[0][~np.isnan(ConList1[0])])
            yr, xr = perform4plot(yr, np.arange(0,1.05,0.05))
            # 第2軸のプロット
            ax23 = ax8.twinx();lines = ax23.plot(xr, yr, ls='-', color='red', marker='o',label='Cumulative ratio',markersize=msize,linewidth = lwidth)
            hold(True)
                     
            ax23.set_yticks(np.linspace(0,1.0,0.2))
            #ホルモンの各binに含まれる数を足しあげる。
            n, bins, patches =ax7.hist(ConList1[1][~np.isnan(ConList1[1])],  normed=F,cumulative=F,bins = 20,color='magenta',alpha=al,histtype=hype)        
            ym = np.add.accumulate(n) / n.sum();xm = np.arange(0,1.05,0.05)#np.convolve(bins, np.ones(2) / 2, mode="same")[1:]

            ym = CumsumCountList(ConList1[1][~np.isnan(ConList1[1])])
            ym, xm = perform4plot(ym, np.arange(0,1.05,0.05))
            ax23 = ax8.twinx();lines = ax23.plot(xm, ym, ls='-', color='magenta', marker='o',label='Cumulative ratio',markersize=msize,linewidth = lwidth)
            hold(True)
            
            ax23.set_yticks(np.linspace(0,1.0,0.2))
            #脂質の各binに含まれる数を足しあげる。
            n, bins, patches =ax7.hist(ConList1[2][~np.isnan(ConList1[2])],  normed=F,cumulative=F,bins = 20,color='green',alpha=al,histtype=hype)        
            yg = np.add.accumulate(n) / n.sum();xg = np.arange(0,1.05,0.05)#np.convolve(bins, np.ones(2) / 2, mode="same")[1:]

            yg = CumsumCountList(ConList1[2][~np.isnan(ConList1[2])])
            yg, xg = perform4plot(yg, np.arange(0,1.05,0.05))
            ax23 = ax8.twinx();lines = ax23.plot(xg, yg, ls='-', color='green', marker='o',label='Cumulative ratio',markersize=msize,linewidth = lwidth)
            hold(True)
            ax23.set_yticks(np.linspace(0,1.0,0.2))
            
            #アミノ酸の各binに含まれる数を足しあげる。
            n, bins, patches =ax7.hist(ConList1[3][~np.isnan(ConList1[3])],  normed=F,cumulative=F,bins = 20,color='blue',alpha=al,histtype=hype)        
            yb = np.add.accumulate(n) / n.sum();xb = np.arange(0,1.05,0.05)#np.convolve(bins, np.ones(2) / 2, mode="same")[1:]

            yb = CumsumCountList(ConList1[3][~np.isnan(ConList1[3])])
            yb, xb = perform4plot(yb, np.arange(0,1.05,0.05))            
            ax23 = ax8.twinx();lines = ax23.plot(xb, yb, ls='-', color='blue', marker='o',label='Cumulative ratio',markersize=msize,linewidth = lwidth)
            hold(True)
            ax23.set_yticks(np.linspace(0,1.0,0.2))
            
            #イオンの各binに含まれる数を足しあげる。
            n, bins, patches =ax7.hist(ConList1[4][~np.isnan(ConList1[4])],  normed=F,cumulative=F,bins = 20,color='purple',alpha=al,histtype=hype)        
            ybl = np.add.accumulate(n) / n.sum();xbl = np.arange(0,1.05,0.05)#np.convolve(bins, np.ones(2) / 2, mode="same")[1:]

            ybl = CumsumCountList(ConList1[4][~np.isnan(ConList1[4])])
            ybl, xbl = perform4plot(ybl, np.arange(0,1.05,0.05))            
            ax23 = ax8.twinx();lines = ax23.plot(xbl, ybl, ls='-', color='purple', marker='o',label='Cumulative ratio',markersize=msize,linewidth = lwidth)
            hold(True)
            ax23.set_yticks(np.linspace(0,1.0,0.2))
            
            #その他の各binに含まれる数を足しあげる。
            n, bins, patches =ax7.hist(ConList1[5][~np.isnan(ConList1[5])],  normed=F,cumulative=F,bins = 20,color='black',alpha=al,histtype=hype)        
            ybl = np.add.accumulate(n) / n.sum();xbl = np.arange(0,1.05,0.05)#np.convolve(bins, np.ones(2) / 2, mode="same")[1:]

            ybl = CumsumCountList(ConList1[5][~np.isnan(ConList1[5])])
            ybl, xbl = perform4plot(ybl, np.arange(0,1.05,0.05))            
            ax23 = ax8.twinx();lines = ax23.plot(xbl, ybl, ls='-', color='black', marker='o',label='Cumulative ratio',markersize=msize,linewidth = lwidth)
            hold(True)
            ax23.set_yticks(np.linspace(0,1.0,0.2))
            ax23.tick_params(labelbottom=False, labelleft=False, labelright='off', labeltop=False)
            #ax10.tick_params(labelright="off")#,left="off")
            ax23.axis('off')          
            
            #異なる代謝同士の各binに含まれる数を足しあげる。
            n, bins, patches =ax7.hist(ConList1[6][~np.isnan(ConList1[6])],  normed=F,cumulative=F,bins = 20,color='gray',alpha=al,histtype=hype)        
            ygr = np.add.accumulate(n) / n.sum();xgr = np.arange(0,1.05,0.05)#np.convolve(bins, np.ones(2) / 2, mode="same")[1:]

            ygr = CumsumCountList(ConList1[6][~np.isnan(ConList1[6])])
            ygr, xgr = perform4plot(ygr, np.arange(0,1.05,0.05))
            ax23 = ax8.twinx();linygres = ax23.plot(xgr, ygr, ls='-', color='gray', marker='o',label='Cumulative ratio', markersize=msize,linewidth = lwidth)
            
            """
            #ax7.set_ylim([0,1.0]);ax8.set_ylim([0,1.0])
            #ax23.axis('off'); ax7.axis('off'); ax8.axis('off')

            #ax23.set_yticks(np.linspace(0,1.0,0.2));ax7.set_yticks(np.linspace(0,1.0,0.2))
            plt.savefig(save_dir + 'AbsHistAveSameMol'+FileNameTag+'.pdf') 
            plt.close()
            
            #######################################33#正だけ_同じ代謝物:累積分布
            fig7, ax7 = plt.subplots(1,1)
            ConListPos1 = [listposr1, listposm1, listposg1, listposb1, listpospl1, listposbl1, listposgr1];         ConList2 = [listr2, listm2, listg2, listb2, listpl2, listbl2, listgr2];
    
            ConListPos1 = [np.array(ConListPos1[x]) for x in range(len(ConListPos1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]

            fig7, ax7 = plt.subplots(1,1)
            x = np.sort(ConListPos1[0][~np.isnan(ConListPos1[0])]); y=np.arange(1,len(x)+1) / len(x)
            x=np.insert(x,0,0);y=np.insert(y,0,0)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[0])
            x = np.sort(ConListPos1[1][~np.isnan(ConListPos1[1])]); y=np.arange(1,len(x)+1) / len(x)
            x=np.insert(x,0,0);y=np.insert(y,0,0)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[1])

            x = np.sort(ConListPos1[2][~np.isnan(ConListPos1[2])]); y=np.arange(1,len(x)+1) / len(x)
            x=np.insert(x,0,0);y=np.insert(y,0,0)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[2])

            x = np.sort(ConListPos1[3][~np.isnan(ConListPos1[3])]); y=np.arange(1,len(x)+1) / len(x)
            x=np.insert(x,0,0);y=np.insert(y,0,0)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[3])

            x = np.sort(ConListPos1[4][~np.isnan(ConListPos1[4])]); y=np.arange(1,len(x)+1) / len(x)
            x=np.insert(x,0,0);y=np.insert(y,0,0)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[4])
            x = np.sort(ConListPos1[5][~np.isnan(ConListPos1[5])]); y=np.arange(1,len(x)+1) / len(x)
            x=np.insert(x,0,0);y=np.insert(y,0,0)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[5])
            x = np.sort(ConListPos1[6][~np.isnan(ConListPos1[6])]); y=np.arange(1,len(x)+1) / len(x)
            x=np.insert(x,0,0);y=np.insert(y,0,0)
            ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=CList[6])
            plt.savefig(save_dir + 'AbsHistAveSameMol'+FileNameTag+'_Pos.pdf')  
    
            ##################################3 ##################################3 ##################################3 ##################################3 ##################################3
            fig10, ax9 = plt.subplots(1,1)
            fig10, ax10 = plt.subplots(1,1)
            
            n, bins, patches =ax9.hist(ConList2[0][~np.isnan(ConList2[0])], density=F, cumulative=F,bins = 20,color='red',alpha=al,histtype=hype)        
            yr = np.add.accumulate(n) / n.sum();xr = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            # 第2軸のプロット
            ax23 = ax10.twinx();lines = ax23.plot(xr, yr, ls='-', color='red', marker='o',label='Cumulative ratio')
 
            n, bins, patches =ax9.hist(ConList2[1][~np.isnan(ConList2[1])],  density=F,cumulative=F,bins = 20,color='magenta',alpha=al,histtype=hype)        
            ym = np.add.accumulate(n) / n.sum();xm = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            ax23 = ax10.twinx();lines = ax23.plot(xm, ym, ls='-', color='magenta', marker='o',label='Cumulative ratio')

            n, bins, patches =ax9.hist(ConList2[2][~np.isnan(ConList2[2])],  density=F,cumulative=F,bins = 20,color='green',alpha=al,histtype=hype)        
            yg = np.add.accumulate(n) / n.sum();xg = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            ax23 = ax10.twinx();lines = ax23.plot(xg, yg, ls='-', color='green', marker='o',label='Cumulative ratio')

            n, bins, patches =ax9.hist(ConList2[3][~np.isnan(ConList2[3])],  density=F,cumulative=F,bins = 20,color='blue',alpha=al,histtype=hype)        
            yb = np.add.accumulate(n) / n.sum();xb = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            ax23 = ax10.twinx();lines = ax23.plot(xb, yb, ls='-', color='blue', marker='o',label='Cumulative ratio')

            n, bins, patches =ax9.hist(ConList2[4][~np.isnan(ConList2[4])],  density=F,cumulative=F,bins = 20,color='black',alpha=al,histtype=hype)        
            ybl = np.add.accumulate(n) / n.sum();xbl = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            ax23 = ax10.twinx();lines = ax23.plot(xbl, ybl, ls='-', color='purple', marker='o',label='Cumulative ratio')
           
            n, bins, patches =ax9.hist(ConList2[5][~np.isnan(ConList2[5])],  density=F,cumulative=F,bins = 20,color='black',alpha=al,histtype=hype)        
            ybl = np.add.accumulate(n) / n.sum();xbl = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            ax23 = ax10.twinx();lines = ax23.plot(xbl, ybl, ls='-', color='black', marker='o',label='Cumulative ratio')           
    
            n, bins, patches =ax9.hist(ConList2[6][~np.isnan(ConList2[6])],  density=F,cumulative=F,bins = 20,color='black',alpha=al,histtype=hype)        
            ybl = np.add.accumulate(n) / n.sum();xbl = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            ax23 = ax10.twinx();lines = ax23.plot(xbl, ybl, ls='-', color='gray', marker='o',label='Cumulative ratio')

            ax23.tick_params(labelbottom=False, labelleft=False, labelright='off', labeltop=False)
           # ax10.tick_params(labelright="off")#,left="off")
            ax23.axis('off')
            #n, bins, patches = plt.hist(list(tempnp1[~np.isnan(tempnp1)]),  normed=False,cumulative=False,bins = 20,color='red',alpha=al,histtype=hype)  
            # 第2軸用値の算出
            #yr = np.add.accumulate(n) / n.sum();xr = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            # 第2軸のプロット
            #ax22 = ax2.twinx();lines = ax22.plot(xr, yr, ls='--', color='r', marker='o',label='Cumulative ratio')
           
            #n, bins, patches = ax2.hist(list(tempnp2[~np.isnan(tempnp2)]),  normed=F,cumulative=F,bins = 20,color='magenta',alpha=al,histtype=hype)        
            #ym = np.add.accumulate(n) / n.sum();xm = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            #ax23 = ax2.twinx();lines = ax23.plot(xm, ym, ls='--', color='magenta', marker='o',label='Cumulative ratio')
            
            #n, bins, patches = ax2.hist(list(tempnp3[~np.isnan(tempnp3)]),  normed=F,cumulative=F,bins = 20,color='green',alpha=al,histtype=hype)        
            #yg = np.add.accumulate(n) / n.sum();xg = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            #ax23 = ax2.twinx();lines = ax23.plot(xg, yg, ls='--', color='green', marker='o',label='Cumulative ratio')
            
            #n, bins, patches = ax2.hist(list(tempnp4[~np.isnan(tempnp4)]),  normed=True,cumulative=True,bins = 20,color='blue',alpha=al,histtype=hype)        
            #yb = np.add.accumulate(n) / n.sum();xb = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            #ax23 = ax2.twinx();lines = ax23.plot(xb, yb, ls='--', color='blue', marker='o',label='Cumulative ratio')

            #n, bins, patches = ax2.hist(list(tempnp5[~np.isnan(tempnp5)]),  normed=True,cumulative=True,bins = 20,color='black',alpha=al,histtype=hype)        
            #ybl = np.add.accumulate(n) / n.sum();xbl = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]
            #ax23 = ax2.twinx();lines = ax23.plot(xbl, ybl, ls='--', color='black', marker='o',label='Cumulative ratio')
            ax9.axis('off'); ax10.axis('off')
            plt.savefig(save_dir + 'AbsHistStdSameMol'+FileNameTag+'.pdf') 
            plt.close()

            
        elif NormlSwitch == 2:#累積分布にする    
            FileNameTag = 'Cumulative'
            fig2, ax2 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = abs(SubjRmean[Col[0]:Col[(c['red'])]]).values.flatten(); tempnp4 = abs(SubjRmean[Col[green]:Col[blue]]).values.flatten(); tempnp2 = abs(SubjRmean[Col[(c['red'])]:Col[magenta]]).values.flatten(); tempnp3 = abs(SubjRmean[Col[magenta]:Col[green]]).values.flatten();  tempnp5 = abs(SubjRmean[Col[blue]:Col[black]]).values.flatten()
            ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)])], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')        
            plt.savefig(save_dir + 'AbsHistAve_'+FileNameTag+'.pdf')        

 #weights = np.ones(len(data))/float(len(data))
        #絶対値ではない           
            fig2, ax2 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None);
            tempnp1 = SubjRmean[Col[0]:Col[(c['red'])]].values.flatten(); tempnp4 = SubjRmean[Col[green]:Col[blue]].values.flatten(); tempnp2 = SubjRmean[Col[(c['red'])]:Col[magenta]].values.flatten(); tempnp3 = SubjRmean[Col[magenta]:Col[green]].values.flatten();  tempnp5 = SubjRmean[Col[blue]:Col[black]].values.flatten()
 #           w= np.ones(len(tempnp1[~np.isnan(tempnp1)]))/float(len(tempnp1[~np.isnan(tempnp1)]));w = np.insert(w,-1,np.ones(len(tempnp2[~np.isnan(tempnp2)]))/float(len(tempnp2[~np.isnan(tempnp2)])));w= np.insert(w,-1,np.ones(len(tempnp3[~np.isnan(tempnp3)]))/float(len(tempnp3[~np.isnan(tempnp3)])));w= np.insert(w,-1,np.ones(len(tempnp4[~np.isnan(tempnp4)]))/float(len(tempnp4[~np.isnan(tempnp4)])));w= np.insert(w,-1,np.ones(len(tempnp5[~np.isnan(tempnp5)]))/float(len(tempnp5[~np.isnan(tempnp5)])))
            w= list(np.ones(len(tempnp1[~np.isnan(tempnp1)]))/float(len(tempnp1[~np.isnan(tempnp1)])));w += [list(np.ones(len(tempnp2[~np.isnan(tempnp2)]))/float(len(tempnp2[~np.isnan(tempnp2)])))];w+= [list(np.ones(len(tempnp3[~np.isnan(tempnp3)]))/float(len(tempnp3[~np.isnan(tempnp3)])))];w+= [list(np.ones(len(tempnp4[~np.isnan(tempnp4)]))/float(len(tempnp4[~np.isnan(tempnp4)])))];w+= [list(np.ones(len(tempnp5[~np.isnan(tempnp5)]))/float(len(tempnp5[~np.isnan(tempnp5)])))]

            print(w);print(tempnp1[~np.isnan(tempnp1)],tempnp2[~np.isnan(tempnp2)],tempnp3[~np.isnan(tempnp3)],tempnp4[~np.isnan(tempnp4)],tempnp5[~np.isnan(tempnp5)])
            ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)])], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')   #,weights=w     
            plt.savefig(save_dir + 'HistAve_'+FileNameTag+'.pdf')        
        #Stdもやる
            fig3, ax3 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = SubjRstd[Col[0]:Col[(c['red'])]].values.flatten();  tempnp4 = SubjRstd[Col[green]:Col[blue]].values.flatten(); tempnp2 = SubjRstd[Col[(c['red'])]:Col[magenta]].values.flatten();  tempnp3 = SubjRstd[Col[magenta]:Col[green]].values.flatten();   tempnp5 = SubjRstd[Col[blue]:Col[black]].values.flatten()
            ax3 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)])], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')        
            plt.savefig(save_dir + 'HistStd_'+FileNameTag+'.pdf')  
            
            """
            #両方がその代謝系で、それ以外はグレー
            fig4, ax4 = plt.subplots(1,1)
            ConList1 = [listr1, listm1, listg1, listb1, listbl1, listgr1];         ConList2 = [listr2, listm2, listg2, listb2, listbl2, listgr2];
    
            ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
    
            ax4 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'],cumulative=True,histtype='stepfilled')        
            plt.savefig(save_dir + 'HistAveSameMol_'+FileNameTag+'.pdf')  
    
            fig5, ax5 = plt.subplots(1,1)
            ax4 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'],cumulative=True,histtype='stepfilled')        
    
            #ax5 = plt.hist([listr2, listm2, listg2, listb2, listbl2, listgr2], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'])        
            plt.savefig(save_dir + 'HistStdSameMol_'+FileNameTag+'.pdf')  
    
            fig6, ax6 = plt.subplots(1,1)
            ax6 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])]], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')        
            plt.savefig(save_dir + 'HistAveSameMolWOGray_'+FileNameTag+'.pdf') 
    
            fig7, ax7 = plt.subplots(1,1)
            ax7 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])]], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')        
            plt.savefig(save_dir + 'HistStdSameMolWOGray_'+FileNameTag+'.pdf')  
            
            fig8, ax8 = plt.subplots(1,1)
            ConList1 = [listrabs1, listmabs1, listgabs1, listbabs1, listblabs1, listgrabs1];         ConList2 = [listrabs2, listmabs2, listgabs2, listbabs2, listblabs2, listgrabs2];
    
            ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
    
            ax8 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'],cumulative=True,histtype='stepfilled'        )
            plt.savefig(save_dir + 'AbsHistAveSameMol_'+FileNameTag+'.pdf')  
    
            fig9, ax9 = plt.subplots(1,1)
            ax9 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])],ConList2[5][~np.isnan(ConList2[5])]], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'],cumulative=True,histtype='stepfilled')        
    
            #ax5 = plt.hist([listr2, listm2, listg2, listb2, listbl2, listgr2], stacked=True, bins = 20,color=['red','magenta','green','blue','black','gray'])        
            plt.savefig(save_dir + 'AbsHistStdSameMol_'+FileNameTag+'.pdf')
            
            fig10, ax10 = plt.subplots(1,1)
            ax10 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])]], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')        
            plt.savefig(save_dir + 'AbsHistAveSameMolWOGray_'+FileNameTag+'.pdf') 
    
            fig11, ax11 = plt.subplots(1,1)
            ax11 = plt.hist([ConList2[0][~np.isnan(ConList2[0])],ConList2[1][~np.isnan(ConList2[1])],ConList2[2][~np.isnan(ConList2[2])],ConList2[3][~np.isnan(ConList2[3])],ConList2[4][~np.isnan(ConList2[4])]], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')        
            plt.savefig(save_dir + 'AbsHistStdSameMolWOGray_'+FileNameTag+'.pdf') 
            """
    
    elif SwitchDict['Incretin']  == 1:#Incretin枠なら
        MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolColor_Eng.xlsx',header=0,encoding = "ISO-8859-1") # 日本語
        
        fig1, ax1 = plt.subplots(1,1)#abs
        fig2, ax2 = plt.subplots(1,1)#normal

        Col = SubjRmean.columns;ColorList = MolColorDF.loc[Col,'MolColor'].dropna()
        #MolColorDF.loc[Col,:].isnull()[0]
        idxbool,colbool = [np.logical_not(np.array(MolColorDF.loc[Col,:].isnull()[0]))]*2
        SubjRmean = SubjRmean.iloc[idxbool,colbool];SubjRstd = SubjRstd.iloc[idxbool,colbool]
        #ColorList = MolColorDF.loc[Col,'MolColor'];
        Col = SubjRmean.columns
        print(Col)
        c = collections.Counter(ColorList);listr1=[];listr2=[];listm1=[];listm2=[];listg1=[];listg2=[];listb1=[];listb2=[];listbl1=[];listpl1=[];listpl2=[];listbl2=[];listgr1=[];listgr2=[]
        listrabs1=[];listrabs2=[];listmabs1=[];listmabs2=[];listgabs1=[];listgabs2=[];listbabs1=[];listbabs2=[];listblabs1=[];listplabs1=[];listplabs2=[];listblabs2=[];listgrabs1=[];listgrabs2=[]
        for i in range(0,len(SubjRmean.columns)):#列を1~84まで
            if i < (c['red']):
                listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[(c['red']-1)]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(c['red']-1)]])
                list1=list(SubjRmean.loc[Col[i],Col[i]:Col[(c['red']-1)]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(c['red']-1)]])

                #list1.remove(1);list2.remove(0);#Std.remove(1)
                ax1.scatter(listabs1,listabs2,color=ColorList[i])
                
                ax2.scatter(list1,list2,color=ColorList[i])

                listrabs1 += listabs1; listrabs2 += listabs2;  listr1 += list1; listr2 += list2
                #GH.ScatterWHeatmap(np.array(listrabs1)[~np.isnan(listrabs1)],np.array(listrabs2)[~np.isnan(listrabs2)],AbsAvemin,AbsAvemax,Stdmin,Stdmax,save_dir+'_red',bin=30)#散布図をヒートマップで
                #GH.ScatterWColorbar(np.array(listrabs1)[~np.isnan(listrabs1)],np.array(listrabs2)[~np.isnan(listrabs2)],0,0,1,save_dir+'red')
                
                #カーネル密度推定して色分け
                plt.close()
                
                stloc = Col[(c['red'])]
                #blue = c['red']+c['blue']-1
                #print(Col[i],ColorList[i])
                
                #plt.close()
                #stloc = Col[magenta]
                green = c['red'] + c['green']
                #print(Col[i],ColorList[i])
            elif i < green:
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[green]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[green]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
    
                    #list1.remove(1);list2.remove(0);#Std.remove(1)
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i]) 
                    listg1 += list1; listg2 += list2;  listgabs1 += listabs1; listgabs2 += listabs2
                    #GH.ScatterWHeatmap(np.array(listgabs1)[~np.isnan(listgabs1)],np.array(listgabs2)[~np.isnan(listgabs2)],AbsAvemin,AbsAvemax,Stdmin,Stdmax,save_dir+'_green',bin=30)#散布図をヒートマップで
                    #GH.ScatterWColorbar(np.array(listgabs1)[~np.isnan(listgabs1)],np.array(listgabs2)[~np.isnan(listgabs2)],0,save_dir+'green')
                    
                    plt.close()
                
                    stloc = Col[green]
                    blue = green + c['blue'] -1
            elif i < (blue):
                list1=list(SubjRmean.loc[Col[i],Col[i]:Col[blue]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[blue]])
                listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[blue]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[blue]])

                #list1.remove(1);list2.remove(0);#Std.remove(1)
                ax1.scatter(listabs1,listabs2,color=ColorList[i])
                ax2.scatter(list1,list2,color=ColorList[i]) 
                listb1 += list1; listb2 += list2;  listbabs1 += listabs1; listbabs2 += listabs2
                #GH.ScatterWHeatmap(np.array(listbabs1)[~np.isnan(listbabs1)],np.array(listbabs2)[~np.isnan(listbabs2)],AbsAvemin,AbsAvemax,Stdmin,Stdmax,save_dir+'_blue',bin=30)#散布図をヒートマップで
                #GH.ScatterWColorbar(np.array(listbabs1)[~np.isnan(listbabs1)],np.array(listbabs2)[~np.isnan(listbabs2)],0,save_dir+'blue')
                
                plt.close()                
                stloc = Col[blue]
                #purple = blue + c['purple']
                #print(Col[i],ColorList[i])

                plt.close()
            list1=list(SubjRmean.loc[Col[i],stloc:]);list2 =list(SubjRstd.loc[Col[i],stloc:])
            listabs1=list(abs(SubjRmean.loc[Col[i],stloc:]));listabs2 =list(SubjRstd.loc[Col[i],stloc:])

            listgr1 += list1;listgr2 += list2;    listgrabs1 += listabs1;listgrabs2 += listabs2
#            GH.ScatterWColorbar(np.array(listgrabs1)[~np.isnan(listgrabs1)],np.array(listgrabs2)[~np.isnan(listgrabs2)],0,save_dir+'gray')

            #listgr1.remove(1);#listgr2.remove(1);#Std.remove(1)
            
            ax1.scatter(listabs1,listabs2,color='gray',alpha=0.5)
            ax2.scatter(list1,list2,color='gray',alpha=0.5)  
        ColDict={}
        #ColDict['red']=c[red; 
        ColDict['green'] = green; ColDict['blue']=blue
        AbsAve = list(abs(SubjRmean.values.flatten()[~np.isnan(SubjRmean.values.flatten())]))
        Std = list(SubjRstd.values.flatten()[~np.isnan(SubjRstd.values.flatten())])
        absr, absp=pearsonr(AbsAve,Std)
        p='{:e}'.format(absp)
        SaveName='_WMolColor_Abs'
            
        ax1.set_title('R=' + str(round(absr,3)) + ', p=' + p,fontsize=10);ax2.set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp,fontsize=10)
        ax1.set_xlabel('Abs(Ave)',fontsize=10);ax2.set_xlabel('Ave',fontsize=10)              
        ax1.set_ylabel('Std',fontsize=10);ax2.set_ylabel('Std',fontsize=10)#,fontproperties=prop)     fig.tight_layout()
        plot_axis = plt.axis()
        fig1.savefig(save_dir +'AbsAveStdScatter' + SaveName +'.pdf',format='pdf')
        fig2.savefig(save_dir +'AveStdScatter' + SaveName +'.pdf',format='pdf')
        plt.close()                 
    #0で　:#とりあえず簡単に片方の代謝系で色分け__ヒストグラム描画　1で代謝ごとの累積分布、2で分布の割合
    if SwitchDict['Incretin_hist'] == 'Anywayhist':
        print('YYYY')
        PlotHist(SubjRmean,Col,c, ColDict, NormlSwitch)
        


def mkPneltoCalcAveVar(AdjustSwitch,FastingDatadict,SubjectName,save_dir):#j被験者20人をまとめて、
    if AdjustSwitch == 0:#全時系列なら
        Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/RDFiwashi.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRPanel = pd.Panel({'iwana':pd.DataFrame(data=None,index=list(Kpc.index),columns=list(Kpc.columns))})
        for i in range(0,len(SubjectName)):
            SubjRPanel[SubjectName[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/RDF' + SubjectName[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRmean = SubjRPanel.mean(axis=0)
        SubjRstd = SubjRPanel.std(axis=0, ddof=0)
        
        SubjRmean.to_excel(save_dir + 'Rmean.xlsx')
        SubjRstd.to_excel(save_dir + 'Rstd.xlsx')  
        SubjRmeanrev = MatHel.adjustMatrlower(SubjRmean); SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
        SubjRmeanrev.to_excel(save_dir + 'Rmean_rev.xlsx')
        SubjRstdrev.to_excel(save_dir + 'Rstd_rev.xlsx')
        

    #生データのパネルを作り
    elif AdjustSwitch == 1:#全時系列なら
        TimeList=[120]#0とTimeListの相関
        Kraw = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/Data/TimeCourse/SubjTmCs_kamasuRaw.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRPanel = pd.Panel({'kanpachi':pd.DataFrame(data=None,index=list(Kraw.index),columns=list(Kraw.columns))})
        SubjTmCsDfNew = pd.DataFrame(data=None,columns=list(Kraw.columns));
        for i in range(0,len(SubjectName)):
            SubjRPanel[SubjectName[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/Data/TimeCourse/SubjTmCs_' + SubjectName[i] + 'Raw.xlsx',header=0,encoding = "ISO-8859-1")
            SubjTmCsDfNew.loc[i,:]=[FastingDatadict[Label[ik]][i] for ik in range(0,84)] ; SubjTmCsDfNew.loc[i+TimeList[0],:] = list(SubjRPanel[SubjectName[i]].loc[120])
        #列はラベル、行に[Fasting,120min]*被験者ぶん
        
        RDF,PDF =  AmT.CorrMolEachSubjAdjustHelper(SubjTmCsDfNew,MolDict,Label)
        RDF.to_excel(save_dir + 'RDF.xlsx')
        PDF.to_excel(save_dir + 'PDF.xlsx')
    elif AdjustSwitch == str(2)+'_old':#全時系列なら
        Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/SubjRPDF_Eng/PDFexp_SubjectNo1_25_B.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRPanel = pd.Panel({'kanpachi':pd.DataFrame(data=None,index=list(Kpc.index),columns=list(Kpc.columns))})
        for i in range(0,len(SubjectName)):
            SubjRPanel[SubjectName[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/SubjRPDF_Eng/RDF' + SubjectName[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRmean = SubjRPanel.mean(axis=0)
        SubjRstd = SubjRPanel.std(axis=0, ddof=0)
        
        SubjRmean.to_excel(save_dir + 'Rmean.xlsx')
        SubjRstd.to_excel(save_dir + 'Rstd.xlsx') 
        
    elif AdjustSwitch == 2:#Incretin枠のデータ各被験者・各実験条件分
        Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Delta/Each3Subj6Cond/RDFexp_SubjectNo1_25_B.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRPanel = pd.Panel({'No1_25_B':pd.DataFrame(data=None,index=list(Kpc.index),columns=list(Kpc.columns))})#パネルの枠作り
        if FastingDatadict['Target'] == '25_Bolus':#各条件で平均取るときのターゲット：'25_Bolus','50_Bolus','75_Bolus'.'25_Continuous','50_Continuous','75_Continuous'
            tempList =['1_25_B','2_25_B','3_25_B']
            
        elif FastingDatadict['Target'] == '50_Bolus':
            tempList =['1_50_B','2_50_B','3_50_B']            
        elif FastingDatadict['Target'] == '75_Bolus':
            tempList =['1_75_B','2_75_B','3_75_B']    
        elif FastingDatadict['Target'] == '25_Continuous':
            tempList =['1_25_R_2h','2_25_R_2h','3_25_R_2h']    
        elif FastingDatadict['Target'] == '50_Continuous':
            tempList =['1_50_R_2h','2_50_R_2h','3_50_R_2h']    
        elif FastingDatadict['Target'] == '75_Continuous':
            tempList =['1_75_R_2h','2_75_R_2h','3_75_R_2h']                
        else:#全員
            tempList =['1_25_B','1_25_R_2h', '1_50_B', '1_50_R_2h', '1_75_B', '1_75_R_2h','2_25_B', '2_25_R_2h', '2_50_B', '2_50_R_2h', '2_75_B', '2_75_R_2h','3_25_B', '3_25_R_2h', '3_50_B', '3_50_R_2h', '3_75_B', '3_75_R_2h' ]

        if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
                SubjRPanel[tempList[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Each3Subj6Cond/RDFexp_SubjectNo' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
        else:#差分なら        
            for i in range(0,len(tempList)):
                SubjRPanel[tempList[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Delta/Each3Subj6Cond/RDFexp_SubjectNo' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
                SubjRPanel[tempList[i]] = pd.read_excel(save_dir+'/RDFexp_SubjectNo' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")

                #SubjRmeanrev = MatHel.adjustMatrlower(SubjRPanel[tempList[i]]); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
                #MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','Each3Subj6Cond/' + tempList[i])
                #MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','Each3Subj6Cond/' + tempList[i])
            
            
            SubjRmean = SubjRPanel.mean(axis=0)
            SubjRstd = SubjRPanel.std(axis=0, ddof=0)
            SubjRmean.to_excel(save_dir +'Each3Subj6Cond/'+ 'Rmean_'+FastingDatadict['Target']+'.xlsx')
            SubjRstd.to_excel(save_dir + 'Each3Subj6Cond/'+'Rstd_'+FastingDatadict['Target']+'.xlsx') 
            
            SubjRmeanrev = MatHel.adjustMatrlower(SubjRmean); SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
    
            MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','Each3Subj6Cond/AllR_'+FastingDatadict['Target'] )
            MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','Each3Subj6Cond/Allabs(R)_'+FastingDatadict['Target'] )
            
            

            SubjRmeanrev.to_excel(save_dir +'Each3Subj6Cond/'+ 'Rmean_Upper_'+FastingDatadict['Target']+'.xlsx')
            SubjRstdrev.to_excel(save_dir +'Each3Subj6Cond/'+ 'Rstd_Upper_'+FastingDatadict['Target']+'.xlsx')
            
            
    elif AdjustSwitch == 4:#全被験者繋いだ、40タイムコースで、各分子間の相関は？(81×40) が6条件 
        Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Comb3Subj/RDFexp_SubjectNo1_25_B.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRPanel = pd.Panel({'RDFexp_SubjectNo1_25_B':pd.DataFrame(data=None,index=list(Kpc.index),columns=list(Kpc.columns))})
        tempList =['RDFexp_SubjectNo1_25_B', 'RDFexp_SubjectNo1_25_R_2h', 'RDFexp_SubjectNo1_50_B', 'RDFexp_SubjectNo1_50_R_2h', 'RDFexp_SubjectNo1_75_B', 'RDFexp_SubjectNo1_75_R_2h' ]
        
        if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
                SubjRPanel[tempList[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Comb6Cond/' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
        else:#差分なら
            for i in range(0,len(tempList)):
                SubjRPanel[tempList[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Delta/Comb3Subj/' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
            SubjRmean = SubjRPanel.mean(axis=0)
            SubjRstd = SubjRPanel.std(axis=0, ddof=0)
            SubjRmeanrev = MatHel.adjustMatrlower(SubjRmean); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
            #SubjRmeanrev  = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180913/Eng/6Cond/Rmeanrev.xlsx',header=0,encoding = "ISO-8859-1")
            MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','Comb3Subj/AllR' )
            MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','Comb3Subj/Allabs(R)' )
            
            SubjRmean.to_excel(save_dir + 'Comb3Subj/'+'Rmean.xlsx')
            SubjRstd.to_excel(save_dir + 'Comb3Subj/' + 'Rstd.xlsx') 

    elif AdjustSwitch ==  5:#Incretin枠各条件平均全被験者平均の2分子内相関：Ave6Cond
        Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Ave6Cond/RDFexp_25_B.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRPanel = pd.Panel({'_25_B':pd.DataFrame(data=None,index=list(Kpc.index),columns=list(Kpc.columns))})

        tempList =['_25_B','_25_R_2h', '_50_B', '_50_R_2h', '_75_B', '_75_R_2h']
        for i in range(0,len(tempList)):
            SubjRPanel[tempList[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Ave6Cond/RDFexp' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
            SubjRmeanrev = MatHel.adjustMatrlower(SubjRPanel[tempList[i]]); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
            MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','Ave6Cond/' + tempList[i])
            MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','Ave6Cond/' + tempList[i])
    
        SubjRmean = SubjRPanel.mean(axis=0)
        SubjRstd = SubjRPanel.std(axis=0, ddof=0)
        SubjRmeanrev = MatHel.adjustMatrlower(SubjRmean); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに

        MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','Ave6Cond/AllR' )
        MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','Ave6Cond/Allabs(R)' )

        SubjRmean.to_excel(save_dir + 'Ave6Cond/' +  'Rmean.xlsx')
        SubjRstd.to_excel(save_dir + 'Ave6Cond/' +  'Rstd.xlsx') 

    elif AdjustSwitch ==  6:#Incretin枠各条件平均全条件繋げて2分子内相関、   CombAve6Cond     
        Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/CombAve6Cond/RDFexp_SubjectNo3_75_R_2hWoDel.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRmeanrev = MatHel.adjustMatrlower(Kpc); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
        MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','CombAve6Cond/AllR')
        MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','CombAve6Cond/Allabs(R)')
        
    elif AdjustSwitch ==  7: #全被験者繋いで、6条件求める：Comb6Cond
        Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Comb6Cond/exp_SubjectNo3WoDel.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRPanel = pd.Panel({'exp_SubjectNo1WoDel':pd.DataFrame(data=None,index=list(Kpc.index),columns=list(Kpc.columns))})
        tempList =['exp_SubjectNo1WoDel','exp_SubjectNo2WoDel', 'exp_SubjectNo3WoDel']

        if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
                SubjRPanel[tempList[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Comb6Cond/' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
        else:#差分なら
            
            for i in range(0,len(tempList)):
                SubjRPanel[tempList[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Delta/Comb6Cond/' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
                SubjRmeanrev = MatHel.adjustMatrlower(SubjRPanel[tempList[i]]); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
                MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','Comb6Cond/' + tempList[i])
                MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','Comb6Cond/' + tempList[i])
            SubjRmean = SubjRPanel.mean(axis=0)
            SubjRstd = SubjRPanel.std(axis=0, ddof=0)
            SubjRmeanrev = MatHel.adjustMatrlower(SubjRmean); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
    
            MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','Comb6Cond/AllR' )
            MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','Comb6Cond/Allabs(R)' )
    
            SubjRmean.to_excel(save_dir + 'Comb6Cond/' +  'Rmean.xlsx')
            SubjRstd.to_excel(save_dir + 'Comb6Cond/' +  'Rstd.xlsx') 
    
    elif AdjustSwitch ==  8: #6条件×全被験者繋いだ、Comb3SubjComb6Cond
    
        Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Comb3SubjComb6Cond/RDFexp_SubjectNo3_75_R_2h.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRmeanrev = MatHel.adjustMatrlower(Kpc); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
        MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','Comb3SubjComb6Cond/AllR')
        MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','Comb3SubjComb6Cond/Allabs(R)')
        
    elif AdjustSwitch ==  'Continuous': #条件がContinuousなら
        Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Continuous/InterMolCorr/RDFisaki.xlsx',header=0,encoding = "ISO-8859-1")
        SubjRPanel = pd.Panel({'isaki':pd.DataFrame(data=None,index=list(Kpc.index),columns=list(Kpc.columns))})
        tempList =['iwashi','karei','shimaaji','unagi']

        if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
                SubjRPanel[tempList[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Comb6Cond/' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
        else:#差分なら
            
            for i in range(0,len(tempList)):
                SubjRPanel[tempList[i]] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Continuous/InterMolCorr/RDF' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
                SubjRmeanrev = MatHel.adjustMatrlower(SubjRPanel[tempList[i]]); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
                MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient', tempList[i]+'/Hist')
                MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient', tempList[i] +'/Hist')
            SubjRmean = SubjRPanel.mean(axis=0)
            SubjRstd = SubjRPanel.std(axis=0, ddof=0)
            SubjRmeanrev = MatHel.adjustMatrlower(SubjRmean); #SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
    
            MSMSH.mkhist(SubjRmeanrev,save_dir,'Correlation coefficient','All/AllR' )
            MSMSH.mkhistabs(SubjRmeanrev,save_dir,'Correlation coefficient','All/Allabs(R)' )
    
            SubjRmean.to_excel(save_dir + 'All/' +  'Rmean.xlsx')
            SubjRstd.to_excel(save_dir + 'All/' +  'Rstd.xlsx')     
    
    return(SubjRPanel,SubjRmean,SubjRstd)


def CalcCorrBetw2MolEachSubjContinuous(SubjTimeSeriesDF,SubjectName,MolDict,Label,AdjustSwitch,FastingDatadict,save_dir):

        if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
                SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Comb6Cond/' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
        else:#差分なら 
            for i in range(0,len(SubjectName)):#被験者の数だけ
                SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Continuous/TimeCourse/English/Delta/' + SubjectName[i] + '_Continuous2h_Delta.xlsx',header=0,encoding = "ISO-8859-1")
                SubjTmCsDf=SubjTmCsDf.drop('time(min)',axis=1)
                Label = list(SubjTmCsDf.columns)
                RDF,PDF =  AmT.CorrMolEachSubjHelper(SubjTmCsDf,Label,Label)
                RDF.to_excel(save_dir + 'RDF'+SubjectName[i] + '.xlsx')
                PDF.to_excel(save_dir + 'PDF'+SubjectName[i] + '.xlsx') 
    
def CalcCorrBetw2MolEachSubj(SubjTimeSeriesDF,SubjectName,MolDict,Label,AdjustSwitch,FastingDatadict,save_dir):
    #ある被験者のみだす
    #Label.remove('GLP-1');Label.remove('Aspartic acid');Label.remove('1-methyl-histidine');Label.remove('Acetoacetate')
    ##Label.remove('Free Fatty Acid');Label.remove('M ethanolamine');Label+=['Free fatty acid']; Label+=['Monoethanolamine']
    #MolDict.pop('GLP-1');MolDict.pop('Aspartic acid');MolDict.pop('1-methyl-histidine');MolDict.pop('Acetoacetate');
    
    if AdjustSwitch == 0:#全時系列なら
        for i in range(0,len(SubjectName)):
            SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/20180427/SubjTmcs_' + SubjectName[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
            RDF,PDF =  AmT.CorrMolEachSubjHelper(SubjTmCsDf,MolDict,Label)
            RDF.to_excel(save_dir + 'RDF'+SubjectName[i] + '.xlsx')
            PDF.to_excel(save_dir + 'PDF'+SubjectName[i] + '.xlsx')
        
    elif AdjustSwitch ==  1:#ある時点だけなら
        
        for i in range(0,len(SubjectName)):
            SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/20180427/SubjTmcs_' + SubjectName[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
            SubjTmCsDfNew = pd.DataFrame(data=None,columns=list(SubjTmCsDf.columns));SubjTmCsDfNew.loc[0,:]=[FastingDatadict[Label[ik]][i] for ik in range(0,84)] ; SubjTmCsDfNew.loc[TimeList[0],:] = list(SubjTmCsDf.loc[120])
            RDF,PDF =  AmT.CorrMolEachSubjAdjustHelper(SubjTmCsDfNew,MolDict,Label)
            RDF.to_excel(save_dir + 'RDF'+SubjectName[i] + '.xlsx')
            PDF.to_excel(save_dir + 'PDF'+SubjectName[i] + '.xlsx')
        
    elif AdjustSwitch ==  2:#Incretin枠のデータ各被験者・各実験条件分
        if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
                SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Comb6Cond/' + tempList[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
        else:#差分なら             
            for i in range(0,len(SubjectName)):#被験者の数だけ
                SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/TimeCourse/Delta/' + SubjectName[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
                if 'Unit' in list(SubjTmCsDf.index):
                    SubjTmCsDf = SubjTmCsDf.drop('Unit')
                imptime = [ -5,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240]
                DelList = list ( set(list(SubjTmCsDf['time(min)'])) - set( imptime ) ) 
                if 2 in DelList:
                    for jj in DelList:
                        SubjTmCsDf = SubjTmCsDf.drop( SubjTmCsDf.index[SubjTmCsDf['time(min)'] == jj][0] )
                SubjTmCsDf.columns = [ 'Free fatty acid'  if 'Free Fatty Acid' == list(SubjTmCsDf.columns)[i] else list(SubjTmCsDf.columns)[i] for i in range(len(SubjTmCsDf.columns))]
                SubjTmCsDf.columns = [ 'Monoethanolamine'  if 'M ethanolamine' == list(SubjTmCsDf.columns)[i] else list(SubjTmCsDf.columns)[i] for i in range(len(SubjTmCsDf.columns))]
                SubjTmCsDf = SubjTmCsDf.drop(['GLP-1','Aspartic acid','1-methyl-histidine','Acetoacetate'],axis=1)
                RDF,PDF =  AmT.CorrMolEachSubjHelper(SubjTmCsDf,MolDict,Label)
                RDF.to_excel(save_dir + 'RDF'+SubjectName[i] + '.xlsx')
                PDF.to_excel(save_dir + 'PDF'+SubjectName[i] + '.xlsx') 
            
    elif AdjustSwitch ==  3:#全被験者2分子内相関なら
        #差分と生を場合分ける
        NewLabel=list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,index_col=0).index)[:83]
        if FastingDatadict['RawDelta']=='Delta':
            filename = 'Delta'
        else:
            filename = 'Raw'
        if  FastingDatadict['Target'] ==   'Bolus':
            for i in range(0,len(SubjectName)):
                SubjTmCsFast = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/'+filename+'_Eng/New/SubjTmcs_'+ SubjectName[i] +filename+'.xlsx',header=0,index_col=0)
                SubjTmCsFast.columns = NewLabel
                SubjTmCsFast.to_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/'+filename+'_Eng/New/SubjTmcs_'+ SubjectName[i] +filename+'.xlsx')
                #SubjTmCs = SubjTmCs.drop('3-Hydroxybutyrate',axis=1); SubjTmCsFast = SubjTmCs.drop(-10,axis=0); 
                #SubjTmCsFast.columns= NewLabel
                #SubjTmCsFast.loc[0] = SubjTmCs.loc[-10:0].mean()#差分なら全て0になるはず
                #20人分を縦にくっつける
                if i > 0:
                    SubjTmCsNew = pd.concat([SubjTmCsNew,SubjTmCsFast],axis=0)#, join_axes=[SubjTmCsFast.columns])
                else:
                    SubjTmCsNew = SubjTmCsFast
        elif FastingDatadict['Target'] ==   'Continuous':#Continuous*3なら
            for i in range(0,len(SubjectName)):
                SubjTmCsFast = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Continuous/TimeCourse/English/'+filename+'/New/'+ SubjectName[i] +'_Continuous2h_'+filename+'.xlsx',header=0,encoding = "ISO-8859-1")
                #SubjTmCs = SubjTmCs.drop('3-Hydroxybutyrate',axis=1); SubjTmCsFast = SubjTmCs.drop(-10,axis=0); 
                SubjTmCsFast.columns= NewLabel
                #SubjTmCsFast.loc[0] = SubjTmCs.loc[-10:0].mean()#差分なら全て0になるはず
                #20人分を縦にくっつける
                if i > 0:
                    SubjTmCsNew = pd.concat([SubjTmCsNew,SubjTmCsFast],axis=0, join_axes=[SubjTmCsFast.columns])
                else:
                    SubjTmCsNew = SubjTmCsFast  
        elif FastingDatadict['Target'] == 'Both':#BC両方なら
            for i in range(0,len(SubjectName)):#Bolusやって、Continuouもやる
                SubjTmCsFastBolus = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/'+filename+'_Eng/New/SubjTmcs_'+ SubjectName[i] +filename+'.xlsx',header=0,encoding = "ISO-8859-1")
                SubjTmCsFastBolus.columns= [NewLabel[i]+'_B' for i in range(len(NewLabel))]
                #20人分を縦にくっつける
                if i > 0:
                    SubjTmCsBolus = pd.concat([SubjTmCsBolus,SubjTmCsFastBolus],axis=0, join_axes=[SubjTmCsFastBolus.columns])
                else:
                    SubjTmCsBolus = SubjTmCsFastBolus  
                SubjTmCsFastContinuous = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Continuous/TimeCourse/English/'+filename+'/New/'+ SubjectName[i] +'_Continuous2h_'+filename+'.xlsx',header=0,encoding = "ISO-8859-1")
                SubjTmCsFastContinuous.columns= [NewLabel[i]+'_C' for i in range(len(NewLabel))]
                #20人分を縦にくっつける
                if i > 0:
                    SubjTmCsContinuous = pd.concat([SubjTmCsContinuous,SubjTmCsFastContinuous],axis=0, join_axes=[SubjTmCsFastContinuous.columns])
                else:
                    SubjTmCsContinuous = SubjTmCsFastContinuous
                    
                if i == 4:
                    SubjTmCsNew = pd.concat([SubjTmCsBolus,SubjTmCsContinuous],axis=1, join_axes=[SubjTmCsBolus.index])#分子名方向にくっつける
                    NewLabel = list(SubjTmCsNew.columns)
        save_dir= save_dir + filename
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)                   
        SubjTmCsNew.to_excel(save_dir + '/AllSubjTmCsDiff.xlsx')
        #SubjTmCsDfNewEng = CNE.ChnageNametoEnglish(SubjTmCsDfNew,1);
        SubjTmCsNew.to_excel(save_dir + '/AllSubjTmCsDiff_20190228.xlsx')
### 相関係数or分散共分散
        if FastingDatadict['relation']=='pearson':
            RDF,PDF =  AmT.CorrMolEachSubjHelper(SubjTmCsNew,NewLabel,NewLabel)             
        elif FastingDatadict['relation']=='spearman':
            RDF,PDF =  AmT.CorrMolEachSubjHelperSpearman(SubjTmCsNew,NewLabel,NewLabel)  
        elif FastingDatadict['relation']=='cov':
            RDF,PDF =  AmT.CorrMolEachSubjHelperCov(SubjTmCsNew,NewLabel,NewLabel)  

        #RDF,PDF =  AmT.CorrMolEachSubjAdjustHelper(SubjTmCsNew,NewLabel,NewLabel)  
        RDF.to_excel(save_dir + '/RDF.xlsx');RDFrev = RDF.copy();RDFrev =MatHel.adjustMatrlower(RDFrev);RDFrev.to_excel(save_dir + '/RDFrev.xlsx');
        PDF.to_excel(save_dir + '/PDF.xlsx');PDFrev = PDF.copy();PDFrev =MatHel.adjustMatrlower(PDFrev);PDFrev.to_excel(save_dir + '/PDFrev.xlsx');   
        SC.UseR(PDFrev,{'EngLabel':'TPSMqval'})        

        GH.mkhist(RDFrev,save_dir,'RDFrev',20 ,'RDFrev');GH.mkhist(PDFrev,save_dir,'PDFrev',20 ,'PDFrev',)

        GH.mkhist(np.abs(RDFrev),save_dir,'AbsRDFrev',20 ,'AbsRDFrev')
        
        print(RDF)
        return(SubjTmCsNew,RDF,PDF)

    elif AdjustSwitch ==  'MultDF':#全条件繋げて全被験者2分子内相関
        NewLabel = list(SubjTimeSeriesDF)
        RDF,PDF =  AmT.CorrMolEachSubjAdjustHelper(SubjTimeSeriesDF,NewLabel,NewLabel)  
        RDF.to_excel(save_dir + 'RDF.xlsx')
        PDF.to_excel(save_dir + 'PDF.xlsx')
        print('Yes')
        return(SubjTimeSeriesDF,RDF,PDF)

    elif AdjustSwitch ==  4:#Incretin枠全条件繋げて全被験者2分子内相関
        if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
            pass
        elif FastingDatadict['RawDelta']=='Delta':#差分と生を場合分ける
            tempDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/' + SubjectName[0] + '.xlsx',header=0,encoding = "ISO-8859-1")
            NewCounter = [SubjectName[0:6], SubjectName[6:12], SubjectName[12:] ]
            
            for i in range(0,len(NewCounter)):#3人分
                rootDF = pd.DataFrame(data=None, index= None, columns = ( tempDF.columns ))
                for ii in range(len(NewCounter[i])):#6条件
                    SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/' + NewCounter[i][ii] + '.xlsx',header=0,encoding = "ISO-8859-1")
                    if 'Unit' in list(SubjTmCsDf.index):
                        SubjTmCsDf = SubjTmCsDf.drop('Unit')
                    imptime = [ -10,-5,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240]
                    DelList = list ( set(list(SubjTmCsDf['time(min)'])) - set( imptime ) ) 
                    if 2 in DelList:
                        for jj in DelList:
                            SubjTmCsDf = SubjTmCsDf.drop( SubjTmCsDf.index[SubjTmCsDf['time(min)'] == jj][0] )
                        #SubjTmCsDf = SubjTmCsDf[ SubjTmCsDf[ 'time(min)' ] != DelList]
                    rootDF = pd.concat([ rootDF, SubjTmCsDf ],axis=0, join_axes=[ rootDF.columns ])
                rootDF.to_excel(save_dir+NewCounter[i][ii]+'WoDel.xlsx')                    
                RDF,PDF =  AmT.CorrMolEachSubjHelper(rootDF,MolDict,Label)
                RDF.to_excel(save_dir + 'RDF'+SubjectName[i] + '.xlsx')
                PDF.to_excel(save_dir + 'PDF'+SubjectName[i] + '.xlsx') 
    elif AdjustSwitch ==  5:#Incretin枠各条件平均全被験者平均の2分子内相関：Ave6Cond
        #つなげた
        if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
            pass
        elif FastingDatadict['RawDelta']=='Delta':#差分と生を場合分ける
            tempDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Delta/' + SubjectName[0] + '.xlsx',header=0,encoding = "ISO-8859-1")
            NewCounter = [[ SubjectName[0],SubjectName[6],SubjectName[12] ], [ SubjectName[1],SubjectName[7],SubjectName[13] ], [ SubjectName[2],SubjectName[8],SubjectName[14] ], [ SubjectName[3],SubjectName[9],SubjectName[15] ] ,[ SubjectName[4],SubjectName[10],SubjectName[16] ] ,[ SubjectName[5],SubjectName[11],SubjectName[17] ]  ]
            
            for i in range(0,len(NewCounter)):#6条件人分
                Kpc = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Delta/' + NewCounter[i][0] + '.xlsx',header=0,encoding = "ISO-8859-1")                
                if 'Unit' in list(Kpc.index):
                        Kpc = Kpc.drop('Unit')

                imptime = [ -5,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240]
                DelList = list ( set(list(Kpc['time(min)'])) - set( imptime ) ) 
                if 2 in DelList:
                    for jj in DelList:
                        Kpc = Kpc.drop( Kpc.index[Kpc['time(min)'] == jj][0] )
                SubjRPanel = pd.Panel({NewCounter[i][0]:pd.DataFrame(data=None,index=list(Kpc['time(min)']),columns=list(Kpc.columns))})

                for ii in range(len(NewCounter[i])):#3人
                    SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/Delta/' + NewCounter[i][ii] + '.xlsx',header=0,encoding = "ISO-8859-1")
                    if 'Unit' in list(SubjTmCsDf.index):
                        SubjTmCsDf = SubjTmCsDf.drop('Unit')
                    imptime = [ -5,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240]
                    DelList = list ( set(list(SubjTmCsDf['time(min)'])) - set( imptime ) ) 
                    if 2 in DelList:
                        for jj in DelList:
                            SubjTmCsDf = SubjTmCsDf.drop( SubjTmCsDf.index[SubjTmCsDf['time(min)'] == jj][0] )
                    SubjTmCsDf.to_excel(save_dir+NewCounter[i][ii] + '_EachData.xlsx')
                        #SubjTmCsDf = SubjTmCsDf[ SubjTmCsDf[ 'time(min)' ] != DelList]
                    SubjTmCsDf.index = SubjTmCsDf['time(min)']
                    SubjRPanel[ NewCounter[i][ii] ] = SubjTmCsDf
                SubjRmean = SubjRPanel.mean(axis=0)
                SubjRstd = SubjRPanel.std(axis=0, ddof=0)                    
                    #rootDF = pd.concat([ rootDF, SubjTmCsDf ],axis=0, join_axes=[ rootDF.columns ])
                SubjRmean.to_excel(save_dir+NewCounter[i][ii]+'WoDel.xlsx')                    
                RDF,PDF =  AmT.CorrMolEachSubjHelper(SubjRmean,MolDict,Label)
                RDF.to_excel(save_dir + 'RDF'+NewCounter[i][ii] + '.xlsx')
                PDF.to_excel(save_dir + 'PDF'+NewCounter[i][ii] + '.xlsx') 
                
    elif AdjustSwitch ==  6:#Incretin枠各条件平均全条件繋げて2分子内相関
        SubjectName = [ 'exp_SubjectNo3_25_BWoDel', 'exp_SubjectNo3_25_R_2hWoDel', 'exp_SubjectNo3_50_BWoDel', 'exp_SubjectNo3_50_R_2hWoDel', 'exp_SubjectNo3_75_BWoDel', 'exp_SubjectNo3_75_R_2hWoDel' ]
        if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
            pass
        elif FastingDatadict['RawDelta']=='Delta':#差分と生を場合分ける
            #平均タイムコースを取ってくる
            tempDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/TimeCourse/Ave6Cond/' + SubjectName[0] + '.xlsx',header=0,encoding = "ISO-8859-1")
            rootDF = pd.DataFrame(data=None, index= None, columns = ( tempDF.columns ))
            
            for i in range(0,len(SubjectName)):#6条件
                SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/TimeCourse/Ave6Cond/' + SubjectName[0] + '.xlsx',header=0,encoding = "ISO-8859-1")
                    #SubjTmCsDf = SubjTmCsDf[ SubjTmCsDf[ 'time(min)' ] != DelList]
                rootDF = pd.concat([ rootDF, SubjTmCsDf ],axis=0, join_axes=[ rootDF.columns ])
            rootDF.to_excel(save_dir+SubjectName[i]+'WoDel.xlsx')                    
            RDF,PDF =  AmT.CorrMolEachSubjHelper(rootDF,MolDict,Label)
            RDF.to_excel(save_dir + 'RDF'+SubjectName[i] + '.xlsx')
            PDF.to_excel(save_dir + 'PDF'+SubjectName[i] + '.xlsx')
            
    elif AdjustSwitch ==  7: #全被験者繋いで、6条件求める
         if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
             SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/' + NewCounter[i][ii] + '.xlsx',header=0,encoding = "ISO-8859-1")
         
         elif FastingDatadict['RawDelta']=='Delta':#差分と生を場合分ける
            tempDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/' + SubjectName[0] + '.xlsx',header=0,encoding = "ISO-8859-1")
            NewCounter = [[ SubjectName[0],SubjectName[6],SubjectName[12] ], [ SubjectName[1],SubjectName[7],SubjectName[13] ], [ SubjectName[2],SubjectName[8],SubjectName[14] ], [ SubjectName[3],SubjectName[9],SubjectName[15] ] ,[ SubjectName[4],SubjectName[10],SubjectName[16] ] ,[ SubjectName[5],SubjectName[11],SubjectName[17] ]  ]
            for i in range(0,len(NewCounter)):#6条件
                rootDF = pd.DataFrame(data=None, index= None, columns = ( tempDF.columns ))
                for ii in range(len(NewCounter[i])):#3人分
                    SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/TimeCourse/Delta/' + NewCounter[i][ii] + '.xlsx',header=0,encoding = "ISO-8859-1")
                    if 'Unit' in list(SubjTmCsDf.index):
                        SubjTmCsDf = SubjTmCsDf.drop('Unit')
                    imptime = [ -5,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240]
                    DelList = list ( set(list(SubjTmCsDf['time(min)'])) - set( imptime ) ) 
                    if 2 in DelList:
                        for jj in DelList:
                            SubjTmCsDf = SubjTmCsDf.drop( SubjTmCsDf.index[SubjTmCsDf['time(min)'] == jj][0] )
                        #SubjTmCsDf = SubjTmCsDf[ SubjTmCsDf[ 'time(min)' ] != DelList]
                    rootDF = pd.concat([ rootDF, SubjTmCsDf ],axis=0, join_axes=[ rootDF.columns ])
                rootDF.to_excel(save_dir+NewCounter[i][ii]+'WoDel.xlsx')                    
                RDF,PDF =  AmT.CorrMolEachSubjHelper(rootDF,MolDict,Label)
                RDF.to_excel(save_dir + 'RDF'+SubjectName[i] + '.xlsx')
                PDF.to_excel(save_dir + 'PDF'+SubjectName[i] + '.xlsx')    

    elif AdjustSwitch ==  8: #6条件×全被験者繋いだ、
         if FastingDatadict['RawDelta']=='Raw':#差分と生を場合分ける
            pass
         elif FastingDatadict['RawDelta']=='Delta':#差分と生を場合分ける
            tempDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/' + SubjectName[0] + '.xlsx',header=0,encoding = "ISO-8859-1")
            NewCounter = [[ SubjectName[0],SubjectName[6],SubjectName[12] ], [ SubjectName[1],SubjectName[7],SubjectName[13] ], [ SubjectName[2],SubjectName[8],SubjectName[14] ], [ SubjectName[3],SubjectName[9],SubjectName[15] ] ,[ SubjectName[4],SubjectName[10],SubjectName[16] ] ,[ SubjectName[5],SubjectName[11],SubjectName[17] ]  ]
            rootDF = pd.DataFrame(data=None, index= None, columns = ( tempDF.columns ))
            for i in range(0,len(SubjectName)):#6条件
                SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Incretin/' + SubjectName[i] + '.xlsx',header=0,encoding = "ISO-8859-1")
                if 'Unit' in list(SubjTmCsDf.index):
                    SubjTmCsDf = SubjTmCsDf.drop('Unit')
                imptime = [ -5,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240]
                DelList = list ( set(list(SubjTmCsDf['time(min)'])) - set( imptime ) ) 
                if 2 in DelList:
                    for jj in DelList:
                        SubjTmCsDf = SubjTmCsDf.drop( SubjTmCsDf.index[SubjTmCsDf['time(min)'] == jj][0] )
                    #SubjTmCsDf = SubjTmCsDf[ SubjTmCsDf[ 'time(min)' ] != DelList]
                rootDF = pd.concat([ rootDF, SubjTmCsDf ],axis=0, join_axes=[ rootDF.columns ])
            rootDF.to_excel(save_dir+SubjectName[i]+'WoDel.xlsx')                    
            RDF,PDF =  AmT.CorrMolEachSubjHelper(rootDF,MolDict,Label)
            RDF.to_excel(save_dir + 'RDF'+SubjectName[i] + '.xlsx')
            PDF.to_excel(save_dir + 'PDF'+SubjectName[i] + '.xlsx') 
                
                
def AllSubjRPDF_Hist(save_dir):

        SubjTmCsDfNew  = pd.read_excel(save_dir + 'AllSubjTmCsDiff.xlsx')
        RDF,PDF =  AmT.CorrMolEachSubjHelper(SubjTmCsDfNew,MolDict,Label)
        RDF.to_excel(save_dir + 'RDF.xlsx');RDFrev = MatHel.adjustMatrlower(RDF);RDFrev.to_excel(save_dir + 'RDFrev.xlsx');
        PDF.to_excel(save_dir + 'PDF.xlsx');PDFrev = MatHel.adjustMatrlower(PDF);PDFrev.to_excel(save_dir + 'PDFrev.xlsx');   
        RDFEng = CNE.ChnageNametoEnglish(RDF,2);PDFEng = CNE.ChnageNametoEnglish(PDF,2);RDFrevEng = CNE.ChnageNametoEnglish(RDFrev,2);PDFrevEng = CNE.ChnageNametoEnglish(PDFrev,2);
        RDFEng.to_excel(save_dir + 'RDF_Eng.xlsx');PDFEng.to_excel(save_dir + 'PDF_Eng.xlsx');RDFrevEng.to_excel(save_dir + 'RDFrev_Eng.xlsx');PDFrevEng.to_excel(save_dir + 'PDFrev_Eng.xlsx')
        #GH.mkHist(RDFrev,'RDFrev',save_dir);GH.mkHist(PDFrev,'PDFrev'+save_dir)
        
        SwitchDict={}
        SwitchDict['Label']='AllSubj2Mol'#list(PDFrev.columns)
        #SC.UseR(PDFrev,SwitchDict)
        #fig1, ax1 = plt.subplots(1,1)
        fig = plt.figure()
        list1 = np.array(RDFrevEng).flatten()[~np.isnan(RDFrevEng.values.flatten().astype(np.float32))]
        list2 = np.array(PDFrevEng).flatten()[~np.isnan(PDFrevEng.values.flatten().astype(np.float32))]        
        ax1 = plt.hist(list1.astype(np.float32))#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)        
        plt.savefig(save_dir + 'Hist_RDFrev.pdf')  ;plt.close()#print(ConList2[4][~np.isnan(ConList2[4])])
        fig = plt.figure()
        ax1 = plt.hist(list2.astype(np.float32))
        plt.savefig(save_dir + 'Hist_PDFrev.pdf')  ;plt.close()
        SubjTmCsDfNewEng  = pd.read_excel(save_dir + 'AllSubjTmCsDiff_Eng.xlsx')

        return(SubjTmCsDfNewEng,RDFrev,PDFrev)
        
def makExcelAnyDF(DF,QvalueDF,UndrThsh,R,Pvalue,Label,file_dir,save_dir):#DF受け取ってQ値閾値以下の組み合わせのみエクセル表示に
    NewIndex= list(set(QvalueDF.index[UndrThsh[0]]))
    Newcolumns= list(set(QvalueDF.columns[UndrThsh[1]]))
    NewDF = pd.DataFrame(data=[],columns=Newcolumns,index=NewIndex)
    #[R[QvalueDF.columns[UndrThsh[1][i]]][QvalueDF.index[UndrThsh[0][i]]] for i in range(0,len(UndrThsh[0]))]
    for i in range(0, len(NewIndex)):
        ColName = QvalueDF.columns[UndrThsh[1][i]]
        IndName = QvalueDF.index[UndrThsh[0][i]]
        NewDF[ColName][IndName] = R[QvalueDF.index[UndrThsh[1][i]]][QvalueDF.columns[UndrThsh[0][i]]]
    return(NewDF)

def plotScattterAnyDF(DF,QvalueDF,UndrThsh,R,Pvalue,Label,file_dir,save_dir):#被験者ごとプロパティと分子種×時間のスキャッタ_色分け
    SubjPropDF = DF#pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjTimeSeriesProp.xlsx',header=0,encoding = "ISO-8859-1")
    FastPropDF = QvalueDF#pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180326/QvalueFuncQvalueTimeSeriesProp.xlsx',header=0,encoding = "ISO-8859-1")
    
    #UndrThsh = np.where(FastPropDF<0.1)
    #FastPropDF.columns[UndrThsh[0][0]]
    #FastPropDF.index[UndrThsh[1][0]]
    mpl.rcParams['font.family'] = 'Arial Unicode MS'
    #for i in range(0,len(UndrThsh[1])):
     #   SubjPropDF[FastPropDF.columns[UndrThsh[0][i]]],SubjPropDF[FastPropDF.index[UndrThsh[0][i]]]
    num_cols=6#len(DataType)*colmol*avestd #1行あたりのグラフ数
    num_rows=6
    num_page = (len(UndrThsh[1]) // 30) + 1
    #mod_mol = num_mol*len(DataType)%(num_cols*num_rows)
    colorlist = []
    for ii in range(0,14):
        colorlist.append(cm.hsv(ii/14.0))
    count=0
    maxlength = len(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]])
    for k in range(0,num_page+1):
        fig, host = plt.subplots(num_rows,num_cols,figsize=(50,50))
        #plt.rcParams["font.size"] = 20
        plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05)
        for j in range(0,num_rows): #行番号,一行につき1被験者   
            for i in range(0,num_cols):#最終列にlegendつけるならnum_cols-1
                if count < len(UndrThsh[1]) and abs(round(R[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]],3)) > 0.5:
                    y=[[]]*len(BolusNoList)
                    if j < 100:
                        for ii in range(0,int(maxlength/len(BolusNoList))):##iiは-10,0,10のブンだけ回る
                            y[ii] = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]][(ii*20):(ii*20+19)],SubjPropDF[FastPropDF.columns[UndrThsh[0][count]]][(ii*20):(ii*20+19)],s=600,c=colorlist[ii])
                            #fig.legend(np.ravel(y),list(SubjPropDF.index),prop={'size':20},bbox_to_anchor=(0.95, 0.95))
                            r, p=pearsonr(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.columns[UndrThsh[0][count]]])
                            #if r <0:
                             #   y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c=['b']*20)
                            #else:
                             #   y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c='r')
                            host[j,i].set_xlabel(FastPropDF.columns[UndrThsh[1][count]],fontsize=50)                
                            host[j,i].set_ylabel(FastPropDF.columns[UndrThsh[0][count]],fontsize=50)#,fontproperties=prop)                 
                            Roundp = round(p,7)
                            b = '%.2e'%Roundp
                            TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                            p = Pvalue[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]]
                            p='{:e}'.format(p)
                            
                            TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
                            q = FastPropDF[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]]
                            q='{:e}'.format(q)
                            TenRoundq = str(q)[0:4] + '$\it{×10^{-'+str(q)[str(q).find('e')+2:]+'}}$'
                            
                            host[j,i].set_title('R=' + str(round(R[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]],3)) + ', p=' + TenRoundp + '\nq=' + TenRoundq,fontsize=40)
                            host[j,i].tick_params(labelsize=40)
                        
                    else:
                        y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.columns[UndrThsh[0][count]]],s=600,c=colorlist)
                        r, p=pearsonr(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.columns[UndrThsh[0][count]]])
                        #if r <0:
                         #   y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c=['b']*20)
                        #else:
                         #   y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c='r')
                        host[j,i].set_xlabel(FastPropDF.columns[UndrThsh[1][count]],fontsize=50)                
                        host[j,i].set_ylabel(FastPropDF.columns[UndrThsh[0][count]],fontsize=50)#,fontproperties=prop)                 
                        Roundp = round(p,7)
                        b = '%.2e'%Roundp
                        TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                        p = Pvalue[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]]
                        p='{:e}'.format(p)
                        
                        TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
                        q = FastPropDF[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]]
                        q='{:e}'.format(q)
                        TenRoundq = str(q)[0:4] + '$\it{×10^{-'+str(q)[str(q).find('e')+2:]+'}}$'
                        
                        host[j,i].set_title('R=' + str(round(R[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]],3)) + ', p=' + TenRoundp + '\nq=' + TenRoundq,fontsize=40)
                        host[j,i].tick_params(labelsize=40)
                        count+=1
                count+=1
        #fig.legend(np.ravel(y),list(SubjPropDF.index),prop={'size':10},bbox_to_anchor=(0.95, 0.95))
        fig.tight_layout()
        plot_axis = plt.axis()
        plt.savefig(save_dir + str(k) + '.pdf',format='pdf')
    print('output:'+ save_dir)
    
   


def plotScatterBtwTimeSeries(save_dir,Switch,Label,MolDict,PlotSwitch):#空腹値間、-10同士、0同士、10同士,,,
    #SubjPropTimeSeriesDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjTimeSeriesDF.xlsx',header=0,encoding = "ISO-8859-1")#生データ(尖り後)
    SubjPropTimeSeriesDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjTimeSeriesDeltaWONaN.xlsx',header=0,encoding = "ISO-8859-1")#差分データ(尖り後)
    SubjFastDFPre = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjPropFasting.xlsx',header=0,encoding = "ISO-8859-1")
    SubjFastDF = SubjFastDFPre[list(SubjFastDFPre.columns[0:84])]
    mpl.rcParams['font.family'] = 'Arial Unicode MS'
    TimeList = ['-10','0','10','20','30','45','60','75','90','120','150','180','210','240']    
    if Switch == 0:#0は空腹値
        print('空腹値')

    elif 'All' in Switch :
        save_dir = save_dir + 'ScatterBtw_AllTimePoint/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)        

        TargetMol = ['イソロイシン','バリン']
        print(TargetMol[0]+'vs'+TargetMol[1])
        SubjTargetXDF = pd.DataFrame(data=None,columns=[],index=[])#まずDF雛形作る
        SubjTargetYDF = pd.DataFrame(data=None,columns=[],index=[])#まずDF雛形作る
        if Switch == 'AllMolZScore':
            save_dir = save_dir + 'ScatterBtw_AllMolTimePoint/'
            if not os.path.isdir(save_dir):
                os.makedirs(save_dir)    
            SubjTargetTempDF = pd.DataFrame(data=None,columns=list(Label),index=[])#まずDF雛形作る
            #このDFに分子種ごとのデータを入れてく
            Count=0
            xTimeList=[]
            xlist=[]
            NaNDF = pd.DataFrame(data=[np.nan]*20)
            if PlotSwitch == 'ON':
                #SubjPropTimeSeriesDF = pd.DataFrame(data=None)#新しく作り直す
                QvalueDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180506/MolTime/QvalueFuncQvalue.xlsx',header=0,encoding = "ISO-8859-1")
                UndrThsh = np.where(QvalueDF<0.1)
                NeedLabelX = [Label[UndrThsh[0][x]] for x in UndrThsh[0]]
                NeedLabelY = [Label[UndrThsh[1][x]] for x in UndrThsh[1]]

            
        for ij in list(SubjPropTimeSeriesDF.columns):           
            r = re.compile("(.*)(_)(.*)") 
            d = r.search(ij)   
            if d.group(3) == Label[Count]:#条件に合う分子のラベルだったら
                xlist.append(ij);  xTimeList.append(d.group(1))
            if len(xlist) == MolDict[Label[Count]]:#そのラベルの数に達したら、
                difList = list(set(TimeList) - set(xTimeList))
                TimeCount =0 
                templist=[]
                for jj in TimeList:                
                    if jj in difList:#条件に合う分子のラベルだったら
                        templist +=list(NaNDF[0])
                    else:
                        templist += list(sp.zscore(np.array(SubjPropTimeSeriesDF[xlist[TimeCount]]),axis=0))
                        TimeCount += 1
                SubjTargetTempDF[Label[Count]] =  templist     
                Count+=1
                if Count < 83:
                    if Label[Count] == '高感度CRP':
                        Count+=1
                xlist=[]
                xTimeList=[]
        SubjTargetTempDF.index = list(SubjPropTimeSeriesDF.index)*14        
        xlist=[]
        ylist = []
        Count = 0
        for ii in list(SubjPropTimeSeriesDF.columns):
            #re.search(r'[^1]+',XlabelProList)
            r = re.compile("(.*)(_)(.*)") 
            d = r.search(ii)         
            if d.group(3) == TargetMol[0]:#条件に合う分子のラベルだったら
                xlist.append(ii)
            if d.group(3) == TargetMol[1]:
                ylist.append(ii)

        if Switch == 'AllZScore':#ZScore化する
            if len(xlist) == len(ylist):#下から繋げてく
                for j in range(len(xlist)):
                    SubjPropTimeSeriesDF[xlist[j]] = zscore(np.array(SubjPropTimeSeriesDF[xlist[j]]),axis=0)
                    SubjTargetXDF = pd.concat([SubjTargetXDF, SubjPropTimeSeriesDF[xlist[j]]],axis=0)
                    
                for k in range(len(ylist)):
                    SubjPropTimeSeriesDF[ylist[k]] = zscore(np.array(SubjPropTimeSeriesDF[ylist[k]]),axis=0)
                    SubjTargetYDF = pd.concat([SubjTargetYDF, SubjPropTimeSeriesDF[ylist[k]]],axis=0)            
                SubjFastDF = pd.concat([SubjTargetXDF,SubjTargetYDF],axis=1)#新しく作る 
                SubjFastDF.columns = TargetMol            
        else:
            if len(xlist) == len(ylist):#下から繋げてく
                for j in range(len(xlist)):
        
                    SubjTargetXDF = pd.concat([SubjTargetXDF, SubjPropTimeSeriesDF[xlist[j]]],axis=0)
                    
                for k in range(len(ylist)):
                    SubjTargetYDF = pd.concat([SubjTargetYDF, SubjPropTimeSeriesDF[ylist[k]]],axis=0)            
                SubjFastDF = pd.concat([SubjTargetXDF,SubjTargetYDF],axis=1)#新しく作る 
                #SubjFastDF.columns = TargetMol
            
    else:
        save_dir = save_dir + 'ScatterBtw' + str(Switch) +'/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        print(str(Switch)+'分時点')
        TimeIdx = Switch
        newlist=[]
        for ii in list(SubjPropTimeSeriesDF.columns):
            #re.search(r'[^1]+',XlabelProList)
            r = re.compile("(.*)(_)(.*)") 
            d = r.search(ii) 
            if d.group(1) == str(TimeIdx):#条件に合う時点のラベルだったら
                newlist.append(ii)
        SubjFastDF = pd.DataFrame(data=None,columns=newlist,index=list(SubjPropTimeSeriesDF.index))#新しく作る               
        for sk in range(len(newlist)):
            SubjFastDF[newlist[sk]] = SubjPropTimeSeriesDF[newlist[sk]]
            
    ###全分子、全時点間の網羅的相関
    MolTimeRDF = pd.DataFrame(data=[],index=Label,columns=Label)
    MolTimePvalueDF = pd.DataFrame(data=[],index=Label,columns=Label)        
        
    SubjFastDF = SubjTargetTempDF
    #SubjFastDF= SubjFastDF.dropna(how='any',axis=0)
    num_cols=6#len(DataType)*colmol*avestd #1行あたりのグラフ数
    num_rows=6
    num_Mols = len(SubjFastDF.columns)
    num_graph = int(num_Mols*num_Mols/2)
    num_page = 0#(num_graph // (num_cols*(num_rows-1))) + 1
    #mod_mol = num_mol*len(DataType)%(num_cols*num_rows)
    colorlist = []
    for ii in range(0,len(SubjFastDF.index)):
        colorlist.append(cm.hsv(ii/20.0))
    countx=0
    county=countx+1
    #pd.concat([SubjFastDF[SubjFastDF.columns[countx]],SubjFastDF[SubjFastDF.columns[county]]])
    for k in range(0,num_page+1):
        fig, host = plt.subplots(num_rows,num_cols,figsize=(50,50))
        #plt.rcParams["font.size"] = 20
        plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05)
        for j in range(0,num_rows): #行番号,一行につき1被験者   
            for i in range(0,num_cols-1):
                y=[[]]*len(SubjFastDF.index)
                if countx<num_Mols-1:
                    #for ii in range(0,len(SubjFastDF.index)):
                        #y[ii] = host[j,i].scatter(SubjFastDF[SubjFastDF.columns[countx]][ii],SubjFastDF[SubjFastDF.columns[county]][ii],s=600,c=colorlist[ii])
                    #fig.legend(np.ravel(y),list(SubjFastDF.index),prop={'size':20},bbox_to_anchor=(0.95, 0.95))
                    if Switch == 'AllMolZScore':
                      SubjFast4WOnanDF = pd.concat([SubjFastDF[SubjFastDF.columns[countx]],SubjFastDF[SubjFastDF.columns[county]]],axis=1)
                      #SubjFastWOnanDF = SubjFast4WOnanDF.dropna(how='any',axis=0)
                      ##r, p=pearsonr(SubjFastWOnanDF[SubjFastDF.columns[countx]],SubjFastWOnanDF[SubjFastDF.columns[county]])
                      #MolTimeRDF.loc[SubjFastDF.columns[countx],SubjFastDF.columns[county]] = r
                      #MolTimePvalueDF.loc[SubjFastDF.columns[countx],SubjFastDF.columns[county]] = p
                    #else:
                    
                    #r, p=pearsonr(SubjFastDF[SubjFastDF.columns[countx]],SubjFastDF[SubjFastDF.columns[county]])
                    #MolTimeRDF.loc[SubjFastDF.columns[countx],SubjFastDF.columns[county]] = r
                    #MolTimePvalueDF.loc[SubjFastDF.columns[countx],SubjFastDF.columns[county]] = p
                    
                    ##host[j,i].set_xlabel(SubjFastDF.columns[countx],fontsize=50)                
                    ##host[j,i].set_ylabel(SubjFastDF.columns[county],fontsize=50)#,fontproperties=prop)                 
                    #Roundp = round(p,7)
                    #b = '%.2e'%Roundp
                    #TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                    
                    #TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
                    ##host[j,i].set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp,fontsize=40)# + '\nq=' + str(round(SubjFastDF[SubjFastDF.columns[countx]][SubjFastDF.columns[county]],3)),fontsize=40)
                    ##host[j,i].tick_params(labelsize=40)
                    county+=1
                    if county>num_Mols-1:
                        countx+=1
                        county=countx+1    

        #fig.tight_layout()
        #plot_axis = plt.axis()
        #plt.savefig(save_dir + str(k) + '.pdf',format='pdf')
        plt.close()
        
    MolTimeRDF.to_excel(save_dir + 'MolTimeR.xlsx')  
    MolTimePvalueDF.to_excel(save_dir + 'MolTimeP.xlsx') 
    print('output:'+ save_dir)                            
    return(UndrThsh,SubjFastDF,QvalueDF)
    
def plotScattterTimeSeriesvsPropWSubjColor(file_dir,save_dir):#被験者ごとプロパティと分子種×時間のスキャッタ_色分け
    SubjPropDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjTimeSeriesProp.xlsx',header=0,encoding = "ISO-8859-1")
    FastPropDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180326/QvalueFuncQvalueTimeSeriesProp.xlsx',header=0,encoding = "ISO-8859-1")
    UndrThsh = np.where(FastPropDF<0.1)
    FastPropDF.columns[UndrThsh[0][0]]
    FastPropDF.index[UndrThsh[1][0]]
    mpl.rcParams['font.family'] = 'Arial Unicode MS'
    #for i in range(0,len(UndrThsh[1])):
     #   SubjPropDF[FastPropDF.columns[UndrThsh[0][i]]],SubjPropDF[FastPropDF.index[UndrThsh[0][i]]]
    num_cols=6#len(DataType)*colmol*avestd #1行あたりのグラフ数
    num_rows=6
    num_page = (len(UndrThsh[1]) // 30) + 1
    #mod_mol = num_mol*len(DataType)%(num_cols*num_rows)
    colorlist = []
    for ii in range(0,20):
        colorlist.append(cm.hsv(ii/20.0))
    count=0
    for k in range(0,num_page+1):
        fig, host = plt.subplots(num_rows,num_cols,figsize=(50,50))
        #plt.rcParams["font.size"] = 20
        plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05)
        for j in range(0,num_rows): #行番号,一行につき1被験者   
            for i in range(0,num_cols-1):
                if count < len(UndrThsh[1]):
                    y=[[]]*len(BolusNoList)
                    if j == 0:
                        for ii in range(0,len(BolusNoList)):
                            y[ii] = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]][ii],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]][ii],s=600,c=colorlist[ii])
                            fig.legend(np.ravel(y),list(SubjPropDF.index),prop={'size':20},bbox_to_anchor=(0.95, 0.95))
                            r, p=pearsonr(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]])
                            #if r <0:
                             #   y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c=['b']*20)
                            #else:
                             #   y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c='r')
                            host[j,i].set_xlabel(FastPropDF.columns[UndrThsh[1][count]],fontsize=50)                
                            host[j,i].set_ylabel(FastPropDF.index[UndrThsh[0][count]],fontsize=50)#,fontproperties=prop)                 
                            Roundp = round(p,7)
                            b = '%.2e'%Roundp
                            TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                            host[j,i].set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp + '\nq=' + str(round(FastPropDF[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]],3)),fontsize=40)
                            host[j,i].tick_params(labelsize=40)
                        count+=1
                    else:
                        y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c=colorlist)
                        r, p=pearsonr(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]])
                        #if r <0:
                         #   y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c=['b']*20)
                        #else:
                         #   y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c='r')
                        host[j,i].set_xlabel(FastPropDF.columns[UndrThsh[1][count]],fontsize=50)                
                        host[j,i].set_ylabel(FastPropDF.index[UndrThsh[0][count]],fontsize=50)#,fontproperties=prop)                 
                        Roundp = round(p,7)
                        b = '%.2e'%Roundp
                        TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                        host[j,i].set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp + '\nq=' + str(round(FastPropDF[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]],3)),fontsize=40)
                        host[j,i].tick_params(labelsize=40)
                        count+=1
        #fig.legend(np.ravel(y),list(SubjPropDF.index),prop={'size':10},bbox_to_anchor=(0.95, 0.95))
        fig.tight_layout()
        plot_axis = plt.axis()
        plt.savefig(save_dir + str(k) + '.pdf',format='pdf')
    print('output:'+ save_dir)
    
    
    
def plotScattterTimeSeriesvsProp(file_dir,save_dir):#被験者ごとプロパティと分子種×時間のスキャッタ
    SubjPropDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjTimeSeriesProp.xlsx',header=0,encoding = "ISO-8859-1")
    FastPropDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180326/QvalueFuncQvalueTimeSeriesProp.xlsx',header=0,encoding = "ISO-8859-1")
    SubjPropDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/SubjTimeSeriesProp_Eng.xlsx',header=0,encoding = "ISO-8859-1")
    FastPropDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/QvalueFuncQvalueTimeSeriesProp_Eng.xlsx',header=0,encoding = "ISO-8859-1")

    UndrThsh = np.where(FastPropDF<0.1)
    FastPropDF.columns[UndrThsh[0][0]]
    FastPropDF.index[UndrThsh[1][0]]
    mpl.rcParams['font.family'] = 'Arial Unicode MS'
    #for i in range(0,len(UndrThsh[1])):
     #   SubjPropDF[FastPropDF.columns[UndrThsh[0][i]]],SubjPropDF[FastPropDF.index[UndrThsh[0][i]]]
    num_cols=6#len(DataType)*colmol*avestd #1行あたりのグラフ数
    num_rows=6
    num_page = (len(UndrThsh[1]) // 36) + 1
    #mod_mol = num_mol*len(DataType)%(num_cols*num_rows)
    count=0
    for k in range(0,num_page+1):
        fig, host = plt.subplots(num_rows,num_cols,figsize=(50,50))
        #plt.rcParams["font.size"] = 20
        plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05)
        for j in range(0,num_rows): #行番号,一行につき1被験者   
            for i in range(0,num_cols):
                if count < len(UndrThsh[1]):
                    y=[[]]*len(BolusNoList)
                    r, p=pearsonr(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]])
                    if r <0:
                        y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c=['b']*20)
                    else:
                        y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c='r')
                    host[j,i].set_xlabel(FastPropDF.columns[UndrThsh[1][count]],fontsize=50)                
                    host[j,i].set_ylabel(FastPropDF.index[UndrThsh[0][count]],fontsize=50)#,fontproperties=prop)                 
                    Roundp = round(p,7)
                    b = '%.2e'%Roundp
                    TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                    host[j,i].set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp + '\nq=' + str(round(FastPropDF[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]],3)),fontsize=40)
                    host[j,i].tick_params(labelsize=40)
                    count+=1
        fig.tight_layout()
        plot_axis = plt.axis()
        plt.savefig(save_dir + str(k) + '.pdf',format='pdf')
    print('output:'+ save_dir)

#空腹値とプロパティ与えてスキャッター吐き出す
def plotScattterFastvsProp(BolusData,BolusLabel):
    #prop = matplotlib.font_manager.FontProperties(fname=r'C:\Windows\Fonts\meiryo.ttc', size=10)
    #フォントファイルを指定してあげる
    #fp = FontProperties(fname=r'C:\WINDOWS\Fonts\yugothic.ttf', size=14)
    #BolusLabel = BolusData[0,BolusNoList[0]]['Label'][0][0][0]
    #行：被験者、列：空腹値のﾃﾞｰﾀﾌﾚｰﾑを作成
    a=[]
    for i in range(10,30):
        if i == 10:
            a = pd.DataFrame([BolusData[0,10]['RawData_Basal'][0]])
        else:
            series=pd.Series(BolusData[0,i]['RawData_Basal'][0])
            a = a.append(series, ignore_index = True)    
    a.columns= BolusLabel
#aの中から取り出したいときはa.loc[行数、a.columns[ｊ](ｊ：列数)]
    for ii in PropertyList:
        tempseries=[]
        for i in range(10,30):        
            tempseries.append(BolusData[0,i][ii][0][0])        
        series=pd.Series(tempseries)
        #series.columns= [BolusLabeldf.values.flatten()]
        a[ii] =series
        """
#被験者と各プロパティだけまとめる
 SubNameidx = list(SubjectName[0:20].values.flatten())
　#a = pd.DataFrame(index=SubNameidx,columns=[])
 a = pd.DataFrame(data=None,index=[SubNameidx],columns=None)
 count=0
 for ii in PropertyListg:
    #print("'"+ii+"'")
    tempseries=[]
    for i in range(10,30):
        if ii==PropertyListg[0]:
        
            tempseries.append(BolusData[0,i][ii][0][0][0])
        else:
            tempseries.append(BolusData[0,i][ii][0][0])
        #print(tempseries)
    series=pd.Series(tempseries)
    #series.columns= [BolusLabeldf.values.flatten()]
    a[count] =tempseries
    count+=1
a.columns=[PropertyListg]
a
a.to_excel(save_dir + 'human_Property2.xlsx')
        """
    #for k in range(0,len(BolusLabel[0].T)):
        #print(a[a.columns[k]])   
    a.to_excel(save_dir + 'subvsFastProp.xlsx')
    """
    subplot_row = 2
    subplot_col = 9            
    num_mol =len(BolusLabel)
    num_graph = np.empty(num_mol)
    num_page = num_mol//subplot_row 
    count=0
    #mod_nutrient = num_nutrient%2
    #どこかでUrateとCreatinineを飛ばす処理をする。
    for i in range(1,num_page+1):#ページ分まわす
        fig, host = plt.subplots(subplot_row,subplot_col,figsize=(20,5))
        fig.subplots_adjust(left=0.3)
        
        for k in range(0,subplot_row):#行分まわす
            
            p1 = host[k,0].scatter(a[PropertyList[0]],a[a.columns[count]])#Raw
            p2 = host[k,1].scatter(a[PropertyList[1]],a[a.columns[count]])
            p3 = host[k,2].scatter(a[PropertyList[2]],a[a.columns[count]])
            p4 = host[k,3].scatter(a[PropertyList[3]],a[a.columns[count]])
            p5 = host[k,4].scatter(a[PropertyList[4]],a[a.columns[count]])
            p6 = host[k,5].scatter(a[PropertyList[5]],a[a.columns[count]])
            p7 = host[k,6].scatter(a[PropertyList[6]],a[a.columns[count]])
            p8 = host[k,7].scatter(a[PropertyList[7]],a[a.columns[count]])
            p9 = host[k,8].scatter(a[PropertyList[8]],a[a.columns[count]])
            for ii in range(0,len(PropertyList)):                  
                #legend
                host[k-1,ii].set_xlabel(PropertyList[ii])                
                host[k-1,ii].set_ylabel(BolusLabel[count])#,fontproperties=prop)                
                r, p=pearsonr(a[PropertyList[ii]],a[a.columns[count]])
                host[k-1,ii].set_title('R=' + str(round(r,3)) + ', p=' + str(round(p,3)))
                #lines = [p1]
                #plt.title(BolusData[0,10]['Label'][i][0][0], loc='center')
                fig.tight_layout() 
            count+=1
                
        plot_axis = plt.axis()
        #plt.legend(lines, p1.get_label())
        
        plt.savefig(save_dir + str(i),format='pdf')
    print('output:'+ save_dir)
    """
"""    

plt.scatter(a[a.columns[0]],a[PropertyList[0]])
plt.xlabel(PropertyList[0])
plt.ylabel(BolusLabel[0][0][0])


#
series=pd.Series([BolusData[0,10]['RawData_Basal'][0]],index=a.columns)
series
a.append(series)


series
series=pd.Series([BolusData[0,10]['RawData_Basal'][0]])
series
a.append(series)
data_frame.columns
a.columns
"""
def pngmerge(save_dir):
    from PIL import Image as img
    canvas = {}
    count=0
    for i in range(0,6):
        canvas[i] = img.new("RGB",(1440,360*3),(255,255,255))
    
        for ii in range(0,3): #1ページ行分
            #vsProperty
            #try:
            a_jpg = img.open(save_dir +  str(count) + ".pdf") # size(100,100)
            canvas[i] .paste(a_jpg,(0,360*ii))
            count+=1
        canvas[i] .save(save_dir + str(i) + '.pdf', 'pdf', quality=100, optimize=True)
#余白は小さく
        
def mkexcelFastingProp(AllData,Label,BolusNoList,FastingDatadict,save_dir,PropertyList):
    SubjPropDF = pd.DataFrame(data=None,index=list(SubjectName[0]),columns=list(FastingDatadict.keys()))
    for j in range(0,len(Label)):
        SubjPropDF[Label[j]] = FastingDatadict[Label[j]]    
    for jj in range(0,len(PropertyList)):
        SubjPropDF[PropertyList[jj]] = list(itertools.chain.from_iterable((list(itertools.chain.from_iterable(AllData[0,BolusNoList[0]:BolusNoList[19]+1][PropertyList[jj]])))))
        SubjPropDF.to_excel(save_dir + 'SubjPropFasting.xlsx')

    print(save_dir + 'SubjPropFasting.xlsx') 
    
def plotScattterFastvsProp(file_dir,save_dir):#被験者ごとプロパティと空腹値のスキャッタ
    SubjPropDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjPropFasting.xlsx',header=0,encoding = "ISO-8859-1")
    FastPropDF = pd.read_excel(file_dir + '/QvalueFuncQvalue.xlsx',header=0,encoding = "ISO-8859-1")
    UndrThsh = np.where(FastPropDF<0.1)
    FastPropDF.columns[UndrThsh[0][0]]
    FastPropDF.index[UndrThsh[1][0]]
    mpl.rcParams['font.family'] = 'Arial Unicode MS'
    #for i in range(0,len(UndrThsh[1])):
     #   SubjPropDF[FastPropDF.columns[UndrThsh[0][i]]],SubjPropDF[FastPropDF.index[UndrThsh[0][i]]]
    num_cols=6#len(DataType)*colmol*avestd #1行あたりのグラフ数
    num_rows=6
    num_page = 1
    #mod_mol = num_mol*len(DataType)%(num_cols*num_rows)
    count=0
    for k in range(0,num_page+1):
        fig, host = plt.subplots(num_rows,num_cols,figsize=(50,50))
        #plt.rcParams["font.size"] = 20
        plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05)
        for j in range(0,num_rows): #行番号,一行につき1被験者   
            for i in range(0,num_cols):
                if count < len(UndrThsh[1]):
                    y=[[]]*len(BolusNoList)
                    r, p=pearsonr(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]])
                    if r <0:
                        y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c='b')
                    else:
                        y = host[j,i].scatter(SubjPropDF[FastPropDF.columns[UndrThsh[1][count]]],SubjPropDF[FastPropDF.index[UndrThsh[0][count]]],s=600,c='r')
                    host[j,i].set_xlabel(FastPropDF.columns[UndrThsh[1][count]],fontsize=50)                
                    host[j,i].set_ylabel(FastPropDF.index[UndrThsh[0][count]],fontsize=50)#,fontproperties=prop)                 
                    Roundp = round(p,7)
                    b = '%.2e'%Roundp
                    TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                    host[j,i].set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp + '\nq=' + str(round(FastPropDF[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]],3)),fontsize=40)
                    host[j,i].tick_params(labelsize=40)
                    count+=1
        fig.tight_layout()
        plot_axis = plt.axis()
        plt.savefig(save_dir + str(k),format='pdf')
    print('output:'+ save_dir)


def plotScatterProp(save_dir):
    import numpy as np
    import pandas as pd
    SubjectPropertyDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjPropWgender.xlsx',header=0,encoding = "ISO-8859-1")
#行方向に10こ、プロパティの数分、列も同じk。列に数は１、２、３、、、、
    plt.rcParams["font.size"] = 10
    colnames = SubjectPropertyDF.columns
    index = SubjectPropertyDF.index
    CorrDF = pd.DataFrame(data=None,index=list(SubjectPropertyDF.columns),columns=list(SubjectPropertyDF.columns))
    PDF = pd.DataFrame(data=None,index=list(SubjectPropertyDF.columns),columns=list(SubjectPropertyDF.columns))

    fig, host = plt.subplots(len(colnames),len(colnames),figsize=(20,20))
    #fig.subplots_adjust(left=0.3)      
    for i in range(0,len(colnames)):#列
        for j in range(0+i,len(colnames)):#行
            if i==j:#同じプロパティのところではヒストグラム
                host[j,i].hist(SubjectPropertyDF[colnames[j]])
                if i == 0: 
                   host[j,0].set_ylabel(colnames[j],fontsize=20)  
                host[9,i].set_xlabel(colnames[i],fontsize=20)
                host[9,i].xaxis.set_label_coords(0.5, -0.3)
                host[j,i].tick_params(labelsize=15)
            #cor=np.corrcoef(PropValue[i],PropValue[1+j])
            #corr[j,i] = cor[0,1]
            else:
                r, p = pearsonr(SubjectPropertyDF[colnames[j]],SubjectPropertyDF[colnames[i]])
                host[j,i].scatter(SubjectPropertyDF[colnames[i]],SubjectPropertyDF[colnames[j]],s=70)
                host[j,i].tick_params(labelsize=15)
                if i == 0: 
                    host[j,0].set_ylabel(colnames[j],fontsize=20)
                host[9,i].set_xlabel(colnames[i],fontsize=20)
                Roundp = round(p,7)
                b = '%.2e'%Roundp
                TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                if colnames[j]=='Gender' or colnames[i]=='Gender':
                    if p<0.05:
                        host[j,i].set_title('p=' + TenRoundp, color='red',fontsize=12) #"%.6f" %(p), '%e' %p
                    else:
                        host[j,i].set_title('p=' + TenRoundp,fontsize=12)
                        #host[j,i].set_title("Drain Current $\it{V_{DS}}$ (mA)",fontsize=10)
                else:
                    if p<0.05:
                        host[j,i].set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp, color='red',fontsize=12)
                    else:
                        host[j,i].set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp,fontsize=12)
                CorrDF.loc[colnames[i],colnames[j]] = r
                PDF.loc[colnames[i],colnames[j]] = p
    CorrDF.to_excel(save_dir+'CorrDF_Prop.xlsx')
    PDF.to_excel(save_dir+'PDF_Prop.xlsx')
    
    fig.tight_layout()
    #plt.tick_params(labelbottom='off')
    plt.savefig(save_dir + '4.pdf')
def _main():
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    sns.set_style('whitegrid')
    #matplotlib inline
    SubjectPropertyDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjPropWgender.xlsx',header=0,encoding = "ISO-8859-1")
    sns.pairplot(SubjectPropertyDF)
    plt.suptitle('This is the main title')

    plt.show()
#plt.savefig(filename, bbox_inches='tight')
#plotScattterFastvsPropSumitomo(BolusData,BolusLabel)
#ADjustAllData用にLabel,BolusNoList,NaNTogariLabel,AllPointLabel,ThHldCondition,DataTypeを作るスクリプト
"""
Label,NaNTogariLabel,AllPointLabel,ThHldCondition,DataType,AllLabel = aD.PreProcess(AllData,BolusLabel)
AllData,dupdict,DirLabel,FileLabel = aD.AdjustAllData(AllData,Label,AllLabel,BolusNoList,NaNTogariLabel,AllPointLabel,ThHldCondition,DataType)
FastingDatadict = aD.mkFastingData(AllData,dupdict,Label,DirLabel,FileLabel,BolusNoList)
print(PropertyListg)
"""
if __name__ == '__main__':
    
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.io
    import pandas as pd
    from Helper import general as ge
    from DataLoader import mkBolusIdx as BI
    from pylab import *
    import matplotlib.font_manager
    from scipy.stats import pearsonr
    import itertools
    from DataLoader import AdjustData011518 as aD
    import matplotlib.cm as cm
    import re
    import MolTimeSeriesCV as MTSCV
    import scipy.stats as sp
    import Helper.AnalMolTime as AmT
    import MatrHelper as MatHel
    import collections
    import Helper.GraphHelper as GH
    import Helper.MolCorrSubjmeanStdHelper as MSMSH
    import Helper.ISHelper as ISH
    from scipy.stats import zscore
    import DataLoader.ChangeNameEng as CNE
    import Helper.AnalSubjTmCsHelper as ASTH
    import Helper.StatCal as SC
    import Helper.AnalPropParamHelper as APPHel


    #########諸々スイッチ
    ColorSwitch = 3#平均と標準偏差の散布図0つけない、1ー色つける 2でAveを絶対値に3でAveを絶対値、同一代謝に色付け   (基本的に3)
    ModelParamSwitch = 0#1でモデルパラメータの相関？
    MolCorrSubjmeanStdSwitch = 1#0で何もしない、1で被験者内2分子網羅的相関係数解析、2で全被験者2分子網羅的相関係数解析なら
    AdjustSwitch = 0 #ある時点に絞る
    if MolCorrSubjmeanStdSwitch == 1 or MolCorrSubjmeanStdSwitch == 2:
        NormlSwitch = 0#0でヒストグラム描画時に各代謝系ごとに正規化しない、1:各代謝ごとに総数で正規化して累積分布 2で累積分布、3で何もしない(1か2が基本)
    
    InsSensSwitch = 0#インスリンに関する感受性checkするなら、  
    NmlzSwitch = 0#1で列方向に標準化
    RowDeltaSwitch=0# 0でRaw, 1でDelta

    
     ###########################################3データ取得
   
    Today = ge.getToday()
    #macとWinで変える
    data_dir ='/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/' 
    #data_dir = 'C:/Users/fujita/Google ドライブ/Kuroda lab c = collections.Counter(ColorList)/Research/Metabolome/file'
    ref_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/To_analize/'
    #ref_dir = 'C:/Users/owner/Documents/fujita/Data/Metabolome/file/To_analize/'
    #name_relation = pd.read_excel(ref_dir + 'LSI測定代謝物_名称対応表.xlsx',header=None,encoding = "ISO-8859-1")
    #macとWinで変える
    matBolus = scipy.io.loadmat('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/AllData_BolusOnly_WoUC_UsableTimecourse.mat')
    #matBolus = scipy.io.loadmat('C:/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/AllData_BolusOnly_WoUC_UsableTimecourse.mat')
    #AllData = mat['AllData']
    BolusData = matBolus['AllData']
    mat = scipy.io.loadmat('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/AllData.mat')
    
    AllData = mat['AllData']
    
    BolusNoListPre = itertools.chain.from_iterable(BI.mkBolusData(BolusData, data_dir, ref_dir)-1)
    BolusNoList=[]
    for i in BolusNoListPre:
        BolusNoList.append(i[0])
        
    AllLabelPre = AllData[0,BolusNoList[0]]['Label'][0]
    AllLabel=[]
    for i in AllLabelPre:
        AllLabel.append(i[0])
    DataType = ['RawData', 'DeltaData']
    BolusLabel = BolusData[0,BolusNoList[0]]['Label'][0]
    BolusLabel = list(itertools.chain.from_iterable(BolusLabel))
    plt.rc("font", family="Arial") #ギリシャ文字など出力のためフォント指定
    Timepoint = pd.read_excel(data_dir + '/Timepoint.xls',header=None,encoding = "ISO-8859-1")
    
    SubjectName = pd.read_excel(data_dir + '/SubjectNameEng.xls',header=None,encoding = "ISO-8859-1")

    SubjTimeSeriesDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjTimeSeriesDF.xlsx',header=0,encoding = "ISO-8859-1")#生データ(尖り後)
    CorrDF = pd.read_excel("/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Rmean_Eng.xlsx",header=0,encoding = "ISO-8859-1")

    save_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/' + Today  + '/'

    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
        
    Label,NaNTogariLabel,AllPointLabel,ThHldCondition,DataType = aD.PreProcess(AllData,BolusLabel)
    #プロパティをまとめる(gender以外9種)
    #閾値以上の尖りをNanにする
    AllData,dupdict,DirLabel,FileLabel = aD.AdjustAllData(AllData,Label,AllLabel,BolusNoList,NaNTogariLabel,AllPointLabel,ThHldCondition,DataType,SubjectName)
    #DealtaDataを作成しAllDataに入れる    
    
    FastingDatadict = aD.mkFastingData(AllData,dupdict,Label,DirLabel,FileLabel,BolusNoList,Timepoint)
    #AllData = aD.mkDeltaData(AllData,dupdict,Label,DirLabel,FileLabel,BolusNoList,FastingDatadict)

    SubjectName = list(pd.read_excel(data_dir + '/SubjectNameEng.xls',header=None,encoding = "ISO-8859-1")[0])

    PropertyList = ['BMI', 'Age', 'Height', 'Weight', 'Waist', 'FastingGlucose', 'TwoHGlucose', 'FastingInsulin', 'TwoHInsulin'];
    #BolusPropertyList = ['BMI', 'Age', 'Height', 'Weight', 'Waist', 'FastingGlucose', 'TwoHGlucose', 'FastingInsulin', 'TwoHInsulin'];
    
    PropertyListg = ['Gender', 'BMI', 'Age', 'Height', 'Weight', 'AC', 'FastingGlucose', 'TwoHGlucose', 'FastingInsulin', 'TwoHInsulin'];
    numPropertyList = len(PropertyList);
    plotScatterProp(save_dir)#プロパティ同士の相関をスキャッター

    #mkexcelFastingProp(AllData,Label,BolusNoList,FastingDatadict,save_dir,PropertyListg)
    file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180125/Exception/'
    #plotScattterFastvsProp(file_dir,save_dir)#被験者ごとプロパティと空腹値のスキャッター
    #plotScattterTimeSeriesvsProp(file_dir,save_dir)#被験者ごとプロパティと分子種×時間のスキャッタ
    #plotScattterTimeSeriesvsPropWSubjColor(file_dir,save_dir)#被験者ごとプロパティと分子種×時間のスキャッタ被験者ごと色分け

    
    MolDict, MolCVAveDict,MolCVStdDict,MolCVCVDict, SubjTimeSeriesCVDF,MolCVVarDict= MTSCV.CalcCV(SubjTimeSeriesDF,Label)#各時点のCVを求めている（今欲しいのではない）
    TimeList = ['AllMolZScore']#-10,0,10..時点で指定できる、or'All'で全時点（関数の中で分子種を設定できる）、'AllZScore'で各時点でZscore化したもの同士の相関
    
    for Switch in TimeList:
        PlotSwitch = 'ON'#'ON'
        #UndrThsh,SubjFastDF,QvalueDF = plotScatterBtwTimeSeries(save_dir,Switch,Label,MolDict,PlotSwitch)#被験者間空腹値の網羅的相関
    #SubjFastDFPre = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjPropFasting.xlsx',header=0,encoding = "ISO-8859-1")
    #SubjFastDF = SubjFastDFPre[list(SubjFastDFPre.columns[0:84])]
    #QvalueDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180508/QvalueFuncQvalue.xlsx',header=0,encoding = "ISO-8859-1")
    #UndrThsh = np.where(QvalueDF<0.1)
    
    ###差分の時
    #Pvalue = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180506/MolTime/ref/MolTimeP.xlsx',header=0,encoding = "ISO-8859-1")
    #R= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180506/MolTime/ref/MolTimeR.xlsx',header=0,encoding = "ISO-8859-1")
    ####空腹値のとき
    #Pvalue = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180508/MolTimeP.xlsx',header=0,encoding = "ISO-8859-1")
    #R = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180508/MolTimeR.xlsx',header=0,encoding = "ISO-8859-1")

    
    if ModelParamSwitch == 1:#1でモデルパラメータの相関？
        NewLabel = ge.ConvertParamLabelRtoPython(Label)
        SubjFastDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/fig/ParamTauAnal_20180511Subject_All/ParamACDE_SubjMollog10.xlsx',header=0,encoding = "ISO-8859-1")
        SubjFastDF.columns = NewLabel + PropertyList
        QvalueDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180516/ParamACDE_ParamACDE/log10/QvalueFuncQvalueParamACDE_ParamACDE.xlsx',header=0,encoding = "ISO-8859-1")
        QvalueDF.columns = NewLabel; QvalueDF.index = NewLabel
        UndrThsh = np.where(QvalueDF<0.1)
        Pvalue = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180516/ParamACDE_ParamACDE/log10/PvalueaaParamACDE_ParamACDE.xlsx',header=0,encoding = "ISO-8859-1")
        Pvalue.columns = NewLabel; Pvalue.index = NewLabel
        R= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180516/ParamACDE_ParamACDE/log10/CorrParamACDE_ParamACDE.xlsx',header=0,encoding = "ISO-8859-1")
        R.columns = NewLabel; R.index = NewLabel
        NewDF = makExcelAnyDF(SubjFastDF,QvalueDF,UndrThsh,R,Pvalue,Label,file_dir,save_dir)
        NewDF.to_excel(save_dir + 'RBetwModelPAramByFDRQ.xlsx')
        
    #plotScattterAnyDF(SubjFastDF,QvalueDF,UndrThsh,R,Pvalue,Label,file_dir,save_dir)#Q値閾値以下のものだけ描画
    
    
    if MolCorrSubjmeanStdSwitch == 1:#被験者内2分子網羅的相関係数解析なら
        ####リスト2つにすればピアソンの相関係数を求められる
        #1.被験者一人選ぶ、分子選ぶ、全組み合わせ相関計算する
        #CalcCorrBetw2MolEachSubj(SubjTimeSeriesDF,SubjectName,MolDict,Label,AdjustSwitch,FastingDatadict,save_dir)
        #2.被験者20人の同じ組み合わせのRの平均、分散
        SubjRPanel,SubjRmean,SubjRstd = mkPneltoCalcAveVar(AdjustSwitch,FastingDatadict,SubjectName,save_dir)
        
        #SubjRmeanrev = MatHel.adjustMatrupper(SubjRmean); SubjRstdrev = MatHel.adjustMatrupper(SubjRstd)#対角をNaNに、上三角をNaNに
        SubjRmeanrev = MatHel.adjustMatrlower(SubjRmean); SubjRstdrev = MatHel.adjustMatrlower(SubjRstd)#対角をNaNに、下三角をNaNに
        SubjRmeanrev = CNE.ChnageNametoEnglish(SubjRmeanrev,2); SubjRstdrev = CNE.ChnageNametoEnglish(SubjRstdrev,2)
        #PlotScatter(save_dir,SubjRmean,SubjRstd,ColorSwitch,NormlSwitch,[])
    
        SubjRmeanrev.to_excel(save_dir + 'RBetwModelPAr.xlsx')
        #GH.mkHist(SubjRmeanrev,save_dir)
        
        #CombDF = MSMSH.mkExcelupper(SubjRmeanrev,SubjRstdrev,0.8)#ある値以上の分子名組み合わせをエクセルにまとめる           
        #CombDF.to_excel(save_dir + 'CombDf.xlsx')
        
        CombDF = MSMSH.mkExcellower(SubjRmeanrev,SubjRstdrev,0.8)#DF入れてある値以下の行、列名取ってきてエクセルにする
        CombDF.to_excel(save_dir + 'CombDf_Corr<0.8.xlsx')  
        
        fig = plt.figure()
        fig.set_size_inches(np.array(CorrDF).shape[1]/5, np.array(CorrDF).shape[0]/6)
        ASTH.heatmap(np.array(CorrDF),list(CorrDF.columns))
        fig.tight_layout() 
        plt.savefig(save_dir+'HeatMap.pdf')
        
        #被験者20人で全部の解析をやる
        #EachSubjAnal(SubjRPanel,SubjectName)
        
    if MolCorrSubjmeanStdSwitch == 2:#全被験者2分子網羅的相関係数解析なら
        CalcCorrBetw2MolEachSubj(SubjTimeSeriesDF,SubjectName,MolDict,Label,AdjustSwitch,FastingDatadict,save_dir)
        SubjTmCsDfNewEng,RDFrev,PDFrev = AllSubjRPDF_Hist(save_dir)
        #Qval = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180826/ParamACDE_ParamACDE__AICClst7_pearson/QvalueStoreyParamACDE_ParamACDE.xlsx')
        #CombDF = APPHel.mkExcellower2(RDFrev,RDFrev,Qval,0.1)#DF入れてある値未満の行、列名取ってきてエクセルにする
        #CombDF.to_excel(save_dir + 'CombDf_Qval<0.1.xlsx')
        #20180827
        #CombDF = MSMSH.mkExcelupper(RDFrev,RDFrev,0.8)#ある値以上の分子名組み合わせをエクセルにまとめる           
        #CombDF.to_excel(save_dir + 'CombDf_Corr>=0.8.xlsx')   

        CombDF = MSMSH.mkExcellower(RDFrev,RDFrev,0.8)#DF入れてある値以下の行、列名取ってきてエクセルにする
        CombDF.to_excel(save_dir + 'CombDf_Corr<0.8.xlsx')   
            
        SwitchDict={};SwitchDict['EngLabel'] = []; SwitchDict['Time'] = [];#SwitchDict['Target']=TargetComb; SwitchDict['Label']=LabelUniq;SwitchDict['MolType'] = 'Param'#str(TimepointList[i])
        #SwitchDict['method'] = 'pearson'# 'pearson' or 'spearman'
        #Corr, Pvalue, save_dir = ASTH.TmSrCorr(TryList[ii],SwitchDict,save_dir)
        #SC.UseR(Pvalue,SwitchDict)
        
        #RDFrev = RDFrev.astype(np.float32);SwitchDict={};SwitchDict['MolColorScatter']=1
        PlotScatter(save_dir,RDFrev,RDFrev,ColorSwitch,NormlSwitch,SwitchDict)
        
        #APPHel.PlotScatterDF(SubjTmCsDfNewEng,CombDF,save_dir ,3)#1で全部赤、3で#散布図：生を赤に、負を青に
        #ASTH.PlotScatter(SwitchDict,save_dir) 
       
    if InsSensSwitch == 1:#インスリンに関する感受性checkするなら、
        NmlzLabel=[]
        SensIdx = 'インスリン'#何に対する感受性か？
        for ii in SubjectName:
            tempsubjTmCs = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw/SubjTmCs_'+ ii + 'Raw.xlsx',header=0,encoding = "ISO-8859-1")
            IdxList = tempsubjTmCs[SensIdx];
            if NmlzSwitch == 1:#1で列方向に標準化
                NmlzLabel+=['Zscore']
                IdxList = zscore(tempsubjTmCs[SensIdx])
                tempsubjTmCs = tempsubjTmCs.drop(SensIdx,axis=1)
                tempsubjTmCs = tempsubjTmCs.apply(sp.zscore, axis=0)
            else:
                NmlzLabel=['Normal']
            NmlzLabel+=['Raw'] if RowDeltaSwitch==0 else ['Delta']# 0でRaw, 1でDelta

                            
            ISH.InsulinSensitivity(tempsubjTmCs,IdxList,SensIdx,NmlzLabel,save_dir+ii+'/')
            
        
    #pngmerge(save_dir)
"""    By20181106
    if __name__ == '__main__':
        _main()
                import matplotlib.pyplot as plt

# Start with one
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot([1,2,3])

# Now later you get a new subplot; change the geometry of the existing
n = len(fig.axes)
for i in range(n):
    fig.axes[i].change_geometry(n+1, 1, i+1)

# Add the new
ax = fig.add_subplot(n+1, 1, n+1)
ax.plot([4,5,6])

plt.show() 
"""