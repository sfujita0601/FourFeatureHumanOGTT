# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 11:00:04 2017

@author: fujita

"""
    
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import pandas as pd
from pylab import *
import matplotlib.font_manager
from scipy.stats import pearsonr
import itertools
import matplotlib.cm as cm
import re
import scipy.stats as sp
import collections
from scipy.stats import zscore

def mksixList(SubjRmean,Col,CnumList):
    tempnp1 = SubjRmean[Col[CnumList[0]]:Col[CnumList[1]-1]].values.flatten(); 
    tempnp4 = SubjRmean[Col[CnumList[3]]:Col[CnumList[4]-1]].values.flatten(); 
    tempnp2 = SubjRmean[Col[CnumList[1]]:Col[CnumList[2]-1]].values.flatten(); 
    tempnp3 = SubjRmean[Col[CnumList[2]]:Col[CnumList[3]-1]].values.flatten();  
    tempnp5 = SubjRmean[Col[CnumList[4]]:Col[CnumList[5]-1]].values.flatten(); 
    tempnp6 = SubjRmean[Col[CnumList[5]]:Col[CnumList[6]-1]].values.flatten()  
    return(tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6)

def mksixgrayList(SubjRmean,Col,CnumList):
    list1=[];list2=[];list3=[];list4=[];list5=[];list6=[];list7=[]
    listr1=[];listr2=[];listm1=[];listm2=[];listg1=[];listg2=[];listb1=[];listb2=[];listbl1=[];listpl1=[];listpl2=[];listbl2=[];listgr1=[];listgr2=[]
    listrabs1=[];listrabs2=[];listmabs1=[];listmabs2=[];listgabs1=[];listgabs2=[];listbabs1=[];listbabs2=[];listblabs1=[];listplabs1=[];listplabs2=[];listblabs2=[];listgrabs1=[];listgrabs2=[]
    listposr1=[];listposr2=[];listposm1=[];listposm2=[];listposg1=[];listposg2=[];listposb1=[];listposb2=[];listposbl1=[];listpospl1=[];listpospl2=[];listposbl2=[];listposgr1=[];listposgr2=[]
    for i in range(0,len(SubjRmean.columns)):
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

def PlotMolList(Col,c,MolColor):
    alist = [[0] for i in range(len(c)+1)]
    for i in range(len(alist)):
        if i==0:
            alist[i] = 0
        else:
            alist[i] = c[MolColor[i-1]] + alist[i-1]
    return(alist)
def plothistEachMetabo(SubjRmean,Col,c,MolColor,save_dir):
            CnumList = PlotMolList(Col,c,MolColor)
            fig2, ax2 = plt.subplots(1,1)
            if len(CnumList) == 7:
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6 = mksixabsList(SubjRmean,Col,CnumList)
                ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=MolColor)        
     
                plt.close()

         
                fig2, ax2 = plt.subplots(1,1)
                tempDF = pd.DataFrame(data=None)
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6 = mksixList(SubjRmean,Col,CnumList)
                ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)])], stacked=True, bins = 20,color=MolColor)        
                plt.close()

                fig4, ax4 = plt.subplots(1,1)
                ConList1 = mksixgrayList(SubjRmean,Col,CnumList)
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         #ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
                MolColorGray = MolColor+['gray']
                ax4 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=MolColorGray)        
                plt.close()    
                fig6, ax6 = plt.subplots(1,1)
                ax6 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])]], stacked=True, bins = 20,color=MolColor)        
                plt.close()        
                fig8, ax8 = plt.subplots(1,1)
                ConList1 = [[abs(ConList1[j][i]) for i in range(len(ConList1[j]))] for j in range(len(ConList1))]
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         #ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]       
                ax8 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=MolColorGray)        
                plt.savefig(save_dir + 'Fig6_b1.pdf') 
                plt.close()
                fig10, ax10 = plt.subplots(1,1)
                ax10 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])]], stacked=True, bins = 20,color=MolColor)        
                xmin, xmax, ymin, ymax = plt.axis()
                plt.ylim([ymin,ymax])
                plt.savefig(save_dir + 'Fig6_b2.pdf')             
                plt.close()

                ax8 = plt.gca()
                ax8.spines['bottom'].set_linewidth(1.5); ax8.spines['left'].set_linewidth(1.5); ax8.spines['right'].set_linewidth(1.5); ax8.spines['top'].set_linewidth(1.5); 
                plt.close()

            elif len(CnumList) == 8:
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6,tempnp7 = mksevenabsList(SubjRmean,Col,CnumList)
                ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)]),list(tempnp7[~np.isnan(tempnp7)])], stacked=True, bins = 20,color=MolColor)                
                plt.close()
                fig2, ax2 = plt.subplots(1,1)
                tempDF = pd.DataFrame(data=None)
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6,tempnp7 = mksevenList(SubjRmean,Col,CnumList)
                ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)]),list(tempnp6[~np.isnan(tempnp6)]),list(tempnp7[~np.isnan(tempnp7)])], stacked=True, bins = 20,color=MolColor)        
                plt.close()                
                fig4, ax4 = plt.subplots(1,1)
                ConList1 = mksevengrayList(SubjRmean,Col,CnumList)
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         #ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
                MolColorGray = MolColor+['gray']
                ax4 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])],ConList1[7][~np.isnan(ConList1[7])]], stacked=True, bins = 20,color=MolColorGray)        
                plt.close()
                fig6, ax6 = plt.subplots(1,1)
                ax6 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=MolColor)        
                plt.close()
            
                fig8, ax8 = plt.subplots(1,1)
                ConList1 = [[abs(ConList1[j][i]) for i in range(len(ConList1[j]))] for j in range(len(ConList1))]
                
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         #ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
        
                ax8 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])],ConList1[7][~np.isnan(ConList1[7])]], stacked=True, bins = 20,color=MolColorGray)        
                plt.savefig(save_dir + 'Fig6_b1.pdf') 
                plt.close()
                
                fig10, ax10 = plt.subplots(1,1)
                ax10 = plt.hist([ConList1[0][~np.isnan(ConList1[0])],ConList1[1][~np.isnan(ConList1[1])],ConList1[2][~np.isnan(ConList1[2])],ConList1[3][~np.isnan(ConList1[3])],ConList1[4][~np.isnan(ConList1[4])],ConList1[5][~np.isnan(ConList1[5])],ConList1[6][~np.isnan(ConList1[6])]], stacked=True, bins = 20,color=MolColor)        
    
                xmin, xmax, ymin, ymax = plt.axis()
                plt.ylim([ymin,ymax])
                plt.savefig(save_dir + 'Fig6_b2.pdf')  
                plt.close()                
def plotculmitivehistEachMetabo(SubjRmean,Col,c,MolColor,save_dir):
            msize=3; lwidth=1.5
            CnumList = PlotMolList(Col,c,MolColor)
            if len(CnumList) == 7:
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6 = mksixabsList(SubjRmean,Col,CnumList)
                ConList1 = mksixgrayList(SubjRmean,Col,CnumList)

                ConList1 = [[abs(ConList1[j][i]) for i in range(len(ConList1[j]))] for j in range(len(ConList1))]
                ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];
                fig7, ax7 = plt.subplots(1,1)
                x = np.sort(ConList1[0][~np.isnan(ConList1[0])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0)

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

                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[3])
    
                x = np.sort(ConList1[4][~np.isnan(ConList1[4])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0);

                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)

                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[4])
                
                x = np.sort(ConList1[5][~np.isnan(ConList1[5])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0);
                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c=MolColor[5])
                
                x = np.sort(ConList1[6][~np.isnan(ConList1[6])]); y=np.arange(1,len(x)+1) / len(x)
                x=np.insert(x,0,0);
                xnew = [[x[i]]+[x[i]] for i in range(len(x)-1)]               
                y = [[y[i]]+[y[i]] for i in range(len(y))]
                xnew = [e for inner_list in xnew for e in inner_list]
                y = [e for inner_list in y for e in inner_list]
                xnew.append(x[-1])
                x=xnew
                y.insert(0,0.0)
                ax7 = plt.plot(x,y,marker='.',markersize=3,linestyle='-',linewidth = 1,c='gray')
                plt.close()
            
            elif len(CnumList) == 8:
                tempnp1,tempnp2,tempnp3,tempnp4,tempnp5,tempnp6,tempnp7 = mksevenabsList(SubjRmean,Col,CnumList)            
                ConList1 = mksevengrayList(SubjRmean,Col,CnumList)
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
                plt.close()
###         
def PlotScatter(save_dir,SubjRmean,SubjRstd,ColorSwitch,NormlSwitch,SwitchDict):
    Ave = list(SubjRmean.values.flatten()[~np.isnan(SubjRmean.values.flatten().astype(float))]);        
    AbsAve = list(abs(SubjRmean.values.flatten()[~np.isnan(SubjRmean.values.flatten().astype(float))]))
    AbsAvemin = min(AbsAve); AbsAvemax = max(AbsAve)
    Std = list(SubjRstd.values.flatten()[~np.isnan(SubjRstd.values.flatten().astype(float))])
    Stdmin = min(Std); Stdmax = max(Std)
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams["font.size"] = 20

    r, p=pearsonr(Ave,Std)
    p='{:e}'.format(p)

    TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
    if ColorSwitch == 0:

        fig, ax = plt.subplots(1,1)
        ax.scatter(Ave,Std)
        SaveName = []
        SaveName='_WNmlColor'
        fig.savefig(save_dir +'AveStdScatter' + SaveName +'.pdf',format='pdf')
    
    elif ColorSwitch == 1:
        MolColorDF = SwitchDict['MolColor']
        fig, ax = plt.subplots(1,1)
        ColorList =[]
        Col = SubjRmean.columns
        for i in range(0,len(SubjRmean.columns)):
            list1=list(SubjRmean.loc[Col[i],Col[i]:]);list2 =list(SubjRstd.loc[Col[i],Col[i]:])
            list1.remove(1);list2.remove(0);
            
            ax.scatter(list1,list2,color='Blue')     
        SaveName='_WBlueColor'
        fig.savefig(save_dir +'AveStdScatter' + SaveName +'.pdf',format='pdf')
        
    elif ColorSwitch == 2:
        MolColorDF = SwitchDict['MolColor']
        fig, ax = plt.subplots(1,1)
        ColorList = []
        Col = SubjRmean.columns
        for i in range(0,len(SubjRmean.columns)):
            list1=list(abs(SubjRmean.loc[Col[i],Col[i]:]));list2 =list(SubjRstd.loc[Col[i],Col[i]:])
            list1.remove(1);list2.remove(0);
            
            ax.scatter(list1,list2,color='blue')
        SaveName='_WBlueColor'
        fig.savefig(save_dir +'abs(Ave)StdScatter' + SaveName +'.pdf',format='pdf')
            
    elif ColorSwitch == 3:#
        MolColorDF = SwitchDict['MolColor']
        
        fig1, ax1 = plt.subplots(1,1)#abs
        fig2, ax2 = plt.subplots(1,1)#normal

        Col = list(SubjRmean.columns)
        try:
            MolColor = SwitchDict['MolColor']
            if SwitchDict['Check'] == 'Amino':
                MolColorDF,a = LH.AACheck(MolColor,MolColor,'EAA')#'protein','ketogenic','EAA','SemiEAA'
            elif SwitchDict['Check'] == 'Glc':
                MolColorDF,a = LH.TCACheck(MolColor,MolColor,'TCA')#'TCA'
        except:
            pass
        ColorList = MolColorDF['MolColor']
        ColorList = MolColorDF.loc[Col,'MolColor']
        c = collections.Counter(ColorList);listr1=[];listr2=[];listm1=[];listm2=[];listg1=[];listg2=[];listb1=[];listb2=[];listbl1=[];listpl1=[];listpl2=[];listbl2=[];listgr1=[];listgr2=[]
        listrabs1=[];listrabs2=[];listmabs1=[];listmabs2=[];listgabs1=[];listgabs2=[];listbabs1=[];listbabs2=[];listblabs1=[];listplabs1=[];listplabs2=[];listblabs2=[];listgrabs1=[];listgrabs2=[]
        listposr1=[];listposr2=[];listposm1=[];listposm2=[];listposg1=[];listposg2=[];listposb1=[];listposb2=[];listposbl1=[];listpospl1=[];listpospl2=[];listposbl2=[];listposgr1=[];listposgr2=[]
        if SwitchDict['MolColorScatter']  == 1:
            for i in range(0,len(SubjRmean.columns)):
                if i < (c['red']):
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[(c['red']-1)]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(c['red']-1)]])
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[(c['red']-1)]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(c['red']-1)]])
                    redNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[(c['red']-1)]]);
                    listpos1=list(redNP[redNP>0]);listpos2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(c['red']-1)]])
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i])
    
                    listrabs1 += listabs1; listrabs2 += listabs2;  listr1 += list1; listr2 += list2
                    listposr1 += listpos1; listposr2 += listpos2
                    plt.close()
                    
                    stloc = Col[(c['red'])]
                    magenta = c['red']+c['magenta']
                elif i < magenta:
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[magenta]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(magenta)]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[magenta]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(magenta)]])
                    magentaNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[magenta]]);
                    listpos1=list(magentaNP[magentaNP>0]);listpos2 =list(SubjRstd.loc[Col[i],Col[i]:Col[(magenta)]])
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i])
    
                    listm1 += list1; listm2 += list2; listmabs1 += listabs1; listmabs2 += listabs2;
                    listposm1 += listpos1; listposm2 += listpos2;

                    plt.close()
                    stloc = Col[magenta]
                    green = magenta + c['green']
                elif i < green:
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[green]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[green]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[green]])
                    greenNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[green]])
                    listpos1=list(greenNP[greenNP>0]);
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i]) 
                    listg1 += list1; listg2 += list2;  listgabs1 += listabs1; listgabs2 += listabs2
                    listposg1 += listpos1; listg2 += list2; 
                    plt.close()
                    
                    stloc = Col[green]
                    blue = green + c['blue']

                elif i < (blue):
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[blue]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[blue]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[blue]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[blue]])
                    blueNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[blue]])
                    listpos1=list(blueNP[blueNP>0]);
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i]) 
                    listb1 += list1; listb2 += list2;  listbabs1 += listabs1; listbabs2 += listabs2
                    listposb1 += listpos1; listb2 += list2; 
                    plt.close()                
                    stloc = Col[blue]
                    purple = blue + c['purple']

                elif i < (purple):
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[purple]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[purple]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[purple]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[purple]])
                    purpleNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[purple]])
                    listpos1=list(purpleNP[purpleNP>0]);
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i]) 
                    listpl1 += list1; listpl2 += list2;  listplabs1 += listabs1; listplabs2 += listabs2
                    listpospl1 += listpos1; listpl2 += list2;
                    plt.close()                
                    stloc = Col[purple]
                    black = purple + c['black']-1

                elif i < (black):
                    list1=list(SubjRmean.loc[Col[i],Col[i]:Col[black]]);list2 =list(SubjRstd.loc[Col[i],Col[i]:Col[black]])
                    listabs1=list(abs(SubjRmean.loc[Col[i],Col[i]:Col[black]]));listabs2 =list(SubjRstd.loc[Col[i],Col[i]:Col[black]])
                    blackNP = np.array(SubjRmean.loc[Col[i],Col[i]:Col[black]])
                    listpos1=list(blackNP[blackNP>0]);
                    ax1.scatter(listabs1,listabs2,color=ColorList[i])
                    ax2.scatter(list1,list2,color=ColorList[i])  
                    listbl1 += list1; listbl2 += list2;  listblabs1 += listabs1; listblabs2 += listabs2
                    listposbl1 += listpos1; listbl2 += list2; 
                    plt.close()
                    stloc = Col[black]
                list1=list(SubjRmean.loc[Col[i],stloc:]);list2 =list(SubjRstd.loc[Col[i],stloc:])
                listabs1=list(abs(SubjRmean.loc[Col[i],stloc:]));listabs2 =list(SubjRstd.loc[Col[i],stloc:])
                greyNP = np.array(SubjRmean.loc[Col[i],stloc:])
                listpos1=list(greyNP[greyNP>0]);
                listgr1 += list1;listgr2 += list2;    listgrabs1 += listabs1;listgrabs2 += listabs2
                listposgr1 += listpos1;listgr2 += list2;
                
                ax1.scatter(listabs1,listabs2,color='gray',alpha=0.5)
                ax2.scatter(list1,list2,color='gray',alpha=0.5)  
                plt.close()
                
            Std = list(SubjRstd.values.flatten()[~np.isnan(SubjRstd.values.flatten())])
            absr, absp=pearsonr(AbsAve,Std)
            p='{:e}'.format(absp)
            SaveName='_WMolColor_Abs'
                
            ax1.set_title('R=' + str(round(absr,3)) + ', p=' + p,fontsize=10);ax2.set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp,fontsize=10)
            ax1.set_xlabel('Abs(Ave)',fontsize=10);ax2.set_xlabel('Ave',fontsize=10)              
            ax1.set_ylabel('Std',fontsize=10);ax2.set_ylabel('Std',fontsize=10)#,fontproperties=prop)     fig.tight_layout()
            plot_axis = plt.axis()
            #fig1.savefig(save_dir +'AbsAveStdScatter' + SaveName +'.pdf',format='pdf')
            plt.close()
            #fig2.savefig(save_dir +'AveStdScatter' + SaveName +'.pdf',format='pdf')
            plt.close()       
            
        if NormlSwitch == 0:
            plothistEachMetabo(SubjRmean,Col,c,list(c),save_dir)
        elif NormlSwitch == 1:
            al=1
            hype='bar'
            FileNameTag = 'Nmlz'
            T=True;F=False
            
            ########################################
            plotculmitivehistEachMetabo(SubjRmean,Col,c,list(c),save_dir)
            CList=['red','magenta','green','blue','purple','black','gray']

            msize=3; lwidth=1.5
            ConList1 = [listrabs1, listmabs1, listgabs1, listbabs1,listplabs1, listblabs1, listgrabs1];         ConList2 = [listrabs2, listmabs2, listgabs2, listbabs2, listplabs2, listblabs2, listgrabs2];
            ConList1 = [np.array(ConList1[x]) for x in range(len(ConList1))];         ConList2 = [np.array(ConList2[x]) for x in range(len(ConList2))]
            plt.savefig(save_dir + 'Fig6_c'+'.pdf') 
            plt.close()
            
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
            ##################################3 ##################################3 ##################################3 ##################################3 ##################################3
            fig10, ax9 = plt.subplots(1,1)
            fig10, ax10 = plt.subplots(1,1)
            
            n, bins, patches =ax9.hist(ConList2[0][~np.isnan(ConList2[0])], density=F, cumulative=F,bins = 20,color='red',alpha=al,histtype=hype)        
            yr = np.add.accumulate(n) / n.sum();xr = np.convolve(bins, np.ones(2) / 2, mode="same")[1:]

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

            ax23.axis('off')
            ax9.axis('off'); ax10.axis('off')
            plt.savefig(save_dir + 'Fig6_c'+'.pdf') 
            plt.close()

            
        elif NormlSwitch == 2:  
            FileNameTag = 'Cumulative'
            fig2, ax2 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = abs(SubjRmean[Col[0]:Col[(c['red'])]]).values.flatten(); tempnp4 = abs(SubjRmean[Col[green]:Col[blue]]).values.flatten(); tempnp2 = abs(SubjRmean[Col[(c['red'])]:Col[magenta]]).values.flatten(); tempnp3 = abs(SubjRmean[Col[magenta]:Col[green]]).values.flatten();  tempnp5 = abs(SubjRmean[Col[blue]:Col[black]]).values.flatten()
            ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)])], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')        
            plt.savefig(save_dir + 'AbsHistAve_'+FileNameTag+'.pdf')        


        
            fig2, ax2 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None);
            tempnp1 = SubjRmean[Col[0]:Col[(c['red'])]].values.flatten(); tempnp4 = SubjRmean[Col[green]:Col[blue]].values.flatten(); tempnp2 = SubjRmean[Col[(c['red'])]:Col[magenta]].values.flatten(); tempnp3 = SubjRmean[Col[magenta]:Col[green]].values.flatten();  tempnp5 = SubjRmean[Col[blue]:Col[black]].values.flatten()
            w= list(np.ones(len(tempnp1[~np.isnan(tempnp1)]))/float(len(tempnp1[~np.isnan(tempnp1)])));w += [list(np.ones(len(tempnp2[~np.isnan(tempnp2)]))/float(len(tempnp2[~np.isnan(tempnp2)])))];w+= [list(np.ones(len(tempnp3[~np.isnan(tempnp3)]))/float(len(tempnp3[~np.isnan(tempnp3)])))];w+= [list(np.ones(len(tempnp4[~np.isnan(tempnp4)]))/float(len(tempnp4[~np.isnan(tempnp4)])))];w+= [list(np.ones(len(tempnp5[~np.isnan(tempnp5)]))/float(len(tempnp5[~np.isnan(tempnp5)])))]

            ax2 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)])], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')   #,weights=w     
            plt.savefig(save_dir + 'HistAve_'+FileNameTag+'.pdf')        

            fig3, ax3 = plt.subplots(1,1)
            tempDF = pd.DataFrame(data=None)
            tempnp1 = SubjRstd[Col[0]:Col[(c['red'])]].values.flatten();  tempnp4 = SubjRstd[Col[green]:Col[blue]].values.flatten(); tempnp2 = SubjRstd[Col[(c['red'])]:Col[magenta]].values.flatten();  tempnp3 = SubjRstd[Col[magenta]:Col[green]].values.flatten();   tempnp5 = SubjRstd[Col[blue]:Col[black]].values.flatten()
            ax3 = plt.hist([list(tempnp1[~np.isnan(tempnp1)]), list(tempnp2[~np.isnan(tempnp2)]), list(tempnp3[~np.isnan(tempnp3)]), list(tempnp4[~np.isnan(tempnp4)]), list(tempnp5[~np.isnan(tempnp5)])], stacked=True, bins = 20,color=['red','magenta','green','blue','black'],cumulative=True,histtype='stepfilled')        
            plt.savefig(save_dir + 'HistStd_'+FileNameTag+'.pdf')  
            


        



