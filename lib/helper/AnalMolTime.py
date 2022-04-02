#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 11:16:05 2018

@author: fujita
"""
import numpy as np
import pandas as pd
import pdb; 
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import matplotlib.cm as cm
import warnings 
from scipy.stats import zscore
from scipy import stats



###
def CalcAUC(list1,Timelist):
    AUC=0
    if len(list1) == len(Timelist):
        length = len(list1)
        for i in range(0,length-1):
            base = abs(Timelist[i+1] - Timelist[i])
            height = (list1[i] + list1[i+1]) / 2.
            
            if height != height:
                height = 0
            AUC += base * height
    else:
        pdb.set_trace()
    return(AUC)
### 
def CalcAUChalf(list1,Timelist,ResponseIdxDict,ij):
    AUCList=[]
    if len(list1) == len(Timelist):
        length = len(list1)
        for i in range(0,length-1):
            base = abs(Timelist[i+1] - Timelist[i])
            height = (list1[i] + list1[i+1]) / 2.
            if height != height:
                height = 0
            if i == 0:
                AUCList.append(base * height)
            else:
                AUCList.append( AUCList[i-1] + base * height)
    else:
        pdb.set_trace()
    
    if len(AUCList)==0:
        print(list1)

    if ij in ResponseIdxDict['DownBoth']:

        AUCaccum = AUCList[-1]
        AUCstar = AUCaccum / 2
        tl = np.array(AUCList) <= AUCstar
    else:
        AUCaccum = AUCList[-1]
        AUCstar = AUCaccum / 2
        tl = np.array(AUCList) <= AUCstar
    TCIdx = [];TCIdxPre = []
    count=0
    try:
        for jj in range(len(tl)): 
            if (tl[jj] == 0) and (count==0):
                TCIdxPre=Timelist[jj+1]   
                TCIdxbefore = Timelist[jj]
                val1 = list1[jj+1] 
                val2 = list1[jj]    
                AUCPre =AUCList[jj-1] 
                count+=1
    
        x=[TCIdxbefore,TCIdxPre]
        y=[val2,val1]
        slope, intercept, r_value, _, _ = stats.linregress(x, y)
        a = slope
        b = (intercept + val2 - slope*TCIdxbefore)
        c = -(val2*TCIdxbefore + TCIdxbefore*intercept + 2*AUCstar - 2*AUCPre)
             
        if slope==0:
            pass
        try:
            x1, x2 = solv_quadratic_equation(float(a), float(b), float(c))
            if (x1>=0) and (x2>=0) :
                if (x1>=TCIdxbefore) and ( x1<=TCIdxPre):
                    tstar = x1
                else:
                    tstar=x2
                
            elif x1 >= 0:
                tstar = x1
            elif x2 >= 0:
                tstar = x2

        except:
            TCIdx+=[240]
    except:  
                AUCaccum=np.nan
                tstar=np.nan      
    return(AUCaccum,tstar)
    
def solv_quadratic_equation(a, b, c):
    if a==0:
        x_1= -c / b
        x_2 = -c/b
    else:
        D = (b**2 - 4*a*c) ** (1/2)
        x_1 = (-b + D) / (2 * a)
        x_2 = (-b - D) / (2 * a)

    return x_1, x_2


###
def CorrMolEachSubjHelper(DF,MolDict,Label):
    RDF = pd.DataFrame(data=None, index=Label, columns=Label)
    PDF = pd.DataFrame(data=None, index=Label, columns=Label)
    dellist = lambda items, indexes: [item for index, item in enumerate(items) if index not in indexes]


    st = 0
    for ij in range(len(MolDict)):
        list1 = list(DF[Label[ij]])
        for j in range(len(MolDict)):
            list2 = list(DF[Label[j]])
            if np.isnan(list2).any() or np.isnan(list1).any(): 
                       print(str(Label)+'nan')
                       aaa = list(np.where(np.isnan(list2))[0])
                       bbb = list(np.where(np.isnan(list1))[0])
                       aaa = list(set(aaa+bbb))
                       list1temp=dellist(list1,aaa)
                       list2temp = dellist(list2,aaa)
                       r, p = pearsonr(list1temp,list2temp)
                       RDF.loc[Label[ij],Label[j]] = r
                       PDF.loc[Label[ij],Label[j]] = p
                       
            else:           
                r, p = pearsonr(list1,list2)
                RDF.loc[Label[ij],Label[j]] = r
                PDF.loc[Label[ij],Label[j]] = p
    return(RDF,PDF)  
     


    
    
