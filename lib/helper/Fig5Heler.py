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
from  .StatCal import UseR, calcpeasonr
from .mkHeatmapHelper import draw_heatmapclustering
import collections
from scipy.stats import zscore







class TmPtCorr:
    def __init__(self):
        self.label=[]
        self.subjectName=[]
        self.optiondict=dict()
        self.timepointlist=[]
        self.save_dir=''
        self.MolColorDF=pd.read_excel("./Data/LabelSummary.xlsx",header=0,index_col=0)
    
    def AnalTmPtCorrEachMol(self,DF):
            if self.optiondict['Target'] =='FoldChnage':
                self.AnalFoldChangeEachMol(DF)
            else:
                self.AnalRawEachMol(DF)
                

    
    def AdjustFasting(self,DF):
            ColList0 = [ '0_' +  self.label[jj] for jj in range(len(self.label))]
            ColList10 = ['-10_' +  self.label[jj] for jj in range(len(self.label))]
            AllSubjTmCsZero = DF[ColList0]; AllSubjTmCsminus = DF[ColList10]; 
            AllSubjTmCsZero.columns=self.label;AllSubjTmCsminus.columns=self.label
            FastingDF = pd.DataFrame(data=None,index=self.subjectName,columns=self.label)
            for i in self.label:
                FastingDF[i] = np.nanmean([np.array(AllSubjTmCsminus[i].astype(float)),np.array(AllSubjTmCsZero[i].astype(float))],axis=0)
            return(FastingDF)
            

    def AnalRawEachMol(self,DF): 
            FastingDF=self.AdjustFasting(DF)
            timepointlist=self.timepointlist;
            timepointlist.remove(0);
            print(timepointlist)
            MolCorrDF = pd.DataFrame(data=None,index=self.label,columns=timepointlist)
            MolPvalDF = pd.DataFrame(data=None,index=self.label,columns=timepointlist)
            try:
                MolCorrDF = pd.read_excel(self.save_dir+'Corr_Wpearson.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
                MolPvalDF = pd.read_excel(self.save_dir+'Pvalue_Wpearson.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
            except:
                for jj in self.label:

                    for ii in timepointlist:
                        list1 = list(FastingDF[jj]); list2 = list(DF[str(ii)+'_'+jj])
                        r,p = calcpeasonr(list1,list2)
                        MolCorrDF[ii][jj] = r; MolPvalDF[ii][jj] = p


            MolCorrDF.to_excel(self.save_dir+'/Corr_Wspearman.xlsx');  MolPvalDF.to_excel(self.save_dir+'/Pvalue_Wspearman.xlsx')      
            plt.figure();plt.hist(list(MolCorrDF.values[~pd.isnull(MolCorrDF.values)]),bins=15);plt.gca().tick_params(axis='both',labelsize=20)#plt.xticks(np.arange(-0.2, 0.69, 0.1))
            plt.title('Corr_Wspearman'); plt.savefig(self.save_dir +'Q<0.1_DistOfCorr_Wspearman.pdf');plt.close()
            plt.figure();plt.hist(list(MolPvalDF.values[~pd.isnull(MolPvalDF.values)]),bins=15);plt.gca().tick_params(axis='both',labelsize=20)#plt.xticks(np.arange(-0.2, 0.69, 0.1))
            plt.title('Pvalue_Wspearman'); plt.savefig(self.save_dir +'Q<0.1_DistOfPvalue_Wspearman.pdf');plt.close()
            
            self.optiondict['EngLabel']='RawEachMol'
            UseR(MolPvalDF,self.optiondict)
            self.optiondict['method']='spearman';self.optiondict['Data']=DF
            if len(self.optiondict['SignLabel'])>1:
                MolCorrDF=MolCorrDF.T[self.optiondict['SignLabel']].T;
                self.label=self.optiondict['SignLabel']
            Optiondict={'MolColor':self.MolColorDF,
                    'Annotate' : 1,
                    'Label':self.label,
                    'title':'Correlation vs Fasting',
                    'Threshold':1.5,
                    'numcluster':4,
                    'metric' : 'euclidean'
                    }

            ClstMol,ClstMolDF,Colordict,ClstDF,ClstColorDF,ColorDF,ClstNoColorDF = draw_heatmapclustering(MolCorrDF,1,'MolColor', Optiondict,self.save_dir+'Correlation vs Fasting_')
    def AnalFoldChangeEachMol(self,DF):
            FastingDF=self.AdjustFasting(DF)
            self.timepointlist.remove(-10);self.timepointlist.remove(0)
            for i in range(len(self.timepointlist)):
                self.optiondict['EngLabel'] = self.label;
                self.optiondict['Time'] = str(self.timepointlist[i])
                ColList = [self.optiondict['Time'] + '_' +  self.optiondict['EngLabel'][jj] for jj in range(len(self.optiondict['EngLabel']))]
                tempDF = DF[ColList]
                FastingDF.columns=ColList
                DF[ColList] =tempDF / FastingDF;# DF.columns=ColList  
                
            
            MolCorrDF = pd.DataFrame(data=None,index=self.label,columns=self.timepointlist)
            MolPvalDF = pd.DataFrame(data=None,index=self.label,columns=self.timepointlist)
            try:
                MolCorrDF = pd.read_excel(self.save_dir+'Corr_Wpearson.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
                MolPvalDF = pd.read_excel(self.save_dir+'Pvalue_Wpearson.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
            except:
                for jj in self.label:
                    #MolCorrDF[0][jj] =1 
                    for ii in self.timepointlist:
                        list1 = list(FastingDF[jj]); list2 = list(AllSubjTmCs[str(ii)+'_'+jj])
                        r,p = SC.calcpeasonr(list1,list2)
                        MolCorrDF[ii][jj] = r; MolPvalDF[ii][jj] = p

            MolCorrDF.to_excel(self.save_dir+'/Corr_Wpearson.xlsx');  MolPvalDF.to_excel(self.save_dir+'/Pvalue_Wpearson.xlsx')      
            plt.figure();plt.hist(list(MolCorrDF.values[~pd.isnull(MolCorrDF.values)]),bins=15);#plt.xticks(np.arange(-0.2, 0.69, 0.1))
            plt.title('Corr_Wpearson'); plt.savefig(self.save_dir +'Q<0.1_DistOfCorr_Wpearson.pdf');plt.close()
            plt.figure();plt.hist(list(MolPvalDF.values[~pd.isnull(MolPvalDF.values)]),bins=15);#plt.xticks(np.arange(-0.2, 0.69, 0.1))
            plt.title('Pvalue_Wpearson'); plt.savefig(self.save_dir +'Q<0.1_DistOfPvalue_Wpearson.pdf');plt.close()
            
            self.optiondict['EngLabel']='RawEachMol'
            SC.UseR(MolPvalDF,self.optiondict)
            self.optiondict['method']='pearson';self.optiondict['Data']=DF
            #PlotScatter(self.optiondict,self.save_dir);plt.close()
            if len(self.optiondict['SignLabel'])!=0:MolCorrDF=MolCorrDF.T[self.optiondict['SignLabel']].T;self.label=self.optiondict['SignLabel']
            Optiondict={'Color':LH.MolColor(self.label),
                    'Annotate' : 1,
                    'Label':self.label,
                    'Title':'Correlation vs Fasting'
                    }
            plotTmCs(MolCorrDF.T,Optiondict,self.save_dir)
            ClstMol,Colordict,ClstDF,ClstColorDF,ColorDF = mD.draw_heatmap(MolCorrDF,1,'MolColor',{'title':'Correlation vs Fasting'},self.save_dir+'Correlation vs Fasting_',cmap='bwr')
        



            