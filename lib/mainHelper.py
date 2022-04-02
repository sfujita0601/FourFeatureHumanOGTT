#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 22:11:02 2021

@author: fujita
"""


import numpy as np
import pandas as pd
import scipy.io
import itertools
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import scipy.stats as sts
from .helper.StatCal import UseR, calcpeasonr
from .helper  import VolcanoPlotHelper as VP
from .helper  import mkHeatmapHelper as mHH
from .helper  import PCAPre as PCAPre
from .helper  import GraphHelper as GH
from .helper  import Fig5Heler as ASTH
from .helper import Fig6Helper as pSFP
from .helper import AnalMolTime as AmT


class OGTThuman(object):
    
    
    def __init__(self,TwodData,save_dir):
        self.label=[]
        self.optiondict=dict()
        self.timepointlist=[0,10,20,30,45,60,75,90,120,150,180,210,240]
        self.save_dir=save_dir 
        self.TwodData = TwodData
        self.FastingData = self.mkFasting()
        self.ThreedData = self.mk3dData()
        self.MolLabel = list(self.FastingData.columns)
        self.TwodSubjTimeMolData= self.mk2dSubjTimeMolData()
        self.TwodConcatTimeSubjMolData= self.TwodConcatTimeSubjMolData()
        self.YDF = self.mkNormalizeByStd(self.ThreedData,self.save_dir)
        self.MolColorDF=pd.read_excel("./Data/LabelSummary.xlsx",header=0,index_col=0)
        self.SignLabel = ['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
                               'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
                               'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate','Growth hormone']
        
    def mkFasting(self):
        MolLabel=list(self.TwodData.index)
        FastingDF=pd.DataFrame(data=None,index=range(20),columns=MolLabel)
        for j in range(1,21):
                FastingDF.iloc[j-1,:]=list(np.nanmean(self.TwodData.iloc[:,[14*(j-1),14*(j-1)+1] ],axis=1))    
        return(FastingDF)
    def mk3dData(self):
            for i in range(1,21):
                AddDF = self.TwodData.iloc[:,(1+14*(i-1)):(14+14*(i-1))].T     
                AddDF.iloc[0,:]=list(self.FastingData.iloc[i-1,:]);AddDF.index=self.timepointlist.copy()
                if i== 1: 
                    TmCsDF = self.TwodData.iloc[:,(1+14*(i-1)):(14+14*(i-1))].T
                    
                    TmCsDF.iloc[0,:]=list(self.FastingData.iloc[i-1,:]);TmCsDF.index=self.timepointlist.copy()
                else:
                    TmCsDF=pd.concat([TmCsDF,AddDF],axis=0) 
            SubjTmCs3dData = pd.DataFrame(data=np.array(TmCsDF),index=pd.MultiIndex.from_product([range(20), list(AddDF.index)]),columns=list(TmCsDF.columns))
            SubjTmCs3dData=SubjTmCs3dData.rename_axis(['Subject','Time'], axis=0)
            return(SubjTmCs3dData)
    def mk2dSubjTimeMolData(self):
            timepointlist=self.timepointlist.copy()
            timepointlist.append(-10)
            for i in timepointlist:
                if i==0:#
                    tempSubjTimeMol =    self.ThreedData.swaplevel(0).loc[i] 
                    tempSubjTimeMol.columns= [str(i) +'_'+ jj for jj in list(tempSubjTimeMol.columns)]
                elif i==-10:
                    SubjTimeMol =    self.ThreedData.swaplevel(0).loc[0] 
                    SubjTimeMol.columns= [str(i) +'_'+ jj for jj in list(SubjTimeMol.columns)]
                    tempSubjTimeMol=pd.concat([tempSubjTimeMol,SubjTimeMol],axis=1) 
                
                else:
                    SubjTimeMol =    self.ThreedData.swaplevel(0).loc[i]
                    SubjTimeMol.columns= [str(i) +'_'+ jj for jj in list(SubjTimeMol.columns)]
                    tempSubjTimeMol=pd.concat([tempSubjTimeMol,SubjTimeMol],axis=1) 
                    #tempDelta = SubjTmCsDelta / ThreedData.fillna(0).std(level=1,ddof=0).mean() 

            return(tempSubjTimeMol)   
    
    def TwodConcatTimeSubjMolData(self):
            timepointlist=self.timepointlist.copy()
            for i in range(20):
                if i==0:#
                    tempSubjTimeMol =    self.ThreedData.loc[i] - self.FastingData.iloc[i,:]
                else:
                    SubjTimeMol =    self.ThreedData.loc[i]- self.FastingData.iloc[i,:]
                    tempSubjTimeMol=pd.concat([tempSubjTimeMol,SubjTimeMol],axis=0)         
            return(tempSubjTimeMol)   


        
    def AnalFig4(self):
        save_dir=self.save_dir+'Fig4/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir) 
        TPSI=self.AnalTPSI(save_dir)
        TVRI=self.AnalTVRI(save_dir)
        MolColorDF_Eng = self.MolColorDF
        Optiondict={'calcR':'','xlabel':'TVRI','ylabel':'TPSI','Annotate':0,'Label':list(MolColorDF_Eng.index),'title':'TVRI vs TPSI','x_err':'','y_err':''}             

        GH.mkScatterWHist(list(TVRI),list(TPSI),save_dir,list(MolColorDF_Eng['MolColor']),Optiondict)#
    def mkNormalizeByStd(self,ThreedData,save_dir):
        for i in range(20):
            if i==0:
                SubjDelta =    ThreedData.loc[i] - ThreedData.loc[i].iloc[0,:] 
                tempDelta = SubjDelta / ThreedData.fillna(0).std(level=1,ddof=0).mean() 
            else:
                SubjTmCsDelta = ThreedData.loc[i] - ThreedData.loc[i].iloc[0,:]  
                SubjDelta=pd.concat([SubjDelta,SubjTmCsDelta],axis=0) 
                tempDelta = SubjTmCsDelta / ThreedData.fillna(0).std(level=1,ddof=0).mean() 
            ##tempDelta.to_excel(save_dir + 'SubjTmCs_'+SubjectName[i]+'DivBymeanRawStdDelta.xlsx')
        ThreedData_rev = pd.DataFrame(data=np.array(SubjDelta),index=pd.MultiIndex.from_product([range(20), list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))

        Subjmean = ThreedData_rev.fillna(0).mean(level=1)/ ThreedData.fillna(0).std(level=1,ddof=0).mean() 
      
        return(Subjmean)
    def AnalFig5(self):
        save_dir=self.save_dir+'Fig5/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir) 
        EngLabel=['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
                   'Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
                   'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate',
               'Growth hormone','Glu + threo-beta-methylaspartate']#'Glu + threo-beta-methylaspartate','Citrate',
         
        TmPtCorr_class = ASTH.TmPtCorr(); 
        SwitchDict=dict({'Target':'Raw'});   timepointlist=self.timepointlist.copy(); #timepointlist.append(-10)
        SwitchDict['SignLabel']=EngLabel
        TmPtCorr_class.save_dir=save_dir;TmPtCorr_class.label=self.MolLabel;TmPtCorr_class.optiondict=SwitchDict;TmPtCorr_class.timepointlist=timepointlist
        TmPtCorr_class.AnalTmPtCorrEachMol(self.TwodSubjTimeMolData)

    def AnalFig6(self):
        save_dir=self.save_dir+'Fig6/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir) 
        Optiondict=dict({'relation':'pearson'})
        Label=list(self.TwodConcatTimeSubjMolData.columns)
        if Optiondict['relation']=='pearson':
            RDF,PDF =  AmT.CorrMolEachSubjHelper(self.TwodConcatTimeSubjMolData,Label,Label)  
        ##self.TwodConcatTimeSubjMolData.to_excel(save_dir+'tempData.xlsx') 
        ##RDF.to_excel(save_dir + '/RDF.xlsx');
        RDFrev = RDF.copy();RDFrev =self.adjustMatrlower(RDFrev);##RDFrev.to_excel(save_dir + '/RDFrev.xlsx');
        aa=RDFrev.astype(float).values.flatten()
        aab=aa[~np.isnan(aa)]
        
        ##PDF.to_excel(save_dir + '/PDF.xlsx');
        PDFrev = PDF.copy(); 
        UseR(PDFrev,{'EngLabel':'TPSMqval'})        
        mHH.draw_heatmap(RDF,'both','MolColor',{'title':''},save_dir,cmap='PuOr_r')

        SwitchDict = {}
        SwitchDict['MolColorScatter']  = 1
        SwitchDict['Check']  = 'Amino'
        SwitchDict['MolColor']=self.MolColorDF

        SubjRmeanrev = self.adjustMatrlower(RDF); SubjRstdrev = self.adjustMatrlower(RDF)

        SubjRmeanrev.to_excel(save_dir+'/SubjRrev.xlsx'); 
        for i in range(2):
                    NormlSwitch = i
                    pSFP.PlotScatter(save_dir,SubjRmeanrev.astype(float),SubjRstdrev.astype(float),3,NormlSwitch,SwitchDict)      

    def AnalFig7(self):
        from .helper import mkNetWorkGraphHelper as GraphHel
        save_dir=self.save_dir+'Fig7SFig6_8_9/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir) 
        LabelSum =  self.MolColorDF
        Optiondict=dict({'relation':'pearson'})
        Label=list(self.TwodConcatTimeSubjMolData.columns)
        if Optiondict['relation']=='pearson':
            RDF,PDF =  AmT.CorrMolEachSubjHelper(self.TwodConcatTimeSubjMolData,Label,Label)  
        RDFrev = self.adjustMatrlower(RDF); 
        ColorSwitchDict  = {1:'ClstColor',2:'MolColor',3:'TimeVarColor'}#
        OptionDict = {'Edge':'CorrCoef', 
                      'EdgeColor':'Color_posneg'}  
        OptionDict['mkglobalDF_Time']='' 
        OptionDict['mkEdge'] = 'Thresh'
        OptionDict['mkEdge_double'] = ''
        OptionDict['pos']='positive' 
        OptionDict['pos_Global']='LabelLoop'
        OptionDict['Node_size'] = ''
        OptionDict['Check'] = ''
        OptionDict['AminoCheck'] = 'EAA'
        OptionDict['MolColor'] =self.MolColorDF
    
        AbsSwitch = 1
        if AbsSwitch == 1:
            RDFrev = np.abs(RDFrev)
        fiveIdx = 0.8
        for i in range(2,3):
            ClstList=[];CorrList=[];VarnumEachClstnode=[]; AvenumEachClstnode =[] ;templist=[]; NodeList =[]
            for Thresh in [0.100+0.05*x for x in range( 1 ) ]:#18
                print('Threshold: %02.2f' %Thresh)
       
                templist = GraphHel.mkNetWorkGraph(RDFrev,LabelSum,Thresh,ColorSwitchDict[i],OptionDict,[],save_dir+'_'+str(Thresh)+'/')
                ClstList += [templist[0]]; AvenumEachClstnode += [templist[1]] ; VarnumEachClstnode += [templist[2]]; NodeList += [templist[3]]
                CorrList += [Thresh]   
            fig1, ax1 = plt.subplots()

            ln1 = ax1.plot(CorrList,ClstList,label='# of Components',marker = 'o',lw=1.5);ax11 = ax1.twinx() ;
            ln11 = ax11.plot(CorrList,NodeList,label='# of Nodes',c='orange',marker = 'o',lw=1.5); 
            params = {'mathtext.default': 'regular' }   
            plt.rcParams.update(params) 
            ax1.set_xlabel('Average of $TPSM^{Abs}$',size=20); 
            ax11.set_ylabel('# of Nodes',size=20);ax1.set_ylabel('# of Components ',size=20); 
    
    
            h1, l1 = ax1.get_legend_handles_labels();
            h11, l11 = ax11.get_legend_handles_labels(); 
            ax1.legend(h1+h11, l1+l11, loc='center left'); plt.tick_params(labelsize=20) ;ax11.xaxis.set_tick_params(labelsize=20)
            xmin, xmax, ymin, ymax = ax11.axis()
    
            ax1.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
            ax11.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    
            ax11.tick_params(axis='both',labelsize=20);  ax1.tick_params(axis='both',labelsize=20)

            fig1.savefig(save_dir + 'Coef vs Clst or Node_' + ColorSwitchDict[i] + '.pdf',bbox_inches="tight")
    
            fig2, ax2 = plt.subplots()
            ln2 = ax2.plot(CorrList,AvenumEachClstnode,label='Mean',marker = 'o',lw=1.5); 
            ax3 = ax2.twinx() ;
            ln3 = ax3.plot(CorrList,VarnumEachClstnode,label='Variance',c='orange',marker = 'o',lw=1.5);
    
            ax2.set_xlabel('Average of $TPSM^{Abs}$',size=20); ax2.set_ylabel('Mean of Component size',size=20);ax3.set_ylabel('Variance of Component size',size=20); 
            h2, l2 = ax2.get_legend_handles_labels();h3, l3 = ax3.get_legend_handles_labels(); 
            ax2.legend(h2+h3, l2+l3, loc='center'); plt.tick_params(labelsize=20) ; ax3.xaxis.set_tick_params(labelsize=20)
            xmin, xmax, ymin, ymax = ax3.axis()
            ax2.legend( loc='center left', borderaxespad=0, fontsize=18)
            ax2.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
            ax3.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    
            ax2.legend_ = None
            ax2.set_yticks([0, 20, 40, 60, 80])
            ax3.set_yticks([0, 20, 40, 60, 80])
    
            ax2.tick_params(axis='both',labelsize=20);  
    
            ax3.set_yticks([0, 500,1000,1500,2000,2500,3000])
            ax3.tick_params(axis='both',labelsize=20)
            fig2.savefig(save_dir + 'Coef vs Ave or Var size of Clst_' + ColorSwitchDict[i] + '.pdf',bbox_inches="tight")
            
    def AnalFig8FigS9(self):
        save_dir=self.save_dir+'Fig8FigS9/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)         
        def plot3D(DF,OptionDict,save_dir):
            try:
                MolColorDF= OptionDict['MolColorDF']
            except:
                MolColorDF = self.MolColorDF
           
        
            ColorList = list(MolColorDF['MolColor'][list(DF.index)])
            ClstColorList = list(MolColorDF['ClstColor'])

            ColLabel = list(DF.columns)
            MolAllLabel = list(DF.index)
            elevList=[30,30,45,45,45,0,0,0]
            azimList=[148,315,45,135,225,45,135,225]
            count=-1
            for i,j,k in  list(itertools.permutations([0,1,2], 3)):
                    count+=1
                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')#ax = Axes3D(fig,projection='3d')
        
                    XLabel=ColLabel[1]
                    YLabel=ColLabel[2]
                    ZLabel=ColLabel[0]
                    X = DF[XLabel]
                    Y = DF[YLabel]
                    Z = DF[ZLabel]
                    ax.scatter(X, Y, Z,'o',c=ColorList)
                    
                    AnDict=dict({'Glucose':'Glc','Insulin':'Ins','C-peptide':'CRP','GIP(Active)':'GIP','Pyruvate':'Pyr','Total bile acid':'TBA',
                           'Citrate':'Cit','Cortisol':'Cor','Free fatty acid':'FFA','Total ketone body':'Ketone','Glutamic acid':'Glu',
                           'Citrulline':'Citr','Methionine':'Met','Isoleucine':'Ile','Leucine':'Leu','Tyrosine':'Tyr','4-Methyl-2-oxopentanoate':'4M2O','Glu + threo-beta-methylaspartate':'Glu+TBM','Growth hormone':'GH'})
        
                    ax.set_xlabel(XLabel,size=10)
                    ax.set_ylabel(YLabel,size=10,rotation='vertical')
                    ax.set_zlabel(ZLabel,size=10,rotation=270)
                    
                    ax.set_xlim(np.min(X), np.max(X))
                    ax.set_ylim(np.min(Y), np.max(Y))
                    ax.set_zlim(np.min(Z), np.max(Z))
                    
                    ax.tick_params(labelsize=10)
                
                    save_dir_figs = save_dir + '/figs'
                    if not os.path.isdir(save_dir_figs):
                        os.makedirs(save_dir_figs)
                    ax.view_init(elev=elevList[count], azim=azimList[count])
                    plt.savefig(save_dir+'3D'+'_'+str(i)+str(j)+str(k)+'.pdf')
            
            pngflag=0
            if pngflag == 1:
                for angle in range(0, 360):
                    ax.view_init(30, angle)
                    plt.savefig(save_dir_figs +'/{0}_{1:03d}.jpg'.format('3D', angle))           

        def calcAUC_AUChalf(self):
            timepointlist=self.timepointlist.copy()
            AUCDF = pd.DataFrame(data=None, index = range(20),columns=list(self.ThreedData.columns))
            TAUChalfDF = pd.DataFrame(data=None, index = range(20),columns=list(self.ThreedData.columns))
            OptionDict=dict()
            OptionDict['DownBoth']=''
            for jj in range(20):      
                for ii in range(len(self.ThreedData.columns)):       
                    list1 = list(self.ThreedData.loc[jj].iloc[:,ii] - self.FastingData.loc[jj].iloc[ii]); 
                    AUC= AmT.CalcAUC(list1,timepointlist)
                    AUCDF.iloc[jj,ii]=AUC
                    list1 = list(np.abs(self.ThreedData.loc[jj].iloc[:,ii] - self.FastingData.loc[jj].iloc[ii]));                    
                    AUC, tstar = AmT.CalcAUChalf(list1,timepointlist,OptionDict,'')
                    TAUChalfDF.iloc[jj,ii] = tstar
            return(AUCDF,TAUChalfDF)
            
        def calcCorrAUC_TAUChalfEach(self,AUCDF,TAUChalfDF,save_dir): 
            def calculate_CorrWpearson(df,OptionDict):  
                OtherList=self.Label;
                if OptionDict['Target'] in OtherList:
                    OtherList.remove(OptionDict['Target'])
                pval_mat = np.zeros((1, len(OtherList)))
                corr_mat = np.zeros((1, len(OtherList)))
                counti=-1;countj=-1
                for i in [OptionDict['Target']]:
                    counti+=1
                    for j in OtherList:   
                        countj+=1
                        corrtest = calcpeasonr(list(df[i]), list(df[j]))
                        corr_mat[counti,countj] = corrtest[0] 
                        pval_mat[counti,countj] = corrtest[1] 
                df_corr = pd.DataFrame(corr_mat, columns=OtherList, index=[0])
                df_pval = pd.DataFrame(pval_mat, columns=OtherList, index=[0])
                return(df_corr,df_pval)
            OptionDict=dict();OptionDict['Target']='Glucose'
            CorrDF,PvalDF = calculate_CorrWpearson(AUCDF,OptionDict)
            CorrDF.to_excel(save_dir+'CorrelationCoefficient_AUC_Glucose.xlsx');PvalDF.to_excel(save_dir+'Pvales_AUC_Glucose.xlsx');
            UseR(PvalDF,{'EngLabel':'GlucoseAUC'})

            CorrDF,PvalDF = calculate_CorrWpearson(TAUChalfDF,OptionDict)
            CorrDF.to_excel(save_dir+'CorrelationCoefficient_TAUChalf_Glucose.xlsx');PvalDF.to_excel(save_dir+'Pvales_TAUChalf_Glucose.xlsx');
            UseR(PvalDF,{'EngLabel':'GlucoseTAUChalf'})            

        def calcAUC_AUChalf_YDF(self,YDF):
            timepointlist=self.timepointlist.copy()
            PropertyDF = pd.DataFrame(data=None, index = list(YDF.columns),columns=['AUC','TAUChalf'])
            TAUChalfDF = pd.DataFrame(data=None, index =list(YDF.columns),columns=['AUC','TAUChalf'])
            OptionDict=dict()
            OptionDict['DownBoth']=''
   
            for jj in list(YDF.columns):     
                    list1 = list(YDF[jj]); 
                    AUC= AmT.CalcAUC(list1,timepointlist)
                    PropertyDF.loc[jj,'AUC']=AUC
                    list1 = list(np.abs(YDF[jj]));                    
                    AUC, tstar = AmT.CalcAUChalf(list1,timepointlist,OptionDict,'')
                    PropertyDF.loc[jj,'TAUChalf'] = tstar
            return(PropertyDF)
        
        def calcCVAUC_TAUChalf(self,AUCDF,TAUChalfDF,PropertyDF,save_dir):
            try:
                Score =pd.read_excel(self.save_dir+'Fig2_3FigS4_5/Score.xlsx',header=0,index_col=0)
            except:
                print('Please Analyse Fig2_3 first')
                
            MolColorDF_Eng = self.MolColorDF
            Optiondict={'calcR':'','xlabel':'CV of AUC','ylabel':'CV of TAUChalf','Annotate':1,'Label':list(MolColorDF_Eng.index),'title':'AUCvsAUChalf','x_err':'','y_err':''}             

            OtherList=self.SignLabel
            OptionDict=dict();OptionDict['DelMolList']=['Total bile acid','Pyruvate',
                               'Cortisol','Glutamic acid','Citrate',
                               'Glu + threo-beta-methylaspartate','Growth hormone']

            if OptionDict['DelMolList'][0] in OtherList:
                [OtherList.remove(OptionDict['DelMolList'][i]) for i in range(len(OptionDict['DelMolList']))  ]  
            AUCCV = AUCDF.abs()[OtherList].std(ddof=1) / AUCDF.abs()[OtherList].mean()
            TAUhalfCCV = TAUChalfDF.abs()[OtherList].std(ddof=1) / TAUChalfDF.abs()[OtherList].mean()
            Optiondict['Label']=OtherList
            GH.mkScatterWHist(list(AUCCV ),list(TAUhalfCCV),save_dir,list(MolColorDF_Eng['MolColor'].loc[OtherList]),Optiondict)#2つのリストの散布図+ヒストグラム
            Score.index=list(MolColorDF_Eng.index) 
            PC1 = Score.loc[OtherList][0]; PC2 = Score.loc[OtherList][1]; 
            Optiondict={'calcR':'','xlabel':'AUC','ylabel':'Score of PC1','Annotate':1,'Label':list(MolColorDF_Eng.index),'title':'AUCvsPC1','x_err':'','y_err':''}             
            Optiondict['Label']=OtherList;Optiondict['calcR']='pearson'

            GH.mkScatterWHist(list(PropertyDF.loc[OtherList]['AUC']),list(PC1),save_dir,list(MolColorDF_Eng['MolColor'].loc[OtherList]),Optiondict)
            Optiondict={'calcR':'','xlabel':'TAUChalf','ylabel':'Score of PC2/ Score of PC1','Annotate':1,'Label':list(MolColorDF_Eng.index),'title':'TAUChalfvsPC1Divby2','x_err':'','y_err':''}             
            Optiondict['Label']=OtherList;Optiondict['calcR']='pearson'

            GH.mkScatterWHist(list(PropertyDF.loc[OtherList]['TAUChalf']),list(PC2/PC1),save_dir,list(MolColorDF_Eng['MolColor'].loc[OtherList]),Optiondict)

            
        DegDF = pd.read_excel('./Result/Fig7SFig6_8_9/_0.6/DegNormalized.xlsx',header=0,index_col=0)#,encoding = "ISO-8859-1",index_col=0)
        TPSI=self.AnalTPSI(save_dir)
        TVRI=self.AnalTVRI(save_dir) 
         
        VarDF = pd.DataFrame(data=None,index=list(DegDF.index))
        VarDF['TPSI'] = TPSI
        VarDF['TVRI'] = list(TVRI)

        VarDF['numEdge']=list(DegDF['numEdge'])
        VarDF.to_excel(save_dir+'3Dplot.xlsx')
       
        plot3D(VarDF,dict(),save_dir)
        
        AUCDF,TAUChalfDF = calcAUC_AUChalf(self)
        self.Label=self.SignLabel.copy()
        calcCorrAUC_TAUChalfEach(self,AUCDF,TAUChalfDF,save_dir)
        PropertyDF=calcAUC_AUChalf_YDF(self,self.YDF)
        calcCVAUC_TAUChalf(self,AUCDF,TAUChalfDF,PropertyDF,save_dir)

    def AnalFigS2(self,file_dir ):
        save_dir=self.save_dir+'FigS2/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)  
        def PreprocessOWTT(file_dir):
            OWTTData=pd.read_excel(file_dir+'Metabolites_and_hormones_data.xlsx',header=0,index_col=0,sheet_name='Raw (Water validation)')    
            DelList=['Unit']
            OWTTData=OWTTData.drop(DelList,axis=1)
            OWTTData=OWTTData.loc[['Glucose','Insulin']]
            Data=pd.DataFrame(data=None,index=[0,10,20,30,60,90,120],columns=['Glucose_Mean','Glucose_Std','Insulin_Mean','Insulin_Std'])
            for i in range(len(Data.index)):
                Data.iloc[i,0] = OWTTData.iloc[0,(i*3):((i*3+3))].mean()
                Data.iloc[i,1] = OWTTData.iloc[0,(i*3):((i*3+3))].std(ddof=1)/np.sqrt(len(OWTTData.iloc[0,(i*3):((i*3+3))]))
                Data.iloc[i,2] = OWTTData.iloc[1,(i*3):((i*3+3))].mean()
                Data.iloc[i,3] = OWTTData.iloc[1,(i*3):((i*3+3))].std(ddof=1)/np.sqrt(len(OWTTData.iloc[1,(i*3):((i*3+3))])) 
        
            return(Data)
        def PreprocessOGTT(file_dir):
            OGTTData=pd.read_excel(file_dir+'Metabolites_and_hormones_data.xlsx',header=0,index_col=0,sheet_name='Raw (glucose validation)')    
            DelList=['Unit']
            OGTTData=OGTTData.drop(DelList,axis=1)
            OGTTData=OGTTData.loc[['Glucose','Insulin']]
            Data=pd.DataFrame(data=None,index=[0,10,20,30,60,90,120],columns=['Glucose_Mean','Glucose_Std','Insulin_Mean','Insulin_Std'])
            for i in range(len(Data.index)):
                Data.iloc[i,0] = OGTTData.iloc[0,(i*3):((i*3+3))].mean()
                Data.iloc[i,1] = OGTTData.iloc[0,(i*3):((i*3+3))].std(ddof=1) /np.sqrt(len(OGTTData.iloc[0,(i*3):((i*3+3))]))
                Data.iloc[i,2] = OGTTData.iloc[1,(i*3):((i*3+3))].mean()
                Data.iloc[i,3] = OGTTData.iloc[1,(i*3):((i*3+3))].std(ddof=1)  /np.sqrt(len(OGTTData.iloc[1,(i*3):((i*3+3))]))     
            return(Data)
        OWTTData = PreprocessOWTT(file_dir)
        OGTTData = PreprocessOGTT(file_dir)

        plt.figure(figsize=(12,10))
        plt.errorbar([0,10,20,30,60,90,120],list(OWTTData.iloc[:,0]),yerr =list(OWTTData.iloc[:,1]),c='cyan', fmt='o',ecolor='cyan',linestyle = 'solid',linewidth=2.0,markersize=5.0)
        plt.errorbar([0,10,20,30,60,90,120],OGTTData.iloc[:,0],OGTTData.iloc[:,1],c='red', fmt='o',ecolor='red',linestyle = 'solid',linewidth=2.0,markersize=5.0)
        plt.tick_params(labelsize=20)
        plt.ylabel('mg/dL',size=20);plt.xlabel('Time (min)',size=20)
        plt.savefig(save_dir + '/TimeCourse_Glucose'  +'.pdf',bbox_inches="tight")
        ##plt.savefig(save_dir + '/TimeCourse_Glucose'  +'.png',bbox_inches="tight")
        plt.close()

        plt.figure(figsize=(12,10))
        plt.errorbar([0,10,20,30,60,90,120],list(OWTTData.iloc[:,2]),yerr =list(OWTTData.iloc[:,3]),c='cyan', fmt='o',ecolor='cyan',linestyle = 'solid',linewidth=2.0,markersize=5.0)
        plt.errorbar([0,10,20,30,60,90,120],OGTTData.iloc[:,2],OGTTData.iloc[:,3],c='red', fmt='o',ecolor='red',linestyle = 'solid',linewidth=2.0,markersize=5.0)
        plt.tick_params(labelsize=20)
        plt.ylabel('μU/mL',size=20);plt.xlabel('Time (min)',size=20)
        
        plt.savefig(save_dir + '/TimeCourse_Insulin'  +'.pdf',bbox_inches="tight")
        ##plt.savefig(save_dir + '/TimeCourse_Insulin'  +'.png',bbox_inches="tight")
        plt.close()   
        
        
    def AnalAddFig1(self,file_dir):
        save_dir=self.save_dir+'AddFig1/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)  
        def Preprocess(file_dir):
            ErrorData=pd.read_excel('./Data/MeasurementError.xlsx',header=0,index_col=0)    
            MolColor =   self.MolColorDF
            FastingCV = self.FastingData.std(ddof=1,axis=0)*100 /self.FastingData.mean(axis=0)
            DF=pd.concat([ErrorData,FastingCV,MolColor],join='inner',axis=1)
            DF=DF.rename(columns={0: 'Between different sample CV(%)'})            
            DF=DF.rename(columns={'Within-run reproducibility(%)': 'Within same sample CV(%)'})
            DF=DF.drop(['Citrate'],axis=0)
            return(DF)
        Data=Preprocess(file_dir)
        fig,ax = plt.subplots(figsize=(6,6))
        list1=list(Data['Within same sample CV(%)'])
        list2=list(Data['Between different sample CV(%)'])
        Data.to_excel(save_dir+'MeasurementErrorvsCVbetweenIndividual.xlsx')
        ax.scatter(Data['Within same sample CV(%)'],Data['Between different sample CV(%)'],c=Data['MolColor'],s=20)
        for i in range(len(Data.index)):
            ax.annotate(list(Data.index)[i],fontsize=10 , xy=(list1[i],list2[i]))

        xmin, xmax, ymin, ymax = ax.axis() 
        x = np.arange(xmin, xmax, 0.1)
        ax.plot(x,x,c='black',linewidth=1)
        ax.tick_params(axis="both", labelsize=20)
        plt.ylabel('Between different sample CV',size=30);plt.xlabel('Within same sample CV',size=30)
        
        plt.savefig(save_dir + '/MeasurementErrorvsFastingCV'  +'.pdf',bbox_inches="tight")
        ##plt.savefig(save_dir + '/MeasurementErrorvsFastingCV'  +'.png',bbox_inches="tight")
        plt.close()          





    def AnalTPSI(self,save_dir):
            def calcEachMolInterSubjCorr(self):  
                    for i in range(20):
                        if i==0:#
                            SubjDelta =    self.ThreedData.loc[i] - self.ThreedData.loc[i].iloc[0,:] 
                        else:
                            SubjTmCsDelta = self.ThreedData.loc[i] - self.ThreedData.loc[i].iloc[0,:]  
                            SubjDelta=pd.concat([SubjDelta,SubjTmCsDelta],axis=0) 
                    ThreedData_rev = pd.DataFrame(data=np.array(SubjDelta),index=pd.MultiIndex.from_product([range(20), list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))
 
                    NewLabel =list(ThreedData_rev.columns) 
                    MolColorDF_Eng = self.MolColorDF

                    MolDict = {}; EachMolVarCorrDF = pd.DataFrame(data=None,index=['R','p'],columns = NewLabel)
                    for ii in NewLabel:
                        SubjTmCs1=[]; SubjTmCs2=[]
                        for jj in range(20):#
                            for kk in range(20 - (jj+1)):
                                SubjTmCs1 += list(ThreedData_rev.loc[jj][ii])
                                SubjTmCs2 += list(ThreedData_rev.loc[kk+(jj+1)][ii])
                        R,P=calcpeasonr(SubjTmCs1,SubjTmCs2)#
                        MolDict.update({ii:[R,P]})
                        EachMolVarCorrDF[ ii ]['R'] = R
                        EachMolVarCorrDF[ ii ]['p']= P
                    MolAveList = [ list(MolDict.values())[j][0] for j in range(len( MolDict) ) ] ; 
                    MolLabel = list(MolDict.keys())
                    x = np.linspace(0, len(MolAveList),len(MolAveList) )

                    MolColorDF_Eng['AveCorr'] = MolAveList 
                    MolColorDF_Eng = MolColorDF_Eng.sort_values(by='AveCorr') ; MolLabel=list(MolColorDF_Eng.index); #DFColor = pd.concat([DFNew,ColorDF], axis=1, join_axes=[DFNew.index])
                    TPSI = list(MolColorDF_Eng['AveCorr']); filename = '_Sorted'
                        
                    fig = plt.figure(figsize=(4,9.6))
                    plt.barh(range(len(TPSI)), TPSI, color=MolColorDF_Eng['MolColor'],tick_label=MolLabel); plt.gca().tick_params(axis='both',labelsize=20) ; #plt.rcParams['axes.linewidth'] = 1.5
                    plt.xticks([-0.2,0.0,0.2,0.4,0.6],['-0.2','0.0','0.2','0.4','0.6']);   plt.tick_params(axis='y',labelsize=7)
        
                    plt.xticks(fontsize=20);plt.savefig(save_dir+'TPSI.pdf',bbox_inches="tight")
                    plt.close()
                    return(MolAveList)
                    
            TPSI=calcEachMolInterSubjCorr(self)
            return(TPSI)
        
    def AnalTVRI(self,save_dir):
        def mkNormalizeZscore(ThreedData,save_dir):
            timepointlist=self.timepointlist.copy()
            for i in timepointlist:
                if i==0:#
                    SubjCenter =    ThreedData.swaplevel(0).loc[i] - ThreedData.swaplevel(0).loc[i].fillna(0).mean()
                    tempZscore= SubjCenter / ThreedData.swaplevel(0).loc[0].fillna(0).std(ddof=1)
                elif i != -10:
                    SubjCenter =    ThreedData.swaplevel(0).loc[i] - ThreedData.swaplevel(0).loc[i].fillna(0).mean()
                    AddZscore= SubjCenter / ThreedData.swaplevel(0).loc[i].fillna(0).std(ddof=1)                   
                    tempZscore=pd.concat([tempZscore,AddZscore],axis=0) 
            ThreedZscore = pd.DataFrame(data=np.array(tempZscore),index=pd.MultiIndex.from_product([range(13), list(SubjCenter.index)]),columns=list(tempZscore.columns))

            tempZscore.to_excel(save_dir + 'SubjTmCs_zscored.xlsx')
            TVRI = 1-ThreedZscore.fillna(0).std(level=1,ddof=1).mean()
            TVRI.to_excel(save_dir+'TVRI.xlsx')
            return(TVRI)
        TVRI=mkNormalizeZscore(self.ThreedData,save_dir)
        self.Drawbar(TVRI,save_dir)
        return(TVRI)
        
        
    def AnalFig2_3(self):
        save_dir=self.save_dir+'Fig2_3FigS4_5/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir) 
        def plotTimeCourse(ClstMol,Colordict,TimeSignDF,EachSwitch,save_dir):#
            from matplotlib import pyplot as plt
            cluster_total = len(ClstMol)

            max_TimeSign = max(TimeSignDF.max(axis=1))
            min_TimeSign = min(TimeSignDF.min(axis=1))
            
            for ii in range(cluster_total):
        
                fig = plt.figure()
                fig.set_size_inches(6,4)
                plt.rcParams['axes.linewidth'] = 2.0
                
                target_dataname = ClstMol['cluster_' + str(ii+1)]
                num_target = len(target_dataname)
                plt.axhline(color="k",linewidth=2)
                handle_list = []
                ForClusterMean=np.empty((0,TimeSignDF.shape[0]),int)
                for ij in range(num_target):#
                  pl = plt.plot(list(TimeSignDF.index),TimeSignDF.loc[:,target_dataname[ij]],'-',lw=4,mfc='none',color = Colordict['cluster_' + str(ii+1)][0])
                  pl[0].set_mec(pl[0].get_color()) #
                  handle_list.append(pl[0])
                  if EachSwitch == 1:
                      fig = plt.figure()
                      fig.set_size_inches(6,4)  
                      pl = plt.plot(list(TimeSignDF.index),TimeSignDF.loc[:,target_dataname[ij]],'-',lw=4,mfc='none',color = 'k')
                      pl = plt.plot(list(TimeSignDF.index),[0]*len(list(TimeSignDF.index)),'-',lw=1,mfc='none',color = 'k')
        
                      pl[0].set_mec(pl[0].get_color()) 
                      handle_list.append(pl[0])
                      plt.xlabel('Time (min.)',size=18)
    
                      plt.ylim([min_TimeSign-0.1,max_TimeSign+0.1])
                      plt.show()
                      fig.savefig(save_dir + 'InClusterTimeCourse_' + target_dataname[ij] + '.pdf')
                      
                  if ij ==0:
                      ForClusterMean = np.append(ForClusterMean,np.array([TimeSignDF.loc[:,target_dataname[ij]]]),axis=0)
                  else:
                      ForClusterMean = np.append(ForClusterMean,np.array([TimeSignDF.loc[:,target_dataname[ij]]]),axis=0)
                ClusterMean=[]
                for i in range(0,len(ForClusterMean[0,:])):
                    ClusterMean.append(np.mean(ForClusterMean[:,i]))
                plt.plot(list(TimeSignDF.index),ClusterMean,color='k',lw=4)#axis=0
               
                plt.tick_params(labelsize=0)
        
                plt.xticks(rotation=270)
                plt.annotate('Cluster'+str(ii+1),xy=(45,0.75),size=50)

        
                plt.tick_params(labelbottom=False,
                        labelleft=False,
                        labelright=False,
                        labeltop=False)
                plt.tick_params(bottom=False,
                        left=False,
                        right=False,
                        top=False)
                if ii == 13:#
                    plt.xlabel('Time (min.)',size=100)
                    plt.xticks([0,60,120,240],['0','60','120','240'],size=100)
                    plt.tick_params(bottom=True)
                    plt.tick_params(labelbottom=True)

                plt.ylim([min_TimeSign-0.1,max_TimeSign+0.1])
                plt.show()

                fig.savefig(save_dir + 'ClusterTimeCourse' + str(ii+1) + '.pdf',bbox_inches="tight")
                plt.close()
                
        def timeCourseAvePCA(XDF,ClstAveTimeSeriesDF,ClstColorDF,ColorDF,ClstMolDF,ClstNoColorDF,save_dir):#PCAする

            AnalSwitch={'Biplot':1,'DrawEllipse':1,'PlotPCACovRatio':1,'LoadingPlotWCluster':1,
                        'ScoreHeatMap':1 ,'FactorLoadingHeatMap':1,
                        'LoadingPCsimulation':1,
                        'calcinner_outer_product' :0, 
                        'LengVector': 0,
                        'WHist':0
                }
            
            LoadingOpt = 0 
            BiplotSwitch={'DotColor':'ClstColor',
                          'Label':'Annotate',
                          'EachMolColor':[],
                          'Check':'',
                          'AminoCheck' :'',
                          'ColorDF':ColorDF,
                          'LoadingPCsimulation':'PC1_2_represen'
                          }#
            OptionSwitch={'Target':''}#
            OptionSwitch['MolColor']= pd.read_excel("./Data/LabelSummary.xlsx",header=0,index_col=0)
            CompareBolus={'TargetDetail' : ''}

            ##############################################
            MolLabel = list(pd.read_excel('./Data/LabelSummary.xlsx',header=0, index_col=0).index)
            
            labelProp =  XDF.columns
            MolLabel=XDF.index
        
            PCAPre.PCAprep(XDF,ClstAveTimeSeriesDF,ClstColorDF,MolLabel,labelProp,ClstMolDF,AnalSwitch,LoadingOpt,BiplotSwitch,OptionSwitch,CompareBolus,ClstNoColorDF,save_dir)
         
        Timepoint=[jj for jj in range(2,14)]        
        #Normalized timecourse
        YDF= self.mkNormalizeByStd(self.ThreedData,save_dir)
        YDF.to_excel(save_dir+'YDF.xlsx')
        Optindict={}
        Optindict['title']=''
        Optindict['Check'] = ''
        Optindict['AminoCheck'] = ''
        Optindict['Threshold']=3.2
        Optindict['numcluster']=13
        Optindict['metric']='euclidean'
        Optindict['MolColor']= self.MolColorDF

        ClstMol,ClstMolDF,Colordict,ClstDF,ClstColorDF,ColorDF,ClstNoColorDF = mHH.draw_heatmapclustering(YDF.T,1,'MolColor',Optindict,save_dir+'/',cmap='PuOr_r')#'bwr' ('YlOrRd'),'Reds'    
        ClstAveTimeSeriesDF=[]
        ClstAveTimeSeriesDF =mHH.mkClstAveTimeSign(save_dir,YDF.T,ClstMol)
        plotTimeCourse(ClstMol,Colordict,YDF,0,save_dir)

        timeCourseAvePCA(YDF.T,ClstAveTimeSeriesDF,ClstColorDF,ColorDF,ClstMolDF,ClstNoColorDF,save_dir+'/')#PCA

        Optiondict={'Color':'MolColor',
                'Annotate' : 1,
                }


    def AnalSTable7(self):
        save_dir=self.save_dir+'STable7/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir) 
        def mk3dDataSTable7(self):
            templist=self.timepointlist.copy()
            timepointlist=[-10]+templist
            for i in range(1,21):
                
                AddDF = self.TwodData.iloc[:,(0+14*(i-1)):(14+14*(i-1))].T     
                AddDF.index=timepointlist
                if i== 1: 
                    TmCsDF = self.TwodData.iloc[:,(0+14*(i-1)):(14+14*(i-1))].T
                    TmCsDF.index=timepointlist
                else:
                    TmCsDF=pd.concat([TmCsDF,AddDF],axis=0) 
            SubjTmCs3dData = pd.DataFrame(data=np.array(TmCsDF),index=pd.MultiIndex.from_product([range(20), list(AddDF.index)]),columns=list(TmCsDF.columns))
            SubjTmCs3dData=SubjTmCs3dData.rename_axis(['Subject','Time'], axis=0)
            return(SubjTmCs3dData)
        ThreedData=mk3dDataSTable7(self)
        MissingvalueDF = pd.DataFrame(index=[jj for jj in range(len(ThreedData.isnull().sum()))])
        MissingvalueDF['Molecule']=list(ThreedData.isnull().sum().index)
        MissingvalueDF['The percentage of missing points (%)']=list(round(ThreedData.isnull().sum()/(20*14)*100,1))

        DF = MissingvalueDF[MissingvalueDF!=0].dropna().sort_values(by='The percentage of missing points (%)')

        DF.index=[i+1 for i in range(len(DF.index))]
        DF.to_excel(save_dir+'Missingvalueratio.xlsx')
        
    def AnalFigS3(self,Data,file_dir,save_dir): 
        save_dir=self.save_dir+'FigS3/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir) 
        Timepoint=[jj for jj in range(2,14)]
        
        ### make Fasting values matrix
        MolLabel=list(Data.index)
        FastingDF=pd.DataFrame(data=None,index=range(20),columns=MolLabel)
        for j in range(1,21):
                FastingDF.iloc[j-1,:]=list(np.nanmean(Data.iloc[:,[14*(j-1),14*(j-1)+1] ],axis=1))
        
        NewtDF = pd.DataFrame(data=None,index=Timepoint,columns=MolLabel)
        NewpDF = pd.DataFrame(data=None,index=Timepoint,columns=MolLabel)
        NewFCDF = pd.DataFrame(data=None,index=Timepoint,columns=MolLabel)#FoldChange
        NewForVolDF = pd.DataFrame(data=None,index=None,columns=['label','p_val','ratio','s_val'])
        LabelList = []; PvalList=[]; FCList = []; s_valList=[];QvalList  =[]      
    
    
        for ii in Timepoint:
            TmPointRaw=Data.iloc[:,[ ii+ 14*(j-1) for j in range(1,21) ]].T
            for jj in MolLabel:
                if TmPointRaw[jj].isnull().any()==1:
                    tempFast = FastingDF.copy()
                    tempFast[jj][list(TmPointRaw[jj].isnull())] = np.nan
                    t, p = sts.ttest_rel(tempFast[jj].dropna(),TmPointRaw[jj].dropna())
                    FC = np.nanmean(list(TmPointRaw[jj].dropna())) / np.nanmean(tempFast[jj].dropna())
                else:
                    t, p = sts.ttest_rel(FastingDF[jj],TmPointRaw[jj])
                    FC = np.nanmean(list(TmPointRaw[jj])) / np.nanmean(FastingDF[jj])
                NewtDF.loc[ii,jj]=t
                NewpDF.loc[ii,jj]=p  
                NewFCDF.loc[ii,jj]=FC
                LabelList.append(str(ii)+'_'+jj)
                PvalList.append(p)
                FCList.append(FC)
                s_valList.append(0)
        QvalueStorey,QvalueBH=UseR(pd.DataFrame(PvalList),{'EngLabel':'FCEachTimevsFasting'})
        NewForVolDF['label']=LabelList; NewForVolDF['p_val']=PvalList;
        NewForVolDF['ratio']=np.log2(FCList); NewForVolDF['s_val'] = -1*np.log10(QvalueStorey)
            
        ##NewForVolDF.to_excel(save_dir+'DFForVolcano.xlsx')
        ##NewtDF.to_excel(save_dir+'tvalue_EachTmPoint.xlsx')
        ##NewpDF.to_excel(save_dir+'pvalue_EachTmPoint.xlsx')
        ##NewFCDF.to_excel(save_dir+'FCvalue_EachTmPoint.xlsx')
        self.AnalVolcanoPlot(NewForVolDF,file_dir,save_dir)
        return(NewForVolDF)

    def AnalVolcanoPlot(self,df,file_dir,save_dir):
            save_dirVP = save_dir+'VolcanoPlot/'
            if not os.path.isdir(save_dirVP):
                os.makedirs(save_dirVP)    
            plt.figure()
            volc = VP.Volcano(df["ratio"], df["p_val"], df["s_val"],df["label"], s_curve_x_axis_overplot=0.5, s_curve_y_axis_overplot=0.5)
            fig, DF = volc.get_fig(self.MolColorDF)
            plt.rcParams['axes.linewidth'] = 2.0
            MolColor=pd.read_excel(file_dir+'LabelSummary.xlsx',header=0,index_col=0)
            Volcano = pd.concat([MolColor,DF],axis=1,join='inner')
            LabelcountDF = Volcano[['Up','Down']]

            plt.savefig(save_dirVP+'/VolcanoPlot.pdf')
            ##LabelcountDF.to_excel(save_dirVP+'Labelcount.xlsx')
            plt.figure()
            volc = VP.Volcano(df["ratio"], df["p_val"], df["s_val"],df["label"], s_curve_x_axis_overplot=0.5, s_curve_y_axis_overplot=0.5)
            fig = volc.get_DistrFC()
            plt.rcParams['axes.linewidth'] = 2.0

            plt.savefig(save_dirVP+'/DistributionOfFoldChange.pdf')
            plt.figure()
            
            VP.DrawPieChart(LabelcountDF,file_dir,save_dirVP)
    def Drawbar(self,DrawDF,save_dir):       
        MolColor = self.MolColorDF
        DrawDF = pd.concat([DrawDF,MolColor],join='inner',axis=1)
        List1 = np.array(list(DrawDF.sort_values(by=0)[0]))
        
        XtickLabel=list(DrawDF.sort_values(by=0).index)
        MolColor=list(DrawDF.sort_values(by=0)['MolColor'])
        fig,ax = plt.subplots(figsize=(4,9.6))
        ax.barh(np.linspace(0, len(List1), len(List1)),List1,tick_label=XtickLabel,color=MolColor)
        plt.gca().tick_params(axis='both',labelsize=20)
        ax.tick_params(axis='y',labelsize=7)
        ax.set_xticks([0, 0.25, 0.5, 0.75])
        ax.set_xticklabels(['0.0', '0.25', '0.50', '0.75'])
    
        plt.savefig(save_dir+'TVRI.pdf',bbox_inches="tight") 
        ##plt.savefig(save_dir+'TVRI.png',bbox_inches="tight") 
    
        plt.close()    

    def adjustMatrupper(self,matr):
        collength =len(matr.columns)
        rowlength = len(matr.index)
        for i in range(0,collength):
            matr.iloc[i,i]=np.nan
            for j in range(i+1,rowlength):
                matr.iloc[i,j] = np.nan
                
        return(matr)
        
    def adjustMatrlower(self, matr):
        collength =len(matr.columns)
        rowlength = len(matr.index)
        for i in range(0,collength):
            matr.iloc[i,i]=np.nan
            for j in range(i+1,rowlength):
                matr.iloc[j,i] = np.nan
                
        return(matr)
        
def Preprocess(file_dir):
    Data=pd.read_excel(file_dir+'Metabolites_and_hormones_data.xlsx',header=0,index_col=0,sheet_name='Raw (glucose) ')    
    DelList=['ID','Unit']
    Data=Data.drop(DelList,axis=1)
    Data=Data[Data['Analysed molecule']==1]
    Data=Data.drop('Analysed molecule',axis=1)

    return(Data)

    
