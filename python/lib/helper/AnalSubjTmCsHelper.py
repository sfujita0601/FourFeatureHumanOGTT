#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 21:39:14 2018

AnalSubjTmCsのHelper
@author: fujita
"""


import numpy as np
import pandas as pd
import scipy.io
import itertools
import sys
import matplotlib as mpl
import os
from Helper import general as ge
import matplotlib.pyplot as plt
from DataLoader import mkMatForModel as MM4Model
from Helper import StatCal as SC
import matplotlib.cm as cm
import Helper.AnalPropParamHelper as APPHel
import Helper.PCA_ChgColor_After030118 as HPCA
import mkdendrogram as Clst 
from scipy.stats import pearsonr, spearmanr
import Helper.AnalMolTime as AmT
import Helper.PCAPre as PCAPre
import matplotlib.cm as cm
import DataLoader.mkBolusIdx as BI
import warnings 
from DataLoader import AdjustData as aD
from scipy.stats import zscore
import Helper.AnalMolTime as AmT
import MatrHelper as MH
from matplotlib.colors import LinearSegmentedColormap
import Helper.IncretinWorkHelper as IWHel
import Helper.MolCorrSubjmeanStdHelper as MSMSH
import re
import Helper.GraphHelper as GH
import Helper.mkNetWorkGraphHelper as NWGraph
import scipy.stats as sts
import DataLoader.ChangeNameEng as CNE
import scipy
import Helper.mkHeatmapHelper as mHH
import AnalPropParam as APP
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import mkdendrogram as mD
import Helper.AdjustClstDF as ACD
from Helper import ContinuousWorkHelper as CWH
import Helper.AdjustClstDF as ACD
import Helper.LabelHeler as LH


color_map = LinearSegmentedColormap('color_map',
    {'red': [(0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,1.0,1.0)], #(x,y0,y1) y1→y0y1→y0...
   'green': [(0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,0.0,0.0)],
    'blue': [(0.0,1.0,1.0),(0.5,1.0,1.0),(1.0,0.0,0.0)]})




def Premkvenn(save_dir):#ベン図を描画する準備
    
    from matplotlib_venn import venn2, venn3
    from matplotlib import pyplot

    #plt.figure(facecolor="azure", edgecolor="coral", linewidth=10)    
    pyplot.rcParams['axes.linewidth'] = 0.0# 軸の線幅edge linewidth。囲みの太さ
    
    # 重なる部分の領域の割合を0にする   
    v = venn2(subsets=(0, 13, 23), set_labels = ('', ''))
    
    # ベン図の色を変更する   
    v.get_patch_by_id('10').set_color('skyblue')
    v.get_patch_by_id('10').set_edgecolor('black')
    v.get_label_by_id('100').set_text('')
    #v.get_patch_by_id('11').set_color('white')
    v.get_patch_by_id('11').set_edgecolor('black')
    v.get_label_by_id('11').set_text('')
    #v.get_patch_by_id('01').set_color('white')
    v.get_label_by_id('01').set_text('')
    v.get_patch_by_id('01').set_edgecolor('black')

    # 背景色を変更する
    pyplot.gca().set_axis_on()
    pyplot.gca().set_facecolor('white')

    plt.savefig(save_dir+'venn.pdf')
    pyplot.show()
 

def PlotScatterDF(DF1,DF2,save_dir,logSwitch):#DF2の1,2列目の組みをDF1の描画
    DF2Col = list(DF2.columns);numfig = len(DF2[DF2Col[0]])
    print(DF2[DF2Col[0]][0]);print(DF2[DF2Col[0]][0])
    num_cols=6#len(DataType)*colmol*avestd #1行あたりのグラフ数
    num_rows=6
    num_page = int(numfig / (num_cols*num_rows))
    #mod_mol = num_mol*len(DataType)%(num_cols*num_rows)
    #DF1 = [DF1.drop(DF1.index[DF1[DF2[DF2Col[0]][x]] > 900000]) for x in range(0,len(DF2[DF2Col[0]]))]    #count=0
    #外れ値ぽいのを抜く？
    #idx =['kinnmedai','buri','iwana','isaki']; DF1 = DF1.drop(idx,axis=0)
    count=0
    for k in range(0,num_page+1):
        fig, host = plt.subplots(num_rows,num_cols,figsize=(50,50))
        plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05)
        for j in range(0,num_rows): #行番号,一行につき1被験者   
            for i in range(0,num_cols):
                    if count < numfig:#                y=[[]]*len(BolusNoList)
                        if logSwitch==0:#普通の
                            r, p=pearsonr(list(DF1[DF2[DF2Col[0]][count]]),list(DF1[DF2[DF2Col[1]][count]]))
                            host[j,i].scatter(list(DF1[DF2[DF2Col[0]][count]]),list(DF1[DF2[DF2Col[1]][count]]),s=600,c='b')
                        elif logSwitch==1: #Xをlog10に
                            r, p=pearsonr(list(DF1[DF2[DF2Col[0]][count]]),list(DF1[DF2[DF2Col[1]][count]]))
                            host[j,i].scatter(list(DF1[DF2[DF2Col[0]][count]]),list(DF1[DF2[DF2Col[1]][count]]),s=600,c='r')
                        elif logSwitch==2:#Xが6を取らないように
                            X = np.log10(list(DF1[DF2[DF2Col[0]][count]]))
                        elif logSwitch==3: #散布図：生を赤に、負を青に
                            r, p=pearsonr(list(DF1[DF2[DF2Col[0]][count]]),list( DF1[DF2[DF2Col[1]][count]] ) ) 
                            r, p=spearmanr(list(DF1[DF2[DF2Col[0]][count]]),list( DF1[DF2[DF2Col[1]][count]] ) ) 
                            
                            
                            host[j,i].scatter(  list(DF1[DF2[DF2Col[0]][count]])  ,list(DF1[DF2[DF2Col[1]][count]]) ,s=600,c='r')                           
                            if r <0:
                                host[j,i].scatter( list(DF1[DF2[DF2Col[0]][count]]) ,list(DF1[DF2[DF2Col[1]][count]]) ,s=600,c='b')
                        ##   host[j,i].scatter(list(DF1[DF2[DF2Col[0]][count]]),list(DF1[DF2[DF2Col[1]][count]]),s=600,c='r')
               
                        ###分子名だけ取り出す                    
                        r = re.compile("(.*)(_)(.*)") 
                        d = r.search(DF2[DF2Col[0]][count]) 
                        host[j,i].set_xlabel(d.group(3),fontsize=40)                
                        host[j,i].set_ylabel(DF2[DF2Col[1]][count],fontsize=40)#,fontproperties=prop)  
                        
                        p = DF2['PVal'][count]
                        r = DF2['Corr'][count]
                        Roundp = round(p,10)
                        b = '%.2e'%Roundp
                        TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                        Roundq = round(DF2['QVal'][count],100)
                        bq = '%.2e'%Roundq
                        TenRoundq = bq[0:bq.find('e')] + '$\it{×10^{-'+str(bq[len(bq)-1])+'}}$'
                        #p_q$\mathsf{ hoge }$
                        host[j,i].set_title('$\it{R}$=' + str(round(r,3)) + ', $\it{p}$=' + TenRoundp + '\n$it{q}$=' +TenRoundq ,fontsize=40)
                        #q値のみ
                        host[j,i].set_title('$\it{R}$=' + str(round(r,3)) + ', $\it{q}$=' +TenRoundq ,fontsize=40)

                        host[j,i].tick_params(labelsize=30)
                        count+=1
                    else:
                        host[j,i].spines['right'].set_visible(False)
                        host[j,i].spines['top'].set_visible(False)
                        host[j,i].spines['left'].set_visible(False)
                        host[j,i].spines['bottom'].set_visible(False)
                        host[j,i].tick_params(labelleft="off",left="off") # y軸の削除
                        host[j,i].tick_params(labelbottom="off",bottom="off") # y軸の削除

        fig.tight_layout()
        plot_axis = plt.axis()
        plt.savefig(save_dir + str(k)+'.pdf',format='pdf')
    print('output:'+ save_dir)
    
def PreFastResponsePlotScatterDF(save_dir):
    
     CombDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/RespIdx/CombDF_rev_sorted.xlsx' ,header=0,encoding = "ISO-8859-1")
     #CombDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/RespIdx/Comb_PropProp.xlsx',header=0,encoding = "ISO-8859-1")
     DF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/RespIdx/SubjMolRespIdx_All_Prop.xlsx' ,header=0,encoding = "ISO-8859-1")
     #Comb_PropProp
     PlotScatterDF(DF,CombDF,save_dir ,3)#1で全部赤、3で#散布図：生を赤に、負を青に

def ChangeCombDF(Label,save_dir):#CombDFの順番を整形
     r = re.compile("(.*)(_)(.*)")     
     CombDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/RespIdx/Comb_Fasting_AUC_Eachother_Prop_ref付き_latest_relateProp.xlsx' ,header=0,encoding = "ISO-8859-1")
     #並びの点数をつけていいく：Fasting(1), AUC(2), Label順、その後Prop順
     orderList=[]
     PropLabelDict = dict({'Age':1, 'Height':2, 'Weight':3, 'Waist':4, 'BMI':5, 'FastingGlucose':6,'TwoHGlucose':7, 'FastingInsulin':8, 'TwoHInsulin':9})
     for i in range(len(list(CombDF['Mol1']))):  
         d = r.search(CombDF['Mol1'][i]) 
         if 'Fasting' == d.group(1):   #空腹値ならはじめ 
             orderIdx = 1
         else:#応答量なら後
             orderIdx = 100
         LabelIdx = Label.index(d.group(3))#分子ごとに順番振る
         PropIdx = PropLabelDict[CombDF['Mol2'][i]]#プロパティの順番
         orderList.append(orderIdx**2*(LabelIdx+1)**1.5+PropIdx)
     CombDF['order'] = orderList; print(CombDF['Mol1']); CombDF.to_excel(save_dir+'CombDF_rev.xlsx')
     CombDF_rev = CombDF.sort_values('order')
     CombDF_rev.to_excel(save_dir+'CombDF_rev_sorted.xlsx')
    
def calcp_FCval(SubjectName,save_dir):#空腹値とのt検定_各時点で
    
    file_dir ='/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Delta_Eng/New/'
    kanpachi = pd.read_excel(file_dir +'SubjTmCs_'+SubjectName[0]+'Delta.xlsx'  ,header=0,encoding = "ISO-8859-1",index_col=0)

    SubjFasting_Eng = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/SubjFasting_Eng.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    #QvalDF  = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181108/QvalueStoreyParamACDE_ParamACDE.xlsx', header=0,encoding = "ISO-8859-1")
    #QvalDF  = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20190323/TimePointSubjMol/AllTimePoint/ParamACDE_ParamACDE__2Mol_pearson/QvalueStoreyParamACDE_ParamACDE.xlsx', header=0,encoding = "ISO-8859-1")
    #QvalDF  = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20190425/TimePointSubjMol_old/DFForVolcano_qlamda.xlsx', header=0,encoding = "ISO-8859-1")
    Timepoint = kanpachi.index.tolist()[1:]
    #3-HBを除く
    Mol = kanpachi.columns.tolist()
#    Mol.remove('3-Hydroxybutyrate')
    NewtDF = pd.DataFrame(data=None,index=Timepoint,columns=Mol)#t検定のt値
    NewpDF = pd.DataFrame(data=None,index=Timepoint,columns=Mol)#t検定のp値
    NewFCDF = pd.DataFrame(data=None,index=Timepoint,columns=Mol)#FoldChange
    NewSignChangeDF =  pd.DataFrame(data=np.zeros([len(Timepoint),len(Mol)]),index=Timepoint,columns=Mol)#有意に変動したポイント
    
    NewForVolDF = pd.DataFrame(data=None,index=None,columns=['label','p_val','ratio','s_val'])
    LabelList = []; PvalList=[]; FCList = []; s_valList=[];QvalList  =[] 
    r = re.compile("(.*)(_)(.*)")     
    #dlist =[r.search(QvalDF.index[i]) for i in range(len(QvalDF.index))]
    
    for ii in Timepoint:#save_dirから各時点を取ってくる

        TmPointDelta = pd.read_excel(save_dir+str(ii)+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        TmPointRaw = SubjFasting_Eng + TmPointDelta#Delta後のDF
        for jj in Mol:#分子で回す
            if TmPointRaw[jj].isnull().any()==1:#Series内にnanがあるなら
                tempFast = SubjFasting_Eng.copy()
                tempFast[jj][TmPointRaw[jj].isnull()] = np.nan#ここでは対象時点と合わせるために
                t, p = sts.ttest_rel(tempFast[jj].dropna(),TmPointRaw[jj].dropna())
                
                
            else:
                t, p = sts.ttest_rel(SubjFasting_Eng[jj],TmPointRaw[jj])
            """%%%片側検定？
            pval3 = p
            pval2 = p / 2.0
            pval1 = 1.0 - pval2
            if t < 0.0:
                w = pval2
                pval2 = pval1
                pval1 = w
            FC = np.nanmean(TmPointRaw[jj]) / np.nanmean(SubjFasting_Eng[jj])
            """ 


                        

            NewtDF.loc[ii,jj]=t
            NewpDF.loc[ii,jj]=p  
            NewFCDF.loc[ii,jj]=FC
            
            
            #d = r.search(ii) 
            #if MolLabelList[jj] in ii:            
            #if (ii==dlist[0].group(1) && (jj==dlist[0].group(1)
            
            #q=QvalDF.loc[str(ii)+'_'+jj]['p_val']
            #if ((q<0.1) and (FC > 1.5)) or ((q<0.1) and (FC < 0.67)):
            #    NewSignChangeDF.loc[ii,jj] = 1
            
            #QvalList.append(q)
            LabelList.append(str(ii)+'_'+jj)
            PvalList.append(p)
            FCList.append(FC)
            s_valList.append(0)
    #####20181108=Q値に変えた
    
    NewForVolDF['label']=LabelList; NewForVolDF['p_val']=PvalList;#QvalList;#
    NewForVolDF['ratio']=FCList; NewForVolDF['s_val'] = s_valList
    #APP.mkHist(NewForVolDF['p_val'],'pvalue',save_dir,xticks=[0,0.2,0.4,0.6,0.8,1.0],xticklabels=['0','0.2','0.4','0.6','0.8','1.0'])#ヒストグラムを作る
        
    NewForVolDF.to_excel(save_dir+'DFForVolcano.xlsx')
    NewtDF.to_excel(save_dir+'tvalue_EachTmPoint.xlsx')
    NewpDF.to_excel(save_dir+'pvalue_EachTmPoint.xlsx')
    NewFCDF.to_excel(save_dir+'FCvalue_EachTmPoint.xlsx')
    NewSignChangeDF.to_excel(save_dir+'SignChange_EachTmPoint.xlsx')
    
def mkTimePointSubjMolxlsx(SubjectName,save_dir):#各時点で、被験者×分子のDF作る
    
    file_dir ='/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Delta_Eng/New/'
    kanpachi = pd.read_excel(file_dir +'SubjTmCs_'+SubjectName[0]+'Delta.xlsx'  ,header=0,index_col=0)
    SubjTmCsDeltaPanel = pd.Panel({SubjectName[0]:pd.DataFrame(data=kanpachi,index=list(kanpachi.index),columns=list(kanpachi.columns))})
    #NewDF = pd.DataFrame(data=None,index=SubjectName,columns=SubjTmCstemp.columns.tolist())
    for i in SubjectName:
           SubjTmCsDelta = pd.read_excel(file_dir +'SubjTmCs_'+i+'Delta.xlsx'  ,header=0,index_col=0)
           SubjTmCsDeltaPanel[i] = SubjTmCsDelta


    for ii in SubjTmCsDelta.index.tolist():
       NewDF = SubjTmCsDeltaPanel.major_xs(ii)        
       NewDF.T.to_excel(save_dir+str(ii)+'.xlsx')
         
            

def givetoGraph(Corr,Qvalue,save_dir):
    ColorSwitchDict  = {1:'ClstColor',2:'MolColor',3:'TimeVarColor'}#1でクラスタ、2で分子,3で時間ばらつき
    GraphDict={}
    GraphDict['pos']='Left'#座標を陽に与える、'positive'、計算する：negative:#AAのみサークルで描画：'AA'、#サークルで描画：'circle'、:#斥力で離す'Left'、　:#PropMetaboMolTimeなら：'Global'
    GraphDict['TargetProp']=''#BMIしか描かない
    GraphDict['QvalAnal'] = 1#Q値を元にした解析を行う：1
    GraphDict['GaraphSwitch'] = '' # 'MolTimePCProp'：MolTimeのPCとPropでの相関結果を元にグラフ描画  'MoltimeProp':分子 ×時間 vs property
    GraphDict['GaraphSwitch_MolTime_PCA_Prop'] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/MolTimePCA_vs_Propert/Corr__W_pearson_CombCoef.xlsx' ,header=0,encoding = "ISO-8859-1")

    GraphDict['GaraphSwitch_MolTime_PCA_Prop_Q'] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/MolTimePCA_vs_Propert/QvalueStorey_CombCoef.xlsx' ,header=0,encoding = "ISO-8859-1")
    GraphDict['GaraphSwitch_MolTime_PCA_Prop_CoefThresh']=0.70
    GraphDict['GaraphSwitch_MolTime_PCA_Prop_Abs'] = 'posneg'
    #GraphOptionDict['Edge']=''#'Subject'}    #''で相関係数の大きさで線の太さを変える、'Subject'で被験者の数だけ、'Global：
    #GraphOptionDict['mkEdge'] = ''#閾値を基準にグラフを描画する場合：'Thresh'、'Comb':#何種類かのファイルを元に書かれた組み合わせの線をひく
    #GraphOptionDict['pos']='' #座標を陽に与える、'positive'、計算する：negative:#AAのみサークルで描画：'AA'、#サークルで描画：'circle'、:#斥力で離す'Left'、　:#PropMetaboMolTimeなら：'Global'
    #GraphOptionDict['Adjacency_matrix']= ''#pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180930/temp_PropMetabo_rev.xlsx' ,header=0,encoding = "ISO-8859-1")
    #QValDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/MolTimePCA_vs_Propert/QvalueStorey_.xlsx' ,header=0,encoding = "ISO-8859-1")
    #CorrDF_Eng = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/MolTimePCA_vs_Propert/Corr__W_pearson.xlsx' ,header=0,encoding = "ISO-8859-1")
    #NWGraph.mkGraph(CorrDF,0.75,4,ColorSwitch,save_dir)#エッジの数を計算して、グラフを描画して,次数分布と次数中心性のグラフも、中心性の高いとこだけ描画も
    NWGraph.mkGraphQval(Corr,Qvalue,0.1,2,ColorSwitchDict[2],GraphDict,save_dir)#Qvalue<0.1イカのやつ

def setpValCorr(SubjectName,ResponseIdxDict, save_dir):
    r = re.compile("(.*)(_)(.*)") 

    #'Opeak', 'Tpeak', 'Ttimeconctant','FC','Gain', 'TPI','Adaptation Precision', 'AUC','Fasting'
    #特定の特徴量のみ
    DelIdx=[ 'Fasting','AUC' ]
    name = DelIdx[0]
    if len(DelIdx) > 1:
        name = []
        for s in DelIdx:
            name+=s
        name=name[0]
    #もしファイルが別にあるなら
    #pearson:Pval=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181023/RespIdx/Fasting_AUC_Gain/Pvalue_Wpearson_NoSeparate.xlsx',header=0,encoding = "ISO-8859-1")
    #pearson:Corr=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181023/RespIdx/Fasting_AUC_Gain/Corr_Wpearson_NoSeparate.xlsx',header=0,encoding = "ISO-8859-1")
    #spearman
    #Pval=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181205/RespIdx/Pvalue_Wpearson_NoSeparate_FastingAUCProp.xlsx',header=0,encoding = "ISO-8859-1")
    #Corr=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181205/RespIdx/Corr_Wpearson_NoSeparate_FastingAUCProp.xlsx',header=0,encoding = "ISO-8859-1")
    
    #List1 = [list(Pval.index)[k]  for k in range(len(list(Pval.index))) if r.search(list(Pval.index)[k]).group(1) in DelIdx ] 
    #PvalNew = Pval.loc[List1]
    #CorrNew = Corr.loc[List1]
    
    #PvalNew.to_excel(save_dir+'Pvalue_Wpearson_'+'W_'+name+'.xlsx')
    #CorrNew.to_excel(save_dir+'Corr_Wpearson_'+'W_'+name+'.xlsx')

    #pearsonで：Qvalue = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181023/Fasting_AUC_Gain/QvalueStoreyParamACDE_ParamACDE.xlsx',header=0,encoding = "ISO-8859-1")
    #spermanで
    #Qvalue = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181205/RespIdx/QvalueStoreyParamACDE_ParamACDE.xlsx',header=0,encoding = "ISO-8859-1")
    #RespIdxProp = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181205/RespIdx/SubjMolRespIdx_FastingAUCProp.xlsx',header=0,encoding = "ISO-8859-1")
    Pval=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191211/RespIdx/Pvalue_Wpearson_NoSeparate.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    Corr=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191211/RespIdx/Corr_Wpearson_NoSeparate.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    Qvalue = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191211/Opeak_Fasting/QvalueStoreyOpeak_Fasting_c.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    RespIdxProp = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191211/RespIdx/SubjMolRespIdx.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    Qvalue.columns = list(Corr.columns)
    CombDF = APPHel.mkExcellower2(Corr,Pval,Qvalue,1)
    CombDF.to_excel(save_dir+'CombDF2.xlsx')
    CombDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191211/RespIdx/CombDF.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)

    APPHel.PlotScatterDF(RespIdxProp,CombDF,save_dir ,3)#1で全部赤、3で#散布図：生を赤に、負を青に
    #givetoGraph(Corr,Qvalue,save_dir)    
    
def calcCorrBetwAUCVartvsRespIdx(save_dir):#特徴量と変動量指標のAUCの相関
    VarIdx=list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181025/RespIdx/SubjIdx_Z.xlsx',header=0,encoding = "ISO-8859-1")['Ttimeconctant'])
    RespIdx=list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181023/RespIdx/AllResp/SubjIdx_Z.xlsx',header=0,encoding = "ISO-8859-1")['Ttimeconctant'])
    r,p = sts.pearsonr(VarIdx,RespIdx)#ピアソン
    #t.to_excel(save_dir + 'Corr'+filename+'.xlsx'); Pval.to_excel(save_dir + 'Pvalue'+filename+'.xlsx');
    #tUpper = MH.adjustMatrlower(t);  PvalUpper = MH.adjustMatrlower(Pval) #対角と下三角をNaNにする
    #tUpper.to_excel(save_dir + 'CorrUpper'+filename+'.xlsx'); PvalUpper.to_excel(save_dir + 'PvalueUpper'+filename+'.xlsx');
    fig=plt.figure(figsize=(5,5))
    plt.scatter(VarIdx,RespIdx,s=8)#, color=MolColor)
    Roundp = round(p,200)
    b = '%.2e'%Roundp
    TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                        
    #plt.plot([0.835,1.0],[0.95, 0.95], c='black')
    #plt.xlim([0.835,1.0])
    
    plt.title('R=' + str(round(r,3)) + ', p<0.01' ,fontsize=10)
    plt.savefig(save_dir +  'TimeConctant_VartvsResp.pdf',format='pdf')

    
def calcpearsonDF(SubjectName,ResponseIdxDict, save_dir):
    RespIdxMolDict={}; NameList=['Fasintg','AUC','AUChalf','AUCPeak','AUCPeakhalf']#'Opeak', 'Tpeak', 'Ttimeconctant','FC','Gain', 'TPI','Adaptation Precision', 'AUC','Fasting' ]
    RespIdxDF = pd.read_excel(save_dir+'SubjMolRespIdx.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    LabelSum =   pd.read_excel("//Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx",header=0,encoding = "ISO-8859-1",index_col=0)

    if ResponseIdxDict['Target_Calcpearson'] in NameList:
        RespIdxDF = pd.read_excel(save_dir+ResponseIdxDict['Target_Calcpearcon']+'.xlsx',header=0,encoding = "ISO-8859-1")
        
    r = re.compile("(.*)(_)(.*)") 
    if ResponseIdxDict['Target_Calcpearson'] == 'Clst6':#変動のおおきかったやつだけにする
        Clst6= list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolRepresentupperClst6_Eng.xlsx',header=0,encoding = "ISO-8859-1").index)

        NewDF = pd.DataFrame(data=None,index=SubjectName,columns= [])

        for i in Clst6:
            List1 = [list(RespIdxDF.columns)[k]  for k in range(len(list(RespIdxDF.columns))) if i == r.search(list(RespIdxDF.columns)[k]).group(3) ] 
            temparray = np.array( RespIdxDF[ List1 ] ); AddDF = pd.DataFrame(data=temparray,index=SubjectName,columns= List1)
            NewDF=pd.concat([NewDF,AddDF],axis=1, join_axes=[NewDF.index])
        NewDF.to_excel(save_dir+'SubjIdx_Clst6.xlsx') 
        RespIdxDF =   NewDF 
    #特定の特徴量を覗く
    #DelIdx='Opeak'
    #for i in list( RespIdxDF.columns ):
     #   #for j in NameList: #特徴量で回す   
      #  if r.search(i).group(1) ==DelIdx:
       #     RespIdxDF.drop(i,axis=1)     
            
    if ResponseIdxDict['WOGender']=='WOgender':#性別除く
        SubjProp = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/SubjPropWOgender.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        GenderProp = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjPropWgender_fujita.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)   
    else:
        SubjProp = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjPropWgender_fujita.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    RespIdxProp = pd.concat([RespIdxDF,SubjProp],axis=1); 
    countIdx=1
    filename='_NoSeparate'
    if ResponseIdxDict['WOGender_Corr']=='Separate': #男女でわけるなら
        RespIdxProp = pd.concat([RespIdxDF,SubjProp],axis=1); 
        RespIdxPropWomen = RespIdxProp[GenderProp['Gender']==0]
        RespIdxProp = RespIdxProp[GenderProp['Gender']==1]
        RespIdxProp.to_excel(save_dir+'SubjMolRespIdx_man.xlsx')
        RespIdxPropWomen.to_excel(save_dir+'SubjMolRespIdx_woman.xlsx')
        
        filename='_man'
        countIdx=2
    for ii in range(countIdx):#男女で分けたら2回行う
        if ii ==1:
            RespIdxProp=RespIdxPropWomen
            filename='_woman'
        RespIdxProp = RespIdxProp.drop(['BMI','Height','Weight','Age','Waist','FastingGlucose','FastingInsulin','TwoHGlucose','TwoHInsulin','Gender'],axis=1)
     #   MolName= list(RespIdxDF.index)
        NewDF = pd.DataFrame(data=None,index=SubjectName)
        print('Ongoing')
    
        """
        t,Pval = SC.calculate_CorrWpearson(RespIdxProp)#ピアソン
    
        t.to_excel(save_dir + 'Corr'+filename+'.xlsx'); Pval.to_excel(save_dir + 'Pvalue'+filename+'.xlsx');
        tUpper = MH.adjustMatrlower(t);  PvalUpper = MH.adjustMatrlower(Pval) #対角と下三角をNaNにする
        tUpper.to_excel(save_dir + 'CorrUpper'+filename+'.xlsx'); PvalUpper.to_excel(save_dir + 'PvalueUpper'+filename+'.xlsx');
       
        """
        if ResponseIdxDict['Aim']=='Calcpearson':#ピアソンで
            Corr,Pval = SC.calculate_CorrWpearson(RespIdxProp)#ピアソン
        elif ResponseIdxDict['Aim']=='Calcspearman':#スピアマンで
            Corr,Pval = SC.calculate_CorrWspeaman(RespIdxProp)#スピアマン
        
        #LastName = list(Corr.index)[ np.where(np.array([Corr.index=='BMI']))[1][0]-1]       
        #Corr = Corr.loc[:LastName,'BMI':];Pval = Pval.loc[:LastName,'BMI':];#NameList[j] += '_Separate'
        
        Corr.to_excel(save_dir + 'Corr_Wpearson'+filename+'.xlsx'); Pval.to_excel(save_dir + 'Pvalue_Wpearson'+filename+'.xlsx');
        CorrUpper=Corr.copy(); PvalUpper=Pval.copy()
        CorrUpper = MH.adjustMatrlower(CorrUpper);  PvalUpper = MH.adjustMatrlower(PvalUpper) #対角と下三角をNaNにする
        CorrUpper.to_excel(save_dir + 'CorrUpper'+filename+'.xlsx'); PvalUpper.to_excel(save_dir + 'PvalueUpper'+filename+'.xlsx');
        plt.hist(Corr.values[~np.isnan(CorrUpper.values)]);plt.title('_Wpearson'); plt.savefig(save_dir + '_DistOfCorr_Wpearson'+filename+'.pdf');plt.close()
        plt.hist(Pval.values[~np.isnan(PvalUpper.values)]);plt.title( '_Wpearson'); plt.savefig(save_dir + '_DistOfPval_Wpearson'+filename+'.pdf');plt.close()
    
        #if 'Prop' in NameList[j]: #and 'CutOff' not in NameList[j]:#プロパティを含んだ行列の時に列方向をプロパティだけにしたい
        #print(list(Corr.index)[ np.where(np.array([Corr.index=='AUC_Glu + threo-beta-methylaspartate']))])
        
        
        #もしファイルが別にあるなら
        #Pval=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181010/RespIdx/Pvalue_Wpearson.xlsx',header=0,encoding = "ISO-8859-1")
        #Corr=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181010/RespIdx/Corr_Wpearson.xlsx',header=0,encoding = "ISO-8859-1")                
        plt.rcParams["font.size"] = 10
### temp_20191211
        def temp20191211(Pval,Corr):        
            Pval = Pval.iloc[0:19,19:]            
            #Corr.to_excel(save_dir + 'Corr_Wpearson'+filename+'.xlsx'); Pval.to_excel(save_dir + 'Pvalue_Wpearson'+filename+'.xlsx');
            #plt.hist(Corr.values[~np.isnan(Corr.values)]);plt.title('_Wpearson'); plt.savefig(save_dir + '_DistOfCorr_Wpearson'+filename+'.pdf');plt.close()
            #plt.hist(Pval.values[~np.isnan(Pval.values)]);plt.title( '_Wpearson'); plt.savefig(save_dir + '_DistOfPval_Wpearson'+filename+'.pdf');plt.close()    
            Pval2 = pd.DataFrame(data=None,index=list(Pval.index),columns=list(Pval.columns))
            for i in range(len(Pval.columns)):
                Pval2.iloc[i,i] = np.diag(Pval)[i] 
            Corr2 = pd.DataFrame(data=None,index=list(Pval.index),columns=list(Pval.columns))
            for i in range(len(Pval.columns)):
                Corr2.iloc[i,i] = np.diag(Corr)[i] 
            Corr2.to_excel(save_dir + 'Corr_Wpearson_20191211'+filename+'.xlsx'); Pval2.to_excel(save_dir + 'Pvalue_Wpearson_20191211'+filename+'.xlsx');
        def temp20200511_Ins(Pval,Corr):        
            Pval2 = Pval[['AUC_Insulin','AUChalf_Insulin']].iloc[13:]   #Insulinの2つの特徴量と対応する他の分子x特徴量行列         
            Pval2['AUC_Insulin']['AUC_Insulin']=np.nan;Pval2['AUChalf_Insulin']['AUChalf_Insulin']=np.nan;
            Corr2 = Corr[['AUC_Insulin','AUChalf_Insulin']].iloc[13:]   #Insulinの2つの特徴量と対応する他の分子x特徴量行列         
            Corr2['AUC_Insulin']['AUC_Insulin']=np.nan;Corr2['AUChalf_Insulin']['AUChalf_Insulin']=np.nan;
            SC.UseR(Pval2,{'EngLabel':'Insulin'+filename})
            #csvをxlsxへ
            df_new = pd.read_csv("/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200524/Insulin"+filename+"/QvalueStoreyInsulin"+filename+"_c.csv", header=0,encoding = "ISO-8859-1",index_col=0)            
            df_new.to_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200524/Insulin'+filename+'/QvalueStoreyInsulin'+filename+'_c.xlsx')            
            Qvalue = df_new#pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200511/Insulin/QvalueStoreyInsulin_c.xlsx', header=0,encoding = "ISO-8859-1",index_col=0)
            return(Pval2,Corr2,Qvalue)  
        def temp20200511_Glc(Pval,Corr):        
            Pval2 = Pval[['AUC_Glucose','AUChalf_Glucose']].iloc[19:]   #Insulinの2つの特徴量と対応する他の分子x特徴量行列         
            Pval2['AUC_Glucose']['AUC_Glucose']=np.nan;Pval2['AUChalf_Glucose']['AUChalf_Glucose']=np.nan;
            Corr2 = Corr[['AUC_Glucose','AUChalf_Glucose']].iloc[19:]   #Insulinの2つの特徴量と対応する他の分子x特徴量行列         
            Corr2['AUC_Glucose']['AUC_Glucose']=np.nan;Corr2['AUChalf_Glucose']['AUChalf_Glucose']=np.nan;
            SC.UseR(Pval2,{'EngLabel':'Glucose'})
            #csvをxlsxへ
            df_new = pd.read_csv("/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200511/Glucose"+filename+"/QvalueStoreyGlucose_c.csv", header=0,encoding = "ISO-8859-1",index_col=0)            
            df_new.to_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200511/Glucose'+filename+'/QvalueStoreyGlucose_c.xlsx')            
            Qvalue = df_new#pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200511/Insulin/QvalueStoreyInsulin_c.xlsx', header=0,encoding = "ISO-8859-1",index_col=0)
            return(Pval2,Corr2,Qvalue)  
        def temp20200524_InsGlc(Pval,Corr,ResponseIdxDict):        
            Pval2 = Pval[['AUC_Glucose','AUChalf_Glucose','AUC_Insulin','AUChalf_Insulin']].iloc[13:]   #Insulinの2つの特徴量と対応する他の分子x特徴量行列         
            Pval2['AUC_Glucose']['AUC_Glucose']=np.nan;Pval2['AUChalf_Glucose']['AUChalf_Glucose']=np.nan;Pval2['AUC_Glucose']['AUChalf_Glucose']=np.nan;
            Pval2['AUC_Insulin']['AUC_Insulin']=np.nan;Pval2['AUChalf_Insulin']['AUChalf_Insulin']=np.nan;Pval2['AUC_Insulin']['AUChalf_Insulin']=np.nan;

            Pval2['AUC_Glucose']['AUC_Insulin']=np.nan; Pval2['AUChalf_Glucose']['AUChalf_Insulin']=np.nan;
            Pval2['AUC_Glucose']['AUChalf_Insulin']=np.nan; Pval2['AUChalf_Glucose']['AUC_Insulin']=np.nan;

            Corr2 = Corr[['AUC_Glucose','AUChalf_Glucose','AUC_Insulin','AUChalf_Insulin']].iloc[13:];   #Insulinの2つの特徴量と対応する他の分子x特徴量行列         
            Corr2['AUC_Glucose']['AUC_Glucose']=np.nan;Corr2['AUChalf_Glucose']['AUChalf_Glucose']=np.nan; Corr2['AUC_Glucose']['AUChalf_Glucose']=np.nan; 
            Corr2['AUC_Insulin']['AUC_Insulin']=np.nan;Corr2['AUChalf_Insulin']['AUChalf_Insulin']=np.nan;Corr2['AUC_Insulin']['AUChalf_Insulin']=np.nan;
            Corr2['AUC_Glucose']['AUChalf_Insulin']=np.nan; Corr2['AUChalf_Glucose']['AUC_Insulin']=np.nan;
            
                
            SC.UseR(Pval2,{'EngLabel':'GlucoseInsulin'+filename})
            #csvをxlsxへ
            df_new = pd.read_csv("/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/"+ResponseIdxDict['Today']+"/GlucoseInsulin"+filename+"/QvalueStoreyGlucoseInsulin"+filename+"_c.csv", header=0,encoding = "ISO-8859-1",index_col=0)            
            df_new.to_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/'+ResponseIdxDict['Today']+'/GlucoseInsulin'+filename+'/QvalueStoreyGlucoseInsulin'+filename+'_c.xlsx')            
            Qvalue = df_new#pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200511/Insulin/QvalueStoreyInsulin_c.xlsx', header=0,encoding = "ISO-8859-1",index_col=0)
            return(Pval2,Corr2,Qvalue)  
        def temp20200524_InsGlcEAchIdx(Pval,Corr,ResponseIdxDict):        
            Pval2 = Pval[['AUC_Glucose','AUChalf_Glucose','AUC_Insulin','AUChalf_Insulin']].iloc[13:]   #Insulinの2つの特徴量と対応する他の分子x特徴量行列         
            Pval2['AUC_Glucose']['AUC_Glucose']=np.nan;Pval2['AUChalf_Glucose']['AUChalf_Glucose']=np.nan;Pval2['AUC_Glucose']['AUChalf_Glucose']=np.nan;
            Pval2['AUC_Insulin']['AUC_Insulin']=np.nan;Pval2['AUChalf_Insulin']['AUChalf_Insulin']=np.nan;Pval2['AUC_Insulin']['AUChalf_Insulin']=np.nan;

            Pval2['AUC_Glucose']['AUC_Insulin']=np.nan; Pval2['AUChalf_Glucose']['AUChalf_Insulin']=np.nan;
            Pval2['AUC_Glucose']['AUChalf_Insulin']=np.nan; Pval2['AUChalf_Glucose']['AUC_Insulin']=np.nan;

            Pval2['AUC_Glucose'].iloc[13:]= [np.nan]*13; Pval2['AUC_Insulin'].iloc[13:]= [np.nan]*13;  Pval2['AUChalf_Glucose'].iloc[0:13] = [np.nan]*13;  Pval2['AUChalf_Insulin'].iloc[0:13] = [np.nan]*13; 
            
            Corr2 = Corr[['AUC_Glucose','AUChalf_Glucose','AUC_Insulin','AUChalf_Insulin']].iloc[13:];   #Insulinの2つの特徴量と対応する他の分子x特徴量行列         
            Corr2['AUC_Glucose']['AUC_Glucose']=np.nan;Corr2['AUChalf_Glucose']['AUChalf_Glucose']=np.nan; Corr2['AUC_Glucose']['AUChalf_Glucose']=np.nan; 
            Corr2['AUC_Insulin']['AUC_Insulin']=np.nan;Corr2['AUChalf_Insulin']['AUChalf_Insulin']=np.nan;Corr2['AUC_Insulin']['AUChalf_Insulin']=np.nan;
            Corr2['AUC_Glucose']['AUChalf_Insulin']=np.nan; Corr2['AUChalf_Glucose']['AUC_Insulin']=np.nan;
 
            Corr2['AUC_Glucose'].iloc[13:]= [np.nan]*13; Corr2['AUC_Insulin'].iloc[13:]= [np.nan]*13;  Corr2['AUChalf_Glucose'].iloc[0:13] = [np.nan]*13;  Corr2['AUChalf_Insulin'].iloc[0:13] = [np.nan]*13; 
            plt.hist(Corr2.values[~np.isnan(Corr2.values)]);plt.title('_Wpearson'); plt.savefig(save_dir + '_DistOfCorr_Wpearson'+filename+'_GlcInsEach.pdf');plt.close()
            plt.hist(Pval2.values[~np.isnan(Pval2.values)]);plt.title( '_Wpearson'); plt.savefig(save_dir + '_DistOfPval_Wpearson'+filename+'_GlcInsEach.pdf');plt.close()
            
            SC.UseR(Pval2,{'EngLabel':'GlucoseInsulin'+filename})
            #csvをxlsxへ
            df_new = pd.read_csv("/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/"+ResponseIdxDict['Today']+"/GlucoseInsulin"+filename+"/QvalueBHGlucoseInsulin"+filename+"_c.csv", header=0,encoding = "ISO-8859-1",index_col=0)            
            df_new.to_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/'+ResponseIdxDict['Today']+'/GlucoseInsulin'+filename+'/QvalueBHGlucoseInsulin'+filename+'_c.xlsx')            
            Qvalue = df_new#pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200511/Insulin/QvalueStoreyInsulin_c.xlsx', header=0,encoding = "ISO-8859-1",index_col=0)
            return(Pval2,Corr2,Qvalue)  
        def NormalCase(Pval2,Corr2,ResponseIdxDict):
            plt.hist(Corr2.values[~np.isnan(Corr2.values)]);plt.title('_Wpearson'); plt.savefig(save_dir + '_DistOfCorr_Wpearson'+filename+'_GlcInsEach.pdf');plt.close()
            plt.hist(Pval2.values[~np.isnan(Pval2.values)]);plt.title( '_Wpearson'); plt.savefig(save_dir + '_DistOfPval_Wpearson'+filename+'_GlcInsEach.pdf');plt.close()
            
            SC.UseR(Pval2,{'EngLabel':'GlucoseInsulin'+filename})
            #csvをxlsxへ
            df_new = pd.read_csv("/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/"+ResponseIdxDict['Today']+"/GlucoseInsulin"+filename+"/QvalueBHGlucoseInsulin"+filename+"_c.csv", header=0,encoding = "ISO-8859-1",index_col=0)            
            df_new.to_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/'+ResponseIdxDict['Today']+'/GlucoseInsulin'+filename+'/QvalueBHGlucoseInsulin'+filename+'_c.xlsx')            
            Qvalue = df_new#pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200511/Insulin/QvalueStoreyInsulin_c.xlsx', header=0,encoding = "ISO-8859-1",index_col=0)
            return(Pval2,Corr2,Qvalue)  


        #Pval2 = temp20191211()
        #Pval2,Corr2,Qvalue = temp20200511_Ins(Pval,Corr)
        #Pval2,Corr2,Qvalue = temp20200511_Glc(Pval,Corr)
        #Pval2,Corr2,Qvalue = temp20200524_InsGlc(Pval,Corr,ResponseIdxDict)
        #Pval2,Corr2,Qvalue = temp20200524_InsGlcEAchIdx(Pval,Corr,ResponseIdxDict)
        Pval2,Corr2,Qvalue = NormalCase(Pval,Corr,ResponseIdxDict)

        #if 'Prop' not in NameList[j]:
        #CorrUpper = MH.adjustMatrlower(Corr);  PvalUpper = MH.adjustMatrlower(Pval) #対角と下三角をNaNにする
        #Corr.to_excel(save_dir+'Corr_Wspeaman_Upper.xlsx')    ; Pval.to_excel(save_dir+'Pvalue_Wspeaman_Upper.xlsx')        
        #SC.UseR(Pval2,{'EngLabel':'Opeak_Fasting'})
        #Qvalue  = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191211/Opeak_Fasting/QvalueStoreyOpeak_Fasting_c.xlsx', header=0,encoding = "ISO-8859-1",index_col=0)
        Qvalue.columns  = list(Pval2.columns);Qvalue.index  = list(Pval2.index);
        #Qvalue = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181023/Fasting_AUC_Gain/QvalueStoreyParamACDE_ParamACDE.xlsx',header=0,encoding = "ISO-8859-1")
        r = re.compile("(.*)(_)(.*)") 
        CombDF = APPHel.mkExcellower2(Corr2,Pval2,Qvalue,0.1)
        ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(list(CombDF['Mol'])[x]).group(3),'MolColor'] if re.compile("(.*)(_)(.*)").search(list(CombDF['Mol'])[x]).group(3) in list(LabelSum.index) else 'white' for x in range(len(CombDF['Mol']))  ]
        CombDF['Color'] = ColorList
        CombDF = CombDF.sort_values(by='Property')
        CombDF.to_excel(save_dir+'CombDF'+filename+'.xlsx')
        APPHel.PlotScatterDF(RespIdxProp,CombDF,save_dir ,3)#filename)#1で全部赤、3で#散布図：生を赤に、負を青に 'malefemale': 男女色分け 'male':男性のみ
        
        givetoGraph(Corr2,Qvalue,save_dir+'/')    

    
def calcCorrBetwMolRespIdx(SubjectName,ResponseIdxDict, save_dir):#特徴量同士に相関が歩かないか
        RespIdxMolDict={}; NameList=['Opeak', 'Fasing']#'Gain','Fasing', 'Tpeak', 'Ttimeconctant','FC', 'TPI','Adaptation Precision', 'AUC' ]
        if ResponseIdxDict['Target_AnaResponse'] == 'SubjSign':#SubjSign ：各被験者の変動量指標の計算
            tempDF = pd.read_excel(save_dir +'AllRespIdx_'+ SubjectName[0] + '__'+ ResponseIdxDict['AbsResp']+'_Eng.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        else:
            tempDF = pd.read_excel(save_dir +'AllRespIdx_'+ SubjectName[0]  + '__'+ ResponseIdxDict['AbsResp']+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)#.drop('time(min)')
        MolName= list(tempDF.index)
        RespIdxDF = pd.read_excel(save_dir+'/SubjMolRespIdx.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        RespIdxDF_Z = RespIdxDF.apply(zscore, axis=0)#全ての列を標準化
        #同じ分子で同じ列にまとめる
        NewDF = pd.DataFrame(data=None,index=[],columns= NameList)
        r = re.compile("(.*)(_)(.*)") 
        #d = r.search(ii) 
                    #if MolLabelList[jj] == d.group(3):
        for i in MolName:
            #for j in NameList: #特徴量で回す                
                List1 = [list(RespIdxDF.columns)[k]  for k in range(len(list(RespIdxDF.columns))) if i == r.search(list(RespIdxDF.columns)[k]).group(3) ]
                temparray = np.array( RespIdxDF_Z[ List1 ] ); AddDF = pd.DataFrame(data=temparray,index=SubjectName,columns= NameList)
                NewDF=pd.concat([NewDF,AddDF],axis=0, join_axes=[NewDF.columns])#各指標×分子で標準化
                #if i in list(RespIdxDF.columns)
            
        NewDF.to_excel(save_dir+'SubjIdx_Z.xlsx')
        RespIdxDF_Z .to_excel(save_dir+'SubjIdx_Z_All.xlsx')
        NewDF=NewDF.dropna(axis=0)
        Corr,Pval = SC.calculate_CorrWspeaman(NewDF)
        Corr.to_excel(save_dir + 'SubjIdx_Z_Corr_Wpearson1.xlsx'); Pval.to_excel(save_dir + 'SubjIdx_Z_Pvalue_Wpearson1.xlsx');
        CorrUpper = MH.adjustMatrlower(Corr);  PvalUpper = MH.adjustMatrlower(Pval) #対角と下三角をNaNにする
        Corr.to_excel(save_dir+'Corr_Wpearson_Upper.xlsx')    ; Pval.to_excel(save_dir+'Pvalue_Wpearson_Upper.xlsx')        

        plt.hist(Corr.values[~np.isnan(Corr.values)]);plt.title('SubjIdx_Z_Wpearson'); plt.savefig(save_dir + 'SubjIdx_Z_DistOfCorr_Wpearson.pdf');plt.close()
        plt.hist(Pval.values[~np.isnan(Pval.values)]);plt.title( 'SubjIdx_Z_Wpearson'); plt.savefig(save_dir + 'SubjIdx_Z_DistOfPval_Wpearson.pdf');plt.close()
        
        GH.DFScatter(NewDF,[],9,9,10,[],save_dir)#DFで受け取って列の分だけ散布図描画

        ##NewDF.to_excel(save_dir + 'NewDF.xlsx')     
        #Namecount=0

        #RemakeDF = pd.DataFrame(data=None,index=SubjectName,columns=['0'])

        #for jj in NameList: #特徴量で回す
         #   Molcount=0            
          #  for ii in MolName:#各分子で回す  
           #     if ResponseIdxDict['Target_Anal']=='mkHist':#ヒストグラム作る
            #        TargetList = [x for x in list(NewDF.loc[ii][jj]) if x == x and x not in  [float("inf"),float("-inf") ] ]
             #       print(TargetList)                
              #      MSMSH.mkhistList( TargetList,save_dir,NameList[Namecount]+MolName[Molcount],NameList[Namecount]+MolName[Molcount],10)#1つのリストのヒストグラム                
               # else:# DF作り直す
                #    List1=list(NewDF.loc[ii][jj])[:-1]
                 #   RemakeDF[ jj + '_' + ii ] = List1
                #Molcount+=1
            #Namecount+=1
        #RemakeDF.to_excel(save_dir+'SubjMolRespIdx.xlsx')
        
        
def calccorrpval(DF,TargetMolComb):
    
    #for jj in range(len(TargetMolComb)):
    if TargetMolComb['Exp'] == '2Mol':
        TargetComb = TargetMolComb['Target']
        TargetCombLabel = TargetMolComb['Label']
        
        DFNewCorr  = pd.DataFrame(data=None,columns=TargetCombLabel,index=TargetCombLabel);    DFNewPval  = pd.DataFrame(data=None,columns=TargetCombLabel,index=TargetCombLabel)
        DFNewScat = pd.DataFrame(data=None,columns=TargetCombLabel,index=None)
        Col = DF.columns
        for ii in range(len(TargetComb)):#全組み合わせ回す
            Target1 =[]; Target2 =[]
            Target1 += [Col[jj]  for jj in range(len(Col)) if TargetComb[ii][0] == Col[jj][:-2]] ; Target2 += [Col[jj]  for jj in range(len(Col)) if TargetComb[ii][1] == Col[jj][:-2]]
            DF1Normz = DF[Target1].apply(zscore, axis=0); DF2Normz = DF[Target2].apply(zscore, axis=0)
            #list1 = DF1Normz.T.values.tolist(); list2 = DF2Normz.T.values.tolist()
            list1 = list(DF1Normz.values.flatten()); list2 = list(DF2Normz.values.flatten())
            if TargetMolComb['method'] == 'pearson':
                if len(list1) != len(list2):
                    print(TargetComb[ii])
                else:
                    r, p = pearsonr(list1,list2)
            elif TargetMolComb['method'] == 'spearman':
                r, p = spearmanr(list1,list2)
            DFNewCorr.loc[TargetComb[ii][0],TargetComb[ii][1]] = r; DFNewPval.loc[TargetComb[ii][0],TargetComb[ii][1]] = p
            DFNewCorr.loc[TargetComb[ii][1],TargetComb[ii][0]] = r; DFNewPval.loc[TargetComb[ii][1],TargetComb[ii][0]] = p
    
            DFNewScat.loc[:,TargetComb[ii][0]] = list1; DFNewScat.loc[:,TargetComb[ii][1]] = list2;
    elif TargetMolComb['Exp'] == 'MolTimeProp':
        for ii in range(len(TargetComb)):#全組み合わせ回す
            Target1 =[]; Target2 =[]
            Target1 += [Col[jj]  for jj in range(len(Col)) if TargetComb[ii][0] in Col[jj]] ; Target2 += [Col[jj]  for jj in range(len(Col)) if TargetComb[ii][1] in Col[jj]]
            DF1Normz = DF[Target1].apply(zscore, axis=0); DF2Normz = DF[Target2].apply(zscore, axis=0) 
        
        
    return(DFNewCorr,DFNewPval,DFNewScat)
    

def DictUpdate(Dict,Name,List):
    Dict.update({Name:List})
    return(Dict)


def FeatureClustering(SubjectName,ResponseIdxDict, save_dir):#特徴量のクラスタリング
    if ResponseIdxDict['Target_Hist']=='RespIdxAll':#特徴量全部なら
        pass
    elif ResponseIdxDict['Target_Hist']=='RespIdxMol':#特徴量×分子20人文ずつなら   
        RespIdxMolDict={}; NameList=['AUChalf']#'Gain', 'Fasting','Tpeak','Ttimeconctant','FC', 'TPI','Adaptation Precision', 'AUC','AUChalf','Fasting' ]#'Opeak','Gain',
        if ResponseIdxDict['Target_AnaResponse'] == 'SubjSign':#SubjSign ：各被験者の変動量指標の計算
            tempDF = pd.read_excel(save_dir +'AllRespIdx_'+ SubjectName[0] + '__'+ ResponseIdxDict['AbsResp']+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        else:
            tempDF = pd.read_excel(save_dir +'AllRespIdx_'+ SubjectName[0] + '__'+ ResponseIdxDict['AbsResp']+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)#.drop('time(min)')
        NewDF = pd.DataFrame(data=None,index=list(tempDF.index),columns=list(tempDF.columns))
        MolName= list(tempDF.index)
        for i in range(0,len(SubjectName)):
            SubjIdxDF = pd.read_excel(save_dir +'AllRespIdx_'+ SubjectName[i] + '__'+ ResponseIdxDict['AbsResp']+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)#.drop('time(min)')
            #SubjIdxDF = CNE.ChnageNametoEnglish(SubjIdxDF ,0)  
            #SubjIdxDF.to_excel(save_dir+'AllRespIdx_'+ SubjectName[0] + '__'+ ResponseIdxDict['AbsResp']+'_Eng.xlsx')
    
            NewDF = pd.concat([NewDF,SubjIdxDF],axis=1, join_axes=[SubjIdxDF.index])
        NewDF.to_excel(save_dir + 'temp_AllSubjRespIdxDF4Clustering.xlsx') 
        Namecount=0

        if ResponseIdxDict['Gender_Hist'] ==0:#性別で分けない
            for jj in NameList: #特徴量で回す
                RemakeDF = pd.DataFrame(data=None,index=MolName,columns=SubjectName)#特徴量ごとに作る
                Molcount=0            
                for ii in MolName:#各分子で回す  
                    if ResponseIdxDict['Target_Anal']=='mkHist':#ヒストグラム作る
                        TargetList = [x for x in list(NewDF.loc[ii][jj]) if x == x and x not in  [float("inf"),float("-inf") ] ]
                        MSMSH.mkhistList( TargetList, save_dir ,NameList[Namecount]+MolName[Molcount], NameList[Namecount]+MolName[Molcount], 10)#1つのリストのヒストグラム                
                    else:# DF作り直す
                        List1=list(NewDF.loc[ii][jj])[1:]
                        if ResponseIdxDict['Nmlz'] == 'Zscore': #ノーマライズするしない、時間方向zscore：'Zscore'
                            List1 = zscore(List1)#scipyの方ではNaN処理できない？
                        elif ResponseIdxDict['Nmlz'] =='DivAve' :#集団平均で割る
                             List1 = list(np.array(List1)/np.mean(np.array(List1)))
                        RemakeDF.loc[ii] = List1
                    Molcount+=1
                Namecount+=1
                #RemakeDF=RemakeDF.drop('0')
                RemakeDF.to_excel(save_dir+jj+'_SubjMolRespIdx4Clustering'+ResponseIdxDict['Nmlz']+'.xlsx')  
                Optindict={}
                Optindict['title']=jj#'Bolus_MolFeature'
                Optindict['Check'] = ''#ヒートマップ描画時AAのラベルなど変える：'Amino','Glc','AminoGlc': どちらも
                Optindict['AminoCheck'] = ''#'protein','ketogenic','EAA','SemiEAA'
                
                #ClstMol,Colordict,ClstDF,ClstColorDF,ColorDF = mD.draw_heatmap(RemakeDF,1,'MolColor',Optindict,save_dir+jj+'_'+ResponseIdxDict['Nmlz'],cmap='bwr')
                
                RemakeDF.columns = [i+'_'+jj for i in list(RemakeDF.columns)]  
                if Namecount==1:#初めの特徴量なら
                    tempDF = RemakeDF
                else:
                    tempDF = pd.concat([tempDF,RemakeDF],axis=1, join_axes=[tempDF.index])

            ClstMol,Colordict,ClstDF,ClstColorDF,ColorDF = mD.draw_heatmap(tempDF,1,'MolColor',Optindict,save_dir+jj+'_'+ResponseIdxDict['Nmlz'] + '_All',cmap='bwr')
            
    #各Cluster, 分子の時系列
    OptionSwitch={'Target':''}#'EachCluster' : 各Cluster のloading時系列を描画
    labelProp =  tempDF.columns
    MolLabel=tempDF.index
    CompareBolus=dict({}); ClstMolDF=pd.DataFrame(data=None); ClstColorDF=pd.DataFrame(data=None);ClstAveTimeSeriesDF=pd.DataFrame(data=None)
     #PC平面プロっト：'Biplot'　＃楕円フィッティング：'DrawEllipse' #各PCの寄与率：'PlotPCACovRatio'
    #各分子のLoading時系列をプロット、クラスター色分け：'LoadingPlotWCluster'　
    AnalSwitch={'Biplot':0,'DrawEllipse':0,'PlotPCACovRatio':0,'LoadingPlotWCluster':0,
                'ScoreHeatMap':0 ,'FactorLoadingHeatMap':0,#Score, FactorLoadingのヒートマップ：'ScoreHeatMap','FactorLoadingHeatMap'
                'ScoreVariance':0,'LoadingPCsimulation':0,#'ScoreVariance':各分子、被験者間のスコアの分散  'LoadingPCsimulation' : #PC1,2固定して片方動かした時の時系列描画 
                'VectorDiagram' : 0,#Bolus->Continuouのベクトル図}
                'calcinner_outer_product' :0, #各平面でのBC間の内積外積を算出、
                'LengVector': 0}#PC1,2score間のベクトル長を算出する。}
    
    LoadingOpt = 0 #FactorLoadingのヒートマップ：#0：データ（縦）×主成分（横）, #1：主成分（縦）×データ（横）
    BiplotSwitch={'DotColor':'MolColor',#'MolColor'：各分子の色、 'ClsterColor'：クラスターの色、'Black':黒、EachMol：各分子の色、'EachMolCluster':分子x○○の時のCluster色'BC':#PCA時に色分けBCにする、'MolColor_AminoGlc':AAのラベルなど変える
                  'Label':'Annotate',#'Annotate':DotにLabel併記
                  'EachMolColor':[],
                  'Check':'',#'Amino','Glc','AminoGlc': どちらも
                  'AminoCheck' :'' #'protein','ketogenic','EAA','SemiEAA' アミノ酸や糖代謝色変え
                  }#'EachMolColor':分子x被験者PCAにおける、各分子でのいろわえk

    PCAPre.PCAprep(tempDF,ClstAveTimeSeriesDF,ClstColorDF,MolLabel,labelProp,ClstMolDF,AnalSwitch,LoadingOpt,BiplotSwitch,OptionSwitch,CompareBolus,save_dir)

    #ClstAveTimeSeriesDF=mD.plotTimeCourse(ClstMol,Colordict,XDF,EachSwitch,save_dir)#変動指標値、クラスターごと個別+平均値タイムコース描画
    #ClstAveTimeSeriesDF = ACD.mkClstAveTimeSign(save_dir,DF.T,ClstMol)#クラスタごとに分けた平均をだす
    #ClstAveTimeSeriesDF.to_excel(save_dir + 'ClstAveTimeSeriesDF.xlsx')
    #BolusContinuousPCA(DF,ClstAveTimeSeriesDF,ClstColorDF,CompareBolus,ClstMol,save_dir)#BolusとContinuous連結してPCA


def mkHistEachIdx(SubjectName,ResponseIdxDict, save_dir):#各特徴量指標の分布
    SubjPropGender = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjPropWgender_fujita.xlsx',header=0,encoding = "ISO-8859-1")['Gender'])
    if ResponseIdxDict['Target_Hist']=='RespIdxAll':#特徴量全部なら
        RespIdxMolDict={}; NameList=[ 'Opeak', 'Tpeak','Ttimeconctant', 'Gain', 'FC', 'TPI','Adaptation Precision', 'AUC' ]
        RespIdxManDict={}; RespIdxWomanDict={};        
        for i in range(0,len(SubjectName)):
            SubjIdxDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181005/AllRespIdx_' + SubjectName[i] + '__Abs.xlsx',header=0,encoding = "ISO-8859-1").drop('time(min)')
            OpeakList = [x for x in list( SubjIdxDF['Opeak'] ) if x == x]
            TpeakList = [x for x in list( SubjIdxDF['Tpeak'] ) if x == x]
            Ttimeconctant = [x for x in list( SubjIdxDF['Ttimeconctant'] ) if x == x]
            GainList = [x for x in list( SubjIdxDF['Gain'] ) if x == x]
            FCList = [x for x in list( SubjIdxDF['FC'] ) if x == x]   
            TPIList = [x for x in list( SubjIdxDF['TPI'] ) if x == x]
            APList = [x for x in list( SubjIdxDF['Adaptation Precision'] )if x == x]
            AUCList = [x for x in list( SubjIdxDF['AUC'] ) if x == x]
            
            RespIdxMolDict.update({'Opeak':OpeakList})
            RespIdxMolDict.update({'Tpeak':TpeakList})            
            RespIdxMolDict.update({'Ttimeconctant':Ttimeconctant})        
            RespIdxMolDict.update({'Gain':GainList})        
            RespIdxMolDict.update({'FC':FCList})        
            RespIdxMolDict.update({'TPI':TPIList})        
            RespIdxMolDict.update({'Adaptation Precision':APList})        
            RespIdxMolDict.update({'AUC':AUCList})        
            #RespIdxMolDict.update({'Opeak':list( SubjIdxDF['Opeak'] )}) 
            List=[ OpeakList, TpeakList,Ttimeconctant,  GainList,FCList, TPIList,APList, AUCList  ]
            if ResponseIdxDict['Gender_Hist'] ==1:#性別で分ける
                if SubjPropGender[i] == 1:#男性なら
                    for ii in range(len(NameList)):
                       RespIdxManDict = DictUpdate(RespIdxManDict,NameList[ii],List[ii]) #関数化したただ辞書に詰めるだけ
                else:#女性なら
                    for ii in range(len(NameList)):                   
                       RespIdxWomanDict = DictUpdate(RespIdxWomanDict,NameList[ii],List[ii]) #関数化したただ辞書に詰めるだけ
        count=0
        for ii in list( RespIdxMolDict.keys() ): 
                if ResponseIdxDict['Gender_Hist'] ==1:#性別で分ける
                    MSMSH.mkhist2List(RespIdxManDict[ii],RespIdxWomanDict[ii], save_dir,NameList[count],NameList[count],10)#1つのリストのヒストグラム
                    #MSMSH.mkhistList(RespIdxWomanDict[ii],save_dir,NameList[count],NameList[count],10)#1つのリストのヒストグラム                
                else:#男女混ぜる
                    MSMSH.mkhistList(RespIdxMolDict[ii],save_dir,NameList[count],NameList[count],10)#1つのリストのヒストグラム
                count+=1
            
    elif ResponseIdxDict['Target_Hist']=='RespIdxMol':#特徴量×分子20人文ずつなら
        RespIdxMolDict={}; NameList=['FC','Slope120']#'Gain', 'Fasting','Tpeak','Ttimeconctant','FC', 'TPI','Adaptation Precision', 'AUC','Fasting' ]#'Opeak','Gain',
        if ResponseIdxDict['Target_AnaResponse'] == 'SubjSign':#SubjSign ：各被験者の変動量指標の計算
            tempDF = pd.read_excel(save_dir +'AllRespIdx_'+ SubjectName[0] + '__'+ ResponseIdxDict['AbsResp']+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        else:
            tempDF = pd.read_excel(save_dir +'AllRespIdx_'+ SubjectName[0] + '__'+ ResponseIdxDict['AbsResp']+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)#.drop('time(min)')
        NewDF = pd.DataFrame(data=None,index=list(tempDF.index),columns=list(tempDF.columns))
        MolName= list(tempDF.index)
        for i in range(0,len(SubjectName)):
            SubjIdxDF = pd.read_excel(save_dir +'AllRespIdx_'+ SubjectName[i] + '__'+ ResponseIdxDict['AbsResp']+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)#.drop('time(min)')
            #SubjIdxDF = CNE.ChnageNametoEnglish(SubjIdxDF ,0)  
            #SubjIdxDF.to_excel(save_dir+'AllRespIdx_'+ SubjectName[0] + '__'+ ResponseIdxDict['AbsResp']+'_Eng.xlsx')

            NewDF = pd.concat([NewDF,SubjIdxDF],axis=1, join_axes=[SubjIdxDF.index])
        NewDF.to_excel(save_dir + 'temp_AllSubjRespIdxDF.xlsx')     
        Namecount=0
        RemakeDF = pd.DataFrame(data=None,index=SubjectName,columns=['0'])

        if ResponseIdxDict['Gender_Hist'] ==0:#性別で分けない
            for jj in NameList: #特徴量で回す
                Molcount=0            
                for ii in MolName:#各分子で回す  
                    if ResponseIdxDict['Target_Anal']=='mkHist':#ヒストグラム作る
                        TargetList = [x for x in list(NewDF.loc[ii][jj]) if x == x and x not in  [float("inf"),float("-inf") ] ]
                        #print(TargetList)                                    
                        MSMSH.mkhistList( TargetList, save_dir ,NameList[Namecount]+MolName[Molcount], NameList[Namecount]+MolName[Molcount], 10)#1つのリストのヒストグラム                
                    else:# DF作り直す
                        List1=list(NewDF.loc[ii][jj])[1:]
                        RemakeDF[ jj + '_' + ii ] = List1
                    Molcount+=1
                Namecount+=1
            RemakeDF=RemakeDF.drop('0',axis=1)
            RemakeDF.to_excel(save_dir+'SubjMolRespIdx.xlsx')
        else:#性別でわけるなら            
            RemakeDF = pd.read_excel(save_dir +'/SubjMolRespIdx.xlsx')
            RemakeDF.index = SubjPropGender
            
            tDF = pd.DataFrame(data=None,index=list(RemakeDF.columns),columns=['Gender'])
            PvalDF = pd.DataFrame(data=None,index=list(RemakeDF.columns),columns=['Gender'])
            
            for jj in NameList: #特徴量で回す
                Molcount=0            
                for ii in MolName:#各分子で回す  
                    TargetList1 = list( RemakeDF.loc[0][ jj + '_' + ii ] )
                    print(TargetList1)
                    TargetList2 = list( RemakeDF.loc[1][ jj + '_' + ii ] )                    
                    if ResponseIdxDict['Target_Anal']=='mkHist':#ヒストグラム作る                        
                        [x for x in list(NewDF.loc[ii][jj]) if x == x and x not in  [float("inf"),float("-inf") ] ]
                        MSMSH.mkhist2List( TargetList1,TargetList2,save_dir,NameList[Namecount]+MolName[Molcount],NameList[Namecount]+MolName[Molcount],5)#1つのリストのヒストグラム                                       
                        MSMSH.mkhistList( TargetList2,save_dir,NameList[Namecount]+MolName[Molcount],NameList[Namecount]+MolName[Molcount],2)#1つのリストのヒストグラム                                    
                    elif ResponseIdxDict['Target_Anal']=='CalcTtest':#男女でt検定                  
                        t, Pval = sts.ttest_ind(TargetList1,TargetList2)                                               
                        #t,Pval = SC.calculate_CorrWpearson(RespIdxProp)#ピアソン
                        tDF.loc[jj+'_'+ii]['Gender']=t
                        PvalDF.loc[jj+'_'+ii]['Gender']=Pval
                        filename='_Ttest'
                    elif ResponseIdxDict['Target_Anal']=='CalcManU':#男女でマンホイットニー検定                 
                        t, Pval = sts.mannwhitneyu(TargetList1,TargetList2)                                               
                        #t,Pval = SC.calculate_CorrWpearson(RespIdxProp)#ピアソン
                        tDF.loc[jj+'_'+ii]['Gender']=t
                        PvalDF.loc[jj+'_'+ii]['Gender']=Pval
                        filename='_ManU'
                    elif ResponseIdxDict['Target_Anal']=='CalcEffectSize':#男女で効果量                          
                        Pval= SC.cohend(TargetList1,TargetList2)
                        tDF.loc[jj+'_'+ii]['Gender']=0
                        PvalDF.loc[jj+'_'+ii]['Gender']=Pval
                        filename='_Cohend'
                    Molcount+=1
                Namecount+=1 
            tDF.to_excel(save_dir + 'StatValue'+filename+'.xlsx'); PvalDF.to_excel(save_dir + 'Pvalue'+filename+'.xlsx');
                       
            plt.hist(tDF.values[~np.isnan(tDF.values.astype(float))].astype(float) );plt.title('_Wpearsonrev'); plt.savefig(save_dir + '_DistOfCorr_Wpearsonrev'+filename+'.pdf');plt.close()
            plt.hist(PvalDF.values[~np.isnan(PvalDF.values.astype(float))].astype(float) );plt.title( '_Wpearsonrev'); plt.savefig(save_dir + '_DistOfPval_Wpearsonrev.'+filename+'.pdf');plt.close()
          
                    
def AnalResponse(SubjectName,ResponseIdxDict, save_dir):#応答の特徴量を計算
    if ResponseIdxDict['ClusteringWOtherIndex']==1: #ばらつきを他の指標と合わせてClustering
        Clustering(ResponseIdxDict,save_dir)#くっつけてクラスタリング
    if ResponseIdxDict['CalcCorrBetwFeature_CVEach']==1: #特徴量ばらつき間の相関_各分子で  
        calcorrbetwMolFetureEach([],ResponseIdxDict,save_dir)#特徴量ばらつき間の相関_各分子で 
    for i in range(0,len(SubjectName)):   #各被験者で特徴量を算出する
        if ResponseIdxDict['Target_AnaResponse'] == 'SubjSign':#SubjSign ：各被験者の変動量指標の計算
            #Delta??#SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/DivBymeanrawStdDelta_Eng/New/SubjTmCs_'+SubjectName[i] + 'DivBymeanRawStdDelta.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
            SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191127/SubjTmCs_' + SubjectName[i]+'DivBymeanRawStdRaw.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)       

            #2019以前#SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/SubjSign/SubjTmCs_' + SubjectName[i] + 'Delta.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        elif ResponseIdxDict['Target_AnaResponse'] == 'SubjB':#BC比較でB
                #SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw_Eng/SubjTmcs_'+SubjectName[i]+'Raw.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
                SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw_Eng/'+SubjectName[i]+'_all.xlsx',header=0,index_col=0)

        elif ResponseIdxDict['Target_AnaResponse'] == 'SubjC':#BC比較でC
                #SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Continuous/TimeCourse/English/Raw/New/'+SubjectName[i]+'_Continuous2h_Raw.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
                SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Continuous/TimeCourse/English/Raw/New/'+SubjectName[i]+'_all.xlsx',header=0,index_col=0)
        elif ResponseIdxDict['Target_AnaResponse'] == 'Average':#被験者平均で  
### temp
                SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20190823/temp_ketoneafter.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        else:
            SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw_Eng/ケトン検出限界含む/SubjTmCs_' + SubjectName[i] + 'Raw.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)#[['Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid','Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate']]#[['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid','Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid','Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate']]#,'Growth hormone']]
        if len(SubjTmCsDf.index) == 14:#-10minあるなら
            Col = list(SubjTmCsDf.columns)
            for ii in Col:
                SubjTmCsDf[ii]=list(pd.to_numeric(SubjTmCsDf[ii],errors='coerce'))
                SubjTmCsDf.index = SubjTmCsDf['time(min)']
            SubjTmCsDf=SubjTmCsDf.fillna(0)    
            FastingDF = np.nansum(SubjTmCsDf.iloc[[0,1]],axis=0) / 2#空腹値を計算する
            #SubjTmCsDf = SubjTmCsDf.drop(-10,axis=0)
            SubjTmCsDf = SubjTmCsDf[1:]
            SubjTmCsDf.iloc[0,:]  =  FastingDF
            if len(SubjTmCsDf.columns) == 84:#3-HB消えてないなら
                SubjTmCsDf = SubjTmCsDf.drop('Total bile acid',axis=1)
            EngLabel = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0)['English'])
            #SubjTmCsDf.columns = EngLabel
            NewDF = pd.DataFrame(data=None, index=list( SubjTmCsDf.index ), columns = ['time(min)'] + list( SubjTmCsDf.columns )[0:-1])
            NewDF=SubjTmCsDf.copy()
            NewDF['time(min)']=list( SubjTmCsDf.index )
            NewDF.to_excel(save_dir+'SubjTmCs_' + SubjectName[i] + 'Raw.xlsx')

        NewDF=SubjTmCsDf.copy()
        NewDF['time(min)']=list( SubjTmCsDf.index )
        SubjTmCsDf=NewDF.copy()

### 解析対象：
        if any(ResponseIdxDict['SignLabel'])==1: SubjTmCsDf=SubjTmCsDf[ResponseIdxDict['SignLabel']+['time(min)']]
        Col = list(SubjTmCsDf.columns)

        if ResponseIdxDict['Nmlz'] == 'Zscore': #ノーマライズするしない、時間方向zscore：'Zscore'
            #SubjTmCsDf.iloc[:,1:] = scipy.stats.zscore(SubjTmCsDf.iloc[:,1:])#time以外：こっちだとNaN処理できてない
            for kj in range(len(SubjTmCsDf)-1):
                SubjTmCsDf.iloc[:,kj] = zscore(SubjTmCsDf.iloc[:,kj+1])
### temp_20190909
        #SubjTmCsDf = SubjTmCsDf[['Growth hormone','SM-C IGF-1','Free fatty acid','time(min)']      ]  
        Gain,TPI,Fastingmean,FastingDF,Opeak,Tpeak,Ttimeconctant,TPrecision, Diff, TP,AUCPeakList,AUCPeakhalfList,RawDiff,Slope120List= IWHel.CalcPeak(SubjTmCsDf,ResponseIdxDict)#時系列内増加max - 空腹値
        SlopeList=[]

        for ij in Col:
            try:
                Gain_dump,Tc_dump,TcThree = AmT.calcOPeak(list(SubjTmCsDf[ij] - SubjTmCsDf[ij][0]),list(SubjTmCsDf.index))#差分を与える：そしてgain
                SlopeList.append(np.abs(TcThree))#傾きは正負をもつ。CV計算するときに厄介なのでabs
            except:
                print(SubjectName[i]+ ij)
                SlopeList.append(np.nan)
    ###############################   AUC,AUChalf
        AUCSeries = Gain.copy()
        AUChalf = Gain.copy()

        for ii in range(len( Diff.columns)-1): #分子の数だけ
            print(Col[ii])

            list1 = list(RawDiff.iloc[:,ii]); TimeList = list(SubjTmCsDf['time(min)']);
            #print(Diff.columns[ii])
            AUCSeries.iloc[ii] = AmT.CalcAUC(list1,TimeList)
           # if Col[ii] not in ResponseIdxDict['DownBoth']:#減少を伴うものは減少したぶんの最大AUC        
            #        list1 = list(Diff.iloc[:,ii]); TimeList = list(SubjTmCsDf['time(min)']);
### 20200510現在、AUC1/2は絶対値のAUCで計算            
            list1 = list(Diff.iloc[:,ii]); TimeList = list(SubjTmCsDf['time(min)']);
            
            AUC, tstar = AmT.CalcAUChalf(list1,TimeList,ResponseIdxDict,Col[ii])
            #AUCSeries.iloc[ii] = AUC
            AUChalf.iloc[ii] = tstar
        #Opeak = np.max(SubjTmCsDf)
    ###############################  FoldChange
        FC = Opeak / FastingDF
        log2FC = np.log2(Opeak / FastingDF)
        #Sensitivity = Gain / float(file_name[i][0:2])
        Adaptation_Precision = SubjTmCsDf.iloc[-5:].mean() / FastingDF
        CalcRespDF = pd.DataFrame(data=None, index = Col , columns=[#'Fasting',
                                                                    #'Opeak', 
                                                                    #'Gain',
                                                                    #'FC',
                                                                    #'log2FC',
                                                                    #'TPI',
                                                                    #'Adaptation Precision',
                                                                     #'TPrecision',
                                                                    # 'AUC',
                                                                    #'Slope3',
                                                                    #'Tpeak',
                                                                    #'Ttimeconctant',
                                                                    # 'AUChalf'
                                                                     ])
       # CalcRespDF['Opeak'] = Opeak; #CalcRespDF['Tpeak'] = Tpeak;  
        #CalcRespDF['Ttimeconctant'] = Ttimeconctant; CalcRespDF['Gain'] = Gain; 
        #CalcRespDF['FC'] = FC;
        #CalcRespDF['TPI'] = TPI;  #CalacRespDF['Sensitivity'] = Sensitivity
        #CalcRespDF['Adaptation Precision'] = Adaptation_Precision; 
        #CalcRespDF['log2FC'] = log2FC;
        #CalcRespDF['Slope120'] = Slope120List;
        CalcRespDF['AUChalf']= AUChalf
        CalcRespDF['AUC'] = AUCSeries; #CalcRespDF['Fasting'] = FastingDF; 
        #CalcRespDF['TPrecision'] = TPrecision;  #
       # CalcRespDF['AUCPeak'] = AUCPeakList;CalcRespDF['AUCPeakhalf'] =AUCPeakhalfList
        
        #CalcRespDF['AUCDivByFasting'] = AUCSeries / FastingDF
        #CalcRespDF['AUCDivByFasting'] = np.log2(AUCSeries / FastingDF)

        CalcRespDF = CalcRespDF.drop('time(min)',axis=0)
        CalcRespDF.to_excel(save_dir + 'AllRespIdx_'+SubjectName[i] +'_'+ResponseIdxDict['Nmlz']+'_' +ResponseIdxDict['AbsResp']+'.xlsx')
        #ここからPanel化していく
        if i == 0:
            #CalcRespPanel = pd.Panel({SubjectName[i]:pd.DataFrame(data=None,index=list(CalcRespDF.index),columns=list(CalcRespDF.columns))})
            SubjDF = CalcRespDF
            #SubjDF = np.abs(CalcRespDF)
        else:
            tempDF = CalcRespDF         
            #tempDF = np.abs(CalcRespDF)       
            SubjDF=pd.concat([SubjDF,tempDF],axis=0)
            #SubjDF=pd.concat([SubjDF,tempDF],axis=0,join_axes=[SubjDF.columns])         
    CalcRespPanel = pd.DataFrame(data=np.array(SubjDF),index=pd.MultiIndex.from_product([SubjectName, list(tempDF.index)]),columns=list(SubjDF.columns))

    if ResponseIdxDict['calcCV']==1: #特徴量ごとにCV計算+plot
        MolFeatureDF = calcCV(CalcRespPanel,ResponseIdxDict,save_dir)
    if ResponseIdxDict['calcVar']==1: #特徴量ごとにVar計算+plot
        MolFeatureDF = calcVar(CalcRespPanel,ResponseIdxDict,save_dir)
        ResponseIdxDict['denominator'] = 'Fasting'
        #calcVarRatio(MolFeatureDF,ResponseIdxDict,save_dir)
    if ResponseIdxDict['calcAveStd']==1: #特徴量ごとにAve.Std計算+plot
        NewStdDF,NewAveDF = calcAveStd(CalcRespPanel,save_dir)#各DFのAve,Std
    if ResponseIdxDict['CalcCorrAveAve']==1: #特徴量ごとのAveとAve波形の特徴量間の相関
        AveFeatureDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20190911/RespIdx/AllRespIdx_SubjectAll__FastReset_Abs_Raw.xlsx',header=0,encoding = "ISO-8859-1")
        NewStdDF,NewAveDF = calcCorrAveAve(NewAveDF,AveFeatureDF,save_dir)      
    if ResponseIdxDict['calcorrbetwMol']==1: #各特徴量、分子間相関
        calcorrbetwMol(CalcRespPanel,save_dir)    
    if ResponseIdxDict['AveStdFeature']==1: #特徴量間の平均と標準へんさの関係性
        plotAveStd(CalcRespPanel,save_dir)
    if ResponseIdxDict['CalcCorrBetwFeature_Draw']==1: #特徴量間の相関_各分子で
        calcorrbetwIdx(CalcRespPanel,save_dir)#特徴量間の相関_各分子で
    if ResponseIdxDict['CalcCorrBetwFeatureEachMol_Draw']==1: #特徴量間の相関_かく分子ごとに
        calcorrbetwIdxAllMol(CalcRespPanel,save_dir)#特徴量間の相関_かく分子ごとに           
    if ResponseIdxDict['CalcCorrBetwFeature_hist']==1: #特徴量間の分布_各分子で
        calcorrbetwIdx_hist(CalcRespPanel,save_dir)#特徴量間の相関_各分子で集めて分布
    if ResponseIdxDict['CalcCorrBetwFeature_CVEach']==1: #特徴量ばらつき間の相関_各分子で  
        calcorrbetwMolFetureEach(MolFeatureDF,ResponseIdxDict,save_dir)#特徴量ばらつき間の相関_各分子で 
    if ResponseIdxDict['CalcCorrBetwFeature_AveEach']==1: #特徴量平均間の相関_分子間で      
        ResponseIdxDict['std']=NewStdDF
        calcorrbetwMolFetureAveEach(NewAveDF,ResponseIdxDict,save_dir)#特徴量平均間の相関_分子間で 
    if ResponseIdxDict['MolFeaturePCA']==1: #特徴量ばらつきのPCA_CalcCVで求めてから
        #あとで消す 201906くらい？
        #MolFeatureDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20190702/RespIdx/MolFeature_temp.xlsx',header=0,encoding = "ISO-8859-1")
        MolFeaturePCA(MolFeatureDF,ResponseIdxDict,save_dir)

def MolFeaturePCA(XDF,ResponseIdxDict,save_dir):#PCAする
     #PC平面プロっト：'Biplot'　＃楕円フィッティング：'DrawEllipse' #各PCの寄与率：'PlotPCACovRatio'
    #各分子のLoading時系列をプロット、クラスター色分け：'LoadingPlotWCluster'　
    AnalSwitch={'Biplot':1,'DrawEllipse':0,'PlotPCACovRatio':1,'LoadingPlotWCluster':0,
                'ScoreHeatMap':1 ,'FactorLoadingHeatMap':1,#Score, FactorLoadingのヒートマップ：'ScoreHeatMap','FactorLoadingHeatMap'
                'ScoreVariance':0,'LoadingPCsimulation':0,#'ScoreVariance':各分子、被験者間のスコアの分散  'LoadingPCsimulation' : #PC1,2固定して片方動かした時の時系列描画 
                'VectorDiagram' : 0,#Bolus->Continuouのベクトル図}
                'calcinner_outer_product' :0, #各平面でのBC間の内積外積を算出、
                'LengVector': 0}#PC1,2score間のベクトル長を算出する。}
    
    LoadingOpt = 0 #FactorLoadingのヒートマップ：#0：データ（縦）×主成分（横）, #1：主成分（縦）×データ（横）
    BiplotSwitch={'DotColor':'ClstColor_direct',#'MolColor'：各分子の色、 'ClsterColor'：クラスターの色、'Black':黒、EachMol：各分子の色、'EachMolCluster':分子x○○の時のCluster色'BC':#PCA時に色分けBCにする、'MolColor_AminoGlc':AAのラベルなど変える
                  'Label':'Annotate',#'Annotate':DotにLabel併記
                  'EachMolColor':[],
                  'Check':'',#'Amino','Glc','AminoGlc': どちらも
                  'AminoCheck' :'EAA' #'protein','ketogenic','EAA','SemiEAA' アミノ酸や糖代謝色変え
                  }#'EachMolColor':分子x被験者PCAにおける、各分子でのいろわえk
    OptionSwitch={'Target':''}#'EachCluster' : 各Cluster のloading時系列を描画
    CompareBolus=dict({}); ClstMolDF=pd.DataFrame(data=None); ClstColorDF=pd.DataFrame(data=None);ClstAveTimeSeriesDF=pd.DataFrame(data=None)
    ##############################################
    MolLabel = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1").index)
    
    #XDF = XDF.T
    ISSIDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/EachMolInterSubjConnectCorr.xlsx',header=0,encoding = "ISO-8859-1").loc['ISSI']
    IDIDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjMolTimeCVAve_Eng.xlsx',header=0,encoding = "ISO-8859-1").loc['IDI'].drop('3-Hydroxybutyrate')
    DF  = pd.concat([XDF,IDIDF],axis=1, join_axes=[XDF.index]) #分子x特徴量
    DF = pd.concat([DF,ISSIDF],axis=1, join_axes=[DF.index]) #分子x特徴量
    DF.to_excel(save_dir+'/MolFeature_IDI_ISSI.xlsx')
    #DF=DF.drop('Opeak',axis=1);DF=DF.drop('Gain',axis=1);DF=DF.drop('log2FC',axis=1);DF=DF.drop('Slope3',axis=1);DF=DF.drop('Tpeak',axis=1);DF=DF.drop('Ttimeconctant',axis=1);
    
    if ResponseIdxDict['Target_PCA_Variance'] == 'Minus': #対象をばらつきのなさに変換する(1からひく。IPSI以外)

        for i in ['Fasting','Opeak','Gain','AUC','log2FC','Slope3','Tpeak','Ttimeconctant', 'AUChalf','IDI']:
            DF[i] = 1 - DF[i]

    if ResponseIdxDict['Target_PCA'] == 'Selected':#選択した特徴量だけなら
        DF=DF.drop('Opeak',axis=1)
        DF=DF.drop('Gain',axis=1)
    
        DF=DF.drop('log2FC',axis=1)
        DF=DF.drop('Slope3',axis=1)
        DF=DF.drop('Tpeak',axis=1)
        DF=DF.drop('Ttimeconctant',axis=1)
    XDF = DF.drop('hs-CRP')#IDIがないため、
    labelProp =  XDF.columns
    MolLabel=XDF.index

    PCAPre.PCAprep(XDF,ClstAveTimeSeriesDF,ClstColorDF,MolLabel,labelProp,ClstMolDF,AnalSwitch,LoadingOpt,BiplotSwitch,OptionSwitch,CompareBolus,save_dir)

def calcorrbetwMolFetureAveEach(MolFeatureDF,ResponseIdxDict,save_dir):#特徴量ばらつき間の相関_分子間
    MolList = list(MolFeatureDF.index)
    IdxList = list(MolFeatureDF.columns)#特徴量
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0).loc[MolList]
    MolColor = MolColor['MolColor']
    MolFeatureDF_Corr=pd.DataFrame(data=None,index=IdxList,columns=IdxList)   
    MolFeatureDF_Pval=pd.DataFrame(data=None,index=IdxList,columns=IdxList)   
    
    nanIdx=0; nanSubj=0
    DF=MolFeatureDF
    if any(DF.T.isnull()) == 1: #DF内に1つでもnanがあれば、
        nanIdx = list(DF[DF.T.isnull().any()].index)
        nanSubj = list(DF.T[DF.T.isnull().any(axis=1)].index)
        DF = DF.dropna(axis=0)
    for j in range(len(IdxList)-1):
        for k in range(len(IdxList)-1):
            Optiondict={'calcR':'pearson','xlabel':IdxList[j],'ylabel':IdxList[k+1],'Annotate':1,'Label':MolList,'title':str(j)+'_'+str(k)+'_'+IdxList[j]+'_'+IdxList[k+1],'errorbar':1,'x_err':list(ResponseIdxDict['std'][IdxList[j]]),'y_err':list(ResponseIdxDict['std'][IdxList[k+1]])}             
            Corr, Pval = GH.mkScatterWHist(list(DF[IdxList[j]]),list(DF[IdxList[k+1]]),save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム   
            try:
                MolFeatureDF_Corr[ IdxList[j]][IdxList[k+1]] = Corr
                MolFeatureDF_Pval[ IdxList[j]][IdxList[k+1]] = Pval
            except:
                pass
        #pg = sns.pairplot(DF.T)
        #pairplot(DF.T)
        #plt.savefig(save_dir  +i +'_'+ str(nanIdx)+'_'+str(nanSubj)+'.png')
        #plt.savefig(save_dir  +i +'_'+ str(nanIdx)+'_'+str(nanSubj)+ '.pdf')
                plt.close()
    MolFeatureDF_Corr.to_excel(save_dir+'MolFeatureDF_Corr.xlsx')
    MolFeatureDF_Pval.to_excel(save_dir+'MolFeatureDF_Pval.xlsx')

    

def calcorrbetwMolFetureEach(MolFeatureDF,ResponseIdxDict,save_dir):#特徴量ばらつき間の相関_各分子で
    #MolList = list(CalcRespPanel.major_axis)
    #IdxList = list(CalcRespPanel.minor_axis)     

    #for i in MolList:     
    #nanIdx=0; nanSubj=0
    ISSIDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/EachMolInterSubjConnectCorr.xlsx',header=0,encoding = "ISO-8859-1").loc['ISSI']
    IDIDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjMolTimeCVAve_Eng.xlsx',header=0,encoding = "ISO-8859-1").loc['IDI'].drop('3-Hydroxybutyrate')

    #tempで
    ###?# MolFeatureDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20190703/Continuous/MolFeature_temp.xlsx',header=0,encoding = "ISO-8859-1")
    MolFeatureDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20190702/RespIdx/MolFeature_temp.xlsx',header=0,encoding = "ISO-8859-1")


    DF  = pd.concat([MolFeatureDF,IDIDF],axis=1, join_axes=[MolFeatureDF.index]) #分子x特徴量
    DF = pd.concat([DF,ISSIDF],axis=1, join_axes=[DF.index]) #分子x特徴量
    if ResponseIdxDict['Target_PCA_Variance'] == 'Minus': #対象をばらつきのなさに変換する(1からひく。IPSI以外)
        for i in ['Fasting','Opeak','Gain','AUC','log2FC','Slope3','Tpeak','Ttimeconctant', 'AUChalf','IDI']:
            DF[i] = 1 - DF[i]
    if any(DF.T.isnull()) == 1: #DF内に1つでもnanがあれば、
        nanIdx = list(DF[DF.T.isnull().any()].index)
        nanSubj = list(DF.T[DF.T.isnull().any(axis=1)].index)
        DF = DF.dropna(axis=0)
    #pg = sns.pairplot(DF.T)
    pairplot(DF)
    plt.savefig(save_dir +'MolFetureEach.png')
    plt.savefig(save_dir + 'MolFetureEach.pdf')
    plt.close()
        
def mkZscore(CorrFastPropDF,label,labelProp):
    meanNP = np.array(CorrFastPropDF) -np.nanmean(np.array(CorrFastPropDF),axis=0,keepdims=True)
    stdNP = np.nanstd(np.array(CorrFastPropDF),axis=0)
    signDF = pd.DataFrame(index=label,columns=labelProp)

    for i in range(len(labelProp)):
        #if i == 0:
         #   signNP = np.array([0]*len(label))
        #else:    
        signNP = meanNP[:,i]/stdNP[i]
        signDF[labelProp[i]] = signNP
    return(signDF)
    
def Clustering(ResponseIdxDict,save_dir):#クラスタリングする
    metab=['Glucose','Pyruvate','Total bile acid',
                                       'Citrate','Free fatty acid','Total ketone body','Glutamic acid',
                  'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate']
    Insmetab=['Glucose','Pyruvate','Free fatty acid','Total ketone body','Citrulline','Isoleucine','Leucine','Tyrosine']

    #ISSIDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/EachMolInterSubjConnectCorr.xlsx',header=0,encoding = "ISO-8859-1").loc['ISSI']
    #IDIDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjMolTimeCVAve_Eng.xlsx',header=0,encoding = "ISO-8859-1").loc['IDI'].drop('3-Hydroxybutyrate')
    MolFeatureDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20190619/RespIdx/MolFeature_temp.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    MolFeatureDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191217/RespIdx/_MolFeature_temp_CV_Delta.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    MolFeaturemeanDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191217/RespIdx/MolFeature_Ave.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    MolFeatureDF = MolFeatureDF.loc[Insmetab]; MolFeaturemeanDF =MolFeaturemeanDF.loc[Insmetab]
    DF=MolFeatureDF[['FC','Ttimeconctant','TPI']];DF.columns=['FC_Variance','Thalf_Variance','TPI_Variance']
    MolFeaturemeanDF = MolFeaturemeanDF[['FC','Ttimeconctant','TPI']]; MolFeaturemeanDF.columns = ['FC_Mean','Thalf_Mean','TPI_Mean']
    DF = pd.concat([DF, MolFeaturemeanDF],axis=1,join_axes=[DF.index])
    #MolFeatureDF=MolFeatureDF.drop('Opeak',axis=1)
    #MolFeatureDF=MolFeatureDF.drop('Gain',axis=1)

    #MolFeatureDF=MolFeatureDF.drop('log2FC',axis=1)
    #MolFeatureDF=MolFeatureDF.drop('Slope3',axis=1)
    #MolFeatureDF=MolFeatureDF.drop('Tpeak',axis=1)
    #MolFeatureDF=MolFeatureDF.drop('Ttimeconctant',axis=1)
    
    #DF  = pd.concat([MolFeatureDF,IDIDF],axis=1, join_axes=[MolFeatureDF.index]) #分子x特徴量
    #DF = pd.concat([DF,ISSIDF],axis=1, join_axes=[DF.index]) #分子x特徴量
    #DF = DF.drop('hs-CRP')#IDIがないため、
    DF = mkZscore(DF,list(DF.index),list(DF.columns))#正規化する#列方向たてに正規化する=相関行列が返ってくる
    if ResponseIdxDict['Target_PCA_Variance'] == 'Minus': #対象をばらつきのなさに変換する(1からひく。IPSI以外)
        for i in ['Fasting','Opeak','Gain','AUC','log2FC','Slope3','Tpeak','Ttimeconctant', 'AUChalf','IDI']:
            DF[i] = 1 - DF[i]

    if ResponseIdxDict['Target_PCA'] == 'Selected':#選択した特徴量だけなら
        DF=DF.drop('Opeak',axis=1)
        DF=DF.drop('Gain',axis=1)
    
        DF=DF.drop('log2FC',axis=1)
        DF=DF.drop('Slope3',axis=1)
        DF=DF.drop('Tpeak',axis=1)
        DF=DF.drop('Ttimeconctant',axis=1)

    Optindict={}
    Optindict['title']='Bolus_MolFeature'
    Optindict['Check'] = ''#ヒートマップ描画時AAのラベルなど変える：'Amino','Glc','AminoGlc': どちらも
    Optindict['AminoCheck'] = 'EAA'#'protein','ketogenic','EAA','SemiEAA'
    
    ClstMol,Colordict,ClstDF,ClstColorDF,ColorDF = mD.draw_heatmap(DF,1,'MolColor',Optindict,save_dir+'Bolus+Continuous_',cmap='bwr')
    #mHH.draw_heatmap(DF,1,'MolColor',Optindict,save_dir+'Bolus+Continuous_',cmap='bwr')
    #各Cluster, 分子の時系列
    OptionSwitch={'Target':''}#'EachCluster' : 各Cluster のloading時系列を描画
    labelProp =  DF.columns
    MolLabel=DF.index
    CompareBolus=dict({}); ClstMolDF=pd.DataFrame(data=None); ClstColorDF=pd.DataFrame(data=None);ClstAveTimeSeriesDF=pd.DataFrame(data=None)
     #PC平面プロっト：'Biplot'　＃楕円フィッティング：'DrawEllipse' #各PCの寄与率：'PlotPCACovRatio'
    #各分子のLoading時系列をプロット、クラスター色分け：'LoadingPlotWCluster'　
    AnalSwitch={'Biplot':1,'DrawEllipse':0,'PlotPCACovRatio':1,'LoadingPlotWCluster':0,
                'ScoreHeatMap':1 ,'FactorLoadingHeatMap':1,#Score, FactorLoadingのヒートマップ：'ScoreHeatMap','FactorLoadingHeatMap'
                'ScoreVariance':0,'LoadingPCsimulation':0,#'ScoreVariance':各分子、被験者間のスコアの分散  'LoadingPCsimulation' : #PC1,2固定して片方動かした時の時系列描画 
                'VectorDiagram' : 0,#Bolus->Continuouのベクトル図}
                'calcinner_outer_product' :0, #各平面でのBC間の内積外積を算出、
                'LengVector': 0}#PC1,2score間のベクトル長を算出する。}
    
    LoadingOpt = 0 #FactorLoadingのヒートマップ：#0：データ（縦）×主成分（横）, #1：主成分（縦）×データ（横）
    BiplotSwitch={'DotColor':'MolColor',#'MolColor'：各分子の色、 'ClsterColor'：クラスターの色、'Black':黒、EachMol：各分子の色、'EachMolCluster':分子x○○の時のCluster色'BC':#PCA時に色分けBCにする、'MolColor_AminoGlc':AAのラベルなど変える
                  'Label':'Annotate',#'Annotate':DotにLabel併記
                  'EachMolColor':[],
                  'Check':'',#'Amino','Glc','AminoGlc': どちらも
                  'AminoCheck' :'EAA' #'protein','ketogenic','EAA','SemiEAA' アミノ酸や糖代謝色変え
                  }#'EachMolColor':分子x被験者PCAにおける、各分子でのいろわえk

    PCAPre.PCAprep(DF,ClstAveTimeSeriesDF,ClstColorDF,MolLabel,labelProp,ClstMolDF,AnalSwitch,LoadingOpt,BiplotSwitch,OptionSwitch,CompareBolus,save_dir)

    #ClstAveTimeSeriesDF=mD.plotTimeCourse(ClstMol,Colordict,XDF,EachSwitch,save_dir)#変動指標値、クラスターごと個別+平均値タイムコース描画
    #ClstAveTimeSeriesDF = ACD.mkClstAveTimeSign(save_dir,DF.T,ClstMol)#クラスタごとに分けた平均をだす
    #ClstAveTimeSeriesDF.to_excel(save_dir + 'ClstAveTimeSeriesDF.xlsx')
    #BolusContinuousPCA(DF,ClstAveTimeSeriesDF,ClstColorDF,CompareBolus,ClstMol,save_dir)#BolusとContinuous連結してPCA

            
            
def  plotAveStd(CalcRespPanel,save_dir):#各特徴量で平均vs標準偏差
    # sklearn.linear_model.LinearRegression クラスを読み込み
    from sklearn import linear_model
    clf = linear_model.LinearRegression()
    #fig,host=plt.subplots(numFeature,numFeature,figsize=(20,20))
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1")
    MolColorList = MolColor['MolColor']

    MolList = list(CalcRespPanel.mean(level=1,axis=0).index)
    IdxList = list(CalcRespPanel.columns)  
    
    AveList=[];StdList=[]
    MolList=['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid','Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid','Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate']#,'Growth hormone']]

    for j in IdxList:#各指標
        for i in MolList: #各分子   
            AveList.append(CalcRespPanel.mean(level=1,axis=0).loc[i][j])
            StdList.append(CalcRespPanel.std(level=1,axis=0,ddof=0).loc[i][j])
        # 散布図
        X=np.log(np.array(AveList)).reshape(len(AveList),1)
        Y=np.log(np.array(StdList)).reshape(len(StdList),1)
        #################################################### logとらない
        mask = ~np.logical_or(np.isnan(AveList), np.isnan(StdList))
        AveList, StdList = np.array(AveList)[mask], np.array(StdList)[mask]
    
    
        X=np.array(AveList).reshape(len(AveList),1)
        Y=np.array(StdList).reshape(len(StdList),1)
        plt.scatter(AveList,StdList )
        
    
        # 予測モデルを作成
        clf.fit(X, Y)
        # 回帰直線
        #plt.plot(np.exp(X), np.exp(clf.predict(X)),lw=1)
       
        #################################################### logとらない
        plt.plot(X, clf.predict(X),lw=1)
        
        #for ii in range(len(AveList)):
         #   plt.text(AveList[ii],StdList[ii],MolList[ii],size=2)
    
        #ついでに相関係数
        mask = ~np.logical_or(np.isnan(X), np.isnan(Y))
        x, y = X[mask], Y[mask]
        r, p = stats.pearsonr(x, y)   
        
        #plt.yscale('log')  # メイン: y軸をlogスケールで描く
        #plt.xscale('log')
        #################################################### logとらない
    
        plt.yticks(size=20);plt.xticks(size=20)
        plt.title("回帰係数="+ str(clf.coef_)+ '\n切片：'+str(clf.intercept_)+'\nR^2='+str(clf.score(X, Y))+'\nR='+str(r)+'\np='+str(p),size=5)
    
        plt.savefig(save_dir+str(j)+'.png')
        plt.savefig(save_dir+str(j)+'.pdf')
        plt.close()
        
        AveList=[];StdList=[]
        res = Y - clf.predict(X)
        # 残差を図示
        plt.scatter(Y , res)
        for ii in range(len(AveList)):
            plt.text(Y[ii],res[ii],MolList[ii],size=2)        
        plt.xlabel('Predictecd values') 
        plt.ylabel('Residuals')
        plt.hlines(y=0, xmin=np.min(Y), xmax=np.max(Y), color='red',lw=1);plt.yticks(size=20);plt.xticks(size=20)
        plt.savefig(save_dir+str(j)+'_res.png')
        plt.savefig(save_dir+str(j)+'_res.pdf')
        plt.close()
        resValue = np.sqrt(res**2)
        
        #NewDF[j] = resValue
        #List1 = list(NewDF[j]); Title = str(j)+'_res'; xlabel = ''; ylabel = j; Color = MolColorList; xticks = MolList; size=20; Titlesize=20; xsize=10;
        #GH.mkBar(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize,save_dir)
    
        #MolColor[j]=List1
        #List1 = list(MolColor[j].sort_values()); xticks = list(MolColor[j].sort_values().index); Color = list(MolColor.sort_values(j)['MolColor'])
        
        #GH.mkSortedBarWHist(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)

        
def calcorrbetwIdx_hist(CalcRespPanel,save_dir):#特徴量間の相関_各分子で集めて分布
    MolList = list(CalcRespPanel.major_axis)
    IdxList = list(CalcRespPanel.minor_axis)  
    for i in MolList:     
        nanIdx=0; nanSubj=0
        DF = CalcRespPanel.major_xs(i)
        if any(DF.T.isnull()) == 1: #DF内に1つでもnanがあれば、
            nanIdx = list(DF[DF.T.isnull().any()].index)
            nanSubj = list(DF.T[DF.T.isnull().any(axis=1)].index)
            DF = DF.dropna(axis=0)
        #ある分子の相関係数行列を作成->Panel化
        CorrDF = calcCorr(DF)
        if i == 'Glucose':#1爪つめら
            Panel = pd.Panel({i:pd.DataFrame(data=None,index=list(CorrDF.index),columns=list(CorrDF.columns))})
            Panel[i] = CorrDF
            #SumarryDF = CombDF 
            #if ij == 'speaman':
             #   SumarryPDF = CombDF
        else:
            Panel[i] = CorrDF
            #Panel.minor_xs('AUC').loc['FC']
            #Panel.minor_xs(i)
    plothist(Panel,save_dir)
    #g = sns.PairGrid(df, dropna=False)#height=1.6, 
    #g.map_diag(draw_hist)

    #Panel = 

def plothist(Panel,save_dir):
    numFeature = len(Panel.minor_axis)
    Feature = list(Panel.minor_axis)
    #fig,host = plt.subplots(numrow,numcol,figsize=(15,15))
    
    fig,host=plt.subplots(numFeature,numFeature,figsize=(20,20))
    #Feature.reverse()
    count=1
    for i in range(numFeature):#0-8
        for j in range(numFeature):#0-9
            if j > i:
                List = np.array(Panel.minor_xs(Feature[i]).loc[Feature[j]])
                #fig.add_subplot(numFeature,numFeature,count)
                #plt.subplot(numFeature,numFeature,count)
                List = List.astype(np.float32)
                host[i,j].hist(List[~np.isnan(List)])
                host[i,j].set_xlim([-1,1])
                count+=1
                host[i,j].set_title(Feature[i]+' vs '+Feature[j])

            if j-i == 1:#最終行
                host[i,j].set_xlabel(Feature[j])
                #if j == 0:#左端
                host[i,j].set_ylabel(Feature[i]) 
            #for jj in len()
            if i >= j:
                host[i,j].tick_params(labelbottom="off",bottom="off") # x軸の削除
                host[i,j].tick_params(labelleft="off",left="off") # y軸の削除
                host[i,j].set_xticklabels([]) 
                host[i,j].spines["right"].set_color("none")  # 右消し
                host[i,j].spines["left"].set_color("none")   # 左消し
                host[i,j].spines["top"].set_color("none")    # 上消し
                host[i,j].spines["bottom"].set_color("none") # 下消し
               #box("off") #枠線の削除
    fig.tight_layout() #ラベルが重ならないように

    plt.savefig(save_dir +'histFeature.png')
    plt.savefig(save_dir  +'histFeature.pdf')
    plt.close()        
    #
def calcCorr(DF):### ある行と異なる行間の網羅的相関
    Col = list(DF.index)
    NewDF = pd.DataFrame(data=None,index=list(DF.index),columns=list(DF.index))
    for i in range(len(Col)-1):
        for j in range(i+1,len(Col)):
            x=np.array(DF.T[Col[i]]);y=np.array(DF.T[Col[j]])
            mask = ~np.logical_or(np.isnan(x.astype(float)), np.isnan(y.astype(float)))
            x, y = x[mask], y[mask]
            r, p = stats.pearsonr(x, y)     
            NewDF.loc[Col[i],Col[j]]=r
            NewDF.loc[Col[j],Col[i]]=r
    return(NewDF)

def calcorrbetwIdxAllMol(CalcRespPanel,save_dir):#特徴量間の相関_かく分子ごとに
    MolList = list(CalcRespPanel.major_axis)
    IdxList = list(CalcRespPanel.minor_axis)     
    MolFeatureDF_Corr = pd.DataFrame(data=None,index=MolList,columns= [str(IdxList[j]) +'_' +str(IdxList[jj]) for j in range(len(IdxList)-1) for jj in np.arange(j+1,len(IdxList))])##各分子の
    MolFeatureDF_Pval = pd.DataFrame(data=None,index=MolList,columns=[str(IdxList[j]) +'_' +str(IdxList[jj]) for j in range(len(IdxList)-1) for jj in np.arange(j+1,len(IdxList))])##各分子の

    for i in MolList:  #ある分子に関して
        
        nanIdx=0; nanSubj=0
        DF = CalcRespPanel.major_xs(i)##各分子のIdx x subject
        if any(DF.T.isnull()) == 1: #DF内に1つでもnanがあれば、
            nanIdx = list(DF[DF.T.isnull().any()].index)
            nanSubj = list(DF.T[DF.T.isnull().any(axis=1)].index)
            DF = DF.dropna(axis=1)
        for j in range(len(IdxList)-1):#各Indexで回す
            for k in range(len(IdxList)-1):#各Indexで回す
                Optiondict={'calcR':'pearson','xlabel':IdxList[j],'ylabel':IdxList[k+1],'Annotate':0,'Label':MolList,'title':str(j)+'_'+str(k)+'_'+IdxList[j]+'_'+IdxList[k+1]+'_'+i+'_'}
                Corr, Pval = GH.mkScatterWHist(list(DF.T[IdxList[j]]),list(DF.T[IdxList[k+1]]),save_dir,['k']*len(DF.T[IdxList[j]]),Optiondict)#2つのリストの散布図+ヒストグラム   
                try:
                    MolFeatureDF_Corr[ IdxList[j] +'_'+ IdxList[k+1] ][i] = Corr
                    MolFeatureDF_Pval[ IdxList[j] +'_'+ IdxList[k+1] ][i] = Pval
                except:
                    pass
        #pg = sns.pairplot(DF.T)
        #pairplot(DF.T)
        #plt.savefig(save_dir  +i +'_'+ str(nanIdx)+'_'+str(nanSubj)+'.png')
        #plt.savefig(save_dir  +i +'_'+ str(nanIdx)+'_'+str(nanSubj)+ '.pdf')
                plt.close()
    MolFeatureDF_Corr.to_excel(save_dir+'MolFeatureDF_Corr.xlsx')
    MolFeatureDF_Pval.to_excel(save_dir+'MolFeatureDF_Pval.xlsx')
    
def calcorrbetwIdx(CalcRespPanel,save_dir):#特徴量間の相関_各分子で
    MolList = list(CalcRespPanel.major_axis)
    IdxList = list(CalcRespPanel.minor_axis)     

    for i in MolList:     
        nanIdx=0; nanSubj=0
        DF = CalcRespPanel.major_xs(i)
        if any(DF.T.isnull()) == 1: #DF内に1つでもnanがあれば、
            nanIdx = list(DF[DF.T.isnull().any()].index)
            nanSubj = list(DF.T[DF.T.isnull().any(axis=1)].index)
            DF = DF.dropna(axis=1)
        #pg = sns.pairplot(DF.T)
        pairplot(DF.T)
        plt.savefig(save_dir  +i +'_'+ str(nanIdx)+'_'+str(nanSubj)+'.png')
        plt.savefig(save_dir  +i +'_'+ str(nanIdx)+'_'+str(nanSubj)+ '.pdf')
        plt.close()

def calcorrbetwMol(CalcRespPanel,save_dir):#各特徴量、分子間相関     
    MolList = list(CalcRespPanel.major_axis)
    IdxList = list(CalcRespPanel.minor_axis)     

    for i in IdxList:                 
        nanIdx=0; nanSubj=0
        DF = CalcRespPanel.minor_xs(i)
        if  i == 'TPrecision':

            DF =   DF.T[['Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid', 'Citrulline', 'Methionine', 'Isoleucine', 'Leucine', 'Tyrosine', '4-Methyl-2-oxopentanoate', 'Glu + threo-beta-methylaspartate']].T#], 'Growth hormone'] #減少or 両方変動分子

        if any(DF.T.isnull()) == 1: #DF内に1つでもnanがあれば、
            nanIdx = list(DF[DF.T.isnull().any()].index)
            nanSubj = list(DF.T[DF.T.isnull().any(axis=1)].index)
            DF = DF.dropna(axis=1)
        #pg = sns.pairplot(DF.T)
        pairplot(DF.T)
        plt.savefig(save_dir  +i +'_'+ str(nanIdx)+'_'+str(nanSubj)+'.png')
        plt.savefig(save_dir  +i +'_'+ str(nanIdx)+'_'+str(nanSubj)+ '.pdf')
        plt.close()

def draw_hist(x, **kws):
    plt.hist(x[~np.isnan(x)])
    plt.ylim([0,5])

def corr_func(x, y, **kws):
    mask = ~np.logical_or(np.isnan(x), np.isnan(y))
    x, y = x[mask], y[mask]
    r, p = stats.pearsonr(x, y)
    ax = plt.gca()
    Roundp = round(p,100)
    Roundr = round(r,100)
    b = '%.2e'%Roundp
    c = '%.2e'%Roundr
    
    TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
    TenRoundr = c[0:c.find('e')] + '$\it{×10^{-'+str(c[len(c)-1])+'}}$'

    #p='{:e}'.format(p)
                                        #Roundq = round(DF2['QVal'][count],100)
                        #bq = '%.2e'%Roundq
                        #TenRoundq = bq[0:bq.find('e')] + '$\it{×10^{-'+str(bq[len(bq)-1])+'}}$'        
    #TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
    #ax.set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp,fontsize=10)    
    ax.annotate(#"r = {:.3f}".format(r),
            'r = '+str(TenRoundr) +'\n'+    #round(r,3)
            'p = '+str(TenRoundp),               
               xy=(.2, .5), 
               xycoords=ax.transAxes,
               size=16)

def pairplot(df):
    g = sns.PairGrid(df, dropna=False)#height=1.6, 
    g.map_diag(draw_hist)
    g.map_upper(sns.regplot, scatter_kws={"s": 8}, fit_reg=False)#, line_kws={"color":  "r"})
    g.map_lower(corr_func)

def calcCV(CalcRespPanel,ResponseIdxDict,save_dir):#各DFの

    MolList = list(CalcRespPanel.mean(level=1,axis=0).index)
    IdxList = list(CalcRespPanel.columns)      
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0).loc[MolList]
    MolColor = MolColor['MolColor']
    NewDF = pd.DataFrame(data=None, index=MolList,columns=IdxList)
    for i in MolList: #各分子   
        for j in IdxList:#各指標,'
            if  j in  ['Tpeak', 'Ttimeconctant','AUChalf','TPrecision','TPI','log2FC']:
                            NewDF.loc[i,j] = CalcRespPanel.swaplevel().loc[i][j].std(ddof=0) / CalcRespPanel.swaplevel().loc[i][j].mean()
                            #NewDF.loc[i,j] = CalcRespPanel.major_xs(i).loc[j].std(ddof=0) / CalcRespPanel.major_xs(i).loc[j].mean()
###上のCVはtemp
            else:
                if CalcRespPanel.swaplevel().loc[i][j].mean() != CalcRespPanel.swaplevel().loc[i][j].mean():#nanなら
                    NewDF.loc[i,j] = np.nan
                else:
                    NewDF.loc[i,j] =  CalcRespPanel.swaplevel().loc[i][j].std(ddof=0) / np.abs(CalcRespPanel.swaplevel().loc[i][j].mean())

                    #NewDF.loc[i,j] = CalcRespPanel.major_xs(i).loc[j].std(ddof=0) / CalcRespPanel.major_xs(i).loc[j].mean()

    CVbar(NewDF,save_dir)
    NewDF.to_excel(save_dir+ '_MolFeature_temp_CV_'+ResponseIdxDict['Target']+'.xlsx')

    if     ResponseIdxDict['calcCV_RelationBC']==1: #BC間で、特徴量のCVに相関関係があるか
        BDF = pd.read_excel(save_dir+'_MolFeature_temp_CV_SubjB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)            
        CDF = pd.read_excel(save_dir+'_MolFeature_temp_CV_SubjC.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)            

        for i in list(BDF.columns):#各特徴量で回す      
            list1=list(BDF[i]); list2=list(CDF[i]); 
            Optiondict={'calcR':'spearman','xlabel':'Bolus','ylabel':'Continuous','Annotate':1,'Label':MolList,'title':'FeatureCVBvsC_'+i}
            GH.mkScatterWHist(list1,list2,save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム        

    return(NewDF)

def calcVar(CalcRespPanel,ResponseIdxDict,save_dir):#各DFの
    MolList = list(CalcRespPanel.mean(level=1,axis=0).index)
    IdxList = list(CalcRespPanel.columns)    
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0).loc[MolList]
    MolColor = MolColor['MolColor']

    NewDF = pd.DataFrame(data=None, index=MolList,columns=IdxList)
    for i in MolList: #各分子   
        for j in IdxList:#['Ttimeconctant','AUChalf','AUC','Fasting','FC']:#IdxList:#各指標,'
            if  j in  ['Fasting']:# 'Ttimeconctant','AUChalf','FC']:
                            NewDF.loc[i,j] = CalcRespPanel.swaplevel().loc[i][j].std(ddof=0) #/ CalcRespPanel.swaplevel().loc[i][j].mean()

                            #NewDF.loc[i,j] = CalcRespPanel.major_xs(i).loc[j].std(ddof=0) / CalcRespPanel.major_xs(i).loc[j].mean()
            else:
                if CalcRespPanel.swaplevel().loc[i][j].mean() != CalcRespPanel.swaplevel().loc[i][j].mean():#nanならCalcRespPanel.major_xs(i).loc[j].mean() != CalcRespPanel.major_xs(i).loc[j].mean():#nanなら
                    NewDF.loc[i,j] = np.nan
                else:
                    NewDF.loc[i,j] = CalcRespPanel.swaplevel().loc[i][j].std(ddof=0)
                    #NewDF.loc[i,j] = CalcRespPanel.major_xs(i).loc[j].var(ddof=0) #/ CalcRespPanel.major_xs(i).loc[j].mean()
            NewDF.loc[i,'AUC_divbyFasting']=CalcRespPanel.swaplevel().loc[i]['AUC'].std(ddof=0) / CalcRespPanel.swaplevel().loc[i]['Fasting'].std(ddof=0)
            NewDF.loc[i,'AUChalf_divbyFasting']=CalcRespPanel.swaplevel().loc[i]['AUChalf'].std(ddof=0) / CalcRespPanel.swaplevel().loc[i]['Fasting'].std(ddof=0)
            
            
            
    CVbar(NewDF,save_dir)
    NewDF.to_excel(save_dir+ '_MolFeature_temp_variance_Fasting_Var_'+ResponseIdxDict['Target']+'.xlsx')

    if     'calcCV_RelationBC' in ResponseIdxDict:
        if ResponseIdxDict['calcCV_RelationBC']==1: #BC間で、特徴量のCVに相関関係があるか
            BDF = pd.read_excel(save_dir+'_MolFeature_temp_CV_SubjB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)            
            CDF = pd.read_excel(save_dir+'_MolFeature_temp_CV_SubjC.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)            
    
            for i in list(BDF.columns):#各特徴量で回す      
                list1=list(BDF[i]); list2=list(CDF[i]); 
                Optiondict={'calcR':'spearman','xlabel':'Bolus','ylabel':'Continuous','Annotate':1,'Label':MolList,'title':'FeatureCVBvsC_'+i}
                GH.mkScatterWHist(list1,list2,save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム        

    return(NewDF)

def calcVarRatio(MolFeatureDF,ResponseIdxDict,save_dir):#分散の比を計算、
    MolList = list(MolFeatureDF.index)

    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0).loc[MolList]
    MolColor = MolColor['MolColor']
    denominator = ResponseIdxDict['denominator']
    MolFeature =list(MolFeatureDF.columns)
    NewDF = MolFeatureDF.copy()
    for i in MolFeature:
        NewDF.loc[:,i] = MolFeatureDF[i] / MolFeatureDF[denominator]
    CVbar(NewDF,save_dir)        
    NewDF.to_excel(save_dir+ '_MolFeature_temp_variance_'+ResponseIdxDict['Target']+'_'+denominator+'(Rario).xlsx')
        
    
    
def plotAveStdBar(NewStdDF,NewAveDF,save_dir):#平均+-標準偏差の棒グラフ   
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0).loc[list(NewStdDF.index)]
    MolColorList = MolColor['MolColor']

    Col = list(NewStdDF.columns); Idx = list(NewStdDF.index)
    
    for i in Col:
        List1 = list(NewAveDF[i]); Title = i; xlabel = ''; ylabel = i; Color = MolColorList; xticks = Idx; size=20; Titlesize=20; xsize=20;
        StdList = list(NewStdDF[i])

        GH.mkBarAveStd(List1, StdList ,Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize,save_dir)
        #GH.mkSortedBarWHist(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)
        MolColor[i] = list(NewAveDF[i])
        MolColor[i+'std'] = list(NewStdDF[i])
        List1 = list(MolColor[i].sort_values()); xticks = list(MolColor[i].sort_values().index); Color = list(MolColor.sort_values(i)['MolColor'])
        #GH.mkSortedBar(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)#任意のリストを代入して棒グラフを描画する
        StdList = list(MolColor.sort_values(i)[i+'std'])
        GH.mkSortedBarAveStdWHist(List1, StdList, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)
    
def calcCorrAveAve(FeatureAveDF,AveFeatureDF,save_dir):#特徴量のAveとAve波形の特徴量相関
    
    FAveCol = list(FeatureAveDF.columns); AveFCol = list(AveFeatureDF.columns); print(FAveCol == AveFCol);MolLabel = list(FeatureAveDF.index);
    MolColor = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1").loc[MolLabel]['MolColor'])

    for i in FAveCol:#各特徴量で回す              
        list1=list(FeatureAveDF[i]); list2=list(AveFeatureDF[i]); 
        Optiondict={'calcR':'pearson','xlabel':'FeatureAveRage','ylabel':'AverageFeature','Annotate':1,'Label':MolLabel,'title':'FAve_vs_AveF'+i}
        GH.mkScatterWHist(list1,list2,save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム
        Optiondict={'calcR':'spearman','xlabel':'FeatureAveRage','ylabel':'AverageFeature','Annotate':1,'Label':MolLabel,'title':'FAve_vs_AveF'+i}
        GH.mkScatterWHist(list1,list2,save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム        
    

def calcAveStd(CalcRespPanel,save_dir):#各DFのAve,Std
    MolList = list(CalcRespPanel.mean(level=1,axis=0).index)
    IdxList = list(CalcRespPanel.columns)    
    NewAveDF = pd.DataFrame(data=None, index=MolList,columns=IdxList)
    NewStdDF = pd.DataFrame(data=None, index=MolList,columns=IdxList)

    for i in MolList: #各分子   
        for j in IdxList:#各指標
            #############################  tempchange! 修正する
            NewStdDF.loc[i,j] = CalcRespPanel.swaplevel().loc[i][j].std(ddof=0)# -CalcRespPanel.swaplevel().loc[i][j].std(ddof=0)#/ CalcRespPanel.swaplevel().loc[i][j].mean()
            NewAveDF.loc[i,j] = CalcRespPanel.swaplevel().loc[i][j].mean()
            
            #NewStdDF.loc[i,j] = CalcRespPanel.major_xs(i).loc[j].std(ddof=0)- CalcRespPanel.major_xs(i).loc[j].std(ddof=0) 
            #NewAveDF.loc[i,j] = CalcRespPanel.major_xs(i).loc[j].mean()
            
    plotAveStdBar(NewStdDF,NewAveDF,save_dir)

    NewStdDF.to_excel(save_dir+ 'MolFeature_Std.xlsx')
    NewAveDF.to_excel(save_dir+ 'MolFeature_Ave.xlsx')

    """
    Optindict=dict({'title':'MolFearture'});ColorSwitch='MolColor'
    RedDF =  pd.DataFrame(data=np.zeros([len(MolList),len(IdxList)]), index=MolList,columns=IdxList)
    #各分子ごと、最小のみ残して、0入れる
    #for i in MolList:    
        #RedDF.loc[i,NewDF.loc[i].idxmin()] = 1 #各分子内で
    for j in IdxList:
        RedDF.loc[NewDF[j].idxmax(),j] = 1 #各分子内で
    RedDF = RedDF.drop('3-Hydroxybutyrate',axis=0)

    RedDF.index=CNE.RefferNewName()
    mHH.draw_redwhitemap(RedDF.T,1,ColorSwitch,Optindict,save_dir,cmap='bwr')#['Reds','Blues']):(a)    
    """
    return(NewStdDF,NewAveDF)

def CVbar(DF,save_dir):#各分子の各特徴量CVの棒グラフ
    #DF = DF.drop('3-Hydroxybutyrate',axis=0)
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0).loc[list(DF.index)]
    MolColorList = MolColor['MolColor'][list(DF.index)]

    Col = list(DF.columns); Idx = list(DF.index)
    
    for i in Col:#['Ttimeconctant','AUChalf','AUC','Fasting','FC']:#Col:
        nanIdx =  list(DF.T[DF.isnull().any()].index)
        #Idx = list( DF.T[DF.isnull().any()].dropna(axis=1).columns)
        try:
            if i != 'hoge':#昔はSubjectName[i]　２０１９１１１７
                if i in nanIdx:#もしnanある分子があるなら、
                    Idx=list(DF[i].dropna().index)
                    MolColorList = MolColor['MolColor'][Idx]
                    List1 = list(DF[i].dropna())
                else:
                    List1 = list(DF[i]);
                Title = i; xlabel = ''; ylabel = i; Color = MolColorList; xticks = Idx; size=40; Titlesize=40; xsize=20;
                GH.mkBar(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize,save_dir+i+'_CV_')
                GH.mkSortedBarWHist(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, dict(),save_dir+i+'_CV_Whist')
                MolColor[i] = list(DF[i])
                List1 = list(MolColor[i].sort_values().dropna()); xticks = list(MolColor[i].sort_values().dropna().index); Color = list(MolColor.sort_values(i).loc[list(MolColor[i].sort_values().dropna().index)]['MolColor'])
                GH.mkSortedBar(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir+i+'_CV_')#任意のリストを代入して棒グラフを描画する
                GH.mkSortedBarWHist(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, dict(),save_dir+i+'_CV_Whist')
        except:
                pass

def DrawPicture(TimeVar, MolCorr, MolLabelDistDict, save_dir):#ある指標を元に図を描画する
    file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/20180712/'
    MolColorIdx  = pd.read_excel(file_dir+ 'MolColor_Eng.xlsx')
    ClstColorIdx  = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng.xlsx')['ClstColor']
    if MolLabelDistDict['MollCorr'] == 'Corr' or MolLabelDistDict['MollCorr'] == 'Corr_Cluster':#濃度変動を元に
        MolCorrAve = list( MolCorr.loc['20C2Ave'] ); MolCorrAveThresh = 0.62
        MolCorrAveOverThresh = list( MolCorr.loc['20C2Ave'][MolCorr.loc['20C2Ave'] > MolCorrAveThresh].index ) 
        MolCorrAveUnderThresh = list( MolCorr.loc['20C2Ave'][MolCorr.loc['20C2Ave'] <= MolCorrAveThresh].index ) 
        
        TimeVarList = list( TimeVar.loc['ZstdAve'] ); TimeVarThresh = 0.4
        TimeVarOverThresh = list(TimeVar.loc['ZstdAve'][TimeVar.loc['ZstdAve'] > TimeVarThresh].index)
        TimeVarUnderThresh = list(TimeVar.loc['ZstdAve'][TimeVar.loc['ZstdAve'] <= TimeVarThresh].index)
        
        CorrAveOvTimeVarOv = set(MolCorrAveOverThresh) & set(TimeVarOverThresh)#波形同じだが順番変わる（右上）
        CorrAveOvTimeVarUnd = set(MolCorrAveOverThresh) - set(CorrAveOvTimeVarOv)#波形同じかつ順番変わらない（左上）
        CorrAveUndTimeVarOv= set( MolCorrAveUnderThresh) & set(TimeVarOverThresh)#波形違うかつ順番変わる（右下）
        CorrAveUndTimeVarUnd = set( MolCorrAveUnderThresh) - set(CorrAveUndTimeVarOv)#波形違うだが順番変わらない（左下）
        
        Row = max( (len(CorrAveOvTimeVarOv) + len(CorrAveUndTimeVarOv) ), (len( CorrAveOvTimeVarUnd ) + len( CorrAveUndTimeVarUnd )) )
        Col = max( (len(CorrAveOvTimeVarOv) + len(CorrAveOvTimeVarUnd)), ( len(CorrAveUndTimeVarOv)+len(CorrAveUndTimeVarUnd) ) )

        tempList = [CorrAveOvTimeVarOv, CorrAveOvTimeVarUnd, CorrAveUndTimeVarOv, CorrAveUndTimeVarUnd ]
        for ii in range(4):
            fig,ax = plt.subplots(1,1)
            
            #quotient = len(tempList[ii]) // 3 
            q, mod = divmod(len(tempList[ii]), 3)
            
            set1=list( np.arange(0,0.9,0.4) )
            x_list = set1 *q + set1[0:mod] #x軸は0,0.4,0.8が分子数だけ

            y_listPre = [[i]*3 for i in np.arange(0.1,q,0.05)] + [[q]*mod]#y軸は1,1,1 2,2,2 3,3,3...など        
            y_list= [item for sublist in y_listPre for item in sublist]
            print(len(x_list))
            print(x_list)
            print(y_list)
            print(len(y_list))
            print(tempList[ii] ) 
            if MolLabelDistDict['MollCorr'] == 'Corr_Cluster':#濃度変動を元に変動量指標クラスタリングも使って            
                for i,(x,y) in enumerate(zip(x_list,y_list)):
                    MolColorIdx['MolColor'][list(tempList[ii])[i]] = 'white'  if ClstColorIdx[list(tempList[ii])[i]] not in ['red','black','blue','cyan'] else MolColorIdx['MolColor'][list(tempList[ii])[i]]
                    ax.annotate(list(tempList[ii])[i],(x,y),size=8, color =  MolColorIdx['MolColor'][list(tempList[ii])[i]] )
            else:
                for i,(x,y) in enumerate(zip(x_list,y_list)):
                    ax.annotate(list(tempList[ii])[i],(x,y),size=8, color =  MolColorIdx['MolColor'][list(tempList[ii])[i]] )
                
            plt.xticks(color="None")
            plt.yticks(color="None")
            plt.tick_params(length=0)
            ax = plt.gca() # get current axis
            ax.spines["right"].set_color("none")  # 右消し
            ax.spines["left"].set_color("none")   # 左消し
            ax.spines["top"].set_color("none")    # 上消し
            ax.spines["bottom"].set_color("none") # 下消し  
              
            plt.savefig(save_dir + str(ii) + '.pdf')
                
    

def mkCountList4His(MolRepsen2MolCombPre,CorrThresh):
    CountDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/Eng/RCount_0.7945723056980873.xlsx',header=0,encoding = "ISO-8859-1")
    Mol1 = list( MolRepsen2MolCombPre['Mol1'][MolRepsen2MolCombPre['Corr']>CorrThresh] )
    Mol2 = list( MolRepsen2MolCombPre['Mol2'][MolRepsen2MolCombPre['Corr']>CorrThresh] )
    CountList = []
    for i in range(len(Mol1)):
              CountList  += [CountDF.loc[Mol1[i],Mol2[i]] / 20 ] 
    return(CountList)
    
def calcpearson(list1,list2,FileName,xlabel,ylabel, save_dir):
    r, p = pearsonr(list1,list2)
    fig,ax = plt.subplots(1,1)
    ax.tick_params(labelsize=10)
    ax.scatter(list1,list2)
    plot_axis = plt.axis()
    ax.set_xlabel(xlabel,fontsize=10)                
    ax.set_ylabel(ylabel,fontsize=10)                
    Roundp = round(p,7)
    b = '%.2e'%Roundp
    TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
    p='{:e}'.format(p)
                        
    plt.plot([0.78, 1],[0.13618187995504369, 0.13618187995504369], c='black')
    plt.xlim([0.78 , 1])
    plt.ylim([-0.1,1.1])    
    TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
    ax.set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp,fontsize=10)
    plt.savefig(save_dir +  FileName + 'Wline.pdf',format='pdf')

def calcpearsonWcountColor(list1,list2,FileName,xlabel,ylabel, save_dir):

    r, p = pearsonr(list1,list2)
    fig,ax = plt.subplots(1,1)
    ax.tick_params(labelsize=10)
    ax.scatter(list1,list2)
    plot_axis = plt.axis()
    ax.set_xlabel(xlabel,fontsize=10)                
    ax.set_ylabel(ylabel,fontsize=10)                
    Roundp = round(p,7)
    b = '%.2e'%Roundp
    TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
    p='{:e}'.format(p)
                        
    #plt.plot([0.78, 1],[0.13618187995504369, 0.13618187995504369], c='black')
    plt.xlim([0.78 , 1])
    plt.ylim([-0.1,1.1])    
    TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
    ax.set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp,fontsize=10)
    RStdCountDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/RmeanStdCountratio.xlsx',header=0,encoding = "ISO-8859-1") 
    
    Corr = list(RStdCountDF[RStdCountDF['Count_ratio'] >= 0.95]['Corr'])
    Std = list(RStdCountDF[RStdCountDF['Count_ratio'] >= 0.95]['Std'])
    ax.scatter(Corr, Std, c='r')
    
    plt.savefig(save_dir +  FileName + 'WColor.pdf',format='pdf')

def mkAveList4His(MolRepsen2MolCombPre,CorrThresh):
    MolRepsen2MolCombAve = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Rmean_Eng.xlsx',header=0,encoding = "ISO-8859-1")
    Mol1 = list( MolRepsen2MolCombPre['Mol1'][MolRepsen2MolCombPre['Corr']>CorrThresh] )
    Mol2 = list( MolRepsen2MolCombPre['Mol2'][MolRepsen2MolCombPre['Corr']>CorrThresh] )
    AveList = []
    for i in range(len(Mol1)):

              AveList  += [MolRepsen2MolCombAve.loc[Mol1[i],Mol2[i]]]
    return(AveList)

def mkStdList4His(MolRepsen2MolCombPre,CorrThresh):
    MolRepsen2MolCombStd = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Rstd_Eng.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    Mol1 = list( MolRepsen2MolCombPre['Mol1'][MolRepsen2MolCombPre['Corr']>CorrThresh] )
    Mol2 = list( MolRepsen2MolCombPre['Mol2'][MolRepsen2MolCombPre['Corr']>CorrThresh] )
    StdList = []
    for i in range(len(Mol1)):

              StdList  += [MolRepsen2MolCombStd.loc[Mol1[i],Mol2[i]]]
    return(StdList)
    
    
def AdjustCVTmCs(SubjAveCorrDict,SubjectName,save_dir):#CVノ時系列吐く
    RawCV = SubjAveCorrDict['Raw'];DeltaCV=SubjAveCorrDict['Delta'];FCCV = SubjAveCorrDict['FC']
    AllCV=pd.concat([RawCV,DeltaCV,FCCV],axis=0,join_axes=[RawCV.columns]) 

    AllCVMulIdxDF = pd.DataFrame(data=np.array(AllCV),index=pd.MultiIndex.from_product([ ['Raw','Delta','FC'],list(RawCV.index)]),columns=list(RawCV.columns))

###  ばらつきの時系列を描画
    DrawTmCsDict={'Draw':'NormalizeFirst_Single',#'NormalizeFirst_Single': #1つの枠に、3条件同時描画なら,'NormalizeFirst_Triple':1つの枠に、1条件x3
                  'BothDelta' : 'RawDiffPctC'}#差分を描画する'RawDeltaFC':#3種類を描画する}　'RawDiffPctC':#3種類を描画する 前後差分、pct変化
    Label = list(AllCVMulIdxDF.columns)

            #Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            #Subjmean_Z = Subjmean.apply(zscore, axis=0)
            #SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean 

    if any(SubjAveCorrDict['SubjTmCsPanel']) ==1:
        if DrawTmCsDict['BothDelta']=='RawDiffPctC':#前点からの差分ならそのばらつきをここで計算する
            DeltaCV=SubjAveCorrDict['SubjTmCsPanel'].groupby(level=0).diff().std(level=1)
            FCCV = SubjAveCorrDict['SubjTmCsPanel'].groupby(level=0).pct_change().std(level=1)
            AllCV=pd.concat([RawCV,DeltaCV,FCCV],axis=0,join_axes=[RawCV.columns]) 
            AllCVMulIdxDF = pd.DataFrame(data=np.array(AllCV),index=pd.MultiIndex.from_product([ ['Raw','Delta','FC'],list(RawCV.index)]),columns=list(RawCV.columns))

        DrawTmCsDict['Draw'] = 'EachSubject_Variance_Triple'
        CWH.drawTmCsMultiIdxDF(SubjAveCorrDict['SubjTmCsPanel'], AllCVMulIdxDF ,Label,DrawTmCsDict,SubjectName,save_dir)#MultiIdx形式のSubjectxTimex項目の時系列を描画する。
    else:
        CWH.drawTmCsMultiIdxDF(AllCVMulIdxDF ,AllCVMulIdxDF , Label,DrawTmCsDict,SubjectName,save_dir)#MultiIdx形式のSubjectxTimex項目の時系列を描画する。
    GH.plotAveStdMultiIndex(AllCVMulIdxDF,save_dir)#各量で平均vs標準偏差_MultiIndexver.

### それぞれの差分のクラスタリング
    RDCV = RawCV - DeltaCV; RFCV = RawCV-FCCV; DFCV = DeltaCV - FCCV
    #RawDeltaFCDiffCV = pd.concat([RDCV,RFCV,DFCV],axis=0,join_axes=[RDCV.columns])
    RawDeltaFCDiffCV = pd.concat([RDCV,RFCV],axis=0,join_axes=[RDCV.columns])

    #RawDeltaFCDiffCV = np.log10(RawDeltaFCDiffCV+np.abs(np.min(np.min(RawDeltaFCDiffCV)))+1)
    #RawDeltaFCDiffCV = np.log10(RawDeltaFCDiffCV+np.abs(np.min(np.min(RawDeltaFCDiffCV)))+1)
    RawDeltaFCDiffCV = RawDeltaFCDiffCV.drop('Growth hormone',axis=1)

    RawDeltaFCDiffCV.to_excel(save_dir+'RawDeltaFCDiffCV.xlsx')

    #ClstMol,Colordict,ClstDF,ClstColorDF = mD.draw_heatmap(RawDeltaFCDiffCV.T,1,'MolColor',{'title':'DifferenceCV'},save_dir,cmap='bwr')#'bwr' ('YlOrRd'),'Reds'
    
    
def calcTmCsAve(OptionDict,SubjectName,save_dir): ## IDIHelperと対応
    print(SubjectName)
    MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)
    MolColor = list(MolColorDF['MolColor'])
    axislabel=2 #0で時間方向、1で分子方向、2で何もしない,3で0点削ってnp.log10  4でRawFCなら全部でで標準化　、7で各時系列の差分をつなげる,

    if axislabel==0:
        axislabelname='Time'
    elif axislabel==1:
        axislabelname='Mol'
    elif axislabel==3:      
        axislabelname='log10'        
    elif axislabel==4:      
        axislabelname='AllTime'   
    elif axislabel==5:      
        axislabelname='AllTime_Wozero'
    elif axislabel ==6:#6で条件の時点間で正規化した時系列をつなげる、
        axislabelname='EachTimeBetwCond'
    elif   axislabel ==7:#7で各時系列の差分をつなげる
        axislabelname='DefferenceBetwCond'
    elif   axislabel ==8:#7で各時系列の差分のつなげて分子内標準化
         axislabelname='DefferenceBetwCondNormalize'
    elif   axislabel ==9:#9で各時系列の差分をつなげる2つだけ       
         axislabelname='DefferenceBetwCond2'
    else:
        axislabelname='No'      
    cmap='Reds'

    if OptionDict['DataFlag'] == 'Raw':
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw_Eng/ケトン検出限界含む/'

        for i in range(len(list(SubjectName))):
            SubjTmCsDelta = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)       
            if i== 0: 
                kanpachi = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
            else:
                kanpachi=pd.concat([kanpachi,SubjTmCsDelta],axis=0,join_axes=[kanpachi.columns])

            #SubjTmCsDeltaPanel[SubjectName[i]] = SubjTmCsDelta
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))

        Subjmean = SubjTmCsDeltaPanel.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)
        SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) / Subjmean

    elif OptionDict['DataFlag'] == 'Delta':
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Delta_Eng/ケトン検出限界含む/'
        #SubjTmCsDeltaPanel = pd.Panel({SubjectName[0]:pd.DataFrame(data=kanpachi,index=list(kanpachi.index),columns=list(kanpachi.columns))})
        for i in range(len(list(SubjectName))):
            SubjTmCsDelta = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Delta.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)       
            #SubjTmCsDeltaPanel[SubjectName[i]] = SubjTmCsDelta
            if i== 0: 
                kanpachi = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Delta.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
            else:
                kanpachi=pd.concat([kanpachi,SubjTmCsDelta],axis=0,join_axes=[kanpachi.columns])            
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)
### 20200114_Deltaのstd or var or CV
        SubjCV = SubjTmCsDeltaPanel.var(level=1,ddof=0) #/ Subjmean
        
    elif  OptionDict['DataFlag'] == 'FC':        
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw_Eng/ケトン検出限界含む/'
        i=0
        #SubjTmCsDeltaPanel = pd.Panel({SubjectName[0]:pd.DataFrame(data=kanpachi,index=list(kanpachi.index),columns=list(kanpachi.columns))})
        for i in range(len(list(SubjectName))):
            SubjTmCsDelta = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)       
            #SubjTmCsDeltaPanel[SubjectName[i]] = SubjTmCsDelta / SubjTmCsDelta.loc[0]
            if i== 0: 
                kanpachi = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                kanpachi = kanpachi / kanpachi.loc[0]
            else:
                kanpachi=pd.concat([kanpachi ,(SubjTmCsDelta/ SubjTmCsDelta.loc[0])],axis=0,join_axes=[kanpachi.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)
### 20200114_FCのstd or var or CV
        SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean   #FCは正規化してあるので平均で割らない？
    elif  OptionDict['DataFlag'] == 'log2FC':        
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw_Eng/ケトン検出限界含む/'
        i=0
        #SubjTmCsDeltaPanel = pd.Panel({SubjectName[0]:pd.DataFrame(data=kanpachi,index=list(kanpachi.index),columns=list(kanpachi.columns))})
        for i in range(len(list(SubjectName))):
            SubjTmCsDelta = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)       
            #SubjTmCsDeltaPanel[SubjectName[i]] = SubjTmCsDelta / SubjTmCsDelta.loc[0]
            if i== 0: 
                kanpachi = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                kanpachi = np.log2(kanpachi / kanpachi.loc[0])
            else:
                kanpachi=pd.concat([kanpachi ,np.log2(SubjTmCsDelta/ SubjTmCsDelta.loc[0])],axis=0,join_axes=[kanpachi.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)
### 20200114_log2FCのstd or var or CV
        SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean   #FCは正規化してあるので平均で割らない？


    elif OptionDict['DataFlag'] == 'RawDeltaFC':#3つの分散をつなげて時系列描画
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw_Eng/ケトン検出限界含む/'

        for i in range(len(list(SubjectName))):
            SubjTmCsDelta = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)       
            if i== 0: 
                kanpachi = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
            else:
                kanpachi=pd.concat([kanpachi,SubjTmCsDelta],axis=0,join_axes=[kanpachi.columns])

            #SubjTmCsDeltaPanel[SubjectName[i]] = SubjTmCsDelta
        SubjTmCsDeltaPanelRaw = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))

        RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191009/Raw_NormalizationFirst_No/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191009/Delta_NormalizationFirst_No/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        FCSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191009/FC_NormalizationFirst_No/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        kanpachi=pd.concat([RawSubjVar,DeltaSubjVar,FCSubjVar],axis=0,join_axes=[RawSubjVar.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(RawSubjVar.index)]),columns=list(RawSubjVar.columns))
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)#ナンデモイイ
        Subjmean_Z = Subjmean.apply(zscore, axis=0)#nanndemoii
        SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) / Subjmean   #何でもいい
        cmap='bwr'
        DrawTmCsDict={'Draw': 'NormalizeFirst_Single_double',#'_Triple':3種類を描画する, 'EachSubject_Triple': #1つの条件で、被験者個別なら　'NormalizeFirst_Single': #1つの枠に、3条件同時描画なら 'NormalizeFirst_Single_double':RawとDeltaのみ
                      'BothDelta' : 'RawDeltaFCWAveStd'}#差分を描画する'RawDeltaFC':#3種類を描画する}
    elif OptionDict['DataFlag'] == 'RawDeltaCVDelta':#RawとDeltaのCVタイムコースの差分
            #SubjTmCsDeltaPanel[SubjectName[i]] = SubjTmCsDelta

        RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191009/Raw_NormalizationFirst_No/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191009/Delta_NormalizationFirst_No/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        FCSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191009/FC_NormalizationFirst_No/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        if any(OptionDict['SignLabel'])==1:RawSubjVar=RawSubjVar[OptionDict['SignLabel']];DeltaSubjVar=DeltaSubjVar[OptionDict['SignLabel']]
        SubjCV=RawSubjVar-DeltaSubjVar;#SubjCV=SubjCV.apply(zscore)
        cmap='bwr'
        ClstMol,Colordict,ClstDF,ClstColorDF = mD.draw_heatmap(SubjCV.T,1,'MolColor',{'title':'CV'+OptionDict['DataFlag']},save_dir+OptionDict['DataFlag']+'/',cmap='bwr')#'bwr' ('YlOrRd'),'Reds'


            
        #Label=list(SubjTmCsDeltaPanel.columns)
        #CWH.drawTmCsMultiIdxDF(SubjTmCsDeltaPanel,SubjTmCsDeltaPanel, Label,DrawTmCsDict,SubjectName,save_dir+OptionDict['DataFlag']+'/')#MultiIdx形式のSubjectxTimex項目の時系列を描画する。
    elif OptionDict['DataFlag'] == 'RawDeltaCVCVDelta':#RawとDeltaのCVCVのDelta
        RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191006/Raw_NormalizationFirst_1/SubjCVCV.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191006/Delta_NormalizationFirst_1/SubjCVCV.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)


        ###Raw-Delta BarGraph
        DeltaRD = RawSubjVar - DeltaSubjVar;    DeltaRDSort = DeltaRD.copy()
        DeltaRDSort['MolColor'] = MolColor ;DeltaRDSort=DeltaRDSort.sort_values(by=0)
        if any(OptionDict['SignLabel'])==1:DeltaRDSort=DeltaRDSort.T[OptionDict['SignLabel']].T.sort_values(by=0)

        GH.mkSortedBarWHist(list(DeltaRDSort[0]), 'R-D', 'Molecule', 'R-D', list(DeltaRDSort['MolColor']), list(DeltaRDSort.index), 10, 15,10, dict(), save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
        Optiondict={'calcR':'spearman',#'spearman'
        'xlabel':'Raw',
            'ylabel':'Delta',
            'Annotate':1,
            'Label':list(RawSubjVar.index),
            'title':'',
            'y=x' : 1#y=xの線を足す

            }
        #if any(BolusDeltaRDSort['SignLabel'])==1:ContinuousDeltaRD=ContinuousDeltaRD.T[OptionDict['SignLabel']].T;BolusDeltaRD=BolusDeltaRD.T[OptionDict['SignLabel']].T
        GH.mkScatterWHist(list(RawSubjVar[0]),list(DeltaSubjVar[0]),save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム


    elif OptionDict['DataFlag'] == 'RawDeltaFC_var':#3つの分散をつなげて時系列描画
        RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Raw/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Delta_all/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        FCSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/FC/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        kanpachi=pd.concat([RawSubjVar,DeltaSubjVar,FCSubjVar],axis=0,join_axes=[RawSubjVar.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(RawSubjVar.index)]),columns=list(RawSubjVar.columns))
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)#ナンデモイイ
        Subjmean_Z = Subjmean.apply(zscore, axis=0)#nanndemoii
        SubjCV = DeltaSubjVar / RawSubjVar
        #SubjTmCsDeltaPanel.std(level=1,ddof=0) / Subjmean   #何でもいい
### 20200114_Rawばらつき vs Deltaばらつき
        SubjMean_Delta=DeltaSubjVar.mean();SubjMean_Raw=RawSubjVar.mean()
        SubjMean_Delta=SubjMean_Delta[['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
                       'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
                       'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate','Growth hormone']]
        SubjMean_Raw=SubjMean_Raw[['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
                       'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
                       'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate','Growth hormone']]
    ### 20200114_ばらつき時系列の時間平均棒グラフ
        OptionDict['save_dir'] =save_dir
        SubjCVMean = SubjMean_Delta / SubjMean_Raw
        GH.mkSortedBarWHist(list((SubjCVMean)), 'VarTimeMean', 'Molecular Name', 'VarTimeMean', LH.MolColor(list(SubjCVMean.index)),list(SubjCVMean.index), 15, 15,15, OptionDict, OptionDict['save_dir'])#ヒストグラム月の任意のリストを代入して棒グラフを描画する
        SubjCVMeanSort = SubjCVMean.sort_values(0)
        GH.mkSortedBarWHist(list((SubjCVMeanSort)), 'VarTimeMean_Sort', 'Molecular Name', 'VarTimeMean', LH.MolColor(list(SubjCVMeanSort.index)),list(SubjCVMeanSort.index), 15, 15,15, OptionDict, OptionDict['save_dir'])#ヒストグラム月の任意のリストを代入して棒グラフを描画する

        cmap='bwr'                 
    elif OptionDict['DataFlag'] == 'RawDeltaFC_CV':#3つのCVをつなげて時系列描画
        RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Raw/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Delta_all/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        #DeltaSubjVar = pd.read_excel('RawDeltaFC_NormalizationFirst_Time_19_zeroあり
        #DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Delta_all/SubjCV_Square_root.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        DeltaSubjVar.loc[0]=0
        FCSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/FC/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        kanpachi=pd.concat([RawSubjVar,DeltaSubjVar,FCSubjVar],axis=0,join_axes=[RawSubjVar.columns])
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(RawSubjVar.index)]),columns=list(RawSubjVar.columns))
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)#ナンデモイイ
        Subjmean_Z = Subjmean.apply(zscore, axis=0)#nanndemoii
        SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) / Subjmean   #何でもいい
        cmap='bwr' 
    elif OptionDict['DataFlag'] == 'Rawvar':#Rawの分散の時間ばらつき
        SubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Raw/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Raw/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Delta_all/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        FCSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/FC/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        kanpachi=pd.concat([RawSubjVar,DeltaSubjVar,FCSubjVar],axis=0,join_axes=[RawSubjVar.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(RawSubjVar.index)]),columns=list(RawSubjVar.columns))
        cmap='bwr'                 
    elif OptionDict['DataFlag'] == 'Deltavar':#Rawの分散の時間ばらつき
        SubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Delta_all/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Raw/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Delta_all/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        FCSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/FC/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        kanpachi=pd.concat([RawSubjVar,DeltaSubjVar,FCSubjVar],axis=0,join_axes=[RawSubjVar.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(RawSubjVar.index)]),columns=list(RawSubjVar.columns))
        cmap='bwr' 
    elif OptionDict['DataFlag'] == 'FCvar':#Rawの分散の時間ばらつき
        SubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/FC/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Raw/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Delta_all/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        FCSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/FC/temp_var.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
        kanpachi=pd.concat([RawSubjVar,DeltaSubjVar,FCSubjVar],axis=0,join_axes=[RawSubjVar.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(RawSubjVar.index)]),columns=list(RawSubjVar.columns))
        cmap='bwr' 
    try:
        Subjmean =  SubjVar#ナンデモイイ
        Subjmean_Z =  SubjVar#nanndemoii
        SubjCV =  SubjVar  
    except:
            pass
#はじめに時系列を平均で割って、正規化する
    if OptionDict['NormalizeFist']==1:
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw_Eng/ケトン検出限界含む/'
        #i=0   
        #kanpachi = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=None) 
        #SubjTmCsDeltaPaneltemp = pd.Panel({SubjectName[0]:pd.DataFrame(data=kanpachi,index=list(kanpachi.index),columns=list(kanpachi.columns))})
        #SubjTmCsDeltaPanel = pd.Panel({SubjectName[0]:pd.DataFrame(data=kanpachi,index=list(kanpachi.index),columns=list(kanpachi.columns))})
        for i in range(len(list(SubjectName))):
            SubjTmCsDelta = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)       
            #SubjTmCsDeltaPaneltemp[SubjectName[i]] = SubjTmCsDelta
            if i== 0: 
                kanpachi = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
            else:
                kanpachi=pd.concat([kanpachi,SubjTmCsDelta],axis=0,join_axes=[kanpachi.columns]) 
        SubjTmCsDeltaPaneltemp = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))

        PanelOptionDict={'Data':'Raw'#Raw,Delta,FC ##Bolus, ContinuousのPanel作るときの正規化
            };#DrawTmCsDict['Data']=PanelOptionDict['Data']
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/'
        if OptionDict['CondFlag'] == 'Bolus':#BolusContinuousならここでMultIndexDF受け取る
            SubjectName=['isaki','iwashi','karei','shimaaji','unagi']
            SubjTmCsDeltaPaneltemp = CWH.mkMultBolus(SubjectName,PanelOptionDict)
        if OptionDict['CondFlag'] == 'Continuous':#BolusContinuousならここでMultIndexDF受け取る
            SubjectName=['isaki','iwashi','karei','shimaaji','unagi']
            SubjTmCsDeltaPaneltemp = CWH.mkMultContinuous(SubjectName,file_dir,PanelOptionDict)

####################### ここでRawを元に正規化する
        for i in range(len(list(SubjectName))):
            TempDF = SubjTmCsDeltaPaneltemp.loc[SubjectName[i]] / SubjTmCsDeltaPaneltemp.mean(level=1)      
            if OptionDict['DataFlag'] == 'Delta':
                TempDF = TempDF - TempDF.loc[0]            
            if (OptionDict['DataFlag'] == 'FC')  or (OptionDict['DataFlag'] == 'RawFC'):
                TempDF = TempDF / TempDF.loc[0]             
            if i== 0: 
                kanpachi = TempDF
            else:
                kanpachi=pd.concat([kanpachi,TempDF],axis=0,join_axes=[kanpachi.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))
################### Bolus Continuous      
        if (OptionDict['CondFlag'] == 'Bolus') or ( OptionDict['CondFlag'] == 'Continuous'):
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean   
            if OptionDict['DataFlag'] == 'RawDeltaFC':#3しゅ合体multDF
                FCSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/'+OptionDict['CondFlag']+'_FC_NormalizationFirst_AllTime_Wozero/SubjCVzero.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/'+OptionDict['CondFlag']+'_Raw_NormalizationFirst_AllTime_Wozero/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/'+OptionDict['CondFlag']+'_Delta_NormalizationFirst_AllTime_Wozero/SubjCVzero.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                kanpachi=pd.concat([RawSubjVar,DeltaSubjVar,FCSubjVar],axis=0,join_axes=[RawSubjVar.columns]) 
                SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(RawSubjVar.index)]),columns=list(RawSubjVar.columns))
                SubjCV = kanpachi    ;cmap='Reds'
            elif OptionDict['DataFlag'] == 'RawDeltaDelta':#RawとDeltaのDelta
                RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Bolus_Raw_NormalizationFirst_AllTime_Wozero/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Bolus_Delta_NormalizationFirst_AllTime_Wozero/SubjCVzero.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                BolusDeltaRD = RawSubjVar - DeltaSubjVar

                RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Continuous_Raw_NormalizationFirst_AllTime_Wozero/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191010/Continuous_Delta_NormalizationFirst_AllTime_Wozero/SubjCVzero.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                
                ContinuousDeltaRD = RawSubjVar - DeltaSubjVar                
            elif OptionDict['DataFlag'] == 'RawDeltaCVCVDelta':#RawとDeltaのCVCVのDelta
                BRawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Progress_report/20191115_invivohuman/Bolus_Continuous/ばらつき時系列/Bolus/SubjCVCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                BDeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Progress_report/20191115_invivohuman/Bolus_Continuous/ばらつき時系列/Bolus/SubjCVCV_Delta.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                BolusDeltaRD = BRawSubjVar - BDeltaSubjVar;BolusDeltaRDSort=BolusDeltaRD.copy()
                BolusDeltaRDSort['MolColor'] = MolColor;BolusDeltaRDSort=BolusDeltaRDSort.sort_values(by=0)
                if any(OptionDict['SignLabel'])==1:BolusDeltaRDSort=BolusDeltaRDSort.T[OptionDict['SignLabel']].T.sort_values(by=0)
                GH.mkSortedBarWHist(list(BolusDeltaRDSort[0]), 'Bolus_R-D', 'Molecule', 'Bolus_R-D', list(BolusDeltaRDSort['MolColor']), list(BolusDeltaRDSort.index), 10, 15,10, dict(), save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
                #BのR vs D
                Optiondict={'calcR':'spearman','xlabel':'Bolus_Raw', 'ylabel':'Bolus_Delta','Annotate':1,'Label':list(BRawSubjVar.index),
                    'title':'Bolus_RawvsDelta',
                    'y=x' : 1#y=xの線を足す
                    }
                #if any(BolusDeltaRDSort['SignLabel'])==1:ContinuousDeltaRD=ContinuousDeltaRD.T[OptionDict['SignLabel']].T;BolusDeltaRD=BolusDeltaRD.T[OptionDict['SignLabel']].T
                GH.mkScatterWHist(list(BRawSubjVar[0]),list(BDeltaSubjVar[0]),save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム

                RawSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Progress_report/20191115_invivohuman/Bolus_Continuous/ばらつき時系列/Continuous/SubjCVCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                DeltaSubjVar = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Progress_report/20191115_invivohuman/Bolus_Continuous/ばらつき時系列/Continuous/SubjCVCV_Delta.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
                
                ContinuousDeltaRD = RawSubjVar - DeltaSubjVar;    ContinuousDeltaRDSort = ContinuousDeltaRD.copy()
                ContinuousDeltaRDSort['MolColor'] = MolColor ;ContinuousDeltaRDSort=ContinuousDeltaRDSort.sort_values(by=0)
                if any(OptionDict['SignLabel'])==1:ContinuousDeltaRDSort=ContinuousDeltaRDSort.T[OptionDict['SignLabel']].T.sort_values(by=0)
                GH.mkSortedBarWHist(list(ContinuousDeltaRDSort[0]), 'Continuous_R-D', 'Molecule', 'Continuous_R-D', list(ContinuousDeltaRDSort['MolColor']), list(ContinuousDeltaRDSort.index), 10, 15,10, dict(), save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
                #BのR vs D
                Optiondict={'calcR':'spearman','xlabel':'Continuous_Raw', 'ylabel':'Continuous_Delta','Annotate':1,'Label':list(BRawSubjVar.index),
                    'title':'continuous_RawvsDelta',
                    'y=x' : 1#y=xの線を足す
                    }
                #if any(BolusDeltaRDSort['SignLabel'])==1:ContinuousDeltaRD=ContinuousDeltaRD.T[OptionDict['SignLabel']].T;BolusDeltaRD=BolusDeltaRD.T[OptionDict['SignLabel']].T
                GH.mkScatterWHist(list(RawSubjVar[0]),list(DeltaSubjVar[0]),save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム

    #BR-CR, BD-CD
                BolusContinuousDeltaRR = BRawSubjVar - RawSubjVar;BolusContinuousDeltaRR =BolusContinuousDeltaRR.copy()
                BolusContinuousDeltaRR['MolColor'] = MolColor;BolusContinuousDeltaRDSort=BolusContinuousDeltaRR.sort_values(by=0)
                if any(OptionDict['SignLabel'])==1:BolusContinuousDeltaRDSort=BolusContinuousDeltaRDSort.T[OptionDict['SignLabel']].T.sort_values(by=0)
                GH.mkSortedBarWHist(list(BolusContinuousDeltaRDSort[0]), 'BolusContinuous_R-R', 'Molecule', 'BolusContinuous_R-R', list(BolusContinuousDeltaRDSort['MolColor']), list(BolusContinuousDeltaRDSort.index), 10, 15,10, dict(), save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する

                BolusContinuousDeltaRD = BDeltaSubjVar - DeltaSubjVar;BolusContinuousDeltaRD =BolusContinuousDeltaRD .copy()
                BolusContinuousDeltaRD['MolColor'] = MolColor;BolusContinuousDeltaRDSort=BolusContinuousDeltaRD.sort_values(by=0)
                if any(OptionDict['SignLabel'])==1:BolusContinuousDeltaRDSort=BolusContinuousDeltaRDSort.T[OptionDict['SignLabel']].T.sort_values(by=0)
                GH.mkSortedBarWHist(list(BolusContinuousDeltaRDSort[0]), 'BolusContinuous_D-D', 'Molecule', 'BolusContinuous_D-D', list(BolusContinuousDeltaRDSort['MolColor']), list(BolusContinuousDeltaRDSort.index), 10, 15,10, dict(), save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
    
    #R-DのBvsC
                Optiondict={'calcR':'spearman',#'spearman'
                'xlabel':'Bolus_R-D',
                    'ylabel':'Continuous_R-D',
                    'Annotate':1,
                    'Label':list(BolusDeltaRD.index),
                    'title':'',
                    'y=x' : 1#y=xの線を足す

                    }
                #if any(BolusDeltaRDSort['SignLabel'])==1:ContinuousDeltaRD=ContinuousDeltaRD.T[OptionDict['SignLabel']].T;BolusDeltaRD=BolusDeltaRD.T[OptionDict['SignLabel']].T
                GH.mkScatterWHist(list(BolusDeltaRD[0]),list(ContinuousDeltaRD[0]),save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム

                #B,CのR-R vs D-D
                Optiondict={'calcR':'spearman',#'spearman'
                'xlabel':'R-R',
                    'ylabel':'D-D',
                    'Annotate':1,
                    'Label':list(BolusContinuousDeltaRD.index),
                    'title':'RR-DD',
                    'y=x' : 1#y=xの線を足す

                    }
                #if any(BolusDeltaRDSort['SignLabel'])==1:ContinuousDeltaRD=ContinuousDeltaRD.T[OptionDict['SignLabel']].T;BolusDeltaRD=BolusDeltaRD.T[OptionDict['SignLabel']].T
                GH.mkScatterWHist(list(BolusContinuousDeltaRR[0]),list(BolusContinuousDeltaRD[0]),save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム

                #B,CのR-R vs　FastingCV-FastingCV
                BDF = pd.read_excel(save_dir+'Continuous/_MolFeature_temp_CV_SubjB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)            
                CDF = pd.read_excel(save_dir+'Continuous/_MolFeature_temp_CV_SubjC.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)            
                Optiondict={'calcR':'spearman',#'spearman'
                'xlabel':'CVCV_B-C',
                    'ylabel':'FastingCV_B-C',
                    'Annotate':1,
                    'Label':list(BolusContinuousDeltaRD.index),
                    'title':'CVCV-FastingCV',
                    'y=x' : 1#y=xの線を足す
                    }
                #if any(BolusDeltaRDSort['SignLabel'])==1:ContinuousDeltaRD=ContinuousDeltaRD.T[OptionDict['SignLabel']].T;BolusDeltaRD=BolusDeltaRD.T[OptionDict['SignLabel']].T
                GH.mkScatterWHist(list(BolusContinuousDeltaRR[0]),list(BDF['Fasting']-CDF['Fasting']),save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム

                #B,のCVCV vs　FastingCV
                BDF = pd.read_excel(save_dir+'Continuous/_MolFeature_temp_CV_SubjB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)            
                Optiondict={'calcR':'spearman',#'spearman'
                'xlabel':'CVCV_B',
                    'ylabel':'FastingCV_B',
                    'Annotate':1,
                    'Label':list(BolusContinuousDeltaRD.index),
                    'title':'CVCV-FastingCV_B',
                    'y=x' : 1#y=xの線を足す
                    }
                #if any(BolusDeltaRDSort['SignLabel'])==1:ContinuousDeltaRD=ContinuousDeltaRD.T[OptionDict['SignLabel']].T;BolusDeltaRD=BolusDeltaRD.T[OptionDict['SignLabel']].T
                GH.mkScatterWHist(list(BRawSubjVar[0]),list(BDF['Fasting']),save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム
                #C,のCVCV vs　FastingCV
                CDF = pd.read_excel(save_dir+'Continuous/_MolFeature_temp_CV_SubjC.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)            
                Optiondict={'calcR':'spearman',#'spearman'
                'xlabel':'CVCV_C',
                    'ylabel':'FastingCV_C',
                    'Annotate':1,
                    'Label':list(BolusContinuousDeltaRD.index),
                    'title':'CVCV-FastingCV_C',
                    'y=x' : 1#y=xの線を足す
                    }
                #if any(BolusDeltaRDSort['SignLabel'])==1:ContinuousDeltaRD=ContinuousDeltaRD.T[OptionDict['SignLabel']].T;BolusDeltaRD=BolusDeltaRD.T[OptionDict['SignLabel']].T
                GH.mkScatterWHist(list(RawSubjVar[0]),list(CDF['Fasting']),save_dir,MolColor,Optiondict)#2つのリストの散布図+ヒストグラム

                if  OptionDict['BCFlag'] == 'ZEachTime':#各時点で標準化も
                    BolusDeltaRD = BolusDeltaRD.apply(zscore,axis=0);ContinuousDeltaRD= ContinuousDeltaRD.apply(zscore,axis=0)
                kanpachi=pd.concat([BolusDeltaRD,ContinuousDeltaRD],axis=0,join_axes=[BolusDeltaRD.columns]) 
                SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([['Bolus','Continuous'], list(BolusDeltaRD.index)]),columns=list(BolusDeltaRD.columns))
                SubjCV = kanpachi    ;cmap='bwr'                
                
            OptionDict['DataFlag'] = OptionDict['CondFlag']+'_'+OptionDict['DataFlag']+'_NormalizationFirst_'+axislabelname

        if (OptionDict['DataFlag'] == 'Raw') or ( OptionDict['DataFlag'] == 'RawDeltaFC'):
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean   
            OptionDict['DataFlag'] = OptionDict['DataFlag']+'_NormalizationFirst_'+axislabelname
        elif OptionDict['DataFlag'] == 'Delta':    
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean
            if axislabel ==3:
                SubjCV = SubjCV.drop(0)
                SubjCV = np.log10(SubjCV)
            OptionDict['DataFlag'] = 'Delta_NormalizationFirst_'+axislabelname            
        elif (OptionDict['DataFlag'] =='FC') or (OptionDict['DataFlag'] == 'RawFC'): #FCを作ってのちにRawと結合   
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean   
            if axislabel ==3:
                SubjCV = SubjCV.drop(0)
                SubjCV = np.log10(SubjCV)
            OptionDict['DataFlag'] = OptionDict['DataFlag'] +'_NormalizationFirst_' +axislabelname
            

        elif OptionDict['DataFlag'] =='RawDeltaFC':    
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean 
            #if axislabel ==3:
             #   SubjCV = SubjCV.drop(0)
              #  SubjCV = np.log10(SubjCV)
    elif OptionDict['CondFlag'] == 'Bolus':#BolusContinuousならここでMultIndexDF受け取る
        PanelOptionDict={'Data':'Raw'#Raw,Delta,FC ##Bolus, ContinuousのPanel作るときの正規化
            };#DrawTmCsDict['Data']=PanelOptionDict['Data']
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/'
        SubjectName=['isaki','iwashi','karei','shimaaji','unagi']
        SubjTmCsDeltaPanel = CWH.mkMultBolus(SubjectName,PanelOptionDict)
    elif OptionDict['CondFlag'] == 'Continuous':#BolusContinuousならここでMultIndexDF受け取る
        PanelOptionDict={'Data':'Raw'#Raw,Delta,FC ##Bolus, ContinuousのPanel作るときの正規化
            };#DrawTmCsDict['Data']=PanelOptionDict['Data']
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/'
        SubjectName=['isaki','iwashi','karei','shimaaji','unagi']
        SubjTmCsDeltaPanel = CWH.mkMultContinuous(SubjectName,file_dir,PanelOptionDict)
    if  OptionDict['NormalizeFist_DF']==1:#はじめに時系列を正規化するとき、DeltaとFCの0点を消す
        SubjCV = SubjCV.drop(0)
### 20191008_平均の平均で割る
    if OptionDict['NormalizeFistDouble']==1:#はじめに時系列を正規化dounble
        TempAdjustNFD(OptionDict,axislabel,axislabelname,SubjectName,save_dir)
 
### 20191008_標準偏差の平均で割る（DeltaならFig5になるはず）
    elif  OptionDict['NormalizeByStd']==1: # #各種正規化あとに、標準偏差の時間平均で正規化する
        OptionDict,SubjTmCsDeltaPanel,Subjmean,Subjmean_Z,SubjCV = mkNormalizeByStd(OptionDict,SubjTmCsDeltaPanel,axislabel,axislabelname,SubjectName,save_dir)
        if OptionDict['ThirdFlag']=='EachSubj':#RawのStd平均で各個人を割ったもの
            SubjCV = OptionDict['RawDelta_NormalizationByStd_EachSubj']
        
    if axislabel ==0:#0で時間方向、1で分子方向、2で何もしない,3で0点消してlog10
        SubjCV=SubjCV.apply(zscore,axis=0);        cmap='bwr' 
    elif axislabel ==1:
        SubjCV = SubjCV.T.apply(zscore,axis=0).T;        cmap='bwr' 
    if not os.path.isdir(save_dir+OptionDict['DataFlag']):
        os.makedirs(save_dir+OptionDict['DataFlag'])  
    SubjTmCsDeltaPanel.std(level=1,ddof=0).to_excel(save_dir+OptionDict['DataFlag']+'/temp_std.xlsx')  
    SubjTmCsDeltaPanel.var(level=1,ddof=0).to_excel(save_dir+OptionDict['DataFlag']+'/temp_var.xlsx') 
    try:
        Subjmean.to_excel(save_dir+OptionDict['DataFlag']+'/temp_mean.xlsx')
        Subjmean_Z.to_excel(save_dir+OptionDict['DataFlag']+'/temp_meanZ.xlsx')
    except:
        pass
    SubjCV.to_excel(save_dir+OptionDict['DataFlag']+'/SubjCV.xlsx')
    SubjCV.to_excel(save_dir+OptionDict['DataFlag']+'/SubjCV.xlsx')

    
    #SubjCV=SubjCV[['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate',#'Total bile acid',
     #             'Citrate','Cortisol','Free fatty acid','Total ketone body',#'Glutamic acid',
      #             'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate','Growth hormone']]

    
###  各正規化の時系列を描画_Bolus Continuous
    list(SubjTmCsDeltaPanel.columns)
    DrawTmCsDict={'Draw': 'NormalizeFirst_Single_double_FC',#'_Triple':3種類を描画する, 'EachSubject_Triple': #1つの条件で、被験者個別なら　'NormalizeFirst_Single': #1つの枠に、3条件同時描画なら 'EachSubject_Variance_Triple': #Raw, Delta, FCの濃度変化とばらつき変化 NormalizeFirst_Single_double : RawとDeltaだけ 'NormalizeFirst_Single_double_FC': #1つの枠に、2条件同時描画なら_FC
                  'BothDelta' : 'RawDeltaFC'}#差分を描画する'RawDeltaFC':#3種類を描画する}
    #Label =['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
     #              'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
      #             'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate','Growth hormone']
    if any(OptionDict['SignLabel'])==1:
        Label=OptionDict['SignLabel']
        SubjTmCsDeltaPanel=SubjTmCsDeltaPanel[Label]
    else:
        Label=list(SubjTmCsDeltaPanel.columns)
    #CWH.drawTmCsMultiIdxDF(SubjTmCsDeltaPanel,SubjTmCsDeltaPanel, Label,DrawTmCsDict,SubjectName,save_dir+OptionDict['DataFlag']+'/')#MultiIdx形式のSubjectxTimex項目の時系列を描画する。
    
###  濃度変化とばらつき変化    
    OptionDict['SubjTmCsPanel'] = SubjTmCsDeltaPanel
    OptionDict['Raw'] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191006/Raw_NormalizationFirst_1/SubjCV.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    OptionDict['Delta'] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191006/Delta_NormalizationFirst_1/SubjCV.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    OptionDict['FC'] = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191006/FC_NormalizationFirst_1/SubjCV.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    #AdjustCVTmCs(OptionDict,SubjectName,save_dir)
    
    
    #SubjmeanFasting = pd.read_excel(save_dir+'temp_fasting.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
    #SubjmeanFasting_Z = SubjmeanFasting.apply(zscore, axis=0)
    #SubjmeanFasting_Z.to_excel(save_dir+'SubjmeanFasting_Z.xlsx')
    try:
        if SubjCV.loc[0].isnull().all() == 1:   #nanの行も削る(Time=0?)
            SubjCV=SubjCV.drop(0)
        #SubjCV=SubjCV.drop('Threonate',axis=1)
    except:
        pass
    try:
        SubjCV=SubjCV.replace([np.inf, -np.inf], np.nan).dropna(axis=1) ## infがある分子は削る
    except:
        ('error_l1686')
        pass
    if (OptionDict['DataFlag'] == 'Delta_old'): ### DeltaかFCならsqrt{CV^2}
        SubjCV=np.sqrt(SubjCV*SubjCV)
        SubjCV.to_excel(save_dir+'SubjCV_Square_root.xlsx')
        
        
### temp_20190916 20人揃っていなくても、CVナラいいのでは？ CRPも復活_SeveralPropertyHelper L312と対応
     
    if OptionDict['DataFlag'] == 'FC':#一応平均と標準偏差の関係を調べておく
        if SubjCV.loc[0].isnull().all() == 1:   #nanの行も削る(Time=0?)
            SubjCV=SubjCV.drop(0)
        #SubjTmCsDeltaPanel=SubjTmCsDeltaPanel.transpose('items','minor', 'major') #被験者、分子、時点になるはず？
        #plotAveStd(SubjTmCsDeltaPanel, save_dir)#平均と標準偏差の関係
            SubjCV.to_excel(save_dir+OptionDict['DataFlag']+'/SubjCV_Z_'+axislabelname+'.xlsx')
    if OptionDict['DataFlag'] == 'RawFC_NormalizationFirst_' +axislabelname: #RawFC
            RawSubjCV = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20190926/Raw_No/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
            if     axislabel in [0,1]:#0で時間方向、1で分子方向、2で何もしない,3で0点消してlog10
                RawSubjCV=RawSubjCV.apply(zscore,axis=axislabel)
            RawFCSubjCV = pd.concat([RawSubjCV,SubjCV],axis=0,join_axes=[RawSubjCV.columns])
            if axislabel==4:#2つ合わせて時間方向に正規化
                RawFCSubjCV=RawFCSubjCV.apply(zscore,axis=0);cmap='bwr' 
            RawFCSubjCV.to_excel(save_dir+'RawFCSubjCV_'+axislabelname+'.xlsx')
### RawFCのAUC比を計算する
            calcRawFCAUCratio(RawSubjCV,SubjCV,axislabelname,OptionDict,save_dir)
            SubjCV = RawFCSubjCV


    if OptionDict['DataFlag'] =='RawDeltaFC_NormalizationFirst_'+axislabelname:
            cmap='Reds'
        ### Delta   
            DeltaSubjCV = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191006/Delta_NormalizationFirst_1/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)#[list(SubjCV.columns)] 
            #DeltaSubjCV = DeltaSubjCV.drop(0)
            if   axislabel == 0:#0で時間方向、1で分子方向、2で何もしない,3で0点消してlog10
                DeltaSubjCV=DeltaSubjCV.apply(zscore,axis=axislabel); cmap='bwr'
                if  OptionDict['NormalizeFist_DF']==1:#はじめに時系列を正規化するとき、DeltaとFCの0点を消す
                    DeltaSubjCV = DeltaSubjCV.drop(0) 
            elif axislabel ==1:
                DeltaSubjCV = DeltaSubjCV.T.apply(zscore,axis=0).T;        cmap='bwr' 
                if  OptionDict['NormalizeFist_DF']==1:#はじめに時系列を正規化するとき、DeltaとFCの0点を消す
                    DeltaSubjCV = DeltaSubjCV.drop(0,axis=0)  
            elif axislabel ==3:
                DeltaSubjCV = DeltaSubjCV.drop(0); cmap='bwr'
                DeltaSubjCV = np.log10(DeltaSubjCV)                
            elif axislabel ==5:
                DeltaSubjCV = DeltaSubjCV.drop(0); cmap='bwr'
 
            RawDeltaSubjCV = pd.concat([SubjCV,DeltaSubjCV],axis=0,join_axes=[SubjCV.columns])
        ### FC
            FCSubjCV = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20191006/FC_NormalizationFirst_1/SubjCV.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)#[list(SubjCV.columns)] 
            if axislabel ==0:#0で時間方向、1で分子方向、2で何もしない,3で0点消してlog10
                FCSubjCV=FCSubjCV.apply(zscore,axis=0);        cmap='bwr' 
                if  OptionDict['NormalizeFist_DF']==1:#はじめに時系列を正規化するとき、DeltaとFCの0点を消す
                    FCSubjCV = FCSubjCV.drop(0) 
            elif axislabel ==1:
                FCSubjCV = FCSubjCV.T.apply(zscore,axis=0).T;        cmap='bwr' 
                if  OptionDict['NormalizeFist_DF']==1:#はじめに時系列を正規化するとき、DeltaとFCの0点を消す
                    FCSubjCV = FCSubjCV.drop(0,axis=0) 
            elif axislabel ==3:
                FCSubjCV = FCSubjCV.drop(0); cmap='bwr'
                FCSubjCV=np.log10(FCSubjCV)
            elif axislabel ==5:
                 FCSubjCV = FCSubjCV.drop(0); cmap='bwr'
              
            RawDeltaFCSubjCV = pd.concat([RawDeltaSubjCV,FCSubjCV],axis=0,join_axes=[RawDeltaSubjCV.columns])
               
            if   axislabel in [4,5]:#0で時間方向、1で分子方向、2で何もしない,3で0点消してlog10, 4で全て繋げてから時間方向に標準化
                #RawDeltaFCSubjCV = RawDeltaFCSubjCV.drop(0)#Rawの0も消えちゃう
                RawDeltaFCSubjCV=RawDeltaFCSubjCV.apply(zscore,axis=0) ; cmap='bwr'
            if   axislabel ==6:#6で条件の時点間で正規化した時系列をつなげる、
                SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(RawDeltaFCSubjCV),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(SubjCV.index)]),columns=list(RawDeltaFCSubjCV.columns))
                RawDF = SubjTmCsDeltaPanel.unstack(level=1).apply(zscore,axis=0,ddof=0).stack(level=1).loc['Raw']
                DeltaDF = SubjTmCsDeltaPanel.unstack(level=1).apply(zscore,axis=0,ddof=0).stack(level=1).loc['Delta']
                FCDF = SubjTmCsDeltaPanel.unstack(level=1).apply(zscore,axis=0,ddof=0).stack(level=1).loc['FC']
                RawDeltaFCSubjCV = pd.concat([RawDF,DeltaDF,FCDF],axis=0,join_axes=[RawDF.columns]) ; cmap='bwr'
            if   axislabel ==7:#7で各時系列の差分をつなげる
                SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(RawDeltaFCSubjCV),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(SubjCV.index)]),columns=list(RawDeltaFCSubjCV.columns))

                RawDF = SubjTmCsDeltaPanel.loc['Raw'];  DeltaDF = SubjTmCsDeltaPanel.loc['Delta']; FCDF =SubjTmCsDeltaPanel.loc['FC']
                RD = RawDF-DeltaDF;RF=RawDF-FCDF; DF=DeltaDF-FCDF
                RawDeltaFCSubjCV = pd.concat([RD,RF,DF],axis=0,join_axes=[RD.columns]) ; cmap='bwr'
            if   axislabel ==8:#7で各時系列の差分のつなげて分子内標準化
                SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(RawDeltaFCSubjCV),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(SubjCV.index)]),columns=list(RawDeltaFCSubjCV.columns))

                RawDF = SubjTmCsDeltaPanel.loc['Raw'];  DeltaDF = SubjTmCsDeltaPanel.loc['Delta']; FCDF =SubjTmCsDeltaPanel.loc['FC']
                RD = RawDF-DeltaDF;RF=RawDF-FCDF; DF=DeltaDF-FCDF
                RawDeltaFCSubjCV = pd.concat([RD,RF,DF],axis=0,join_axes=[RD.columns]) ; cmap='bwr'
                RawDeltaFCSubjCV = RawDeltaFCSubjCV.apply(zscore,axis=0,ddof=0)
            if   axislabel ==9:#9で各時系列の差分をつなげる2つだけ
                SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(RawDeltaFCSubjCV),index=pd.MultiIndex.from_product([['Raw','Delta','FC'], list(SubjCV.index)]),columns=list(RawDeltaFCSubjCV.columns))

                RawDF = SubjTmCsDeltaPanel.loc['Raw'];  DeltaDF = SubjTmCsDeltaPanel.loc['Delta']; FCDF =SubjTmCsDeltaPanel.loc['FC']
                RD = RawDF-DeltaDF;RF=RawDF-FCDF; #DF=DeltaDF-FCDF
                RawDeltaFCSubjCV = pd.concat([RD,RF],axis=0,join_axes=[RD.columns]) ; cmap='bwr'
                
            OptionDict['DataFlag'] = 'RawDeltaFC_NormalizationFirst_'   +axislabelname

            RawDeltaFCSubjCV.to_excel(save_dir+OptionDict['DataFlag']+'/RawDeltaFCSubjCV.xlsx')
            SubjCV = RawDeltaFCSubjCV

    if not os.path.isdir(save_dir+OptionDict['DataFlag']):
        os.makedirs(save_dir+OptionDict['DataFlag'])  

### クラスタリング  描画分子選択  
    #SubjCV=SubjCV.drop('hs-CRP',axis=1)   
    #SubjCV=SubjCV[['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate',
     #              'Free fatty acid',
      #             'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate',
       #        ]]
            
           # ['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
            #      'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
             #      'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate','Growth hormone']]
    #SubjCV = SubjCV[['Total ketone body','Free fatty acid','Malate','Lactate','cis-Aconitate','Citrate','Pyruvate','Succinate','Glucose',
     #      'Alanine','Arginine', 'Asparagine', 'Citrulline','Leucine','Glutamic acid','Glutamine','Histidine','Lysine','Valine',
      #     'Methionine','Ornithine','Phenylalanine','Proline','Serine','Threonine','Tryptophan','Tyrosine','Taurine','Creatine','Creatinine','Betaine','Choline','Glycerophosphate','Urate']]
    
    #0,30,60,90,120dale
    #SubjCV = SubjCV.loc[[0,30,60,90,120]]
    #metab=['Triglyceride','Total cholesterol','LDL cholesterol','HDL cholesterol','Insulin','Glucose','Lactate','Alanine','Arginine', 
    #'Citrulline','Cystine','Glutamic acid','Glutamine','Glycine','Histidine','Hydroxyproline', 'Asparagine',
    #'Isoleucine','Leucine','Lysine','Ornithine','Phenylalanine','Proline','Serine','Threonine','Taurine',
    #'Tryptophan','Tyrosine','Valine','3-Methylhistidine']
    #SubjCV = SubjCV[metab]
    
### CVCVを作成
    SubjCVStd = SubjCV.std(ddof=0,axis=0); SubjCVvar = SubjCV.var(ddof=0,axis=0);
    SubjCVMean = SubjCV.mean(axis=0)
### 20200114_ばらつき時系列の時間平均棒グラフ
    OptionDict['save_dir'] =save_dir
    #GH.mkSortedBarWHist(list((SubjCVMean)), 'log10(FCVarTimeMean)', 'Molecular Name', 'log10(FCVarTimeMean)', LH.MolColor(list(SubjCVMean.index)),list(SubjCVMean.index), 15, 15,15, OptionDict, OptionDict['save_dir'])#ヒストグラム月の任意のリストを代入して棒グラフを描画する
    SubjCVMeanSort = SubjCVMean.sort_values(0)
    #GH.mkSortedBarWHist(list((SubjCVMeanSort)), 'log10(FCVarTimeMean)_Sort', 'Molecular Name', 'log10(FCVarTimeMean)', LH.MolColor(list(SubjCVMeanSort.index)),list(SubjCVMeanSort.index), 15, 15,15, OptionDict, OptionDict['save_dir'])#ヒストグラム月の任意のリストを代入して棒グラフを描画する

    SubjCVStd.to_excel(save_dir+OptionDict['DataFlag']+'/SubjCVStd.xlsx'); SubjCVMean.to_excel(save_dir+OptionDict['DataFlag']+'/SubjCVMean.xlsx');SubjCVvar.to_excel(save_dir+OptionDict['DataFlag']+'/SubjCVvar.xlsx')
    SubjCVCV = SubjCVStd / SubjCVMean
    SubjCVCV.to_excel(save_dir+OptionDict['DataFlag']+'/SubjCVCV.xlsx')    
### temp_Delta用_20190926
    #SubjCV=SubjCV.drop('Total bile acid',axis=1)
    #SubjCV=SubjCV.drop('Glutamic acid',axis=1)
    #SubjCV = SubjCV.drop('Growth hormone',axis=1);
    #SubjCV = SubjCV.drop('hs-CRP',a xis=1);
    ClstMol,Colordict,ClstDF,ClstColorDF,ColorDF = mD.draw_heatmap(SubjCV.T,1,'MolColor',{'title':'mean'+OptionDict['DataFlag']},save_dir+OptionDict['DataFlag']+'/',cmap='PuOr_r')#'bwr' ('YlOrRd'),'Reds'    
    ClstAveTimeSeriesDF = ACD.mkClstAveTimeSign(save_dir,SubjCV.T,ClstMol)#クラスタごとに分けた平均をだす

    if any(OptionDict['SignLabel'])==1:SubjCV=SubjCV[OptionDict['SignLabel']]

    if OptionDict['ThirdFlag']!='EachSubj':#RawのStd平均で各個人を割ったもの
        ClstMol,Colordict,ClstDF,ClstCVColorDF,ColorDF = mD.draw_heatmap(SubjCV.T,1,'MolColor',{'title':'CV'+OptionDict['DataFlag']},save_dir+OptionDict['DataFlag']+'/',cmap='bwr')#'bwr' ('YlOrRd'),'Reds'
        ClstAveTimeSeriesDF = ACD.mkClstAveTimeSign(save_dir,SubjCV.T,ClstMol)#クラスタごとに分けた平均をだす
    else:
        print('l_2255')
        #ClstAveTimeSeriesDF=pd.DataFrame(); ClstColorDF=pd.DataFrame(); ColorDF=pd.DataFrame();
    EachSwitch=0#個別のタイムコース作成するなら1
    #各Cluster, 分子の時系列
    #mD.plotTimeCourse(ClstMol,Colordict,SubjCV,EachSwitch,save_dir)#変動指標値、クラスターごと個別+平均値タイムコース描画
    #CompareClstMol = mD.mkCompareTableFile(ClstMol)#住友さんのクラスターごとの分子と、本解析のクラスター毎の分子で重複があるならファイルにする
    #mD.mkCompareTableLtx('','')
    #mD.mkClusterMolText(CompareClstMol)#クラスター数だけテキストファイルつくる
    #mD.ForCompareTableLtx(CompareClstMol)#クラスター数だけ＝行毎テキストファイルつくる

### PCA
    TimeCourseAvePCA(SubjCV.T,ClstAveTimeSeriesDF,ClstColorDF,ColorDF,save_dir+OptionDict['DataFlag']+'/')#PCAする

    Optiondict={'Color':'MolColor',
                'Annotate' : 1,#分子名などをAnnotateするなら
                'Label':list(SubjCV.columns)
                }
### CV時系列
    #SubjCV= SubjCV[['Glucose','LDL cholesterol']]
    #plotTmCs(SubjCV,Optiondict,save_dir+OptionDict['DataFlag']+'/'+axislabelname)#CVノタイムコース受け取って1つのグラフに重ねて描画

### 各時点CVの時間平均と標準偏差のplot
    #AdjustPlot(SubjCVMean,SubjCVStd,OptionDict,save_dir+OptionDict['DataFlag']+'/')
### ばらつきのなさの指標
    xlabel = 'SubjCVCV'#'ISCO'or'SubjStd'
    #AdjustISCO(SubjCVCV,xlabel,OptionDict,save_dir)#xlabelでファイル名変更可
    #1-SubjCVCV
def calcRawFCAUCratio(RawCV,FCCV,axislabelname,OptionDict,save_dir):#RawとFCのばらつき時系列のAUCを計算する
    if any(OptionDict['SignLabel'])==1:RawCV=RawCV[OptionDict['SignLabel']];FCCV=FCCV[OptionDict['SignLabel']]

    MolName = list(RawCV.columns)
    AUCDF = pd.DataFrame(data=None,index=['Raw','FC'],columns=MolName)
    for i in ['Raw','FC']:
        for j in MolName:
            if i=='Raw':                
                CVList = list(RawCV[j])
            else:
                CVList=list(FCCV[j])
            TimeList = list(RawCV.index)
            AUCDF.loc[i,j] = AmT.CalcAUC(CVList,TimeList)
    AUCDF.to_excel(save_dir+'AUCratiotoRaw.xlsx')
    RawtoFCratio = list(np.array(AUCDF.loc['FC']) / np.array(AUCDF.loc['Raw']))
    RawtoFCratioDF = pd.DataFrame(data=None,index=[0],columns=MolName);RawtoFCratioDF.loc[0]=RawtoFCratio
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)

    ### ソートして棒グラフを棒グラフを描画する
    RawtoFCratioColor=pd.concat([RawtoFCratioDF,MolColor.T],axis=0,join_axes=[RawtoFCratioDF.columns])
    RawtoFCratioColorSort = RawtoFCratioColor.T.sort_values(by=0).T
    RawtoFCratioColorSort.to_excel(save_dir+'AUCratiotoRawWColor_dump.xlsx')
    List1= list(RawtoFCratioColorSort.loc[0]); Title='AUCratio'; MolColor = list(RawtoFCratioColorSort.loc['MolColor']); MolNameList = list(RawtoFCratioColorSort.columns)
    GH.mkSortedBarWHist(List1, Title, '', 'AUCratio (Fold change/Raw)', MolColor, MolNameList, 60, 20,0, dict(),save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する

def AdjustPlot(list1,list2,OptionDict,save_dir):
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)
    ColorList = list(MolColor['MolColor'])
        
    OptionDict['Annotate'] = 1#分子名などをAnnotateするなら
    OptionDict['Label']= list(list1.index)
    OptionDict['calcR']='pearson';OptionDict['xlabel']='CVMean';OptionDict['ylabel']='CVStd';
    OptionDict['title']='CVAve_CVStd'
    GH.mkScatterWHist(list1,list2,save_dir,ColorList,OptionDict)#2つのリストの散布図+ヒストグラム

    
def AdjustISCO(SubjCVCV,xlabel,OptionDict,save_dir):#ばらつきのなさの指標の棒グラフ描画
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)
    
    MolCVCVDictSortWColor=pd.concat([SubjCVCV.T,MolColor],axis=1,join_axes=[SubjCVCV.T.index])
    NewDF4CVCV = MolCVCVDictSortWColor.sort_values(by=0)
    ColorsCV = list(NewDF4CVCV['ClstColor'])
    
    MolColorsCV = list(NewDF4CVCV['MolColor'])
    
    #OptionDict={}
    OptionDict['Annotate'] = 1#分子名などをAnnotateするなら
    OptionDict['calcR']='pearson';OptionDict['xlabel']=xlabel;OptionDict['ylabel']='ISPS';Title='ISCO'
    #GH.mkSortedBarWHist(List1, Title, '', 'ISCO', MolColorsCV, XtickLabel, 60, 30,0, dict(),save_dir+OptionDict['DataFlag']+'/')#ヒストグラム月の任意のリストを代入して棒グラフを描画する

    ##### IDIの小さい順#np.reciprocal(
    #List1 = np.sort(np.array(list(MolCVCVDict.values()))) 
    List1 = list(NewDF4CVCV[0])

    fig,ax = plt.subplots(figsize=(4.8,6.4))#np.linspace(0, 1, len(MolCVAveDict)),
    ax.barh(np.linspace(0, len(List1), len(List1)),List1,tick_label=list(NewDF4CVCV.index),color=MolColorsCV)
    #for tick in ax.get_xticklabels():
        #tick.set_rotation(270)#plot.tick_params(axis='both', which='major', labelsize=10)
    plt.gca().tick_params(axis='both',labelsize=20)
    ax.tick_params(axis='y',labelsize=5)
    ax.set_ylabel(OptionDict['xlabel'],size=30)
    #ax.set_ylim(-20,25)

    plt.savefig(save_dir+OptionDict['DataFlag']+'/MolTimeCVCV_WMolColorSorted' + OptionDict['xlabel'] + '.pdf',bbox_inches="tight") 
    plt.savefig(save_dir+OptionDict['DataFlag']+'/MolTimeCVCV_WMolColorSorted' + OptionDict['xlabel'] + '.png',bbox_inches="tight") 

    
def TimeCourseAvePCA(XDF,ClstAveTimeSeriesDF,ClstColorDF,ColorDF,save_dir):#PCAする
    #PC平面プロっト：'Biplot'　＃楕円フィッティング：'DrawEllipse' #各PCの寄与率：'PlotPCACovRatio'
    #各分子のLoading時系列をプロット、クラスター色分け：'LoadingPlotWCluster'　
    AnalSwitch={'Biplot':1,'DrawEllipse':1,'PlotPCACovRatio':0,'LoadingPlotWCluster':1,
                'ScoreHeatMap':0 ,'FactorLoadingHeatMap':0,#Score, FactorLoadingのヒートマップ：'ScoreHeatMap','FactorLoadingHeatMap'
                'ScoreVariance':0,'LoadingPCsimulation':0,#'ScoreVariance':各分子、被験者間のスコアの分散  'LoadingPCsimulation' : #PC1,2固定して片方動かした時の時系列描画 
                'VectorDiagram' : 0,#Bolus->Continuouのベクトル図}
                'calcinner_outer_product' :0, #各平面でのBC間の内積外積を算出、
                'LengVector': 0,
                'WHist':0#ヒストグラムつきにする
        }#PC1,2score間のベクトル長を算出する。}
    
    LoadingOpt = 0 #FactorLoadingのヒートマップ：#0：データ（縦）×主成分（横）, #1：主成分（縦）×データ（横）
    BiplotSwitch={'DotColor':'ClstColor_direct',#'MolColor'：各分子の色、 'ClsterColor'：クラスターの色、'Black':黒、EachMol：各分子の色、'EachMolCluster':分子x○○の時のCluster色'BC':#PCA時に色分けBCにする、'MolColor_AminoGlc':AAのラベルなど変える 'ClstColor_direct'クラスタリングした直後のクラスターの色（昔のClsterの色を使う可能性があるから）,'PC1_2_represen', 'PC1_2'　　　　　　　　　　　　　
                  'Label':'',#'Annotate':DotにLabel併記
                  'EachMolColor':[],
                  'Check':'',#'Amino','Glc','AminoGlc': どちらも
                  'AminoCheck' :'',#'EAA', #'protein','ketogenic','EAA','SemiEAA' アミノ酸や糖代謝色変え
                  'ColorDF':ColorDF
                  }#'EachMolColor':分子x被験者PCAにおける、各分子でのいろわえk
    OptionSwitch={'Target':''}#'EachCluster' : 各Cluster のloading時系列を描画
    CompareBolus={'TargetDetail' : ''}#BC各被験者なら);'EachSubject'
    ClstMolDF=pd.DataFrame(data=None);#ClstColorDF=pd.DataFrame(data=None);#ClstAveTimeSeriesDF=pd.DataFrame(data=None) 
    ##############################################
    MolLabel = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1", index_col=0).index)
    
    #XDF = XDF.T
    #ISSIDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/EachMolInterSubjConnectCorr.xlsx',header=0,encoding = "ISO-8859-1", index_col=0).loc['ISSI']
    #IDIDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjMolTimeCVAve_Eng.xlsx',header=0,encoding = "ISO-8859-1", index_col=0).loc['IDI'].drop('3-Hydroxybutyrate')
    #DF  = pd.concat([XDF,IDIDF],axis=1, join_axes=[XDF.index]) #分子x特徴量
    #DF = pd.concat([DF,ISSIDF],axis=1, join_axes=[DF.index]) #分子x特徴量
    #DF.to_excel(save_dir+'/MolFeature_IDI_ISSI.xlsx')
    #DF=DF.drop('Opeak',axis=1);DF=DF.drop('Gain',axis=1);DF=DF.drop('log2FC',axis=1);DF=DF.drop('Slope3',axis=1);DF=DF.drop('Tpeak',axis=1);DF=DF.drop('Ttimeconctant',axis=1);
    
    labelProp =  XDF.columns
    MolLabel=XDF.index

    PCAPre.PCAprep(XDF,ClstAveTimeSeriesDF,ClstColorDF,MolLabel,labelProp,ClstMolDF,AnalSwitch,LoadingOpt,BiplotSwitch,OptionSwitch,CompareBolus,save_dir)

        
def plotTmCs(TmCsDF,Optiondict,save_dir):#CVのタイムコース受け取って1つのグラフに重ねて描画
    (num_rows,num_cols)=(1,1)
    fig, host = plt.subplots(num_rows,num_cols,figsize=(11.25,7.5))#20,10))
    plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05);molcount=0
    Col=list(TmCsDF.columns)
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,index_col=0)

    if Optiondict['Color'] == 'MolColor':
        MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,index_col=0)
    elif Optiondict['Color'] == 'Blue':
        MolColor = pd.DataFrame(data=None,index=Col,columns=['MolColor'])
        MolColor['MolColor'] = ['blue']*len(Col)
    elif Optiondict['Color'] == 'Green':
        MolColor = pd.DataFrame(data=None,index=Col,columns=['MolColor'])
        MolColor['MolColor'] = ['green']*len(Col)
    elif Optiondict['Color'] == 'Red':
        MolColor = pd.DataFrame(data=None,index=Col,columns=['MolColor'])
        MolColor['MolColor'] = ['red']*len(Col)        
    elif type(Optiondict['Color']) == list:
        MolColor = pd.DataFrame(data=None,index=Col,columns=['MolColor'])
        MolColor['MolColor'] = Optiondict['Color']         
    for j in range(0,num_rows): #行番号,一行につき1被験者 
        for i in range(0,num_cols):#列で回す
            for molcount in range(len(Col)):#len(Col)):
            
                host.plot(list(TmCsDF.index), list(TmCsDF[Col[molcount]]),c=MolColor['MolColor'][Col[molcount]], marker="o",linestyle = 'solid',linewidth=2.0,markersize=5.0)
                if Optiondict['Annotate'] == 1:#分子名などをAnnotateするなら
                        AnDict=dict({'Glucose':'Glc','Insulin':'Ins','C-peptide':'CRP','GIP(Active)':'GIP','Pyruvate':'Pyr','Total bile acid':'TBA',
                               'Citrate':'Cit','Cortisol':'Cor','Free fatty acid':'FFA','Total ketone body':'Ketone','Glutamic acid':'Glu',
                               'Citrulline':'Citr','Methionine':'Met','Isoleucine':'Ile','Leucine':'Leu','Tyrosine':'Tyr','4-Methyl-2-oxopentanoate':'4M2O','Glu + threo-beta-methylaspartate':'Glu+TBM','Growth hormone':'GH'})
                        if Optiondict['Label'][molcount] in ['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
                               'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
                               'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate','Growth hormone']:
                    #時系列の最終点に
                            try:#横軸が時間なら
                                host.annotate(AnDict[Optiondict['Label'][molcount]],fontsize=10, xy=(list(TmCsDF.index)[-1]+5,list(TmCsDF[Col[molcount]])[-1]))
                            except:#項目なら
                                host.annotate(AnDict[Optiondict['Label'][molcount]],fontsize=10, xy=(len(TmCsDF.index)-0.8,list(TmCsDF[Col[molcount]])[-1]))
                                
                        else:
                            try:
                                host.annotate(Optiondict['Label'][molcount],fontsize=1, xy=(list(TmCsDF.index)[-1]+5,list(TmCsDF[Col[molcount]])[-1]))
                            except:
                                host.annotate(Optiondict['Label'][molcount],fontsize=1, xy=(len(TmCsDF.index)-0.8,list(TmCsDF[Col[molcount]])[-1]))
                                
                molcount += 1
        try:#もし時間間隔がズレるなら
            #host.set_xticks([0,100,200])
            #host.set_xticklabels(['0','100','200'],size=20)
            host.set_xticks(list(TmCsDF.index))
            host.set_xticklabels([int(list(TmCsDF.index)[i]) for i in range(len(TmCsDF.index))])
        except:    
            host.set_xticks(list(TmCsDF.index))
            host.set_xticklabels(list(TmCsDF.index))            
            
            #host.set_xticks([0,60,120,180,240])
            #host.set_xticklabels(['0','60','120','180','240'])
        plt.rcParams["font.size"] = 25
        host.set_ylabel('Correlation',fontsize=40)#濃度
        xmin, xmax, ymin, ymax = plt.axis()
        #plt.hlines(0.6, xmin, xmax, "black", linestyles='dashed',linewidth=1)
        #plt.hlines(-0.6, xmin, xmax, "black", linestyles='dashed',linewidth=1)

        #host.set_ylabel('CV',fontsize=40)#濃度
        host.set_xlabel("Time (min.)",fontsize=30)#時系列
        host.tick_params(labelsize=15, axis='both')

        plot_axis = plt.axis()
        #host[3,2].axis('off')#いらないところはこれで消す
        #host[3,3].axis('off')
        try:
            plt.title(Optiondict['Title'])
        except:
            pass
        plt.savefig(save_dir + 'Correlation_TimeCourse.pdf',bbox_inches="tight")
        plt.savefig(save_dir + 'Correlation_TimeCourse.png',bbox_inches="tight")

        #plt.savefig(save_dir + 'CV_TimeCourse.pdf',bbox_inches="tight")
        #plt.savefig(save_dir + 'CV_TimeCourse.png',bbox_inches="tight")

        plt.close()

def mkNormalizeByStd(OptionDict,SubjTmCsDeltaPanel,axislabel,axislabelname,SubjectName,save_dir):
    #SubjectName=list(set(SubjTmCsDeltaPanel.index.get_level_values(level=0)))
    if OptionDict['DataFlag'] == 'Raw':
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            #SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean
            SubjCV = Subjmean / SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() #/ Subjmean   
            OptionDict['DataFlag'] = 'Raw_NormalizationByStd'
    elif OptionDict['DataFlag'] == 'Delta':    
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)
        SubjCV = Subjmean /SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() #/ Subjmean
        if axislabel ==3:
            SubjCV = SubjCV.drop(0)
            SubjCV = np.log10(SubjCV)
        OptionDict['DataFlag'] = 'Delta_NormalizationByStd'            
    elif OptionDict['DataFlag'] =='FC':    
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)
        SubjCV = Subjmean /SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() #/ Subjmean   
        if axislabel ==3:
            SubjCV = SubjCV.drop(0)
            SubjCV = np.log10(SubjCV)
        OptionDict['DataFlag'] = 'FC_NormalizationByStd' 
    elif OptionDict['DataFlag'] =='RawDeltaFC':    
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)
        SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean 
        if axislabel ==3:
            SubjCV = SubjCV.drop(0)
            SubjCV = np.log10(SubjCV)
    if OptionDict['SecondFlag']=='RawRaw': #正規化方法、RawRaw:RawのStd平均で除してRawをPanelに                    
        SubjTmCsDeltaPanel = SubjTmCsDeltaPanel / SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() 

        for i in range(len(SubjectName)):
            tempRaw=SubjTmCsDeltaPanel.loc[SubjectName[i]]
            tempRaw.to_excel(save_dir + 'SubjTmCs_'+SubjectName[i]+'DivBymeanRawStdRaw.xlsx')
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)
        SubjCV = Subjmean / SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() #/ Subjmean   
        OptionDict['DataFlag'] = 'RawRaw_NormalizationByStd'
    if OptionDict['SecondFlag']=='RawDelta': #正規化方法、RawDelta:RawのStd平均から差分とる                    
        for i in range(len(SubjectName)):
            if i==0:#差分データ構築
                SubjDelta =    SubjTmCsDeltaPanel.loc[SubjectName[i]] - SubjTmCsDeltaPanel.loc[SubjectName[i],0]    
                tempDelta = SubjDelta / SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() 
            else:
                SubjTmCsDelta = SubjTmCsDeltaPanel.loc[SubjectName[i]] - SubjTmCsDeltaPanel.loc[SubjectName[i],0]   
                SubjDelta=pd.concat([SubjDelta,SubjTmCsDelta],axis=0,join_axes=[SubjDelta.columns]) 
                tempDelta = SubjTmCsDelta / SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() 
            tempDelta.to_excel(save_dir + 'SubjTmCs_'+SubjectName[i]+'DivBymeanRawStdDelta.xlsx')
        SubjTmCsDeltaPanel_rev = pd.DataFrame(data=np.array(SubjDelta),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))
        SubjCV = SubjTmCsDeltaPanel_rev.mean(level=1) / SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() 
        OptionDict['DataFlag'] = 'RawDelta_NormalizationByStd'
        Subjmean = SubjTmCsDeltaPanel_rev.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)     
        if OptionDict['ThirdFlag']=='EachSubj':#RawのStd平均で各個人を割ったもの
            MolName=list(SubjCV.columns); NewCol = [i +'_' + jj  for jj in SubjectName for i in MolName]
            NewDF = pd.DataFrame(data=None,columns=NewCol,index = list(SubjCV.index))
            for ii in  range(len(SubjectName)):           
                NewDF.iloc[:,0+83*ii:83+83*ii]=np.array(SubjTmCsDeltaPanel_rev.loc[SubjectName[ii]])
            NewDF.to_excel(save_dir + 'RawDelta_NormalizationByStd_EachSubj_WNan.xlsx')
            NewDF=NewDF.dropna(axis=1)
            NewDF.to_excel(save_dir + 'RawDelta_NormalizationByStd_EachSubj_WoNan.xlsx')

            OptionDict['RawDelta_NormalizationByStd_EachSubj']  = NewDF
            #SubjTmCsDeltaPanel_rev SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() 

    if OptionDict['SecondFlag']=='RawDeltaEach': #正規化方法、RawDeltaEach:Rawの各時点のStdで除してから差分とる                    
        for i in range(len(SubjectName)):
            if i==0:#差分データ構築
                SubjDelta = SubjTmCsDeltaPanel.loc[SubjectName[i]] - SubjTmCsDeltaPanel.loc[SubjectName[i],0]    
            else:
                SubjTmCsDelta = SubjTmCsDeltaPanel.loc[SubjectName[i]] - SubjTmCsDeltaPanel.loc[SubjectName[i],0]   
                SubjDelta=pd.concat([SubjDelta,SubjTmCsDelta],axis=0,join_axes=[SubjDelta.columns]) 
        SubjTmCsDeltaPanel_rev = pd.DataFrame(data=np.array(SubjDelta),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))

        SubjCV = SubjTmCsDeltaPanel_rev.mean(level=1) / SubjTmCsDeltaPanel.std(level=1,ddof=0)
        OptionDict['DataFlag'] = 'RawDelta_NormalizationByStdEach'
        Subjmean = SubjTmCsDeltaPanel_rev.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)  

    if OptionDict['SecondFlag']=='RawFCEach': #正規化方法、RawFCEach:Rawの各時点のStdで除してからFCとる                    
        for i in range(len(SubjectName)):
            if i==0:#差分データ構築
                SubjDelta =    SubjTmCsDeltaPanel.loc[SubjectName[i]] / SubjTmCsDeltaPanel.loc[SubjectName[i],0]    
            else:
                SubjTmCsDelta = SubjTmCsDeltaPanel.loc[SubjectName[i]] / SubjTmCsDeltaPanel.loc[SubjectName[i],0]   
                SubjDelta=pd.concat([SubjDelta,SubjTmCsDelta],axis=0,join_axes=[SubjDelta.columns]) 
        SubjTmCsDeltaPanel_rev = pd.DataFrame(data=np.array(SubjDelta),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))

        SubjCV = SubjTmCsDeltaPanel_rev.mean(level=1) / SubjTmCsDeltaPanel.std(level=1,ddof=0)
        OptionDict['DataFlag'] = 'RawFC_NormalizationByStdEach'
        Subjmean = SubjTmCsDeltaPanel_rev.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)  
        
    if OptionDict['SecondFlag']=='RawFC': #正規化方法、RawFC:RawのStd平均で除してからFCに                  
        for i in range(len(SubjectName)):
            if i==0:#FCデータ構築
                SubjDelta =    SubjTmCsDeltaPanel.loc[SubjectName[i]] / SubjTmCsDeltaPanel.loc[SubjectName[i],0]    
            else:
                SubjTmCsDelta = SubjTmCsDeltaPanel.loc[SubjectName[i]] / SubjTmCsDeltaPanel.loc[SubjectName[i],0]   
                SubjDelta=pd.concat([SubjDelta,SubjTmCsDelta],axis=0,join_axes=[SubjDelta.columns]) 
        SubjTmCsDeltaPanel_rev = pd.DataFrame(data=np.array(SubjDelta),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))
        SubjCV = SubjTmCsDeltaPanel_rev.mean(level=1) / SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() 
        OptionDict['DataFlag'] = 'RawFC_NormalizationByStd'
        Subjmean = SubjTmCsDeltaPanel_rev.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)   

    if OptionDict['SecondFlag']=='RawRawEach': #正規化方法、RawRawEach:Rawの各時点のStdで除してRawをPanelに                    
        for i in range(len(SubjectName)):
            if i==0:#FCデータ構築
                SubjDelta = SubjTmCsDeltaPanel.loc[SubjectName[i]] / SubjTmCsDeltaPanel.std(level=1,ddof=0)
            else:
                SubjTmCsDelta = SubjTmCsDeltaPanel.loc[SubjectName[i]] / SubjTmCsDeltaPanel.std(level=1,ddof=0)  
                SubjDelta=pd.concat([SubjDelta,SubjTmCsDelta],axis=0,join_axes=[SubjDelta.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(SubjDelta),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))
        Subjmean = SubjTmCsDeltaPanel.mean(level=1)
        Subjmean_Z = Subjmean.apply(zscore, axis=0)
        SubjCV = Subjmean / SubjTmCsDeltaPanel.std(level=1,ddof=0).mean() #/ Subjmean   
        # SubjTmCsDeltaPanel.std(level=1,ddof=0)#各時点でstdは1のはず
        OptionDict['DataFlag'] = 'RawRaw_NormalizationByStdEach'
        
    return(OptionDict,SubjTmCsDeltaPanel,Subjmean,Subjmean_Z,SubjCV)
def heatmap(data_2d,label_name): 
  import mpl_toolkits.axes_grid1
  num_x = data_2d.shape[1]
  num_y = data_2d.shape[0]

  c_max = max([data_2d.max(), abs(data_2d.min())])
  colors = np.array([plt.cm.hsv(i/len(data_2d)) for i in range(len(data_2d))])
  ax = plt.imshow(data_2d,aspect='auto',interpolation='nearest',
             cmap=color_map,vmin=-c_max,vmax=c_max)
  plt.gca().xaxis.tick_top()
  divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
  #cax = divider.append_axes('right', '5%', pad='3%')
  char = plt.colorbar(shrink=1)
  char.ax.tick_params(labelsize=10)
  plt.xticks(range(num_x),label_name,rotation=90,size=7)#np.arange(1,num_x+1).astype('<U'))
  plt.yticks(range(num_y),label_name,size=7)
  #plt.title('Number of principal components')
  bx = plt.gca()
  #label = bx.set_xlabel('Number of principal components', fontsize = 15)
  bx.xaxis.set_label_coords(0.5,1.03)
  #if ChangeColorSwitch == 1:#ラベルの色変え、ラベルが代謝順になっていないと使えないが
  MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolColor_Eng.xlsx')
  [t.set_color(i) for (i,t) in
     zip(list(MolColorDF['MolColor']),plt.gca().yaxis.get_ticklabels())]
  [t.set_color(i) for (i,t) in
     zip(list(MolColorDF['MolColor']),plt.gca().xaxis.get_ticklabels())]
  #plt.ylabel(')
  plt.gca().get_xaxis().set_ticks_position('none')
  plt.gca().get_yaxis().set_ticks_position('none')


def DrawWaterDelta(DF,DF2,Label,save_dir):
    mpl.rcParams['font.family'] = 'Arial Unicode MS'
    num_mol =len(Label)#//len(DataType)#分子数 ,plotStyledict,Titledict
    colmol = 2#1行あたりに書く分子数
    avestd = 1#平均波形のみ別に書くなら2それ以外は1
    num_cols=4#len(DataType)*colmol*avestd #1行あたりのグラフ数
    num_rows=7
    num_page = num_mol*2//(num_cols*num_rows)
    mod_mol = num_mol%(num_cols*num_rows)
    molcount=-1
    Label2=[]
    for i in range(len(Label)):
        Label2.append(Label[i])
        Label2.append(Label[i]) 
    Label = Label2#Labelは2分子ずつ並ぶ
    print(Label)
    for k in range(0, num_page +1):#ページ作る
        fig, host = plt.subplots(num_rows,num_cols,figsize=(30,40))
        plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05)
        for j in range(0,num_rows): #行番号,一行につき1被験者 
            for i in range(0,num_cols):#列で回す
                if molcount < len(Label):
                    molcount += 1
                    DF.loc[Label[molcount],'Incretin':] = [np.nan if type(list(DF.loc[Label[molcount]])[ii]) == str else list(DF.loc[Label[molcount]])[ii] for ii in range(1,len(DF.loc[Label[molcount]])) ]
                    DF2.loc[Label[molcount],:] = [np.nan if type(list(DF2.loc[Label[molcount]])[ii]) == str else list(DF2.loc[Label[molcount]])[ii] for ii in range(len(DF2.loc[Label[molcount]])) ]
                    
                    if i%2 == 0:#偶数なら
                        host[j,i].plot(list(DF.loc['Time','Incretin':]), list(DF.loc[Label[molcount],'Incretin':]),c='blue', marker="o",linestyle = 'solid',linewidth=2.0,markersize=10.0)
                        host[j,i].plot(list(DF2.loc['Time','Incretin':][1:]), list(DF2.loc[Label[molcount],'Incretin':][1:]),c='black', marker="o",linestyle = '-',linewidth=2.0,markersize=10.0)
                        host[j,i].set_ylabel( DF.loc[Label[molcount],'Unit'],fontsize=40)#濃度
                    else:#奇数なら
                        DFlen = len(list(DF.loc['Time','Incretin':])); DF2len = len(list(DF2.loc['Time','Incretin':]))
                        host[j,i].plot(list(DF.loc['Time','Incretin':]), list(DF.loc[Label[molcount],'Incretin':] - DF.loc[Label[molcount],'Incretin':][0]) ,c='blue', marker="o",linestyle = 'solid',linewidth=2.0,markersize=10.0)
                        host[j,i].plot(list(DF2.loc['Time','Incretin':][1:]), list(DF2.loc[Label[molcount],'Incretin':][1:] - DF2.loc[Label[molcount],'Incretin':][1] ) ,c='black', marker="o",linestyle = '-',linewidth=2.0,markersize=10.0)
                        host[j,i].set_ylabel('Delta',fontsize=40)#濃度

                    host[j,i].set_title(Label[molcount],loc='center')#BolusLabel[molcount]
                    plt.rcParams["font.size"] = 25
                    host[j,i].set_xlabel("time (min.)",fontsize=30)#時系列


        plot_axis = plt.axis()
        #host[3,2].axis('off')#いらないところはこれで消す
        #host[3,3].axis('off')
        if not os.path.isdir(save_dir + '/TimeCourseWater_Delta/'):
            os.makedirs(save_dir + '/TimeCourseWater_Delta/')
        plt.savefig(save_dir + '/TimeCourseWater_Delta/' + str(k) +'.pdf',bbox_inches="tight")
        plt.close()
        print(save_dir + '/TimeCourseWate_Deltar/' + str(k) +'.pdf')
        
        
def DrawWater(DF,DF2,Label,save_dir):
    mpl.rcParams['font.family'] = 'Arial Unicode MS'
    num_mol =len(Label)#//len(DataType)#分子数 ,plotStyledict,Titledict
    colmol = 2#1行あたりに書く分子数
    avestd = 1#平均波形のみ別に書くなら2それ以外は1
    num_cols=4#len(DataType)*colmol*avestd #1行あたりのグラフ数
    num_rows=7
    num_page = num_mol//(num_cols*num_rows)
    mod_mol = num_mol%(num_cols*num_rows)
    molcount=-1
    for k in range(0, num_page +1):#ページ作る
        fig, host = plt.subplots(num_rows,num_cols,figsize=(30,40))
        plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05)
        for j in range(0,num_rows): #行番号,一行につき1被験者 
            for i in range(0,num_cols):#列で回す
                if molcount < len(Label):
                    molcount += 1
                    DF.loc[Label[molcount],'Incretin':] = [np.nan if type(list(DF.loc[Label[molcount]])[ii]) == str else list(DF.loc[Label[molcount]])[ii] for ii in range(1,len(DF.loc[Label[molcount]])) ]
                    DF2.loc[Label[molcount],:] = [np.nan if type(list(DF2.loc[Label[molcount]])[ii]) == str else list(DF2.loc[Label[molcount]])[ii] for ii in range(len(DF2.loc[Label[molcount]])) ]
                    
                    host[j,i].plot(list(DF.loc['Time','Incretin':]), list(DF.loc[Label[molcount],'Incretin':]),c='blue', marker="o",linestyle = 'solid',linewidth=2.0,markersize=10.0)
                    host[j,i].plot(list(DF2.loc['Time','Incretin':]), list(DF2.loc[Label[molcount],'Incretin':]),c='black', marker="o",linestyle = '-',linewidth=2.0,markersize=10.0)
                    
                    host[j,i].set_title(Label[molcount],loc='center')#BolusLabel[molcount]
                    plt.rcParams["font.size"] = 25
                    host[j,i].set_ylabel( DF.loc[Label[molcount],'Unit'],fontsize=40)#濃度
                    host[j,i].set_xlabel("time (min.)",fontsize=30)#時系列


        plot_axis = plt.axis()
        #host[3,2].axis('off')#いらないところはこれで消す
        #host[3,3].axis('off')
        if not os.path.isdir(save_dir + '/TimeCourseWater/'):
            os.makedirs(save_dir + '/TimeCourseWater/')
        plt.savefig(save_dir + '/TimeCourseWater/' + str(k) +'.pdf',bbox_inches="tight")
        plt.close()
        print(save_dir + '/TimeCourseWater/' + str(k) +'.pdf')

#列同士の相関を網羅的に求めていく
def TmSrCorr(DF,SwitchDict,save_dir):
    #被験者×パラメタ_分子なら
    if SwitchDict['MolType'] == 'MolTimeProp':
        print('yes')
        #for j in range(len(TryList)): 
        if SwitchDict['method'] == 'pearson':            
            Corr,Pval = SC.calculate_CorrWpearson(DF); Name = '_pearson'
        elif SwitchDict['method'] == 'spearman':
            Corr,Pval = SC.calculate_CorrWspeaman(DF); Name = '_spearman'
        plt.rcParams["font.size"] = 10
        LastName = list(Corr.index)[ np.where(np.array([Corr.index=='BMI']))[1][0]-1]            
        Corr = Corr.loc[:LastName,'BMI':'TwoHInsulin'];Pval = Pval.loc[:LastName,'BMI':'TwoHInsulin'];
        CorrUpper = Corr; PvalUpper =  Pval
        Corr.to_excel(save_dir + 'Corr_' +'_W'+Name+'.xlsx'); Pval.to_excel(save_dir + '/Pvalue_' + '_W' + Name+ '.xlsx');
        plt.hist(Corr.values[~np.isnan(Corr.values)]);plt.title('_W' + Name ); plt.savefig(save_dir +'_DistOfCorr_W' + Name + '.pdf');plt.close()
        plt.hist(Pval.values[~np.isnan(Pval.values)]);plt.title( '_W' + Name); plt.savefig(save_dir +'_DistOfPval_W' + Name + '.pdf');plt.close()


    elif SwitchDict['MolType'] == 'Param':
        ColList =list(DF.columns) #[SwitchDict['Time'] + '_' + SwitchDict['EngLabel'][jj] for jj in range(len(SwitchDict['EngLabel']))]
        Corr,Pval,DFNewScat = APPHel.calccorrpval(DF,SwitchDict)
        plt.rcParams["font.size"] = 10
        DFNewScat.to_excel(save_dir+'DFNewScat_'+SwitchDict['method']+'.xlsx')
        Corr.to_excel(save_dir+ '/Corr_' +'_W'+SwitchDict['method']+'.xlsx'); Pval.to_excel(save_dir + '/Pvalue_' + '_W'+SwitchDict['method']+'.xlsx');
        CorrUpper = MH.adjustMatrlower(Corr);  PvalUpper = MH.adjustMatrlower(Pval) #対角と下三角をNaNにする
        Corr.to_excel(save_dir+'/Corr_'  + '_W'+SwitchDict['method']+'_Upper.xlsx')    ; Pval.to_excel(save_dir+'/Pvalue_'+'_W'+SwitchDict['method']+'_Upper.xlsx')        

        plt.figure();plt.hist(list(Corr.values[~pd.isnull(Corr.values)]));plt.title('_W'+SwitchDict['method']); plt.savefig(save_dir +'_DistOfCorr_'+SwitchDict['method']+'.pdf');plt.close()
        plt.figure();plt.hist(list(Pval.values[~pd.isnull(Pval.values)]));plt.title( '_W'+SwitchDict['method']); plt.savefig(save_dir +'_DistOfPval_W'+SwitchDict['method']+'.pdf');plt.close()
        
    ####被験者×各時点分子濃度なら
    elif SwitchDict['MolType'] == 'Different':
        #同じ時間の奴だけ取り出してその数だけCorrとか出せば良い
        if SwitchDict['Corr_Mol_Time']==1:#相関計算時に分子x時間で指定するなら 
            ColList = [SwitchDict['Time'] + '_' +  SwitchDict['EngLabel'][jj] for jj in range(len(SwitchDict['EngLabel']))]
        else:
            ColList=list(DF.columns)
        #DFCol = list(DF.columns)
        #MolList = [DF.loc[re.compile("(.*)(_)(.*)").search(DFCol[x]).group(1),'ClstColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
        DF=DF[ColList]
        Corr,Pval = SC.calculate_CorrWpearson(DF.loc[:,ColList])
        plt.rcParams["font.size"] = 10
        if not os.path.isdir(save_dir+SwitchDict['Time']):
          os.makedirs(save_dir+SwitchDict['Time'])  
        if SwitchDict['Corr_Mol_Time']==1:#相関計算時に分子x時間で指定するなら 

            NewColList = [SwitchDict['EngLabel'][jj] + '_' + SwitchDict['Time'] for jj in range(len(SwitchDict['EngLabel']))]
        else:
            NewColList =SwitchDict['EngLabel']
        Corr=Corr[ColList]
        Corr.columns=NewColList;Corr.index=NewColList;
        Corr.to_excel(save_dir +SwitchDict['Time']+ '/Corr_' +'_Wpearson.xlsx'); Pval.to_excel(save_dir +SwitchDict['Time']+ '/Pvalue_' + '_Wpearson.xlsx');
        AllCorr = Corr.copy(); AllPval = Pval.copy()
        CorrUpper = MH.adjustMatrlower(Corr);  PvalUpper = MH.adjustMatrlower(Pval) #対角と下三角をNaNにする
        Corr.to_excel(save_dir+SwitchDict['Time']+'/Corr_'  + '_Wpearson_Upper.xlsx')    ; Pval.to_excel(save_dir+SwitchDict['Time']+'/Pvalue_'+'_Wpearson_Upper.xlsx')        
#Upperじゃなくていいのか？
        plt.xlim([-1,1])
        plt.hist(CorrUpper.values[~np.isnan(Corr.values)]);plt.title('_Wpearson'); plt.savefig(save_dir +SwitchDict['Time']+'/_DistOfCorr_Wpearson.pdf');plt.close()
        plt.hist(PvalUpper.values[~np.isnan(Pval.values)]);plt.title( '_Wpearson'); plt.savefig(save_dir + SwitchDict['Time']+'/_DistOfPval_Wpearson.pdf');plt.close()

        
    elif SwitchDict['MolType'] == 'All':
        #for j in range(len(TryList)):             
        Corr,Pval = SC.calculate_CorrWpearson(DF)
        plt.rcParams["font.size"] = 10
        Corr.to_excel(save_dir + 'Corr_' +'_Wpeason.xlsx'); Pval.to_excel(save_dir + '/Pvalue_' + '_Wpeason.xlsx');
        CorrUpper = MH.adjustMatrlower(Corr);  PvalUpper = MH.adjustMatrlower(Pval) #対角と下三角をNaNにする
        CorrUpper.to_excel(save_dir+'/Corr_'  + '_Wpeason_Upper.xlsx')    ; PvalUpper.to_excel(save_dir+'/Pvalue_'+'_Wpeason_Upper.xlsx')        
        plt.hist(CorrUpper.values[~np.isnan(Corr.values)]);plt.title('_Wpeason'); plt.savefig(save_dir +'_DistOfCorr_Wpeason_Upper.pdf');plt.close()
        plt.hist(PvalUpper.values[~np.isnan(Pval.values)]);plt.title( '_Wpeason'); plt.savefig(save_dir +'_DistOfPval_Wpeason_Upper.pdf');plt.close()


    return(CorrUpper, PvalUpper,AllCorr,AllPval, save_dir)
    

    
def PlotScatter(SwitchDict,save_dir):
     #file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180806'#2分子濃度変動
     #file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180827'
     file_dir =save_dir
     
     #sp =ij;#Peason or speaman#####PArama?Param?  Corr_CutOffTimePropRepesen_Separate_Clst7_WPeason
     #Corr= pd.read_excel(file_dir+'/'+file_label+'/Corr_CutOffTimeProp'+Rep+'_Separate_'+file_label2+'W'+sp+'.xlsx',header=0,encoding = "ISO-8859-1");
     #Pvalue=pd.read_excel(file_dir+'/'+file_label+'/Pvalue_CutOffTimeProp'+Rep+'_Separate_'+file_label2+'W'+sp+'.xlsx',header=0,encoding = "ISO-8859-1");
     #Qvalue=pd.read_excel(file_dir+'/'+file_label+'/QvalueStoreyCutOff_Prop.xlsx',header=0,encoding = "ISO-8859-1")
     
     #Corr= pd.read_excel(file_dir+'/'+file_label+'/Corr_ParamACDElog10SubjTimeVar_Separate_WPeason.xlsx',header=0,encoding = "ISO-8859-1");
     Corr= pd.read_excel(file_dir+'/'+ 'Corr__W'+SwitchDict['method']+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0);
     Pvalue=pd.read_excel(file_dir+'/'+ 'Pvalue__W'+SwitchDict['method']+'.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
     
     Qvalue=pd.read_excel(file_dir+'/'+SwitchDict['EngLabel']+'/'+'QvalueStorey'+SwitchDict['EngLabel']+'_c.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
     #Qvalue=pd.read_excel(file_dir+'/QvalueBH_.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
     #Qvalue=pd.read_excel(file_dir+'/'+ SwitchDict['Time']+ '/QvalueStorey'+ SwitchDict['Time']+'_'+ SwitchDict['Time']+ '.xlsx',header=0,encoding = "ISO-8859-1")
     
     #Qvalue.index =list(Pvalue.index); Qvalue.columns = list(Pvalue.columns)
     #Qvalue.to_excel(file_dir+'/QvalueStorey_.xlsx')
     
     Qvalue.index =list(Pvalue.index); Qvalue.columns = list(Pvalue.columns)
     Qvalue.to_excel(file_dir+'/QvalueBH_.xlsx')
     
     #CutOffRep = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Progress_report/20180621_invivohumanmtg/ModelParam/代表分子/ParamaACDElog10Repesen.xlsx',header=0,encoding = "ISO-8859-1");
     
     #Corr= pd.read_excel(file_dir+'/'+file_label+'/Corr_ParamACDElog10SubjTimeVar_Separate_WPeason.xlsx',header=0,encoding = "ISO-8859-1");
     Corr.columns=list(Pvalue.columns);Corr.index=list(Pvalue.index) # Corrのラベルを修正する
     CombDF = APPHel.mkExcellower2(Corr,Pvalue,Qvalue,0.1)
#     if not os.path.isdir(save_dir +file_label ):
#      os.makedirs(save_dir + file_label)
     CombDF.to_excel(file_dir+ '/CombDF_CorrPvalQval_Param.xlsx')
     plt.figure();plt.hist(list(CombDF['Corr'].values[~pd.isnull(CombDF['Corr'].values)]),bins=15);plt.xticks(np.arange(-0.2, 0.69, 0.1))
     plt.title('Corr_W'+SwitchDict['method']+'_Q<0.1'); plt.savefig(save_dir +'Q<0.1_DistOfCorr_W'+SwitchDict['method']+'.pdf');plt.close()
     #DF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180706/時定数/CutOffProDF.xlsx',header=0,encoding = "ISO-8859-1")
     #DF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180827/DFNewScat_pearson.xlsx',header=0,encoding = "ISO-8859-1")
     #DF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180806/Score+Prop.xlsx')
     #DF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180802/DFNewScat.xlsx')
     #DF = pd.read_excel(file_dir+'/'+file_label+'/_'+file_label+'_ParamaACDECutOffTimeProp'+Rep+'.xlsx',header=0,encoding = "ISO-8859-1")
     DF=SwitchDict['Data']
     #ParamaACDElog10SubjTimeVar#ParamaACDElog10#OParamaACDElog10Prop
     APPHel.PlotScatterDF(DF,CombDF,save_dir ,3)#1で全部赤、3で#散布図：生を赤に、負を青に
     plt.close()

def TempAdjustNFD(OptionDict,axislabel,axislabelname,SubjectName,save_dir):#平均の平均でわる
        file_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/Raw_Eng/ケトン検出限界含む/'
        #i=0   
        #kanpachi = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=None) 
        #SubjTmCsDeltaPaneltemp = pd.Panel({SubjectName[0]:pd.DataFrame(data=kanpachi,index=list(kanpachi.index),columns=list(kanpachi.columns))})
        #SubjTmCsDeltaPanel = pd.Panel({SubjectName[0]:pd.DataFrame(data=kanpachi,index=list(kanpachi.index),columns=list(kanpachi.columns))})
        for i in range(len(list(SubjectName))):
            SubjTmCsDelta = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0)       
            #SubjTmCsDeltaPaneltemp[SubjectName[i]] = SubjTmCsDelta
            if i== 0: 
                kanpachi = pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
            else:
                kanpachi=pd.concat([kanpachi,SubjTmCsDelta],axis=0,join_axes=[kanpachi.columns]) 
        SubjTmCsDeltaPaneltemp = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))

        for i in range(len(list(SubjectName))):#ここでRawを元に正規化する
            TempDF = SubjTmCsDeltaPaneltemp.loc[SubjectName[i]] / SubjTmCsDeltaPaneltemp.mean(level=1).mean()     
            if OptionDict['DataFlag'] == 'Delta':
                TempDF = TempDF - TempDF.loc[0]            
            if OptionDict['DataFlag'] == 'FC':
                TempDF = TempDF / TempDF.loc[0]             
            if i== 0: 
                kanpachi = TempDF
            else:
                kanpachi=pd.concat([kanpachi,TempDF],axis=0,join_axes=[kanpachi.columns]) 
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDelta.index)]),columns=list(SubjTmCsDelta.columns))
                 
        if OptionDict['DataFlag'] == 'Raw':
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean   
            OptionDict['DataFlag'] = 'Raw_NormalizationFirst_'+axislabelname
        elif OptionDict['DataFlag'] == 'Delta':    
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean
            if axislabel ==3:
                SubjCV = SubjCV.drop(0)
                SubjCV = np.log10(SubjCV)
            OptionDict['DataFlag'] = 'Delta_NormalizationFirst_'+axislabelname            
        elif OptionDict['DataFlag'] =='FC':    
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean   
            if axislabel ==3:
                SubjCV = SubjCV.drop(0)
                SubjCV = np.log10(SubjCV)
            OptionDict['DataFlag'] = 'FC_NormalizationFirst_'+axislabelname 
        elif OptionDict['DataFlag'] =='RawDeltaFC':    
            Subjmean = SubjTmCsDeltaPanel.mean(level=1)
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean 

            if axislabel ==3:
                SubjCV = SubjCV.drop(0)
                SubjCV = np.log10(SubjCV)
        elif OptionDict['DataFlag'] == 'DeltaStd':    
            Subjmean = SubjTmCsDeltaPanel.mean(level=1) / SubjTmCsDeltaPanel.std(level=1,ddof=0).mean()
            Subjmean_Z = Subjmean.apply(zscore, axis=0)
            SubjCV = SubjTmCsDeltaPanel.std(level=1,ddof=0) #/ Subjmean
            if axislabel ==3:
                SubjCV = SubjCV.drop(0)
                SubjCV = np.log10(SubjCV)
            OptionDict['DataFlag'] = 'Delta_NormalizationFirst_'     +axislabelname            
class AnalSLE:
    def __init__(self):
        self.label=[]
        self.subjectName=[]
        self.optiondict=dict()
        self.timepointlist=[]
        self.save_dir=''
        
    def Anal(self,DF):
            if self.optiondict['Target'] =='FoldChnage':
                #self.AnalFoldChangeSLE(DF)#FCに対する各時点の相関解析 
                #self.AnalFCSLEIndivConcat(DF)#FCに対する定常（t=-10,0）の全被験者に対するある被験者
                pass
            else:
                #self.AnalRawSLE(DF)#Rawに対する各時点の相関解析 
                #self.AnalRawSLEIndiv(DF)#Rawに対する各時点の相関解析 
                #self.AnalRawSLEIndivConcat(DF)#Rawに対する定常（t=-10,0）の全被験者に対するある被験者
                pass
            self.DrawSLETimeCourse()
            self.DrawSLETimeCourseMol()#各個人各分子SLEの時間変化を描画
    def DrawSLETimeCourseMol(self):#各個人各分子SLEの時間変化を描画
        (num_rows,num_cols)=(1,1)

        
        for s in self.subjectName:#各被験者ごとにH(k,t)_delta, H(t)_deltaを得る
            fig, host = plt.subplots(num_rows,num_cols,figsize=(11.25,7.5))#20,10))
            plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05);molcount=0
            count=0
            for k in self.label:#ある被験者の全分子
                TmCsDF=pd.read_excel(self.save_dir+s+'_H_kt_delta_DF.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
                Col=list(TmCsDF.columns)#Time
                
                host.plot(Col,list(TmCsDF.loc[k]), color=cm.hsv(count/len(self.label)), marker="o",linestyle = 'solid',linewidth=2.0,markersize=5.0)
                #host.legend(self.label,prop={'size':20},bbox_to_anchor=(0.95, 0.9))
                                
                    #host.set_xticks([0,60,120,180,240])
                    #host.set_xticklabels(['0','60','120','180','240'])
                host.annotate(k,fontsize=1, xy=(Col[-1],list(TmCsDF.loc[k])[-1]))
                host.set_ylim([-0.05,1.0])
        
                plot_axis = plt.axis()
                #host[3,2].axis('off')#いらないところはこれで消す
                #host[3,3].axis('off')
                try:
                    plt.title(Optiondict['Title'])
                except:
                    pass
                count+=1
            plt.rcParams["font.size"] = 25
            host.set_ylabel('',fontsize=40)#濃度
            xmin, xmax, ymin, ymax = plt.axis()
            #plt.hlines(0.6, xmin, xmax, "black", linestyles='dashed',linewidth=1)
            #plt.hlines(-0.6, xmin, xmax, "black", linestyles='dashed',linewidth=1)
    
            #host.set_ylabel('CV',fontsize=40)#濃度
            host.set_xlabel("Time (min.)",fontsize=30)#時系列
            host.tick_params(labelsize=15, axis='both')
            host.set_xticks(list(TmCsDF.columns))
            host.set_xticklabels(list(TmCsDF.columns)) 
            plt.savefig(self.save_dir + s+'_H_delta_k_TimeCourse.pdf',bbox_inches="tight")
            plt.savefig(self.save_dir + s+'_H_delta_k_TimeCourse.png',bbox_inches="tight")
            plt.close()
    def DrawSLETimeCourse(self):#各個人の平均SLEの時間変化を描画
        (num_rows,num_cols)=(1,1)
        fig, host = plt.subplots(num_rows,num_cols,figsize=(11.25,7.5))#20,10))
        plt.subplots_adjust(wspace=0.4, hspace=0.6,right=0.85,left=0.05);molcount=0
        count=0
        for s in self.subjectName:#各被験者ごとにH(k,t)_delta, H(t)_deltaを得る
            TmCsDF=pd.read_excel(self.save_dir+s+'_H_t_delta_DF.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
            Col=list(TmCsDF.columns)#Time
            
            host.plot(Col,list(TmCsDF.iloc[0]), color=cm.hsv(count/len(self.subjectName)), marker="o",linestyle = 'solid',linewidth=2.0,markersize=5.0)
            host.set_xticks(list(TmCsDF.columns))
            host.set_xticklabels(list(TmCsDF.columns))            
                            
                #host.set_xticks([0,60,120,180,240])
                #host.set_xticklabels(['0','60','120','180','240'])
            host.annotate(s,fontsize=1, xy=(Col[-1],list(TmCsDF.iloc[0])[-1]))

            plt.rcParams["font.size"] = 25
            host.set_ylabel('',fontsize=40)#濃度
            xmin, xmax, ymin, ymax = plt.axis()
            #plt.hlines(0.6, xmin, xmax, "black", linestyles='dashed',linewidth=1)
            #plt.hlines(-0.6, xmin, xmax, "black", linestyles='dashed',linewidth=1)
    
            #host.set_ylabel('CV',fontsize=40)#濃度
            host.set_xlabel("Time (min.)",fontsize=30)#時系列
            host.tick_params(labelsize=15, axis='both')
    
            plot_axis = plt.axis()
            #host[3,2].axis('off')#いらないところはこれで消す
            #host[3,3].axis('off')
            try:
                plt.title(Optiondict['Title'])
            except:
                pass
            count+=1
        plt.savefig(self.save_dir + 'H_delta_TimeCourse.pdf',bbox_inches="tight")
        plt.savefig(self.save_dir + 'H_delta_TimeCourse.png',bbox_inches="tight")
        plt.close()
        
    def AnalFoldChangeSLE(self,DF):
        pass
        
        
    def AnalRawSLE(self,DF):#Rawに対するSLE解析
        for s in self.subjectName:#各被験者ごとにH(k,t)_delta, H(t)_deltaを得る
            H_kt_delta_DF = pd.DataFrame(data=None,index=self.label,columns=self.timepointlist)
            H_t_delta_DF = pd.DataFrame(data=None,index=[0],columns=self.timepointlist)

        #被験者19人で作成して、1人ずつ当てはめるようにする
            for k in self.label:
                for t in self.timepointlist:
                    self.optiondict['EngLabel'] = self.label;
                    self.optiondict['Time'] = str(t)
                    ColList = [self.optiondict['Time'] + '_' +  self.optiondict['EngLabel'][jj] for jj in range(len(self.optiondict['EngLabel']))]
                    tempDF = DF[ColList]
                    #tempDF.columns=self.label
                    #ある分子kを中心にした時の19人のPCC、Stdを算出
                    tempDF_19 = tempDF.drop(s,axis=0)
                    P_itDF_19,Std_19,Corr_19 = self.calcCorr(tempDF_19,k)
                    #エントロピーの算出
                    H_kt_19 = self.calcH(P_itDF_19)
                    #1人追加した時の20人のPCC,Stdを算出
                    P_itDF_20,Std_20,Corr_20 = self.calcCorr(tempDF,k)

                    #1人追加した時のエントロピー算出
                    H_kt_20 = self.calcH(P_itDF_20)

                    #1人追加した時の差分エントロピーの算出
                    Std_delta = np.abs(Std_20 - Std_19)
                    H_kt_delta = Std_delta * (np.abs(H_kt_20 - H_kt_19))
                    #差分エントロピーの平均の算出
                    H_kt_delta_DF.loc[:,t]=list(H_kt_delta)
                    print(s+'_'+k+'_'+str(t))
                    #ある分子kのある被験者sを抜いた時のCorr
                    Corr_19.to_excel(self.save_dir+s+'_'+k+'_'+str(t)+'_Corr19DF.xlsx')
                    Std_19.to_excel(self.save_dir+s+'_'+k+'_'+str(t)+'_Std19DF.xlsx')
            H_t_delta_DF.loc[0]=list(H_kt_delta_DF.mean(axis=0))
            H_kt_delta_DF.to_excel(self.save_dir+s+'_H_kt_delta_DF.xlsx')
            H_t_delta_DF.to_excel(self.save_dir+s+'_H_t_delta_DF.xlsx')
    def AnalFCSLEIndivConcat(self,DF):#Rawに対するSLE解析20人でreference作るver.
        #-10minと0minの平均した空腹値DFを作成
        FastingDF=self.AdjustFasting(DF)
        # Reference作る
        #たて、分子、横、4時点分？
        self.optiondict['EngLabel'] = self.label;

        ColListformer=[]
        Colm = []
        for p in ['-10','0']:
            for o in self.subjectName:
                Colm += [o+p]
        tempDF_Ref= pd.DataFrame(data=None,index=self.label,columns=Colm)#分子 x 時間_被験者

        for jj in self.optiondict['EngLabel']:
            #for kk in ['-10','0']:#初めの2点
            #ColListformer = [kk + '_' +  self.optiondict['EngLabel'][jj] for jj in range(len(self.optiondict['EngLabel']))]
            tempFastPreDF=DF['-10_' + jj]
            tempFastPreDF.columns=jj
            FastPreDF=tempFastPreDF / FastingDF[jj]#ラベルを分子に戻し、空腹値で割った
            
            #FastPreDF.columns='-10_' + jj
            tempFastPostDF=DF['0_' + jj]
            tempFastPostDF.columns=jj            
            FastPostDF  = tempFastPostDF / FastingDF[jj]#ラベルを分子に戻し、空腹値で割った
            #FastPostDF.columns='0_' + jj#そしてまたラベルを変える

            tempDF_Ref.loc[jj]=list( FastPreDF )+list(FastPostDF )#ある分子の行に時間_被験者ぶん入れる

##############b     ここから

        tempDF=tempDF_Ref.copy()
        
        for s in self.subjectName[19:20]:#各被験者ごとにH(k,t)_delta, H(t)_deltaを得る
            H_kt_delta_DF = pd.DataFrame(data=None,index=self.label,columns=self.timepointlist)
            H_t_delta_DF = pd.DataFrame(data=None,index=[0],columns=self.timepointlist)

        #被験者19人で作成して、1人ずつ当てはめるようにする
            for k in self.label:
                self.optiondict['EngLabel'] = self.label;

                for t in self.timepointlist:
                    self.optiondict['Time'] = str(t)
                    ColList = [self.optiondict['Time'] + '_' +  self.optiondict['EngLabel'][jj] for jj in range(len(self.optiondict['EngLabel']))]
                    tempDF[s+'_'+str(t)] = [ (DF.loc[s][ColList])[ik] / (FastingDF.loc[s][list(self.optiondict['EngLabel'])])[ik] for ik in range(len(ColList)) ]#Referenceにサンプルを足す
                    #tempDF.columns=self.label
                    #ある分子kを中心にした時の19人のPCC、Stdを算出
                    #tempDF_19 = tempDF.drop(s,axis=0)
                    P_itDF_Ref,Std_Ref,Corr_Ref = self.calcCorr(tempDF_Ref.T,k)
                    #エントロピーの算出
                    H_kt_Ref = self.calcH(P_itDF_Ref)
                    #1人追加した時の20人のPCC,Stdを算出
                    P_itDF_20,Std_20,Corr_20 = self.calcCorr(tempDF.T,k)

                    #1人追加した時のエントロピー算出
                    H_kt_20 = self.calcH(P_itDF_20)

                    #1人追加した時の差分エントロピーの算出
                    Std_delta = np.abs(Std_20 - Std_Ref)
                    H_kt_delta = Std_delta * np.abs(H_kt_20 - H_kt_Ref)
                    #差分エントロピーの平均の算出
                    H_kt_delta_DF.loc[:,t]=list(H_kt_delta)
                    print(s+'_'+k+'_'+str(t))
                    tempDF = tempDF.drop(s+'_'+str(t),axis=1)#Referenceから足したのを引く
                #ある分子kのある被験者の、初めの4点+1点を使った時のCorr
                    
                    Corr_20.to_excel(self.save_dir+s+'_'+k+'_'+str(t)+'_CorrDF.xlsx')
                    Std_20.to_excel(self.save_dir+s+'_'+k+'_'+str(t)+'_StdDF.xlsx')
                #ある分子kのある被験者の、初めの4点を使った時のCorr
                Corr_Ref.to_excel(self.save_dir+s+'_'+k+'_CorrRefDF.xlsx')
                Std_Ref.to_excel(self.save_dir+s+'_'+k+'_StdRefDF.xlsx')
            H_t_delta_DF.loc[0]=list(H_kt_delta_DF.mean(axis=0))
            H_kt_delta_DF.to_excel(self.save_dir+s+'_H_kt_delta_DF.xlsx')
            H_t_delta_DF.to_excel(self.save_dir+s+'_H_t_delta_DF.xlsx') 
            
    def AnalRawSLEIndivConcat(self,DF):#Rawに対するSLE解析20人でreference作るver.
        # Reference作る
        #たて、分子、横、4時点分？
        self.optiondict['EngLabel'] = self.label;

        ColListformer=[]
        Colm = []
        for p in ['-10','0']:
            for o in self.subjectName:
                Colm += [o+p]
        tempDF_Ref= pd.DataFrame(data=None,index=self.label,columns=Colm)#分子 x 時間_被験者

        for jj in self.optiondict['EngLabel']:
            #for kk in ['-10','0']:#初めの2点
            #ColListformer = [kk + '_' +  self.optiondict['EngLabel'][jj] for jj in range(len(self.optiondict['EngLabel']))]
            
            tempDF_Ref.loc[jj]=list(DF['-10_' + jj])+list(DF['0_' + jj] )#ある分子の行に時間_被験者ぶん入れる
        tempDF=tempDF_Ref.copy()
        
        for s in self.subjectName[19:20]:#各被験者ごとにH(k,t)_delta, H(t)_deltaを得る
            H_kt_delta_DF = pd.DataFrame(data=None,index=self.label,columns=self.timepointlist)
            H_t_delta_DF = pd.DataFrame(data=None,index=[0],columns=self.timepointlist)

        #被験者19人で作成して、1人ずつ当てはめるようにする
            for k in self.label:
                self.optiondict['EngLabel'] = self.label;


                for t in self.timepointlist:
                    self.optiondict['Time'] = str(t)
                    ColList = [self.optiondict['Time'] + '_' +  self.optiondict['EngLabel'][jj] for jj in range(len(self.optiondict['EngLabel']))]
                    tempDF[s+'_'+str(t)] = list(DF.loc[s][ColList])#Referenceにサンプルを足す
                    #tempDF.columns=self.label
                    #ある分子kを中心にした時の19人のPCC、Stdを算出
                    #tempDF_19 = tempDF.drop(s,axis=0)
                    P_itDF_Ref,Std_Ref,Corr_Ref = self.calcCorr(tempDF_Ref.T,k)
                    #エントロピーの算出
                    H_kt_Ref = self.calcH(P_itDF_Ref)
                    #1人追加した時の20人のPCC,Stdを算出
                    P_itDF_20,Std_20,Corr_20 = self.calcCorr(tempDF.T,k)

                    #1人追加した時のエントロピー算出
                    H_kt_20 = self.calcH(P_itDF_20)

                    #1人追加した時の差分エントロピーの算出
                    Std_delta = np.abs(Std_20 - Std_Ref)
                    H_kt_delta = Std_delta * np.abs(H_kt_20 - H_kt_Ref)
                    #差分エントロピーの平均の算出
                    H_kt_delta_DF.loc[:,t]=list(H_kt_delta)
                    print(s+'_'+k+'_'+str(t))
                    tempDF = tempDF.drop(s+'_'+str(t),axis=1)#Referenceから足したのを引く
                #ある分子kのある被験者の、初めの4点+1点を使った時のCorr
                    
                    Corr_20.to_excel(self.save_dir+s+'_'+k+'_'+str(t)+'_CorrDF.xlsx')
                    Std_20.to_excel(self.save_dir+s+'_'+k+'_'+str(t)+'_StdDF.xlsx')
                #ある分子kのある被験者の、初めの4点を使った時のCorr
                Corr_Ref.to_excel(self.save_dir+s+'_'+k+'_CorrRefDF.xlsx')
                Std_Ref.to_excel(self.save_dir+s+'_'+k+'_StdRefDF.xlsx')
            H_t_delta_DF.loc[0]=list(H_kt_delta_DF.mean(axis=0))
            H_kt_delta_DF.to_excel(self.save_dir+s+'_H_kt_delta_DF.xlsx')
            H_t_delta_DF.to_excel(self.save_dir+s+'_H_t_delta_DF.xlsx') 
                    
    def AnalRawSLEIndiv(self,DF):#Rawに対するSLE解析個々人でreference作るver.
        for s in self.subjectName:#各被験者ごとにH(k,t)_delta, H(t)_deltaを得る
            H_kt_delta_DF = pd.DataFrame(data=None,index=self.label,columns=self.timepointlist)
            H_t_delta_DF = pd.DataFrame(data=None,index=[0],columns=self.timepointlist)

        #被験者19人で作成して、1人ずつ当てはめるようにする
            for k in self.label:
                self.optiondict['EngLabel'] = self.label;

                # Reference作る
                #たて、分子、横、4時点分？
                ColListformer=[]
                tempDF_Ref= pd.DataFrame(data=None,index=self.label,columns=['-10','0','10','20'])

                for kk in ['-10','0','10','20']:#初めの4点
                    ColListformer = [kk + '_' +  self.optiondict['EngLabel'][jj] for jj in range(len(self.optiondict['EngLabel']))]

                    tempDF_Ref[kk]=list(DF.loc[s][ColListformer])
                    tempDF=tempDF_Ref.copy()
                for t in self.timepointlist:
                    self.optiondict['Time'] = str(t)
                    ColList = [self.optiondict['Time'] + '_' +  self.optiondict['EngLabel'][jj] for jj in range(len(self.optiondict['EngLabel']))]
                    tempDF[t] = list(DF.loc[s][ColList])#Referenceにサンプルを足す
                    #tempDF.columns=self.label
                    #ある分子kを中心にした時の19人のPCC、Stdを算出
                    #tempDF_19 = tempDF.drop(s,axis=0)
                    P_itDF_Ref,Std_Ref,Corr_Ref = self.calcCorr(tempDF_Ref.T,k)
                    #エントロピーの算出
                    H_kt_Ref = self.calcH(P_itDF_Ref)
                    #1人追加した時の20人のPCC,Stdを算出
                    P_itDF_20,Std_20,Corr_20 = self.calcCorr(tempDF.T,k)

                    #1人追加した時のエントロピー算出
                    H_kt_20 = self.calcH(P_itDF_20)

                    #1人追加した時の差分エントロピーの算出
                    Std_delta = np.abs(Std_20 - Std_Ref)
                    H_kt_delta = Std_delta * np.abs(H_kt_20 - H_kt_Ref)
                    #差分エントロピーの平均の算出
                    H_kt_delta_DF.loc[:,t]=list(H_kt_delta)
                    print(s+'_'+k+'_'+str(t))
                #ある分子kのある被験者の、初めの4点+1点を使った時のCorr
                    
                    Corr_20.to_excel(self.save_dir+s+'_'+k+'_'+str(t)+'_CorrDF.xlsx')
                    Std_20.to_excel(self.save_dir+s+'_'+k+'_'+str(t)+'_StdDF.xlsx')
                #ある分子kのある被験者の、初めの4点を使った時のCorr
                Corr_Ref.to_excel(self.save_dir+s+'_'+k+'_CorrRefDF.xlsx')
                Std_Ref.to_excel(self.save_dir+s+'_'+k+'_StdRefDF.xlsx')
            H_t_delta_DF.loc[0]=list(H_kt_delta_DF.mean(axis=0))
            H_kt_delta_DF.to_excel(self.save_dir+s+'_H_kt_delta_DF.xlsx')
            H_t_delta_DF.to_excel(self.save_dir+s+'_H_t_delta_DF.xlsx')                    
    def calcCorr(self,DF,k):#83!/(83-2)!2!=3403通りのPCCに対する、分子k-i間のPCCがP_i
        templabel=self.label.copy()
        templabel.remove(k)

        P_itDF=pd.DataFrame(data=None,index=[0],columns=templabel)
        CorrUpper, PvalueUpper, Corr, Pval, self.save_dir = TmSrCorr(DF,self.optiondict,self.save_dir)
        CorrAbs = np.abs(Corr)#絶対値にする
        if self.optiondict['Corr_Mol_Time']==1:#相関計算時に分子x時間で指定するなら 

            CorrD = CorrAbs.drop(k+ '_' + self.optiondict['Time'],axis=1)
            P_itDF.iloc[0] = [list(CorrD.loc[k+ '_' + self.optiondict['Time']])[ii] / np.sum( list(CorrD.loc[k+ '_' + self.optiondict['Time']]) ) for ii in range(len(templabel))]
        else:
            CorrD = CorrAbs.drop(k,axis=1)
            P_itDF.iloc[0] = [list(CorrD.loc[k])[ii] / np.sum( list(CorrD.loc[k]) ) for ii in range(len(templabel))]
            
        Std=DF.std(axis=0)
        Std.to_excel(self.save_dir+k+'_StdDF.xlsx')

        return(P_itDF,Std,Corr)

    def calcH(self,P_itDF):#エントロピー算出
        M=len(P_itDF.columns)#82? 
        H_kt = -1 / np.log2(M) * np.sum( [P_itDF[jj] * np.log2(P_itDF.astype(np.float64))[jj] for jj in list(P_itDF.columns)] )
        return(H_kt)
    def AdjustFasting(self,DF):
            ColList0 = [ '0_' +  self.label[jj] for jj in range(len(self.label))]
            ColList10 = ['-10_' +  self.label[jj] for jj in range(len(self.label))]
            AllSubjTmCsZero = DF[ColList0]; AllSubjTmCsminus = DF[ColList10]; 
            AllSubjTmCsZero.columns=self.label;AllSubjTmCsminus.columns=self.label
            FastingDF = pd.DataFrame(data=None,index=self.subjectName,columns=self.label)
            for i in self.label:
                FastingDF[i] = np.nanmean([np.array(AllSubjTmCsminus[i]),np.array(AllSubjTmCsZero[i])],axis=0)
            return(FastingDF)        
        
class TmPtCorr:
    def __init__(self):
        self.label=[]
        self.subjectName=[]
        self.optiondict=dict()
        self.timepointlist=[]
        self.save_dir=''
    
    def AnalTmPtCorrEachMol(self,DF):#ある分子内で空腹値vsその他の時点の相関
            if self.optiondict['Target'] =='FoldChnage':
                self.AnalFoldChangeEachMol(DF)
            else:
                self.AnalRawEachMol(DF)#Rawに対する各時点の相関解析 
                
    def mkContinuousDF(self,BolusDF):#Continuousの各時点x各分子行列作る
        SubjectName=['isaki','iwashi','karei','shimaaji','unagi']
        for i in range(0,len(SubjectName)):   #各被験者で特徴量を算出する
            SubjTmCsDf = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Continuous/TimeCourse/English/Raw/New/'+SubjectName[i]+'_Continuous2h_Raw.xlsx',header=0,encoding = "ISO-8859-1",index_col=0).drop('time(min)',axis=1)
            EngLabel = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1")['English'])
            SubjTmCsDf.columns=EngLabel;SubjTmCsDf.to_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Continuous/TimeCourse/English/Raw/New/'+SubjectName[i]+'_Continuous2h_Raw_20191223.xlsx')
            if len(SubjTmCsDf.index) == 14:#-10minあるなら
                FastingDF = np.nansum(SubjTmCsDf.iloc[[0,1]],axis=0) / 2#空腹値を計算する
                SubjTmCsDf = SubjTmCsDf.drop(-10,axis=0)
                SubjTmCsDf.iloc[0,:]  =  FastingDF
                if len(SubjTmCsDf.columns) == 84:#3-HB消えてないなら
                    SubjTmCsDf = SubjTmCsDf.drop('Total bile acid',axis=1)
                EngLabel = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1")['English'])
                SubjTmCsDf.columns = EngLabel
                NewDF = pd.DataFrame(data=None, index=list( SubjTmCsDf.index ), columns = ['time(min)'] + list( SubjTmCsDf.columns )[0:-1])
                NewDF=SubjTmCsDf.copy()
                NewDF['time(min)']=list( SubjTmCsDf.index )
                NewDF.to_excel(save_dir+'SubjTmCs_' + SubjectName[i] + 'Raw.xlsx')

        #NewDF=SubjTmCsDf.copy()
        #NewDF['time(min)']=list( SubjTmCsDf.index )
        
        #
            #SubjTmCsDf=NewDF.copy()

            if i== 0: 
                kanpachi = SubjTmCsDf#pd.read_excel(file_dir + 'SubjTmCs_' + SubjectName[i] +'Raw.xlsx',header=0,encoding = "ISO-8859-1", index_col=0) 
            else:
                kanpachi=pd.concat([kanpachi,SubjTmCsDf],axis=0,join_axes=[kanpachi.columns])
    
                #SubjTmCsDeltaPanel[SubjectName[i]] = SubjTmCsDelta
        SubjTmCsDeltaPanel = pd.DataFrame(data=np.array(kanpachi),index=pd.MultiIndex.from_product([SubjectName, list(SubjTmCsDf.index)]),columns=list(SubjTmCsDf.columns))
        TimeList = list(SubjTmCsDf.index);MolList = list(SubjTmCsDf.columns);ColList=[ str(i) + '_' + j for j in MolList for i in [-10] + TimeList ]
        NewDF = pd.DataFrame(data=None,index=SubjectName,columns=ColList)
### -10minの取り扱いはあとで修正する!!! 20191223
        for jj in MolList:#分子で回す
            for ii in [-10] + TimeList :#時点で回す
                if  ii == -10:
                    NewDF[str(ii) + '_' + jj] = list(SubjTmCsDeltaPanel.swaplevel().loc[0][jj])
                else:
                    NewDF[str(ii) + '_' + jj] = list(SubjTmCsDeltaPanel.swaplevel().loc[ii][jj])
        NewDF.to_excel(self.save_dir + 'SubjTimeSeriesDFWnan_Continuous_Eng_Ketoneplus.xlsx')
        print('end')

    def AnalTmPtCorr(self,DF):
            if self.optiondict['Target'] =='FoldChnage':
                self.AnalFoldChange(DF)
            else:
                self.AnalRaw(DF)#Rawに対する各時点の相関解析 
    
    def AdjustFasting(self,DF):
            ColList0 = [ '0_' +  self.label[jj] for jj in range(len(self.label))]
            ColList10 = ['-10_' +  self.label[jj] for jj in range(len(self.label))]
            AllSubjTmCsZero = DF[ColList0]; AllSubjTmCsminus = DF[ColList10]; 
            AllSubjTmCsZero.columns=self.label;AllSubjTmCsminus.columns=self.label
            FastingDF = pd.DataFrame(data=None,index=self.subjectName,columns=self.label)
            for i in self.label:
                FastingDF[i] = np.nanmean([np.array(AllSubjTmCsminus[i].astype(float)),np.array(AllSubjTmCsZero[i].astype(float))],axis=0)
            return(FastingDF)
            
    def AnalFoldChange(self,DF):#FoldChangeに対する各時点の相関解析
            FastingDF=self.AdjustFasting(DF)
            for i in range(len(self.timepointlist)):
                self.optiondict['EngLabel'] = self.label;
                self.optiondict['Time'] = str(self.timepointlist[i])
                ColList = [self.optiondict['Time'] + '_' +  self.optiondict['EngLabel'][jj] for jj in range(len(self.optiondict['EngLabel']))]
                tempDF = DF[ColList]
                tempDF.columns=self.label
                DF =tempDF / FastingDF; DF.columns=ColList
                #時間だけ取り出したDFを渡す
                CorrUpper, PvalueUpper, Corr, Pval, self.save_dir = TmSrCorr(DF,self.optiondict,self.save_dir)
                mD.draw_heatmapWXY(Corr,self.save_dir+self.optiondict['Time']+'/',1,'MolColor','bwr')#bwr  Reds
                mHH.draw_heatmapWXY(Corr,self.save_dir+self.optiondict['Time']+'/',1,'MolColor',dict(),cmap='bwr')
                plt.close()
                self.optiondict['EngLabel']=self.optiondict['Time']
                SC.UseR(PvalueUpper,self.optiondict)
                self.optiondict['method']='pearson';self.optiondict['Data']=DF
                ASTH.PlotScatter(SwitchDict,self.save_dir+self.optiondict['Time']+'/');plt.close()
            count=99
            for jj in self.label:#分子ごとに76ほんの時系列
                count+=1
                tempIdx=self.label.copy();tempIdx.remove(jj)
                TmCsDF = pd.DataFrame(data=None,index=tempIdx,columns=self.timepointlist)
                for k in self.timepointlist:#時点で回す
                    CorrDF = pd.read_excel(self.save_dir+str(k)+'/Corr__Wpeason.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
                    CorrDF.index=self.label;CorrDF.columns=self.label
                    for m in tempIdx:#jj以外の分子
                        TmCsDF.loc[m,k]=CorrDF[m][jj]
### 0.6以上にするなら                        
                ThrList=list(set(list(TmCsDF.T[TmCsDF.T>0.6].dropna(how='all',axis=1).columns)+list(TmCsDF.T[TmCsDF.T<-0.6].dropna(how='all',axis=1).columns)))
                TmCsDF.loc[ThrList].to_excel(self.save_dir+jj+'.xlsx')
                Optiondict={'Color':LH.MolColor(tempIdx),#分子名を与えれば色のリストを吐く,
                        'Annotate' : 1,#分子名などをAnnotateするなら
                        'Label':tempIdx,
                        'Title':jj
                        }
                ASTH.plotTmCs(TmCsDF.T[tempIdx],Optiondict,self.save_dir+str(count)+'_')#CVのタイムコース受け取って1つのグラフに重ねて描画
            
    def AnalRaw(self,DF):#Rawに対する各時点の相関解析  
            for i in range(len(self.timepointlist)):
                AllSubjTmCs = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/SubjTimeSeriesDFWnan_Eng_Ketoneplus.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)#最新のを使う
    
                self.optiondict['EngLabel'] = self.label;
                self.optiondict['Time'] = str(self.timepointlist[i])
                CorrUpper, PvalueUpper, Corr, Pval, self.save_dir = TmSrCorr(DF,self.optiondict,self.save_dir)
                mD.draw_heatmapWXY(Corr,self.save_dir+self.optiondict['Time']+'/',1,'MolColor','bwr')#bwr  Reds
                mHH.draw_heatmapWXY(Corr,self.save_dir+self.optiondict['Time']+'/',1,'MolColor',dict(),cmap='bwr')
                plt.close()
                self.optiondict['EngLabel']=self.optiondict['Time']
                SC.UseR(PvalueUpper,Sself.optiondict)
                self.optiondict['method']='pearson';self.optiondict['Data']=DF
                ASTH.PlotScatter(self.optiondict,self.save_dir+self.optiondict['Time']+'/');plt.close()
            count=99
            for jj in self.label:#分子ごとに76ほんの時系列
                count+=1
                tempIdx=self.label.copy();tempIdx.remove(jj)
                TmCsDF = pd.DataFrame(data=None,index=tempIdx,columns=self.timepointlist)
                for k in self.timepointlist:#時点で回す
                    CorrDF = pd.read_excel(self.save_dir+str(k)+'/Corr__Wpearson.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
                    CorrDF.index=self.label;CorrDF.columns=self.label
                    for m in tempIdx:#jj以外の分子
                        TmCsDF.loc[m,k]=CorrDF[m][jj]
                ThrList=list(set(list(TmCsDF.T[TmCsDF.T>0.6].dropna(how='all',axis=1).columns)+list(TmCsDF.T[TmCsDF.T<-0.6].dropna(how='all',axis=1).columns)))
                TmCsDF.loc[ThrList].to_excel(self.save_dir+jj+'.xlsx')
                Optiondict={'Color':LH.MolColor(tempIdx),#分子名を与えれば色のリストを吐く,
                        'Annotate' : 1,#分子名などをAnnotateするなら
                        'Label':tempIdx,
                        'Title':jj
                        }
                ASTH.plotTmCs(TmCsDF.T[tempIdx],Optiondict,self.save_dir+str(count)+'_')#CVのタイムコース受け取って1つのグラフに重ねて描画

    def AnalRawEachMol(self,DF):#Rawに対する各時点の相関解析  
            FastingDF=self.AdjustFasting(DF)
            timepointlist=self.timepointlist;
            timepointlist.remove(0);
            timepointlist.remove(-10)
            print(timepointlist)
            MolCorrDF = pd.DataFrame(data=None,index=self.label,columns=timepointlist)
            MolPvalDF = pd.DataFrame(data=None,index=self.label,columns=timepointlist)
            try:
                MolCorrDF = pd.read_excel(self.save_dir+'Corr_Wpearson.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
                MolPvalDF = pd.read_excel(self.save_dir+'Pvalue_Wpearson.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
            except:
                print('here I am')
                for jj in self.label:#分子ごとに76ほんの時系列
                    #MolCorrDF[0][jj] =1 #便宜的にtime=0にCorr=1を代入
                    for ii in timepointlist:
                        #AllSubjTmCs = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/SubjTimeSeriesDFWnan_Eng_Ketoneplus.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)#最新のを使う
                        list1 = list(FastingDF[jj]); list2 = list(DF[str(ii)+'_'+jj])
                        r,p = SC.calcpeasonr(list1,list2)#SC.calcspearmanr(list1,list2)##
                        MolCorrDF[ii][jj] = r; MolPvalDF[ii][jj] = p
                        #print(str(ii)+'_'+jj)

            MolCorrDF.to_excel(self.save_dir+'/Corr_Wspearman.xlsx');  MolPvalDF.to_excel(self.save_dir+'/Pvalue_Wspearman.xlsx')      
            plt.figure();plt.hist(list(MolCorrDF.values[~pd.isnull(MolCorrDF.values)]),bins=15);#plt.xticks(np.arange(-0.2, 0.69, 0.1))
            plt.title('Corr_Wspearman'); plt.savefig(self.save_dir +'Q<0.1_DistOfCorr_Wspearman.pdf');plt.close()
            plt.figure();plt.hist(list(MolPvalDF.values[~pd.isnull(MolPvalDF.values)]),bins=15);#plt.xticks(np.arange(-0.2, 0.69, 0.1))
            plt.title('Pvalue_Wspearman'); plt.savefig(self.save_dir +'Q<0.1_DistOfPvalue_Wspearman.pdf');plt.close()
            
            self.optiondict['EngLabel']='RawEachMol'
            SC.UseR(MolPvalDF,self.optiondict)
            self.optiondict['method']='spearman';self.optiondict['Data']=DF
            #PlotScatter(self.optiondict,self.save_dir);plt.close()
            if any(self.optiondict['SignLabel'])==1:MolCorrDF=MolCorrDF.T[self.optiondict['SignLabel']].T;self.label=self.optiondict['SignLabel']
            Optiondict={'Color':LH.MolColor(self.label),#分子名を与えれば色のリストを吐く,
                    'Annotate' : 1,#分子名などをAnnotateするなら
                    'Label':self.label,
                    'Title':'Correlation vs Fasting'
                    }
            plotTmCs(MolCorrDF.T,Optiondict,self.save_dir)#CVのタイムコース受け取って1つのグラフに重ねて描画
            ClstMol,Colordict,ClstDF,ClstColorDF,ColorDF = mD.draw_heatmap(MolCorrDF,1,'MolColor',{'title':'Correlation vs Fasting'},self.save_dir+'Correlation vs Fasting_',cmap='PuOr_r')
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
                for jj in self.label:#分子ごとに76ほんの時系列
                    #MolCorrDF[0][jj] =1 #便宜的にtime=0にCorr=1を代入
                    for ii in self.timepointlist:
                        #AllSubjTmCs = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Timecourse/SubjTimeSeriesDFWnan_Eng_Ketoneplus.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)#最新のを使う
                        list1 = list(FastingDF[jj]); list2 = list(AllSubjTmCs[str(ii)+'_'+jj])
                        r,p = SC.calcpeasonr(list1,list2)
                        MolCorrDF[ii][jj] = r; MolPvalDF[ii][jj] = p
                        print(str(ii)+'_'+jj)

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
            Optiondict={'Color':LH.MolColor(self.label),#分子名を与えれば色のリストを吐く,
                    'Annotate' : 1,#分子名などをAnnotateするなら
                    'Label':self.label,
                    'Title':'Correlation vs Fasting'
                    }
            plotTmCs(MolCorrDF.T,Optiondict,self.save_dir)#CVのタイムコース受け取って1つのグラフに重ねて描画
            ClstMol,Colordict,ClstDF,ClstColorDF,ColorDF = mD.draw_heatmap(MolCorrDF,1,'MolColor',{'title':'Correlation vs Fasting'},self.save_dir+'Correlation vs Fasting_',cmap='bwr')
