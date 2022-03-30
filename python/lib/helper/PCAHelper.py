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
#import Helper.GraphHelper as GH
#import MolPlot3D as ThreeD
import os
#import Helper.LabelHeler as LH

color_map = LinearSegmentedColormap('color_map',
    {'red': [(0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,1.0,1.0)], #(x,y0,y1) y1→y0y1→y0...
   'green': [(0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,0.0,0.0)],
    'blue': [(0.0,1.0,1.0),(0.5,1.0,1.0),(1.0,0.0,0.0)]})

def DelLabem(num_row1,MolLabel):
 
    r = re.compile("(.*)(_)(.*)") 
    MolLabelNew = []    

    for i in range(num_row1):

        d = r.search(MolLabel[i])
        MolLabelNew.append(d.group(1)) 

    return(MolLabelNew)
    
def calclenScore(pca_score,BiplotSwitch,CompareBolus,MolLabel,save_dir):#B,Cそれぞれの原点からの距離も算出、棒グラフと3Dplot?
    num_col1 = pca_score.shape[1]//2
    num_row1 = pca_score.shape[0]//2#83になる(metabolome)    
    pca_score1 = pca_score[0:num_row1,[0,1]]#B
    pca_score2 = pca_score[num_row1:,[0,1]]   #C    
    PC1norm_array = np.empty(num_row1)
    PC2norm_array = np.empty(num_row1)

    for i in range(num_row1):
        PC1norm_array[i] = np.linalg.norm(pca_score1[i,:])
        PC2norm_array[i] = np.linalg.norm(pca_score2[i,:])
    #try:    
    MolLabelNew=DelLabem(num_row1,MolLabel)
    MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1") # 日本語

    BCNormDF =CompareBolus['VINorm'].T; BCNormDF['Bolus']= list(PC1norm_array); BCNormDF['Continuous']= list(PC2norm_array); BCNormDF['MolColor']= list(MolColorDF['MolColor']); BCNormDF.index = MolLabelNew
    Optiondict={'xlabel':'VINorm','ylabel':'PC12Norm','Annotate':1,'Label':MolLabelNew,'calcR':'pearson','title':'VINorm vs BolusNorm'}
    GH.mkScatterWHist(list(BCNormDF['Norm']),list(PC1norm_array),save_dir,BiplotSwitch['EachMolColor'],Optiondict)#2つのリストの散布図+ヒストグラム
    Optiondict['title'] = 'VINorm vs ContinuousNorm'
    GH.mkScatterWHist(list(BCNormDF['Norm']),list(PC2norm_array),save_dir,BiplotSwitch['EachMolColor'],Optiondict)#2つのリストの散布図+ヒストグラム
    ############# BC間距離とB,Cそれぞれの距離、それぞれの距離間の散布図
    ### Bolus vs BCNorm
    BCNormDF =CompareBolus['VINorm'].T; BCNormDF['Bolus']= list(PC1norm_array); BCNormDF['Continuous']= list(PC2norm_array); BCNormDF['MolColor']= list(MolColorDF['MolColor']); BCNormDF.index = MolLabelNew
    BCNormDF['lengthPC12'] = list(CompareBolus['lengthPC12']);
    Optiondict={'xlabel':'lengthPC12','ylabel':'BolusNorm','Annotate':1,'Label':MolLabelNew,'calcR':'pearson','title':'lengthPC12 vs BolusNorm'}
    Optiondict['y=x'] = 1#y=xの線を足す
    GH.mkScatterWHist(list(BCNormDF['lengthPC12']),list(PC1norm_array),save_dir,BiplotSwitch['EachMolColor'],Optiondict)#2つのリストの散布図+ヒストグラム

    ###  Continuous vs BCNorm
    Optiondict={'xlabel':'lengthPC12','ylabel':'ContinuousNorm','Annotate':1,'Label':MolLabelNew,'calcR':'pearson'}
    Optiondict['title'] = 'lengthPC12 vs ContinuousNorm';Optiondict['y=x'] = 1#y=xの線を足す 
    GH.mkScatterWHist(list(BCNormDF['lengthPC12']),list(PC2norm_array),save_dir,BiplotSwitch['EachMolColor'],Optiondict)#2つのリストの散布図+ヒストグラム

    ### Bolus vs Continuous
    Optiondict={'xlabel':'BolusNorm','ylabel':'ContinuousNorm','Annotate':1,'Label':MolLabelNew,'calcR':'pearson'}
    Optiondict['title'] = 'BolusNorm vs ContinuousNorm';Optiondict['y=x'] = 1#y=xの線を足す
    GH.mkScatterWHist(list(PC1norm_array),list(PC2norm_array),save_dir,BiplotSwitch['EachMolColor'],Optiondict)#2つのリストの散布図+ヒストグラム

    ### Bolus vs Continuous w BCNorm
    Optiondict={'xlabel':'BolusNorm','ylabel':'ContinuousNorm','Annotate':1,'Label':MolLabelNew,'calcR':'pearson'}
    Optiondict['title'] = 'BolusNorm vs ContinuousNorm w BCNorm';Optiondict['y=x'] = 1#y=xの線を足す
    Optiondict['markersize'] = [i*50 for i in CompareBolus['lengthPC12'] ]
    GH.mkScatterWHist(list(PC1norm_array),list(PC2norm_array),save_dir,BiplotSwitch['EachMolColor'],Optiondict)#2つのリストの散布図+ヒストグラム


    ############# 3種の棒グラフ
    PC1SortDF = BCNormDF.sort_values(by='Bolus');List1 = list(PC1SortDF['Bolus']); Title= 'PC1Norm'; xlabel = ''; ylabel = 'Bolus'; Color = PC1SortDF['MolColor']; 
    xticks=list(PC1SortDF.index); size=20; xsize=5; Titlesize=0
    GH.mkSortedBarWHist(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)#任意のリストを代入して棒グラフを描画する
    PC1SortDF = BCNormDF.sort_values(by='Continuous');List1 = list(PC1SortDF['Continuous']); Title= 'PC2Norm'; xlabel = ''; ylabel = 'Continuous'; Color = PC1SortDF['MolColor']; 
    xticks=list(PC1SortDF.index); size=20; xsize=5; Titlesize=0
    GH.mkSortedBarWHist(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir)#任意のリストを代入して棒グラフを描画する
    
    ### 3次元
    #ThreeD.plot3D(BCNormDF,save_dir)

    #except:
     #   pass
    
def calclenbetwScore(pca_score,BiplotSwitch,CompareBolus,MolLabel,save_dir):#PC1,2score間のベクトル長を算出する。
    num_col1 = pca_score.shape[1]//2
    num_row1 = pca_score.shape[0]//2#83になる(metabolome)    
    pca_score1 = pca_score[0:num_row1,[0,1]]#B
    pca_score2 = pca_score[num_row1:,[0,1]]   #C    
    PC12norm_array = np.empty(num_row1)

    for i in range(num_row1):
        PC12norm_array[i] = np.linalg.norm(pca_score1[i,:]-pca_score2[i,:])
        
    try:    
        MolLabelNew=DelLabem(num_row1,MolLabel)
        BCNormDF =CompareBolus['VINorm'] 
        Optiondict={'xlabel':'VINorm','ylabel':'PC12Norm','Annotate':1,'Label':MolLabelNew,'calcR':'pearson','title':'VINorm vs PC12Norm'}
        GH.mkScatterWHist(list(BCNormDF.loc['Norm']),list(PC12norm_array),save_dir,BiplotSwitch['EachMolColor'],Optiondict)#2つのリストの散布図+ヒストグラム
    except:
        pass
    
    return(list(PC12norm_array))

def plotPCAResult(pca_score,facecolor='b',marker='o',cov_ratio=[],label_name=[],fontsize=1, markersize=20):
  #観測点・Loadingをプロット
  plt.axhline(color="k",lw=1)
  plt.axvline(color="k",lw=1)
  
  if isinstance(facecolor,(tuple,list)):
    facecolor = convertColor(facecolor)

  #if markersize==20:#デフォルトならそのまま     
  h = plt.scatter(pca_score[:,0], pca_score[:,1],c=facecolor, marker=marker, s= markersize, edgecolor="none")
  #else:#デフォルトでないなら一つずつ変える
   #   for i in range(len(markersize)):
    #      h = plt.scatter(pca_score[:,0], pca_score[i,1],c=facecolor, marker=marker, s= markersize[i], edgecolor="none")


  ## 説明率配列があれば軸名に追記
  if len(cov_ratio)>0:
    plt.xlabel('PC1 ({0:.0f}%)'.format(cov_ratio[0]*100))
    plt.ylabel('PC2 ({0:.0f}%)'.format(cov_ratio[1]*100))
  ## ラベル配列があれば点のそばに記入
  hori_align = 'left'
  vert_align = 'bottom'
  if type(label_name) == list:#リストなら
      for ii,label in enumerate(label_name):        
        if pca_score[ii,0]<0:
          hori_align = 'right'
        if pca_score[ii,1]<0:
          vert_align = 'top'
        """
        if random.random()<0.5:
          hori_align = 'right'
        if random.random()<0.5:
          vert_align = 'top'
        """
        plt.annotate(label,(pca_score[ii,0],pca_score[ii,1]), size=fontsize,
                     horizontalalignment=hori_align, verticalalignment=vert_align)
  else:
        if pca_score[0,0]<0:
          hori_align = 'right'
        if pca_score[0,1]<0:
          vert_align = 'top'
        plt.annotate(label_name,(pca_score[0,0],pca_score[0,1]), size=fontsize,
                     horizontalalignment=hori_align, verticalalignment=vert_align)      
  return h

def PlotAngle(pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF,BiplotSwitch,CompareBolus):  ## Bolus->Rampの変化角θと仰角Φのscatter
  import math
  num_col1 = pca_score.shape[1]//2
  num_row1 = pca_score.shape[0]//2
  MolLabel=DelLabem(num_row1,MolLabel)

  def vector2angle(v1,v2):
    cal_angle = math.acos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
    if (v1[0]*v2[1]-v1[1]*v2[0] < 0):
      cal_angle = -cal_angle
    return cal_angle
  
  fig = plt.figure()
  fig.set_size_inches(12,9)
  plt.rc("font", family="Arial") #ギリシャ文字など出力のためフォント指定
  
  num_id = num_row1
  theta_array = np.empty(num_id)
  phi_array = np.empty(num_id)
  pca_score1 = pca_score[0:num_row1,[0,1]]
  pca_score2 = pca_score[num_row1:,[0,1]]
  
  for ii in range(num_id):
    theta_array[ii] = vector2angle(pca_score1[ii,:], pca_score2[ii,:])
    phi_array[ii] = vector2angle(-pca_score1[ii,:], pca_score2[ii,:]-pca_score1[ii,:])
  
  #cluster_total = cluster_number1.max()
  #Bolusのクラスタ色
  #cmap = plt.cm.get_cmap('jet', cluster_total)
  #cmap_list = [cmap(i) for i in range(cmap.N)]
  
  for ii in range(num_row1):
    target_pos = ii#(cluster_number1==(ii+1))
    plotPCAResult(np.vstack((theta_array[target_pos],phi_array[target_pos])).T,BiplotSwitch['EachMolColor'][ii],
                      cov_ratio=[], label_name=MolLabel[target_pos],markersize=CompareBolus['lengthPC12'][ii]*10)
  for ii in range(num_row1): # θを2π進める
    target_pos = ii#(cluster_number1==(ii+1))
    plotPCAResult(np.vstack((theta_array[target_pos]+2*math.pi,phi_array[target_pos])).T,BiplotSwitch['EachMolColor'][ii],
                      cov_ratio=[], label_name=MolLabel[target_pos],markersize=CompareBolus['lengthPC12'][ii]*10)
  for ii in range(num_row1): # Φを2π進める
    target_pos = ii#(cluster_number1==(ii+1))
    plotPCAResult(np.vstack((theta_array[target_pos],phi_array[target_pos]+2*math.pi)).T,BiplotSwitch['EachMolColor'][ii],
                      cov_ratio=[], label_name=MolLabel[target_pos],markersize=CompareBolus['lengthPC12'][ii]*10)
  
  plt.axis((-math.pi,2*math.pi,-math.pi,2*math.pi))
  plt.xticks(np.arange(-1,2.1,0.5)*math.pi,
             [r"$-\pi$",r"$-\frac{1}{2}\pi$",0,r"$\frac{1}{2}\pi$",r"$\pi$",
              r"$\frac{3}{2}\pi$",r"$2\pi$"])
  plt.yticks(np.arange(-1,2.1,0.5)*math.pi,
             [r"$-\pi$",r"$-\frac{1}{2}\pi$",0,r"$\frac{1}{2}\pi$",r"$\pi$",
              r"$\frac{3}{2}\pi$",r"$2\pi$"])
  plt.xlabel(r"$\theta$ (Bolus -> O -> Ramp2h)")
  plt.ylabel(r"$\phi$ (O -> Bolus -> Ramp2h)")
  plt.grid(linestyle=':')
  plt.title("Bolus -> Ramp2h angle")

def calcinner_outer_product(pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF):
    #BC間の2点の内積、外積
    def vector2angle(v1,v2):
        cal_angle = math.acos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
        if (v1[0]*v2[1]-v1[1]*v2[0] < 0):
            cal_angle = -cal_angle
        return(cal_angle)
    
    num_col1 = pca_score.shape[1]//2#C
    num_row1 = pca_score.shape[0]//2#B

    innerList=[];outerList=[]

    innerList = [np.dot(pca_score[i,[0,1]], pca_score[i+num_row1,[0,1]]) for i in range(num_row1)]#BのPC1    
    outerList = [np.cross(pca_score[i,[0,1]], pca_score[i+num_row1,[0,1]]) for i in range(num_row1)]#BのPC1
    AngleList = [vector2angle(pca_score[i,[0,1]], pca_score[i+num_row1,[0,1]]) for i in range(num_row1)]#BのPC1



    return(innerList,outerList,AngleList)
    
def VectorDiagram(pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF):
    ## Bolus->Rampをベクトル場で表現今開いているPC平面に追加する形で。
    num_col1 = pca_score.shape[1]//2
    num_row1 = pca_score.shape[0]//2
    
    #cluster_total = cluster_number1.max()
      #Bolusのクラスタ色
    #cmap = plt.cm.get_cmap('jet', cluster_total)
    #cmap_list = [cmap(i) for i in range(cmap.N)]
      
    #for ii in range(cluster_total):
        #target_pos = (cluster_number1==(ii+1))
       # pch.plotPCAResult(pca_score[:num_row1,:][target_pos,:],cmap_list[ii],cov_ratio)
    #Rampのクラスタ色
    #cluster_total = cluster_number2.max()
    #cmap = plt.cm.get_cmap('jet', cluster_total)
    #cmap_list = [cmap(i) for i in range(cmap.N)]
      
    #for ii in range(cluster_total):
       # target_pos = (cluster_number2==(ii+1))
      #  pch.plotPCAResult(pca_score[num_row1:,:][target_pos,:],cmap_list[ii],cov_ratio)
      
     #pch.plotPCAcoef(pca_coef, data_dic['time_str'])
    plt.quiver(pca_score[0:num_row1,[0]],
         pca_score[0:num_row1,[1]],
         pca_score[num_row1:,[0]]-pca_score[0:num_row1,[0]],
         pca_score[num_row1:,[1]]-pca_score[0:num_row1,[1]],
         angles='xy',scale_units='xy',scale=1,width=0.003)
    plt.title("Bolus -> Ramp2h")

     #plt.close()

def PC1_2(scoreDF):#PC1-2で分類
    from sklearn.metrics import confusion_matrix
    PC1_p_num = len(scoreDF[scoreDF[1] >= 0].index)
    PC1_p = list(scoreDF[scoreDF[1] >= 0].index)
    PC1_m_num =len(scoreDF[scoreDF[1] < 0].index)
    PC1_m =list(scoreDF[scoreDF[1] < 0].index)
    PC2_p_num = len(scoreDF[scoreDF[2] >= 0].index)
    PC2_p = list(scoreDF[scoreDF[2] >= 0].index)
    PC2_m_num = len(scoreDF[scoreDF[2] < 0].index)
    PC2_m = list(scoreDF[scoreDF[2] < 0].index)
    
    PC1 = list(scoreDF[1] >= 0)#[0]*len(PC1_m)+[1]*len(PC1_)#0:マイナス、1:プラス
    PC2 = list(scoreDF[2] < 0)#[0]*len(set(PC1_m)  & set(PC2_p))+#0:プラス、1:マイナス
    ForDF = confusion_matrix(PC1, PC2)
    array([[26, 18],
       [13, 26]])
    NewDF = pd.DataFrame(data=None,index=['PC2_+']*6+['PC2_-']*5,columns=['PC1_-']*5+['PC1_+']*6)
    count12px=5;count12py=0;
    for i in list(scoreDF.index):#各分子
        if (i in PC1_p) and (i in PC2_p):#PC1,2+なら
            NewDF.iloc[count12py,count12px]=i
            count12p+=1
            #NewDF[]

def calcVarPC1_2(tempDF):#PC1,2,1+2のばらつきを算出
    #行方向の分散を計算する
    tempVar=[]
    for ii in range(len(tempDF.index)):#行の数だけ
        tempVar.append(np.nanvar(tempDF.loc[ii+1,:]))
    
    return(tempVar[0],tempVar[1],tempVar[0]+tempVar[1])    
def calcScoreVarBC(ScoreDF,label,save_dir):#BC5データで、各分子x条件x各PCとPC1+2でのばらつき描画
    r = re.compile("(.*)(-)(.*)"); 
    MolDF =        pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1") 
    EngLabel = list(MolDF['English']);MolColor=list(MolDF['MolColor'])
    PC1VarList=[];PC2VarList=[];PC12VarList=[];xticks=[];ColorList=[]
    for ii in range(len(EngLabel)):
        #MolLabel = re.compile("(.*)(_)(.*)").search(ii).group(1)
        Mollabel = [s for s in list(ScoreDF.columns) if EngLabel[ii] in s]#まずは分子でとる
        Blist = [ij for ij in Mollabel if r.search(ij).group(3) == 'B']
        Clist = [ij for ij in Mollabel if r.search(ij).group(3) == 'C']
        
        tempDF = ScoreDF[Blist]
        #ある分子の全被験者、PCのスコア
        #tempDF.to_excel(save_dir+EngLabel[ii]+'_Score_dump.xlsx')
        ####各PCにおける分散を計算
        PC1Var,PC2Var,PC12Var = calcVarPC1_2(tempDF)
        PC1VarList.append(PC1Var);PC2VarList.append(PC2Var);PC12VarList.append(PC12Var);xticks.append(EngLabel[ii]+'_B')
        tempDF = ScoreDF[Clist]
        PC1Var,PC2Var,PC12Var = calcVarPC1_2(tempDF)
        PC1VarList.append(PC1Var);PC2VarList.append(PC2Var);PC12VarList.append(PC12Var);xticks.append(EngLabel[ii]+'_C')
        ColorList.append(MolColor[ii]);ColorList.append(MolColor[ii]);
 
    SortDF=pd.DataFrame(data=None,columns=['PC1','PC2','PC12'],index=xticks)
    SortDF['PC1']=PC1VarList;SortDF['PC2']=PC2VarList;SortDF['PC12']=PC12VarList;SortDF['MolColor']=ColorList;SortDF.to_excel(save_dir+'/1stDF.xlsx')
    SortDF=SortDF.dropna(axis=0)
    PC1DF=SortDF.sort_values(by='PC1');    PC2DF=SortDF.sort_values(by='PC2');    PC12DF=SortDF.sort_values(by='PC12')
    GH.mkSortedBarWHist(list(PC1DF['PC1']), 'PC1', '', 'PC1', list(PC1DF['MolColor']), list(PC1DF.index), 20, 5,0, save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
    GH.mkSortedBarWHist(list(PC2DF['PC2']), 'PC2', '', 'PC2', list(PC2DF['MolColor']), list(PC2DF.index), 20, 5,0, save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
    GH.mkSortedBarWHist(list(PC12DF['PC12']), 'PC1+PC2', '', 'PC1+PC2', list(PC12DF['MolColor']), list(PC12DF.index), 20, 5,0, save_dir)#ヒストグラム月の任意のリストを代入して棒グラフを描画する
    Optiondict=dict({'calcR':'spearman','xlabel':'PC1','ylabel':'PC2','title':'PC1vsPC2','Label':list(PC1DF.index),'Annotate':1})
    GH.mkScatterWHist(list(PC1DF['PC1']),list(PC1DF['PC2']),save_dir,list(PC1DF['MolColor']),Optiondict)#2つのリストの散布図+ヒストグラム
def calcVar(tempDF,save_dir):#各PCにおける分散を計算
    #行方向の分散を計算して、一番大きなPCを吐く
    tempVar=[]
    for ii in range(len(tempDF.index)):#行の数だけ
        tempVar.append(np.nanvar(tempDF.loc[ii+1,:]))
    tempDF['Variance'] = tempVar
    PCmax = max(tempVar);PCmaxIdx = tempVar.index(max(tempVar))+1
    return(tempDF,PCmax,PCmaxIdx)

def calcScoreVar(ScoreDF,Englabel,save_dir):#各分子、各PCの被験者間スコア分散算出
    EngLabel = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng.xlsx',header=0,encoding = "ISO-8859-1")['English'])
    EngLabel.remove('3-Hydroxybutyrate')
    MolPCDict=dict({'Mol':[],'PC':[]})
    PCmaxList=[];    PCmaxIdxList=[]
    #ある分子名を含む列を抽出
    for ii in range(len(EngLabel)):
        #MolLabel = re.compile("(.*)(_)(.*)").search(ii).group(1)
        Mollabel = [s for s in list(ScoreDF.columns) if EngLabel[ii] in s]
        tempDF = ScoreDF[Mollabel]
        #ある分子の全被験者、PCのスコア
        #tempDF.to_excel(save_dir+EngLabel[ii]+'_Score_dump.xlsx')
        ####各PCにおける分散を計算
        tempDF, PCmax,PCmaxIdx = calcVar(tempDF,save_dir)
        #tempDF.to_excel(save_dir +EngLabel[ii]+ '_ScoreWVar_dump.xlsx')
        
        PCmaxList.append(PCmax);    PCmaxIdxList.append(PCmaxIdx)
    df = pd.DataFrame(data=[PCmaxList,PCmaxIdxList], columns=EngLabel,index=['maxVar','maxPCIdx'])
    df.to_excel(save_dir + 'maxPCEachMol.xlsx')    
    ######13x83のTableでmax持つPCのみ赤くする
    dataframe = pd.DataFrame(data=None, index=['1','2','3','4','5','6','7','8','9','10','11','12','13'], columns = EngLabel)
    for i in range(len(EngLabel)):
        loc = df[EngLabel[i]]['maxPCIdx']
        dataframe[EngLabel[i]][str(int(loc))] = 1
    dataframe.to_excel(save_dir + 'macPXIdx.xlsx')
    
    #fig, ax = plt.subplots()
        #pca_score = np.log10(pca_score)
    #fig.set_size_inches(pca_score.shape[1], pca_score.shape[0]/5.0)
    
#PC平面に楕円を描く
def drawEllipse(XX,pca_score,ClstMolDF,cov_ratioIdx,MolLabel,ClstColorDF):#楕円用の長軸、短軸を算出する。
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from scipy.stats import chi2
    from scipy.sparse.linalg import svds

    #ClstColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/20180317/ClstColor.xlsx')
    #ClstColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/20190226/Delta_AveByStd_2/ClstColor.xlsx')

    #ClstColorDF=ClstMolDF
    #原点とスコアの重心を中心とした固有値、固有ベクトル
    #fig = plt.figure()
    #ax = plt.axes()
    a = plt.subplot(111, aspect='equal')
    dimention=2
    pca_score2,pca_coef,cov_ratio,pca_components_,W, v = PCA(XX,dimention)
    ScoreThetaMolNamedict = {}
    MolNameList=[]
    EucDistDict = {}
    RadianDict = {}
    #クラスターに含まれるラベル分のスコアを平均すれば良い
    #スコアはラベル順に並んでる
    #for i in range(0,2):#len(pca_score)):#PCの数だけやる
    for j in range(len(ClstMolDF.columns)):#クラスターの数だけやる
        pca_score[j,cov_ratioIdx]#cov_ratioIdx軸のスコア
        MolLabel[j]#上のスコアに対応するラベル
        Con = list(ClstMolDF.columns)
        #Con.reverse()
        TargetMolList = ClstMolDF[Con[j]]
        TargetMolList = TargetMolList.dropna(how='any')
        #print(TargetMolList)
        #ScoreThetaMolNameList=[]
        #TargetMolListのラベルがLabel中何番目か
        PcaScoreAve = np.empty((0,2), int)
        for jj in range(len(TargetMolList)):
            labelloc = np.where(np.array(MolLabel)==TargetMolList[jj])#クラスターj中の分子種たち
            PcaScoreForAve = np.append(np.array([pca_score[labelloc[0][0],cov_ratioIdx]]),np.array([pca_score[labelloc[0][0],cov_ratioIdx+1]]), axis=0)#2次元のPCスコア行列(x,y)を縦に足してく
            PcaScoreAve = np.append(PcaScoreAve,[PcaScoreForAve],axis=0)
            ScoreThetaMolNamedict.update({TargetMolList[jj]:PcaScoreForAve})
        #ScoreThetaMolNamedict.update({ScoreThetaMolNameList})
        x,y = np.mean(PcaScoreAve,axis=0) #これがクラスターjの楕円の中心座標
        
        #クラスタの座標を特異値分解して単位円にかける
        ########################################住友さんは平均0だけ
        PcaScoreAve = colCentering(PcaScoreAve)
        PcaScoreAve = np.array(PcaScoreAve)
        
        #分散共分散行列を特異値分解しても
        #U, s, VT = np.linalg.svd(np.dot(PcaScoreAve.T,PcaScoreAve), full_matrices=True)
        #元の行列を特異値分解しても
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
        
def PCA(X,dimention):#PCAしたい行列Xを入れる（numpyで）
    from sklearn.decomposition import PCA
    #pca = PCA(n_components=dimention)
    #pca.fit(X)
    pca = PCA()
    X[np.isnan(X)] = 0
    pca_score = pca.fit_transform(X)#主成分スコア
    pca.components_, pca.mean_, pca.get_covariance() #主成分Loading(列方向に並んでる)、PCAの平均、共分散
    cv = pca.get_covariance()
    W, v = np.linalg.eig(cv)#固有値、固有ベクトルv[:,0]が第1固有ベクトル

    cov_ratio = pca.explained_variance_ratio_#寄与率
    #print(pca.components_[0,:])    
    pca_coef= pca.components_.T* np.sqrt(pca.explained_variance_)#FactorLoading
    Xd = pca.transform(X)



    return(pca_score,pca_coef,cov_ratio,pca.components_,W, v)
def plot3DPCA(labelProp,pca_score,pca_coef,Label,BiplotSwitch,save_dir, linecolor='r',fontsize=1):
    save_dir_figs = save_dir + '/figs'    
    ScattColorSwitch=2
    if not os.path.isdir(save_dir_figs):
            os.makedirs(save_dir_figs)
    for angle in [0,30,60,90,120,150,180,210,240,270,300,330]:#range(0, 180):
        fig = plt.figure()
        axes = fig.add_subplot(111, projection="3d")
        axes.view_init(30, angle)      
        ############# 散布図

        Color=BiplotSwitch['DotColor']
    
        MolColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',index_col=0)
        #MolColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/20190227/Delta_AveByStd_2/MolColorDFDelta_AveByStd_2.xlsx')
    
        hori_align = 'right'
        vert_align = 'top'
        max_xscore = max(pca_score[:,0])*1.1
        max_yscore = max(pca_score[:,1])*1.2
        min_xscore = min(pca_score[:,0])
        min_yscore = min(pca_score[:,1])
        min_zscore = min(pca_score[:,2])
        max_zscore = max(pca_score[:,2])   
        
        XIdx=0;YIdx=1;ZIdx= 2
        for i in range(len(pca_score)):
            if Color=='Black':#平凡ver
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx],pca_score[i,ZIdx], c='k', label=labelProp[i],s =2)
            elif Color=='EachMol':#被験者x分子PCAにおける分子ごとの色分け
                
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx],pca_score[i,ZIdx], c=BiplotSwitch['EachMolColor'][i], label=labelProp[i],s =2)
            elif Color=='EachMolCluster':#○○x分子PCAにおける分子ごとの色分け            
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx],pca_score[i,ZIdx], c=BiplotSwitch['EachMolColor'][i], label=labelProp[i],s =2)            
            elif Color in ['BC','BWoC','CWoB']:#○○x分子PCAにおけるBCの色分け            
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx],pca_score[i,ZIdx], c=BiplotSwitch['EachMolColor'][i], label=labelProp[i],s =2)            
            elif  Color=='MolColor_AminoGlc':#AAのラベルなど変える：'Amino','Glc','AminoGlc': どちらも
                try:
                    if BiplotSwitch['Check'] == 'Amino':#AAのラベルなど変える
                        MolColor,a = LH.AACheck(MolColorDF,MolColorDF,BiplotSwitch['AminoCheck'])#'protein','ketogenic','EAA','SemiEAA'
                    elif BiplotSwitch['Check'] == 'Glc':#糖代謝系のラベルなど変える
                        MolColor,a = LH.TCACheck(MolColorDF,MolColorDF,'TCA')#'TCA'
                    elif BiplotSwitch['Check'] == 'AminoGlc':#糖代謝系、AAのラベルなど変える
                        MolColor,a = LH.AACheck(MolColorDF,MolColorDF,BiplotSwitch['AminoCheck'])
                        MolColor,a = LH.TCACheck(MolColorDF,MolColorDF,'TCA')#'TCA'
                except:
                    pass
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx],pca_score[i,ZIdx], c=MolColor.loc[labelProp[i]]['MolColor'], label=labelProp[i],s =2)            

            else:
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx],pca_score[i,ZIdx], c=MolColorDF.loc[labelProp[i]][Color], label=labelProp[i],s =3)
            if BiplotSwitch['Label']=='Annotate':
                axes.text(pca_score[i,XIdx],pca_score[i,YIdx],pca_score[i,ZIdx],labelProp[i],size=1,zorder=100)## 各要素にラベルを表示ラベル、x座標、y座標
            axes.set_xlabel('PC1',fontsize=10)
            axes.set_ylabel('PC2',fontsize=10)
            axes.set_zlabel('PC3',fontsize=10)
            axes.set_xlim([min_xscore,max_xscore])
            axes.set_ylim([min_yscore,max_yscore])
            axes.set_zlim([min_zscore,max_zscore])

        scale_ratio = 10#(np.max(np.sqrt(np.sum(np.power(pca_score[:,:3],3),axis=1, keepdims=True)))/
                     #np.max(np.sqrt(np.sum(np.power(pca_coef[:,:3],3),axis=1, keepdims=True))))
        
        for ii in range(pca_coef.shape[0]):#列の数
            #print([0,pca_coef[ii,0]*max_norm],[0,pca_coef[ii,1]*max_norm])
            if ScattColorSwitch==2:
                axes.plot([0,pca_coef[ii,XIdx]*scale_ratio*0.5],[0,pca_coef[ii,YIdx]*scale_ratio*0.5],[0,pca_coef[ii,ZIdx]*scale_ratio*0.5],'r',linewidth=0.5)
            else:
                axes.plot([0,pca_coef[ii,XIdx]*scale_ratio]*0.7,[0,pca_coef[ii,YIdx]*scale_ratio]*0.7,[0,pca_coef[ii,ZIdx]*scale_ratio]*0.7,MolColorDF.loc[Label[ii]][Color],linewidth=0.5)
            if BiplotSwitch['Label']=='Annotate':
                axes.text(pca_coef[ii,XIdx]*scale_ratio*0.5,pca_coef[ii,YIdx]*scale_ratio*0.5,pca_coef[ii,ZIdx]*scale_ratio*0.5,Label[ii],size=5,zorder=100)## 各要素にラベルを表示ラベル、x座標、y座標
     
    
        #pngflag=0
        #if pngflag == 1:
    
            #for angle in range(0, 360):
            #    ax.view_init(30, angle)
        plt.savefig(save_dir_figs +'/{0}_{1:03d}.pdf'.format('3D', angle))
        plt.close()

def mkEachMolColor(pca_score,TimeSignDF):#:分子x被験者PCAにおける、各分子でのいろわえk
    MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    
    if len(pca_score[:,0]) > 0:#雑だが、被験者x分子とかの時
        #if re.compile("(.*)(_)(.*)").search(jj).group(3) in list(PosDF.columns) :
         #   ii=re.compile("(.*)(_)(.*)").search(jj).group(3)            
        alist = []
        for i in list(TimeSignDF.columns):
            #print(i.get_text())
            nowlabel = re.compile("(.*)(_)(.*)").search(i).group(1)
            alist.append(MolColorDF['MolColor'][nowlabel])
        return(alist)

def mkEachMolClsterColor(pca_score,TimeSignDF,ClstColorDF,ClstMolDF):#:分子x○○PCAにおける、各分子のClusterでの色で色分け
    MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolColor_Eng.xlsx')
    
    if len(pca_score[:,0]) > 0:#雑だが、被験者x分子とかの時
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
                print(i)
            alist.append(ClstColorDF[Cluster][0])
        return(alist)
        
def mkEachBCColor(pca_score,TimeSignDF):#:分子x被験者PCAにおける、BCでの色分け
    MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolColor_Eng.xlsx')
    
    if len(pca_score[:,0]) > 0:#雑だが、被験者x分子とかの時
        #if re.compile("(.*)(_)(.*)").search(jj).group(3) in list(PosDF.columns) :
         #   ii=re.compile("(.*)(_)(.*)").search(jj).group(3)            
        alist = []
        for i in list(TimeSignDF.columns):
            #print(i.get_text())
            nowlabel = re.compile("(.*)(_)(.*)").search(i).group(3)
            if nowlabel == 'B':#Bolusは青
                alist.append('blue')
            else:#Continuousは緑
                alist.append('green')
                
        return(alist)   
        
def mkEachBColorWOC(pca_score,TimeSignDF):#:分子x被験者PCAにおける、Bのみ色付き
    MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolColor_Eng.xlsx')
    
    if len(pca_score[:,0]) > 0:#雑だが、被験者x分子とかの時
        #if re.compile("(.*)(_)(.*)").search(jj).group(3) in list(PosDF.columns) :
         #   ii=re.compile("(.*)(_)(.*)").search(jj).group(3)            
        alist = []
        for i in list(TimeSignDF.columns):
            #print(i.get_text())
            nowlabel = re.compile("(.*)(_)(.*)").search(i).group(3)
            nowmollabel = re.compile("(.*)(_)(.*)").search(i).group(1)
            if nowlabel == 'B':#Bolusは青
                alist.append(MolColorDF['MolColor'][nowmollabel])
            else:#Continuousは緑
                alist.append('white')
                
        return(alist)   
        
def mkEachCColorWOB(pca_score,TimeSignDF):#:分子x被験者PCAにおける、Bのみ色付き
    MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolColor_Eng.xlsx',index_col=0)
    
    if len(pca_score[:,0]) > 0:#雑だが、被験者x分子とかの時
        #if re.compile("(.*)(_)(.*)").search(jj).group(3) in list(PosDF.columns) :
         #   ii=re.compile("(.*)(_)(.*)").search(jj).group(3)            
        alist = []
        for i in list(TimeSignDF.columns):
            #print(i.get_text())
            nowlabel = re.compile("(.*)(_)(.*)").search(i).group(3)
            nowmollabel = re.compile("(.*)(_)(.*)").search(i).group(1)

            if nowlabel == 'B':#Bolusは青
                alist.append('white')
            else:#Continuousは緑
                alist.append(MolColorDF['MolColor'][nowmollabel])

                
        return(alist)   


def PCAScatter(pca_score,label,XIdx,YIdx,ClstMol,MolLabel,BiplotSwitch):#スコア、スコアのラベル、分散説明率、散布図の色、クラスターの色？、分子ラベル、スイッチ
    plt.axhline(color="k")
    plt.axvline(color="k")
    Color=BiplotSwitch['DotColor']

    MolColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',index_col=0)
    #MolColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/20190227/Delta_AveByStd_2/MolColorDFDelta_AveByStd_2.xlsx')

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
        if Color=='Black':#平凡ver
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c='k', label=label[i],s =2)
        elif Color=='EachMol':#被験者x分子PCAにおける分子ごとの色分け
            
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=BiplotSwitch['EachMolColor'][i], label=label[i],s =2)
        elif Color=='EachMolCluster':#○○x分子PCAにおける分子ごとの色分け            
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=BiplotSwitch['EachMolColor'][i], label=label[i],s =2)            
        elif Color in ['BC','BWoC','CWoB']:#○○x分子PCAにおけるBCの色分け            
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=BiplotSwitch['EachMolColor'][i], label=label[i],s =2)            
        elif  Color=='MolColor_AminoGlc':#AAのラベルなど変える：'Amino','Glc','AminoGlc': どちらも
                try:
                    if BiplotSwitch['Check'] == 'Amino':#AAのラベルなど変える
                        MolColor,a = LH.AACheck(MolColorDF,MolColorDF,BiplotSwitch['AminoCheck'])#'protein','ketogenic','EAA','SemiEAA'
                    elif BiplotSwitch['Check'] == 'Glc':#糖代謝系のラベルなど変える
                        MolColor,a = LH.TCACheck(MolColorDF,MolColorDF,'TCA')#'TCA'
                    elif BiplotSwitch['Check'] == 'AminoGlc':#糖代謝系、AAのラベルなど変える
                        MolColor,a = LH.AACheck(MolColorDF,MolColorDF,BiplotSwitch['AminoCheck'])
                        MolColor,a = LH.TCACheck(MolColorDF,MolColorDF,'TCA')#'TCA'
                except:
                    pass
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=MolColor.loc[MolLabel[i]]['MolColor'], label=label[i],s =2)            
        elif  Color=='ClstColor_direct':#AAのラベルなど変える：'Amino','Glc','AminoGlc': どちらも
                MolColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200219/RawDelta_NormalizationByStd/ClstNoColorDFmeanRawDelta_NormalizationByStd.xlsx',index_col=0)
                #MolColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200610/RawDelta_NormalizationByStd/ClstNoColorDFmeanRawDelta_NormalizationByStd.xlsx',index_col=0)
                MolColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20201108/RawDelta_NormalizationByStd/ClstNoColorDFCVRawDelta_NormalizationByStd.xlsx',index_col=0)
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=MolColorDF.loc[MolLabel[i]]['ClsterColor'], label=label[i],s =3)
        elif Color=='ClstColor':
                axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=BiplotSwitch['EachMolColor'].loc[MolLabel[i]]['ClsterColor'], label=label[i],s =3)
        
        else:
            axes.scatter(pca_score[i,XIdx],pca_score[i,YIdx], c=MolColorDF.loc[MolLabel[i]][Color], label=label[i],s =3)
        print(BiplotSwitch['Label'] + 'BiplotSwitch')
        if BiplotSwitch['Label']=='Annotate':
            axes.annotate(label[i], xy=(pca_score[i,XIdx],pca_score[i,YIdx]), size=3,horizontalalignment=hori_align, verticalalignment=vert_align)## 各要素にラベルを表示ラベル、x座標、y座標
    axes.set_xlabel('PC'+str(XIdx+1),fontsize=10)
    axes.set_ylabel('PC'+str(YIdx+1),fontsize=10)
    axes.set_xlim([-max_xscore,max_xscore])
    axes.set_ylim([-max_yscore,max_yscore])
    #axes.set_yticks([-6,-4.5,-3,-1.5,0.0,1.5,3.0,4.5,6])
    #axes.set_yticklabels(['-6.0','-4.5','-3.0','-1.5','0.0','1.5','3.0','4.5','6.0'],size=15)

def plotPCAcoef(pca_coef,LabelProp,pca_score,i,j,MolLabel,BiplotSwitch, linecolor='r',fontsize=1):

    Color=BiplotSwitch['DotColor']
    #ScaleFactor= max(pca_score[:,cov_ratioIdx])
    MolColorDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/TimeSeriesAnalysis/20180214/Cls13/MolColorDF.xlsx',index_col=0)
    ScattColorSwitch=2
    if ScattColorSwitch==0:
        Color='Color'
    else:
        Color='ClsterColor'
    scale_ratio = (np.max(np.sqrt(np.sum(np.power(pca_score[:,:2],2),axis=1, keepdims=True)))/
                 np.max(np.sqrt(np.sum(np.power(pca_coef[:,:2],2),axis=1, keepdims=True))))
    for ii in range(pca_coef.shape[0]):#列の数
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
    #plt.xlim([-max_xscore,max_xscore])
    #plt.ylim([-max_yscore,max_yscore])
    

def PCA(X,dimention):#PCAしたい行列Xを入れる（numpyで）
    from sklearn.decomposition import PCA
    #pca = PCA(n_components=dimention)
    #pca.fit(X)
    pca = PCA()
    X[np.isnan(X.astype(float))] = 0
    pca_score = pca.fit_transform(X)
    pca.components_, pca.mean_, pca.get_covariance() #主成分(列方向に並んでる)、PCAの平均、共分散
    cv = pca.get_covariance()
    W, v = np.linalg.eig(cv)#固有値、固有ベクトルv[:,0]が第1固有ベクトル
    #固有ベクトル = PCAの主成分：共分散行列に主成分をかけるのと同じ
    #cv.dot(v[:,0].reshape(2,1)) == v[:,0]*W[0]
    #Projectionする
    #pca_coef = pca.transform(np.identity(X.shape[1])*max_norm) 
    cov_ratio = pca.explained_variance_ratio_
    #print(pca.components_[0,:])    
    max_norm = np.max(np.sqrt(np.sum(np.power(X.astype(float),2), axis=1))) #最大ノルムでスケーリング
    #pca_coef = pca.transform(np.identity(X.shape[1])*np.sqrt(pca.explained_variance_)) #単位行列を変換→主成分空間上で元の軸を表す
    pca_coef= pca.components_.T* np.sqrt(pca.explained_variance_)#これだとFactorLoading
    #pca_coef = scipy.stats.pearsonr(pca_score,X).T
    #pca_coef = np.transpose(pca.components_)
    Xd = pca.transform(X)
    #print(Xd)
    # covariance matrix x eigenvector
    #固有値の根×固有ベクトルが院試負荷量
    #print(v.size,W.size)
    #print(v)
    # display 

    # 主成分の寄与率を出力する
    #print(pca.components_)
    #print('各次元の寄与率: {0}'.format(pca.explained_variance_ratio_))
    #print('累積寄与率: {0}'.format(sum(pca.explained_variance_ratio_)))
    return(pca_score,pca_coef,cov_ratio,pca.components_.T,W, v)

def colCentering(data_2d):
  row_mean = np.mean(data_2d, axis=0,keepdims=True)
  return data_2d - row_mean

def plotPCACovRatio(cov_ratio):
  # 分散の説明率プロット
  num_var = len(cov_ratio)
  x_tick = np.arange(1,num_var+1)
  plt.bar(x_tick, cov_ratio)
  plt.plot(x_tick, np.cumsum(cov_ratio),'-o', mfc='none',mec='b',mew=8,linewidth=3)
  plt.xticks(x_tick, fontsize=40)#20
  plt.yticks(fontsize=40)#20

  plt.axis([1-0.4, num_var+0.5, 0,1])#左に寄せる
  plt.xlabel("Number of PC", fontsize=40)#10
  plt.ylabel("Explained variance ratio", fontsize=40)#10
  plt.rcParams['axes.linewidth'] = 1.5# 軸の線幅edge linewidth。囲みの太さ
  
def heatPCA(data_2d,label_name): 
  import mpl_toolkits.axes_grid1
  num_x = data_2d.shape[1]
  num_y = data_2d.shape[0]
  print('y')
  print(data_2d)
  c_max = max([data_2d.max(), abs(data_2d.min())])
  colors = np.array([plt.cm.hsv(i/len(data_2d)) for i in range(len(data_2d))])
  ax = plt.imshow(data_2d,aspect='auto',interpolation='nearest',
             cmap='PuOr_r',#color_map, 
### 20200214 ヒートマップの色はここで制御
             vmin=-c_max,vmax=c_max)
  plt.gca().xaxis.tick_top()
  #divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax) ###temp_20191006 なくても問題ない？
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
  plt.rcParams['axes.linewidth'] = 1.5# 軸の線幅edge linewidth。囲みの太さ

  #if ChangeColorSwitch == 1:#ラベルの色変え、ラベルが代謝順になっていないと使えないが
  #MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolColor_Eng.xlsx')
  #[t.set_color(i) for (i,t) in
   #  zip(list(MolColorDF['MolColor']),ax.yaxis.get_ticklabels())]

  #plt.ylabel(')
  plt.gca().get_xaxis().set_ticks_position('none')
  plt.gca().get_yaxis().set_ticks_position('none')