#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 14:58:40 2018

@author: fujita
"""

from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import mpl_toolkits.axes_grid1
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
#import AdjustData011518 as aD
import itertools
#import DataLoader.mkBolusIdx as BI
import matplotlib.cm as cm
import scipy.io
import pandas as pd
import numpy as np
#import Helper.LabelHeler as LH
from scipy.spatial.distance import squareform
#import Helper.ChangeDictToDataFrame as ChDicToDF


def draw_heatmapclustering(a,MolColorSwitch,ColorSwitch,Optindict,save_dir,cmap='bwr'):#bwr#['Reds','Blues']):
    from matplotlib import pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage, dendrogram
    import mpl_toolkits.axes_grid1
    from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
    import itertools
    import matplotlib.cm as cm
    import scipy.io
    import re
    
    metric = 'euclidean'
    method = 'ward'#'centroid'#'ward'#'average'
    threshold_ratio=0.7
    plt.rcParams['font.family'] = 'Arial Unicode MS'
    MolColor = Optindict['MolColor']

    colorlist = ['greenyellow','magenta','deepskyblue','purple','orange','lime','gold','red','blue','green','cyan','black','mediumspringgreen','greenyellow','plum']
    colorlist = ['purple','blue','magenta','deepskyblue','red','orange','gold','cyan','mediumspringgreen','greenyellow','lime']

    vmin = min(a.values.flatten());vmax = max(a.values.flatten())
    v = max(np.abs(vmin), np.abs(vmax)) ## 0真ん中のcm.jetなら　v = max(np.abs(vmin), np.abs(vmax))
    set_link_color_palette(colorlist)
    ih = 70#30 #* len(a.columns) / 13;  40 元は90
    wl = 120#40#ih * 4/3 #* len(a.index)/83 30 元は120
    fig, ax = plt.subplots( figsize=(ih,wl))#元は(3,4)がちょうど良い 論文は30,40
    #
    main_axes = plt.gca()
    #
    divider = make_axes_locatable(main_axes)#(main_axes)(ax)
    
    #xdendro_axes = divider.append_axes("top", 0.5, pad=0)#x軸方向にも分けるなら
    ydendro_axes = divider.append_axes("left",size=5 , pad=0.25)#左側に軸を加える、size=2.5数値で枝の幅決定#size=0.25が良い ih / 6も？長い時は25　size=5,10もある

    plt.sca(ydendro_axes)
    #VI_BolusContinuousは8出切ってる
    threshold,numcluster,distance_sort =  Optindict['Threshold'] ,Optindict['numcluster'],'descending'#'ascending'#'descending'#変動量指標ではこれ使ってた：（基本これにしておく）4.3,13,'descending'#Stdev_Fast:3.9,13:4.5,10:Std_ep:4.5,10  4.7, 11: 4.4,12  :4.2,13  :3.9,14, VI : 3.6,13, BC_VI : 8,13, BC_Delta : 1.5,12, 1.46,13, 1.4,14,1.2,15  BC_FC : 2,13, ### 3.3,13   
    #BC両方_3.6なら34  BC 7.9　#Fujita論文クラスタリングは3.2の13
    plt.rcParams['font.size'] = 3#元は3だった
    #
    ylinkage = linkage(pdist(a, metric=metric), method=method,
                       metric=metric)
    ydendro = dendrogram(ylinkage, orientation='left',leaf_font_size = 0.01,no_labels=1,#truncate_mode = 'level',leaf_font_size = '10',(0.5が良い？, 0.05もある) #no_labels=True, 
                         color_threshold=threshold,above_threshold_color='black',#,link_color_func=lambda k: LeafColors[k]#'descending'#ここでクラスターの切り方変えられる
                         count_sort='False'# distance_sort=distance_sort,#count_sort='False'#
                         #show_leaf_counts=1
                         )#count_sort='False'
    ydendro_axes.invert_yaxis()

    #
    plt.rcParams['lines.linewidth'] = 10#元は1 軸の線幅edge linewidth。囲みの太さ：10なら太め
    plt.rcParams['axes.linewidth'] = 0.4#元は0.4->4
    try:
        if Optindict['Check'] == 'Amino':#AAのラベルなど変える
            MolColor,a = LH.AACheck(MolColor,a,Optindict['AminoCheck'])#'protein','ketogenic','EAA','SemiEAA'
        elif Optindict['Check'] == 'Glc':#糖代謝系のラベルなど変える
            MolColor,a = LH.TCACheck(MolColor,a,'TCA')#'TCA'
        elif Optindict['Check'] == 'AminoGlc':#糖代謝系、AAのラベルなど変える
            MolColor,a = LH.AACheck(MolColor,a,Optindict['AminoCheck'])
            MolColor,a = LH.TCACheck(MolColor,a,'TCA')#'TCA'
    except:
        pass
    plt.axvline(x=threshold, c='k',lw=0.5,ls='dashed')#元は0.5->8 クラスタ分ける棒線たて    
    #plt.gca().set_axis_off()#デンドログラムの枠消す
    plt.gca().tick_params(axis='x', labelsize=100)#* len(a.columns) / 13) labelsize:40,100
    #
    a = a.loc[[a.index[i] for i in ydendro['leaves']]]
    #
    plt.sca(ax)#(main_axes)
    a=a.astype(float)
    if (cmap == 'Reds') or (cmap == 'afmhot_r'):
        aa = ax.imshow(a,  aspect='equal', interpolation='nearest', cmap=cmap, vmin=0, vmax=v)#よく使うのは-4.5, 4.5
    else:
        aa = ax.imshow(a,  aspect='auto', interpolation='nearest', cmap=cmap, vmin=-v, vmax=v)#よく使うのは-4.5, 4.5
    ys, xs = np.meshgrid(range(a.shape[0]),range(a.shape[1]),indexing='ij')
    for (x,y,val) in zip(xs.flatten(), ys.flatten(), a.values.flatten()):
        #plt.text(x,y,'{0:.2f}%'.format(val*100), horizontalalignment='center',verticalalignment='center',size=45, rotation = 90) #小さい時はsize5、大きい時は45
        pass        
    #divider = mpl_toolkits.axes_grid1.make_axes_locatable(main_axes)#これは使わない
    ###
    cluster_assign = fcluster(ylinkage,5,#numcluster,#threshold_ratio*ylinkage[:,2]
                                 criterion='distance')
    # クラスタリング結果の値を取得
    clustered = fcluster(ylinkage, numcluster,#threshold,
                         criterion='maxclust'#'distance'
                         )
    cluster_along_id =clustered[ydendro['leaves']]
    ###
### temp_20190823    

    for ii in range(1,numcluster+1):  #クラスター別れるところに黒線を引く  横線
       cluster_pos = np.where(cluster_along_id==ii)[0][0]#[-1]
       plt.axhline(y=cluster_pos-0.5, lw=15,color='k')  #GHのうえ ：lw=1

### カラーバーの設定
    char = plt.colorbar(aa,pad=0.2,shrink=0.0001)#,orientation ='horizontal')#,shrink=0.6)#元は0.2(pad), 0.25,1(shrink))
    #char.set_clim(-1.,1.2)#色の範囲を手動で
    char.ax.tick_params(labelsize=175)#元は0.5->75 -> 175
    
    plt.gca().yaxis.tick_right()#分子種ラベルを右にやる
    #軸の設定
    #for i in range(a.shape[1]):
    xl = plt.xticks(range(a.shape[1]), a.columns, rotation=270, size=130,color='k')#ize1, 10, 75
    plt.xlabel('Time(min.)',size=10)#論文ではsize=100
    ax.set_xticklabels(a.columns,minor=False)
    #ax.xaxis.label.set_color('r')
    plt.tick_params(axis='both',direction='in',width=0.05)
    #[t.set_color('red') for t in xl.xaxis.label.set_color()]
    bottom, top = ax.get_ylim() # y軸の上限、下限の取得
    plt.yticks(range(a.shape[0]), list(a.index), size=100)#* len(a.index)/20)#元はsize2 小さめは20大きめは35 -> 75
    plt.ylim(top,bottom)
    #文字に色をつける
    b=pd.concat([a,MolColor.reindex(index=a.index)],axis=1)
    b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
#    [t.set_color(i) for (i,t) in
#     zip(list(b['Color']),ax.yaxis.get_ticklabels())] 
### なぜか軸が反対になったら
    #ax.invert_yaxis()   
### 軸ラベル
    if MolColorSwitch == 1:    #分子名に色をつける
        try:
            MolColor = Optindict['MolColor']
        except:
            pass
        if ('Glucose' in list(a.index)) or ('Glucose (mg/dL)' in list(a.index)):#だいたい被験者数の方が少ない   昔：len(a.index) > len(a.columns) or
            if '_' in list(a.index)[0] :#パラメタとかなら
                NewColorList=[]
                for ij in list(a.index):
                    r = re.compile("(.*)(_)(.*)"); d = r.search(ij); NewColorList += [MolColor['MolColor'][d.group(1)] ]              
                [t.set_color(i) for (i,t) in zip(NewColorList,ax.yaxis.get_ticklabels())]
            else:
                b=pd.concat([a,MolColor.reindex(index=a.index)],axis=1)#,join_axes=[a.index])
                b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.yaxis.get_ticklabels())]
        elif len(a.index) < len(a.columns) or 'Glucose' not in list(a.index):#列にパラメタとか分子とか
            if '_' in list(a.columns)[0] :
                NewColorList=[]
                for ij in list(a.columns):
                    r = re.compile("(.*)(_)(.*)"); d = r.search(ij); NewColorList += [MolColor['MolColor'][d.group(1)]]               
                [t.set_color(i) for (i,t) in zip(NewColorList,ax.xaxis.get_ticklabels())]  
            else:
                b=pd.concat([a,MolColor.reindex(index=a.index)],axis=1)
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.xaxis.get_ticklabels())]  
    else:
        ColorSwitch=''
    
    plt.gca().xaxis.set_ticks_position('none')
    plt.gca().yaxis.set_ticks_position('none')
    plt.gca().invert_yaxis()
    tempcolor = ydendro['color_list']
    #plt.title(Optindict['title'])
    
    plt.savefig(save_dir + 'dendrogram'+ColorSwitch+Optindict['title']+'.pdf',bbox_inches="tight")    
    plt.savefig(save_dir + 'dendrogram'+ColorSwitch+Optindict['title']+'.png',bbox_inches="tight")    

    
    #クラスター番号と分子種を辞書に
    ClstMol={}#シトルリンが上
    Colordict={}#色は図の下から枝の組み合わせで入ってる
    ii=0
    for i in range(1,numcluster+1):#クラスター数だけ作る 
        ClstMol.update({'cluster_' + str(i):[]})### numcluster-i
        Colordict.update({'cluster_' + str(i):[]})###
        #clusteriiにまずa.index入れる
        try:#cluster_along_idで連続して番号があれば
            while cluster_along_id[ii]==cluster_along_id[ii+1]:
                ClstMol['cluster_' + str(i)].append(a.index[ii])###
                ii+=1

            ClstMol['cluster_' + str(i)].append(a.index[ii])###

            ii+=1

        except:#cluster_along_idで連続して番号がなければ
            ClstMol['cluster_' + str(i)].append(a.index[ii])###
    #波形の色を指定するためのラベル_一つしかないクラスターは黒って決めちゃう
    rmnum=[]
    for i in range(len(tempcolor)-1):
        if tempcolor[i]=='black':
            rmnum.append(i)
    for i in rmnum:
        tempcolor[i]=tempcolor[i-1]
    #print(tempcolor)
    i=0
    tempcolor2=[]
    for ii in range(numcluster):
        while i<len(tempcolor)-2 and tempcolor[i]==tempcolor[i+1]:#隣が同じ色である限り
                i+=1
        try:
            tempcolor2.append(tempcolor[i])
        except:
            pass
        print(tempcolor2)
        i+=1
    count=0
    #print(tempcolor2)
    for i in range(numcluster):
        if len(ClstMol['cluster_' + str(1+i)])==1:### numcluster-i
            Colordict['cluster_' + str(1+i)].append('black')### numcluster-i
        else:
            try:
                #print(Colordict)
                Colordict['cluster_' + str(1+i)].append(tempcolor2[count])###numcluster-i
                count+=1
            except:
                pass
    
    ClstMolDF = ChDicToDF.ChangeDictToDF(ClstMol)
    #ClstMolDF.columns= list(ClstMolDF.columns).reverse()    
    ClstMolDF.to_excel(save_dir + 'ClstMol.xlsx')
    ClstColorDF = ChDicToDF.ChangeDictToDF(Colordict)
    ClstColorDF.to_excel(save_dir + 'ClstColor.xlsx')
###  いつか調べる20190313 なんのこと？
    
    ColorDF =pd.DataFrame()
    
    ClstNoColorDF = ChDicToDF.ChangeDictToDF_MolColor(ClstMol,Colordict)
    tempDF = MolColor.drop('ClstColor',axis=1)
    NewDF=pd.concat([tempDF,ClstNoColorDF.reindex(index=tempDF.index)],axis=1)#,join_axes=[tempDF.index])
    NewDF.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
    ClstNoColorDF.to_excel(save_dir+'ClstNoColorDF'+Optindict['title']+'.xlsx')
    ColorDF = pd.concat([NewDF,b[['MolColor','Metabo']].reindex(index=NewDF.index)],axis=1)#,join_axes=[NewDF.index])
    #print(Colordict)
      
    return(ClstMol,ClstMolDF,Colordict,b,ClstColorDF,ColorDF,ClstNoColorDF)



def draw_heatmap(a,MolColorSwitch,ColorSwitch,Optindict,save_dir,cmap='bwr'):#['Reds','Blues']):

    metric = 'euclidean'
    method = 'ward'#'average'
    threshold_ratio=0.7
    plt.rcParams['font.family'] = 'Arial Unicode MS'
    data_dir ='/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/' 
    ref_dir = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/To_analize/'
    matBolus = scipy.io.loadmat('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/AllData_BolusOnly_WoUC_UsableTimecourse.mat')
    mat = scipy.io.loadmat('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/AllData.mat')
    BolusData = matBolus['AllData']
    AllData = mat['AllData']
    MolColor = pd.read_excel(ref_dir + '/LabelFujita.xlsx',header=0,index_col=0)    

    EngList = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/MolColor.xlsx',header=0)['English'])
    MolColor.index = EngList
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,index_col=0).iloc[:83,:]
    #MolColor.index = MolColor['MolName']
    max_TimeSign = max(a.max(axis=1))
    min_TimeSign = min(a.min(axis=1))
    #colors = ['deepskyblue','lime','greenyellow','mediumspringgreen','orange','purple','red','blue','magenta','black', 'cyan','k','k','k','k','k']
    colorlist = ['blue','red','purple','orange','mediumspringgreen','greenyellow','lime','deepskyblue','cyan','magenta','gold']
    #030818物理学会用に等代謝とかが強調される感じの色順番にした_Cls13
    colorlist = ['magenta','deepskyblue','purple','orange','lime','green','red','cyan','blue','gold','mediumspringgreen','black','greenyellow']

    vmin = min(a.values.flatten());vmax = max(a.values.flatten())
    v = max(np.abs(vmin), np.abs(vmax))
    set_link_color_palette(colorlist)
    fig, ax = plt.subplots(figsize=(90,60)) #15,10
    #
    #main_axes = plt.gca()
    #
    #divider = make_axes_locatable(main_axes)#(main_axes)(ax)
    
    #xdendro_axes = divider.append_axes("top", 0.5, pad=0)
    #ydendro_axes = divider.append_axes("left",size=0.5, pad=0)#左側に軸を加える、数値で枝の幅決定
    """
    plt.sca(xdendro_axes)
    xlinkage = linkage(pdist(a.T, metric=metric), method=method, metric=metric)
    xdendro = dendrogram(xlinkage, orientation='top', no_labels=True,
                         distance_sort='descending',
                         link_color_func=lambda x: 'blue')
    plt.gca().set_axis_off()
    a = a[[a.columns[i] for i in xdendro['leaves']]]
    """
    
    #plt.sca(ydendro_axes)
    #plt.xlabel('Observation Points', fontsize=55)#
    #
    #plt.xlabel('Distance', fontsize=1)#
    #plt.xscale("log")#
    threshold,numcluster,distance_sort = 5.15,13,'descending'#変動量指標ではこれ使ってた：　4.3,13,'descending'#Stdev_Fast:3.9,13:4.5,10:Std_ep:4.5,10  4.7, 11: 4.4,12  :4.2,13  :3.9,14,
    plt.rcParams['font.size'] = 5#3だった
    #
    #ylinkage = linkage(pdist(a, metric=metric), method=method,
#                       metric=metric)
                       
    #ydendro = dendrogram(ylinkage, orientation='left',leaf_font_size = 10,#truncate_mode = 'level',leaf_font_size = '10', #no_labels=True, 
#                         distance_sort=distance_sort,color_threshold=threshold,above_threshold_color='black'#,link_color_func=lambda k: LeafColors[k]#'descending'#ここでクラスターの切り方変えられる
 #                        )#count_sort='False'
    #
    plt.rcParams['lines.linewidth'] = 5#0.5
    plt.rcParams['axes.linewidth'] = 0.1

    #plt.axvline(x=threshold, c='k',lw=0.5,ls='dashed')#クラスタ分ける棒線
    
    
    #plt.gca().set_axis_off()
    #
    #a = a.ix[[a.index[i] for i in ydendro['leaves']]]
    #
    #plt.sca(ax)#(main_axes)
    try:
        if Optindict['Check'] == 'Amino':#AAのラベルなど変える
            MolColor,ac = LH.AACheck(MolColor,a.T,Optindict['AminoCheck'])#'protein','ketogenic','EAA','SemiEAA'
        elif Optindict['Check'] == 'Glc':#糖代謝系のラベルなど変える
            MolColor,a = LH.TCACheck(MolColor,a.T,'TCA')#'TCA'
        elif Optindict['Check'] == 'AminoGlc':#糖代謝系、AAのラベルなど変える
            MolColor,a = LH.AACheck(MolColor,a.T,Optindict['AminoCheck'])
            MolColor,a = LH.TCACheck(MolColor,a.T,'TCA')#'TCA'
    except:
        pass
    a=a.astype(float)
    aa = ax.imshow(a.T,  aspect='equal', interpolation='nearest', cmap=cmap, vmin=-v, vmax=v)#よく使うのは-4.5, 4.5；aspect='auto';'nearest'
    
    #divider = mpl_toolkits.axes_grid1.make_axes_locatable(main_axes)これは使わなくて良い
    ###
    #cluster_assign = fcluster(ylinkage,5,#numcluster,#threshold_ratio*ylinkage[:,2]
#                                 criterion='distance')
    # クラスタリング結果の値を取得
    #clustered = fcluster(ylinkage, numcluster,#threshold,
     #                    criterion='maxclust'#'distance'
      #                   )
    #cluster_along_id =clustered[ydendro['leaves']]
    ###
    """
    for ii in range(1,numcluster+1):#クラスター別れるところに黒線を引く
        cluster_pos = np.where(cluster_along_id==ii)[0][-1]
        plt.axhline(y=cluster_pos+0.5, lw=1,color='k')
        #print(cluster_pos)
    """    
    #カラーバーの設定
    char = plt.colorbar(aa,pad=0.02)
    char.ax.tick_params(labelsize=120)
    #plt.gca().yaxis.tick_right()#分子種ラベルを右にやる
    #軸の設定
    #for i in range(a.shape[1]):
    plt.xticks(range(a.shape[0]), a.index, rotation=270, size=40,color='k')#15,135
    #ax.set_xticklabels(a.columns,minor=False)
    #ax.xaxis.label.set_color('r')
    plt.tick_params(axis='both',direction='in',width=0.05)
    #[t.set_color('red') for t in xl.xaxis.label.set_color()]
    aindex=list(a.columns)
    #aindex.reverse()
    plt.yticks(range(a.shape[1]),a.columns , size=40)#15,135
    #文字に色をつける
    b=pd.concat([a,MolColor],axis=1)#,join_axes=[a.index])
    b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
#    [t.set_color(i) for (i,t) in
#     zip(list(b['Color']),ax.yaxis.get_ticklabels())]
    
    if MolColorSwitch == 1:    #分子名に色をつける
        if len(a.index) > len(a.columns) or 'Glucose' in list(a.columns):#だいたい被験者数の方が少ない 
            if '_' in list(a.columns)[0] :#パラメタとかなら
                NewColorList=[]
                for ij in list(a.columns):
                    r = re.compile("(.*)(_)(.*)"); d = r.search(ij); NewColorList += [MolColorDict[d.group(1)] ]              
                [t.set_color(i) for (i,t) in zip(NewColorList,ax.yaxis.get_ticklabels())]
            else:
                b=pd.concat([a,MolColor],axis=1)#,join_axes=[a.columns])
                b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.yaxis.get_ticklabels())]
        elif len(a.index) < len(a.columns):#列にパラメタとか分子とか
            if '_' in list(a.columns)[0] :
                NewColorList=[]
                for ij in list(a.columns):
                    r = re.compile("(.*)(_)(.*)"); d = r.search(ij); NewColorList += [MolColorDict[d.group(1)]]               
                [t.set_color(i) for (i,t) in zip(NewColorList,ax.xaxis.get_ticklabels())]  
            else:
                b=pd.concat([a,MolColor],axis=1)#,join_axes=[a.columns])
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.xaxis.get_ticklabels())] 
    elif MolColorSwitch == 'both':#x-y共に色つける
                b=pd.concat([a,MolColor],axis=1)#,join_axes=[a.index])
                b.to_excel(save_dir+'MolColorDF'+Optindict['title']+'.xlsx')
                [t.set_color(i) for (i,t) in zip(list(b[ColorSwitch]),ax.xaxis.get_ticklabels())]
                ColorListrevse = list(b[ColorSwitch]);#ColorListrevse.reverse()
                [t.set_color(i) for (i,t) in zip(ColorListrevse,ax.yaxis.get_ticklabels())]
                xtic = plt.xticks();ytic = plt.yticks()
                turn = sorted(set(list(b['MolColor'])), key=list(b['MolColor']).index)
                pos=0
                #for ii in range(len(turn)):
                 #   pos += list(b['MolColor']).count(turn[ii])
                  #  plt.axvline(x=pos-0.5, c='k',lw=1,ls='solid')#クラスタ分ける棒線
                   # plt.axhline(y=83-(pos+0.5), c='k',lw=1,ls='solid')#クラスタ分ける棒線
                

                
    else:
        ColorSwitch=''
    plt.gca().xaxis.set_ticks_position('none')
    plt.gca().yaxis.set_ticks_position('none')
    #plt.gca().invert_yaxis() ###  なぜかy軸を逆転させる
#    tempcolor = ydendro['color_list']
    #plt.title(Optindict['title'])
    
    plt.savefig(save_dir + 'heatmap'+ColorSwitch+Optindict['title']+'.pdf')
    plt.savefig(save_dir + 'heatmap'+ColorSwitch+Optindict['title']+'.png')
    

    
    #クラスター番号と分子種を辞書に
    ClstMol={}#シトルリンが上
    Colordict={}#色は図の下から枝の組み合わせで入ってる
    ii=0
    for i in range(numcluster):#クラスター数だけ作る
        ClstMol.update({'cluster_' + str(numcluster-i):[]})
        Colordict.update({'cluster_' + str(numcluster-i):[]})
        #clusteriiにまずa.index入れる
        try:#cluster_along_idで連続して番号があれば
            while cluster_along_id[ii]==cluster_along_id[ii+1]:
                ClstMol['cluster_' + str(numcluster-i)].append(a.index[ii])
                ii+=1

            ClstMol['cluster_' + str(numcluster-i)].append(a.index[ii])

            ii+=1

        except:#cluster_along_idで連続して番号がなければ
            ClstMol['cluster_' + str(numcluster-i)].append(a.index[ii])
    #波形の色を指定するためのラベル_一つしかないクラスターは黒って決めちゃう
    rmnum=[]
    """ 
    for i in range(len(tempcolor)-1):
        if tempcolor[i]=='black':
            rmnum.append(i)
    for i in rmnum:
        tempcolor[i]=tempcolor[i-1]
    #print(tempcolor)
    i=0
    tempcolor2=[]
    for ii in range(numcluster):
        while i<len(tempcolor)-2 and tempcolor[i]==tempcolor[i+1]:#隣が同じ色である限り
                i+=1
        try:
            tempcolor2.append(tempcolor[i])
        except:
            pass
        print(tempcolor2)
        i+=1
    count=0
    #print(tempcolor2)
    for i in range(numcluster):
        if len(ClstMol['cluster_' + str(numcluster-i)])==1:
            Colordict['cluster_' + str(numcluster-i)].append('black')
        else:
            try:
                #print(Colordict)
                Colordict['cluster_' + str(numcluster-i)].append(tempcolor2[count])
                count+=1
            except:
                pass
        
    #print(Colordict)

    """ 
    return(ClstMol,Colordict)

