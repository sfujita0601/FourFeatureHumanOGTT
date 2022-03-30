#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 11:32:53 2018
#ぶち込まれたデータで大人しくいろんなグラフを作る
@author: fujita
"""

from numpy.random import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import math
from scipy.stats import pearsonr, spearmanr
import pandas as pd
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D

'''
ベン図を書くテスト
参考資料：https://pypi.python.org/pypi/matplotlib-venn
'''

def Piledbargraph(cov_ratio,numlevel):#辞書の形でcov_ratioいれる
  # 分散の説明率プロット
  num_var = len(cov_ratio['1'])
  x_tick = np.arange(1,num_var+1)
  cov_ratiocum=[]
  count=0
  for i in range(numlevel):
      y2=cov_ratio[str(i+1)]
      if count>0:
          plt.bar(x_tick, y2, bottom=y1, color=cm.GnBu_r(count/numlevel))#, width=0.5)
          y1=[ y1[i] +y2[i] for i in range(len(y1))]
      else:
          plt.bar(x_tick, y2, color=cm.GnBu_r(count/numlevel))# width=0.5)
          y1=y2
          # 棒グラフ内に数値を書く
      for x, y in zip(x_tick, y1):
        plt.text(x, y, count+1, ha='center', va='center',size=1)
      cov_ratiocum+=y2
      count+=1
            
  
  #plt.plot(x_tick, np.cumsum(cov_ratiocum),'-o', mfc='none',mec='b',mew=8,linewidth=3)
  plt.xticks(x_tick, fontsize=40)#20
  plt.yticks(fontsize=40)#20

  plt.axis([1-0.4, num_var+0.5, 0,1])#左に寄せる
  plt.xlabel("Number of PC", fontsize=40)#10
  plt.ylabel("Explained variance ratio", fontsize=40)#10
  plt.rcParams['axes.linewidth'] = 1.5# 軸の線幅edge linewidth。囲みの太さ    
    

def move_sn_y(offs=0,host=0, dig=0, side='left', omit_last=False):
    """Move scientific notation exponent from top to the side.
    
    Additionally, one can set the number of digits after the comma
    for the y-ticks, hence if it should state 1, 1.0, 1.00 and so forth.

    Parameters
    ----------
    offs : float, optional; <0>
        Horizontal movement additional to default.
    dig : int, optional; <0>
        Number of decimals after the comma.
    side : string, optional; {<'left'>, 'right'}
        To choose the side of the y-axis notation.
    omit_last : bool, optional; <False>
        If True, the top y-axis-label is omitted.

    Returns
    -------
    locs : list
        List of y-tick locations.

    Note
    ----
    This is kind of a non-satisfying hack, which should be handled more
    properly. But it works. Functions to look at for a better implementation:
    ax.ticklabel_format
    ax.yaxis.major.formatter.set_offset_string
    """

    # Get the ticks
    locs = host.get_yticks()

    # Put the last entry into a string, ensuring it is in scientific notation
    # E.g: 123456789 => '1.235e+08'
    llocs = '%.3e' % locs[-1]

    # Get the magnitude, hence the number after the 'e'
    # E.g: '1.235e+08' => 8
    yoff = int(str(llocs).split('e')[1])

    # If omit_last, remove last entry
    if omit_last:
        slocs = locs[:-1]
    else:
        slocs = locs

    # Set ticks to the requested precision
    form = r'$%.'+str(dig)+'f$'
    plt.yticks(locs, list(map(lambda x: form % x, slocs/(10**yoff))))

    # Define offset depending on the side
    if side == 'left':
        offs = -.18 - offs # Default left: -0.18
    elif side == 'right':
        offs = 1 + offs    # Default right: 1.0
        
    # Plot the exponent
    plt.text(offs, .98, r'$\times10^{%i}$' % yoff, transform =
            plt.gca().transAxes, verticalalignment='top')

    # Return the locs
    return locs
    
def plot3D(DF,OptionDict,save_dir):#3列分あるDFを入れて、plot   
    #MolColorDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary.xlsx',header=0,encoding = "ISO-8859-1")
    #MolColorDF = MolColorDF.drop('高感度CRP'); MolColorDF = MolColorDF.drop('3ハイドロキシ酪酸')   
    try:
        if OptionDict['plot3DColor'] == 'MolColor':
            ColorList = list(DF['MolColor']);
        if OptionDict['plot3DColor'] == 'ClstColor':
            ColorList = list(MolColorDF['ClstColor'])
        if OptionDict['plot3DColor'] == 'PC12':
            ColorList = OptionDict['LoadingColor']
    except:
        pass
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')#ax = Axes3D(fig,projection='3d')    
    
    ColLabel = list(DF.columns)
    IdxLabel = list(DF.index)
    #それぞれの軸で最大1,
    #try: #描画する
    X = list(DF[ColLabel[1]])
    Y = list(DF[ColLabel[0]])
    Z = list(DF[ColLabel[2]])
    if OptionDict['plot3Dseparation'] == 'Discrete':#離散的に
        stcount=0;edcount=13;inum=len(X)//13;count=0
        for i in range(inum):
            if count == 83:
                count=0
            #ax.scatter(X, Y, Z,'o',c=ColorList)
            #ax.scatter(X, Y, Z,'o',c=Y,cmap='hsv')
            #ax.scatter(X[stcount:edcount], Y[stcount:edcount], Z[stcount:edcount],'-',lw=1,c=np.array(Z[stcount:edcount]),cmap='jet')
            ax.plot(X[stcount:edcount], Y[stcount:edcount], Z[stcount:edcount],lw=1,c=cm.jet(count/83))#c=Y,cmap='hsv')
            #for ij# in range(stcount,edcount):
                #ax.plot(X[ij],Y[ij],Z[ij],c=cm.hsv(count/83))
            count+=1
            stcount=edcount;edcount+=13
        #ax.plot_surface(X, Y, Z)
    #except:
     #   pass
    elif OptionDict['plot3Dseparation'] == 'Continuous':#連続的に
        stcount=0;edcount=13;inum=len(X)//13;count=0#時系列1本ずつ
        Loading = OptionDict['Loading'] 
        Lmax = np.max(Loading); Lmin=np.min(Loading)
        for i in range(inum):
            if count == 8300000:#Loadingの値で色を決めたい。
                count=0
            #ax.scatter(X, Y, Z,'o',c=ColorList)
            #ax.scatter(X, Y, Z,'o',c=Y,cmap='hsv')
            #ax.scatter(X[stcount:edcount], Y[stcount:edcount], Z[stcount:edcount],'-',lw=1,c=np.array(Z[stcount:edcount]),cmap='jet')
            ax.plot(X[stcount:edcount], Y[stcount:edcount], Z[stcount:edcount],lw=1,c=cm.jet(((Loading[count]-Lmin)/(Lmax-Lmin))))#c=Y,cmap='hsv')
            #for ij# in range(stcount,edcount):
                #ax.plot(X[ij],Y[ij],Z[ij],c=cm.hsv(count/83))
            count+=1
            stcount=edcount;edcount+=13
    else:
        #ax.plot(X, Y, Z,lw=1)#,c=cm.jet(((Loading[count]-Lmin)/(Lmax-Lmin))))#c=Y,cmap='hsv')
        
        #ax.plot_surface(np.array(X), np.array(Y), np.array(Z), rstride=1, cstride=1, cmap=cm.coolwarm)
        stcount=0;edcount=1;inum=len(X)//13;count=0#時系列1本ずつ
        Loading = OptionDict['Loading'] 
        Lmax = np.max(Loading); Lmin=np.min(Loading)
        for i in range(len(X)):#全店鬱ver.
            if count == 8300000:#Loadingの値で色を決めたい。
                count=0
            #ax.scatter(X, Y, Z,'o',c=ColorList)
            #ax.scatter(X, Y, Z,'o',c=Y,cmap='hsv')
            #ax.scatter(X[stcount:edcount], Y[stcount:edcount], Z[stcount:edcount],'-',lw=1,c=np.array(Z[stcount:edcount]),cmap='jet')
            ax.scatter(X[stcount:edcount], Y[stcount:edcount], Z[stcount:edcount],lw=1,c=cm.jet(((Loading[count]-Lmin)/(Lmax-Lmin))))#c=Y,cmap='hsv')
            #for ij# in range(stcount,edcount):
                #ax.plot(X[ij],Y[ij],Z[ij],c=cm.hsv(count/83))
            count+=1
            stcount=edcount;edcount+=1
            
        
    try:    #Annotaionする
        if OptionDict['plot3DAnnotation'] == 0:
            for ij in range(len(IdxLabel)):#分子名のプロット
                #ax.scatter(PropX[ij], PropY[ij], PropZ[ij])
                ax.text(X[ij], Y[ij], Z[ij],IdxLabel[ij],size=5,zorder=100,color = ColorList[ij] )#, transfoerm=ax.transAxes)     
                #ax.plot(X, Y, Z,'o')#),zorder=100,color = ColorList,ms=4, mew=0.5)#, transfoerm=ax.transAxes)     
    except:
        pass
    #for jj in range(numMol):#分子種×時間のプロット
        #ax.scatter(MolTimeX[jj], MolTimeY[jj], MolTimeZ[jj])
     #   ax.text(MolTimeX[jj], MolTimeY[jj], MolTimeZ[jj], MolTimeLabel[jj],size=5,zorder=100)#, transfoerm=ax.transAxes)    

    #for kl in range(len(UndrThsh[0])): #スキャッター
     #   ax.plot([PropX[UndrThsh[1][kl]], MolTimeX[kl]], [PropY[UndrThsh[1][kl]], MolTimeY[kl]], [PropZ[UndrThsh[1][kl]],MolTimeZ[kl]], linewidth=0.5)    
    #plt.axis("off")
    #ax.spines["right"].set_color("none")  # 右消し
    #ax.spines["left"].set_color("none")   # 左消し
    #ax.spines["top"].set_color("none")    # 上消し
    #ax.spines["bottom"].set_color("none") # 下消し  
    # 軸ラベルの設定
    ax.set_xlabel(ColLabel[1],size=10)
    ax.set_ylabel(ColLabel[0],size=10)
    ax.set_zlabel(ColLabel[2],size=10)
    
    ax.set_xlim(np.min(X), np.max(X))
    ax.set_ylim(np.min(Y), np.max(Y))
    ax.set_zlim(np.min(Z), np.max(Z))
    
    ax.tick_params(labelsize=10)#direction = "inout", length = 5, colors = "blue")

    #save_dir_figs = save_dir + '/figs'
    
    #if not os.path.isdir(save_dir_figs):
     #   os.makedirs(save_dir_figs)
    #for angle in range(0, 360):
     #   ax.view_init(30, angle)
      #  plt.savefig(save_dir_figs +'/{0}_{1:03d}.jpg'.format('3D', angle))        
    plt.savefig(save_dir+'3D.pdf')
    for angle in range(0, 360):
        ax.view_init(30, angle)
        plt.savefig(save_dir + "figs/{0}_{1:03d}.jpg".format('aaall', angle))

    #plt.show()    

def mkvenn():
    from matplotlib_venn import venn2
    from matplotlib import pyplot
    '''
    ベン図を描く
    '''
    print("run venn diagram")
 
    venn2(subsets=(3, 2, 1))
    pyplot.show()
    
    
def colCentering(data_2d):
  row_mean = np.mean(data_2d, axis=0,keepdims=True)
  return data_2d - row_mean
#PC平面に楕円を描く、2次元の行列と、クラスターのラベルを入れる
def drawEllipse(ax1,XX,groups, MolLabel,ColLabel1,ColLabel2):#楕円用の長軸、短軸を算出する。
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from scipy.stats import chi2
    from scipy.sparse.linalg import svds

    #ClstColorDF=ClstMolDF
    #原点とスコアの重心を中心とした固有値、固有ベクトル

    a = ax1#plt.subplot(111, aspect='equal')
    dimention=2
    #pca_score,pca_coef,cov_ratio,pca_components_,W, v = PCA(XX,dimention)
    ScoreThetaMolNamedict = {}
    MolNameList=[]
    EucDistDict = {}
    RadianDict = {}
    #クラスターに含まれるラベル分のスコアを平均すれば良い
    #スコアはラベル順に並んでる
    #for i in range(0,2):#len(pca_score)):#PCの数だけやる
    for status, group in groups:
    #for j in range(len(ClstMolDF.columns)):#クラスターの数だけやる
        TargetMolList = list(group.index)
        Color = group['color'][0]
        #TargetMolList = TargetMolList.dropna(how='any')
        #print(TargetMolList)
        #ScoreThetaMolNameList=[]
        #TargetMolListのラベルがLabel中何番目か
        PcaScoreAve = np.empty((0,2), int)#各クラスターのPCスコアの平均を入れる
        for jj in range(len(TargetMolList)):#あるクラスターに含まれている分子のスコアを足す
            #labelloc = np.where(np.array(MolLabel)==TargetMolList[jj])#クラスターj中の分子種たち
            PcaScoreForAve = np.append(np.array([XX.loc[TargetMolList[jj],ColLabel1]]),np.array([XX.loc[TargetMolList[jj],ColLabel2]]), axis=0)#2次元の行列(x,y)を縦に足してく
            PcaScoreAve = np.append(PcaScoreAve,[PcaScoreForAve],axis=0)#あるクラスターの(x,y)を足す
            ScoreThetaMolNamedict.update({TargetMolList[jj]:PcaScoreForAve})
        #ScoreThetaMolNamedict.update({ScoreThetaMolNameList})
        x,y = np.mean(PcaScoreAve,axis=0) #これがクラスターjの楕円の中心座標
        
        #クラスタの座標を特異値分解して単位円にかける
        ########################################住友さんは平均0だけ
        PcaScoreAve = colCentering(PcaScoreAve)
        PcaScoreAve = np.array(PcaScoreAve)
        

        U, S, V = np.linalg.svd(PcaScoreAve, full_matrices=True)
       
        #内積か分散共分散入れるかスイッチ
        data = np.cov(PcaScoreAve)
        data = np.dot(PcaScoreAve.T,PcaScoreAve)
        W, V_pca = np.linalg.eigh(data)
        index = W.argsort()[::-1]#固有値大きい順になっていないようなのでソート
        V_pca = V_pca[:, index]#ソート順に入れ直す

        #result1 = math.atan2(V[0][1],V[0][0])

        degree = math.degrees(math.atan2(V.T[1][0],V.T[0][0]))#result1,rad

        Theta = degree
        a.plot(x,y, c=Color,markersize =10,marker='+')
        a.set_ylim([-0.1, 1.1])
        #a.set_xlim([0, 0.7])
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
        
        if len(S)>1:
            #print([s[0]-np.abs(x),s[1]-np.abs(y)])
            e = patches.Ellipse(xy=(x,y), width=S[0], height=S[1], angle=Theta,fc='none',ec=Color,linewidth=1)
            a.add_artist(e)
            #temp = np.dot(I,np.diag(s)/2)
            #tempIV = np.dot(temp,V)
            #Ell = tempIV + [x,y]
            #plt.scatter(Ell[:,0],Ell[:,1],s=0.1)
        else:
            pass

def mkScatterYAngleWHist(list1,list2,save_dir,ColorList,Optiondict):#2つのリスト(2つ目は角度）の散布図+ヒストグラム
    
    fig = plt.figure(figsize=(8,4))
    
    try:
        if Optiondict['calcR']=='pearson':
            r, p = pearsonr(list1,list2)#ピアソン
        elif Optiondict['calcR']=='spearman':
            r, p = spearmanr(list1,list2)#スピアマン
    except:
            Optiondict['calcR']='pearson';r, p = pearsonr(list1,list2)#ピアソン
        
    RankTimeVarDF = pd.DataFrame(data=None,columns=['Rank','TimeVar'])
    RankTimeVarDF['Rank'] = list1; RankTimeVarDF['TimeVar']= list2
    # サブプロットを8:2で分割
    """
    ax1 = fig.add_axes((0, 0.75, 0.75, 0.75))
    ax2 = fig.add_axes((0, 0.75, 0.75, 0.25), sharex=ax1)
    ax3 = fig.add_axes((0.75, 0, 0.25, 0.75), sharey=ax1)
    """#fig.add_axes((左下, 左上, x長さ, y長さ)x座標、y座標、幅、高さ
    ax1 = fig.add_axes((0, 0.2, 0.8, 0.8))
    ax2 = fig.add_axes((0, 1, 0.8, 0.3), sharex=ax1)
    ax3 = fig.add_axes((0.8, 0.2, 0.2, 0.8), sharey=ax1)
    
    try:
        markersize= Optiondict['markersize'] 
    except:
        markersize= 50
        
    # 散布図のx軸のラベルとヒストグラムのy軸のラベルを非表示
    #ax1.tick_params(labelbottom="off");#ax2.tick_params(color='white')
    ax1.tick_params(labelsize=15,direction='out')
    ax2.tick_params(labelleft="true",labelbottom="off",direction='out',labelsize=15);ax2.tick_params(axis='x',color='white')
    ax3.tick_params(labelleft="false",labelsize=15);ax3.tick_params(axis='y',color='white',direction='out')
    ax1.scatter(list1,list2,c=ColorList,s=markersize)# = sns.jointplot(x=DFCol[i]+ str(num[jj]) , y='CV of log10(Param'+DF2Col[ii]+')',color=ColorList,space=0, data=SNSDF)
    ax2.hist(list1, bins=20,ec='black')
    ax3.hist(list2, bins=20,orientation="horizontal",ec='black')
    plot_axis = plt.axis()

    ax1.set_xlabel(Optiondict['xlabel'],fontsize=20)
    ax1.set_ylabel(Optiondict['ylabel'],fontsize=20)  
  
    ax1.set_ylim([0,0.8])
    #ax1.set_xticks([0, 0.2, 0.4, 0.6])
    #ax1.set_xticklabels(['0.0', '0.2', '0.4', '0.6'])#, rotation=30, fontsize='small')
    ax1.set_yticks([0, 0.2, 0.4, 0.6])
    ax1.set_yticklabels(['0.0', '0.2', '0.4', '0.6'])#, rotation=30, fontsize='small')
###
    ###### Angle用に
    ax1.set_yticks(np.arange(-1,1.1,0.5)*math.pi)
    ax1.set_yticklabels([r"$-\pi$",r"$-\frac{1}{2}\pi$",0,r"$\frac{1}{2}\pi$",r"$\pi$",
              r"$\frac{3}{2}\pi$",r"$2\pi$"])  
    #y=xの直線描画
    xmin, xmax, ymin, ymax = ax1.axis() 
    x = np.arange(xmin, 1, 0.1)
    #ax1.plot(x,x,c='red',linewidth=1)
    
    #y軸のmax、minを求める
    ymin=min(list2); ymax=max(list2)
    xmin=min(list1); xmax=max(list1)
    
    try:
        if Optiondict['calcR'] in ['pearson','spearman']:
            #ax1.set_ylim([-0.1, 1.1])
            Roundp = round(p,100)
            b = '%.2e'%Roundp
            TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
            p='{:e}'.format(p)
                                
            TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
            #昔：ax1.annotate('$\it{R}$ = ' + str(round(r,3)) + ', $\it{p}$ = ' + TenRoundp,fontsize=12, xy=(0.43,0.9))
            #ax1.annotate('$\it{R}$ = ' + str(round(r,3)) + ', $\it{p}$ =' +str(Roundp),fontsize=12, xy=(0.5,0.7))
            ax1.annotate('$\it{R}$ = ' + str(round(r,3)) + ', $\it{p}$ = ' + TenRoundp,fontsize=12, xy=(xmin+0.1,ymax-0.1))    
    except:
        pass
    try:
        if Optiondict['mkScatterWHist_drawEllipse']==1:#散布図とヒストグラムのために楕円フィットするなら
            NewDF = Optiondict['NewDF'];groups = Optiondict['groups']; MolLabel = Optiondict['MolLabel']
            drawEllipse(ax1,NewDF,groups, MolLabel,'TimeVar','AveCorr')#2つのリストの散布図+ヒストグラム
    except:
        pass
    if Optiondict['Annotate'] == 1:#分子名などをAnnotateするなら
        for i in range(len(list1)):
            ax1.annotate(Optiondict['Label'][i],fontsize=2, xy=(list1[i],list2[i]))
    try:        
        if Optiondict['y=x'] == 1:#y=xの線を足す
                x = np.linspace(xmin,xmax)  # xの値域(0, 1, 2, 3)
                y = x             # 直線の式
                
                ax1.plot(x,y,"r-",linewidth=1)#,label='y=x') 
    except:
        pass
    
    title = Optiondict['title']+Optiondict['calcR']
    plt.savefig(save_dir +'Scatter_'+title+'.pdf',format='pdf',bbox_inches="tight")
    RankTimeVarDF.to_excel(save_dir + 'DF_'+title+'.xlsx')
      
def mkScatterWHist(list1,list2,save_dir,ColorList,Optiondict):#2つのリストの散布図+ヒストグラム
    xs, ys = (6,6)#Fig3 なら12,6 #何用？16.6,9 1対1は 8:8
    fig = plt.figure(figsize=(xs,ys))
    
    try:
        if Optiondict['calcR']=='pearson':
            r, p = pearsonr(list1,list2)#ピアソン
        elif Optiondict['calcR']=='spearman':
            r, p = spearmanr(list1,list2)#スピアマン
    except:
            Optiondict['calcR']='pearson';r, p = pearsonr(list1,list2)#ピアソン
        
    RankTimeVarDF = pd.DataFrame(data=None,columns=[Optiondict['xlabel'],Optiondict['ylabel']])
    RankTimeVarDF[Optiondict['xlabel']] = list1; RankTimeVarDF[Optiondict['ylabel']]= list2
    # サブプロットを8:2で分割
    if (xs==12) and (ys==6):
        ax1 = fig.add_axes((0, 0.25, 0.75, 0.75))
        ax2 = fig.add_axes((0, 1, 0.75, 0.25), sharex=ax1)#上？
        ax3 = fig.add_axes((0.75, 0.25, 0.17, 0.75), sharey=ax1)#右？
    else:
        ax1 = fig.add_axes((0, 0.25, 0.75, 0.75))
        ax2 = fig.add_axes((0, 1, 0.75, 0.25), sharex=ax1)#上？
        ax3 = fig.add_axes((0.75, 0.25, 0.25, 0.75), sharey=ax1)#右？        
    #fig.add_axes((左下, 左上, x長さ, y長さ)x座標、y座標、幅、高さ
    #ax1 = fig.add_axes((0, 0.2, 0.8, 0.8))
    #ax2 = fig.add_axes((0, 1, 0.8, 0.3), sharex=ax1)
    #ax3 = fig.add_axes((0.8, 0.2, 0.2, 0.8), sharey=ax1)
    
    try:
        markersize= Optiondict['markersize'] 
    except:
        markersize= 100 #元は50
        
    # 散布図のx軸のラベルとヒストグラムのy軸のラベルを非表示
    #ax1.tick_params(labelbottom="off");#ax2.tick_params(color='white')
    ax1.tick_params(labelsize=20,direction='out')
    ax3.tick_params(labelleft=False,labelsize=20,axis='y',color='white',left=False,direction='out')#;ax3.tick_params()
    ax3.tick_params(labelleft=False,labelsize=20,axis='x')#;ax3.tick_params()
    #エラーバーを書くなら、
    try:
        if Optiondict['errorbar']==1:
            x_errlist=Optiondict['x_err']
            y_errlist=Optiondict['y_err']
            for x, y, x_err,y_err,c in zip(list1,list2, x_errlist, y_errlist,ColorList):
                ax1.errorbar(x, y, xerr = x_err, yerr = y_err,  fmt=c, elinewidth=0.5,capsize=1, ecolor='black')
            ax1.scatter(list1,list2,c=ColorList,s=markersize, facecolor='white')# = sns.jointplot(x=DFCol[i]+ str(num[jj]) , y='CV of log10(Param'+DF2Col[ii]+')',color=ColorList,space=0, data=SNSDF)

    except:
        ax1.scatter(list1,list2,c=ColorList,s=markersize,edgecolors='black',linewidths=0.1)# = sns.jointplot(x=DFCol[i]+ str(num[jj]) , y='CV of log10(Param'+DF2Col[ii]+')',color=ColorList,space=0, data=SNSDF) , alpha=0.4)
        #ax1.scatter(list1,list2,edgecolors=ColorList,s=markersize,facecolor='None',linewidth=0.2)# = sns.jointplot(x=DFCol[i]+ str(num[jj]) , y='CV of log10(Param'+DF2Col[ii]+')',color=ColorList,space=0, data=SNSDF)

        print('RRRR')        
    ax2.hist(list1, bins=20,ec='black')
    ax3.hist(list2, bins=20,orientation="horizontal",ec='black')
    ax2.tick_params(labelleft="true",labelbottom=False,direction='out',labelsize=5,axis='x',color='white');
    ax2.tick_params(labelleft="true",labelsize=20,axis='y');
    plot_axis = plt.axis()
    #ax1.set_xlim([ 0.0005,0.0025])
### 枠線制御

    axis=['top','bottom','left','right']
    line_width=[1,1,1,1]
    
    for a,w in zip(axis, line_width):  # change axis width
        ax1.spines[a].set_linewidth(w)
        ax2.spines[a].set_linewidth(w)        
        ax3.spines[a].set_linewidth(w)
        
    ax1.set_xlabel(Optiondict['xlabel'],fontsize=20)
    ax1.set_ylabel(Optiondict['ylabel'],fontsize=20)  
#### TPSI_RIOT用  
    #ax1.set_ylim([-0.1,0.8])
    #ax1.set_xlim([75,170])
    #ax1.set_ylim([-0.7,0.7])

    #ax1.set_xticks([ 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]) #[0, 0.2, 0.4, 0.6]
    #ax1.set_xticklabels(['0.25', '0.50', '0.75', '1.00', '1.25', ''])#, rotation=30, fontsize='small') #'0.0', '0.2', '0.4', '0.6'
    #ax1.set_yticks([0, 0.2, 0.4, 0.6])
    #ax1.set_yticklabels(['0.0', '0.2', '0.4', '0.6'])#, rotation=30, fontsize='small')
### テンソル分解用
    #ax1.set_xlim([-0.3,-0.16])
    #ax1.set_xlim([-0.3,-0.016])

    #y=xの直線描画
    xmin, xmax, ymin, ymax = ax1.axis() 
    x = np.arange(xmin, 1, 0.1)
    #ax1.plot(x,x,c='red',linewidth=1)
    
    #y軸のmax、minを求める
    ymin=min(list2); ymax=max(list2)
    xmin=min(list1); xmax=max(list1)
    try:
        if Optiondict['calcR'] in ['pearson','spearman']:
            #ax1.set_ylim([-0.1, 1.1])
            Roundp = round(p,100)
            Roundr = round(r,3)
            b = '%.2e'%Roundp
            TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
            p='{:e}'.format(p)
            r='{:e}'.format(r)
            TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
            TenRoundr = str(r)[0:4] + '$\it{×10^{-'+str(r)[str(r).find('e')+2:]+'}}$'

            #ax1.annotate('$\it{R}$ = ' + str(round(r,3)) + ', $\it{p}$ = ' + TenRoundp,fontsize=12, xy=(120,0.6))
            #ax1.annotate('$\it{r}$ = ' + TenRoundr + ', $\it{p}$ =' +TenRoundp,fontsize=12, xy=(-0.05,0.64))
            ax1.annotate('$\it{r}$ = ' + TenRoundr + ', $\it{p}$ =' +TenRoundp,fontsize=20, xy=(-730,-10.5))

            #ax1.annotate('$\it{R}$ = ' + TenRoundr + ', $\it{p}$ = ' + TenRoundp,fontsize=20, xy=(xmin+(xmin/100)+80,ymax-(ymax/100)+0.02))    
            #plt.title('$\it{R}$ = ' + TenRoundr + ', $\it{p}$ = ' + TenRoundp,fontsize=10,loc='upper')#, xy=(xmin+(xmin/100)+80,ymax-(ymax/100)+0.02))    
            #print('7')
            #ax2.set_title('$\it{R}$ = ' + TenRoundr + ', $\it{p}$ = ' + TenRoundp,fontsize=20)#, xy=(xmin+(xmin/100)+80,ymax-(ymax/100)+0.02))    
            
        else:
            r=0;p=0   

    except:
        print('44')
        r=0;p=0  
       

    try:
        if Optiondict['mkScatterWHist_drawEllipse']==1:#散布図とヒストグラムのために楕円フィットするなら
            NewDF = Optiondict['NewDF'];groups = Optiondict['groups']; MolLabel = Optiondict['MolLabel']
            drawEllipse(ax1,NewDF,groups, MolLabel,'TimeVar','AveCorr')#2つのリストの散布図+ヒストグラム
    except:
        pass
    if Optiondict['Annotate'] == 1:#分子名などをAnnotateするなら
        AnDict=dict({'Glucose':'Glc','Insulin':'Ins','C-peptide':'CRP','GIP(Active)':'GIP','Pyruvate':'Pyr','Total bile acid':'TBA',
                   'Citrate':'Cit','Cortisol':'Cor','Free fatty acid':'FFA','Total ketone body':'Ketone','Glutamic acid':'Glu',
                   'Citrulline':'Citr','Methionine':'Met','Isoleucine':'Ile','Leucine':'Leu','Tyrosine':'Tyr','4-Methyl-2-oxopentanoate':'4M2O','Glu + threo-beta-methylaspartate':'Glu+TBM','Growth hormone':'GH'})
        for i in range(len(list1)):#list(MolCorr.columns) #分子名などをAnnotateするなら
            if Optiondict['Label'][i] in ['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid',
                   'Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid',
                   'Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate','Growth hormone']:
                ax1.annotate(AnDict[Optiondict['Label'][i]],fontsize=0.1, xy=(list1[i],list2[i]))

            else:
                try:
                    ax1.annotate(Optiondict['Label'][i],fontsize=Optiondict['LabelSize'], xy=(list1[i],list2[i]))
                    
                except:
                    ax1.annotate(Optiondict['Label'][i],fontsize=7, xy=(list1[i],list2[i]))
        RankTimeVarDF.index=Optiondict['Label']
    try:        
        if Optiondict['y=x'] == 1:#y=xの線を足す
                x = np.linspace(xmin,xmax)  # xの値域(0, 1, 2, 3)
                y = x             # 直線の式
                
                ax1.plot(x,y,"r-",linewidth=1)#,label='y=x') 
    except:
        pass
    # ax1.set_xlim([-0.28,-0.17])
### temp_20200219    
    #plt.ylim([-5,5])
    title = Optiondict['title']+Optiondict['calcR']

    RankTimeVarDF.to_excel(save_dir + 'DF_'+title+'.xlsx')

    plt.savefig(save_dir +'Scatter_'+title+'.pdf',format='pdf',bbox_inches="tight")
    plt.savefig(save_dir +'Scatter_'+title+'.png',format='png',bbox_inches="tight")

    return(r,p)
def mkSortedBarAveStdWHist(List1, stdList, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir):#ヒストグラム月の任意のリストを代入して棒グラフを描画する
    #縦横などの微調整はその都度行う。
    
    x = np.linspace( 1, len(List1), len(List1) )
    
    #fig.add_axes((左下, 左上, x長さ, y長さ)x座標、y座標、幅、高さ
    """
    #縦長なら
    fig = plt.figure(figsize=(12,16))
    ax1 = fig.add_axes((0, 0.8, 1, 0.8))
    ax2 = fig.add_axes((0, 1.6, 1, 0.2), sharex=ax1)
    ax1.barh( x, List1, color=Color,tick_label=xticks,linewidth=1)
    """
    #横長なら
    fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_axes((0, 1, 0.8, 1))
    ax2 = fig.add_axes((0.8, 1, 0.2, 1), sharey=ax1)
    width=0.8
    
    ax1.bar( x, List1, width, yerr=stdList, color=Color,tick_label=xticks,linewidth=0.5,ec='black')
    ax2.hist(List1, bins=20,orientation="horizontal",ec='black')

    #plt.bar(bars, heights, width, tick_label=label, yerr=std,
     #   align='center', alpha=0.5, ecolor='black', capsize=5)
    plt.title(Title,fontsize=Titlesize)
    ax1.set_xlabel(xlabel,fontsize=xsize)
    ax1.set_ylabel(ylabel,fontsize=size)
    xmin, xmax, ymin, ymax = ax1.axis() 
    #ax1.plot([xmin, xmax],[1, 1], "black", linestyle='dashed') # normal way
    #plt.xticks(x_tick)
    #ax1.set_xticks(rotation=270,fontsize=xsize)
    ax1.tick_params(axis='y',labelsize=40);    ax1.tick_params(axis='x',rotation=270,labelsize=xsize)


    #ax1.set_yticks(fontsize=size)
    ax2.tick_params(labelbottom="off",bottom="off") # x軸の削除
    ax2.tick_params(labelleft="false",labelsize=15);ax2.tick_params(axis='y',color='white',direction='out')

    #plt.gca().tick_params(axis='both',labelsize=20)

    plt.savefig(save_dir+Title+'_Sort.pdf',bbox_inches="tight")
    plt.savefig(save_dir+Title+'_Sort.png',bbox_inches="tight")

    plt.close()    

def mkSortedAngleBarWHist(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir):#ヒストグラム月の任意のリストを代入して棒グラフを描画する
    #縦横などの微調整はその都度行う。
    
    x = np.linspace( 1, len(List1), len(List1) )
    
    #fig.add_axes((左下, 左上, x長さ, y長さ)x座標、y座標、幅、高さ
    """
    #縦長なら
    fig = plt.figure(figsize=(12,16))
    ax1 = fig.add_axes((0, 0.8, 1, 0.8))
    ax2 = fig.add_axes((0, 1.6, 1, 0.2), sharex=ax1)
    ax1.barh( x, List1, color=Color,tick_label=xticks,linewidth=1)
    """
    #横長なら
    fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_axes((0, 1, 0.8, 1))
    ax2 = fig.add_axes((0.8, 1, 0.2, 1), sharey=ax1)
    ax1.bar( x, List1, color=Color,tick_label=xticks,linewidth=1,ec='black')
    ax2.hist(List1, bins=20,orientation="horizontal",ec='black')

    plt.title(Title,fontsize=Titlesize)
    ax1.set_xlabel(xlabel,fontsize=xsize)
    ax1.set_ylabel(ylabel,fontsize=size)
    xmin, xmax, ymin, ymax = ax1.axis() 
    #ax1.plot([xmin, xmax],[1, 1], "black", linestyle='dashed') # normal way
    #plt.xticks(x_tick)
    #ax1.set_xticks(rotation=270,fontsize=xsize)
    ax1.tick_params(axis='y',labelsize=40);    ax1.tick_params(axis='x',rotation=270,labelsize=xsize)


    #ax1.set_yticks(fontsize=size)
    ax2.tick_params(labelbottom="off",bottom="off") # x軸の削除
    ax2.tick_params(labelleft="false",labelsize=15);ax2.tick_params(axis='y',color='white',direction='out')
    ax1.set_yticks(np.arange(-1,1.1,0.5)*math.pi)
    ax1.set_yticklabels([r"$-\pi$",r"$-\frac{1}{2}\pi$",0,r"$\frac{1}{2}\pi$",r"$\pi$",
              r"$\frac{3}{2}\pi$",r"$2\pi$"])
    #plt.gca().tick_params(axis='both',labelsize=20)

    plt.savefig(save_dir+Title+'_Sort_Angle.pdf',bbox_inches="tight")
    plt.savefig(save_dir+Title+'_Sort_Angle.png',bbox_inches="tight")

    plt.close()
    
def mkSortedBarWHist(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, OptionDict, save_dir):#ヒストグラム付きの任意のリストを代入して棒グラフを描画する
    #縦横などの微調整はその都度行う。List1, Title, '', LabelDict['ylabelbar'], MolColorsCV, XtickLabel, 60, 30,0, save_dir
    
    x = np.linspace( 1, len(List1), len(List1) )
    
    #fig.add_axes((左下, 左上, x長さ, y長さ)x座標、y座標、幅、高さ
    
    #縦長なら
    fig = plt.figure(figsize=(12,16))
    ax1 = fig.add_axes((0, 0.8, 1, 0.8))
    ax2 = fig.add_axes((0, 1.6, 1, 0.2), sharex=ax1)
    ax1.barh( x, List1, color=Color,tick_label=xticks,linewidth=1)
    """
    #横長なら
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_axes((0, 1, 0.8, 1))
    ax2 = fig.add_axes((0.8, 1, 0.2, 1), sharey=ax1)
    ax1.bar( x, List1, color=Color,tick_label=xticks,linewidth=1,ec='black')
    ax2.hist(List1, bins=20,orientation="horizontal",ec='black')
    """
    ### xticksの色を変えたい
    #label = plt.ylabel("y-label")
    #label.set_color("red")
    ll=ax1.get_xticklabels() 
    #[i.set_color('blue') if '_B' in i else 'green' for i in ax1.get_xticklabels() ]
    #### 20190822_temp　いらないから?
    try:
        [t.set_color(i) for (i,t) in zip(OptionDict['xtickxolor'],ax1.xaxis.get_ticklabels())]
    except:
        pass
    
    try:
        if type(OptionDict['plotStar'])==list:#特定の分子に☆印をつける
            fs=30
            kwargs = dict(ha='center', va='bottom')
            if fs is not None:
                kwargs['fontsize'] = fs
            for ii in range(len(OptionDict['plotStar'])):
                lindx = xticks.index(OptionDict['plotStar'][ii])
                text='*'
                if 'B-C' in Title: #ばぐ？ 
                    lindy = List1[lindx] -6 #+ (np.max(np.abs(List1)))*0.02*List1[lindx]/np.abs(List1[lindx])
                else:
                    lindy = List1[lindx]
                lindx+=1
                ax1.text(lindx,lindy, text, **kwargs)
    except:
        pass
    plt.title(Title,fontsize=Titlesize)
    ax1.set_xlabel(xlabel,fontsize=xsize)
    ax1.set_ylabel(ylabel,fontsize=size)
    xmin, xmax, ymin, ymax = ax1.axis() 
    #ax1.plot([xmin, xmax],[1, 1], "black", linestyle='dashed') # normal way
    #plt.xticks(x_tick)
    #ax1.set_xticks(rotation=270,fontsize=xsize)
    ax1.tick_params(axis='y',labelsize=40);    ax1.tick_params(axis='x',rotation=270,labelsize=xsize)

    #ax1.set_yticks(fontsize=size)
    ax2.tick_params(labelbottom="off",bottom="off") # x軸の削除
    ax2.tick_params(labelleft="off",labelsize=15);ax2.tick_params(axis='y',color='white',direction='out')

    #plt.gca().tick_params(axis='both',labelsize=20)
### 軸の最大最小
    plt.ylim([0,1])
    
    plt.savefig(save_dir+Title+'_Sort.pdf',bbox_inches="tight")
    plt.savefig(save_dir+Title+'_Sort.png',bbox_inches="tight")

    plt.close()
    
def mkSortedBar(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize, save_dir):#任意のリストを代入して棒グラフを描画する
    
    
    x = np.linspace( 1, len(List1), len(List1) )
    fig = plt.figure(figsize=(16,16))

    plt.bar( x, List1,width=1, color=Color,tick_label=xticks,linewidth=1,ec='black')   
    
    plt.title(Title,fontsize=Titlesize)
    plt.xlabel(xlabel,fontsize=xsize)
    plt.ylabel(ylabel,fontsize=size)
    #plt.xticks(x_tick)
    plt.xticks(rotation=270,fontsize=xsize)
    plt.yticks(fontsize=size)
    plt.savefig(save_dir+Title+'_Sort.pdf',bbox_inches="tight")
    plt.close()

def mkBarAveStd(List1, StdList, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize,save_dir):#任意のリストを代入して棒グラフを描画する
    
    
    x = np.linspace( 1, len(List1), len(List1) )
    fig = plt.figure(figsize=(20,12))
#error_kw=dict(lw=5, capsize=5, capthick=3)
    plt.bar( x, List1,width=1,yerr = StdList, color=Color,tick_label=xticks,linewidth=0.5,ec='black',error_kw=dict(lw=1, capsize=1, capthick=1))
    
    
    plt.title(Title,fontsize=Titlesize)
    plt.xlabel(xlabel,fontsize=xsize)
    plt.ylabel(ylabel,fontsize=size)
    #plt.xticks(x_tick)
    plt.xticks(rotation=270,fontsize=xsize) 
    plt.yticks(fontsize=size)    
    xmin, xmax, ymin, ymax = plt.axis() 
    #plt.plot([xmin, xmax],[1, 1], "black", linestyle='dashed') # normal way
    plt.savefig(save_dir+Title+'_.pdf',bbox_inches="tight")
    plt.close()
    
def mkBar(List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize,save_dir):#任意のリストを代入して棒グラフを描画する
    
    
    x = np.linspace( 1, len(List1), len(List1) )
    
    fig = plt.figure(figsize=(5,3))
    ax1 = fig.add_axes((0, 1, 1, 1))

    ax1.bar( [0+0.1*i for i in range(len(xticks))], List1,width=0.05, color=Color,tick_label=xticks,linewidth=0.5,ec='black')
    
    ax1.set_ylim([-1.0,1.0])
    
    xmin, xmax, ymin, ymax = ax1.axis() 
    
    p = ax1.hlines([0], xmin, xmax, "black", linestyles='solid',linewidth=1)     # hlines
    
    axis=['top','bottom','left','right']
    line_width=[1,1,1,1]
    
    for a,w in zip(axis, line_width):  # change axis width
        ax1.spines[a].set_linewidth(w)

    ax1.set_xlim(xmin, xmax)
    ax1.set_title(Title,fontsize=Titlesize)
    ax1.set_xlabel(xlabel,fontsize=xsize)
    ax1.set_ylabel(ylabel,fontsize=size)
    #plt.xticks(x_tick)
    ax1.set_xticklabels(labels=xticks,rotation=270,fontsize=xsize) 
    #ax1.set_xticklabels(fontsize=size)    
    xmin, xmax, ymin, ymax = ax1.axis() 
    #plt.plot([xmin, xmax],[1, 1], "black", linestyle='dashed') # normal way
    plt.savefig(save_dir+Title+'_.pdf',bbox_inches="tight")
    plt.close()


def DFScatter(DF,OptionDict,numcol,numrow,fsize,MolColor,save_dir):#DFで受け取って列の分だけ散布図描画
     #相関取りやすくするため、リスト化
    #FastCVList = list( FastCV.iloc[0,:] ); MolTimeAveList=list( MolTimeAve.iloc[0,:] ); MolTimeCVAVeList=list( MolTimeAve.iloc[0,:] ); MolTimeAveList=list( MolTimeAve.iloc[0,:] ); MolTimeAveList=list( MolTimeAve.iloc[0,:] )
    #MolName = list(FastCV.columns)
   
    Col = list(DF.columns);Idx=list(DF.index)
    indnames=Col
    fig,ax = plt.subplots(1,1)
    
 
    fig,host = plt.subplots(numrow,numcol,figsize=(15,15))
    
    for i in range(0,len(indnames)):#列
        for j in range(0+i,len(indnames)):#行    
        
            if i==j:#同じプロパティのところではヒストグラム
                    
                    host[j,i].hist(DF[indnames[i]][~np.isnan(DF[indnames[i]])])
                    if i == 0: 
                       host[j,0].set_ylabel(indnames[j],fontsize=fsize)  
                    host[numrow-1,i].set_xlabel(indnames[i],fontsize=fsize)
                    host[numrow-1,i].xaxis.set_label_coords(0.5, -0.3)
                    host[j,i].tick_params(labelsize=fsize)
                    host[j,i].tick_params(width = 0.5, length = 1)
            else:
                mask = ~np.logical_or(np.isnan(DF[indnames[i]]), np.isnan(DF[indnames[j]]))
                x, y = DF[indnames[i]][mask], DF[indnames[j]][mask]
                try:
                    r, p = pearsonr(x, y)  
                except:
                    r, p =0,0
                host[j,i].scatter(DF[indnames[i]],DF[indnames[j]],s=8)#, color=MolColor)
                host[j,i].tick_params(labelsize=fsize)
                if i == 0: 
                    host[j,0].set_ylabel(indnames[j],fontsize=fsize)
                host[numrow-1,i].set_xlabel(indnames[i],fontsize=fsize)
                Roundp = round(p,200)
                b = '%.2e'%Roundp
                TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'        
                host[j,i].set_title('R=' + str(round(r,3)) + ', p=' + TenRoundp,fontsize=fsize)
                host[j,i].tick_params(width = 0.5, length = 1)
        
                 #TenRoundp = b[0:b.find('e')] + '$\it{×10^{-'+str(b[len(b)-1])+'}}$'
                  #      p = Pvalue[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]]
                   #     p='{:e}'.format(p)
                        
                    #    TenRoundp = str(p)[0:4] + '$\it{×10^{-'+str(p)[str(p).find('e')+2:]+'}}$'
                     #   q = FastPropDF[FastPropDF.columns[UndrThsh[1][count]]][FastPropDF.index[UndrThsh[0][count]]]
                      #  q='{:e}'.format(q)
                       # TenRoundq = str(q)[0:4] + '$\it{×10^{-'+str(q)[str(q).find('e')+2:]+'}}$'

    host[0,1].axis("off");host[0,2].axis("off");host[0,3].axis("off");host[0,4].axis("off");host[0,5].axis("off");
    host[1,2].axis("off");host[1,3].axis("off");host[1,4].axis("off");host[1,5].axis("off");
    host[2,3].axis("off");host[2,4].axis("off");host[2,5].axis("off");
    host[3,4].axis("off");host[3,5].axis("off");
    host[4,5].axis("off");
    #Yaxismin = min(MolTimeAveList); Yaxismax=max(MolTimeAveList)
    
    #plt.xticks(MolName, rotation=270)
    #plt.yticks(MolName, rotation=270)    
    #plt.grid()

    #plt.ylim([0,1])
    #plt.yticks( np.arange(0, 1, 0.1) )
    #plt.tight_layout()
    #plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.3f'))#y軸小数点以下3桁表示
    fig.tight_layout()
    plt.savefig(save_dir+'AllCorrSubjVariation.pdf')#,bbox_inches="tight")
    
    
def ScatterWColorbar(x,y,value,ymin,ymax,save_dir):
    import math
    import numpy as np
    from scipy.stats import gaussian_kde        
    #x = np.random.rand(100)
    #y = np.random.rand(100)
    # 乱数を 100 件生成
    """
    bin_num = math.sqrt(len(x))#便の数
    bins = np.linspace(min(x), max(x), bin_num)
    bin_indice = np.digitize(x, bins)
    value = bin_indice#np.random.rand(len(x))
    print(value)
    """
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    #print(z)
    value=(z - min(z))/(max(z) - min(z))
    fig = plt.figure(figsize=(5,5))
    # 散布図を表示
    plt.scatter(x, y, s=10, c=value,vmin=0,vmax=1, cmap='jet')
    #plt.axis([-0.05,1.05,0,0.7])#ここは好きに変えて良い

    #plt.axis([-0.05,1.05,0,0.7])
    plt.ylim([ymin, ymax])
    #plt.xlim([0,0.99])
    # カラーバーを表示
    plt.colorbar()
    plt.savefig(save_dir+'hiatmap_color.pdf')
    plt.close()
    
def mkHist(DF,filename,save_dir):#ヒストグラムを作る
    #DFの各列ごとのヒストグラム
    fig, ax = plt.subplots(1,1,figsize=(6.4,4.8))
    row = list(DF.index)
    col = list(DF.columns)
    ax= DF.plot(y=col, bins=50, alpha=0.5, figsize=(16,4),kind='hist')
    plt.savefig(save_dir + 'Hist' +filename + '.pdf')

def mkhist2List(List1,List2, save_dir,xlabel,filename,bins):#2つのリストのヒストグラム
        fig = plt.figure(figsize=(6.4,4.8))
        plt.hist(List1,bins=bins)#,cumulative=1,normed=True)#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)        
        plt.hist(List2,bins=bins)#,cumulative=1,normed=True)#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)        

        plt.xlabel(xlabel,fontsize=10); plt.ylabel('Frequency',fontsize=10);plt.tick_params(labelsize=10)
        plt.savefig(save_dir +  filename +'Hist_RDFrev.pdf')  ;plt.close()#print(ConList2[4][~np.isnan(ConList2[4])])

    
def mkhistList(List1,save_dir,xlabel,filename,bins,title):#1つのリストのヒストグラム
        fig = plt.figure(figsize=(6.4,4.8))
        plt.hist(List1,bins=bins,ec='black',align='left' )#,cumulative=1,normed=True)#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)  
        #alin 各棒の中心を X 軸目盛上のどの横位置で出力するか。 ‘left‘,‘mid‘,‘right‘ から選択。デフォルト値: ‘mid’
        #plt.xlim([-4,2]) #モデルパラメタのとき
        plt.xlabel(xlabel,fontsize=20); plt.ylabel('Frequency',fontsize=20);plt.tick_params(labelsize=20);plt.title(title,size=20)
### temp 縦軸をlog10に
        plt.yscale('log');#plt.xticks(np.arange(0, 110 + 1, 10)) #np.arange(-110, 20 + 1, 20) np.arange(0, 110 + 1, 10)
        plt.savefig(save_dir +  filename +'Hist.pdf',bbox_inches="tight")  ;plt.close()#print(ConList2[4][~np.isnan(ConList2[4])])

def mkhistListWline(List1,save_dir,xlabel,filename,bins,OptionDict,title):#1つのリストのヒストグラム with 線を引ける
        fig = plt.figure(figsize=(6.4,4.8))
        
        if OptionDict['drawLine']==1:
            Thresh =OptionDict['Thresh']
            ThreshValue = np.percentile(List1,Thresh)
            plt.hist(List1,bins=bins,ec='black',align='left' )#,cumulative=1,normed=True)#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)  
                #xmin, xmax, ymin, ymax = ax1.axis() 
                #x = np.arange(xmin, 1, 0.1)
            plt.axvline(x=ThreshValue, ymin=0, ymax=1,lw=1,c='r')
        else:
            plt.hist(List1,bins=bins,ec='black',align='left' )#,cumulative=1,normed=True)#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)  

        #alin 各棒の中心を X 軸目盛上のどの横位置で出力するか。 ‘left‘,‘mid‘,‘right‘ から選択。デフォルト値: ‘mid’
        plt.xlim([0,1]) #モデルパラメタのとき
        plt.xlabel(xlabel,fontsize=20); plt.ylabel('Frequency',fontsize=20);plt.tick_params(labelsize=20);plt.title(title+'_r(p=0.01)='+str(round(ThreshValue,4)),size=20)
### temp 縦軸をlog10に
        #plt.yscale('log');#plt.xticks(np.arange(0, 110 + 1, 10)) #np.arange(-110, 20 + 1, 20) np.arange(0, 110 + 1, 10)
        plt.savefig(save_dir +  filename +'HistWLine.pdf',bbox_inches="tight")  ;plt.close()#print(ConList2[4][~np.isnan(ConList2[4])])


def mkhistabs(SubjRmeanrev,save_dir,xlabel,bins,filename):#DF中全てのヒストグラム絶対値
        fig = plt.figure(figsize=(6.4,4.8))
        list1 = np.array(SubjRmeanrev).flatten()[~np.isnan(SubjRmeanrev.values.flatten().astype(np.float32))]
        ax1 = plt.hist(np.abs( list1.astype(np.float32) ) ,bins=bins)#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)        
        plt.xlabel(xlabel,fontsize=20); plt.ylabel('Frequency',fontsize=20);plt.tick_params(labelsize=20)
        plt.savefig(save_dir +  filename +'Hist_RDFrev.pdf')  ;plt.close()#print(ConList2[4][~np.isnan(ConList2[4])])

    
def mkhist(SubjRmeanrev,save_dir,xlabel,bins, filename):#DF中全てのヒストグラム
        fig = plt.figure(figsize=(6.4,4.8))
        list1 = np.array(SubjRmeanrev).flatten()[~np.isnan(SubjRmeanrev.values.flatten().astype(np.float32))]
        ax1 = plt.hist(list1.astype(np.float32),bins=bins, align='mid')#,range=(11, 20)  )#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)        
        #plt.xlim([11,20])
        plt.xlabel(xlabel,fontsize=20); plt.ylabel('Frequency',fontsize=20);plt.tick_params(labelsize=20)
        plt.savefig(save_dir +  filename +'Hist_RDFrev.pdf')  ;plt.close()#print(ConList2[4][~np.isnan(ConList2[4])])

    
def ScatterWHeatmap(x,y,xmin,xmax,ymin,ymax,save_dir,bin):#散布図をヒートマップで
    fig,ax = plt.subplots(1,1)
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bin,normed=True)#x ,y があれば色付けできる
    extent = [xmin,xmax,ymin,ymax]#xmin,xmax,ymin,ymaxで描画する    
    plt.clf()
    im=plt.imshow(heatmap.T, extent=extent,interpolation='nearest',origin='lower', cmap=plt.get_cmap('jet'))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size="2%", pad=0.1)#
    #fig.add_axes(cax)
    fig.colorbar(im,cax=cax)
    #fig.tight_layout()
    #plot_axis = plt.axis()
    plt.savefig(save_dir+'.pdf')
    plt.show()
    plt.close()
    
def  plotAveStdMultiIndex(MultiIdxDF,save_dir):#各量で平均vs標準偏差_MultiIndexver.
    # sklearn.linear_model.LinearRegression クラスを読み込み
    from sklearn import linear_model
    clf = linear_model.LinearRegression()
    #fig,host=plt.subplots(numFeature,numFeature,figsize=(20,20))
    MolColor = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
    MolColorList = MolColor['MolColor']

    MolList = list(MultiIdxDF.columns)
    IdxList = list(set(MultiIdxDF.index.get_level_values(level=0)))#Index   
    AveList=[];StdList=[]
    #MolList=['Glucose','Insulin','C-peptide','GIP(Active)','Pyruvate','Total bile acid','Citrate','Cortisol','Free fatty acid','Total ketone body','Glutamic acid','Citrulline','Methionine','Isoleucine','Leucine','Tyrosine','4-Methyl-2-oxopentanoate','Glu + threo-beta-methylaspartate']#,'Growth hormone']]

    for j in IdxList:#各指標
        for i in MolList: #各分子   
            AveList.append(MultiIdxDF.loc[j][i].mean())
            StdList.append(MultiIdxDF.loc[j][i].std(ddof=0) )
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
        r, p = pearsonr(x, y)   
        
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

 
if __name__ == "__main__": 

    
    # Generate some test data
    x = np.random.randn(5000)
    y = np.random.randn(10000)
    
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=50)#x ,y があれば色付けできる
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]#xmin,xmax,ymin,ymaxで描画する
    
    plt.clf()
    plt.imshow(heatmap.T, extent=extent, origin='lower',cmap=plt.get_cmap('jet'))
    plt.show()
"""

    fig,ax = plt.subplots(1,1,figsize=(16,5))
    heatmap, xedges, yedges, Image = ax.hist2d(x, y, bins=bin,normed=True,cmap=cm.jet)#x ,y があれば色付けできる
    extent = [xmin,xmax,ymin,ymax]#xmin,xmax,ymin,ymaxで描画する    
    plt.clf()
    plt.imshow(heatmap.T, extent=extent, cmap=plt.get_cmap('jet'))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    fig.add_axes(cax)
    fig.colorbar(Image,cax=cax)
    fig.tight_layout()
    plot_axis = plt.axis()
    plt.savefig(save_dir+'.pdf')
    plt.show()
 
    
    
    
    import seaborn as sns
    from numpy.random import *
    pts = [] #(x,y)の点を格納する配列　←　これに自分のデータを入れればok
     
    #テストデータとして、点をランダムに生成
    num_points=1000 #ランダムに生成する点の数
    for x in range(num_points):
        X=randn()
        Y=randn()
        pts.append((X,Y)) #正規分布の乱数
 
    monitor_size=(1024,1024) #画面サイズ
    sns.pairplot(pts)
    hm = sns.heatmap(pts)
    img = hm.heatmap(
        points=pts, #点の配列(x,y)
        size=monitor_size,#画面サイズ
        dotsize=30, #ガウシアンのσ(標準偏差)に対応
        #area=((0, 0),(monitor_size[0], monitor_size[1])), #座標の設定
        scheme='classic', #スタイル(classic,fire    ,omg,pbj,pgaitch)
        opacity=200 #透明度
        )
    img.save("classic.png")
"""