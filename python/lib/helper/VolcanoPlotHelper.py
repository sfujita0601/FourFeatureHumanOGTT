#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 14:32:29 2018

@author: fujita
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.style.use('ggplot')


def DrawPieChart(LabelcountDF,save_dir):#増加と減少、またはどちらもの円グラフ描画
    LabelSum =   pd.read_excel("//Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx",header=0,index_col=0).iloc[:83,:]
    Up =[]; Down =[]; Both=[]; 
    for i in list(LabelcountDF.index):
        #if (LabelcountDF.loc[i]['Down'] > 0) and (LabelcountDF.loc[i]['Up'] > 0):
         #   Both.append(i)
        if LabelcountDF.loc[i]['Down'] > 0:
            Down.append(i)
        elif LabelcountDF.loc[i]['Up'] > 0:
            Up.append(i)
    Others = len(LabelSum.index) - len(Up) - len(Down) - len(Both)
    
    #label = list(LabelSum['English'])
    
    #label = [ '増加', '減少', '両方', '変化なし' ]
    #label = [ 'Increased', 'Decreased', 'Both', 'Not Changed' ]
    label = [ 'Increased', 'Decreased',  'Not Changed' ]

    colorList = [ 'red', 'blue', 'purple', 'gray']
    colorList = [ 'red', 'blue',  'gray']
    colorList = [ 'sienna', 'indigo',  'gray']

    #aa,labels, texts = plt.pie(list ( [len(Up) ,len(Down) ,len(Both) ,Others ] ), labels=label,startangle=90,autopct="%.1f%%",counterclock=False, colors=colorList, pctdistance=0.75)#, textprops={'color': "white", 'weight': "bold"})
    aa,labels, texts = plt.pie(list ( [len(Up) ,len(Down)  ,Others ] ), labels=label,startangle=90,autopct="%.1f%%",counterclock=False, colors=colorList, pctdistance=0.75)#, textprops={'color': "white", 'weight': "bold"})

    #plt.axis('equal');
    #plt.title('Ratio of each metabolism',size=15)
    for t in texts:
      t.set_horizontalalignment('center')
      t.set_color('white')
      t.set_weight("bold")
      t.set_size(15)
    for a  in labels: 
      a.set_size(15)
    plt.axis('equal')
    plt.savefig(save_dir + ' piechart_UpDown.pdf' )#,bbox_inches="tight")


class Volcano(object):
    """
    create a Volcano plot from log2(ratios) and corresponding -log10(p_values)
    ToDo: take care of infinite ratios
    
    e.g. usage
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.style.use('ggplot')
    # when using jupyter notebook:
    # %matplotlib inline
    
    
    dict_ = {'p_val': {1: 0.208057931774703, 2: 0.063320586966883294, 3: 0.11424685505629198, 4: 0.46130291511303301, 5: 0.31662387800522196, 6: 0.35821379490648098, 7: 0.0559720700871537, 8: 0.048805611096553, 9: 0.27635717881946303, 336: 7.1364464137191392, 264: 0.50054528789819708, 530: 4.9027283450578603, 83: 4.4128053670050704, 565: 5.7615691096534203}, 'ratio': {1: -0.22443504333496, 2: 0.074678929646808001, 3: 0.27419026692708204, 4: 0.48245075770786605, 5: 0.274400329589845, 6: -0.53869597117106194, 7: 0.13976793819003699, 8: 0.20341746012369896, 9: 0.54393492804633292, 336: 4.2886840820312502, 264: -3.9391211668650299, 530: 4.5683480792575404, 83: 2.5429632663726802, 565: 4.04515149858263}, 's_val': {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0, 336: 11.159643571086351, 264: 0.0, 530: 6.9983089690857891, 83: 1.0852958656755318, 565: 7.2227303405181722}, 'label': {1: 'CRYAB', 2: 'HDLBP', 3: 'TAF15', 4: 'GNAO1', 5: 'KHSRP', 6: 'HSPA4', 7: 'PABPC1;PABPC4', 8: 'MINOS1', 9: 'SPCS2', 336: 'APOL2', 264: 'WISP3', 530: 'WARS', 83: 'HLA-C', 565: 'GBP1'}}
    df = pd.DataFrame(dict_)
    volc = Volcano(df["ratio"], df["p_val"], df["label"], s_curve_x_axis_overplot=0.5, s_curve_y_axis_overplot=0.5)
    fig = volc.get_fig()
    
    for S-curve calculation see Jan Christian Refsgaard's PhD thesis page 90-91
    from Perseus documentation:
    S0: Artificial within groups variance (default: 0). 
    It controls the relative importance of t-test p-value and difference between means. 
    At s0=0 only the p-value matters, while at nonzero s0 also the difference of means 
    plays a role. See Tusher, Tibshirani and Chu (2001) PNAS 98, pp5116-21 for details.
    """
    dict_ = {'p_val': {1: 0.208057931774703, 2: 0.063320586966883294, 3: 0.11424685505629198, 4: 0.46130291511303301, 5: 0.31662387800522196, 6: 0.35821379490648098, 7: 0.0559720700871537, 8: 0.048805611096553, 9: 0.27635717881946303, 336: 7.1364464137191392, 264: 0.50054528789819708, 530: 4.9027283450578603, 83: 4.4128053670050704, 565: 5.7615691096534203}, 'ratio': {1: -0.22443504333496, 2: 0.074678929646808001, 3: 0.27419026692708204, 4: 0.48245075770786605, 5: 0.274400329589845, 6: -0.53869597117106194, 7: 0.13976793819003699, 8: 0.20341746012369896, 9: 0.54393492804633292, 336: 4.2886840820312502, 264: -3.9391211668650299, 530: 4.5683480792575404, 83: 2.5429632663726802, 565: 4.04515149858263}, 's_val': {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0, 336: 11.159643571086351, 264: 0.0, 530: 6.9983089690857891, 83: 1.0852958656755318, 565: 7.2227303405181722}, 'label': {1: 'CRYAB', 2: 'HDLBP', 3: 'TAF15', 4: 'GNAO1', 5: 'KHSRP', 6: 'HSPA4', 7: 'PABPC1;PABPC4', 8: 'MINOS1', 9: 'SPCS2', 336: 'APOL2', 264: 'WISP3', 530: 'WARS', 83: 'HLA-C', 565: 'GBP1'}}
    df = pd.DataFrame(dict_)
    #ラベルとp_val, FCratioを

    #df = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181102/TimePointSubjMol/DFForVolcanoLog.xlsx',header=0,encoding = "ISO-8859-1")
    #volc = Volcano(df["ratio"], df["p_val"], df["s_val"],df["label"], s_curve_x_axis_overplot=0.5, s_curve_y_axis_overplot=0.5)
    #fig = volc.get_fig()
    #plt.savefig('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181102/TimePointSubjMol/VolcanoPlot.pdf')
    def __init__(self, ratio, p_val, s_val,label=None, s_curve_x_axis_overplot=0.5, s_curve_y_axis_overplot=0.5):
        """
        careful use ratio not difference as in Perseus 
        ratio of 0.5 instead of difference -2
        :param ratio: Pandas.Series or Numpy.Array or List of log2(ratios)
        :param p_val: Pandas.Series or Numpy.Array or List of -log10(p-values) 
        :param label: Pandas.Series or Numpy.Array or ListOfString
        """
        assert len(ratio) == len(p_val)
        self.df = pd.DataFrame({"ratio": ratio, "p_val": p_val})
        if label is not None:
            self.df["label"] = label
        self.s_curve_y_axis_overplot = s_curve_y_axis_overplot
        self.p_val_cutoff = self.get_p_val_cutoff()
        self.ratio_cutoff = self.get_ratio_cutoff()
        #self.df["s_val"] = self.df.apply(self.calc_s_from_row, axis=1)
        self.df["s_val"] = s_val
        self.ratio_for_s = pd.Series(np.linspace(self.df["ratio"].min() - s_curve_x_axis_overplot, self.df["ratio"].max() + s_curve_x_axis_overplot, num=1000))
        self.p_for_s_larger_1 = self.ratio_for_s.apply(self.calc_p_for_s_equals_1)

    def get_p_val_cutoff(self):
        """
        p_val_cutoff = 0.05
        pc = 3.5 + median(p_val(50% lowest log2_ratios)) --> is what Jan uses for whatever reason ???
        -log10_pval of 2.0 --> pval of 0.01
        """
        ### hard coded cutoff of 1%
        # return math.log(0.01, 10) * -1
        ### Jan's cutoff, but how to justify???
        quant = self.df["ratio"].quantile(0.5)
        return 2.0 + self.df.loc[self.df["ratio"] < quant, "p_val"].median()

    def get_ratio_cutoff(self):
        """
        log2_ratio_cutoff = 2.0 
        ratio_cutoff_high = 2 + median(ratio(50% lowest log10_p_values))
        ratio_cutoff_low = 0.5 - median(ratio(50% lowest log10_p_values))        
        """
        ### hard coded cutoff of 2 fold enrichment or depletion
        # return math.log(0.5, 2), math.log(2, 2)
        ### Jan's cutoff, how to justify???
        quant = self.df["p_val"].quantile(0.5)
        median_ = self.df.loc[self.df["p_val"] < quant, "ratio"].median()
        ratio_cutoff_high = np.log2(1.5)#2.0 + median_#modified
        ratio_cutoff_low = np.log2(0.67)#-2.0 - median_#modified
        return ratio_cutoff_low, ratio_cutoff_high

    def calc_s_from_row(self, row):
        p_val = row["p_val"]
        ratio = row["ratio"]
        return self.calc_s(p_val, ratio)

    def calc_s(self, p_val, ratio):
        """
        so the algorithmn for finding stuff with s > 1 is:
        discard stuff below the ratio_cutoff
        discard stuff below the p-val cutoff
        do the calcuation for the stuff above BOTH cutoffs and accept all with s > 1
        s = (p_val - p_val_cutoff) * (ratio - ratio_cutoff)
        :param p_val: Float(-log10 p-value)
        :param ratio: Float(log2 ratio)
        :return: Float
        """
        ratio_cutoff_low, ratio_cutoff_high = self.ratio_cutoff
        if ratio > 0:
            ratio_delta = ratio - ratio_cutoff_high
            if ratio_delta < 0:
                return 0
        elif ratio < 0:
            ratio_delta = ratio - ratio_cutoff_low
            if ratio_delta > 0:
                return 0
        else: 
            ratio_delta=0
            return 0
        ratio_delta = abs(ratio_delta)
        p_val_delta = p_val - self.p_val_cutoff
        if p_val_delta < 0:
            return 0
        return p_val_delta * ratio_delta

    def calc_p_for_s_equals_1(self, ratio):
        """
        :param ratio: Float(log2 ratio)
        :return: Float
        """
        ratio_cutoff_low, ratio_cutoff_high = self.ratio_cutoff
        ratio_delta_high = ratio - ratio_cutoff_high
        ratio_delta_low = ratio - ratio_cutoff_low

        if ratio > ratio_cutoff_high:
            return (1.0 / ratio_delta_high) + self.p_val_cutoff
        elif ratio < ratio_cutoff_low:
            return (1.0 / (ratio_delta_low * -1)) + self.p_val_cutoff
        else:
            return np.nan

    def get_fig(self, title="Volcano plot", s_value_cutoff=1000,ColorSwich =2,Labelswitch='labeltable') :#色を変えられる、labeltable'で閾値超えたTable吐く, plotlabel:add mol label on fig
        fig, ax1 = plt.subplots(figsize=(12, 12))
        ax1.set_title(title, fontsize=22, fontweight="bold")

        cond_p_val = self.p_for_s_larger_1 <= (self.df["s_val"].max() + self.s_curve_y_axis_overplot)
        cond_pos = self.ratio_for_s > 0
        x1 = self.ratio_for_s[cond_pos & cond_p_val]
        y1 = self.p_for_s_larger_1[cond_pos & cond_p_val]

        cond_neg = self.ratio_for_s < 0
        x2 = self.ratio_for_s[cond_neg & cond_p_val]
        y2 = self.p_for_s_larger_1[cond_neg & cond_p_val]

        x3 = self.df["ratio"]
        y3 = self.df["s_val"]
        
        
        ####色分け
        ref_file = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/'
        MolColorDF =  pd.read_excel(ref_file+ 'LabelSummary_Eng_No3HB.xlsx',index_col=0)
        NewColorList=[]    
        if ColorSwich ==1:#分子で色分け
            import re
            for label in self.df["label"]:
                r = re.compile("(.*)(_)(.*)"); d = r.search(label); NewColorList += [MolColorDF['MolColor'][d.group(3)] ]           
        elif ColorSwich ==2:#応答の方向で色分け
            cond = ( (self.df["ratio"] <= np.log2(0.67) )|( self.df["ratio"] >= np.log2(1.5)) ) & (self.df["s_val"] > -1*np.log10(0.1) ) #s_value_cutoff#：本来ならば閾値を計算する

            for ratio in range(len(self.df["ratio"])):
                if cond[ratio] ==1 :      
                    if self.df["ratio"][ratio]>0:
                        NewColorList +=['sienna']#red
                    else:
                        NewColorList +=['indigo']#blue
                else:
                     NewColorList +=['black']
        else:
            NewColorList='black'
        ####補助線
        xmin=np.min(x3); xmax=np.max(x3); ymin=np.min(y3); ymax = np.max(y3)
        print(ymin)
        ax1.set_ylim([0,ymax+1])
        ax1.plot([np.log2(0.67), np.log2(0.67)],[0, ymax+1], "black", linestyle='dashed' ,linewidth=1.5)
        ax1.plot([np.log2(1.5), np.log2(1.5)],[0, ymax+1], "black", linestyle='dashed' ,linewidth=1.5)
        ax1.plot([xmin,xmax],[1,1], "black", linestyle='dashed' ,linewidth=1.5)

        ax1.xaxis.set_tick_params(labelsize=20)
        ax1.yaxis.set_tick_params(labelsize=20) 
        plt.rcParams['axes.linewidth'] = 2.0# 軸の線幅edge linewidth。囲みの太さ

        DF=pd.DataFrame(data=None,columns=['Up','Down'])
        if 'plotlabel' in Labelswitch :  #in  self.df.columns:#閾値超えたらLabel表示
            cond = self.df["s_val"] > s_value_cutoff#：本来ならば閾値を計算する
            cond = ( (self.df["ratio"] <= np.log2(0.67) )|( self.df["ratio"] >= np.log2(1.5)) ) & (self.df["s_val"] > -1*np.log10(0.1) ) #s_value_cutoff#：本来ならば閾値を計算する
            
            LabelCountDict={}
            for index_, row in self.df[cond].iterrows():
                
                label = row["label"]               
                x_coord = row["ratio"]
                y_coord = row["s_val"]
                ax1.annotate(label, xy=(x_coord, y_coord), xycoords='data', xytext=(5, 5),
                    textcoords='offset points', arrowprops=dict(arrowstyle="-"), fontsize=3)
                LabelCountDict.update({label:x_coord}) 
                
        if Labelswitch == 'labeltable':#各分子で閾値超えた方向と数
            cond = self.df["s_val"] > s_value_cutoff#：本来ならば閾値を計算する
            cond = ( (self.df["ratio"] <= np.log2(0.67) )|( self.df["ratio"] >= np.log2(1.5)) ) & (self.df["s_val"] > -1*np.log10(0.1) ) #s_value_cutoff#：本来ならば閾値を計算する
            
            LabelCountDict={}
            for index_, row in self.df[cond].iterrows():
                
                label = row["label"]               
                x_coord = row["ratio"]
                y_coord = row["s_val"]

                LabelCountDict.update({label:x_coord}) 
                #別の関数で各分子のupperListとlowerListを作成
                DF = self.MolUpperLowerCount(LabelCountDict)
                #tempDF=pd.concat([MolColorDF,DF],axis=1, join_axes=[MolColorDF.index])
                #DF = tempDF.dropna(axis=0)
                #DF.to_excel(save_dir+'LavelCount.xlsx')
               #'}}$'
               
        #ax1.plot(x1, y1, 'r-', x2, y2, 'r-',
        #'$\it{×10^{-'
        #for ij in range(len(x3)):
        ax1.scatter(x3, y3, color=NewColorList)#, size=10)  # , alpha=0.7)
        ax1.set_xlabel('log$\it{_{2}}$(Fold Change)', fontsize=40)
        ax1.set_ylabel('-log$\it{_{10}}$(q-value)', fontsize=40)
        plt.rcParams['font.family'] = 'Arial' #全体のフォントを設定 
        ax1.tick_params(axis='both',labelsize=30)

        return fig,DF
    
    def MolUpperLowerCount(self,Dict):#辞書を受けっとって、各分子で閾値を超えた数をそれぞれまとめる
        import re
        r = re.compile("(.*)(_)(.*)");
        #ラベルだけのListを作る
        TimetoMolLabel = list(Dict.keys() )
        LabelList = list( set([ r.search(i).group(3) for i in TimetoMolLabel ]) ) 
        NewDF = pd.DataFrame(data=None,index=LabelList,columns=['Up','Down'])

        for j in LabelList:#分子のラベル
            posList=[]
            negList =[]
            for i in TimetoMolLabel:#時点と分子名のラベル
            
                if (r.search(i).group(3) == j) and Dict[i] > 0:#log2FCが正なら
                    posList.append(Dict[i])
                elif (r.search(i).group(3) == j) and Dict[i] < 0:#負なら
                    negList.append(Dict[i])
            NewDF.loc[j,'Up'] = len(posList)
            NewDF.loc[j,'Down'] = len(negList)
        return(NewDF)
        
    def get_DistrFC(self,title='Distribution of Fold Change', bins=100):
        fig, ax1 = plt.subplots(figsize=(12, 12))
        ax1.set_title(title, fontsize=22, fontweight="bold")
        
        ratio = self.df["ratio"]

        xmin=np.min(ratio); xmax=np.max(ratio); #ymin=np.min(y3); ymax = np.max(y3)

        ax1.hist(ratio.values[~np.isnan(ratio.values)], range=(xmin, xmax), bins = bins, color='#1f77b4',#'C0',
                              ec='black');plt.title(title); #plt.savefig(save_dir + '_DistOfCorr_Wpearson'+filename+'.pdf');plt.close()
        ymin, ymax= ax1.get_ylim()
        print(ymax,ymin)
        ax1.plot([np.log2(0.67), np.log2(0.67)],[0, ymax+20], "black", linestyle='dashed' )   
        ax1.plot([np.log2(1.5), np.log2(1.5)],[0, ymax+20], "black", linestyle='dashed' )
        ax1.set_ylim([0,270.95])
        ax1.set_xlabel('log$\it{_{2}}$(FoldChange)', fontsize=40)
        ax1.set_ylabel('Frequency', fontsize=40)

        ax1.tick_params(axis='both',labelsize=30)
        ax1.xaxis.set_tick_params(labelsize=30)
        return(fig)