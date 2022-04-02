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



def DrawPieChart(LabelcountDF,file_dir,save_dir):
    LabelSum =   pd.read_excel(file_dir+"/LabelSummary.xlsx",header=0,index_col=0).iloc[:83,:]
    Up =[]; Down =[]; Both=[]; 
    for i in list(LabelcountDF.index):
        if LabelcountDF.loc[i]['Down'] > 0:
            Down.append(i)
        elif LabelcountDF.loc[i]['Up'] > 0:
            Up.append(i)
    Others = len(LabelSum.index) - len(Up) - len(Down) - len(Both)

    label = [ 'Increased', 'Decreased',  'Not Changed' ]
    colorList = [ 'sienna', 'indigo',  'gray']
    aa,labels, texts = plt.pie(list ( [len(Up) ,len(Down)  ,Others ] ), labels=label,startangle=90,autopct="%.1f%%",counterclock=False, colors=colorList, pctdistance=0.75)#, textprops={'color': "white", 'weight': "bold"})

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
    def __init__(self, ratio, p_val, s_val,label=None, s_curve_x_axis_overplot=0.5, s_curve_y_axis_overplot=0.5):

        assert len(ratio) == len(p_val)
        self.df = pd.DataFrame({"ratio": ratio, "p_val": p_val})
        if label is not None:
            self.df["label"] = label
        self.s_curve_y_axis_overplot = s_curve_y_axis_overplot
        self.p_val_cutoff = self.get_p_val_cutoff()
        self.ratio_cutoff = self.get_ratio_cutoff()
        self.df["s_val"] = s_val
        self.ratio_for_s = pd.Series(np.linspace(self.df["ratio"].min() - s_curve_x_axis_overplot, self.df["ratio"].max() + s_curve_x_axis_overplot, num=1000))
        self.p_for_s_larger_1 = self.ratio_for_s.apply(self.calc_p_for_s_equals_1)

    def get_p_val_cutoff(self):

        quant = self.df["ratio"].quantile(0.5)
        return 2.0 + self.df.loc[self.df["ratio"] < quant, "p_val"].median()

    def get_ratio_cutoff(self):

        quant = self.df["p_val"].quantile(0.5)
        median_ = self.df.loc[self.df["p_val"] < quant, "ratio"].median()
        ratio_cutoff_high = np.log2(1.5)
        ratio_cutoff_low = np.log2(0.67)
        return ratio_cutoff_low, ratio_cutoff_high

    def calc_s_from_row(self, row):
        p_val = row["p_val"]
        ratio = row["ratio"]
        return self.calc_s(p_val, ratio)

    def calc_s(self, p_val, ratio):

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

        ratio_cutoff_low, ratio_cutoff_high = self.ratio_cutoff
        ratio_delta_high = ratio - ratio_cutoff_high
        ratio_delta_low = ratio - ratio_cutoff_low

        if ratio > ratio_cutoff_high:
            return (1.0 / ratio_delta_high) + self.p_val_cutoff
        elif ratio < ratio_cutoff_low:
            return (1.0 / (ratio_delta_low * -1)) + self.p_val_cutoff
        else:
            return np.nan

    def get_fig(self, MolColorDF, title="Volcano plot", s_value_cutoff=1000,ColorSwich =2,Labelswitch='labeltable') :#色を変えられる、labeltable'で閾値超えたTable吐く, plotlabel:add mol label on fig
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
        NewColorList=[]    
        if ColorSwich ==1:
            import re
            for label in self.df["label"]:
                r = re.compile("(.*)(_)(.*)"); d = r.search(label); NewColorList += [MolColorDF['MolColor'][d.group(3)] ]           
        elif ColorSwich ==2:
            cond = ( (self.df["ratio"] <= np.log2(0.67) )|( self.df["ratio"] >= np.log2(1.5)) ) & (self.df["s_val"] > -1*np.log10(0.1) ) 

            for ratio in range(len(self.df["ratio"])):
                if cond[ratio] ==1 :      
                    if self.df["ratio"][ratio]>0:
                        NewColorList +=['sienna']
                    else:
                        NewColorList +=['indigo']
                else:
                     NewColorList +=['black']
        else:
            NewColorList='black'

        xmin=np.min(x3); xmax=np.max(x3); ymin=np.min(y3); ymax = np.max(y3)

        ax1.set_ylim([0,ymax+1])
        ax1.plot([np.log2(0.67), np.log2(0.67)],[0, ymax+1], "black", linestyle='dashed' ,linewidth=1.5)
        ax1.plot([np.log2(1.5), np.log2(1.5)],[0, ymax+1], "black", linestyle='dashed' ,linewidth=1.5)
        ax1.plot([xmin,xmax],[1,1], "black", linestyle='dashed' ,linewidth=1.5)

        ax1.xaxis.set_tick_params(labelsize=20)
        ax1.yaxis.set_tick_params(labelsize=20) 
        plt.rcParams['axes.linewidth'] = 2.0

        DF=pd.DataFrame(data=None,columns=['Up','Down'])
        if 'plotlabel' in Labelswitch :  
            cond = self.df["s_val"] > s_value_cutoff
            cond = ( (self.df["ratio"] <= np.log2(0.67) )|( self.df["ratio"] >= np.log2(1.5)) ) & (self.df["s_val"] > -1*np.log10(0.1) ) 
            
            LabelCountDict={}
            for index_, row in self.df[cond].iterrows():
                
                label = row["label"]               
                x_coord = row["ratio"]
                y_coord = row["s_val"]
                ax1.annotate(label, xy=(x_coord, y_coord), xycoords='data', xytext=(5, 5),
                    textcoords='offset points', arrowprops=dict(arrowstyle="-"), fontsize=3)
                LabelCountDict.update({label:x_coord}) 
                
        if Labelswitch == 'labeltable':
            cond = self.df["s_val"] > s_value_cutoff
            cond = ( (self.df["ratio"] <= np.log2(0.67) )|( self.df["ratio"] >= np.log2(1.5)) ) & (self.df["s_val"] > -1*np.log10(0.1) )
            
            LabelCountDict={}
            for index_, row in self.df[cond].iterrows():
                
                label = row["label"]               
                x_coord = row["ratio"]
                y_coord = row["s_val"]

                LabelCountDict.update({label:x_coord}) 
                DF = self.MolUpperLowerCount(LabelCountDict)

        ax1.scatter(x3, y3, color=NewColorList)
        ax1.set_xlabel('log$\it{_{2}}$(Fold Change)', fontsize=40)
        ax1.set_ylabel('-log$\it{_{10}}$(q-value)', fontsize=40)
        plt.rcParams['font.family'] = 'Arial'
        ax1.tick_params(axis='both',labelsize=30)

        return fig,DF
    
    def MolUpperLowerCount(self,Dict):
        import re
        r = re.compile("(.*)(_)(.*)");

        TimetoMolLabel = list(Dict.keys() )
        LabelList = list( set([ r.search(i).group(3) for i in TimetoMolLabel ]) ) 
        NewDF = pd.DataFrame(data=None,index=LabelList,columns=['Up','Down'])

        for j in LabelList:
            posList=[]
            negList =[]
            for i in TimetoMolLabel:
            
                if (r.search(i).group(3) == j) and Dict[i] > 0:
                    posList.append(Dict[i])
                elif (r.search(i).group(3) == j) and Dict[i] < 0:
                    negList.append(Dict[i])
            NewDF.loc[j,'Up'] = len(posList)
            NewDF.loc[j,'Down'] = len(negList)
        return(NewDF)
        
    def get_DistrFC(self,title='Distribution of Fold Change', bins=100):
        fig, ax1 = plt.subplots(figsize=(12, 12))
        ax1.set_title(title, fontsize=22, fontweight="bold")
        
        ratio = self.df["ratio"]

        xmin=np.min(ratio); xmax=np.max(ratio); 

        ax1.hist(ratio.values[~np.isnan(ratio.values)], range=(xmin, xmax), bins = bins, color='#1f77b4',
                              ec='black');plt.title(title); 
        ymin, ymax= ax1.get_ylim()

        ax1.plot([np.log2(0.67), np.log2(0.67)],[0, ymax+20], "black", linestyle='dashed' )   
        ax1.plot([np.log2(1.5), np.log2(1.5)],[0, ymax+20], "black", linestyle='dashed' )
        ax1.set_ylim([0,270.95])
        ax1.set_xlabel('log$\it{_{2}}$(FoldChange)', fontsize=40)
        ax1.set_ylabel('Frequency', fontsize=40)

        ax1.tick_params(axis='both',labelsize=30)
        ax1.xaxis.set_tick_params(labelsize=30)
        return(fig)