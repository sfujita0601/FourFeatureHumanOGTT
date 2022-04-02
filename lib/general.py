# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:50:52 2017

@author: owner
"""
def getToday():
    import datetime
    now = datetime.datetime.now()
    Today = '{0:%Y%m%d}'.format(now)
    return Today

def ConvertParamLabelRtoPython(Label):#Rで文字化けしたものを統一的な表記に直して打ち込む    
    Labela = [Label[i] + '__a' for i in range(len(Label))]; Labelc = [Label[i] + '__c' for i in range(len(Label))]; Labeld = [Label[i] + '__d' for i in range(len(Label))]; Labele = [Label[i] + '__e' for i in range(len(Label))]
    Labelacde = Labela + Labelc + Labeld + Labele
    #DF.columns = Labelacde; DF.index = Labelacde
    return(Labelacde)
    
def delstrinlist(list1):#リスト内にある文字列str要素を削除する、

     dellist=[ i for i in range(len(list1)) if isinstance(list1[i],type('str'))]
     dellist.reverse()
     #for j in dellist:
      #   list1.pop(j)
     return(dellist)
