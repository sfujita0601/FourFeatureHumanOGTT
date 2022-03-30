#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 02:38:12 2021

@author: fujita
"""
import numpy as np
import scipy.stats as sts
import pandas as pd
from numpy import mean
from numpy import var
from math import sqrt
import six
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols
import itertools
from scipy import linalg
import scipy as sp
from sklearn.covariance import GraphicalLasso
import seaborn as sns
import pydot
import pingouin as pg
def UseR(Pvalue,SwitchDict):
    import pyper
    import pandas as pd
    
    # Python で CSV のデータを読み出す
    #wine = pd.read_csv("wine.csv")
    
    # R のインスタンスを作る
    r = pyper.R(use_numpy = 'True', use_pandas='True')
    
    # Python のオブジェクトを R に渡す
    r.assign("data", Pvalue)
    data=r.get("data")
    r("Pvalue <- data")

    r('library(qvalue)')
    r('library(xlsx)')
    r('Pvalue <- as.matrix((data))')
    r('Pvalue<-as.numeric(Pvalue)')
        
    r('  q.values <- p.adjust(Pvalue, method = "BH")')
    result = r.get("q.values")
    if len(Pvalue.index) < 2:
        r('  QValue <- matrix(q.values,ncol=length(names(data)), nrow=1)')
    #そうでない時
    else:
        r('  QValue <- matrix(q.values,ncol=length(colnames(data)), nrow=length(rownames(data)))')
        r('  colnames(QValue) <- names(data)')
        r('rownames(QValue) <- rownames(data)')
    QvalueBH = r.get("QValue")
    r('qobj <- qvalue(p = Pvalue, fdr.level=0.1)')
    if len(Pvalue.index) < 2:
        r('  colnames(qobj$qvalues) <- names(data)')
        r('rownames(qobj$qvalues) <- rownames(data)')
    else:
        r('  QValue <- matrix(qobj$qvalues,ncol=length(colnames(data)), nrow=length(rownames(data)))')

        r('  colnames(QValue) <- colnames(data)')
        r('rownames(QValue) <- rownames(data)')
    QvalueStorey = r.get("QValue")


    
    return(QvalueStorey,QvalueBH)


def calcpeasonr(list1,list2):#二つのリスト間の重複しているところ計算する
    list1str = delstrinlist(list1)#リスト内にある文字列str要素を削除する、
    list2str = delstrinlist(list2)#リスト内にある文字列str要素を削除する、
    list12str = list(set(list1str+list2str))

    dellist = lambda items, indexes: [item for index, item in enumerate(items) if index not in indexes]

    list1=dellist(list1,list12str)
    list2=dellist(list2,list12str)


    if np.isnan(list2).any() or np.isnan(list1).any():#nanがあったら、両方のリストに存在するとこだけで相関をとる
           #pdb.set_trace() 
           list1nan = list(np.where(np.isnan(list1))[0])
           list2nan = list(np.where(np.isnan(list2))[0])
           list12nan = list(set(list1nan+list2nan))
           list1temp=dellist(list1,list12nan)
           list2temp = dellist(list2,list12nan)
           try:#できたら
               r, p = sts.pearsonr(list1temp,list2temp)
           except:#どっちも空とかなら
               r, p = [np.nan,np.nan]


               
    else:    
        if len(list1)!=len(list2):
            print('1')
        r, p = sts.pearsonr(list1,list2)
    return(r,p)

def delstrinlist(list1):#リスト内にある文字列str要素を削除する、

     dellist=[ i for i in range(len(list1)) if isinstance(list1[i],type('str'))]
     dellist.reverse()
     #for j in dellist:
      #   list1.pop(j)
     return(dellist)