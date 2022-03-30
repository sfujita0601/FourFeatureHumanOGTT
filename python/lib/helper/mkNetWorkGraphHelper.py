#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 30 16:53:59 2018

@author: fujita
"""

import networkx as nx
import numpy as np
import pandas as pd
import scipy.io
import itertools
import sys
import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import re
import collections
import itertools
import Helper.MolCorrSubjmeanStdHelper as MSMSH
import xml.etree.ElementTree as ET 
# coding=utf8
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Author: Aric Hagberg (hagberg@lanl.gov)
"""Betweenness centrality measures."""
from heapq import heappush, heappop
from itertools import count
import random

import networkx as nx

__all__ = ['betweenness_centrality', 'edge_betweenness_centrality',
           'edge_betweenness']


def betweenness_centrality(G, k=None, normalized=True, weight=None,
                           endpoints=False, seed=None):
    r"""Compute the shortest-path betweenness centrality for nodes.

    Betweenness centrality of a node $v$ is the sum of the
    fraction of all-pairs shortest paths that pass through $v$

    .. math::

       c_B(v) =\sum_{s,t \in V} \frac{\sigma(s, t|v)}{\sigma(s, t)}

    where $V$ is the set of nodes, $\sigma(s, t)$ is the number of
    shortest $(s, t)$-paths,  and $\sigma(s, t|v)$ is the number of
    those paths  passing through some  node $v$ other than $s, t$.
    If $s = t$, $\sigma(s, t) = 1$, and if $v \in {s, t}$,
    $\sigma(s, t|v) = 0$ [2]_.

    Parameters
    ----------
    G : graph
      A NetworkX graph.

    k : int, optional (default=None)
      If k is not None use k node samples to estimate betweenness.
      The value of k <= n where n is the number of nodes in the graph.
      Higher values give better approximation.

    normalized : bool, optional
      If True the betweenness values are normalized by `2/((n-1)(n-2))`
      for graphs, and `1/((n-1)(n-2))` for directed graphs where `n`
      is the number of nodes in G.

    weight : None or string, optional (default=None)
      If None, all edge weights are considered equal.
      Otherwise holds the name of the edge attribute used as weight.

    endpoints : bool, optional
      If True include the endpoints in the shortest path counts.

    Returns
    -------
    nodes : dictionary
       Dictionary of nodes with betweenness centrality as the value.

    See Also
    --------
    edge_betweenness_centrality
    load_centrality

    Notes
    -----
    The algorithm is from Ulrik Brandes [1]_.
    See [4]_ for the original first published version and [2]_ for details on
    algorithms for variations and related metrics.

    For approximate betweenness calculations set k=#samples to use
    k nodes ("pivots") to estimate the betweenness values. For an estimate
    of the number of pivots needed see [3]_.

    For weighted graphs the edge weights must be greater than zero.
    Zero edge weights can produce an infinite number of equal length
    paths between pairs of nodes.

    References
    ----------
    .. [1] Ulrik Brandes:
       A Faster Algorithm for Betweenness Centrality.
       Journal of Mathematical Sociology 25(2):163-177, 2001.
       http://www.inf.uni-konstanz.de/algo/publications/b-fabc-01.pdf
    .. [2] Ulrik Brandes:
       On Variants of Shortest-Path Betweenness
       Centrality and their Generic Computation.
       Social Networks 30(2):136-145, 2008.
       http://www.inf.uni-konstanz.de/algo/publications/b-vspbc-08.pdf
    .. [3] Ulrik Brandes and Christian Pich:
       Centrality Estimation in Large Networks.
       International Journal of Bifurcation and Chaos 17(7):2303-2318, 2007.
       http://www.inf.uni-konstanz.de/algo/publications/bp-celn-06.pdf
    .. [4] Linton C. Freeman:
       A set of measures of centrality based on betweenness.
       Sociometry 40: 35–41, 1977
       http://moreno.ss.uci.edu/23.pdf
    """
    betweenness = dict.fromkeys(G, 0.0)  # b[v]=0 for v in G
    if k is None:
        nodes = G
    else:
        random.seed(seed)
        nodes = random.sample(G.nodes(), k)
    for s in nodes:
        # single source shortest paths
        if weight is None:  # use BFS
            S, P, sigma = _single_source_shortest_path_basic(G, s)
        else:  # use Dijkstra's algorithm
            S, P, sigma = _single_source_dijkstra_path_basic(G, s, weight)
        # accumulation
        if endpoints:
            betweenness = _accumulate_endpoints(betweenness, S, P, sigma, s)
        else:
            betweenness = _accumulate_basic(betweenness, S, P, sigma, s)
    # rescaling
    betweenness = _rescale(betweenness, len(G), normalized=normalized,
                           directed=G.is_directed(), k=k)
    return betweenness


def edge_betweenness_centrality(G, k=None, normalized=True, weight=None,
                                seed=None):
    r"""Compute betweenness centrality for edges.

    Betweenness centrality of an edge $e$ is the sum of the
    fraction of all-pairs shortest paths that pass through $e$

    .. math::

       c_B(e) =\sum_{s,t \in V} \frac{\sigma(s, t|e)}{\sigma(s, t)}

    where $V$ is the set of nodes, $\sigma(s, t)$ is the number of
    shortest $(s, t)$-paths, and $\sigma(s, t|e)$ is the number of
    those paths passing through edge $e$ [2]_.

    Parameters
    ----------
    G : graph
      A NetworkX graph.

    k : int, optional (default=None)
      If k is not None use k node samples to estimate betweenness.
      The value of k <= n where n is the number of nodes in the graph.
      Higher values give better approximation.

    normalized : bool, optional
      If True the betweenness values are normalized by $2/(n(n-1))$
      for graphs, and $1/(n(n-1))$ for directed graphs where $n$
      is the number of nodes in G.

    weight : None or string, optional (default=None)
      If None, all edge weights are considered equal.
      Otherwise holds the name of the edge attribute used as weight.

    Returns
    -------
    edges : dictionary
       Dictionary of edges with betweenness centrality as the value.

    See Also
    --------
    betweenness_centrality
    edge_load

    Notes
    -----
    The algorithm is from Ulrik Brandes [1]_.

    For weighted graphs the edge weights must be greater than zero.
    Zero edge weights can produce an infinite number of equal length
    paths between pairs of nodes.

    References
    ----------
    .. [1]  A Faster Algorithm for Betweenness Centrality. Ulrik Brandes,
       Journal of Mathematical Sociology 25(2):163-177, 2001.
       http://www.inf.uni-konstanz.de/algo/publications/b-fabc-01.pdf
    .. [2] Ulrik Brandes: On Variants of Shortest-Path Betweenness
       Centrality and their Generic Computation.
       Social Networks 30(2):136-145, 2008.
       http://www.inf.uni-konstanz.de/algo/publications/b-vspbc-08.pdf
    """
    betweenness = dict.fromkeys(G, 0.0)  # b[v]=0 for v in G
    # b[e]=0 for e in G.edges()
    betweenness.update(dict.fromkeys(G.edges(), 0.0))
    if k is None:
        nodes = G
    else:
        random.seed(seed)
        nodes = random.sample(G.nodes(), k)
    for s in nodes:
        # single source shortest paths
        if weight is None:  # use BFS
            S, P, sigma = _single_source_shortest_path_basic(G, s)
        else:  # use Dijkstra's algorithm
            S, P, sigma = _single_source_dijkstra_path_basic(G, s, weight)
        # accumulation
        betweenness = _accumulate_edges(betweenness, S, P, sigma, s)
    # rescaling
    for n in G:  # remove nodes to only return edges
        del betweenness[n]
    betweenness = _rescale_e(betweenness, len(G), normalized=normalized,
                             directed=G.is_directed())
    return betweenness

def mkNetWorkGraph(CorrDF,LabelSum,Thresh,ColorSwitchDict,OptionDict,SubjectName,save_dir):
    r = re.compile("(.*)(_)(.*)") 
    import networkx as nx
    import collections
    import itertools
    from statistics import mean, median,variance,stdev
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    ID = 0.75
    path = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/Eng/RCount_'+ str(ID) +'.xlsx'
    
    try:#隣接行列を定義しているならば、
        Adjacency_matrix = OptionDict['Adjacency_matrix']
        deledge_lists,Combedge,edge_lists = mkEdge_Color(Adjacency_matrix, OptionDict,save_dir)

    except:#隣接行列を定義していなければ
        if OptionDict['mkEdge'] == 'Thresh':#閾値を基準にグラフを描画する場合
            deledge_lists,Combedge,edge_lists,Adjacency_matrix = mkAdjacency_matrix_Thresh(CorrDF,Thresh,path,OptionDict,save_dir)

        elif OptionDict['mkEdge'] == 'Subject':#被験者何人以上で線をひくみたいな
            ID = 0.85#以下のパスがあれば被験者保存解析
            path = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/Eng/RCount_'+ str(ID) +'.xlsx'
            deledge_lists,Combedge,edge_lists,Adjacency_matrix = mkAdjacency_matrix_Thresh(CorrDF,Thresh,path,OptionDict,save_dir)

        elif OptionDict['mkEdge'] == 'Comb':#何種類かのファイルを元に書かれた組み合わせの線をひく
            deledge_lists,Combedge,edge_lists,Adjacency_matrix = mkAdjacency_matrix_Comb(CorrDF,Thresh,path,OptionDict,save_dir)

    # 描画の準備
    G = nx.Graph()
    G.add_weighted_edges_from(edge_lists)
    #####################################################中心性解析
    between_centers =  Centrality(G,Adjacency_matrix, save_dir)
    #####################################################中心性解析
    
    plt.figure(figsize=(20,22))  #描画対象に合わせて設定する
    plt.rcParams["font.size"] = 5
    
    ##############################################エッジの太さ
    if OptionDict['Edge'] == 'Subject':
        OptionDict['Thresh'] = Thresh
        G,GEdge,EdgeDense,edge_color, sumnummin = mkEdgeWidth(G,list(SubjectName[0]),OptionDict,CorrDF)
        edge_width = [ d["weight"]*(d["weight"]*2)**2.1 for (u,v,d) in GEdge]
        EdgeDense=0.3
    elif OptionDict['Edge'] == 'Subject_rev':  
        OptionDict['Thresh'] = Thresh
        G,GEdge,EdgeDense,edge_color, sumnummin = mkEdgeWidth(G,list(SubjectName[0]),OptionDict,CorrDF)
        edge_width = [ d["weight"]*(d["weight"]*2)**2.1 for (u,v,d) in GEdge]
        #EdgeDense=0.3        
    elif OptionDict['Edge'] == 'Comb':
        G,GEdge,EdgeDense,edge_color, sumnummin = mkEdgeColor(G,OptionDict)
        edge_width = [ d["weight"]*(d["weight"]*2)**2.1 for (u,v,d) in GEdge]
        EdgeDense=0.3     
    elif OptionDict['Edge']=='CorrCoef':  
        G,GEdge,EdgeDense,edge_color, sumnummin = mkEdgeWidth(G,list(SubjectName),OptionDict,CorrDF)
        edge_width = [ d["weight"] for (u,v,d) in G.edges(data=True)]
    elif OptionDict['Edge']=='Black': 
        edge_width = [ d["weight"]*(d["weight"]*2)**2.1 for (u,v,d) in G.edges(data=True)]
        EdgeDense=3        
        edge_color = [ 'black' for (u,v,d) in G.edges(data=True) ]
    elif OptionDict['mkglobalDF_Time']=='Each':   
        edge_width = [ d["weight"]*(d["weight"]*2)**2.1 for (u,v,d) in G.edges(data=True)]
        G,GEdge,EdgeDense,edge_color, sumnummin = mkEdgeColo_GlobalEachr(G,OptionDict)
    else:
        edge_width = [ d["weight"]*(d["weight"]*2)**100 for (u,v,d) in G.edges(data=True)]
        EdgeDense=0.1

        #edge_color = [ 'black' for (u,v,d) in G.edges(data=True) ]
        edge_color = [ 'red'  if (u,v) in deledge_lists or (v,u) in deledge_lists else 'black' for (u,v,d) in G.edges(data=True) ]
        #[ print('red:' + u +v )  if (u,v) in deledge_lists else print('black:'+ u+v ) for (u,v,d) in G.edges(data=True) ]
    edge_color = [ 'black' for (u,v,d) in G.edges(data=True) ]
    #aa =[ [u,v,d] for (u,v,d) in G.edges(data=True)]
###################################エッジ
####################################エッジの本数   
    Gedge = G.edges(data=True)
    MetaboDict=dict(zip(list(LabelSum.index),list(LabelSum['Metabo'])))
    MetabonumDict=dict(zip(['glucose','Hormone','Lipid','AminoAcid','Ion','Others'],[13,3,12,31,6,18]))
    MolEdgeDF = pd.DataFrame(data=np.array([0]*83),index=list(Adjacency_matrix.columns),columns=['numEdge'])
    for (u,v,d) in Gedge:#ある相関係数以上の組みが（グラフ可される組み）入ってる
        #代謝グループのリスト、数のdict
        MolEdgeDF.loc[u]+=1/MetabonumDict[MetaboDict[v]]
        MolEdgeDF.loc[v]+=1/MetabonumDict[MetaboDict[u]]


######################################### エッジの色
    if OptionDict['EdgeColor']=='Color_posneg': #相関係数の正負で決める
        #edge_width = [ d["weight"]*(d["weight"]*2)**2.1 for (u,v,d) in G.edges(data=True)]
        EdgeDense=0.1
        pos = nx.shell_layout(G)
        Idx= list(pos.keys());Col= list(pos.keys());
        edge_color = mkEachEdgeColor(G)   
############################################################################ ノードの大きさ設定   
    tag_list = []
    for (u,v,d) in G.edges(data=True):
        tag_list += [u] +[ v]
        #print([u,v,d])

    tag_count = collections.Counter(tag_list)#.most_common(53)
    aa=[];bb=[]
    for ii in range(len(tag_count)):
        aa+=[list(tag_count.keys())[ii]]
        bb += [{"count":list(tag_count.values())[ii]}]
    G.add_nodes_from((aa[ij],bb[ij]) for ij in range(len(aa)))
   # G.add_nodes_from([(tag, {"count":count}) for tag,count in tag_count])
    node_size = [ d["count"]*100 for (n,d) in G.nodes(data=True)]
    node_size = [2000 for (n,d) in G.nodes(data=True) ] #いつもは3000
    if OptionDict['Node_size'] == 'Centrality':#ノードの大きさを媒介中心性に応じて大きく
        node_size = [500000**(v+0.55) for v in between_centers.values()]
    print([ [n,d["count"] ]for (n,d) in G.nodes(data=True)])
            #G.remove_edge(u, v)
            
#################################################################### ノードポジション設定
    np.random.seed(seed=10000)    #ノードポジションの再現性を確保するためseedを設定する
    #pos = nx.spring_layout(G,k=0.15)    #ノードのポジションの計算
    pos = nx.shell_layout(G)
    
    PosDict={}
    LabelSum =   pd.read_excel("//Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx",header=0,index_col=0).iloc[:83,:] 
    count=0;countb=0;countc=0
    #糖代謝、脂質代謝、アミノ酸代謝など分ける？
    Glu = list(LabelSum['English'])[0:13]
    Hormon = list(LabelSum['English'])[13:16]
    Lipid = list(LabelSum['English'])[16:28]
    AA = list(LabelSum['English'])[29:60]
    ION = list(LabelSum['English'])[61:66]
    Others = list(LabelSum['English'])[61:]
    PropLabel =['BMI', 'Age', 'Weight', 'Waist', 'Height', 'TwoHGlucose', 'TwoHInsulin', 'FastingGlucose', 'FastingInsulin', 'Gender']        
    if OptionDict['pos'] == 'negative': #座標を陽に与える、'positive'
        for jj in list(pos.keys()):
            r = re.compile("(.*)(_)(.*)") ; d = r.search(jj) ;
            if jj in ['TwoHGlucose', 'TwoHInsulin', 'FastingGlucose', 'FastingInsulin']:
                count+=1
                x = 0.15*count; y = 0.15#np.random.uniform(0.1, 0.2)
                PosDict.update({jj:x})
            elif jj in ['BMI', 'Weight', 'Waist', 'Height']:
                countb+=1
                x = 0.6+0.15*countb; y = 0.15#np.random.uniform(0.1, 0.2)
                PosDict.update({'Height':0.75});PosDict.update({'Weight':1.05});PosDict.update({'BMI':0.9});PosDict.update({'Waist':1.2})
            elif jj == 'Age':
                x = 0.01; y = 0.15#np.random.uniform(0.1, 0.2)
                PosDict.update({jj:x})     
            elif jj == 'Methionine':
                x = 0.2;y = 0.16
            elif jj == 'Citrulline':
                x = 0.4;y = 0.16
            elif jj== 'Ester type Cho':
                x = 0.6;y = 0.28
            elif jj =='LDL cholesterol':
                x = 0.65;y = 0.24       
            elif jj =='Carnitine':
                x = 0.22; y = 0.2
            elif  jj in AA :
                x = np.random.uniform(countb, 0.8);y = np.random.uniform(0.007+countc, 0.15)
                #x = np.random.uniform(0, 0.8);y = np.random.uniform(0.05, 0.15)
                countb += 0.01; countc = 0.01
            elif  jj in Glu:
                x = 0.2 + count*0.1; y = 0.002
                #x = np.random.uniform(0.2, 0.5);y = np.random.uniform(0.005, 0.08)
            elif  jj in Lipid or jj == 'Acetoacetate':
                x = np.random.uniform(0.55, 0.8);y = np.random.uniform(0.18, 0.3)
            elif  jj in Others:
                x = np.random.uniform(0,0.3);y = np.random.uniform(0.18, 0.3)
            else :
                x = np.random.uniform(0, 0.2); y = np.random.uniform(0.02, 0.08) #; y=np.where(yy >0.1 , np.random.uniform(0.2, 0.3), yy)
            count += 1
            pos[jj] =(x,y)  #これで各ノーどに座標を入れれる,PosDict[jj]
            
    elif OptionDict['pos']=='positive': #座標を陽に与える、'positive'  
        #PosDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Position_20180912_2アミノ酸がいい感じに配置される.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        PosDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Position_20180912_2アミノ酸がいい感じに配置されるIncretin用.xlsx',header=0,index_col=0)
        #PosDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/RespIdx/Position_20181023_rev.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        labels={}
        for jj in list(pos.keys()):
            #if re.compile("(.*)(_)(.*)").search(jj).group(3) in list(PosDF.columns) :
             #   ii=re.compile("(.*)(_)(.*)").search(jj).group(3)
              #  pos[jj] = (PosDF[ii][0],PosDF[ii][1]) 
            labels[jj]=jj
                      #for idx, node in enumerate(G.nodes()):
                      #   labels[node] = FAGLabel[idx]
            #else:
            pos[jj] = (PosDF[jj][0],PosDF[jj][1])

    elif OptionDict['pos']=='positive_Prop': #Prop座標を陽に与える、'positive_Prop'  
        #PosDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/RespIdx/Position_latest.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        PosDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/RespIdx/Position_latest_rev.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)

        for jj in list(pos.keys()):
            pos[jj] = (PosDF[jj][0],PosDF[jj][1])
        #print(max(pos.values), )
        #print('aaaa')
    elif OptionDict['pos']=='AA':#AAのみサークルで描画
        print( [ i[0]  for i in pos.values()] )
        #print( list(pos.values()) )
        PosDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Position_20180912_2.xlsx',header=0,encoding = "ISO-8859-1",index_col=0)
        for jj in list(pos.keys()):#AA以外は陽に与える
            if  jj in Glu+Lipid+Others:

                pos[jj] = (PosDF[jj][0],PosDF[jj][1])
            else:#AAは正規化する
                pos[jj] = (PosDF[jj][0]*0.4,PosDF[jj][1]*0.4)
                PosDF.loc['x',jj] = PosDF[jj][0]*0.4
                PosDF.loc['y',jj] = PosDF[jj][1]*0.4
                #pos[jj] = pos[jj] *10
                #PosDF[jj][0] = pos[jj][0]
                #PosDF[jj][1] = pos[jj][1]  
        PosDF.to_excel(save_dir + 'Position_20180912.xlsx')      
        #PosDF.to_excel(save_dir + 'tempDF.xlsx')
    elif OptionDict['pos']=='circle':#サークルで描画
        pos = nx.shell_layout(G)  
    elif OptionDict['pos']=='Left':#斥力で離す
        pos = nx.spring_layout(G,k=0.5)  
    elif OptionDict['pos']=='Global':#PropMetaboMolTimeなら
        Metabo = ['glucose', 'Hormone', 'Lipid', 'AminoAcid', 'Ion', 'Others']
        Metabopos = {'glucose':18, 'Hormone':3, 'Lipid':10, 'AminoAcid':6, 'Ion':3, 'Others':0}
        #Metabopos = {'glucose':18, 'Hormone':15, 'Lipid':10, 'AminoAcid':6, 'Ion':3, 'Others':0}

        Metabocolor = {'glucose':'red', 'Hormone':'magenta', 'Lipid':'green',  'AminoAcid':'blue', 'Ion':'purple', 'Others':'black'}
        #MolMetabo = dict(zip(PropLabel, Propposx))
        
        PropLabel =['BMI', 'Age', 'Weight', 'Waist', 'Height', 'TwoHGlucose', 'TwoHInsulin', 'FastingGlucose', 'FastingInsulin', 'Gender']        
        Propposx = [ -0.015, -0.005, -0.0075, -0.0025, -0.0070, -0.0045, -0.0015, -0.0025, -0.015, 0.02]
        Propposy = [ 5.5, 4.5, 12, 8, 6, 16, 14, 19, 18, 12.8]
        Propposx = dict(zip(PropLabel, Propposx)); Propposy = dict(zip(PropLabel, Propposy))
        
        PropFastAfter = [ 'Glu_Fast',	'Glu_After',	'Hormon_Fast',	'Hormon_After',	'Lipid_Fast',	'Lipid_After',	'AA_Fast',	'AA_After',	'Ion_After',	'Others_Fast',	'Others_After']
        PropFastAfter = [ 'Glu_After',	'Hormon_Fast',	'Hormon_After',		'Lipid_After',	'AA_Fast',	'AA_After',		'Others_After']
        PropFastAfterpos = {'Glu_After':20,	'Hormon_Fast':4,	'Hormon_After':2,		'Lipid_After':10,	'AA_Fast':15,	'AA_After':6,		'Others_After':18}
        PropFastAftercolor = {'Glu_After':'red', 'Hormon_Fast':'magenta', 'Hormon_After':'magenta','Lipid_After':'green',  'AA_Fast':'blue', 'AA_After':'blue','Ion':'purple', 'Others_After':'black'}
        
        PropFast2H4H = ['Glu_Fast', 'Glu_2H', 'Glu_4H', 'Hormone_Fast', 'Hormone_2H', 'Hormone_4H', 'Lipid_Fast',  'Lipid_2H', 'Lipid_4H', 'AA_Fast',  'AA_2H', 'AA_4H', 'Ion_Fast', 'Ion_2H', 'Ion_4H', 'Others_Fast', 'Others_2H','Others_4H']
        PropFast2H4Hpos = {'Glu_2H':20,	'Glu_4H':16, 'Hormone_Fast':4,	'Hormone_2H':2,	'Hormone_4H':8,	'Lipid_4H':10,	'Lipid_2H':12,'AA_Fast':15,	'AA_4H':6,'AA_2H':5,		'Others_2H':18}
        PropFast2H4Hcolor = {'Glu_2H':'red', 'Glu_4H':'red','Hormone_Fast':'magenta', 'Hormone_2H':'magenta','Hormone_4H':'magenta','Lipid_2H':'green', 'Lipid_4H':'green',  'AA_Fast':'blue', 'AA_2H':'blue', 'AA_4H':'blue','Ion':'purple', 'Others_2H':'black'}

        #Propposx = [ -0.01, -0.005, -0.005, -0.0025, -0.0025, -0.0075, -0.0025, -0.0025, -0.01]
        #Propposy = [ 8, 4.5, 12, 8, 6, 16, 14, 20, 20]
        RespIdxDict = {'SM-C IGF-1':0,'Cystine':0.45,'Glutarate':0.9,'Glutamic acid':1.35,'TG(Neutral lipid)':1.8,'Histidine':2.25 }
        RespIdxGainDict = {'Serum iron':12,'Lactate':12.45,'Glucose':21.9,'2-Hydroxypentanoate':21.0,'Insulin':21.45}#,'TG(Neutral lipid)':1.8,'Histidine':2.25 }
       
        MolNeed=[]
        MolLabel = list(LabelSum.index)
        for jj in list(pos.keys()):
            if jj in MolLabel:
                MolNeed.append(jj)
                
        Mol = [MolNeed[i] for i in range(len(MolNeed))]
        Mol.reverse()
        Molpos = [-1+0.6*j for j in range(len(MolNeed))]
        Molpos =dict(zip(Mol,Molpos))
        
        Time = ['-10','0','10','20','30','45','60','75','90','120','150','180','210','240']
        Timepos = [0+2*ii for ii in range(len(Time))]
        Timepos =dict(zip(Time,Timepos))    
        
        PosDF = pd.DataFrame(data=None,index=['x','y','Label'],columns=PropLabel)
        pos = nx.shell_layout(G);poscount=4.5#初期値：-0.2
        yFast=2;yGain=20#12.9
        yAUC=13;FAGLabel=[]
        for jj in list(pos.keys()):
            tempLabel = [ j for j in Mol if jj[-5:] in j ] 
            tempLabel.append('')
            #print( tempLabel)
            if jj in Metabo:
                pos[jj] = (0,Metabopos[jj])
            elif jj in PropLabel:
                pos[jj] = (-0.005,poscount)
                #x = np.random.uniform(-0.01,0);y = np.random.uniform(-0.005, 20)
                x =Propposx[jj];  y = Propposy[jj] 
                pos[jj] = ( x - 0.01, y )
                FAGLabel+=[jj]
            elif jj in Mol:
                pos[jj] = (0.005+ 0.005,Molpos[jj])  
            elif jj in PropFastAfter:
                pos[jj] = (0.005+ 0.005,PropFastAfterpos[jj]) 
            elif jj in PropFast2H4H:
                pos[jj] = (0.005+ 0.005,PropFast2H4Hpos[jj] )
            elif jj in Time:
                pos[jj] = (0.01,Timepos[jj]) 
            elif  tempLabel[0] in Mol:
                tempLabel = [ j for j in Mol if jj[-5:] in j ] 
                pos[jj] = (0.020,Molpos[tempLabel[0]]+np.random.uniform(-0.3, 0.3))                 
                #print(jj)                
            elif OptionDict['pos_Global']=='Fasting_AUC_Gain':
                try:
                    if r.search(jj).group(1)=='Fasting':
                        #if r.search(jj).group(3) in list( RespIdxDict.keys() ):
                         #   yFastrev = RespIdxDict[r.search(jj).group(3)]
                          #  yFast+=0.45
                           # y=yFastrev
                        #else:
                            y=yFast
                            yFast+=0.5
                            x=0.04
                            x = np.random.uniform(0.03,0.07)
                            NodeLabel=r.search(jj).group(3)
                        #x = np.random.uniform(0,0.3);y = np.random.uniform(0,9)
                    elif r.search(jj).group(1)=='Gain':
                        #if r.search(jj).group(3) in list( RespIdxGainDict.keys() ):
                         #   yGainrev = RespIdxGainDict[r.search(jj).group(3)]
                            #yGain+=0.45
                          #  y=yGainrev 
                        #else:
                            y=yGain
                            yGain+=0.45
                            x=0.04
                            x = np.random.uniform(0.02,0.06)
                            
                            NodeLabel=r.search(jj).group(3)
                    #elif r.search(jj).group(1)=='AUC':
                        #if r.search(jj).group(3) in list( RespIdxGainDict.keys() ):
                         #   yGainrev = RespIdxGainDict[r.search(jj).group(3)]
                            #yGain+=0.45
                           # y=yGainrev                           
                    else:
                        yAUC+=0.5
                        x = np.random.uniform(0.03,0.07);y = np.random.uniform(10,20)
                        y=yAUC
                        
                    pos[jj] = (x,y)
                    NodeLabel=r.search(jj).group(3)
                    PosDF.loc['x',jj] = pos[jj][0]
                    PosDF.loc['y',jj] = pos[jj][1]
                    #pos[jj] = pos[jj] *10
                    #PosDF[jj][0] = pos[jj][0]
                    #PosDF[jj][1] = pos[jj][1]  
                except:
                    pos[jj] = (0.040,poscount*1.4) 
                    poscount += 0.2  
                    NodeLabel=jj
                
                FAGLabel+=[NodeLabel]

            else:
                try:
                    if r.search(jj).group(1)=='Fasting':
                        if r.search(jj).group(3) in list( RespIdxDict.keys() ):
                            yFastrev = RespIdxDict[r.search(jj).group(3)]
                            yFast+=0.45
                            y=yFastrev
                        else:
                            y=yFast
                            yFast+=0.45
                        x=0.04
                        #x = np.random.uniform(0,0.3);y = np.random.uniform(0,9)
                    elif r.search(jj).group(1)=='Gain':
                        if r.search(jj).group(3) in list( RespIdxGainDict.keys() ):
                            yGainrev = RespIdxGainDict[r.search(jj).group(3)]
                            #yGain+=0.45
                            y=yGainrev 
                        else:
                            y=yGain
                            yGain+=0.45
                        x=0.04
                    #elif r.search(jj).group(1)=='AUC':
                        #if r.search(jj).group(3) in list( RespIdxGainDict.keys() ):
                         #   yGainrev = RespIdxGainDict[r.search(jj).group(3)]
                            #yGain+=0.45
                           # y=yGainrev 
                          
                    else:
                        yGain-=0.45
                        x = np.random.uniform(0.09,0.1);y = np.random.uniform(10,20)
                        y=yGain
                    pos[jj] = (x,y)
                except:
                    pos[jj] = (0.040,poscount*1.4) 
                    poscount += 0.2                
####################################################################
    #陽に与える
    # 座標を設定する．indexがid，代入している値が座標．
    #pos={}
    #pos["nodeA"]=(0,0)
    #pos["nodeB"]=(1,1)
    #pos["nodeC"]=(3,8)
    
    #print(FAGLabel)
    #PosDF.loc['Label']=FAGLabel
    #PosDF.to_excel(save_dir + 'Position_20181023.xlsx')
    #print( list(PosDF.loc['Label']     )  )
    #FAGLabel =  list(PosDF.loc['Label']     )         
    numNode = nx.number_of_nodes(G) 
    NumClst = nx.number_connected_components(G)
    
    #平均クラスターサイズと分散を求める
    EachClstnode = sorted(nx.connected_components(G), key = len, reverse=True)
    numEachClstnode = [];VarnumEachClstnode=[]
    for i in range(len(EachClstnode)):
        numEachClstnode += [len(EachClstnode[i])]
    AvenumEachClstnode = mean(numEachClstnode) ;
    VarnumEachClstnode = variance(numEachClstnode) if len(numEachClstnode)>1 else 0
#################################################################### ノードの色
    ColorList=[]
    if 'Acetoacetate' in list(G.nodes) or 'GLP-1' in list(G.nodes) or 'β- aminoisobutyric acid' in list(G.nodes) or '1-methyl-histidine' in list(G.nodes) or 'Aspartic acid' in list(G.nodes) or 'M ethanolamine' in list(G.nodes):
        LabelSum.loc['Acetoacetate','MolColor'] = 'green'
        LabelSum.loc['GLP-1','MolColor'] = 'red'
        LabelSum.loc['β- aminoisobutyric acid','MolColor'] = 'blue'
        LabelSum.loc['1-methyl-histidine','MolColor'] = 'blue'
        LabelSum.loc['Aspartic acid','MolColor'] = 'blue'
        LabelSum.loc['M ethanolamine','MolColor'] = 'blue'       
    if ColorSwitchDict  == 'ClstColor':#クラスタで色分け
        nodeList = list(G.nodes);ColorList = [LabelSum.loc[NodeList[x],'ClstColor'] for x in range(len(NodeList))]
    elif ColorSwitchDict  == 'MolColor':#分子で色分け
        NodeList = list(G.nodes);
        try:
            if re.compile("(.*)(_)(.*)").search(NodeList[0]).group(1) in list(LabelSum.index) :
                ii=re.compile("(.*)(_)(.*)").search(NodeList[0]).group(1)
                ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1),'MolColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1) in list(LabelSum.index) else 'white' for x in range(len(NodeList))  ]
        except:
            pass
        try:
            if OptionDict['Check'] == 'Amino':#AAのラベルなど変える
                MolColor,a = LH.AACheck(LabelSum,LabelSum,OptionDict['AminoCheck'])#'protein','ketogenic','EAA','SemiEAA'
            elif OptionDict['Check'] == 'Glc':#糖代謝系のラベルなど変える
                MolColor,a = LH.TCACheck(MolColor,a,'TCA')#'TCA'
            elif OptionDict['Check'] == 'AminoGlc':#糖代謝系、AAのラベルなど変える
                MolColor,a = LH.AACheck(LabelSum,a,OptionDict['AminoCheck'])
                LabelSum,a = LH.TCACheck(MolColor,a,'TCA')#'TCA'
                #LabelSum['MolColor'] = MolColor
        except:
            pass
        ColorList = [LabelSum.loc[NodeList[x],'MolColor'] for x in range(len(NodeList))]
        #except:
         #   pass
    elif ColorSwitchDict == 'TimeVarColor':#時間ばらつきで色分け
        NodeList = list(G.nodes);ColorList = [LabelSum.loc[NodeList[x],'TimeVarColor'] for x in range(len(NodeList))]
    elif ColorSwitchDict == 'Global':#PropMetaboMolTimeなら
        print('TYRYRRV')
        NodeList = list(G.nodes)
        for i in NodeList:
            if i in Metabo:
                ColorList.append(Metabocolor[i])
            elif i in PropFastAfter:
                ColorList.append(PropFastAftercolor[i])
            elif i in PropFast2H4H:
                ColorList.append(PropFast2H4Hcolor[i])
            else:
                ColorList.append('white')
    elif ColorSwitchDict == 'RespIdx':#応答特徴量なら
        NodeList = list(G.nodes);
        #print( r.search(NodeList[0]))
        for i in NodeList:
            if i in PropLabel:
                    ColorList.append('white')            
            else:         
                    ColorList.append(LabelSum.loc[ r.search(i).group(3),'MolColor'] )#if  r.search(NodeList[x]).group(3) in list(LabelSum.index)  else 'white'  for x in range(len(NodeList)) ]      
    else:
        ColorList=['white']*len(G.nodes)    
    nx.draw_networkx_nodes(G, pos, node_size=node_size,node_color=ColorList,alpha=0.5)  
    #print(ColorList)

    #print_strongly_connected_components(edge_lists)
    
    if OptionDict['pos_Global'] == 'LabelLoop':#ラベルをループして描く        
        nx.draw_networkx_labels(G,pos=pos,font_size=12,alpha=1,font_weight="bold")#node_color=ColorList,
    elif OptionDict['pos_Global']=='Fasting_AUC_Gain':
        labels={}
        #print(FAGLabel)
        for idx, node in enumerate(G.nodes()):
            labels[node] = FAGLabel[idx]
        #print(labels)
        nx.draw_networkx_labels(G,pos=pos,labels=labels,font_size=15,alpha=1,font_weight="bold")#node_color=ColorList,
    else:
        #pass #ノードのラベル消す時
        nx.draw_networkx_labels(G,pos=pos,labels=labels,font_size=16,alpha=1,font_weight="bold")#node_color=ColorList,        
    nx.draw_networkx_edges(G, pos, alpha=0.5,edge_color=edge_color, width=edge_width)
    plt.axis("off")
    numsub = 20 - len(set(edge_color))
    #plt.title('相関係数平均>'+('%02.2f' %Thresh)+'_島の数：'+str(NumClst) + '   相関係数>'+('%02.2f' %Thresh)+ 'となる最小被験者数：'+str(sumnummin),size='30')
    plt.title('相関係数平均>'+('%02.2f' %Thresh)+'_島の数：'+str(NumClst),size='30')   
    plt.savefig(save_dir+'NewtWork_' +('%02.2f' %Thresh)+'_'+ str(NumClst) + ColorSwitchDict+'.pdf',bbox_inches="tight")
    plt.savefig(save_dir+'NewtWork_' +('%02.2f' %Thresh)+'_'+ str(NumClst) + ColorSwitchDict+'.png',bbox_inches="tight")
##### community
    """
    import community
    #import networkx as nx
    #G = nx.karate_club_graph()
    partition = community.best_partition(G)
    #print(ColorList)
    #labels = [i for i in range(nx.number_of_nodes(G))]#dict([(i, str(i)) for i in range(nx.number_of_nodes(G))])
    ColorDict = dict(zip([aa[i] for i in range(len(aa))],[ColorList[i] for i in range(len(aa))]))
    nx.set_node_attributes(G,  ColorDict,'Color')

    nx.set_node_attributes(G,  labels,'label')
    nx.set_node_attributes(G,  partition,'community')

    nx.write_gml(G, save_dir+"community_4.gml")
    #plt.show()
    """
    Gdeg = dict(nx.degree(G))

    #平均次数  
    #次数の分布
    plt.figure()
    bins = range(1,12)
    plt.rcParams["font.size"] = 15
    plt.hist(Gdeg.values(), bins=15)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.ylabel('Frequency')
    plt.xlabel('Degree')    
    plt.tick_params(labelsize=10) 
    plt.xlim([0,70])
    plt.ylim([0,20])
    plt.savefig(save_dir+'相関係数>'+('%02.2f' %Thresh)+'DistributionOFK.pdf', dpi=300)
    #plt.show()    


    #ブリッジ係数を算出する。
    ##BCDict = GraphHel.calcBrdgingCoefficient(Adjacency_matrix,labels,Gdeg,between_centers,save_dir)
    
    # 次数中心性の分布を見る
    plt.figure()
    plt.hist(nx.degree_centrality(G).values(), bins=20)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel('次数中心性')
    plt.ylabel('分子名数')
    plt.title('分子次数中心性の分布')
    plt.tight_layout()
    plt.savefig(save_dir+'dc.pdf', dpi=300)
    #plt.show()
    
    MolEdgeDF.to_excel(save_dir+'DegNormalized.xlsx')
   
    # 中心性の高いノードを抽出し描く
    """
    central_nodes = [k for k, v in nx.degree_centrality(G).items() if v >= 0.4]
    G_central = G.subgraph(central_nodes)
    if ColorSwitch  == 1:#クラスタで色分け
        NodeList = list(G_central.node);ColorList = [LabelSum.loc[NodeList[x],'ClstColor'] for x in range(len(NodeList))]
    elif ColorSwitch  == 2:#分子で色分け
        NodeList = list(G_central.node);ColorList = [LabelSum.loc[NodeList[x],'MolColor'] for x in range(len(NodeList))]
    elif ColorSwitch == 3:#時間ばらつきで色分け
        NodeList = list(G.node);ColorList = [LabelSum.loc[NodeList[x],'TimeVarColor'] for x in range(len(NodeList))]
    
    #draw_char_graph(G_central, save_dir+'central.pdf',ColorList)
    """
    return([NumClst,AvenumEachClstnode,VarnumEachClstnode, numNode,MolEdgeDF])

def metaboGroup(CorrDF,LabelSum):#代謝グルーぷ同じなら1,違うなら0の行列を作成
    Col = list(CorrDF.columns)
    ColorList = list(LabelSum.loc[Col,'MolColor']    )
    
    MolColorDF = pd.DataFrame(data =np.ones([83,83]),index=ColorList,columns=ColorList)
    ColSet = list(set(MolColorDF.columns))
    for i in range(len(ColSet)):
        MolColorDF.loc[ColSet[i],ColSet[i]]=0
    MolColorDF.index= Col;MolColorDF.columns= Col
        
    return(CorrDF[MolColorDF==1.0])

def calcBrdgingCoefficient(Adjacency_matrix,label,Gdeg,between_centers,save_dir):
    Dv=0
    Di=0
    BCDict=Gdeg.copy()
    CRDict = Gdeg.copy()
    
    for i in label:
        AdLabel = list(Adjacency_matrix.loc[i][Adjacency_matrix.loc[i]==1].index) #分子iが隣接している分子名
        Dv = 1/Gdeg[i]
        for j in  AdLabel:
            Di += 1/Gdeg[j]
        BCDict[i] = Dv / Di
        CRDict[i] = BCDict[i]* between_centers[i]
    BCDict_sorted = sorted(BCDict.items(), key=lambda x: x[1], reverse=True)
    CRDict_sorted = sorted(CRDict.items(), key=lambda x: x[1], reverse=True)

    drawBarCent(BCDict_sorted,'BridgingCoefficient', save_dir)#媒介中心性ソート棒グラフ
    drawBarCent(CRDict_sorted,'BridgingCentrality', save_dir)#媒介中心性ソート棒グラフ


    return(BCDict)
            
            
def Centrality(G, Adjacency_matrix,save_dir):#中心性解析
    between_centers = nx.betweenness_centrality(G)
    #print('媒介中心性')
    #print(sorted(between_centers.items(), key=lambda x: x[1], reverse=True))
### 0のkey削除
    between_centers = {k: v for k, v in between_centers.items() if v!= 0}

    between_centers_sorted = sorted(between_centers.items(), key=lambda x: x[1], reverse=True)
    drawBarCent(between_centers_sorted,'Betweenness Centrality', save_dir)#媒介中心性ソート棒グラフ
    mkDF(between_centers_sorted,'Betweenness Centrality', save_dir)#媒介中心性のエクセルファイル出力
    
    calcnumedge(Adjacency_matrix,between_centers_sorted,save_dir)#エッジ数ソート棒グラフ
    
    close_centers = nx.closeness_centrality(G)
    #print('近接中心性')
    #print(sorted(close_centers.items(), key=lambda x: x[1], reverse=True)[:5])
    close_centers_sorted = sorted(close_centers.items(), key=lambda x: x[1], reverse=True)
    drawBarCent(close_centers_sorted,'Closeness Centrality', save_dir)#媒介中心性ソート棒グラフ
    

    degree_centers = nx.degree_centrality(G)
    #print('次数中心性')
    #print(sorted(degree_centers.items(), key=lambda x: x[1], reverse=True)[:5])
    degree_centers_sorted = sorted(degree_centers.items(), key=lambda x: x[1], reverse=True)
    drawBarCent(degree_centers_sorted,'Degree Centrality', save_dir)#媒介中心性ソート棒グラフ

    ##eigen_centers = nx.eigenvector_centrality_numpy(G)
    #print('固有ベクトル中心性')
    #print(sorted(eigen_centers.items(), key=lambda x: x[1], reverse=True)[:5])
    ##eigen_centers_sorted = sorted(eigen_centers.items(), key=lambda x: x[1], reverse=True)
    ##drawBarCent(eigen_centers_sorted,'Eigen Centrality', save_dir)#媒介中心性ソート棒グラフ
    
    
    #print('クラスター係数')
    #print(nx.average_clustering(G))
    #average_clustering = nx.average_clustering(G)
    #average_clustering_sorted = sorted(average_clustering.items(), key=lambda x: x[1], reverse=True)
    #drawBarCent(average_clustering_sorted,'Average Clustering', save_dir)#媒介中心性ソート棒グラフ

    #print('平均距離')
    #print(nx.average_shortest_path_length(G))
    
    #G_random = nx.connected_watts_strogatz_graph(30, 2, 1)
    #print('ランダム距離と比較して')
    #print(nx.average_shortest_path_length(G_random))
    #G_random_sorted = sorted(G_random.items(), key=lambda x: x[1], reverse=True)
    #drawBarCent(G_random_sorted,'vs Random distance', save_dir)#媒介中心性ソート棒グラフ
    return(between_centers)    

def mkDF(between_centers,file_name,save_dir):#媒介中心性などをエクセルで出力
    valueList = [  between_centers[i][1]  for i in range(len(between_centers) ) ]
    nameList =     [ between_centers[i][0] for i in range(len(between_centers) ) ]

    NewDF = pd.DataFrame(data=valueList,columns=[file_name],index=nameList)
    NewDF.T.to_excel(save_dir+ file_name + '_DF.xlsx')
    
# obsolete name
def drawBarCent(between_centers,file_name,save_dir):#媒介中心性などの大きさをを棒グラフで表す
    #タプルでくるので値を取り出す
    valueList = [  between_centers[i][1]  for i in range(len(between_centers) ) ]
    nameList =     [ between_centers[i][0] for i in range(len(between_centers) ) ]
    x = np.linspace(0, len(between_centers),len(between_centers) )
    fig = plt.figure()
    plt.bar( x, valueList  )
    #plt.yscale("symlog")
    plt.xticks(x,nameList,rotation=270 ,size=8)
    plt.gca().yaxis.set_tick_params(labelsize=20)
    
    plt.title(file_name)
    plt.rcParams['axes.linewidth'] = 1.5 #軸の太さを設定。目盛りは変わらない
    
    plt.savefig(save_dir+ file_name + '.pdf',bbox_inches="tight")
    plt.savefig(save_dir+ file_name + '.png',bbox_inches="tight")


    
    plt.close()
    
def calcnumedge(DF,between_centers,save_dir):#隣接行列から1の数をカウントし、ノーどごとにエッジの数を棒グラフ化
    Col = [ between_centers[i][0] for i in range(len(between_centers) ) ]

    NameList = Col
    ValueList = []

    for j in Col:
        ValueList.append( len( DF[j][ DF[j] == 1 ] ) )
    CombDict = dict( zip( NameList, ValueList ) )      
    
    
    CombDict_sorted = sorted(CombDict.items(), key=lambda x: x[1], reverse=True)
    ValueList = [ CombDict_sorted[i][1] for i in range(len(CombDict_sorted) ) ]
    NameList= [ CombDict_sorted[i][0] for i in range(len(CombDict_sorted) ) ]
    x = np.linspace(0, len(ValueList),len(ValueList) )
    fig = plt.figure()
    plt.bar( x,ValueList )
    plt.xticks(x,NameList,rotation=270,size=5 )
    
    
    plt.title('エッジ数')
    plt.savefig(save_dir+'BarnumEdge.pdf',bbox_inches="tight")  
    plt.close()
    
def edge_betweenness(G, k=None, normalized=True, weight=None, seed=None):
    return edge_betweenness_centrality(G, k, normalized, weight, seed)


# helpers for betweenness centrality

def _single_source_shortest_path_basic(G, s):
    S = []
    P = {}
    for v in G:
        P[v] = []
    sigma = dict.fromkeys(G, 0.0)    # sigma[v]=0 for v in G
    D = {}
    sigma[s] = 1.0
    D[s] = 0
    Q = [s]
    while Q:   # use BFS to find shortest paths
        v = Q.pop(0)
        S.append(v)
        Dv = D[v]
        sigmav = sigma[v]
        for w in G[v]:
            if w not in D:
                Q.append(w)
                D[w] = Dv + 1
            if D[w] == Dv + 1:   # this is a shortest path, count paths
                sigma[w] += sigmav
                P[w].append(v)  # predecessors
    return S, P, sigma


def _single_source_dijkstra_path_basic(G, s, weight):
    # modified from Eppstein
    S = []
    P = {}
    for v in G:
        P[v] = []
    sigma = dict.fromkeys(G, 0.0)    # sigma[v]=0 for v in G
    D = {}
    sigma[s] = 1.0
    push = heappush
    pop = heappop
    seen = {s: 0}
    c = count()
    Q = []   # use Q as heap with (distance,node id) tuples
    push(Q, (0, next(c), s, s))
    while Q:
        (dist, _, pred, v) = pop(Q)
        if v in D:
            continue  # already searched this node.
        sigma[v] += sigma[pred]  # count paths
        S.append(v)
        D[v] = dist
        for w, edgedata in G[v].items():
            vw_dist = dist + edgedata.get(weight, 1)
            if w not in D and (w not in seen or vw_dist < seen[w]):
                seen[w] = vw_dist
                push(Q, (vw_dist, next(c), v, w))
                sigma[w] = 0.0
                P[w] = [v]
            elif vw_dist == seen[w]:  # handle equal paths
                sigma[w] += sigma[v]
                P[w].append(v)
    return S, P, sigma


def _accumulate_basic(betweenness, S, P, sigma, s):#S：クラスター含有ノード名前、P：Sクラスターのs以外のノードが繋がっているノード、s：対象のノー度
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()#Sのケツをとる
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            betweenness[w] += delta[w]# wの媒介性を求める
    return betweenness


def _accumulate_endpoints(betweenness, S, P, sigma, s):
    betweenness[s] += len(S) - 1
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            betweenness[w] += delta[w] + 1
    return betweenness


def _accumulate_edges(betweenness, S, P, sigma, s):
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            c = sigma[v] * coeff
            if (v, w) not in betweenness:
                betweenness[(w, v)] += c
            else:
                betweenness[(v, w)] += c
            delta[v] += c
        if w != s:
            betweenness[w] += delta[w]
    return betweenness


def _rescale(betweenness, n, normalized, directed=False, k=None):
    if normalized:
        if n <= 2:
            scale = None  # no normalization b=0 for all nodes
        else:
            scale = 1.0 / ((n - 1) * (n - 2))
    else:  # rescale by 2 for undirected graphs
        if not directed:
            scale = 0.5
        else:
            scale = None
    if scale is not None:
        if k is not None:
            scale = scale * n / k
        for v in betweenness:
            betweenness[v] *= scale
    return betweenness


def _rescale_e(betweenness, n, normalized, directed=False, k=None):
    if normalized:
        if n <= 1:
            scale = None  # no normalization b=0 for all nodes
        else:
            scale = 1.0 / (n * (n - 1))
    else:  # rescale by 2 for undirected graphs
        if not directed:
            scale = 0.5
        else:
            scale = None
    if scale is not None:
        if k is not None:
            scale = scale * n / k
        for v in betweenness:
            betweenness[v] *= scale
    return betweenness


def mkEdgeColo_GlobalEachr(G,OptionDict):
    
    Gedge = G.edges(data=True)
    EdgeDense=[]
    
    #CountDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/被験者ごと_平均/Adjacency_matrix_BH.xlsx',header=0,encoding = "ISO-8859-1")
    CountDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180920/ColorDF.xlsx',header=0,encoding = "ISO-8859-1")
    edge_colorList=[]
    for (u,v,d) in Gedge:#全ノードの組み合わせが入ってる
        #try:
         #   EdgeWidth =+int(CountDF.loc[u,v])
        #except:
         #   print(u ,v )
        #print(int(CountDF.loc[u,v]))
        edge_color =  CountDF.loc[u,v]#ClstColorList[( 20 - int(CountDF.loc[u,v]) )]#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
        edge_colorList.append(edge_color) 
        sumnummin = []
        EdgeWidth=0.3
        #DrawHistnumEdge(EdgenumSubj)
        #print(edge_color)
        
    return(Gedge,EdgeDense,edge_colorList, sumnummin)
    
def mkEdge_Color(Adjacency_matrix, OptionDict,save_dir):#隣接行列を与えられ、張るエッジと色を決めておく
    # networkxに引き渡すエッジリストを準備
    Idx = Adjacency_matrix
    edge_lists = [];tempedge_list =[]
    deledge_lists = []#次消えるエッジ
    Combedge=[]; Combtemedge =[]     
    for i in Idx:#fileだけ回す
        for j in Idx:#行の数だけ回す
            if Adjacency_matrix.loc[i,j] == 1:#隣接行列で1なら
                tmp = (i,j,0.5) #(from,to,weight)のタプルを作成
                edge_lists.append(tmp)  
                #Adjacency_matrix.loc[ref_file[i].iloc[ii][0],ref_file[i].iloc[ii][1] ] = color
                #Adjacency_matrix.loc[ref_file[i].iloc[ii][1],ref_file[i].iloc[ii][0] ] = color            
                        #組み合わせだけとをり出す
            #Combedge = [edge_lists[j][0:2] for j in range(len(edge_lists))]; Combtemedge = [tempedge_list[j][0:2] for j in range(len(tempedge_list))]
    Adjacency_matrix.to_excel(save_dir+'Color_matrix.xlsx')
    
    return(deledge_lists,Combedge,edge_lists)
    
def setMolTimeColor(Idx,Col,CorrDF):#関数で色を決める     
    #CorrDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/CombDF_Eng.xlsx',header=0,encoding = "ISO-8859-1")

    try:
        try:
            if len(CorrDF.loc[Idx]) == 5:#1つしかない
                if CorrDF.loc[Idx,'Corr'] > 0:
    
                    color = 'red'    
                else:
                    color = 'blue'
            else:
                if any(CorrDF.loc[Idx,'Corr']) > 0 or CorrDF.loc[Idx,'Corr'] > 0:
                    color = 'red'
                else:
                    color = 'blue'     
        except:
            if len(CorrDF.loc[Col]) == 5:#1つしかない
                if  CorrDF.loc[Col,'Corr'] > 0:
                    color = 'red'
                else:
                    color = 'blue'
            else:
                if any(CorrDF.loc[Col,'Corr']) > 0 or CorrDF.loc[Col,'Corr'] > 0:
                    color = 'red'
                else:
                    color = 'blue'       
    except:
        print(Idx, Col)
        color='black'
    return(color)
    
def checkcolor(Idx,Col):#隣接行列から色行列を作成するときに色を割り当てる
    Metabo = ['glucose', 'Hormone', 'Lipid', 'AminoAcid', 'Ion', 'Others']
    Prop =['BMI', 'Age', 'Weight', 'Waist', 'Height', 'TwoHGlucose', 'TwoHInsulin', 'FastingGlucose', 'FastingInsulin']        

    if Idx in Prop or Idx in Metabo:#プロパティは黒
        color = 'black'
    elif Col in Prop  or Col in Metabo:#プロパティは黒
        color = 'black'
    else:#その他は MolTimeを含む
        color = setMolTimeColor(Idx,Col)#関数で色を決める        
    return(color)
    
def setColorResp(Idx,Col,CombDF):#関数で色を決める     
    #CorrDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/CombDF_Eng.xlsx',header=0,encoding = "ISO-8859-1")
    try:
        try:

            if len(CombDF.loc[Idx].columns)>0:#1つだけだったらexceptへ
                if list(CombDF.loc[Idx][ CombDF.loc[Idx]['Mol2']==Col]['Corr'])[0] > 0:
    
                    color = 'red'    
                else:
                    color = 'blue'
    
        except:#Idxが一つしかない or Idxがindexに入ってない
            try:#Idxが一つしかない
                if  list(CombDF.loc[Idx]['Corr'])[0]> 0:
                    color = 'red'
                else:
                    color = 'blue'
            except:#Idxがindexに入ってない
                try:#IdxがMol2にある
                    if  len(CombDF[ CombDF['Mol2']==Idx]['Corr'] )>0:#1つだけだったらexceptへ
                        if list(CombDF[ (CombDF['Mol2']==Idx) & (CombDF['Mol1']==Col)]['Corr'])[0] > 0:            
                            color = 'red'    
                        else:
                            color = 'blue'  
                            print('ここまできた')
                except:
                    try:
                        if list(CombDF[ CombDF['Mol2']==Idx]['Corr'])[0] > 0:
    
                            color = 'red'
                        else:
                            color = 'blue' 
                        print('ここまできたのか？')
                    except:
                        print('何がいけないのか')
                        #print(Idx, Col)
           # else:
            #    if any(CorrDF.loc[Col,'Corr']) > 0 or CorrDF.loc[Col,'Corr'] > 0:
             #       color = 'red'
              #  else:
               #     color = 'blue'       
    except:
        #print(Idx, Col)
        color='black'
    return(color)
    
def mkEachEdgeColor(G):
    Gedge = G.edges(data=True)
    #特徴量とprop
    #CorrDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/RespIdx/Property_Fasting_AUC/Comb_Fasting_AUC_Eachother_Prop.xlsx',header=0,encoding = "ISO-8859-1")
    #最新相関係数
    CorrDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/全被験者_つなげて/RDF.xlsx',header=0,index_col=0)
    #CorrDF.index=Idx
     #Color_matrix作る
    #Idx = list(Adjacency_matrix.columns)
    #Color_matrix=pd.DataFrame(data=None,index=Idx,columns = Idx)
    edge_colorList=[]
    for (u,v,d) in Gedge:#全ノードの組み合わせが入ってる
            #if Adjacency_matrix.loc[jj, ii]==1:
            if CorrDF.loc[u][v] > 0:
                edge_colorList.append('sienna')#'black' 'orange' 'sienna'
            else:
                edge_colorList.append('indigo')#'blue' 'purple' 'indigo'
                
                """
                color = setColorResp(u,v,CorrDF)#関数で色を決める
                #Color_matrix.loc[u,v] =  color
                edge_colorList.append(color)
                """
    #Color_matrix.to_excel(save_dir+'ColorDF.xlsx')
    return(edge_colorList)   


 
    
def mkEachMolTime(Adjacency_matrix,save_dir):     #MolTimeで選ばれるやつだけ別に作る
    MolTimePropQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/QvalueFuncQvalueTimeSeriesProp_Eng.xlsx',header=0,encoding = "ISO-8859-1")
    LabelSum = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng.xlsx',header=0,encoding = "ISO-8859-1")
    Idx = list(Adjacency_matrix.columns)
    MolLabel = list(LabelSum.index)   
    Time = ['-10','0','10','20','30','45','60','75','90','120','150','180','210','240']
    
    #MolTimebのところが1だったら新たなIdxにする
    for i in MolLabel:
        for j in Time:
            if Adjacency_matrix.loc[j, i]==1:
                Adjacency_matrix.loc[j+'_'+i, i] = 1
                Adjacency_matrix.loc[i, j+'_'+i] = 1
                
    Adjacency_matrix.drop(Time, axis=0)               
    Adjacency_matrix.drop(Time, axis=1)  

    #Color_matrix作る
    Idx = list(Adjacency_matrix.columns)
    Color_matrix=pd.DataFrame(data=None,index=Idx,columns = Idx)
    
    for ii in Idx:
        for jj in Idx:
            if Adjacency_matrix.loc[jj, ii]==1:
                color = checkcolor(jj, ii)
                Color_matrix.loc[jj,ii] =  color
                
    Color_matrix.to_excel(save_dir+'ColorDF.xlsx')
    return(Adjacency_matrix)        
        
        
    
def CombMolTime_revMetaboMol(Adjacency_matrix):#MolTimeで選ばれるやつだけPropMolもつなぐ
    MolTimePropQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/QvalueFuncQvalueTimeSeriesProp_Eng.xlsx',header=0,encoding = "ISO-8859-1")
    
    LabelSum = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng.xlsx',header=0,encoding = "ISO-8859-1")

    Idx = list( MolTimePropQvalDF.index ); Col =list( MolTimePropQvalDF.columns )
    MetaboLabel = list(LabelSum['Metabo'])
    MolLabel = list(LabelSum.index)
    MetaboMolDict = dict( zip(MolLabel,MetaboLabel  ) )

    #一度でもそのMolTimeが選ばれたら描画する
    for i in Idx:#行
        for j in Col:#列
            r = re.compile("(.*)(_)(.*)") 
            d = r.search(i) 
            time = d.group(1)   
            Mol = d.group(3)             
            if MolTimePropQvalDF.loc[i, j] < 0.1:
                #時点だけに戻す

                    #if MolLabelList[jj] in ii:
                    #if MolLabelList[jj] == d.group(3):
                        #newlist.append(ii)
                        #print(d.group(3))
                Adjacency_matrix.loc[time,Mol] = 1
                Adjacency_matrix.loc[Mol,time] = 1 
                
            #elif LabelSum.loc[Mol, 'Metabo'] == MetaboMolDict[Mol]:
                Adjacency_matrix.loc[Mol,MetaboMolDict[Mol]] = 1
                Adjacency_matrix.loc[MetaboMolDict[Mol],Mol] = 1   
            elif     Mol == 'Glucose':
                print('yy')
            else:
                pass
                #Adjacency_matrix.loc[Mol,MetaboMolDict[Mol]] = 0                
                #Adjacency_matrix.loc[MetaboMolDict[Mol],Mol] = 0
    return(Adjacency_matrix) 
    
def CombMolTime(Adjacency_matrix) :#分子と時間をつなぐ
    MolTimePropQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/QvalueFuncQvalueTimeSeriesProp_Eng.xlsx',header=0,encoding = "ISO-8859-1")
    Idx = list( MolTimePropQvalDF.index ); Col =list( MolTimePropQvalDF.columns )
    #一度でもそのMolTimeが選ばれたら描画する
    for i in Idx:#行
        for j in Col:#列
            r = re.compile("(.*)(_)(.*)") 
            d = r.search(i) 
            time = d.group(1)   
            Mol = d.group(3)             
            if MolTimePropQvalDF.loc[i, j] < 0.1:
                #時点だけに戻す

                    #if MolLabelList[jj] in ii:
                    #if MolLabelList[jj] == d.group(3):
                        #newlist.append(ii)
                        #print(d.group(3))
                Adjacency_matrix.loc[time,Mol] = 1
                Adjacency_matrix.loc[Mol,time] = 1                
            else:
                pass
                #Adjacency_matrix.loc[time,Mol] = 0                
                #Adjacency_matrix.loc[Mol,time] = 0
    return(Adjacency_matrix)  
    
def CombMetaboMol(Adjacency_matrix) :#代謝と分子をつなぐ
    LabelSum = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng.xlsx',header=0,encoding = "ISO-8859-1")
    MetaboLabel = list(LabelSum['Metabo'])
    MolLabel = list(LabelSum.index)
    for i in MolLabel:#分子ごと回す
        for j in MetaboLabel:
            if LabelSum.loc[i, 'Metabo'] == j:
                Adjacency_matrix.loc[i,j] = 1
                Adjacency_matrix.loc[j,i] = 1                
            else:
                Adjacency_matrix.loc[i,j] = 0                
                Adjacency_matrix.loc[j,i] = 0
    return(Adjacency_matrix)              
    
def CombPropMetaboRespIdx(Adjacency_matrix):#分子の特徴量  ただしPropとつながっている分子特徴量のみ  
    PropMeaboQvalDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181023/Fasting_AUC_Gain_Eachother/QvalueStoreyParamACDE_ParamACDE.xlsx',header=0,encoding = "ISO-8859-1")
    PropMeaboQvalDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181025/Fasting_AUC_Eachother/QvalueStoreyParamACDE_ParamACDE.xlsx',header=0,encoding = "ISO-8859-1")

    pval= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/RespIdx',header=0,encoding = "ISO-8859-1")
    #NewDF=
    for i in list(PropMeaboQvalDF.index):#行だけ回す
        for j in list(PropMeaboQvalDF.columns):#列だけ回す
            if PropMeaboQvalDF.loc[i,j]< 0.1 or PropMeaboQvalDF.loc[j,i]< 0.1 :
                if  Adjacency_matrix.loc[i,'BMI':'TwoHInsulin'].any() == 1 and Adjacency_matrix.loc[j,'BMI':'TwoHInsulin'].any()==1:
                    Adjacency_matrix.loc[i,j] = 1
                    Adjacency_matrix.loc[j,i] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                    Adjacency_matrix.loc[j,i] = 0  
                #tempList
            else:
                Adjacency_matrix.loc[i,j] = 0
                Adjacency_matrix.loc[j,i] = 0                
    return(Adjacency_matrix)  
            
def CombPropMetabo(Adjacency_matrix):#プロパティと代謝をつなぐ
    PropMeaboQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/Fisher/Crossqval_Timedup.xlsx',header=0,encoding = "ISO-8859-1")
    PropMeaboQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180930/Crossqval_Timedup.xlsx',header=0,encoding = "ISO-8859-1")
    PropMeaboQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181012/Crossqval_Timedup_Fasvs2Hvs4H.xlsx',header=0,encoding = "ISO-8859-1")
    PropMeaboQvalDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjPropRespIdx/FastGain/QvalueStoreyParamACDE_ParamACDE_WGender.xlsx' ,header=0,encoding = "ISO-8859-1")
    PropMeaboQvalDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181023/Fasting_AUC_Gain/QvalueStoreyParamACDE_ParamACDE.xlsx' ,header=0,encoding = "ISO-8859-1")
    PropMeaboQvalDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181023/Fasting_AUC/QvalueStoreyParamACDE_ParamACDE.xlsx' ,header=0,encoding = "ISO-8859-1")

    for i in list(PropMeaboQvalDF.index):#行だけ回す
        for j in list(PropMeaboQvalDF.columns):#列だけ回す
            if PropMeaboQvalDF.loc[i,j] < 0.1:
                Adjacency_matrix.loc[i,j] = 1
                Adjacency_matrix.loc[j,i] = 1
            else:
                Adjacency_matrix.loc[i,j] = 0
                Adjacency_matrix.loc[j,i] = 0                
    return(Adjacency_matrix)     


def CombPropProp(Adjacency_matrix,PropSwich):#プロパティどうし相関有意なとこに1を入れてく：とりあえずp<0.05
    
    PropPDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/PDF_Prop_lower.xlsx',header=0,encoding = "ISO-8859-1")
    if PropSwich == 'P':
        for i in list(PropPDF.index):#行だけ回す
            for j in list(PropPDF.columns):#列だけ回す
                if PropPDF.loc[i,j] < 0.05:
                    Adjacency_matrix.loc[i,j] = 1
                    Adjacency_matrix.loc[j,i] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                    Adjacency_matrix.loc[j,i] = 0 
    elif PropSwich ==  'Q':
        PropPDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjPropRespIdx/Prop_PropQ_double.xlsx',header=0,encoding = "ISO-8859-1")
        for i in list(PropPDF.index):#行だけ回す
            for j in list(PropPDF.columns):#列だけ回す
                if PropPDF.loc[i,j] < 0.1:
                    Adjacency_matrix.loc[i,j] = 1
                    Adjacency_matrix.loc[j,i] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                    Adjacency_matrix.loc[j,i] = 0         
    return(Adjacency_matrix)                
    
def mkglobalDF(OptionDict,save_dir):##PropMetaboMolTimeをつなげるための行列を作成する
    PropSwich = 'Q';#'P'ならp値
    print('プロパティは'+PropSwich+'値で処理してます')
    PropCorrDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/CorrDF_Prop_lower.xlsx',header=0,encoding = "ISO-8859-1")#P値
    #PropMeaboQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/Fisher/Crossqval_Timedup.xlsx',header=0,encoding = "ISO-8859-1")
    PropMeaboQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20180930/Crossqval_Timedup.xlsx',header=0,encoding = "ISO-8859-1")
    PropMeaboQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181012/Crossqval_Timedup_Fasvs2Hvs4H.xlsx',header=0,encoding = "ISO-8859-1")
    PropMeaboQvalDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/SubjPropRespIdx/FastGain/QvalueStoreyParamACDE_ParamACDE_WGender.xlsx' ,header=0,encoding = "ISO-8859-1")

    MolTimePropQvalDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/MolTimevsProperty/QvalueFuncQvalueTimeSeriesProp_Eng.xlsx',header=0,encoding = "ISO-8859-1")
    LabelSum = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng.xlsx',header=0,encoding = "ISO-8859-1")

    
    
    Idx = list(PropCorrDF.columns)+ list(PropMeaboQvalDF.index) + list(LabelSum.index) + ['-10','0','10','20','30','45','60','75','90','120','150','180','210','240']
    Adjacency_matrix=pd.DataFrame(data=None,index=Idx,columns = Idx)
    Adjacency_matrix.to_excel(save_dir+'temp.xlsx')
    #隣接行列を埋めていく
    #まずはプロパティ同士
    Adjacency_matrix = CombPropProp(Adjacency_matrix,PropSwich)
    Adjacency_matrix.to_excel(save_dir+'temp_PropProp.xlsx')
    #そしてプロパティと代謝    
    Adjacency_matrix = CombPropMetabo(Adjacency_matrix)    
    Adjacency_matrix.to_excel(save_dir+'temp_PropMetabo.xlsx') 
    
    #そして代謝と分子
    #Adjacency_matrix = CombMetaboMol(Adjacency_matrix)    
    #Adjacency_matrix.to_excel(save_dir+'temp_MetaboMol.xlsx') 
    #最後に分子と時間
    #Adjacency_matrix = CombMolTime(Adjacency_matrix)    
    #Adjacency_matrix = CombMolTime_revMetaboMol(Adjacency_matrix)     #MolTimeで選ばれるやつだけPropMolもつなぐ
    #Adjacency_matrix.to_excel(save_dir+'temp_MolTime.xlsx')   
    
    #特徴量同士
    Adjacency_matrix = CombPropMetaboRespIdx(Adjacency_matrix)    
    Adjacency_matrix.to_excel(save_dir+'temp_PropMetaboRespIdx.xlsx')     
    
    if OptionDict['mkglobalDF_Time']=='mk':#時点も一緒に描画なら：Comb, 各分子の時点も一緒なら：'Each' 'mk'：作る
        Adjacency_matrix = mkEachMolTime(Adjacency_matrix,save_dir)     #MolTimeで選ばれるやつだけ別に作る
        Adjacency_matrix.to_excel(save_dir+'temp_EachMolTime.xlsx')
    
def mkEdgeColor(G,OptionDict):#参照するファイルごとに色をわける
    ref_file = OptionDict['ref_file'] #参照とするfile
    num_file =  len(ref_file)
    Gedge = G.edges(data=True)
    EdgeDense=[]
    
    CountDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/被験者ごと_平均/Adjacency_matrix_BH.xlsx',header=0,encoding = "ISO-8859-1")
    CountDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/被験者ごと_平均/Adjacency_matrix_Storey.xlsx',header=0,encoding = "ISO-8859-1")
    CountDF =  pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20181214/_Raw/CololForEdge.xlsx',header=0,encoding = "ISO-8859-1")


    edge_colorList=[]
    for (u,v,d) in Gedge:#ある相関係数以上の組みが（グラフ可される組み）入ってる
        try:
            EdgeWidth =+int(CountDF.loc[u,v])
        except:
            pass
            #print(u ,v )
        #print(int(CountDF.loc[u,v]))
        edge_color =  CountDF.loc[u,v]#ClstColorList[( 20 - int(CountDF.loc[u,v]) )]#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
        edge_colorList.append(edge_color) 
        sumnummin = []
        #DrawHistnumEdge(EdgenumSubj)
        #print(edge_color)
    return(Gedge,EdgeDense,edge_colorList, sumnummin)

def mkAdjacency_matrix_Comb(CorrDF,Thresh,path,OptionDict,save_dir):#何種類かのファイルを元に書かれた組み合わせの線をひく
    ref_file = OptionDict['ref_file'] #参照とするfile
    num_file =  len(ref_file)
        
    # networkxに引き渡すエッジリストを準備
    edge_lists = [];tempedge_list =[]
    deledge_lists = []#次消えるエッジ
    Combedge=[]; Combtemedge =[]  
    Adjacency_matrix=pd.DataFrame(data=None,index=CorrDF.index,columns = CorrDF.columns)
    for i in range(num_file):#fileだけ回す
        a = 1 if i<100 else 2 if i>100 else 0
        color = 'red' if i == 0  else 'blue' if  i == 1 else 'green'#CorrNoPramがblue
        for ii in range(len(ref_file[i].index)):#行の数だけ回す
            tmp = (ref_file[i].iloc[ii][0],ref_file[i].iloc[ii][1],0.5) #(from,to,weight)のタプルを作成
            edge_lists.append(tmp)  
            Adjacency_matrix.loc[ref_file[i].iloc[ii][0],ref_file[i].iloc[ii][1] ] = color
            Adjacency_matrix.loc[ref_file[i].iloc[ii][1],ref_file[i].iloc[ii][0] ] = color            
                    #組み合わせだけとをり出す
            #Combedge = [edge_lists[j][0:2] for j in range(len(edge_lists))]; Combtemedge = [tempedge_list[j][0:2] for j in range(len(tempedge_list))]
    Adjacency_matrix.to_excel(save_dir+'Adjacency_matrix.xlsx')
    
    return(deledge_lists,Combedge,edge_lists,Adjacency_matrix)
    
    
def mkAdjacency_matrix_Thresh(CorrDF,Thresh,path,OptionDict,save_dir):
    
    if OptionDict['mkEdge_double'] == 'CorrStd':#相関係数と標準偏差を元に線を引く
        SubjRstdrev = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/five/Rstd_rev_Eng.xlsx',header=0,encoding = "ISO-8859-1")
        ThreshStd = np.percentile( (SubjRstdrev.values.flatten()[~np.isnan(SubjRstdrev.values.flatten())] ),5)
        
        # 相関マトリクス作成
        corr_matrix = CorrDF
        Adjacency_matrix = pd.DataFrame(data=None,index=CorrDF.index,columns = CorrDF.columns)
        # networkxに引き渡すエッジリストを準備
        edge_lists = [];tempedge_list =[]
        deledge_lists = []#次消えるエッジ
        Combedge=[]; Combtemedge =[]
        for i in corr_matrix.index.values :
            for j in corr_matrix.index.values :
                #描画したいものだけを抽出する
                if (corr_matrix.loc[i,j] > Thresh) & (corr_matrix.loc[i,j] != 1) & (SubjRstdrev.loc[i,j] <= ThreshStd):#域値以上1ではない
                    tmp = (i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                    edge_lists.append(tmp)
                    #ついでに隣接行列を作成、吐き出す
                    Adjacency_matrix.loc[i,j] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                    
                if (corr_matrix.loc[i,j] > Thresh+0.05): #& (corr_matrix.loc[i,j] != 1) :#域値以上1ではない
                    tmpsec = (i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                    tempedge_list.append(tmpsec) 
                    #組み合わせだけとをり出す
                    Combedge = [edge_lists[j][0:2] for j in range(len(edge_lists))]; Combtemedge = [tempedge_list[j][0:2] for j in range(len(tempedge_list))]
                else:
                    Combedge.append(0); Combtemedge.append(0)
        if all([x == '' for x in Combedge]) & all([y == '' for y in Combtemedge]):
            deledge_lists = []   
        else:    
            deledge_lists = list([Combedge[jj] for jj in range(len(Combedge)) if Combedge[jj] not in  Combtemedge ])#次消えるエッジ
        
        Adjacency_matrix.to_excel(save_dir+'Adjacency_matrix.xlsx')        
    elif os.path.exists(path) == 0:#まだファイルがなければ、陽に与えた閾値を超えたもののみエッジはる
            
        # 相関マトリクス作成
        corr_matrix = CorrDF
        Adjacency_matrix = pd.DataFrame(data=None,index=CorrDF.index,columns = CorrDF.columns)
        # networkxに引き渡すエッジリストを準備
        edge_lists = [];tempedge_list =[]
        deledge_lists = []#次消えるエッジ
        Combedge=[]; Combtemedge =[]
        for i in corr_matrix.index.values :
            for j in corr_matrix.index.values :
                #描画したいものだけを抽出する
                if (corr_matrix.loc[i,j] > Thresh) & (corr_matrix.loc[i,j] != 1) :#域値過多1ではない
                    #print(i,j)
                    tmp = (i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                    edge_lists.append(tmp)
                    #ついでに隣接行列を作成、吐き出す
                    Adjacency_matrix.loc[i,j] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                    
                if (corr_matrix.loc[i,j] > Thresh+0.05): #& (corr_matrix.loc[i,j] != 1) :#域値以上1ではない
                    tmpsec = (i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                    tempedge_list.append(tmpsec) 
                    #組み合わせだけとをり出す
                    Combedge = [edge_lists[j][0:2] for j in range(len(edge_lists))]; Combtemedge = [tempedge_list[j][0:2] for j in range(len(tempedge_list))]
                else:
                    Combedge.append(0); Combtemedge.append(0)
        if all([x == '' for x in Combedge]) & all([y == '' for y in Combtemedge]):
            deledge_lists = []   
        else:
    
            deledge_lists = list([Combedge[jj] for jj in range(len(Combedge)) if Combedge[jj] not in  Combtemedge ])#次消えるエッジ
        
        Adjacency_matrix.to_excel(save_dir+'Adjacency_matrix.xlsx')
    else:
        Thresh=18
        CountDF =  pd.read_excel(path,header=0,encoding = "ISO-8859-1")
        MSMSH.mkhist(CountDF,save_dir,'Number of Subjects','SubjCount')

        # 相関マトリクス作成
        corr_matrix = CountDF 
        Adjacency_matrix = pd.DataFrame(data=None,index=CountDF .index,columns = CountDF .columns)
        # networkxに引き渡すエッジリストを準備
        edge_lists = [];tempedge_list =[]
        deledge_lists = []#次消えるエッジ
        Combedge=[]; Combtemedge =[]
        for i in corr_matrix.index.values :
            for j in corr_matrix.index.values :
                #描画したいものだけを抽出する
                if (corr_matrix.loc[i,j] > Thresh) & (corr_matrix.loc[i,j] != 1) :#域値以上1ではない
                    tmp = (i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                    edge_lists.append(tmp)
                    #ついでに隣接行列を作成、吐き出す
                    Adjacency_matrix.loc[i,j] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                    
                if (corr_matrix.loc[i,j] > Thresh+0.05): #& (corr_matrix.loc[i,j] != 1) :#域値以上1ではない
                    tmpsec = (i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                    tempedge_list.append(tmpsec) 
                    #組み合わせだけとをり出す
                    Combedge = [edge_lists[j][0:2] for j in range(len(edge_lists))]; Combtemedge = [tempedge_list[j][0:2] for j in range(len(tempedge_list))]
                else:
                    Combedge.append(0); Combtemedge.append(0)
        if all([x == '' for x in Combedge]) & all([y == '' for y in Combtemedge]):
            deledge_lists = []   
        else:
    
            deledge_lists = list([Combedge[jj] for jj in range(len(Combedge)) if Combedge[jj] not in  Combtemedge ])#次消えるエッジ            
    return(deledge_lists,Combedge,edge_lists,Adjacency_matrix)

def mkposmtr(DF,LabelName):
    #'Age','臨床指標','身体組成'の数が多いぶん重み付け
    r = re.compile("(.*)(_)(.*)") ; d = r.search(LabelName) ;
    #PosDict.update({'Height':0.75});PosDict.update({'Weight':0.9});PosDict.update({'BMI':1.05});PosDict.update({'Waist':1.2})

    wAge = np.sum(DF.loc[LabelName,'Age'])
    wBody = np.sum(DF.loc[LabelName,['BMI','Height', 'Weight', 'Waist']])
    wMed = np.sum(DF.loc[LabelName,['FastingGlucose','TwoHGlucose', 'FastingInsulin', 'TwoHInsulin']])
    #x = (np.random.uniform(-0.05, 0.05)*wAge + np.random.uniform(0.5, 1.0)*wBody + np.random.uniform(0.1, 0.5)*wMed) / (wAge + wBody + wMed)
    Age = DF.loc[LabelName,'Age']; BMI = DF.loc[LabelName,'BMI']; Height = DF.loc[LabelName,'Height']; Weight = DF.loc[LabelName,'Weight']; Waist = DF.loc[LabelName,'Waist']
    FasGlu = DF.loc[LabelName,'FastingGlucose']; FasIns = DF.loc[LabelName,'FastingInsulin']; THGlu = DF.loc[LabelName,'TwoHGlucose']; THIns = DF.loc[LabelName,'TwoHInsulin']; 
    x = ( 0.9*BMI + 0.75*Height + 1.05*Weight + 1.25*Waist + 0.45*FasGlu + 0.6*FasIns + 0.15*THGlu + 0.3*THIns) / (wAge + wBody + wMed) + np.random.uniform(-0.07, 0.07)
    if x != x:
        x=0
    try:
        if  d.group(3) in ['Growth hormone (GH)','SM-C IGF-1','Cortisol']:
            x = ( -0.1 * Age + 0.85*BMI + 1.2*Height + 1.1*Weight + 1.4*Waist) / (wAge + wBody + wMed)  + np.random.uniform(-0.07, 0.07)
    except:
        print('エッジの重み付けでえらー')
    return(x)        

def mkEdgeWidth(G,SubjectName,OptionDict,CorrDF):#エッジの太さを被験者と相関係数から決める
    Gedge = G.edges(data=True)
    
    ID = 0.7945723056980873 #if len(Gedge) == 122 else 0.8
    print(ID)
    EdgeDense=[]
    count=0
    edge_colorList=[]
    EdgenumSubj=[]
    ClstColorList=list(['darkred','red','sienna','salmon','orange',
                           'gold','goldenrod','yellow','green','olive','lime',
                           'purple','blue','slateblue','cyan','teal',
                           'magenta','pink','gray','black','white'])
    #ClstColorList=list(['darkred','red','purple','magenta','blue','green',
                    #    'pink',
                     #      'lime',                       
                      #     'gray','black'])
    #ClstColorList.reverse()
    
    
    #カウント用の DFつくる    

    path = '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/Eng/RCount_'+ str(ID) +'.xlsx'
    if os.path.exists(path) == 0:#まだファイルがなければ


        SubjR = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/Eng/RDF' + SubjectName[1] + '.xlsx',header=0,encoding = "ISO-8859-1")
        #SubjR = CNE.ChnageNametoEnglish(SubjR,2)
        CountDF = pd.DataFrame(data=None,index=list(SubjR.index),columns=list(SubjR.columns) )
        
        for (u,v,d) in Gedge:#ある相関係数以上の組みが（グラフ可される組み）入ってる
            EdgeWidth =0
            color_count=0
            for i in range(0,len(SubjectName)):
                SubjR = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/Eng/RDF' + SubjectName[i] +  '.xlsx',header=0,encoding = "ISO-8859-1")
                #SubjR = CNE.ChnageNametoEnglish(SubjR,2)
                #SubjR.to_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/Eng/RDF' + SubjectName[i] + '.xlsx')
                if SubjR.loc[u,v] > OptionDict['Thresh']:#その組みがある被験者でも相関係数0.8以上なら太くする
                    EdgeWidth += 0.0018   
                    color_count += 1      
            CountDF.loc[u,v] = color_count 
            CountDF.loc[v, u] = color_count 
            
            #edge_color =  cm.jet(color_count/20)#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
            #print(color_count)
            
            
            edge_color =  ClstColorList[color_count]#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
            
            """
            cmap = get_cmap('jet')
            
            norm = Normalize(vmin=21 - len(ClstColorList), vmax=20)
            mappable = ScalarMappable(cmap=cmap, norm=norm)
            mappable._A = []
    
        
            cmap = mpl.colors.ListedColormap(ClstColorList)
            fig=plt.figure()
            t=np.linspace(0,1,10)
            x = np.linspace(0,0,10)
            yy=np.linspace(0,0,10)
            ax = plt.scatter(x,yy,c=t,cmap=cmap,marker='.',lw=0)#vmin=min(pca_score),vmax=max(pca_score))
            #cbar = plt.colorbar(mappable)
            ticks = np.arange(norm.vmin, norm.vmax+1, 1)
            CF = plt.colorbar(ax)
            CF.set_ticks(ticks)
            #CF.set_ticklabels(ticks)
        
            #CF.ax.set_yticklabels([str(s) for s in ticks])
            plt.clim(( 20 - len(ClstColorList)) ,20)
            plt.savefig(save_dir+'colorbar.pdf')
            """
            
        #colorlist = ['blue','red','purple','orange','mediumspringgreen','greenyellow','lime','deepskyblue','cyan','magenta','gold']
            #RStd = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Rstd_Eng.xlsx',header=0,encoding = "ISO-8859-1")
            #RStd = RStd /RStd.max().max()
            d["weight"] = EdgeWidth
            #d["weight"] += int(0.5*RStd.loc[u,v])#cm.jet(RStd.loc[u,v])
            #edge_color =  cm.jet(RStd.loc[u,v]*2)#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
            #edge_color =  cm.jet(RStd.loc[u,v]*2)#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
    
                #EdgeDense += [EdgeWidth]
                
            edge_colorList.append(edge_color)    
            count+=1
            
            #print(len(Gedge) - count)
        CountDF.to_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Subjects/Eng/RCount'+str(OptionDict['Thresh']) +'.xlsx')
    
    elif OptionDict['Edge'] == 'Subject_rev': #各被験者でその相関係数以上の分子組み合わせがあるんだったら持ってくる
        #19人以上のエッジを赤く
        CountDF =  pd.read_excel(path,header=0,encoding = "ISO-8859-1")
        RstdDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Rstd_Eng.xlsx',header=0,encoding = "ISO-8859-1")
        CountList=[];RStdList=[]
        for (u,v,d) in Gedge:#ある相関係数以上の組みが（グラフ可される組み）入ってる
            #描画される分子における
            try:
                EdgeWidth =+int(CountDF.loc[u,v])
            except:
                print(u ,v )

            edge_color =  'red'  if CountDF.loc[u,v] > 18 else 'black' 
            #edge_color = 'red'  if RstdDF.loc[u,v] < 0.05 else 'black'

            EdgenumSubj.append(CountDF.loc[u,v])
            edge_colorList.append(edge_color) 
            sumnummin = CountDF.min().min()   
            
            CountList += [CountDF.loc[u,v]]
            RStdList += [RstdDF.loc[u,v]] 
            
        fig=plt.figure()
        ax = plt.scatter(CountList, RStdList)#,c=t,cmap=cmap,marker='.',lw=0)#vmin=min(pca_score),vmax=max(pca_score))
        plt.tick_params(labelsize=10)
        plt.xlabel('Count',fontsize=10);plt.ylabel('Std',fontsize=10)
        plt.title('相関係数>'+str(OptionDict['Thresh'])+'となる被験者の数と相関係数平均の標準偏差',fontsize=10)
        plt.savefig(save_dir+'StdVsCount'+str(OptionDict['Thresh']) +'.pdf')   
        plt.close()
        
    elif OptionDict['Edge'] == 'CorrCoef':#相関係数0.2, 0.4, 0.6, 0.8で太さ変える
        CorrDF = np.abs(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/全被験者_つなげて/RDF.xlsx',header=0,index_col=0))
### スピアマン
        CorrDF=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210609/TmCsCorrBetwMol/pearson/Delta/RDF.xlsx',header=0,index_col=0)
        
        for (u,v,d) in Gedge:#ある相関係数以上の組みが（グラフ可される組み）入ってる
            if CorrDF.loc[u,v] < 0.6:#0.5-0.6
                EdgeWidth = 0.25
            elif CorrDF.loc[u,v] < 0.7:#0.6-0.7
                EdgeWidth = 0.5
            elif CorrDF.loc[u,v] < 0.8:#0.7-0.8
                EdgeWidth = 0.25#1
            elif    CorrDF.loc[u,v] < 0.9:#0.8-0.9
                EdgeWidth = 3.5
            elif CorrDF.loc[u,v] < 1.0:#0.9-1.0
                EdgeWidth = 7
            else:
                EdgeWidth = 3  

            d["weight"] = EdgeWidth            
            #print(int(CountDF.loc[u,v]))
        edge_color = 0#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
            #print(CountDF.loc[u,v])
        edge_colorList=0
        sumnummin = 0
    elif OptionDict['Edge'] == 'CorrCoef_Tensor':#相関係数0.2, 0.4, 0.6, 0.8で太さ変える
        if OptionDict['CorrCoef_Tensor']== 'metabotimecorr':
            if OptionDict['Comp']==1:
                IntegDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210610/metabotimecorr/Comp1/tensor_reconstruction_Integration_metabometabotime_time0_metabo1py_Ext.xlsx',header=0,index_col=0) 
                ExtDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210610/metabotimecorr/SignificantHighMolQvaleBH_0_.xlsx',header=0,index_col=0)            
            elif OptionDict['Comp']==2:
                IntegDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210610/metabotimecorr/Comp2/tensor_reconstruction_Integration_metabometabotime_time1_metabo1py_Ext.xlsx',header=0,index_col=0) 
                ExtDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210610/metabotimecorr/SignificantHighMolQvaleBH_10_01.xlsx',header=0,index_col=0)            
            elif OptionDict['Comp']==3:
                IntegDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210610/metabotimecorr/Comp3/tensor_reconstruction_Integration_metabometabotime_time2_metabo1py_Ext.xlsx',header=0,index_col=0) 
                ExtDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210610/metabotimecorr/SignificantHighMolQvaleBH_230_023.xlsx',header=0,index_col=0)            
    
        elif OptionDict['CorrCoef_Tensor']== 'metabocorr':# 各時点の相関          
            IntegDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210607/5%/tensor_reconstruction_Integration_timemetabo_time0_metabo1py_Ext.xlsx',header=0,index_col=0) 
            #IntegDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210607/5%_timecomp2/tensor_reconstruction_Integration_timemetabo_time1_metabo1py_Ext.xlsx',header=0,index_col=0) 
            #IntegDF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210607/5%_timecomp3/tensor_reconstruction_Integration_timemetabo_time2_metabo1py_Ext.xlsx',header=0,index_col=0) 

            Ext12DF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210608/metaboCorr/SignificantHighMolQvaleBH_12.xlsx',header=0,index_col=0) 
            Ext13DF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210608/metaboCorr/SignificantHighMolQvaleBH_13.xlsx',header=0,index_col=0) 
            Ext16DF= pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210608/metaboCorr/SignificantHighMolQvaleBH_16.xlsx',header=0,index_col=0) 
        Target = list(set(ExtDF.index))#list(set(Ext12DF.index)|set(Ext13DF.index)|set(Ext16DF.index))
        IntegDF=IntegDF.loc[Target][Target]
        IntegList=list(IntegDF.index)
        IntegDF.index= ['GIP(Active)' if IntegList[i]=='GIP(Active)(pmol/L)'  else 'SM-C IGF-1' if IntegList[i]=='SM - C IGF - 1(ng/mL)' else 'C-peptide' if IntegList[i]== 'C peptide(ng/mL)' else 'Triglyceride' if IntegList[i]== 'TG (Neutral lipid)(mg/dL)'  else 'α-Amino-n-butyric acid' if IntegList[i]==  'a-ABA(μM)' else 'Growth hormone' if IntegList[i]== 'Growth hormone (GH)(ng/mL)' else delete_brackets(IntegList[i]) for i in range(len(IntegList)) ]    
        IntegDF.columns= ['GIP(Active)' if IntegList[i]=='GIP(Active)(pmol/L)'  else 'SM-C IGF-1' if IntegList[i]=='SM - C IGF - 1(ng/mL)' else 'C-peptide' if IntegList[i]== 'C peptide(ng/mL)' else 'Triglyceride' if IntegList[i]== 'TG (Neutral lipid)(mg/dL)'  else 'α-Amino-n-butyric acid' if IntegList[i]==  'a-ABA(μM)' else 'Growth hormone' if IntegList[i]== 'Growth hormone (GH)(ng/mL)' else delete_brackets(IntegList[i]) for i in range(len(IntegList)) ]    

        for (u,v,d) in Gedge:#ある相関係数以上の組みが（グラフ可される組み）入ってる
            if CorrDF.loc[u,v] < np.percentile(IntegDF.values,90):#0.5-0.6
                EdgeWidth = 0.25
            elif CorrDF.loc[u,v] < np.percentile(IntegDF.values,95):#0.6-0.7
                EdgeWidth = 0.50
            elif CorrDF.loc[u,v] < np.percentile(IntegDF.values,97.5):#0.7-0.8
                EdgeWidth = 1
            elif    CorrDF.loc[u,v] < np.percentile(IntegDF.values,99):#0.8-0.9
                EdgeWidth = 3
            elif CorrDF.loc[u,v] < np.percentile(IntegDF.values,100):#0.9-1.0
                EdgeWidth = 7
            else:
                EdgeWidth = 3  

            d["weight"] = EdgeWidth            
            #print(int(CountDF.loc[u,v]))
        edge_color = 0#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
            #print(CountDF.loc[u,v])
        edge_colorList=0
        sumnummin = 0
    else:    #各被験者でその相関係数以上の分子組み合わせがあるんだったら持ってくる
        CountDF =  pd.read_excel(path,header=0,encoding = "ISO-8859-1")
        
        for (u,v,d) in Gedge:#ある相関係数以上の組みが（グラフ可される組み）入ってる
            try:
                EdgeWidth =+int(CountDF.loc[u,v])
            except:
                print(u ,v )
            #print(int(CountDF.loc[u,v]))
            edge_color =  ClstColorList[( 20 - int(CountDF.loc[u,v]) )]#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
            #print(CountDF.loc[u,v])
            EdgenumSubj.append(CountDF.loc[u,v])
            edge_colorList.append(edge_color) 
            sumnummin = CountDF.min().min()
        #DrawHistnumEdge(EdgenumSubj)
        #print(edge_color)
    return(G,Gedge,EdgeDense,edge_colorList, sumnummin)


def delete_brackets(s):
    """
    括弧と括弧内文字列を削除
    """
    """ brackets to zenkaku """
    table = {
        "(": "（",
        ")": "）",
        "<": "＜",
        ">": "＞",
        "{": "｛",
        "}": "｝",
        "[": "［",
        "]": "］"
    }
    for key in table.keys():
        print(s)
        s = s.replace(key, table[key])
    """ delete zenkaku_brackets """
    l = ['（[^（|^）]*）', '【[^【|^】]*】', '＜[^＜|^＞]*＞', '［[^［|^］]*］',
         '「[^「|^」]*」', '｛[^｛|^｝]*｝', '〔[^〔|^〕]*〕', '〈[^〈|^〉]*〉']
    for l_ in l:
        s = re.sub(l_, "", s)
    """ recursive processing """
    return delete_brackets(s) if sum([1 if re.search(l_, s) else 0 for l_ in l]) > 0 else s   
    
def mkGraphQval(DF,QDF,QThresh,DegThresh,ColorSwitchDict,OptionDict,save_dir):#エッジの数を計算して、グラフを描画して,次数分布と次数中心性のグラフも、中心性の高いとこだけ描画も
    r = re.compile("(.*)(_)(.*)") 

    LabelSum =   pd.read_excel("//Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng_No3HB.xlsx",header=0,encoding = "ISO-8859-1",index_col=0)
    PropLabel =['BMI', 'Age', 'Weight', 'Waist', 'Height', 'TwoHGlucose', 'TwoHInsulin', 'FastingGlucose', 'FastingInsulin']
    for ijk in range(1):
            # 相関マトリクス作成
        corr_matrix = DF#相関行列でないと
        
        if OptionDict['TargetProp']!='':#BMIしか描かない
            corr_matrix[ DF.columns[DF.columns!=OptionDict['TargetProp']] ] = np.nan
            QDF[ QDF.columns[QDF.columns!=OptionDict['TargetProp']] ] = np.nan

        Adjacency_matrix = pd.DataFrame(data=None,index=DF.index,columns = DF.columns)
        PosDict={}
        # networkxに引き渡すエッジリストを準備
        edge_lists = []
        DFInd = list(corr_matrix.index)
        #if len(corr_matrix.index) == len(corr_matrix.columns):
        
        for i in corr_matrix.index.values :
            for j in corr_matrix.columns.values :
                #描画したいものだけを抽出する：Thresh以上
                if (QDF.loc[i,j] < QThresh) & (corr_matrix.loc[i,j] != 1) :
                    tmp = ( i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                    edge_lists.append(tmp)
                                    #ついでに隣接行列を作成、吐き出す
                    Adjacency_matrix.loc[i,j] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                        #座標を決める配列も作っておく
            #PosDict.update({i:mkposmtr(Adjacency_matrix,i)})
         
        #print(edge_lists)
        #else:
         #   None
        # 描画の準備
        G = nx.Graph()
    
        G.add_weighted_edges_from(edge_lists)
        
        plt.figure(figsize=(25,25))  #描画対象に合わせて設定する#MolTimePropは25,25; paramparamは15,15; paramTmeVarは12,9
        
        edge_colorList = ['black' for (u,v,d) in G.edges(data=True)]
        
        if OptionDict['GaraphSwitch'] == 'MolTimeProp': # 'MolTimePCProp'：MolTimeのPCとPropでの相関結果を元にグラフ描画

            #エッジの太さとかの設定
            edge_width = [((d["weight"]+1)*1.5)**2.1 for (u,v,d) in G.edges(data=True)]#[ 2 for (u,v,d) in G.edges(data=True)]#[ d["weight"]*100 for (u,v,d) in G.edges(data=True)]
            edge_width = [np.abs( ((d["weight"]+0)*1) )  for (u,v,d) in G.edges(data=True)]#[ 2 for (u,v,d) in G.edges(data=True)]#[ d["weight"]*100 for (u,v,d) in G.edges(data=True)]
        
            tag_list = []#ノードの大きさを決めるときに使う
            for (u,v,d) in G.edges(data=True):
                tag_list += [u] +[ v]
                #print([u,v,d])
        
            tag_count = collections.Counter(tag_list)#.most_common(50)
            aa=[];bb=[]
            for ii in range(len(tag_count)):
                aa+=[list(tag_count.keys())[ii]]
                bb += [{"count":list(tag_count.values())[ii]}]
            #print([aa,bb])
            G.add_nodes_from((aa[ij],bb[ij]) for ij in range(len(aa)))
            #G.add_nodes_from([(tag, {"count":count}) for tag,count in tag_count])
            #print([(tag, {"count":count}) for tag,count in tag_count])
            #ノードの大きさを設定する
            node_size = [ d["count"]*200 for (n,d) in G.nodes(data=True)]#[ 300 for (n,d) in G.nodes(data=True)]#[ d["count"]*200 for (n,d) in G.nodes(data=True)] # len(d)*2000 for (n,d) in G.nodes(data=True)]
            node_size = [ 1000 for (n,d) in G.nodes(data=True)]#[ 300 for (n,d) in G.nodes(data=True)]#[ d["count"]*200 for (n,d) in G.nodes(data=True)] # len(d)*2000 for (n,d) in G.nodes(data=True)]

           # print([ [n,d["count"] ]for (n,d) in G.nodes(data=True)])
            
        
            #np.random.seed(seed=1)    #ノードポジションの再現性を確保するためseedを設定する: MolTimePropは1, ParamTimeVarは590
        
            pos = nx.spring_layout(G)  if ijk == 0 else nx.shell_layout(G) #ノードのポジションの計算
            count=0;countb=0;
            #糖代謝、脂質代謝、アミノ酸代謝など分ける？
            Glu = list(LabelSum['English'])[0:13]
            Hormon = list(LabelSum['English'])[13:16]
            Lipid = list(LabelSum['English'])[16:29]
            AA = list(LabelSum['English'])[29:60]
            ION = list(LabelSum['English'])[60:66]
            Others = list(LabelSum['English'])[66:]
            
            edge_colorList = ['black' for (u,v,d) in G.edges(data=True)]
            #if 'BMI' in list(pos.keys()):#プロパティなので
            try:
                    for jj in list(pos.keys()):
                        r = re.compile("(.*)(_)(.*)") ; d = r.search(jj) ;
                        if jj in ['TwoHGlucose', 'TwoHInsulin', 'FastingGlucose', 'FastingInsulin']:
                            count+=1
                            x = 0.15*count; y = 0.15#np.random.uniform(0.1, 0.2)
                            #x=0;y=0
                            y=0
                            PosDict.update({'FastingGlucose':1});PosDict.update({'TwoHGlucose':2});PosDict.update({'FastingInsulin':3});PosDict.update({'TwoHInsulin':4})
                            
                            PosDict.update({jj:x})
                        elif jj in ['BMI', 'Weight', 'Waist', 'Height']:
                            countb+=1
                            x = 0.6+0.15*countb; y = 0.15#np.random.uniform(0.1, 0.2)
                            #x = 10;
                            y = 0
                            PosDict.update({'Height':-3});PosDict.update({'Weight':-2});PosDict.update({'BMI':0});PosDict.update({'Waist':-1})
                        elif jj == 'Age':
                            x = -0.05; y = 0.05#np.random.uniform(0.1, 0.2)
                            #x = 0.01; y = 0.15#np.random.uniform(0.1, 0.2)
                            y=0
                            PosDict.update({'Age':-3.5})
                            PosDict.update({jj:x})
                        elif  d.group(3) in AA :
                            x = np.random.uniform(-4, 0);y = np.random.uniform(0.2, 1.5)
                            #if d.group(3) in ['-10','0','10','20','30']:
                             #   pass
                            PosDict.update({jj:x})
                        elif  d.group(3) in Glu:
                            x = np.random.uniform(1.5, 2);y = np.random.uniform(0.2,1.5)
                            PosDict.update({jj:x})
                        elif  d.group(3) in Lipid:
                            x = np.random.uniform(0.1, 2);y = np.random.uniform(-0.2, -1)
                            PosDict.update({jj:x})
                        elif  d.group(3) in Others:
                            x = np.random.uniform(-3, -1);y = np.random.uniform(-1, -2)
                            PosDict.update({jj:x})
                        elif  d.group(3) in ION:
                            x = np.random.uniform(-0.2, -2);y = np.random.uniform(-1, 0)
                            PosDict.update({jj:x})
                        else :
                            x = np.random.uniform(-1.5, -0.1); y = np.random.uniform(-1, -0.1) #; y=np.where(yy >0.1 , np.random.uniform(0.2, 0.3), yy)
                            PosDict.update({jj:x})
                        pos[jj] =(PosDict[jj],y)  #これで各ノーどに座標を入れれる
            except:
                print('ポジションでエラー')
        #ある範囲内で乱数を与えるとか：臨床指標はnp.random.uniform(0, 0.3)：：身体組成は
        
        else:
            if type(OptionDict['GaraphSwitch_MolTime_PCA_Prop']) == 0: #濃度変動とPCとプロパティwoつなぐ
                # 相関マトリクス作成
                CoefDF = OptionDict['GaraphSwitch_MolTime_PCA_Prop']
                QDF = OptionDict['GaraphSwitch_MolTime_PCA_Prop_Q']
                CoefThresh = OptionDict['GaraphSwitch_MolTime_PCA_Prop_CoefThresh']
                corr_matrix = CoefDF#相関行列でないと
                if OptionDict['TargetProp']!='':#BMIしか描かない
                    corr_matrix[ DF.columns[DF.columns!=OptionDict['TargetProp']] ] = np.nan
                    QDF[ QDF.columns[QDF.columns!=OptionDict['TargetProp']] ] = np.nan
        
                QDF = QDF.drop(['PC4','PC5','PC6','PC7','PC8','PC9','PC10'])
                corr_matrix = corr_matrix.drop(['PC4','PC5','PC6','PC7','PC8','PC9','PC10'])
                Adjacency_matrix = pd.DataFrame(data=None,index=DF.index,columns = DF.columns)
                PosDict={}
                # networkxに引き渡すエッジリストを準備
                edge_lists = []
                edge_colorList = []
                DFInd = list(corr_matrix.index)
                #if len(corr_matrix.index) == len(corr_matrix.columns):
                
                for i in corr_matrix.index.values :
                    for j in corr_matrix.columns.values :
                        #描画したいものだけを抽出する：Thresh以上
                        if j in PropLabel :
                            if (QDF.loc[i,j] < QThresh) & (corr_matrix.loc[i,j] != 1):  
                                tmp = ( i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                                edge_lists.append(tmp)
                                                #ついでに隣接行列を作成、吐き出す
                                Adjacency_matrix.loc[i,j] = 1
                                edge_colorList.append('black')
                            else:
                                Adjacency_matrix.loc[i,j] = 0
                        else:
                            if OptionDict['GaraphSwitch_MolTime_PCA_Prop_Abs'] == 'Abs':#絶対値でLopding判断なら
                                 if (np.abs(QDF.loc[i,j] )> CoefThresh) & (corr_matrix.loc[i,j] != 1) :#Q値が０。１noPCor soukann takai PropTime
                                    tmp = ( i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                                    edge_lists.append(tmp)
                                                    #ついでに隣接行列を作成、吐き出す
                                    Adjacency_matrix.loc[i,j] = 1
                                    edge_colorList.append('black')
                                 else:
                                    Adjacency_matrix.loc[i,j] = 0 
                            else:#正負も考えるなら
                                 if (np.abs(QDF.loc[i,j] )> CoefThresh) & (corr_matrix.loc[i,j] != 1) :#Q値が０。１noPCor soukann takai PropTime
                                    tmp = ( i,j,corr_matrix.loc[i,j]*0.025) #(from,to,weight)のタプルを作成
                                    edge_lists.append(tmp)
                                                    #ついでに隣接行列を作成、吐き出す
                                    Adjacency_matrix.loc[i,j] = 1
                                    edge_color = 'red' if QDF.loc[i,j] >0 else 'blue'
                                    edge_colorList.append(edge_color)
                                 else:
                                    Adjacency_matrix.loc[i,j] = 0                                 
                                
                            #edge_colorList += [edge_color]
                            #座標を決める配列も作っておく
                    PosDict.update({i:mkposmtr(Adjacency_matrix,i)})
                    #edge_color = [ 'red'  if (u,v) in deledge_lists else 'black' for (u,v,d) in G.edges(data=True) ]
                    #print(len(edge_colorList) )
                #print(edge_lists)
                #else:
                 #   None
                # 描画の準備
                G = nx.Graph()
            
                G.add_weighted_edges_from(edge_lists)   
            #print(edge_lists)
               
            node_size = [  1000 for (n,d) in G.nodes(data=True)]#[ d["count"]*200 for (n,d) in G.nodes(data=True)] # len(d)*2000 for (n,d) in G.nodes(data=True)]
            edge_width = [2 for (u,v,d) in G.edges(data=True)]#[ d["weight"]*100 for (u,v,d) in G.edges(data=True)]
            font_size = [50 if 'PC' in u else 5 for (u,v,d) in G.edges(data=True) ]
            for (u,v,d) in G.edges(data=True):
                #print(u ,v, d)
                pass
            if OptionDict['pos']=='Left':#斥力で離す
                    pos = nx.spring_layout(G,k=1) 
            elif OptionDict['pos']=='circle':#サークルで描画 
                pos = nx.shell_layout(G);poscount=0
            else:
                poscount=0
                pos = nx.spring_layout(G) 
                for jj in list(pos.keys()):
    
                    poscount += 0.2
                    if 'PC' in jj:
                        pos[jj] = (0,poscount)
                    elif jj in PropLabel:
                        pos[jj] = (-0.005,poscount)
                        Propposx = [ -0.01, -0.005, -0.0075, -0.0025, -0.0070, -0.0045, -0.0020, -0.0025, -0.01]
                        Propposy = [ 8, 4.5, 12, 8, 6, 16, 14, 18, 20]
                        Propposx = dict(zip(PropLabel, Propposx)); Propposy = dict(zip(PropLabel, Propposy))
    
                        pos[jj] = (-0.005,poscount)
                        #x = np.random.uniform(-0.01,0);y = np.random.uniform(-0.005, 20)
                        x =Propposx[jj];  y = Propposy[jj] 
                        pos[jj] = ( x - 0.01, y ) 
                    elif OptionDict['pos']=='positive': #座標を陽に与える、'positive'                     
                        PosDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Position_20180912_2.xlsx',header=0,encoding = "ISO-8859-1")
                        PosDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/Property/2分子濃度変動相関/Position_20180912_2アミノ酸がいい感じに配置されるIncretin用.xlsx',header=0,encoding = "ISO-8859-1")
                        PropLabel =['BMI', 'Age', 'Weight', 'Waist', 'Height', 'TwoHGlucose', 'TwoHInsulin', 'FastingGlucose', 'FastingInsulin']        
                        if r.search(jj).group(3) in list(PosDF.columns):
                            pos[jj] = (PosDF[ r.search(jj).group(3) ][0],PosDF[ r.search(jj).group(3) ][1])
                        #elif jj in PropLabel:
                          
                    else:
                        pos[jj] = (0.005,poscount)
        NodeList = list(G.node);NodeList=[ NodeList[y] if '_' in NodeList[y] else   NodeList[y]+'_0' for y in range(len(NodeList))]
        if ColorSwitchDict  == 'ClstColor':#クラスタで色分け
            ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1),'ClstColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1) in list(LabelSum.index) else 'white' for x in range(len(NodeList))  ]
            FName = '_ClstColor'
        elif ColorSwitchDict  == 'MolColor':#分子で色分け
            #ColorList = [LabelSum.loc[NodeList[x],'MolColor'] if NodeList[x] in list(LabelSum.index) else 'white' for x in range(len(NodeList))  ]

            ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1),'MolColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1) in list(LabelSum.index) else 'white' for x in range(len(NodeList))  ]
            FName = '_MolColor'
        elif ColorSwitchDict == 'TimeVarColor':#時間ばらつきで色分け
             ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1),'TimeVarColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1) in list(LabelSum.index) else 'white' for x in range(len(NodeList))  ]
             FName = '_TimeVarColor' 
        
        #Adjacency_matrix.columns['Color']=ColorList
        Adjacency_matrix.to_excel(save_dir+'Adjacency_matrix.xlsx')
        nx.write_gml(G, save_dir+'network.gml')        
        cclname = '_Nml_' if ijk == 0 else '_Circle_'

  
        
        nx.draw_networkx_nodes(G, pos, node_size=node_size,node_color=ColorList,alpha=0.6)
        #nx.draw_networkx(G,pos=pos,font_size=15,alpha=0.8,edge_cmap=plt.cm.Greys,edge_vmin=-3e4,node_color=ColorList)
        nx.draw_networkx_labels(G,pos=pos,font_size=10,alpha=1,font_weight="bold")#node_color=ColorList,
        nx.draw_networkx_edges(G, pos, alpha=0.4, edge_color=edge_colorList, width=edge_width)
        plt.axis("off")
        #plt.savefig(save_dir+'default'+FName+cclname+'_Var' + OptionDict['GaraphSwitch_MolTime_PCA_Prop_Abs'] + '_' +str(OptionDict['GaraphSwitch_MolTime_PCA_Prop_CoefThresh']) + '.pdf')
        plt.savefig(save_dir+'default'+FName+cclname+'_Var' + '.pdf')

        plt.close()  
        
        Gdeg = dict(nx.degree(G))
        
            #次数の分布
        bins=20
        plt.figure()
        plt.hist(Gdeg.values(), bins=bins)
        #plt.yscale('log')
        #plt.xscale('log')
        plt.ylabel('Frequency')
        plt.xlabel('Degree')
        plt.savefig(save_dir+'DistributionOFK_' + str(bins) + '.pdf', dpi=300,bbox_inches="tight")
        #plt.close()
        
        
            # 次数中心性の分布を見る
        plt.figure()
        plt.hist(nx.degree_centrality(G).values(), bins=10)
        #plt.yscale('log')
        #plt.xscale('log')
        plt.xlabel('次数中心性')
        plt.ylabel('分子名数')
        plt.title('分子次数中心性の分布')
        plt.tight_layout()
        plt.savefig(save_dir+'dc.pdf', dpi=300,bbox_inches="tight")
        #plt.close()
        
        ################################## 媒介中心性
        between_centers = nx.betweenness_centrality(G)
                #print('媒介中心性')
        #print(sorted(between_centers.items(), key=lambda x: x[1], reverse=True))
        between_centers_sorted = sorted(between_centers.items(), key=lambda x: x[1], reverse=True)
        drawBarCent(between_centers_sorted,'Betweenness Centrality',save_dir)
        Centrality(G, Adjacency_matrix,save_dir)#中心性解析
        ################################## エッジの数 
        
        calcnumedge(Adjacency_matrix,between_centers_sorted,save_dir)      
        # 中心性の高いノードを抽出し描く
    """
    central_nodes = [k for k, v in nx.degree_centrality(G).items() if v >= DegThresh]
    G_central = G.subgraph(central_nodes)
    NodeList = list(G.node);NodeList=[ NodeList[y] if '_' in NodeList[y] else '0_'+NodeList[y] for y in range(len(NodeList))]

    if ColorSwitch  == 1:#クラスタで色分け
        ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'ClstColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
        FName = '_ClstColor'
    elif ColorSwitch  == 2:#分子で色分け
        ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'MolColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
        FName = '_MolColor'
    elif ColorSwitch == 3:#時間ばらつきで色分け
         ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'TimeVarColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
         FName = '_TimeVarColor' 
 
    draw_char_graph(G_central, save_dir+'central.pdf',node_color=ColorList)
    plt.axis("off")
    #plt.savefig(save_dir+'default'+FName+'.pdf')
    plt.close()  
    
    Gdeg = dict(nx.degree(G))
    
        #次数の分布
    plt.hist(Gdeg.values(), bins=15)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.ylabel('Frequency')
    plt.xlabel('Degree')
    plt.savefig(save_dir+'DistributionOFK.pdf', dpi=300)
    plt.close()
    
    
        # 次数中心性の分布を見る
    plt.hist(nx.degree_centrality(G).values(), bins=20)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel('次数中心性')
    plt.ylabel('分子名数')
    plt.title('分子次数中心性の分布')
    plt.tight_layout()
    plt.savefig(save_dir+'dc.pdf', dpi=300)
    plt.close()
    
        # 中心性の高いノードを抽出し描く
    central_nodes = [k for k, v in nx.degree_centrality(G).items() if v >= DegThresh]
    G_central = G.subgraph(central_nodes)
    NodeList = list(G.node);NodeList=[ NodeList[y] if '_' in NodeList[y] else '0_'+NodeList[y] for y in range(len(NodeList))]
    if ColorSwitch  == 1:#クラスタで色分け
        ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'ClstColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
        FName = '_ClstColor'

    elif ColorSwitch  == 2:#分子で色分け
        ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'MolColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
        FName = '_MolColor'

    elif ColorSwitch == 3:#時間ばらつきで色分け
         ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'TimeVarColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
         FName = '_TimeVarColor' 

    draw_char_graph(G_central, save_dir+'central'+FName+'.pdf',node_color=ColorList)
    
    
# グラフ可視化ベース関数
def draw_char_graph(G, fname, node_color='red', map_node_edge=True, **kwargs):
    #if map_node_edge:#変な色"
        #kwargs['node_size'] = np.array([d['freq'] for k, d in G.nodes(data=True)]) / 1e3
        #kwargs['edge_color'] = [d['freq'] for u, v, d in G.edges(data=True)]
    plt.figure(figsize=(10, 10))
    nx.draw(G,
            node_color=node_color,
            edge_cmap=plt.cm.Greys,
            edge_vmin=-3e4,
            width=0.5,
            with_labels=True,
            font_family='IPAexGothic',
            font_size=16,
            font_color='black',
            **kwargs)
    plt.savefig(fname, dpi=300)
    plt.show()
    """

def mkGraph(DF,CorrThresh,DegThresh,ColorSwitch,save_dir):#エッジの数を計算して、グラフを描画して,次数分布と次数中心性のグラフも、中心性の高いとこだけ描画も



    LabelSum =   pd.read_excel("//Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/file/FluctuationIndex/LabelSummary_Eng.xlsx",header=0,encoding = "ISO-8859-1")
    
        # 相関マトリクス作成
    corr_matrix = DF
    
    # networkxに引き渡すエッジリストを準備
    edge_lists = []
    
    #if len(corr_matrix.index) == len(corr_matrix.columns):
    
    for i in corr_matrix.index.values :
        for j in corr_matrix.columns.values :
            #描画したいものだけを抽出する：Thresh以上
            if (np.abs(corr_matrix.loc[i,j]) > CorrThresh) & (corr_matrix.loc[i,j] != 1) :
                tmp = (i,j,corr_matrix.loc[i,j]*0.05) #(from,to,weight)のタプルを作成
                edge_lists.append(tmp)
    #else:
     #   None
     
    # 描画の準備
    G = nx.Graph()

    G.add_weighted_edges_from(edge_lists)
    
    plt.figure(figsize=(18,18))  #描画対象に合わせて設定する
        #エッジの太さとかの設定
    edge_width = [ d["weight"]*(d["weight"]*500)**2.1 for (u,v,d) in G.edges(data=True)]

    tag_list = []
    for (u,v,d) in G.edges(data=True):
        tag_list += [u] +[ v]
        #print([u,v,d])

    tag_count = collections.Counter(tag_list).most_common(50)
    G.add_nodes_from([(tag, {"count":count}) for tag,count in tag_count])
    node_size = [ 300 for (n,d) in G.nodes(data=True)]#[ d["count"]*200 for (n,d) in G.nodes(data=True)]
#    print([ [n,d["count"] ]for (n,d) in G.nodes(data=True)])
    np.random.seed(seed=1)    #ノードポジションの再現性を確保するためseedを設定する
    #pos = nx.spring_layout(G)    #ノードのポジションの計算]
    pos = nx.shell_layout(G)
    NodeList = list(G.node);NodeList=[ NodeList[y] if '_' in NodeList[y] else '0_'+NodeList[y] for y in range(len(NodeList))]
    if ColorSwitch  == 1:#クラスタで色分け
        ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'ClstColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
        FName = '_ClstColor'
    elif ColorSwitch  == 2:#分子で色分け
        ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'MolColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
        FName = '_MolColor'
    elif ColorSwitch == 3:#時間ばらつきで色分け
         ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'TimeVarColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
         FName = '_TimeVarColor' 
    nx.draw_networkx(G,pos=pos,font_size=15,alpha=0.8,edge_cmap=plt.cm.Greys,edge_vmin=-3e4,node_color=ColorList)
    plt.axis("off")
    plt.savefig(save_dir+'default'+FName+'.pdf')
    plt.close()
    
    Gdeg = dict(nx.degree(G))
    
        #次数の分布
    plt.hist(Gdeg.values(), bins=5)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.ylabel('Frequency')
    plt.xlabel('Degree')
    plt.savefig(save_dir+'DistributionOFK.pdf', dpi=300)
    plt.close()
    
    
        # 次数中心性の分布を見る
    plt.hist(nx.degree_centrality(G).values(), bins=20)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel('次数中心性')
    plt.ylabel('分子名数')
    plt.title('分子次数中心性の分布')
    plt.tight_layout()
    plt.savefig(save_dir+'dc.pdf', dpi=300)
    plt.close()
    
        # 中心性の高いノードを抽出し描く
    central_nodes = [k for k, v in nx.degree_centrality(G).items() if v >= DegThresh]
    G_central = G.subgraph(central_nodes)
    NodeList = list(G.node);NodeList=[ NodeList[y] if '_' in NodeList[y] else '0_'+NodeList[y] for y in range(len(NodeList))]
    if ColorSwitch  == 1:#クラスタで色分け
        ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'ClstColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
        FName = '_ClstColor'

    elif ColorSwitch  == 2:#分子で色分け
        ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'MolColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
        FName = '_MolColor'

    elif ColorSwitch == 3:#時間ばらつきで色分け
         ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3),'TimeVarColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(3) in list(LabelSum.index) else 'black' for x in range(len(NodeList))  ]
         FName = '_TimeVarColor' 

    draw_char_graph(G_central, save_dir+'central'+FName+'.pdf',node_color=ColorList)