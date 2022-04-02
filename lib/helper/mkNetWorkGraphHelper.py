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

import xml.etree.ElementTree as ET 

"""Betweenness centrality measures."""
from heapq import heappush, heappop
from itertools import count
import random

import networkx as nx



def mkNetWorkGraph(CorrDF,LabelSum,Thresh,ColorSwitchDict,OptionDict,SubjectName,save_dir):
    r = re.compile("(.*)(_)(.*)") 
    import networkx as nx
    import collections
    import itertools
    from statistics import mean, median,variance,stdev
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    ID = 0.75
    path = ''
    
    try:
        Adjacency_matrix = OptionDict['Adjacency_matrix']
        deledge_lists,Combedge,edge_lists = mkEdge_Color(Adjacency_matrix, OptionDict,save_dir)
    except:
        if OptionDict['mkEdge'] == 'Thresh':
            deledge_lists,Combedge,edge_lists,Adjacency_matrix = mkAdjacency_matrix_Thresh(CorrDF,Thresh,path,OptionDict,save_dir)

        elif OptionDict['mkEdge'] == 'Comb':
            deledge_lists,Combedge,edge_lists,Adjacency_matrix = mkAdjacency_matrix_Comb(CorrDF,Thresh,path,OptionDict,save_dir)

    G = nx.Graph()
    G.add_weighted_edges_from(edge_lists)
    #####################################################
    between_centers =  Centrality(G,Adjacency_matrix, save_dir)
    #####################################################
    plt.figure(figsize=(20,22))
    plt.rcParams["font.size"] = 5
    #############################################
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

        edge_color = [ 'red'  if (u,v) in deledge_lists or (v,u) in deledge_lists else 'black' for (u,v,d) in G.edges(data=True) ]
    edge_color = [ 'black' for (u,v,d) in G.edges(data=True) ]

###################################
###################################  
    Gedge = G.edges(data=True)
    MetaboDict=dict(zip(list(LabelSum.index),list(LabelSum['Metabo'])))
    MetabonumDict=dict(zip(['glucose','Hormone','Lipid','AminoAcid','Ion','Others'],[13,3,12,31,6,18]))
    MolEdgeDF = pd.DataFrame(data=np.array([0]*83),index=list(Adjacency_matrix.columns),columns=['numEdge'])
    for (u,v,d) in Gedge:
        MolEdgeDF.loc[u]+=1/MetabonumDict[MetaboDict[v]]
        MolEdgeDF.loc[v]+=1/MetabonumDict[MetaboDict[u]]
#########################################
    if OptionDict['EdgeColor']=='Color_posneg': 
        EdgeDense=0.1
        pos = nx.shell_layout(G)
        Idx= list(pos.keys());Col= list(pos.keys());
        edge_color = mkEachEdgeColor(G,CorrDF)   
############################################################################ 
    tag_list = []
    for (u,v,d) in G.edges(data=True):
        tag_list += [u] +[ v]
    tag_count = collections.Counter(tag_list)
    aa=[];bb=[]
    for ii in range(len(tag_count)):
        aa+=[list(tag_count.keys())[ii]]
        bb += [{"count":list(tag_count.values())[ii]}]
    G.add_nodes_from((aa[ij],bb[ij]) for ij in range(len(aa)))
    node_size = [ d["count"]*100 for (n,d) in G.nodes(data=True)]
    node_size = [2000 for (n,d) in G.nodes(data=True) ]
    if OptionDict['Node_size'] == 'Centrality':#
        node_size = [500000**(v+0.55) for v in between_centers.values()]
#################################################################### 
    np.random.seed(seed=10000)   
    pos = nx.shell_layout(G)
            
    if OptionDict['pos']=='positive':
        PosDF = pd.read_excel("./Data/NetworkPsition.xlsx",header=0,index_col=0)
        labels={}
        for jj in list(pos.keys()):
            labels[jj]=jj
            pos[jj] = (PosDF[jj][0],PosDF[jj][1])
    elif OptionDict['pos']=='positive_Prop': 
        for jj in list(pos.keys()):
            pos[jj] = (PosDF[jj][0],PosDF[jj][1])
    elif OptionDict['pos']=='circle':
        pos = nx.shell_layout(G)  
    elif OptionDict['pos']=='Left':
        pos = nx.spring_layout(G,k=0.5)      
    numNode = nx.number_of_nodes(G) 
    NumClst = nx.number_connected_components(G)
    
    EachClstnode = sorted(nx.connected_components(G), key = len, reverse=True)
    numEachClstnode = [];VarnumEachClstnode=[]
    for i in range(len(EachClstnode)):
        numEachClstnode += [len(EachClstnode[i])]
    AvenumEachClstnode = mean(numEachClstnode) ;
    VarnumEachClstnode = variance(numEachClstnode) if len(numEachClstnode)>1 else 0
####################################################################
    ColorList=[]
    if 'Acetoacetate' in list(G.nodes) or 'GLP-1' in list(G.nodes) or 'β- aminoisobutyric acid' in list(G.nodes) or '1-methyl-histidine' in list(G.nodes) or 'Aspartic acid' in list(G.nodes) or 'M ethanolamine' in list(G.nodes):
        LabelSum.loc['Acetoacetate','MolColor'] = 'green'
        LabelSum.loc['GLP-1','MolColor'] = 'red'
        LabelSum.loc['β- aminoisobutyric acid','MolColor'] = 'blue'
        LabelSum.loc['1-methyl-histidine','MolColor'] = 'blue'
        LabelSum.loc['Aspartic acid','MolColor'] = 'blue'
        LabelSum.loc['M ethanolamine','MolColor'] = 'blue'       
    if ColorSwitchDict  == 'ClstColor':
        nodeList = list(G.nodes);ColorList = [LabelSum.loc[NodeList[x],'ClstColor'] for x in range(len(NodeList))]
    elif ColorSwitchDict  == 'MolColor':
        NodeList = list(G.nodes);
        try:
            if re.compile("(.*)(_)(.*)").search(NodeList[0]).group(1) in list(LabelSum.index) :
                ii=re.compile("(.*)(_)(.*)").search(NodeList[0]).group(1)
                ColorList = [LabelSum.loc[re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1),'MolColor'] if re.compile("(.*)(_)(.*)").search(NodeList[x]).group(1) in list(LabelSum.index) else 'white' for x in range(len(NodeList))  ]
        except:
            pass
        ColorList = [LabelSum.loc[NodeList[x],'MolColor'] for x in range(len(NodeList))]

    elif ColorSwitchDict == 'TimeVarColor':
        NodeList = list(G.nodes);ColorList = [LabelSum.loc[NodeList[x],'TimeVarColor'] for x in range(len(NodeList))]
    elif ColorSwitchDict == 'RespIdx':
        NodeList = list(G.nodes);
        for i in NodeList:
            if i in PropLabel:
                    ColorList.append('white')            
            else:         
                    ColorList.append(LabelSum.loc[ r.search(i).group(3),'MolColor'] )#if  r.search(NodeList[x]).group(3) in list(LabelSum.index)  else 'white'  for x in range(len(NodeList)) ]      
    else:
        ColorList=['white']*len(G.nodes)    
    nx.draw_networkx_nodes(G, pos, node_size=node_size,node_color=ColorList,alpha=0.5)  
  
    if OptionDict['pos_Global'] == 'LabelLoop':     
        nx.draw_networkx_labels(G,pos=pos,font_size=12,alpha=1,font_weight="bold")#node_color=ColorList,
    elif OptionDict['pos_Global']=='Fasting_AUC_Gain':
        labels={}

        for idx, node in enumerate(G.nodes()):
            labels[node] = FAGLabel[idx]
        nx.draw_networkx_labels(G,pos=pos,labels=labels,font_size=15,alpha=1,font_weight="bold")#node_color=ColorList,
    else:
        nx.draw_networkx_labels(G,pos=pos,labels=labels,font_size=16,alpha=1,font_weight="bold")#node_color=ColorList,        
    nx.draw_networkx_edges(G, pos, alpha=0.5,edge_color=edge_color, width=edge_width)
    plt.axis("off")
    numsub = 20 - len(set(edge_color))
    plt.title('Mean(CorrelationCoefficient)>'+('%02.2f' %Thresh)+'_#OfComponents：'+str(NumClst),size='30')   
    plt.savefig(save_dir+'NewtWork_' +('%02.2f' %Thresh)+'_'+ str(NumClst) + ColorSwitchDict+'.pdf',bbox_inches="tight")
    ##plt.savefig(save_dir+'NewtWork_' +('%02.2f' %Thresh)+'_'+ str(NumClst) + ColorSwitchDict+'.png',bbox_inches="tight")

    Gdeg = dict(nx.degree(G))
    MolEdgeDF.to_excel(save_dir+'DegNormalized.xlsx')

    return([NumClst,AvenumEachClstnode,VarnumEachClstnode, numNode,MolEdgeDF])



def Centrality(G, Adjacency_matrix,save_dir):
    between_centers = nx.betweenness_centrality(G)
    between_centers = {k: v for k, v in between_centers.items() if v!= 0}

    between_centers_sorted = sorted(between_centers.items(), key=lambda x: x[1], reverse=True)
    drawBarCent(between_centers_sorted,'Betweenness Centrality', save_dir)
    mkDF(between_centers_sorted,'Betweenness Centrality', save_dir)
    
    calcnumedge(Adjacency_matrix,between_centers_sorted,save_dir)
    
    close_centers = nx.closeness_centrality(G)
    close_centers_sorted = sorted(close_centers.items(), key=lambda x: x[1], reverse=True)
    ##drawBarCent(close_centers_sorted,'Closeness Centrality', save_dir)
    

    degree_centers = nx.degree_centrality(G)
    degree_centers_sorted = sorted(degree_centers.items(), key=lambda x: x[1], reverse=True)
    ##drawBarCent(degree_centers_sorted,'Degree Centrality', save_dir)

    return(between_centers)    

def mkDF(between_centers,file_name,save_dir):
    valueList = [  between_centers[i][1]  for i in range(len(between_centers) ) ]
    nameList =     [ between_centers[i][0] for i in range(len(between_centers) ) ]

    NewDF = pd.DataFrame(data=valueList,columns=[file_name],index=nameList)
    NewDF.T.to_excel(save_dir+ file_name + '_DF.xlsx')
    

def drawBarCent(between_centers,file_name,save_dir):
    valueList = [  between_centers[i][1]  for i in range(len(between_centers) ) ]
    nameList =     [ between_centers[i][0] for i in range(len(between_centers) ) ]
    x = np.linspace(0, len(between_centers),len(between_centers) )
    fig = plt.figure()
    plt.bar( x, valueList  )
    #plt.yscale("symlog")
    plt.xticks(x,nameList,rotation=270 ,size=8)
    plt.gca().yaxis.set_tick_params(labelsize=20)
    
    plt.title(file_name)
    plt.rcParams['axes.linewidth'] = 1.5 
    
    plt.savefig(save_dir+ file_name + '.pdf',bbox_inches="tight")
    ##plt.savefig(save_dir+ file_name + '.png',bbox_inches="tight")
    plt.close()
    
def calcnumedge(DF,between_centers,save_dir):
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
    
    plt.savefig(save_dir+'BarnumEdge.pdf',bbox_inches="tight")  
    plt.close()
    


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


def _accumulate_basic(betweenness, S, P, sigma, s):
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            betweenness[w] += delta[w]
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
    
    edge_colorList=[]
    for (u,v,d) in Gedge:
        edge_color =  CountDF.loc[u,v]#ClstColorList[( 20 - int(CountDF.loc[u,v]) )]#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
        edge_colorList.append(edge_color) 
        sumnummin = []
        EdgeWidth=0.3

        
    return(Gedge,EdgeDense,edge_colorList, sumnummin)
    
def mkEdge_Color(Adjacency_matrix, OptionDict,save_dir):
    Idx = Adjacency_matrix
    edge_lists = [];tempedge_list =[]
    deledge_lists = []
    Combedge=[]; Combtemedge =[]     
    for i in Idx:#
        for j in Idx:
            if Adjacency_matrix.loc[i,j] == 1:
                tmp = (i,j,0.5) 
                edge_lists.append(tmp)  
    ##Adjacency_matrix.to_excel(save_dir+'Color_matrix.xlsx')
    
    return(deledge_lists,Combedge,edge_lists)
    
def mkEachEdgeColor(G,CorrDF):
    Gedge = G.edges(data=True)

    edge_colorList=[]
    for (u,v,d) in Gedge:
            #if Adjacency_matrix.loc[jj, ii]==1:
            if CorrDF.loc[u][v] > 0:
                edge_colorList.append('sienna')#'black' 'orange' 'sienna'
            else:
                edge_colorList.append('indigo')#'blue' 'purple' 'indigo'
    ##Color_matrix.to_excel(save_dir+'ColorDF.xlsx')
    return(edge_colorList)   
   

def mkAdjacency_matrix_Comb(CorrDF,Thresh,path,OptionDict,save_dir):
    ref_file = OptionDict['ref_file'] 
    num_file =  len(ref_file)
    edge_lists = [];tempedge_list =[]
    deledge_lists = []
    Combedge=[]; Combtemedge =[]  
    Adjacency_matrix=pd.DataFrame(data=None,index=CorrDF.index,columns = CorrDF.columns)
    for i in range(num_file):#
        a = 1 if i<100 else 2 if i>100 else 0
        color = 'red' if i == 0  else 'blue' if  i == 1 else 'green'
        for ii in range(len(ref_file[i].index)):
            tmp = (ref_file[i].iloc[ii][0],ref_file[i].iloc[ii][1],0.5) 
            edge_lists.append(tmp)  
            Adjacency_matrix.loc[ref_file[i].iloc[ii][0],ref_file[i].iloc[ii][1] ] = color
            Adjacency_matrix.loc[ref_file[i].iloc[ii][1],ref_file[i].iloc[ii][0] ] = color            
    Adjacency_matrix.to_excel(save_dir+'Adjacency_matrix.xlsx')
    
    return(deledge_lists,Combedge,edge_lists,Adjacency_matrix)
    
    
def mkAdjacency_matrix_Thresh(CorrDF,Thresh,path,OptionDict,save_dir):
    
    if OptionDict['mkEdge_double'] == 'CorrStd':
        ThreshStd = np.percentile( (SubjRstdrev.values.flatten()[~np.isnan(SubjRstdrev.values.flatten())] ),5)
        
        corr_matrix = CorrDF
        Adjacency_matrix = pd.DataFrame(data=None,index=CorrDF.index,columns = CorrDF.columns)

        edge_lists = [];tempedge_list =[]
        deledge_lists = []
        Combedge=[]; Combtemedge =[]
        for i in corr_matrix.index.values :
            for j in corr_matrix.index.values :
                if (corr_matrix.loc[i,j] > Thresh) & (corr_matrix.loc[i,j] != 1) & (SubjRstdrev.loc[i,j] <= ThreshStd):
                    tmp = (i,j,corr_matrix.loc[i,j]*0.025) #(from,
                    edge_lists.append(tmp)

                    Adjacency_matrix.loc[i,j] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                    
                if (corr_matrix.loc[i,j] > Thresh+0.05):
                    tmpsec = (i,j,corr_matrix.loc[i,j]*0.025) 
                    tempedge_list.append(tmpsec) 
                    Combedge = [edge_lists[j][0:2] for j in range(len(edge_lists))]; Combtemedge = [tempedge_list[j][0:2] for j in range(len(tempedge_list))]
                else:
                    Combedge.append(0); Combtemedge.append(0)
        if all([x == '' for x in Combedge]) & all([y == '' for y in Combtemedge]):
            deledge_lists = []   
        else:    
            deledge_lists = list([Combedge[jj] for jj in range(len(Combedge)) if Combedge[jj] not in  Combtemedge ])#次消えるエッジ
        
        Adjacency_matrix.to_excel(save_dir+'Adjacency_matrix.xlsx')        
    elif os.path.exists(path) == 0:
            
        corr_matrix = CorrDF
        Adjacency_matrix = pd.DataFrame(data=None,index=CorrDF.index,columns = CorrDF.columns)
        edge_lists = [];tempedge_list =[]
        deledge_lists = []
        Combedge=[]; Combtemedge =[]
        for i in corr_matrix.index.values :
            for j in corr_matrix.index.values :
                if (corr_matrix.loc[i,j] > Thresh) & (corr_matrix.loc[i,j] != 1) :
                    #print(i,j)
                    tmp = (i,j,corr_matrix.loc[i,j]*0.025) 
                    edge_lists.append(tmp)
                    Adjacency_matrix.loc[i,j] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                    
                if (corr_matrix.loc[i,j] > Thresh+0.05): 
                    tmpsec = (i,j,corr_matrix.loc[i,j]*0.025) 
                    tempedge_list.append(tmpsec) 
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

        corr_matrix = CountDF 
        Adjacency_matrix = pd.DataFrame(data=None,index=CountDF .index,columns = CountDF .columns)
        # 
        edge_lists = [];tempedge_list =[]
        deledge_lists = []
        Combedge=[]; Combtemedge =[]
        for i in corr_matrix.index.values :
            for j in corr_matrix.index.values :
                if (corr_matrix.loc[i,j] > Thresh) & (corr_matrix.loc[i,j] != 1) :
                    tmp = (i,j,corr_matrix.loc[i,j]*0.025) 
                    edge_lists.append(tmp)
                    Adjacency_matrix.loc[i,j] = 1
                else:
                    Adjacency_matrix.loc[i,j] = 0
                    
                if (corr_matrix.loc[i,j] > Thresh+0.05): #& (corr_matrix.loc[i,j] != 1) :
                    tmpsec = (i,j,corr_matrix.loc[i,j]*0.025) 
                    tempedge_list.append(tmpsec) 
                    Combedge = [edge_lists[j][0:2] for j in range(len(edge_lists))]; Combtemedge = [tempedge_list[j][0:2] for j in range(len(tempedge_list))]
                else:
                    Combedge.append(0); Combtemedge.append(0)
        if all([x == '' for x in Combedge]) & all([y == '' for y in Combtemedge]):
            deledge_lists = []   
        else:
    
            deledge_lists = list([Combedge[jj] for jj in range(len(Combedge)) if Combedge[jj] not in  Combtemedge ])#次消えるエッジ            
    return(deledge_lists,Combedge,edge_lists,Adjacency_matrix)

     

def mkEdgeWidth(G,SubjectName,OptionDict,CorrDF):
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
     path = ''

    if OptionDict['Edge'] == 'Subject_rev': 
        CountDF =  pd.read_excel(path,header=0,encoding = "ISO-8859-1")
        CountList=[];RStdList=[]
        for (u,v,d) in Gedge:#
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
        plt.savefig(save_dir+'StdVsCount'+str(OptionDict['Thresh']) +'.pdf')   
        plt.close()
        
    elif OptionDict['Edge'] == 'CorrCoef':
        
        for (u,v,d) in Gedge:
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
 
    else:    
        CountDF =  pd.read_excel(path,header=0,encoding = "ISO-8859-1")
        
        for (u,v,d) in Gedge:
            try:
                EdgeWidth =+int(CountDF.loc[u,v])
            except:
                pass
            #print(int(CountDF.loc[u,v]))
            edge_color =  ClstColorList[( 20 - int(CountDF.loc[u,v]) )]#'red'  if RStd.loc[u,v] > 0.3 else 'black' 
            EdgenumSubj.append(CountDF.loc[u,v])
            edge_colorList.append(edge_color) 
            sumnummin = CountDF.min().min()
    return(G,Gedge,EdgeDense,edge_colorList, sumnummin)


def delete_brackets(s):

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
    

