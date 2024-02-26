#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 06:01:19 2022

@author: as822

goal: 
    - quantify Nll, Nlh and Nhh interactions over time
    - build webs of this connectivity at different timepoints
    - fit exponential function to this (or calculate the critical mass assuming 
                                        that the max for an individual country 
                                        will be the max for others proportional 
                                        to # neurosurgeons)

"""

# =============================================================================
# import libraries
# =============================================================================
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
# from matplotlib.patches import ConnectionStyle
import sys
from itertools import permutations
import seaborn as sns
import matplotlib as mpl
import math
import os
# import statsmodels.stats as stats
import pingouin as pg
np.random.seed(0)
plt.rcParams["font.family"] = "serif"
# =============================================================================
# define functions
# =============================================================================

def fulladj(df, Ucountry): # get the full adjacency matrix for all the data
    if len(df.shape)==2:# have more than 1 paper
        df_year=df.loc[df["Publication year"]>0]
        
        adjacency=np.zeros([len(Ucountry),len(Ucountry)])
        
        countries=df_year[authcont];
        
        for idx in df_year.index:
            uset=countries.loc[idx][:].dropna().unique()
            iset=np.empty(len(uset),);
            for j in range(len(uset)):
                iset[j]=Ucountry.index(uset[j]);
            perm = permutations(iset,2);# connections between two countries
            for i in perm:
                adjacency[int(i[0]),int(i[1])]=adjacency[int(i[0]),int(i[1])]+1;

    else:# only one paper
        df_year=df;
        adjacency=np.zeros([len(Ucountry),len(Ucountry)])
        
        countries=df_year[authcont];
        
        uset=countries.dropna().unique()
        iset=np.empty(len(uset),);
        for j in range(len(uset)):
            iset[j]=Ucountry.index(uset[j]);
        perm = permutations(iset,2);# connections between two countries
        for i in perm:
            adjacency[int(i[0]),int(i[1])]=adjacency[int(i[0]),int(i[1])]+1;

    # print(df_year)

    A = pd.DataFrame(adjacency, index = Ucountry, columns = Ucountry);
    return adjacency, A

def getadj(df, year, Ucountry): # get the adjacency matrix for an individual year
    df_year=df.loc[df["Publication year"]==year]
    adjacency=np.zeros([len(Ucountry),len(Ucountry)])
    
    countries=df_year[authcont];
    
    for idx in df_year.index:
        uset=countries.loc[idx][:].dropna().unique();
        iset=np.empty(len(uset),);
        for j in range(len(uset)):
            iset[j]=Ucountry.index(uset[j]);
        perm = permutations(iset,2);# connections between two countries
        for i in perm:
            adjacency[int(i[0]),int(i[1])]=adjacency[int(i[0]),int(i[1])]+1;

    A = pd.DataFrame(adjacency, index = Ucountry, columns = Ucountry);
    return adjacency, A

def getarcs(A, Ucountry, lmicbins): # get data requisite to plot the arcs
    nLMIC=lmicbins.LMICever.sum();
    nHIC=len(lmicbins)-nLMIC;
    
    # pit HIC on the left, LMIC on the right
    thetaHIC=np.linspace(start = 135, stop = 225, num=nHIC)
    thetaLMIC=np.linspace(start = -65, stop = 65, num=nLMIC)
    
    xH=np.cos(thetaHIC*math.pi/180)
    yH=np.sin(thetaHIC*math.pi/180)
    xL=np.cos(thetaLMIC*math.pi/180)
    yL=np.sin(thetaLMIC*math.pi/180)

    iL=0;
    iH=0;
    G = nx.Graph()
    for i in range(len(lmicbins)):
        if lmicbins.LMICever[i]==1:
            G.add_node(i,
                        # country=Ucountry[i],
                       country=lmicbins.Country[i],
                       color='red',
                       setting='LMIC',
                       pos=(xL[iL],yL[iL]),# pos=(xL[len(xL)-iL-1],yL[len(yL)-iL-1]),
                       ha='left')
            
            iL=iL+1;
            
        if lmicbins.LMICever[i]==0:
            G.add_node(i,
                       country=lmicbins.Country[i],
                       color='blue',
                       setting='HIC',
                       pos=(xH[iH],yH[iH]),
                       ha='right')            
            iH=iH+1;


    for i in range(len(lmicbins)):# go through again and add edges
        icountry=np.where(np.array(Ucountry)==lmicbins.Country[i])[0][0]
        # print(icountry)
        temp=list(A.iloc[icountry][:]>0);
        # print(temp)
        # plt.plot(temp)
        idx=[s for s, x in enumerate(temp) if x];

        val=list(A.iloc[icountry][idx])
        # print('here')
        for j in range(len(idx)):
            G.add_edge(i,idx[j],weight=val[j]);
            
    return G
    
def plotarc(G,fig,ax): #plot the arcs

    fig.patch.set_facecolor('gray')
    
    # get the position of each element
    pos = nx.get_node_attributes(G,'pos')
    nx.draw_networkx_nodes(G,pos,node_size=5,ax=ax)
    weights=list(nx.get_edge_attributes(G,'weight').values())

    # get edges
    # edgeidx=0;
    # normweight=weights/np.max(weights);# how it should be
    normweight=weights/np.double(10);# just for visualization relative to max
    cmap = mpl.cm.get_cmap('magma_r')

    edgeprintorder=np.argsort(weights)
    edgelist=list(G.edges());
    

    for idxval in edgeprintorder:
        edge=edgelist[idxval]
        
        source, target = edge
        if lmicbins.LMICever[edge[1]]==lmicbins.LMICever[edge[0]]: #from same setting
            rad = 0.7
        else:
            rad = 0.05
            
        arrowprops=dict(arrowstyle="-",
                        color=cmap(normweight[idxval]),
                        connectionstyle=f"arc3,rad={rad}",
                        linestyle= '-',
                        linewidth=1,# linewidth=np.log(weights[edgeidx])+0.5,
                        alpha=1,)
        ax.annotate("",
                xy=pos[source],
                xytext=pos[target],
                arrowprops=arrowprops
               )
        # edgeidx=edgeidx+1;
    
    # norm = mpl.colors.Normalize(vmin=1,vmax=np.max(weights))# this is adjusted to plot multiple 
    #temporarily make vmax 10
    norm = mpl.colors.Normalize(vmin=1,vmax=10)# this is adjusted to plot multiple 
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    plt.axis('off')
    
    labels=nx.get_node_attributes(G,'country');
    ha=nx.get_node_attributes(G,'ha');
    nx.draw_networkx_labels(G, pos,labels=labels,font_size=12,font_color='black')
    cax = plt.axes([0.9, 0.1, 0.03, 0.8])
    mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,orientation='vertical')
    
def Ncounter(A,lmicbins):
    iLMIC=np.where(lmicbins.LMICever==1)[0]
    iHIC=np.where(lmicbins.LMICever==0)[0]
    
    Nlh=A.iloc[iLMIC,iHIC].sum().sum()
    Nhh=A.iloc[iHIC,iHIC].sum().sum()/2
    Nll=A.iloc[iLMIC,iLMIC].sum().sum()/2
    
    return Nlh, Nhh, Nll

def Ncounter_norm(A,lmicbins):
    iLMIC=np.where(lmicbins.LMICever==1)[0]
    iHIC=np.where(lmicbins.LMICever==0)[0]
    
    Nlh=A.iloc[iLMIC,iHIC].sum().sum()
    Nhh=A.iloc[iHIC,iHIC].sum().sum()/2
    Nll=A.iloc[iLMIC,iLMIC].sum().sum()/2
    
    # denom=(A.sum()>0).sum();
    # if denom==0:
        # denom=1;# correct zero denom
    
    # Nlh_norm=Nlh/denom
    # Nll_norm=Nll/denom
    # Nhh_norm=Nhh/denom
    
    if (Nlh+Nhh+Nll)>0:
        Nlh_norm=Nlh/(Nlh+Nhh+Nll)
        Nll_norm=Nll/(Nlh+Nhh+Nll)
        Nhh_norm=Nhh/(Nlh+Nhh+Nll)
    else:
        Nlh_norm=np.nan
        Nll_norm=np.nan
        Nhh_norm=np.nan
    
    return Nlh_norm, Nhh_norm, Nll_norm


# =============================================================================
# load the raw data
# =============================================================================
df=pd.read_csv('Final_Global Nsurg Bibliometrics - Global Nsurg Bibliometrics--Final -- 220810_funding_numbers3.csv',header=1)
authcont = [col for col in df if col.endswith('country')]

lmicbins=pd.read_csv('UpdatedLMICBins.csv')

pubyears=np.sort(df['Publication year'].unique())# get unique publication years
yrs=pubyears[:-1];# exclude the last year

df[["Nlh","Nhh","Nll"]]=np.nan
df["F_none"]=np.nan;

fund = [col for col in df if col.startswith('F_')]

# =============================================================================
# Get key elements of the data
# =============================================================================

Ucountry=list((pd.unique(df[authcont].values.ravel('K'))))
Ucountry = [x for x in Ucountry if str(x) != 'nan']
Ucountry=list(np.sort(Ucountry))
# Ucountry.remove("New Zealand")# only has publications by itself with no collab

# =============================================================================
# define save array
# =============================================================================
savearr = [1,1,1,1,1,];
# savearr=[0,0,0,0,0]; #2,3,4,5 save
# savearr=[0,0,0,0,0]; #2,3,4,5 save

# =============================================================================
# do some processing for df
# =============================================================================
for i in range(len(df)):
    adjacency, A=fulladj(df.loc[i,:],Ucountry);
    df.loc[i,"Nlh"],df.loc[i,"Nhh"],df.loc[i,"Nll"]=Ncounter(A,lmicbins);
    

df["F_none"]=(df[fund].sum(axis=1,skipna=True)==0).astype(int)

fund = [col for col in df if col.startswith('F_')]# redefine
fund_true=['Industry Funded','Government Funded','University Funded',
           'Charity Funded','Not Funded']

# =============================================================================
# Get raw # for manuscript
# =============================================================================
# print(df[fund].sum())
# adjacency, A=fulladj(df.loc[df['Publication year']<2021,:], Ucountry);
# Nlh,Nhh,Nll = Ncounter(A,lmicbins);

df_1985_2020=df.loc[df["Publication year"]<2021];
adjacency, A=fulladj(df_1985_2020, Ucountry);
Nlh,Nhh,Nll = Ncounter(A,lmicbins);

df=df_1985_2020;# just to make it clean

# =============================================================================
# Figure 2
# =============================================================================
if savearr[0]==1:
    
    figname="F2_2-N_yrs";
    
    #init df
    N=pd.DataFrame(data=None,index=yrs,columns=["Nlh","Nhh","Nll"])
    
    # get # manuscripts by year
    count=pd.DataFrame(data=None,index=yrs,columns=["count"])
    
    for i in yrs:
        adjacency, A=getadj(df,i, Ucountry);
        N.loc[i,"Nlh"],N.loc[i,"Nhh"],N.loc[i,"Nll"]=Ncounter(A,lmicbins)
        count.loc[i,"count"] = (df["Publication year"]==i).sum()
    
    c_Nll="#4ba173";
    c_Nlh="#9900ff";
    c_Nhh="#980000ff";
    
    fig, ax = plt.subplots(ncols=2, nrows=3, sharex='col',figsize=(10,10),
                           sharey='row',gridspec_kw={'width_ratios': [1.7, 1]})
  
    ### A
    ax[0,0].plot(yrs,count['count'].to_list(),'-o',color='k',linewidth=1.5,label='count');
    ax[0,0].legend(loc='upper center',fontsize=16,ncol=3,)
    
    ax[0,0].set_xticks(np.arange(np.min(yrs),np.max(yrs),1),minor=True,
                     color='k',alpha=0.2,linewidth=0.5);
    ax[0,0].set_xticks(np.arange(np.min(yrs),np.max(yrs)+5,5),minor=False,
                     color='k',alpha=0.5,linewidth=1);
    l,r = ax[0,0].get_xlim()
    b,t = ax[0,0].get_ylim()
    ax[0,0].text(0.02*(r-l)+l,t-0.08*(t-b),'A',fontsize=20,ha='left',va='top',
               backgroundcolor="#d6e1ff")
    
    ax[0,0].set_ylabel('number of manuscripts by year')
    ax[0,0].grid(True, 'major',color='k',alpha=0.3,linewidth=1)
    ax[0,0].grid(True, 'minor',color='k',alpha=0.2,linewidth=0.5)
    
    ### B
    ax[0,1].plot(yrs,count['count'].to_list(),'-o',color='k',linewidth=1.5,label='count');
    # ax[0,1].legend(loc='upper center',fontsize=16,ncol=3,)
    
    ax[0,1].set_xticks(np.arange(np.min(yrs),np.max(yrs),1),minor=True,
                      color='k',alpha=0.2,linewidth=0.5);
    ax[0,1].set_xticks(np.arange(np.min(yrs),np.max(yrs)+5,5),minor=False,
                      color='k',alpha=0.5,linewidth=1);
    ax[0,1].set_xlim(left=2009.5,right=2020.5);
    l,r = ax[0,1].get_xlim()
    b,t = ax[0,1].get_ylim()
    ax[0,1].text(0.04*(r-l)+l,t-0.08*(t-b),'B',fontsize=20,ha='left',va='top',
               backgroundcolor="#d6e1ff")
    
    ax[0,1].grid(True, 'major',color='k',alpha=0.3,linewidth=1)
    ax[0,1].grid(True, 'minor',color='k',alpha=0.2,linewidth=0.5)
    

    ### C
    ax[1,0].plot(N.index,N.Nlh,'-o',color=c_Nlh,linewidth=1.5,label=r'$N_{l,h}$')
    ax[1,0].plot(N.index,N.Nll,'-o',color=c_Nll,linewidth=1.5,label=r'$N_{l,l}$')
    ax[1,0].plot(N.index,N.Nhh,'-o',color=c_Nhh,linewidth=1.5,label=r'$N_{h,h}$')
    ax[1,0].legend(loc='upper center',fontsize=16,ncol=3,)
    ax[1,0].set_xticks(np.arange(np.min(yrs),np.max(yrs),1),minor=True,
                     color='k',alpha=0.2,linewidth=0.5);
    ax[1,0].set_xticks(np.arange(np.min(yrs),np.max(yrs)+5,5),minor=False,
                     color='k',alpha=0.5,linewidth=1);
    l,r = ax[1,0].get_xlim()
    b,t = ax[1,0].get_ylim()
    ax[1,0].text(0.02*(r-l)+l,t-0.08*(t-b),'C',fontsize=20,ha='left',va='top',
               backgroundcolor="#d6e1ff")
    ax[1,0].set_ylabel('global coauthorship connections')
    ax[1,0].grid(True, 'major',color='k',alpha=0.3,linewidth=1)
    ax[1,0].grid(True, 'minor',color='k',alpha=0.2,linewidth=0.5)
    
    
    ### D
    ax[1,1].plot(N.index,N.Nlh,'-o',color=c_Nlh,linewidth=1.5,label=r'$N_{l,h}$')
    ax[1,1].plot(N.index,N.Nll,'-o',color=c_Nll,linewidth=1.5,label=r'$N_{l,l}$')
    ax[1,1].plot(N.index,N.Nhh,'-o',color=c_Nhh,linewidth=1.5,label=r'$N_{h,h}$')
    ax[1,1].set_xticks(np.arange(np.min(yrs),np.max(yrs),1),minor=True,
                     color='k',alpha=0.2,linewidth=0.5);
    ax[1,1].set_xticks(np.arange(np.min(yrs),np.max(yrs)+5,5),minor=False,
                     color='k',alpha=0.5,linewidth=1);
    ax[1,1].set_xlim(left=2009.5,right=2020.5);
    
    l,r = ax[1,1].get_xlim()
    b,t = ax[1,1].get_ylim()
    ax[1,1].text(0.04*(r-l)+l,t-0.08*(t-b),'D',fontsize=20,ha='left',va='top',
               backgroundcolor="#d6e1ff")
    ax[1,1].grid(True, 'major',color='k',alpha=0.3,linewidth=1)
    ax[1,1].grid(True, 'minor',color='k',alpha=0.2,linewidth=0.5)
    
    
    ### E
    melt=pd.melt(df,value_vars=["Nlh","Nhh","Nll"],id_vars=["Publication year"])
    melt=melt.loc[melt["Publication year"]<2021,:]# avoid 2021 cus incomplete
    sns.lineplot(data=melt,x="Publication year",y="value",hue="variable", 
                 ax=ax[2,0],legend=None,palette=[c_Nlh,c_Nhh,c_Nll],)
    
    
    ax[2,0].set_xticks(np.arange(np.min(yrs),np.max(yrs),1),minor=True,
                     color='k',alpha=0.2,linewidth=0.5);
    ax[2,0].set_xticks(np.arange(np.min(yrs),np.max(yrs)+5,5),minor=False,
                     color='k',alpha=0.5,linewidth=1);
    l,r = ax[2,0].get_xlim()
    b,t = ax[2,0].get_ylim()
    ax[2,0].text(0.02*(r-l)+l,t-0.08*(t-b),'E',fontsize=20,ha='left',va='top',
               backgroundcolor="#d6e1ff")
    ax[2,0].set_ylabel('average individual connections')
    ax[2,0].grid(True, 'major',color='k',alpha=0.3,linewidth=1)
    ax[2,0].grid(True, 'minor',color='k',alpha=0.2,linewidth=0.5)
    ax[2,0].set_ylim(b,t)
    
    ### F
    sns.lineplot(data=melt,x="Publication year",y="value",hue="variable", 
                 ax=ax[2,1],legend=None,palette=[c_Nlh,c_Nhh,c_Nll],)

    ax[2,1].set_xticks(np.arange(np.min(yrs),np.max(yrs),1),minor=True,
                     color='k',alpha=0.2,linewidth=0.5);
    ax[2,1].set_xticks(np.arange(np.min(yrs),np.max(yrs)+5,5),minor=False,
                     color='k',alpha=0.5,linewidth=1);
    ax[2,1].set_xlim(left=2009.5,right=2020.5);
    
    ax[2,1].grid(True, 'major',color='k',alpha=0.3,linewidth=1)
    ax[2,1].grid(True, 'minor',color='k',alpha=0.2,linewidth=0.5)
    
    l,r = ax[2,1].get_xlim()
    b,t = ax[2,1].get_ylim()
    ax[2,1].text(0.04*(r-l)+l,t-0.08*(t-b),'F',fontsize=20,ha='left',va='top',
               backgroundcolor="#d6e1ff")
    plt.tight_layout()
    fig.savefig(os.getcwd()+"/figures/"+figname+".png", dpi=600)
    fig.savefig(os.getcwd()+"/figures/"+figname+".svg")
    fig.savefig(os.getcwd()+"/figures/"+figname+".pdf")



# =============================================================================
# Figure 3
# =============================================================================
if savearr[1]==1:
    # melt=pd.melt(df,value_vars=["Nlh","Nhh","Nll"],id_vars=["Publication year"]+fund)
    figname="F3_4-FnoF"
    temp=df.loc[df["Publication year"]<2021,["Publication year","Nlh","Nhh","Nll"]+fund]# avoid 2021 cus incomplete
    
    temp['idx']=np.arange(0,len(temp));
    
    melt=pd.melt(temp,value_vars=["Nlh","Nhh","Nll"],id_vars=["Publication year","idx"]+fund)
    melt=melt.loc[melt["Publication year"]<2021,:]# avoid 2021 cus incomplete
    
    melt.loc[melt['F_none']==1,'F_none']="Not Funded"
    melt.loc[melt['F_none']==0,'F_none']="Funded"
    
    from scipy.stats import bartlett
    
    # check bartlett
    stat,p = bartlett(temp.Nlh,temp.Nll,temp.Nhh)# do bartlett test to check if pop
    # don't have equal variance... low p means unequal variance
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bartlett.html
    file = open("figures/"+figname+'.txt',"w")
    print("bartlett p = {p:08.20f}".format(p=p))
    file.write("\nbartlett test\n")
    file.write("bartlett p = {p:08.20f} \n".format(p=p))
    
   
    
    # do Welch stuff - only compare within Nll, Nlh, Nhh because comparing 
    # outside doesn't really provide good info
    # https://www.geeksforgeeks.org/how-to-perform-welchs-anova-in-python/
    welch_Nll=pg.welch_anova(data=temp,dv='Nll',between='F_none')
    print("\nwelch Nll grouped by binary funding: \n")
    print(welch_Nll.to_string())
    file.write("\nwelch Nll grouped by binary funding: \n")
    file.write(welch_Nll.to_string())    
    
    welch_Nlh=pg.welch_anova(data=temp,dv='Nlh',between='F_none')
    print("\nwelch Nlh grouped by binary funding: \n")
    print(welch_Nlh.to_string())
    file.write("\n\nwelch Nlh grouped by binary funding: \n")
    file.write(welch_Nlh.to_string())
    
    welch_Nhh=pg.welch_anova(data=temp,dv='Nhh',between='F_none')
    print("\nwelch Nhh grouped by binary funding: \n")
    print(welch_Nhh.to_string())
    file.write("\n\nwelch Nhh grouped by binary funding: \n")
    file.write(welch_Nhh.to_string())
    
    welch_p=pd.DataFrame(data=[welch_Nlh["p-unc"][0],welch_Nhh["p-unc"][0],
                            welch_Nll["p-unc"][0]], index=['Nlh','Nhh','Nll'],
                              columns=['pval']);
    
    # df_rm_anova=pg.rm_anova(data=melt,dv='value',within=['F_none','variable'],
    #                         subject='idx')
    
    # df_ancova=pg.ancova(data=temp,dv='Nlh',covar=['Nll','Nhh'],between='F_none',)
    
    # do ANOVA to have it, show cannot do multicomparison without bonferroni correction
    # df_anova=pg.anova(data=melt,dv='value',between=['variable','F_none'],
    #                   detailed=True,ss_type=2)
    # print("\n2-way ANOVA for all factors:\n")
    # print(df_anova.to_string())
    # file.write("\n2-way ANOVA for all factors:\n")
    # file.write(df_anova.to_string())
    
    gs_kw = dict(width_ratios=[1, 1], height_ratios=[2, 1])
    # fig, ax = plt.subplots(figsize=(8,12),ncols=2,nrows=2,gridspec_kw=gs_kw)
    fig, ax = plt.subplot_mosaic([['upper','upper'],
                                  ['ll','lr']],
                                 figsize=(4,6),gridspec_kw=gs_kw,
                                 sharey=False,
                                 layout='constrained')
    ax1=ax['ll']
    ax2=ax['lr']
    ax=ax['upper']
    # ax1=ax[1][0];
    # ax2=ax[1][1];
    # ax=ax[0][:];
    
    sns.barplot(data=melt, x="variable",y="value",hue="F_none",
                ax=ax,n_boot=1000,seed=0)#,errorbar=('ci,95'))
    
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[0:], labels=labels[0:],loc=8,
              bbox_to_anchor=(0.5,-0.20),ncol=2)
    ax.set_xticklabels([r'$\overline{n}_{l,h}$',
                        r'$\overline{n}_{h,h}$',
                        r'$\overline{n}_{l,l}$'])
    
    # ax.set_xticklabels=['$r\bar{n_{l,h}}$','r\n_{h,h}','r\n_{l,l}']
    ax.set_xlabel('')
    ax.set_ylabel('average connectivity')
    
    l,u=ax.get_ylim()
    
    p=0.024;
    n=(welch_p.pval<0.05).sum()
    P=(1-p*(n+2));
    maxy=u;
    u=(P*l+maxy-l)/P-p*(u-l)
    ax.set_ylim(l,u)
    bound=p*(u-l)
    y0=maxy;
    
    for j,pval in enumerate(welch_p.pval):
        if pval<0.05:
  
            numstar=1;
            
            if p<0.001:
                numstar=3;
            elif p<0.01:
                numstar=2;
            
            x0=j-0.2
            x1=j+0.2
            # xmid=(x1+x0)/2;
            # xmid=np.min([x0,x1])+0.5
            xmid=j
            
            # marker_style.update(markeredgecolor="none", markersize=15)
            
            ax.plot([x0,x1],[y0,y0],'-|',color='k',linewidth=1)
            
            dx=0.1
            squaresize=5;
            starsize=5;
            
            if numstar==2:
                ax.plot(xmid-dx/2,y0,marker='o',color='w',markersize=squaresize)
                ax.plot(xmid+dx/2,y0,marker='o',color='w',markersize=squaresize)
                # ax[i_ax,1].plot(xmid,y0,marker='o',color='w',markersize=squaresize)
                ax.plot(xmid-dx/2,y0,marker=(6, 2, 0),color='k',markersize=starsize)
                ax.plot(xmid+dx/2,y0,marker=(6, 2, 0),color='k',markersize=starsize)
            elif numstar==1:
                ax.plot(xmid,y0,marker='o',color='w',markersize=squaresize)
                ax.plot(xmid,y0,marker=(6, 2, 0),color='k',markersize=starsize)              
            else:
                ax.plot(xmid,y0,marker='o',color='w',markersize=squaresize)
                ax.plot(xmid-dx,y0,marker='o',color='w',markersize=squaresize)
                ax.plot(xmid+dx,y0,marker='o',color='w',markersize=squaresize)
                ax.plot(xmid,y0,marker=(6, 2, 0),color='k',markersize=starsize)              
                ax.plot(xmid-dx,y0,marker=(6, 2, 0),color='k',markersize=starsize) 
                ax.plot(xmid+dx,y0,marker=(6, 2, 0),color='k',markersize=starsize)                
            
            GH=pg.pairwise_gameshowell(data=temp,dv=welch_p.index[j],
                                       between='F_none',)
            
            print("\nGames Howell pairwise test: \n")
            print(GH.to_string())
            file.write("\n\nGames Howell pairwise test: \n")
            file.write(GH.to_string())
            
            print(GH)
    

    
    df_anova=pg.anova(data=temp,dv='Nlh',between=['F_Gov','F_Industry','F_Charity'],detailed=True,ss_type=3)
    
    print("\nn-way ANOVA on Nlh:\n")
    print(df_anova.to_string())
    file.write("\n\nn-way ANOVA on Nlh:\n")
    file.write(df_anova.to_string())
    
    
    temp.loc[(temp.F_Gov==0),"F_Gov"]="Not Funded"
    temp.loc[(temp.F_Gov==1),"F_Gov"]="Funded"
    
    temp.loc[(temp.F_Charity==0),"F_Charity"]="Not Funded"
    temp.loc[(temp.F_Charity==1),"F_Charity"]="Funded"
    
    
    # pgov=df_anova.loc[0,'p-unc']
    # pgov=pg.pairwise_gameshowell(data=temp,dv='Nlh',between='F_Gov',).pval
    sns.barplot(data=temp, x="F_Gov",y="Nlh",#hue="F_Gov",
                ax=ax1,n_boot=1000,seed=0,order=['Funded','Not Funded'])#,errorbar=('ci,95'))
    ax1.set_ylabel(r'$\overline{n}_{l,h}$')
    ax1.set_xlabel('Government Funding')
    ax1.set_xticklabels(['',''])
    ax1.legend([],[], frameon=False)
    ax1.set_ylim(0,15)
    y0=14.2
    ax1.plot([0,1],[y0,y0],'-|',color='k',linewidth=1)
    xmid=0.5;
    dx=0.2
    ptemp=pg.pairwise_tukey(data=temp,dv='Nlh',between='F_Gov',)['p-tukey'][0]
    
    ax0=ax
    ax=ax1
    if ptemp<0.001:
        ax.plot(xmid,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid-dx,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid+dx,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid,y0,marker=(6, 2, 0),color='k',markersize=starsize)              
        ax.plot(xmid-dx,y0,marker=(6, 2, 0),color='k',markersize=starsize) 
        ax.plot(xmid+dx,y0,marker=(6, 2, 0),color='k',markersize=starsize)             
        
    elif ptemp<0.01:
        ax.plot(xmid-dx/2,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid+dx/2,y0,marker='o',color='w',markersize=squaresize)
        # ax[i_ax,1].plot(xmid,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid-dx/2,y0,marker=(6, 2, 0),color='k',markersize=starsize)
        ax.plot(xmid+dx/2,y0,marker=(6, 2, 0),color='k',markersize=starsize)
        
    elif ptemp<0.05:
        ax.plot(xmid,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid,y0,marker=(6, 2, 0),color='k',markersize=starsize)                  
    
    ax1=ax
    
    sns.barplot(data=temp, x="F_Charity",y="Nlh",
                ax=ax2,n_boot=1000,seed=0,order=['Funded','Not Funded'])#,errorbar=('ci,95'))
    ax2.set_ylabel('')
    ax2.set_xlabel('Charity Funding')
    ax2.set_xticklabels(['',''])
    ax2.legend([],[], frameon=False)
    ax2.set_ylim(0,15)
    ax2.plot([0,1],[y0,y0],'-|',color='k',linewidth=1)

    ptemp=pg.pairwise_tukey(data=temp,dv='Nlh',between='F_Charity',)['p-tukey'][0]
    ax=ax2
    if ptemp<0.001:
        ax.plot(xmid,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid-dx,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid+dx,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid,y0,marker=(6, 2, 0),color='k',markersize=starsize)              
        ax.plot(xmid-dx,y0,marker=(6, 2, 0),color='k',markersize=starsize) 
        ax.plot(xmid+dx,y0,marker=(6, 2, 0),color='k',markersize=starsize)             
        
    elif ptemp<0.01:
        ax.plot(xmid-dx/2,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid+dx/2,y0,marker='o',color='w',markersize=squaresize)
        # ax[i_ax,1].plot(xmid,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid-dx/2,y0,marker=(6, 2, 0),color='k',markersize=starsize)
        ax.plot(xmid+dx/2,y0,marker=(6, 2, 0),color='k',markersize=starsize)
        
    elif ptemp<0.05:
        ax.plot(xmid,y0,marker='o',color='w',markersize=squaresize)
        ax.plot(xmid,y0,marker=(6, 2, 0),color='k',markersize=starsize)    

    ax2=ax
    
    print("\npost-hoc Tukey for F_Gov: \n")
    THSD=pg.pairwise_tukey(data=temp,dv='Nlh',between='F_Gov',)
    print(THSD.to_string())
    file.write("\n\npost-hoc Tukey for F_Gov: \n")
    file.write(THSD.to_string())
    
    print("\npost-hoc Tukey for F_Charity: \n")
    THSD=pg.pairwise_tukey(data=temp,dv='Nlh',between='F_Charity',)
    print(THSD.to_string())
    file.write("\n\npost-hoc Tukey for F_Charity: \n")
    file.write(THSD.to_string())
    
    
    l,r = ax0.get_xlim()
    b,t = ax0.get_ylim()
    
    ax0.text(0.03*(r-l)+l,t-0.03*(t-b),'A',fontsize=15,ha='left',va='top',
               backgroundcolor="#d6e1ff")
    
    l,r = ax1.get_xlim()
    b,t = ax1.get_ylim()
    ax1.text(0.05*(r-l)+l,t-0.05*(t-b),'B',fontsize=15,ha='left',va='top',
               backgroundcolor="#d6e1ff")
    
    l,r = ax2.get_xlim()
    b,t = ax2.get_ylim()
    ax2.text(0.05*(r-l)+l,t-0.05*(t-b),'C',fontsize=15,ha='left',va='top',
               backgroundcolor="#d6e1ff")
    
    
    # plt.tight_layout()
    file.close()


    fig.savefig(os.getcwd()+"/figures/"+figname+".png", dpi=600)
    fig.savefig(os.getcwd()+"/figures/"+figname+".svg")
    fig.savefig(os.getcwd()+"/figures/"+figname+".pdf")

# =============================================================================
# Figure 4_1
# =============================================================================
if savearr[2]==1:
    
    figname="F4_1_Fyrs";
    c_Nll="#4ba173";
    c_Nlh="#9900ff";
    c_Nhh="#980000ff";
    
    # fig, axd = plt.subplots(ncols=3, nrows=len(fund), sharex='col',
    #                         figsize=(10,20),sharey='col',
    #                         gridspec_kw={'width_ratios': [1, 1]})
    
    fig, axd = plt.subplots(ncols=3, nrows=len(fund), sharex='col',
                            figsize=(20,20),sharey='col',
                            gridspec_kw={'width_ratios': [1, 1, 0.6]})


    for idx,ax in enumerate(axd[:,0]):
        temp=df.loc[df[fund[idx]]==1,:]
        N=pd.DataFrame(data=None,index=yrs,columns=["Nlh","Nhh","Nll"])
        count=pd.DataFrame(data=None,index=yrs,columns=["count"])
        
        for i in yrs:
            adjacency, A=getadj(temp,i, Ucountry);
            N.loc[i,"Nlh"],N.loc[i,"Nhh"],N.loc[i,"Nll"]=Ncounter(A,lmicbins)
            count.loc[i,"count"] = (df["Publication year"]==i).sum()
        
        ax.plot(N.index,N.Nlh,'-o',color=c_Nlh,linewidth=1.5,label=r'$N_{l,h}$')
        ax.plot(N.index,N.Nll,'-o',color=c_Nll,linewidth=1.5,label=r'$N_{l,l}$')
        ax.plot(N.index,N.Nhh,'-o',color=c_Nhh,linewidth=1.5,label=r'$N_{h,h}$')
        # ax.legend(loc='upper center',fontsize=16,ncol=3,)
        ax.set_xticks(np.arange(np.min(yrs),np.max(yrs),1),minor=True,
                         color='k',alpha=0.2,linewidth=0.5);
        ax.set_xticks(np.arange(np.min(yrs),np.max(yrs)+5,5),minor=False,
                         color='k',alpha=0.5,linewidth=1);
        ax.set_ylabel('global coauthorship connections')
        ax.set_xlim(left=2009.5,right=2020.5);
        ax.grid(True, 'major',color='k',alpha=0.3,linewidth=1)
        ax.grid(True, 'minor',color='k',alpha=0.2,linewidth=0.5)
        ax.set_title(fund_true[idx])
        
        if idx==4:
            ax.set_xlabel('publication year')
            ax.legend(loc='upper center',fontsize=16,ncol=3,
                      bbox_to_anchor=(1.1,-0.2))
        
    letters=['A1','B1','C1','D1','E1'];
    
    for idx,ax in enumerate(axd[:,0]):
        l,r = ax.get_xlim()
        b,t = ax.get_ylim()
        ax.text(0.05*(r-l)+l,t-0.08*(t-b),letters[idx],fontsize=20,ha='left',va='top',
                   backgroundcolor="#d6e1ff")
    

    for idx,ax in enumerate(axd[:,1]):
        temp=df.loc[df[fund[idx]]==1,:]
        
        melt=pd.melt(temp,value_vars=["Nlh","Nhh","Nll"],id_vars=["Publication year"])
        melt=melt.loc[melt["Publication year"]<2021,:]# avoid 2021 cus incomplete
        
        sns.lineplot(data=melt,x="Publication year",y="value",hue="variable", 
                     ax=ax,legend=None,palette=[c_Nlh,c_Nhh,c_Nll],)
        ax.set_xticks(np.arange(np.min(yrs),np.max(yrs),1),minor=True,
                         color='k',alpha=0.2,linewidth=0.5);
        ax.set_xticks(np.arange(np.min(yrs),np.max(yrs)+5,5),minor=False,
                         color='k',alpha=0.5,linewidth=1);
        
        ax.set_ylabel('average individual connections')
        if idx==4:
            ax.set_xlabel('publication year')
        ax.grid(True, 'major',color='k',alpha=0.3,linewidth=1)
        ax.grid(True, 'minor',color='k',alpha=0.2,linewidth=0.5)
        ax.set_xlim(left=2009.5,right=2020.5);
        ax.set_title(fund_true[idx])
        
    letters=['A2','B2','C2','D2','E2'];
    for idx,ax in enumerate(axd[:,1]):
        l,r = ax.get_xlim()
        b,t = ax.get_ylim()
        ax.text(0.05*(r-l)+l,t-0.08*(t-b),letters[idx],fontsize=20,ha='left',va='top',
                   backgroundcolor="#d6e1ff")
        
    for idx,ax in enumerate(axd[:,2]):
        temp=df.loc[df[fund[idx]]==1,:]
        
        adjacency, A = fulladj(temp, Ucountry);
        G=getarcs(A, Ucountry, lmicbins);
        
        pos = nx.get_node_attributes(G,'pos')
        nx.draw_networkx_nodes(G,pos,node_size=5,ax=ax,node_color='k')
        weights=list(nx.get_edge_attributes(G,'weight').values())
        edgeidx=0;
        maxwt=17;
        
        edgeprintorder=np.argsort(weights)
        edgelist=list(G.edges());
        
        normweight=weights/np.double(maxwt);# just for visualization relative to max
        cmap = mpl.cm.get_cmap('magma_r')
        # camp=mpl.colors.LinearSegmentedColormap()
        # cmap=sns.color_palette("rocket")
        
        for idxval in edgeprintorder:
            edge=edgelist[idxval]
            
            source, target = edge
            if lmicbins.LMICever[edge[1]]==lmicbins.LMICever[edge[0]]: #from same setting
                rad = 0.7
            else:
                rad = 0.05
                
            arrowprops=dict(arrowstyle="-",
                            color=cmap(normweight[idxval]),
                            connectionstyle=f"arc3,rad={rad}",
                            linestyle= '-',
                            linewidth=1,# linewidth=np.log(weights[edgeidx])+0.5,
                            alpha=1,)
            ax.annotate("",
                    xy=pos[source],
                    xytext=pos[target],
                    arrowprops=arrowprops
                   )
        
        norm = mpl.colors.Normalize(vmin=0,vmax=maxwt)
        
        ax.set_title(fund_true[idx])
        ax.set_aspect('equal');
        ax.patch.set_facecolor('gray')

    letters=['A3','B3','C3','D3','E3'];
    for idx,ax in enumerate(axd[:,2]):
        l,r = ax.get_xlim()
        b,t = ax.get_ylim()
        ax.text(0.03*(r-l)+l,t-0.08*(t-b),letters[idx],fontsize=20,ha='left',va='top',)
        ax.text(0.01*(r-l)+l,t-0.99*(t-b),'HIC',fontsize=16,ha='left',va='bottom',)
        ax.text(0.99*(r-l)+l,t-0.99*(t-b),'LMIC',fontsize=16,ha='right',va='bottom',)
        
        if idx==4:
            cax = plt.axes([0.75, 0.08, 0.15, 0.02])
            mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
            cax.tick_params(labelsize=16)
            cax.set_title('total coauthorship connections')
        
   
    fig.savefig(os.getcwd()+"/figures/"+figname+".png", dpi=600)
    fig.savefig(os.getcwd()+"/figures/"+figname+".svg")
    fig.savefig(os.getcwd()+"/figures/"+figname+".pdf")


# =============================================================================
# Figure 5_2
# =============================================================================

if savearr[3]==1:
    
    figname="F5_2_Connectome";
    
    fig, axd = plt.subplots(ncols=4, nrows=4,figsize=(18,20),
                            gridspec_kw={'width_ratios': [1, 1, 1,0.1]})
    # axd=np.reshape(axd,[6,]);
    connectome_yrs=[2020,2019,2018,2017,2016,2015]
    
    all_yrs=np.arange(yrs.min(),yrs.max()+1)
    adjacency, A = getadj(df,2005, Ucountry);
    hm=pd.DataFrame(data=None,index=A.index,columns=None)
    
    for idx,cyr in enumerate(all_yrs):
        
        adjacency, A = getadj(df,cyr, Ucountry);
        
        hm[str(cyr)]=A.sum(axis=1)

    hm2=hm.copy();
    hm2['Countryname']=hm.index.to_list()
    hm2['LMICever']=np.nan
    hm2['LMICidx']=np.nan
    for idx in range(len(hm2)):
        cntry=hm2.index[idx];
        hm2.loc[cntry,'LMICever']=lmicbins.loc[lmicbins.Country.to_list().index(cntry),'LMICever']
        hm2.loc[cntry,'LMICidx']=lmicbins.Country.to_list().index(cntry);
    
    # hm2['LMICever']=lmicbins.LMICever.to_list();
    # hm2['Countryname']=lmicbins.Country.to_list();
    hm2['sum']=hm.sum(axis=1)
    hm3=hm2.copy();
    hm3=hm3.sort_values(by=['LMICever','sum'],ascending=[True,False])
    
    hic=hm3.loc[hm3.LMICever==0,:].head(10)
    lmic=hm3.loc[hm3.LMICever==1,:].head(10)
    
    
    axcon=np.reshape(axd[:2,:3],[6,]);

    letters=['A1','A2','A3','A4','A5','A6']    

    for idx,cyr in enumerate(connectome_yrs):
        
        ax=axcon[5-idx]
        ax.set_title(cyr,fontsize=20)
        # cyr=2015
        
                
        adjacency, A = getadj(df,cyr, Ucountry);
        G=getarcs(A, Ucountry, lmicbins);
        
        pos = nx.get_node_attributes(G,'pos')
        # nx.draw_networkx_nodes(G,pos,node_size=5,ax=ax,node_color='k')
        weights=list(nx.get_edge_attributes(G,'weight').values())
        edgeidx=0;
        maxwt=10;
        
        edgeprintorder=np.argsort(weights)
        edgelist=list(G.edges());
        
        normweight=weights/np.double(maxwt);# just for visualization relative to max
        cmap = mpl.cm.get_cmap('magma_r')
        # camp=mpl.colors.LinearSegmentedColormap()
        # cmap=sns.color_palette("rocket")
        
        for idxval in edgeprintorder:
            edge=edgelist[idxval]
            
            source, target = edge
            if lmicbins.LMICever[edge[1]]==lmicbins.LMICever[edge[0]]: #from same setting
                rad = 0.7
            else:
                rad = 0.05
                
            arrowprops=dict(arrowstyle="-",
                            color=cmap(normweight[idxval]),
                            connectionstyle=f"arc3,rad={rad}",
                            linestyle= '-',
                            linewidth=1,# linewidth=np.log(weights[edgeidx])+0.5,
                            alpha=1,)
            ax.annotate("",
                    xy=pos[source],
                    xytext=pos[target],
                    arrowprops=arrowprops
                   )
        
        nx.draw_networkx_nodes(G,pos,node_size=5,ax=ax,node_color='k')
        nx.draw_networkx_nodes(G,pos,node_size=50,ax=ax,node_color='b',
                                nodelist=hic.LMICidx.to_list(),
                                node_shape='*')
        nx.draw_networkx_nodes(G,pos,node_size=50,ax=ax,node_color='b',
                                nodelist=lmic.LMICidx.to_list(),
                                node_shape='*')
        
        # for val in hic.LMICidx:
        #     posval=pos[val];
        #     ax.arrow(x=-1.1,y=posval[1],dx=posval[0]+1.1,dy=0,
        #              length_includes_head=True,color='k',linewidth=2)
            
        # for val in lmic.LMICidx:
        #     posval=pos[val];
        #     ax.arrow(x=1.1,y=posval[1],dx=posval[0]-1.1,dy=0,
        #              length_includes_head=True,color='k',linewidth=2)
        
        

        norm = mpl.colors.Normalize(vmin=0,vmax=maxwt)
        
        # ax.set_title(cyr)
        ax.set_aspect('equal');
        # ax.patch.set_facecolor('gray')
        ax.axis('off')
        
        l,r = ax.get_xlim()
        b,t = ax.get_ylim()
        ax.text(0.05*(r-l)+l,t-0.99*(t-b),'HIC',fontsize=16,ha='left',va='bottom',)
        ax.text(0.95*(r-l)+l,t-0.99*(t-b),'LMIC',fontsize=16,ha='right',va='bottom',)
        ax.text(0.0*(r-l)+l,t-0.0*(t-b),letters[5-idx],fontsize=25,ha='left',va='top',)
        
        
        
        
        if idx==0:
            # plt.subplots_adjust(bottom=0.1, right=0.9, top=0.9)
            # cax = plt.axes([0.93, 0.1, 0.03, 0.8])
            gs = axd[0,3].get_gridspec()
            # remove the underlying axes
            for axtemp in axd[:2,3]:
                axtemp.remove()
            cax = fig.add_subplot(gs[:2,3])
            mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,
                                       orientation='vertical',)
            cax.tick_params(labelsize=16)
            cax.set_ylabel('yearly coauthorship connections',fontsize=18)
            
            
    
    fig.patch.set_facecolor('gray')    
    
    gs = axd[0,2].get_gridspec()
    # remove the underlying axes
    for axtemp in axd[2,:3]:
        axtemp.remove()
    row3 = fig.add_subplot(gs[2,:3])
    
    gs = axd[0,3].get_gridspec()
    # remove the underlying axes
    for axtemp in axd[3,:3]:
        axtemp.remove()
    row4 = fig.add_subplot(gs[3,:3])
    
    cmap = mpl.cm.get_cmap('magma_r')
    highest=hm3[hm.columns].max().max()
    norm = mpl.colors.Normalize(vmin=0,vmax=highest)
    
    gs = axd[2,3].get_gridspec()
    # remove the underlying axes
    for axtemp in axd[2:,3]:
        axtemp.remove()
    cax = fig.add_subplot(gs[2:,3])
    mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,
                               orientation='vertical',)
    cax.tick_params(labelsize=16)
    cax.set_ylabel('total yearly coauthorship connections',fontsize=18)
    
    sns.heatmap(hic[hm.columns],yticklabels=hic.index,ax=row3,square=False,
                cbar=False,cmap=cmap,vmin=0,vmax=highest);
    row3.tick_params(axis='x', labelsize=15,rotation=45)
    row3.tick_params(axis='y', labelsize=15)
    
    sns.heatmap(lmic[hm.columns],yticklabels=lmic.index,ax=row4,square=False,
                cbar=False,cmap=cmap,vmin=0,vmax=highest,);
    row4.tick_params(axis='y', labelsize=15)
    row4.tick_params(axis='x', labelsize=15,rotation=45)
    # row4.set_xticklabels(fontsize=12)
    
    # row3.set_title('B')
    # row3.text(-0.05*(r-l)+l,t-0.08*(t-b),'B',fontsize=20,ha='left',va='top',)
    fig.text(0.13, 0.47, 'B', fontsize=25)
    fig.text(0.13, 0.27, 'C', fontsize=25)
    
    fig.savefig(os.getcwd()+"/figures/"+figname+".png", dpi=600)
    fig.savefig(os.getcwd()+"/figures/"+figname+".svg")    
    fig.savefig(os.getcwd()+"/figures/"+figname+".pdf")    
    
    
