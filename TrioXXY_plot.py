#!/usr/bin/python
# -*- coding: UTF-8 -*-

##=====================================
## XXY plot
##=====================================

import os,sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
from adjustText import adjust_text

def depth_ratio_plot(target_chr, df_melt, Depth_ratio_box):

    df_melt_target = df_melt[df_melt['chr']==target_chr]
    df_melt_control = df_melt[~df_melt.chr.isin([target_chr, 'chrX'])]

    fig, ax = plt.subplots(figsize=(4,4))
    sns.boxplot(x="sample", y="depth", data=df_melt_control,
        palette="Set2", zorder=1)
    sns.swarmplot(x="sample", y="depth", data=df_melt_control,
        color="Black", size=3, zorder=2)
    plt.scatter(x="sample", y="depth", data=df_melt_target,
        s=10, color='red', zorder=3)

    plt.xlabel("sample", fontdict={'size': 10})
    plt.ylabel("read depth", fontdict={'size': 10})
    samples = df_melt_control['sample'].unique().tolist()
    xtick_indexs = range(0,len(samples))
    xtick_labels = [s.replace('_depth','') for s in samples]
    plt.xticks(xtick_indexs, xtick_labels, fontsize=8)
    plt.yticks(fontsize=8)

    plt.tight_layout()
    plt.savefig(Depth_ratio_box, dpi=500)
    plt.close(fig)

def depth_chr_plot(df_melt_depth, Depth_sta_per_chr_line):
    fig, ax = plt.subplots(figsize=(12,4))
    ax.set_xlabel('Chr', fontsize=12)
    ax.set_ylabel('Marker mean depth', fontsize=12)
    ax = sns.lineplot(x='chr', y='depth', data = df_melt_depth,
        hue='sample', linewidth = 1.5)

    plt.tight_layout()
    plt.savefig(Depth_sta_per_chr_line, dpi=500)
    plt.close(fig)

def marker_num_chr_plot(df_melt_marker, Marker_sta_per_chr_bar):
    fig, ax = plt.subplots(figsize=(12,4))
    ax.set_xlabel('Chr', fontsize=12)
    ax.set_ylabel('Marker number', fontsize=12)
    ax = sns.barplot(x='chr', y='marker_number', data = df_melt_marker,
        hue='sample')

    plt.tight_layout()
    plt.savefig(Marker_sta_per_chr_bar, dpi=500)
    plt.close(fig)

def Muta_num_plot(target_chr, df, muta_count_bar):
    fig,ax = plt.subplots(1,1,figsize=(4,4))
    # muta nuber
    ax = sns.barplot(x="sample", y=target_chr+'_muta_num',
        palette="Set2", hue="sex", data=df)
    plt.legend(fontsize=8)
    ax.set_title(target_chr+': Muta number')
    ax.set_xlabel("sample", fontsize=10)
    ax.set_ylabel("Muta number", fontsize=10)
    plt.xticks(rotation=30, fontsize=8)
    plt.yticks(fontsize=8)

    plt.tight_layout()
    plt.savefig(muta_count_bar, dpi=500)
    plt.close(fig)

def Muta_het_ratio_plot(target_chr, df, het_ratio_plot):
    fig,ax = plt.subplots(1,1,figsize=(4,4))
    # muta nuber
    ax = sns.scatterplot(x="sex", y=target_chr+'_Het_Hom_ratio', data=df,
        hue="sample", alpha=0.8, palette="Set2")
    plt.legend(fontsize=8)
    ax.set_title(target_chr+': Het/Hom mutation ratio')
    ax.set_xlabel("Sex From Ped", fontsize=10)
    ax.set_ylabel("Het/Hom ratio", fontsize=10)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)

    plt.tight_layout()
    plt.savefig(het_ratio_plot, dpi=500)
    plt.close(fig)

def Chr_karyotyping_plot(muta_count_file, Depth_ratio_file, karyotyping_plot_out):
    df_muta_count = pd.read_csv(muta_count_file, sep='\t')
    df_Depth_ratio = pd.read_csv(Depth_ratio_file, sep='\t')
    df_kary = pd.merge(df_muta_count, df_Depth_ratio)
    df_kary_female = df_kary[df_kary['sex']=='female']
    df_kary_male = df_kary[df_kary['sex']=='male']
    female_num = len(df_kary_female['sample'].tolist())

    fig, ax = plt.subplots(figsize=(4, 4))
    back_color_male = 'skyblue'
    back_color_female = 'pink'
    back_color_KSIII = 'darkgray'
    x_sta = "median_depth_ratio"
    y_sta = "chrX_Het_Hom_ratio"
    samples = df_kary['sample']

    # male region
    xy1 = (0,0)
    width1 = 0.85
    hight1 = 0.3
    # female region
    xy2 = (0.85,0.3)
    # if female_num != 0:
    width2 = max(df_kary[x_sta].tolist()) - 0.8
    hight2 = max(df_kary[y_sta].tolist()) - 0.25
    # else:
    #     width2_1 = 0.5
    #     hight2_1 = 0.5
    #     width2_2 = max(df_kary[x_sta].tolist()) - 0.8
    #     hight2_2 = max(df_kary[y_sta].tolist()) - 0.25
    #     width2 = max([width2_1, width2_2])
    #     hight2 = max([hight2_1, hight2_2])
    # KS-III
    xy3 = (0.85,0)
    ax.add_patch(patches.Rectangle(xy1,width1,hight1, color =back_color_male, alpha=0.9, zorder=1)) # male
    plt.text(0, 0, 'Male', size = 6)
    ax.add_patch(patches.Rectangle(xy2,width2,hight2, color =back_color_female, alpha=0.9, zorder=1)) # female
    plt.text(0.85, 0.3, 'Female/KS-I/KS-II', size = 6)
    ax.add_patch(patches.Rectangle(xy3,width2,hight1, color =back_color_KSIII, alpha=0.9, zorder=1)) # KS-III
    plt.text(0.85, 0, 'KS-III', size = 6)

    sns.scatterplot(x=x_sta, y=y_sta, data=df_kary,
        s=10, hue="sex", alpha=0.8, palette="Set1", zorder=2)

    for i in range(len(samples)):
        texts = [plt.text(df_kary[x_sta][i], df_kary[y_sta][i], \
            df_kary['sample'][i], \
            horizontalalignment='left', fontsize = 3, zorder=2)]
        adjust_text(texts,)

    plt.legend(fontsize=8)
    ax.set_title('chrX karyotyping')
    ax.set_xlabel("depth_ratio", fontsize=10)
    ax.set_ylabel("chrX_Het_Hom_ratio", fontsize=10)

    plt.tight_layout()
    plt.savefig(karyotyping_plot_out, dpi=500)
    plt.close(fig)

def sample_Y_muta_venn(h_sample_Y_muta, sample1, sample2, fig_out_path):
    Y_muta1 = h_sample_Y_muta[sample1]
    Y_muta2 = h_sample_Y_muta[sample2]
    fig_out = os.path.join(fig_out_path, sample1+'_'+sample2+'_Y_muta_compare.png')
    fig, ax1 = plt.subplots(figsize=(5,4))
    venn2(subsets = [set(Y_muta1),set(Y_muta2)],
        set_labels = [sample1,sample2],
        set_colors = ('r','b')
        )

    plt.tight_layout()
    plt.savefig(fig_out, dpi=500)
    plt.close(fig)

def Chr_origin_plot(parent, sta_result, plot_out):
    df = pd.read_csv(sta_result, sep='\t')

    fig, ax = plt.subplots(figsize=(4, 4))
    color_1 = 'skyblue'
    color_2 = 'pink'
    color_3 = 'wheat'
    x_sta = parent+"_info_het_ratio"
    y_sta = parent+"_info_hom_ratio"
    x_num = parent+"_info_het_num"
    y_num = parent+"_info_hom_num"
    all_num = parent+"_unique_num"
    samples = df['sample']

    if parent == 'father':
            xy1 = (0,0)
            width1 = 0.25
            hight1 = 0.25
            xy2 = (0.75,0)
            width2 = 0.25
            hight2 = 0.25
            ax.add_patch(patches.Rectangle(xy1,width1,hight1, color =color_1, alpha=0.3, zorder=1))
            ax.add_patch(patches.Rectangle(xy2,width2,hight2, color =color_2, alpha=0.3, zorder=1))

    if parent == 'mother':
        xy1 = (0,0.25)
        width1 = 0.25
        hight1 = 0.5
        xy2 = (0.25,0)
        width2 = 0.5
        hight2 = 0.25
        xy3 = (0.75,0)
        width3 = 0.25
        hight3 = 0.25
        ax.add_patch(patches.Rectangle(xy1,width1,hight1, color =color_1, alpha=0.3, zorder=1))
        ax.add_patch(patches.Rectangle(xy2,width2,hight2, color =color_2, alpha=0.3, zorder=1))
        ax.add_patch(patches.Rectangle(xy3,width3,hight3, color =color_3, alpha=0.3, zorder=1))

    sns.scatterplot(x=x_sta, y=y_sta, data=df,
        s=7, hue="sample", alpha=0.8, palette="Set1", zorder=2)

    for i in range(len(samples)):
        texts = [plt.text(df[x_sta][i], df[y_sta][i], \
            'Het: '+str(df[x_num][i])+'\n'+'Hom: '+str(df[y_num][i])+'\n'+'Ref: '+str(df[all_num][i]-df[x_num][i]-df[y_num][i]), \
            horizontalalignment='left', fontsize = 3, zorder=2)]
        adjust_text(texts,)

    plt.legend(fontsize=8)
    ax.set_title(parent + '-specific mutation')
    ax.set_xlabel("Het", fontsize=10)
    ax.set_ylabel("Hom", fontsize=10)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)

    plt.tight_layout()
    plt.savefig(plot_out, dpi=500)
    plt.close(fig)
