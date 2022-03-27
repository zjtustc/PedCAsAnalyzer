#!/usr/bin/python
# -*- coding: UTF-8 -*-

##=====================================
## XXY analysis
##=====================================

import os,sys
import pandas as pd
from script import TrioXXY_plot
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-task", action='store', dest='task_ID',
        help="taskID.")
parser.add_argument("-ped", action='store', dest='ped',
        help="pedigree ped file.")
parser.add_argument("-chr", nargs='+', dest='A_chr', default='NA',
        help="Autosome chr analysis, for example: chr21.")
parser.add_argument("-plot", action='store', dest='plot', default='Y',
        help="Y----plot, N----don't plot")
parser.add_argument("-MAF", action='store', dest='MAF', default='0.02',
        help="parent sepcific marker MAF")

paras = parser.parse_args()

################################### Step1: preprocess ###################################

# Step1.1: Preprocess: Sample check
def sampleCheck():
    sampleCheckResult = 'not pass'
    Step1_Note = []
    # samples from ped file
    h_ped = {}
    father = []
    mother = []
    with open(ped_file,'r')as fi:
        for line in fi:
            line_list = line.strip('\r\n').split('\t')
            sample = line_list[1]
            sample_sex = line_list[4]
            if str(sample_sex) == '1':
                sample_sex_new = 'male'
            elif str(sample_sex) == '2':
                sample_sex_new = 'female'
            else:
                sample_sex_new = 'unknow'
            if line_list[2] not in father and str(line_list[2]) != '0':
                father.append(str(line_list[2]))
            if line_list[3] not in mother and str(line_list[3]) != '0':
                mother.append(str(line_list[3]))
            h_ped[sample] = line_list[0:4]+[sample_sex_new]+line_list[5:]

    ped_samples = list(h_ped.keys())
    # samples from vcf files
    vcf_files = [file for file in os.listdir(taskDir) if file.split('.')[-1] == 'vcf']
    vcf_samples = [file.split('.')[0] for file in vcf_files]
    # compare
    ped_unique_samples = [sample for sample in ped_samples if sample not in vcf_samples]
    vcf_unique_samples = [sample for sample in vcf_samples if sample not in ped_samples]

    parents_list = father + mother
    patient_list = [sample for sample in h_ped.keys() if str(h_ped[sample][5])=='2' and sample not in parents_list]
    control_list = [sample for sample in h_ped.keys() if str(h_ped[sample][5])=='1' and sample not in parents_list]
    child_list = sorted(patient_list+control_list)

    if len(ped_unique_samples) != 0:
        Step1_Note.append('%s%s'%('## ',', '.join(ped_unique_samples)+' in ped files but vcf not found!'))
    elif len(vcf_unique_samples) != 0:
        Step1_Note.append('%s%s'%('## ',', '.join(vcf_unique_samples)+' vcf found but not in ped files!'))
    elif len(father) > 1 or len(mother) > 1:
        Step1_Note.append('%s'%('## one trio for one analysis'))
    else:
        sampleCheckResult = 'pass'
    return sampleCheckResult, Step1_Note, father, mother, child_list, h_ped

# Step1.2: Preprocess: Process & Merge VCF files
def get_PAR():
    h_par = {'X':[[60001,2699520],[154931044,155260560]], \
    'Y':[[10001,2649520],[59034050,59363566]]}
    return h_par

def readInMuta(fileIn):
    isfindex = 1
    h_muta_sample = {}
    with open(fileIn, 'r') as fi:
        for line in fi:
            if line.startswith('#'):
                continue
            else:
                line_list = line.strip('\r\n').split('\t')
                chrom = 'chr' + line_list[0].replace('chr', '')
                start = line_list[1]
                ref = line_list[3]
                alt = line_list[4]
                # qual = line_list[5]
                fmat = line_list[8]
                fmat_s = line_list[9].split(':')
                ## format index only need to extract once
                if isfindex:
                    fmat_index = getFormatIndex(fmat)
                    isfindex = 0
                GT = fmat_s[fmat_index['GT']]
                if 'GQ' in fmat_index:
                    GQ = fmat_s[fmat_index['GQ']]
                else:
                    GQ = '90'
                if 'AD' in fmat_index:
                    AD = fmat_s[fmat_index['AD']]
                elif 'RO' in fmat_index and 'AO' in fmat_index:
                    RO = fmat_s[fmat_index['RO']]
                    AO = fmat_s[fmat_index['AO']]
                    AD = '%s,%s' %(RO, AO)
                else:
                    AD = 'NA,NA'
                Ref_DP = AD.split(',')[0]
                Alt_DP = AD.split(',')[1]
                if ',' not in alt:
                    site = '\t'.join([chrom, start, start, ref, alt])
                    info = GT+':'+GQ+':'+Ref_DP+':'+Alt_DP
                    if chrom == 'X':
                        inPar = 'False'
                        for par in h_par['X']:
                            if par[0] <= int(start) and par[1] >= int(start):
                                inPar = 'True'
                            else:
                                continue
                            if inPar == 'False':
                                h_muta_sample[site] = info
                    if chrom == 'Y':
                        inPar = 'False'
                        for par in h_par['Y']:
                            if par[0] <= int(start) and par[1] >= int(start):
                                inPar = 'True'
                            else:
                                continue
                            if inPar == 'False':
                                h_muta_sample[site] = info
                    else:
                        h_muta_sample[site] = info
    return h_muta_sample

def getFormatIndex(fmat):
    h = {}
    fmat = fmat.split(':')
    for i in fmat:
        h[i] = fmat.index(i)
    return h

def merge_ped_vcf(all_sample,all_sample_muta,fileOut):
    all_muta = []
    for h_muta in all_sample_muta:
        all_muta += list(h_muta.keys())
    all_muta = list(set(all_muta))

    miss_geno = ':'.join(['0/0', '0', '0', '0'])

    with open(fileOut,'w')as fo:
        header = '\t'.join(['Chr','Start','End','Ref','Alt'] + all_sample)
        fo.write('%s\n'%(header))

        for muta in all_muta:
            L_muta = [muta]
            for sample in all_sample:
                sample_info = all_sample_muta[all_sample.index(sample)].get(muta,miss_geno)
                L_muta.append(sample_info)
            muta_line = '\t'.join(L_muta)
            fo.write('%s\n'%(muta_line))

def VCF_process_main():
    all_sample_muta = []
    # Process and Merge VCF files !
    for sample in sample_list:
        vcfIn = os.path.join(taskDir, sample+'.vcf')

        h_muta_sample = readInMuta(vcfIn)
        all_sample_muta.append(h_muta_sample)

    merge_ped_vcf(sample_list,all_sample_muta,muta_file)

################################### Step2: Genotype quality ###################################

def quality_site(GT_info):
    GT_info_list = GT_info.split(':')
    GT = GT_info_list[0]
    GQ = GT_info_list[1]
    Ref_DP = GT_info_list[2]
    Alt_DP = GT_info_list[3]

    Qual_marker = 'PASS'

    if GT == '0/0':
        Qual_marker = 'PASS'
    else:
        DP = int(Ref_DP) + int(Alt_DP)
        if DP<15 or int(GQ)<65:
            Qual_marker = 'FAIL'
        if DP != 0:
            alt_percent = int(float('%.2f' % (int(Alt_DP)/(DP)))*100)
            if GT == '0/1' and alt_percent not in range(33,66):
                Qual_marker = 'FAIL'
            elif GT == '1/1' and alt_percent < 75:
                Qual_marker = 'FAIL'
    return Qual_marker

def quality_all_sample(line_list):
    L_sample_seq_marker = []
    sample_GT_infos = line_list[5:]
    for sample_GT_info in sample_GT_infos:
        sample_seq_marker = quality_site(sample_GT_info)
        L_sample_seq_marker.append(sample_seq_marker)
    if 'FAIL' in L_sample_seq_marker:
        return 'FAIL'
    else:
        return 'PASS'

def Genotype_quality_main():
    with open(muta_file,'r')as fi:
        with open(muta_quality_file,'w')as fo:
            header = fi.readline().strip('\r\n')
            fo.write('%s\n'%(header))
            for line in fi:
                line_list = line.strip('\r\n').split('\t')
                Seq_quality_result = quality_all_sample(line_list)
                if Seq_quality_result == 'PASS':
                    fo.write('%s'%(line))
                else:
                    continue

################################### Step3: Depth analysis ###################################
def get_chr_df(_chr):
    df_chr = muta_df[muta_df['Chr']==_chr]
    return df_chr

def depth_count(fileIn):
    h_depth = {}
    with open(fileIn,'r')as fi:
        header_list = fi.readline().strip('\r\n').split('\t')
        for line in fi:
            line_list = line.strip('\r\n').split('\t')
            _chr = line_list[0]
            if _chr in chr_list:
                ped_infos = line_list[5:]
                ped_gts = [i.split(':')[0] for i in ped_infos]
                ped_depths = [int(i.split(':')[2]) + int(i.split(':')[3]) for i in ped_infos]
                if _chr not in h_depth.keys():
                    h_depth[_chr] = [ped_depths]
                else:
                    h_depth[_chr].append(ped_depths)
            else:
                pass

    return h_depth

def depth_sta(target_chr, h_depth):
    h_sta = {}
    for _chr in h_depth.keys():
        maker_nums = []
        mean_depths = []
        for i in range(len(h_ped.keys())):
            sample_depths = [depths[i] for depths in h_depth[_chr]]
            sample_depths = [i for i in sample_depths if i != 0]
            sample_maker_num = len(sample_depths)
            sample_mean_depth = np.mean(sample_depths)
            maker_nums.append(str(sample_maker_num)) # sample muta number
            mean_depths.append(str(sample_mean_depth).split('.')[0]) # sample muta mean depth
        if _chr == 'chrX':
            h_sta[_chr] = ['chrX'] + maker_nums + mean_depths
        else:
            h_sta[_chr] = ['Auto'] + maker_nums + mean_depths

    header2_1 = [s+'_maker_num' for s in sample_list]
    header2_2 = [s+'_depth' for s in sample_list]
    Depth_matrix = os.path.join(X_outputDir, 'Depth_matrix.txt')
    with open(Depth_matrix,'w')as fo:
        header = '\t'.join(['chr','chr_class']+header2_1+header2_2)
        fo.write('%s\n'%(header))
        for _chr in chr_list:
            fo.write('%s\n'%(_chr + '\t' +'\t'.join(h_sta[_chr])))

    df = pd.read_csv(Depth_matrix,sep='\t')

    df_target = df[df['chr']==target_chr]
    df_control = df[~df.chr.isin([target_chr, 'chrX'])]

    Depth_ratio = os.path.join(task_outputDir,target_chr, target_chr+'_Depth_ratio.txt')
    with open(Depth_ratio,'w')as fo2:
        header = '\t'.join(['sample','sex','min_depth_ratio','median_depth_ratio','max_depth_ratio'])
        fo2.write('%s\n'%(header))
        for sample in sample_list:
            sample_sex = h_ped[sample][4]
            max_control_depth = df_control[sample+'_depth'].max()
            median_control_depth = df_control[sample+'_depth'].median()
            min_control_depth = df_control[sample+'_depth'].min()
            targetChr_depth = df_target[sample+'_depth'].tolist()[0]
            note_i_list = '\t'.join(['%.2f'%(targetChr_depth/max_control_depth),'%.2f'%(targetChr_depth/median_control_depth),'%.2f'%(targetChr_depth/min_control_depth)])
            result = sample + '\t' + sample_sex + '\t' + note_i_list
            fo2.write('%s\n'%(result))

def Muta_count(target_chr, df_targetChr, muta_count):
    with open(muta_count,'w')as fo:
        header_list = ['sample','sex',target_chr+'_muta_num',target_chr+'_Het_Hom_ratio']
        header = '\t'.join(header_list)
        fo.write('%s\n'%(header))
        for sample in sample_list:
            sample_sex = h_ped[sample][4]
            sample_het_num = len(df_targetChr[df_targetChr[sample].str.contains('0/1')])
            sample_hom_num = len(df_targetChr[df_targetChr[sample].str.contains('1/1')])
            muta_num = sample_het_num + sample_hom_num
            if sample_hom_num != 0:
                Het_Hom_ratio = round(sample_het_num/sample_hom_num,3)
            else:
                Het_Hom_ratio = 0 # 这里不严谨
            final_result = [sample,sample_sex,muta_num,Het_Hom_ratio]
            final_result = '\t'.join([str(i) for i in final_result])
            fo.write('%s\n'%(final_result))

def Depth_analysis_main(target_chr, df_targetChr, resultDir):
    h_depth = depth_count(muta_quality_file)
    depth_sta(target_chr, h_depth)
    Depth_matrix = os.path.join(X_outputDir,'Depth_matrix.txt')

    if paras.plot == 'Y':
        header2_1 = [s+'_maker_num' for s in sample_list]
        header2_2 = [s+'_depth' for s in sample_list]
        df = pd.read_csv(Depth_matrix,sep='\t')
        df_melt_depth = df.melt(id_vars=["chr","chr_class"]+header2_1,var_name='sample',value_name='depth')
        df_melt_marker = df.melt(id_vars=["chr","chr_class"]+header2_2,var_name='sample',value_name='marker_number')
        Depth_ratio_box = os.path.join(resultDir, target_chr+'_Depth_ratio_box.png')
        Depth_sta_per_chr_line = os.path.join(resultDir, 'Depth_sta_per_chr_line.png')
        Marker_sta_per_chr_bar = os.path.join(resultDir, 'Marker_sta_per_chr_bar.png')
        TrioXXY_plot.depth_ratio_plot(target_chr, df_melt_depth, Depth_ratio_box)
        TrioXXY_plot.depth_chr_plot(df_melt_depth, Depth_sta_per_chr_line)
        TrioXXY_plot.marker_num_chr_plot(df_melt_marker, Marker_sta_per_chr_bar)

    muta_count = os.path.join(resultDir,target_chr+'_muta_count.txt')
    Muta_count(target_chr, df_targetChr, muta_count)

    if paras.plot == 'Y':
        muta_count_df =  pd.read_csv(muta_count,sep='\t')
        muta_count_bar = os.path.join(resultDir, target_chr+'_muta_count_bar.png')
        TrioXXY_plot.Muta_num_plot(target_chr, muta_count_df, muta_count_bar)
        het_ratio_plot = os.path.join(resultDir, target_chr+'_het_ratio.png')
        TrioXXY_plot.Muta_het_ratio_plot(target_chr, muta_count_df, het_ratio_plot)

    Depth_ratio = os.path.join(resultDir, target_chr+'_Depth_ratio.txt')
    if paras.plot == 'Y' and target_chr == 'chrX':
        karyotyping_plot_out  = os.path.join(resultDir, 'chrX_karyotyping.png')
        TrioXXY_plot.Chr_karyotyping_plot(muta_count, Depth_ratio, karyotyping_plot_out)

################################### Step4: Y mutation analysis ###################################

def get_sample_Y_muta(sample):
    L_sample_Y_muta = []
    df_sample_Y_muta = df_Y[df_Y[sample].str.contains('1/1')]
    for index,row in df_sample_Y_muta.iterrows():
        Y_muta = '_'.join([str(row['Start']), str(row['Ref']), str(row['Alt'])])
        L_sample_Y_muta.append(Y_muta)
    return L_sample_Y_muta

def Y_hom_muta_number_sta():
    h_sample_Y_muta = {}
    for sample in sample_list:
        L_sample_Y_muta = get_sample_Y_muta(sample)
        h_sample_Y_muta[sample] = L_sample_Y_muta

    if paras.plot == 'Y':
        i = 0
        while (i < len(sample_list)):
            sample1 = sample_list[i]
            for k in range(i+1, len(sample_list)):
                sample2 = sample_list[k]
                TrioXXY_plot.sample_Y_muta_venn(h_sample_Y_muta, sample1, sample2, X_outputDir)
            i += 1

################################### Step5: X origin analysis ###################################
def Muta_anno(df_target_chr, target_chr, out_path):
    h_Low_MAF_site_1kg = {}
    with open(MAF_site_1kg_file, 'r')as fi:
        for line in fi:
            line_list = line.strip('\r\n').split('\t')
            site_maf = line_list[-1]
            if float(site_maf) <= MAF:
                muta_site = '\t'.join(line_list[0:4])
                h_Low_MAF_site_1kg[muta_site] = site_maf

    muta_file = os.path.join(out_path, target_chr + '_muta.txt')
    df_target_chr.to_csv(muta_file, sep='\t', index=False)

    muta_fre_anno_file = os.path.join(out_path, target_chr + '_' + str(MAF) + '_muta_fre_anno.txt')
    with open(muta_file,'r')as fi:
        with open(muta_fre_anno_file,'w')as fo:
            header = fi.readline().strip('\r\n')
            fo.write('%s\t%s\n'%(header, '1kG_maf'))
            for line in fi:
                line_list = line.strip('\r\n').split('\t')
                muta_site = '\t'.join([line_list[0].lstrip('chr'), line_list[1], line_list[3], line_list[4]])
                if muta_site in h_Low_MAF_site_1kg.keys():
                    _1kG_maf = h_Low_MAF_site_1kg[muta_site]
                    fo.write('%s\t%s\n'%(line.strip('\r\n'), _1kG_maf))

def Chr_origin_analysis(df_target_chr, target_chr, father_unique_df, mother_unique_het_df, sta_result4):
    mother_unique_het_num = len(mother_unique_het_df)
    father_unique_num = len(father_unique_df)

    if len(father) == 1 and len(mother) == 1:
        father_muta_num = len(df_target_chr[~df_target_chr[father[0]].str.contains('0/0')])
        mother_muta_num = len(df_target_chr[~df_target_chr[mother[0]].str.contains('0/0')])
    elif len(father) == 1 and len(mother) == 0:
        father_muta_num = len(df_target_chr[~df_target_chr[father[0]].str.contains('0/0')])
        mother_muta_num = 'NA'
    elif len(mother) == 1 and len(father) == 0:
        mother_muta_num = len(df_target_chr[~df_target_chr[mother[0]].str.contains('0/0')])
        father_muta_num = 'NA'

    if mother_unique_het_num == 0:
        mother_info_het_num = 'NA'
        mother_info_het_ratio = 'NA'
        mother_info_hom_num = 'NA'
        mother_info_hom_ratio = 'NA'
    if father_unique_num == 0:
        father_info_het_num = 'NA'
        father_info_het_ratio = 'NA'
        father_info_hom_num = 'NA'
        father_info_hom_ratio = 'NA'

    with open(sta_result4,'w')as fo:
        header_list = ['sample','sample_sex','father_muta_num','mother_muta_num','sample_muta_num','Het_Hom_ratio', \
        'mother_unique_num','mother_info_het_num','mother_info_het_ratio','mother_info_hom_num','mother_info_hom_ratio', \
        'father_unique_num','father_info_het_num','father_info_het_ratio','father_info_hom_num','father_info_hom_ratio']
        header = '\t'.join(header_list)
        fo.write('%s\n'%(header))
        for sample in child_list:
            sample_muta_num = len(df_target_chr[~df_target_chr[sample].str.contains('0/0')])
            sample_sex = h_ped[sample][4]
            sample_het_num = len(df_target_chr[df_target_chr[sample].str.contains('0/1')])
            sample_hom_num = len(df_target_chr[df_target_chr[sample].str.contains('1/1')])
            if sample_hom_num != 0:
                Het_Hom_ratio = round(sample_het_num/sample_hom_num,3)
            else:
                Het_Hom_ratio = 'NA'

            if mother_unique_het_num != 0:
                sample_mother_inher_het_df = mother_unique_het_df[mother_unique_het_df[sample].str.contains('0/1')]
                mother_info_het_num = len(sample_mother_inher_het_df)
                mother_info_het_ratio = round(mother_info_het_num/mother_unique_het_num,3)
                mother_info_hom_num = len(mother_unique_het_df[mother_unique_het_df[sample].str.contains('1/1')])
                mother_info_hom_ratio = round(mother_info_hom_num/mother_unique_het_num,3)
            if father_unique_num != 0:
                sample_father_inher_het_df = father_unique_df[father_unique_df[sample].str.contains('0/1')]
                father_info_het_num = len(sample_father_inher_het_df)
                father_info_het_ratio = round(father_info_het_num/father_unique_num,3)
                father_info_hom_num = len(father_unique_df[father_unique_df[sample].str.contains('1/1')])
                father_info_hom_ratio = round(father_info_hom_num/father_unique_num,3)

            final_result = [sample,sample_sex,father_muta_num,mother_muta_num,sample_muta_num,Het_Hom_ratio, \
                mother_unique_het_num,mother_info_het_num,mother_info_het_ratio,mother_info_hom_num,mother_info_hom_ratio, 
                father_unique_num,father_info_het_num,father_info_het_ratio,father_info_hom_num,father_info_hom_ratio]
            final_result = '\t'.join([str(i) for i in final_result])
            fo.write('%s\n'%(final_result))

def Chr_origin_analysis_main(target_chr, df_target_chr, df_target_chr_anno, targetChr_outputDir):
    if target_chr == 'chrX':
        father_unique_marker = '1/1'
    else:
        father_unique_marker = '0/1'

    if len(father) == 1 and len(mother) == 1:
        father_unique_df1 = df_target_chr[(df_target_chr[father[0]].str.contains(father_unique_marker)) & (df_target_chr[mother[0]].str.contains('0/0'))]
        mother_unique_df1 = df_target_chr[(df_target_chr[mother[0]].str.contains('0/1')) & (df_target_chr[father[0]].str.contains('0/0'))]
        father_unique_df2 = df_target_chr_anno[df_target_chr_anno[father[0]].str.contains(father_unique_marker)]
        mother_unique_df2 = df_target_chr_anno[df_target_chr_anno[mother[0]].str.contains('0/1')]

        sta_result4_1 = os.path.join(targetChr_outputDir,target_chr+'_origin_trios.txt')
        Chr_origin_analysis(df_target_chr, target_chr, father_unique_df1, mother_unique_df1, sta_result4_1)
        if paras.plot == 'Y':
            father_plot_out_1 = os.path.join(targetChr_outputDir,target_chr+'_father_origin_trios.png')
            TrioXXY_plot.Chr_origin_plot('father', sta_result4_1, father_plot_out_1)
            mother_plot_out_1 = os.path.join(targetChr_outputDir,target_chr+'_mother_origin_trios.png')
            TrioXXY_plot.Chr_origin_plot('mother', sta_result4_1, mother_plot_out_1)

        sta_result4_2 = os.path.join(targetChr_outputDir,target_chr+'_origin_maf.txt')
        Chr_origin_analysis(df_target_chr, target_chr, father_unique_df2, mother_unique_df2, sta_result4_2)
        if paras.plot == 'Y':
            father_plot_out2 = os.path.join(targetChr_outputDir,target_chr+'_father_origin_maf.png')
            TrioXXY_plot.Chr_origin_plot('father', sta_result4_2, father_plot_out2)
            mother_plot_out2 = os.path.join(targetChr_outputDir,target_chr+'_mother_origin_maf.png')
            TrioXXY_plot.Chr_origin_plot('mother', sta_result4_2, mother_plot_out2)

    elif len(father) == 1 and len(mother) == 0:
        mother_unique_df2 = pd.DataFrame()
        father_unique_df2 = df_target_chr_anno[df_target_chr_anno[father[0]].str.contains(father_unique_marker)]

        sta_result4_2 = os.path.join(targetChr_outputDir,target_chr+'_origin_maf.txt')
        Chr_origin_analysis(df_target_chr, target_chr,father_unique_df2, mother_unique_df2, sta_result4_2)
        if paras.plot == 'Y':
            father_plot_out2 = os.path.join(targetChr_outputDir,target_chr+'_father_origin_maf.png')
            TrioXXY_plot.Chr_origin_plot('father', sta_result4_2, father_plot_out2)

    elif len(mother) == 1 and len(father) == 0:
        father_unique_df2 = pd.DataFrame()
        mother_unique_df2 = df_target_chr_anno[df_target_chr_anno[mother[0]].str.contains('0/1')]

        sta_result4_2 = os.path.join(targetChr_outputDir,target_chr+'_origin_maf.txt')
        Chr_origin_analysis(df_target_chr, target_chr,father_unique_df2, mother_unique_df2, sta_result4_2)
        if paras.plot == 'Y':
            mother_plot_out2 = os.path.join(targetChr_outputDir,target_chr+'_mother_origin_maf.png')
            TrioXXY_plot.Chr_origin_plot('mother', sta_result4_2, mother_plot_out2)

def plot_Chr_marker(target_chr, parent, parent_unique_df, _type, sample, parent_Chr_plot, targetChr_outputDir):
    parent_info_hom = parent_unique_df[parent_unique_df[sample].str.contains('1/1')]
    parent_info_het = parent_unique_df[parent_unique_df[sample].str.contains('0/1')]
    parent_info_ref = parent_unique_df[parent_unique_df[sample].str.contains('0/0')]

    parent_info_hom_file = os.path.join(targetChr_outputDir, sample + '_' + parent + '_info_hom_' + _type + '.txt')
    parent_info_hom.to_csv(parent_info_hom_file,index=False,sep='\t')
    parent_info_het_file = os.path.join(targetChr_outputDir, sample + '_' + parent + '_info_het_' + _type + '.txt')
    parent_info_het.to_csv(parent_info_het_file,index=False,sep='\t')
    parent_info_ref_file = os.path.join(targetChr_outputDir, sample + '_' + parent + '_info_ref_' + _type + '.txt')
    parent_info_ref.to_csv(parent_info_ref_file,index=False,sep='\t')

    cmd = "Rscript %s/Chromosome_plot.R %s %s %s %s %s %s" %(
                    './script', target_chr.lstrip('chr'), gbandFile, parent_info_het_file, parent_info_hom_file, parent_info_ref_file, parent_Chr_plot)
    os.system(cmd)

def marker_distribution_main(target_chr, df_target_chr, df_target_chr_anno, targetChr_outputDir):
    if target_chr == 'chrX':
        father_unique_marker = '1/1'
    else:
        father_unique_marker = '0/1'

    if len(father) == 1 and len(mother) == 1:
        father_unique_df1 = df_target_chr[(df_target_chr[father[0]].str.contains(father_unique_marker)) & (df_target_chr[mother[0]].str.contains('0/0'))]
        mother_unique_df1 = df_target_chr[(df_target_chr[mother[0]].str.contains('0/1')) & (df_target_chr[father[0]].str.contains('0/0'))]
        father_unique_df2 = df_target_chr_anno[df_target_chr_anno[father[0]].str.contains(father_unique_marker)]
        mother_unique_df2 = df_target_chr_anno[df_target_chr_anno[mother[0]].str.contains('0/1')]

        for sample in child_list:
            father_Chr_plot1 = os.path.join(targetChr_outputDir, target_chr + '_' +sample + '_father_muta_plot_trios.png')
            plot_Chr_marker(target_chr, 'father', father_unique_df1, 'trio', sample, father_Chr_plot1, targetChr_outputDir)
            mother_Chr_plot1 = os.path.join(targetChr_outputDir, target_chr + '_' +sample + '_mother_muta_plot_trios.png')
            plot_Chr_marker(target_chr, 'mother', mother_unique_df1, 'trio', sample, mother_Chr_plot1, targetChr_outputDir)

            father_Chr_plot2 = os.path.join(targetChr_outputDir, target_chr + '_' +sample + '_father_muta_plot_maf.png')
            plot_Chr_marker(target_chr, 'father', father_unique_df2, 'maf', sample, father_Chr_plot2, targetChr_outputDir)
            mother_Chr_plot2 = os.path.join(targetChr_outputDir, target_chr + '_' +sample + '_mother_muta_plot_maf.png')
            plot_Chr_marker(target_chr, 'mother', mother_unique_df2, 'maf', sample, mother_Chr_plot2, targetChr_outputDir)

    elif len(father) == 1 and len(mother) == 0:
        father_unique_df2 = df_target_chr_anno[df_target_chr_anno[father[0]].str.contains(father_unique_marker)]
        for sample in child_list:
            father_Chr_plot2 = os.path.join(targetChr_outputDir, target_chr + '_' +sample + '_father_muta_plot_maf.png')
            plot_Chr_marker(target_chr, 'father', father_unique_df2, 'maf', sample, father_Chr_plot2, targetChr_outputDir)

    elif len(mother) == 1 and len(father) == 0:
        mother_unique_df2 = df_target_chr_anno[df_target_chr_anno[mother[0]].str.contains('0/1')]
        for sample in child_list:
            mother_Chr_plot2 = os.path.join(targetChr_outputDir, target_chr + '_' +sample + '_mother_muta_plot_maf.png')
            plot_Chr_marker(target_chr, 'mother', mother_unique_df2, 'maf', sample, mother_Chr_plot2, targetChr_outputDir)

if __name__ == "__main__":
    home = os.getcwd()
    taskID = paras.task_ID
    pedFile = paras.ped
    MAF = float(paras.MAF)
    uploadDir = os.path.join(home, 'upload')
    taskDir = os.path.join(uploadDir, taskID)
    os.path.exists(taskDir) or os.mkdir(taskDir)
    ped_file = os.path.join(taskDir, pedFile)
    outputDir = os.path.join(home, 'output')
    os.path.exists(outputDir) or os.mkdir(outputDir)
    task_outputDir = os.path.join(outputDir, taskID)
    os.path.exists(task_outputDir) or os.mkdir(task_outputDir)
    X_outputDir = os.path.join(task_outputDir,'chrX')
    os.path.exists(X_outputDir) or os.mkdir(X_outputDir)
    MAF_site_1kg_file = os.path.join(home,'script','1kg_X_MAF_0.1.txt')
    gbandFile = os.path.join(home,'script','hg19_G-band.txt')

# Step1: Preprocess
    # Step1.1: Sample check
    chr_list = ['chr' + str(i) for i in range(1,23)] + ['chrX']
    h_par = get_PAR()
    sampleCheckResult, Step1_Note, father, mother, child_list, h_ped = sampleCheck()
    if sampleCheckResult == 'not pass':
        print ('%s%s'%(taskID, ' # Step 1: Input error !'))
        for note in Step1_Note:
            print(note)
    else:
        print ('%s%s'%(taskID, ' # Step 1.1: Found vcf file for every sample in ped file !'))
        sample_list = father + mother + child_list
        # Step1.2: Process VCFs
        muta_file = os.path.join(X_outputDir, '01_'+taskID+'_muta.txt')
        VCF_process_main()
        print ('%s%s'%(taskID, ' # Step 1.2: VCF processing over !'))

        # Step2: Genotype quality
        muta_quality_file = os.path.join(X_outputDir, '02_'+taskID+'_muta_quality.txt')
        Genotype_quality_main()
        print ('%s%s'%(taskID, ' # Step 2: Genotype quality filtering over !'))
        muta_df = pd.read_csv(muta_quality_file,sep='\t')
        df_X = get_chr_df("chrX")
        df_Y = get_chr_df("chrY")

        # Step3: Depth analysis
        Depth_analysis_main('chrX', df_X, X_outputDir)
        print ('%s%s'%(taskID, ' # Step 3: Depth analysis over !'))

        # Step4: Y muta analysis
        Y_hom_muta_number_sta()
        print ('%s%s'%(taskID, ' # Step 4: Y mutation analysis over !'))

        # Step5: X origin analysis
        if len(father) == 1 or len(mother) == 1:
            Muta_anno(df_X, 'chrX', X_outputDir)
            X_muta_fre_anno_file = os.path.join(X_outputDir, 'chrX_'+str(MAF)+'_muta_fre_anno.txt')
            df_X_anno = pd.read_csv(X_muta_fre_anno_file, sep='\t')

            Chr_origin_analysis_main('chrX', df_X, df_X_anno, X_outputDir)
            print ('%s%s'%(taskID, ' # Step 5: X origin analysis over !'))
            # Step6: X marker distribution
            if paras.plot == 'Y':
                marker_distribution_main('chrX', df_X, df_X_anno, X_outputDir)
                print ('%s%s'%(taskID, ' # Step 6: X marker plot over !'))

        # repeat depth analysis for autosomes
        target_chr_list = paras.A_chr
        for target_chr in target_chr_list:
            if target_chr in ['chr' + str(i) for i in range(1,23)]:
                targetChr_outputDir = os.path.join(task_outputDir,target_chr)
                os.path.exists(targetChr_outputDir) or os.mkdir(targetChr_outputDir)
                df_targetChr = get_chr_df(target_chr)
                Depth_analysis_main(target_chr, df_targetChr, targetChr_outputDir)
                print ('%s%s%s%s'%(taskID, ', ', target_chr, ' Depth analysis over !'))