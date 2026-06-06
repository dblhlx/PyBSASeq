"""
PyBSASeq, version 3.14
Created on Fri Oct  5 08:22:16 2018
Updated on 
@author: Jianbo Zhang
"""
import os
import sys
import time
import datetime
import argparse
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel
from scipy.signal import savgol_filter
import ast
try:
    from fisher import pvalue_npy
except ImportError:
    try:
        from scipy.stats import fisher_exact
    except ImportError:
        print('Either the module \'Fisher\' (https://github.com/brentp/fishers_exact_test) or \'fisher_exact\' from \'scipy.stat\' is needed to performe Fisher\'s exact tests.')
        sys.exit()


def sm_allelefreq(pop_struc, bulk_size, rep):
    '''
    This function is for the calculation of the ALT allele frequency in the bulks via simulation
    under the null hypothesis (the SV is not associated with the trait).

    An AA, Aa, and aa individual carries 0%, 50%, and 100% of the alt (a) allele, respectively.
    If A/a is not associated with the trait (null hypothesis), the AA:Aa:aa ratios are 0.25:0.5:0.25, 0.5:0:0.5, 
    and 0.5:0.5:0, respectively, in both bulks of an F2 population, an RIL population, or a back crossed population.
    '''
    freq_list = []
    pop = [0.0, 0.5, 1.0]

    if pop_struc == 'F2':
        prob = [0.25, 0.5, 0.25]
    elif pop_struc == 'RIL':
        prob = [0.5, 0.0, 0.5]
    elif pop_struc == 'BC':
        prob = [0.5, 0.5, 0.0]

    for __ in range(rep):
        alt_freq = np.random.choice(pop, bulk_size, p=prob).mean()
        freq_list.append(alt_freq)

    return sum(freq_list)/len(freq_list)


def sort_chrm(l):
    a, b = [], []

    for eml in l:
        if eml.isdigit():
            a.append(eml)
        else:
            b.append(eml)
    a.sort(key=int)
    b.sort()
    a.extend(b)

    return a


def chrm_filtering(df, chromosome_list):
    # Many reference genomes contain unmapped fragments that tend to be small and are not informative for SV-trait association, filtering them out makes the chromosome list more readable.
    all_chrm_ids, all_chrm_range, small_chrm_ids = [], [], []
    for chrm in chromosome_list:
        chrm_size = df[df.CHROM==chrm]['POS'].max().item()

        if chrm_size <= min_frag_size:
            small_chrm_ids.append(chrm)
        else:
            all_chrm_ids.append(chrm)
            all_chrm_range.append([1, chrm_size])       # Startpoint and endpoint

    return [all_chrm_ids, all_chrm_range, small_chrm_ids]


def select_chrms(df, rgn):
    # Create a chromosome list, which can be very long because of the unmapped fragments
    global chrm_prnt_dict
    chrm_prnt_dict = {}
    raw_chrm_ids_ori = df['CHROM'].unique().tolist()
    raw_chrm_ids = sort_chrm(raw_chrm_ids_ori)

    # Filter out chromosomes and unmapped fragments smaller than the sliding window
    # Make the chromosome list more readable and meaningful
    chrm_id_range = chrm_filtering(df, raw_chrm_ids)
    chrm_ids = chrm_id_range[0]
    chrm_range = chrm_id_range[1]
    small_frags = chrm_id_range[2]

    m = 1
    for eml in chrm_ids:
        chrm_prnt_dict[eml] = str(m)
        m += 1

    if small_frags != []:
        print('\nThe chromosomes/fragments below are filtered out because of their small sizes:')
        print(small_frags, '\n')

    dsrd_chrm_ids, dsrd_chrm_range, dsrd_chrm_sizes = [], [], []
    # Chromosomal region arguments are specified
    if rgn[0] == -1:
        print(f'The chromosomes below are greater than {min_frag_size} bp and are suitable for BSA-Seq analysis:')
        print(chrm_ids,'\n')
        print('Although a subset of the above chromosomes can be selected for analysis in next step, it is strongly recommended to have all the chromosomes included when running the script the first time.\n')

        if chrm_order == False:
            dsrd_chrm_ids = chrm_ids
        else:
        # while right_input.lower() != 'yes':
        #     input_string = input('Enter the names of the desired chromosomes for analysis in order and separate each name with a comma, or press the ENTER/RETURN key if the above list is what you want:\n')
        #     if input_string == '':
        #         dsrd_chrm_ids = chrm_ids
        #     else:
        #         dsrd_chrm_ids = [x.strip() for x in input_string.split(',')]

        #     print('Invalid chromosome names, if any, will be removed.\n')

            # Filter out possible invalid chromosome name(s)
            invalid_chrm_list, frags = [], []
            for ch in chrm_order:
                if ch not in chrm_ids:
                    invalid_chrm_list.append(ch)
                elif ch in small_frags:
                    frags.append(ch)
                else:
                    dsrd_chrm_ids.append(ch)

            if invalid_chrm_list != []:
                print('These chromosome names are invalid:', invalid_chrm_list,'\n')

            if frags != []:
                print(f'The chromosomes below are less than {min_frag_size} bp and are not suitable for BSA-Seq analysis:')
                print(frags,'\n')

            print('Chromosomes for analysis:')
            print(dsrd_chrm_ids,'\n')

            # right_input = input('Are the chromosome names in the above list in the right order (yes or no)?\n')
            # print('\n')

        # Create dictionary/list containing the sizes of all the chromosomes
        i = 0
        for ch in chrm_ids:
            if ch in dsrd_chrm_ids:
                temp = chrm_range[i]
                dsrd_chrm_range.append(temp)
                dsrd_chrm_sizes.append(temp[1])
            i += 1

    # Handle the cases in which one or more chromosomal regions are specified
    else:
        l, n = len(rgn), len(rgn) % 3
        if n != 0:
            rgn = rgn[0:l-n]
            print('Each chromosome name should be followed by the starting and ending points of the interested region.')

        i = 0
        while i < l:
            if rgn[i] == 1000:
                ch_id = 'X'
            elif rgn[i] == 1001:
                ch_id = 'Y'
            elif rgn[i] == 1002:
                ch_id = 'Z'
            elif rgn[i] == 1003:
                ch_id = 'W'
            elif rgn[i] == 1004:
                ch_id = 'U'
            elif rgn[i] == 1005:
                ch_id = 'V'
            else:
                ch_id = str(rgn[i])

            if ch_id not in chrm_ids:
                print(ch_id, 'is not a valid chromosome name.')
            else:
                temp_chrm_size = chrm_range[chrm_ids.index(ch_id)][1]
                if rgn[i+1] < 1:
                    rgn[i+1] = 1
                if rgn[i+2] > temp_chrm_size:
                    rgn[i+2] = temp_chrm_size

                if rgn[i+2] - rgn[i+1] > min_frag_size:
                    dsrd_chrm_ids.append(ch_id)
                    dsrd_chrm_range.append([rgn[i+1], rgn[i+2]])
                    dsrd_chrm_sizes.append(rgn[i+2]-rgn[i+1]+1)
                else:
                    print(f'The size of the interested region on chromosome {ch_id} should be greater than ', {min_frag_size}, 'bp.')

            i += 3

    return [dsrd_chrm_ids, dsrd_chrm_range, dsrd_chrm_sizes]


def bulk_names(df):
    global header
    header = df.columns.values.tolist()
    bulks = []

    # Obtain the bulk IDs from the header
    # try:
    for ftr_name in header:
        if ftr_name.endswith('.AD') or ftr_name.endswith('_AD'):
            bulks.append(ftr_name.split('.')[0])

    return bulks

    # except (NameError, IndexError):
    #     print('The allele depth (AD) field is missing. Please include the AD field in the input file.')
    #     sys.exit()


def xticks_property(l):
    # Automatically adjust xticks
    # Calculate distance between xticks and determine the tick unit
    lrgst_chrm_length = max(l)
    div_list = [500000000, 200000000, 100000000, 50000000, 20000000, 10000000, 5000000, 2000000, 1000000]
    max_xticks, min_xticks= 7, 3

    if lrgst_chrm_length/div_list[0] > max_xticks:
        div_unit = 1000000000
        length_unit = 'Gb'
        rmzero = 1e-9
    elif lrgst_chrm_length/div_list[-1] < min_xticks:
        div_unit = 100000
        length_unit = '\u00D7100 kb'
        rmzero = 1e-5
    else:
        for i in div_list:
            if lrgst_chrm_length/i <= max_xticks and lrgst_chrm_length/i >= min_xticks:
                div_unit = i
                if i/1000000 >= 100:
                    length_unit = '\u00D7100 Mb'
                    rmzero = 1e-8
                elif i/1000000 >= 10:
                    length_unit ='\u00D710 Mb'
                    rmzero = 1e-7
                elif i/1000000 >= 1:
                    length_unit = 'Mb'
                    rmzero = 1e-6

                break

    return [div_unit, rmzero, length_unit]


def sv_filtering(df):
    print(f'Perform SV filtering - {sample}')
    df = df.copy()

    # Identify SVs not informative and remove these SVs from the dataframe
    df_ignored = df[~df.CHROM.isin(selected_chrms)]
    df_ignored.to_csv(os.path.join(filtering_path, 'ignored.csv'), index=None)
    df = df.drop(index=df_ignored.index)

    # Identify SVs with an 'NA' value(s) and remove these SVs from the dataframe
    df_na = df[df.isnull().any(axis=1)]
    df_na.to_csv(os.path.join(filtering_path, 'na.csv'), index=None)
    df.dropna(inplace=True)
    misc.append(['Number of SVs after NA drop - '+sample_name, len(df.index)])
    misc.append([''])

    # Identify SVs with more than one ALT allele
    df_m_alts = df[df.ALT.str.contains(',')]

    # Identify SVs with a single ALT allele
    df_1_alt = df.drop(index=df_m_alts.index)

    # Identify one-ALT SVs with zero REF read in both bulks
    df_1_alt_fake = df_1_alt[(df_1_alt[fb_ad].str.startswith('0')) & (df_1_alt[sb_ad].str.startswith('0'))]
    df_1_alt_fake.to_csv(os.path.join(filtering_path, 'fake_1alt.csv'), index=None)

    # Remove one-ALT SVs with zero REF read in both bulks
    df_1_alt_real = df_1_alt.drop(index=df_1_alt_fake.index)
    df_1_alt_real.to_csv(os.path.join(filtering_path, 'real_1alt.csv'), index=None)

    # Using 'str.count' make the code below simpler and easier to understand
    df_3_alts = df_m_alts[df_m_alts.ALT.str.count(',') >= 2]
    df_3_alts.to_csv(os.path.join(filtering_path, '3alts.csv'), index=None)
    df_2_alts = df_m_alts.drop(index=df_3_alts.index)

    # A two-ALT SV is a real SV if the REF read is zero in both bulks
    df_2_alts_real = df_2_alts[(df_2_alts[fb_ad].str.startswith('0')) & (df_2_alts[sb_ad].str.startswith('0'))]

    # The two-ALT SV may be caused by allele heterozygosity if the REF read is not zero
    # Repetitive sequences in the genome or sequencing artifacts are other possibilities
    df_2_alts_het = df_2_alts.drop(index=df_2_alts_real.index)
    df_2_alts_het.to_csv(os.path.join(filtering_path, '2alts_het.csv'), index=None)

    # Making a copy to suppress the warning message. Updating the AD values of these SVs is required
    df_2_alts_real = df_2_alts_real.copy()
    df_2_alts_real.to_csv(os.path.join(filtering_path, 'real_2alts_before.csv'), index=None)

    # Update the AD values of the above SVs by removing the REF read that is zero (i.e. remove '0,' from '0,x,y')
    if not df_2_alts_real.empty:
        df_2_alts_real[fb_ad] = df_2_alts_real[fb_ad].str.slice(start=2)
        df_2_alts_real[sb_ad] = df_2_alts_real[sb_ad].str.slice(start=2)
        df_2_alts_real[['REF', 'ALT']] = df_2_alts_real.ALT.str.split(',', expand=True)
        df_2_alts_real.to_csv(os.path.join(filtering_path, 'real_2alts_after.csv'), index=None)

    # Concatenate 1ALT_Real and 2ALT_Real
    df = pd.concat([df_1_alt_real, df_2_alts_real])
    # df = df_1_alt_real.copy()

    # Filter out 1 bp deletions
    indel_1 = df[(df.REF.str.contains(r'\*')) | (df.ALT.str.contains(r'\*'))]
    # df = df.drop(index=indel_1.index)

    # Remove InDels from the SV dataframe
    indel_2 = df[(df.REF.str.len()>1) | (df.ALT.str.len()>1)]
    # df = df.drop(index=indel_2.index)

    df_indel = pd.concat([indel_1, indel_2])
    df_indel.to_csv(os.path.join(filtering_path, 'indel.csv'), index=None)

    # Obtain REF reads, ALT reads, and locus reads of each SV
    df[[fb_ad_ref, fb_ad_alt]] = df[fb_ad].str.split(',', expand=True).astype(int)
    df[fb_ld] = df[fb_ad_ref] + df[fb_ad_alt]
    df[[sb_ad_ref, sb_ad_alt]] = df[sb_ad].str.split(',', expand=True).astype(int)
    df[sb_ld] = df[sb_ad_ref] + df[sb_ad_alt]

    # Filter out the SVs with zero locus reads in either bulk
    sv_0ld = df[(df[fb_ld] == 0) | (df[sb_ld] == 0)]
    sv_0ld.to_csv(os.path.join(filtering_path, '0ld.csv'), index=None)
    df = df.drop(index=sv_0ld.index)

    # Filter out SVs in which GT and AD are not consistent
    df[[fb_gt_ref, fb_gt_alt]] = df[fb_gt].str.split('/|\\|', expand=True)
    df[[sb_gt_ref, sb_gt_alt]] = df[sb_gt].str.split('/|\\|', expand=True)
    gt_ad = df[(df[fb_gt_ref]==df[fb_gt_alt]) & (df[sb_gt_ref]==df[sb_gt_alt]) & (df[fb_gt_ref]==df[sb_gt_ref])]
    df = df.drop(index=gt_ad.index)
    gt_ad.to_csv(os.path.join(filtering_path, 'gt_ad.csv'), index=None)

    # Filter out SVs in which bulk GT and REF/ALT are not consistent
    gt = df[~(((df[fb_gt_ref]==df.REF) | (df[fb_gt_ref]==df.ALT)) & \
               ((df[fb_gt_alt]==df.REF) | (df[fb_gt_alt]==df.ALT)) & \
               ((df[sb_gt_ref]==df.REF) | (df[sb_gt_ref]==df.ALT)) & \
               ((df[sb_gt_alt]==df.REF) | (df[sb_gt_alt]==df.ALT)))]
    gt.to_csv(os.path.join(filtering_path, 'gt.csv'), index=None)
    df = df.drop(index=gt.index)

    df.sort_values(['ChrmSortID', 'POS'], inplace=True)

    print(f'SV filtering completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return df


def read_filter(df):
    fb_rep_ld, sb_rep_ld = df[fb_ld].mean() + 3*df[fb_ld].std(), df[sb_ld].mean() + 3*df[sb_ld].std() 
    df_rep = df[(df[fb_ld]>fb_rep_ld) | (df[sb_ld]>sb_rep_ld)]
    df_rep.to_csv(os.path.join(filtering_path, 'repetitive_seq.csv'), index=None)
    df = df.drop(index=df_rep.index)

    df_lowread = df[(df[fb_ld]<3) | (df[sb_ld]<3)]
    df_lowread.to_csv(os.path.join(filtering_path, 'lowread.csv'), index=None)
    df = df.drop(index=df_lowread.index)

    return df


def sv_filtering_final(df):
    svs_lowq = df[(df[fb_gq]<gq_value) | (df[sb_gq]<gq_value)]
    svs_lowq.to_csv(os.path.join(filtering_path, 'svs_lowq.csv'), index=None)
    df = df.drop(index=svs_lowq.index)

    sv_highratio = df[(df[fb_ad_alt]/df[fb_ad_ref]>2) & (df[sb_ad_alt]/df[sb_ad_ref]>2)]
    sv_highratio.to_csv(os.path.join(filtering_path, 'highratio.csv'), index=None)
    df = df.drop(index=sv_highratio.index)

    sv_lowratio = df[(df[fb_ad_alt]/df[fb_ad_ref]<0.5) & (df[sb_ad_alt]/df[sb_ad_ref]<0.5)]
    sv_lowratio.to_csv(os.path.join(filtering_path, 'lowratio.csv'), index=None)
    df = df.drop(index=sv_lowratio.index)

    # SVs with very high LD are likely from the repetitive genetic elements in the genome
    # Remove SVs that could be from the repetitive elements
    # print([df[fb_ld].mean(), df[fb_ld].std(), df[sb_ld].mean(), df[sb_ld].std()])
    df = read_filter(df)

    gt_nm = df[df[fb_gt]!=df[sb_gt]]
    gt_nm.to_csv(os.path.join(filtering_path, 'gt_nm.csv'), index=None)
    # df = df.drop(index=gt_nm.index)

    highly_ed = df[((df[fb_ad_alt]/df[fb_ad_ref]>2) & (df[sb_ad_alt]/df[sb_ad_ref]<0.5)) | ((df[fb_ad_alt]/df[fb_ad_ref]<0.5) & (df[sb_ad_alt]/df[sb_ad_ref]>2))]
    highly_ed.to_csv(ed_file, index=None)

    bulk_homo_svs = sv_highratio = df[((df[fb_ad_alt]==0) & (df[sb_ad_ref]==0)) | ((df[fb_ad_ref]==0) & (df[sb_ad_alt]==0))]
    bulk_homo_svs.to_csv(bulk_homo_svs_file, index=None)

    return df


def g_statistic_array(o1, o3, o2, o4):
    '''
    Calculate G-statistic using numpy arrays as input
    o1 - o4 are 4 numpy arrays of observed values that are greater than or equal to zero
    '''
    # Ignore errors caused by 'Divide by zero' or logarithm of zero, and let numpy.where to handle these situations
    np.seterr(all='ignore')

    # Calculate the expected values under the null hypothesis
    e1 = np.where(o1+o2+o3+o4!=0, (o1+o2)*(o1+o3)/(o1+o2+o3+o4), 0)
    e2 = np.where(o1+o2+o3+o4!=0, (o1+o2)*(o2+o4)/(o1+o2+o3+o4), 0)
    e3 = np.where(o1+o2+o3+o4!=0, (o3+o4)*(o1+o3)/(o1+o2+o3+o4), 0)
    e4 = np.where(o1+o2+o3+o4!=0, (o3+o4)*(o2+o4)/(o1+o2+o3+o4), 0)

    # Calculate the log-likelihood ratios
    llr1 = np.where(o1/e1>0, 2*o1*np.log(o1/e1), 0.0)
    llr2 = np.where(o2/e2>0, 2*o2*np.log(o2/e2), 0.0)
    llr3 = np.where(o3/e3>0, 2*o3*np.log(o3/e3), 0.0)
    llr4 = np.where(o4/e4>0, 2*o4*np.log(o4/e4), 0.0)

    return np.where(e1*e2*e3*e4==0, 0.0, llr1+llr2+llr3+llr4)


def sv_ci(row):
    '''
    Estimate the thresholds of Δ(allele frequency) and G-statistic value using its simulated AD values for each SV in the dataset.
    Each row in the dataset represents an SV.
    This function will be used only when the module 'Fisher' is available.
    '''
    # Create an array with 10000 (rep) simulated ALT reads of the SV - first bulk
    yb_ld, eb_ld = row[fb_ld], row[sb_ld]
    yb_alt_array = np.random.binomial(yb_ld, fb_freq, rep)
    # yb_ld_array = np.full(rep, row[fb_ld])
    # Create an array with 10000 (rep) simulated ALT reads of the SV - second bulk
    eb_alt_array = np.random.binomial(eb_ld, sb_freq, rep)
    # eb_ld_array = np.full(rep, row[sb_ld])

    # Create simulated allele frequency and Δ(allele frequency) arrays of the SV
    # yb_af_array = yb_alt_array/yb_ld
    # eb_af_array = eb_alt_array/eb_ld
    daf_array = eb_alt_array/eb_ld - yb_alt_array/yb_ld
    daf_abs_array = np.abs(daf_array)

    # Create a G-statistic array of the SV
    gs_array = g_statistic_array(yb_alt_array, yb_ld-yb_alt_array, eb_alt_array, eb_ld-eb_alt_array)

    # Obtain the percentile of the above arrays
    # ci_yb_af = np.percentile(yb_af_array, percentile_ci)
    # ci_eb_af = np.percentile(eb_af_array, percentile_ci)
    ci_daf = np.percentile(daf_array, percentile_ci)
    ci_gs = np.percentile(gs_array, percentile_th)
    ci_daf_abs = np.percentile(daf_abs_array, percentile_th)

    # return [ci_yb_af, ci_eb_af, ci_daf, ci_gs]
    return [ci_daf, ci_gs, ci_daf_abs]


def sv_ci_scipy(row):
    '''
    Perform Fisher's exact test for each SV in the dataset using its actual AD values in both bulks 
    and estimate the thresholds of Δ(allele frequency) and G-statistic value using its simulated AD values.
    Each row in the dataset represents an SV.
    This function will be used if the 'Fisher' module is not installed. It is slow for large dataset.
    '''
    # Perform Fisher's exact test for each SV using the actual REF/ALT reads
    try:
        __, fe = fisher_exact([[row[fb_ad_ref], row[fb_ad_alt]], [row[sb_ad_ref], row[sb_ad_alt]]])
    except TypeError:
        fe = 'NA'

    # Perform Fisher's exact test for each SV using the simulated REF/ALT reads
    try:
        __, sm_fe = fisher_exact([[row[sm_fb_ad_ref], row[sm_fb_ad_alt]], [row[sm_sb_ad_ref], row[sm_sb_ad_alt]]])
    except TypeError:
        sm_fe = 'NA'

    a = [fe, sm_fe]

    if test == True:
        a.extend(sv_ci(row))

    return a


def sv_ci_gw(df):
    # Using this function to calculate the threshold if 'Fisher' is not available
    print('Estimate the genome-wide thresholds of the sSV/totalSV, G-statistic, and \u0394(allele frequency) at the SV level.')
    ratio_list, gs_list, daf_list, daf_list_abs = [], [], [], []
    sv_smpl = df.sample(sv_per_sw, replace=True)

    sv_smpl['STAT'] = sv_smpl.apply(sv_ci, axis=1)

    # Create new columns containing allele frequency, Δ(allele frequency) confidence intervals, and G-statistic thresholds
    # sv_smpl[[fb_af_ci, sb_af_ci, 'DAF_CI', 'GS_CI']] = pd.DataFrame(sv_smpl.STAT.values.tolist(), index=sv_smpl.index)
    sv_smpl[['DAF_CI', 'GS_Thrshld', 'DAF_abs_Thrshld']] = pd.DataFrame(sv_smpl.STAT.values.tolist(), index=sv_smpl.index)

    sv_smpl['DAF_CI_LB'] = sv_smpl['DAF_CI'].apply(lambda x: x[0]).astype(float)
    sv_smpl['DAF_CI_UB'] = sv_smpl['DAF_CI'].apply(lambda x: x[1]).astype(float)

    return [sv_smpl['DAF_CI_LB'].mean(), sv_smpl['DAF_CI_UB'].mean(), sv_smpl['GS_Thrshld'].mean(), sv_smpl['DAF_abs_Thrshld'].mean()]


def thresholds_gw_approximate(df):
    # Using this function to calculate the threshold if 'Fisher' is not available
    print('Estimate the approximate genome-wide thresholds of the sSV/totalSV, G-statistic, and \u0394(allele frequency) at the sliding window level.')
    ratio_list, gs_list, daf_list, daf_list_abs = [], [], [], []
    for __ in range(rep):
        sv_smpl = df.sample(sv_per_sw, replace=True)
        ssv_smpl = sv_smpl[sv_smpl['sm_FE_P']<sm_alpha]

        ratio_list.append(len(ssv_smpl.index)/sv_per_sw)
        gs_list.append(np.mean(g_statistic_array(sv_smpl[fb_ld]-sv_smpl[sm_fb_ad_alt], sv_smpl[sm_fb_ad_alt], sv_smpl[sb_ld]-sv_smpl[sm_sb_ad_alt], sv_smpl[sm_sb_ad_alt])))
        daf_list.append(np.mean(sv_smpl[sm_sb_ad_alt]/sv_smpl[sb_ld]-sv_smpl[sm_fb_ad_alt]/sv_smpl[fb_ld]))
        daf_list_abs.append(np.mean(np.absolute(sv_smpl[sm_sb_ad_alt]/sv_smpl[sb_ld]-sv_smpl[sm_fb_ad_alt]/sv_smpl[fb_ld])))

    misc.append(['Genome-wide sSV/totalSV ratio threshold', np.percentile(ratio_list, percentile_th)])
    misc.append(['Genome-wide G-statistic threshold', np.percentile(gs_list, percentile_th)])
    misc.append(['Genome-wide Delta(AF) threshold', np.percentile(daf_list, percentile_ci)])
    misc.append(['Genome-wide absolute Delta(AF) threshold', np.percentile(daf_list_abs, percentile_th)])
    print(f'Threshold calculation completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return [np.percentile(ratio_list, percentile_th), np.percentile(gs_list, percentile_th), np.percentile(daf_list, percentile_ci), np.percentile(daf_list_abs, percentile_th)]


def thresholds_gw(df):
    # For the calculation of the genome-wide threshold if 'Fisher' is installed
    print('Estimate the genome-wide thresholds of the sSV/totalSV, G-statistic, and \u0394(allele frequency).')
    gw_ratio_list, gw_gs_list, gw_daf_list, gw_daf_abs_list = [], [], [], []

    for __ in range(rep):
        sv_smpl = df.sample(sv_per_sw, replace=True)

        gw_fb_ad_alt_arr = np.random.binomial(sv_smpl[fb_ld], fb_freq).astype(np.uint32)
        gw_fb_ad_ref_arr = sv_smpl[fb_ld].to_numpy().astype(np.uint32) - gw_fb_ad_alt_arr
        gw_sb_ad_alt_arr = np.random.binomial(sv_smpl[sb_ld], sb_freq).astype(np.uint32)
        gw_sb_ad_ref_arr = sv_smpl[sb_ld].to_numpy().astype(np.uint32) - gw_sb_ad_alt_arr

        __, __, gw_sm_fe_p_arr = pvalue_npy(gw_fb_ad_alt_arr, gw_fb_ad_ref_arr, gw_sb_ad_alt_arr, gw_sb_ad_ref_arr)
        # gw_sm_fe_OR_arr = (gw_fb_ad_alt_arr * gw_sb_ad_ref_arr)/(gw_fb_ad_ref_arr * gw_sb_ad_alt_arr)

        # Calculate sSV/totalSV ratio
        sSV_arr = np.where(gw_sm_fe_p_arr<sm_alpha, 1, 0)
        gw_ratio_list.append(np.mean(sSV_arr))

        # Calculate G-statistic
        gw_gs_arr = g_statistic_array(gw_fb_ad_ref_arr, gw_fb_ad_alt_arr, gw_sb_ad_ref_arr, gw_sb_ad_alt_arr)
        gw_gs_list.append(np.mean(gw_gs_arr))

        # Calculate allele frequency
        gw_daf_arr = gw_sb_ad_alt_arr/sv_smpl[sb_ld] - gw_fb_ad_alt_arr/sv_smpl[fb_ld]
        gw_daf_list.append(np.mean(gw_daf_arr))
        gw_daf_abs_list.append(np.mean(np.absolute(gw_daf_arr)))

    misc.append(['Genome-wide sSV/totalSV ratio threshold', np.percentile(gw_ratio_list, percentile_th)])
    misc.append(['Genome-wide G-statistic threshold', np.percentile(gw_gs_list, percentile_th)])
    misc.append(['Genome-wide Delta(AF) ratio threshold', np.percentile(gw_daf_list, percentile_ci)])
    misc.append(['Genome-wide absolute Delta(AF) ratio threshold', np.percentile(gw_daf_abs_list, percentile_th)])
    print(f'Threshold calculation completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return [np.percentile(gw_ratio_list, percentile_th), np.percentile(gw_gs_list, percentile_th), np.percentile(gw_daf_list, percentile_ci), np.percentile(gw_daf_abs_list, percentile_th)]


def thresholds_sw(df):
    df = df.copy()
    df['DI'] = di(df)
    df = df[df.DI==1]
    # For the calculation of the sliding window-specific sSV/totalSV threshold
    sw_ratio_list, sw_gs_list, sw_daf_list, sw_daf_abs_list = [], [], [], []

    # Convert the LD column to a numpy array
    sw_fb_ld_arr = df[fb_ld].to_numpy().astype(np.uint32)
    sw_sb_ld_arr = df[sb_ld].to_numpy().astype(np.uint32)

    for __ in range(rep):
        # Create new columns with simulated AD values based on the LD values
        sw_fb_ad_alt_arr = np.random.binomial(df[fb_ld], fb_freq).astype(np.uint32)
        sw_fb_ad_ref_arr = sw_fb_ld_arr - sw_fb_ad_alt_arr
        sw_sb_ad_alt_arr = np.random.binomial(df[sb_ld], sb_freq).astype(np.uint32)
        sw_sb_ad_ref_arr = sw_sb_ld_arr - sw_sb_ad_alt_arr

        # Calculate the P-value via Fisher's Exact test
        __, __, sw_sm_fe_p_arr = pvalue_npy(sw_fb_ad_alt_arr, sw_fb_ad_ref_arr, sw_sb_ad_alt_arr, sw_sb_ad_ref_arr)
        # # Calculate the odd ratio, not needed for BSA-Seq analysis
        # sw_sm_fe_OR_arr = (sw_fb_ad_alt_arr * sw_sb_ad_ref_arr)/(sw_fb_ad_ref_arr*sw_sb_ad_alt_arr)

        # Calculate sSV/totalSV
        sSV_arr = np.where(sw_sm_fe_p_arr<sm_alpha, 1, 0)
        sw_ratio_list.append(np.mean(sSV_arr))

        # Calculate G-statistic
        sw_gs_arr = g_statistic_array(sw_fb_ad_ref_arr, sw_fb_ad_alt_arr, sw_sb_ad_ref_arr, sw_sb_ad_alt_arr)
        sw_gs_list.append(np.mean(sw_gs_arr))

        # Calculate allele frequency
        sw_daf_arr = sw_sb_ad_alt_arr/sw_sb_ld_arr - sw_fb_ad_alt_arr/sw_fb_ld_arr
        sw_daf_list.append(np.mean(sw_daf_arr))
        sw_daf_abs_list.append(np.mean(np.absolute(sw_daf_arr)))

    return [np.percentile(sw_ratio_list, percentile_th), np.percentile(sw_gs_list, percentile_th), np.percentile(sw_daf_list, percentile_ci), np.percentile(sw_daf_abs_list, percentile_th)]


def zeroSV(li):
    # Replace 'divide by zero' with the nearest value.
    if li != []:
        li.append(li[-1])   # Assign the previous value to the empty sliding window if the list is not empty
    else:
        li.append('empty')  # Assign 'empty' to the first sliding windows that is empty


def replace_zero(li):
    # Replace the 'empty' placeholders at the beginning of the list with the nearest non-empty value
    i = 0
    while li[i]=='empty':
        i += 1

    j = 0
    while j < i:
        li[j] = li[i]
        j += 1


def di(df):
    # https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
    # https://www.dataquest.io/blog/settingwithcopywarning/
    if not os.path.exists(dgns_path):
        os.makedirs(dgns_path)

    chromosomes = df.CHROM.unique().tolist()
    di_list = []
    for chrm in chromosomes:
        ch = df[df.CHROM==chrm].copy()
        first_pos_value = ch['POS'].iloc[0]                 # ch.at[0, 'POS'] may not work for dataframe subsets
        pos_l1 = ch.POS.tolist()
        pos_l1.insert(0, first_pos_value-seg_size)          # ensure the first SV is included even if it located at the position 1
        pos_l1.pop()

        ch['POS1'] = pos_l1
        ch['DSTNC'] = ch.POS - ch.POS1
        l = ch.DSTNC.tolist()

        i, a = 0, []

        while i < len(l):
            if l[i] >= seg_size:
                a.append(1)
                i += 1
            else:
                a.append(0)
                temp = l[i]
                i += 1
                m = i
                while i < len(l)-1 and temp + l[i] < seg_size:
                    a.append(0)
                    temp += l[i]
                    i += 1
                if m != i:
                    if temp < seg_size and i >= len(l)-1:
                        a.append(0)
                    else:
                        a.append(1)
                    i += 1

        di_list = di_list + a

        ch['DI'] = a
        temp_df = ch[['POS', 'POS1', 'DSTNC', 'DI']]
        if type(chrm) is int:
            temp_df.to_csv(os.path.join(dgns_path, str(chrm)+'_'+str(first_pos_value)+'.csv'), index=None)
        else:
            temp_df.to_csv(os.path.join(dgns_path, chrm+'_'+str(first_pos_value)+'.csv'), index=None)

    return di_list


def bsaseq_plot(df):
    '''
    wm_list: list of warning messages
    sw_dict: a dictionary with the chromosome ID as its keys; the value of each key is a list containing the chromosome ID, the sSV/totalSV ratio in each sliding window, and the startpoint of the sliding window
    misc: miscellaneous information
    sv_region: genomic region above the threshold
    sw_data_frame: DataFrame containing sliding windows
    '''
    print('Prepare SV data for plotting via the sliding window algorithm')

    tl = [i for i in selected_chrms]
    for x in tl:
        y = df[df.CHROM==x]
        # print(x, len(y.index))
        if len(y.index) < min_sv:
            selected_chrms.remove(x)

    # print(selected_chrms)
    if selected_chrms == []:
        print('Your dataset contains too few SVs. No analysis is performed')
        sys.exit()

    # Plot layout setup
    height_ratio = [1,0.8,0.8,0.8]
    xt_pro = xticks_property(chrmSzL)
    fig, axs = plt.subplots(nrows=len(height_ratio), ncols=len(selected_chrms), figsize=(20, 12.8), sharex='col', sharey='row', gridspec_kw={'width_ratios': chrmSzL, 'height_ratios': height_ratio})

    global misc
    global sv_region, sw_data_frame
    sg_y_ratio_list = []     # Smoothed sSV/totalSV ratios of a chromosome or a selected genomic region
    sg_gs_list = []
    sg_t_gs_list = []
    sg_daf_list = []
    sg_nt_daf_list = []
    sg_pt_daf_list = []
    sw_rows = []
    wm_list, sw_dict, sv_region = [], {}, []
    hdr = df.columns.values.tolist()

    # Analyze each chromosome separately
    num_sv_on_chr = []      # List containing the sSV/totalSV, # of sSVs, and # of totalSVs of a chromosome
    ratio_peak_list = []         # List containing the peak of each chromosome 
    i = 1
    for chrm in selected_chrms:
        plot_sp = chrmSzD[i-1][0]               # The start point of the selected genomic region in the plot
        plot_ep = chrmSzD[i-1][1]

        ch = df[(df.CHROM==chrm) & (df.POS>=plot_sp) & (df.POS<=plot_ep)]

        num_sv_on_chr.append([chrm, ch['sSV'].sum(), len(ch.index), ch['sSV'].sum()/len(ch.index)])

        sw_str = plot_sp                # The beginning of the window
        sw_end = plot_sp+sw_size        # The end of the window 
        icrs = incremental_step         # The incremental step
        x = []                          # Sliding window position on the chromosome
        y = []                          # Number of sSVs in the sliding window
        y_total = []                    # Number of total SVs in the sliding window
        y_ratio = []                    # sSV/totalSV ratio of the sliding window
        y5 = []                         # G-statistic
        y6 = []                         # The threshold of the G-statistic
        y7 = []                         # y7: Δ(allele frequency)
        y8, y9 = [], []                 # The confidence interval of the Δ(allele frequency)

        # chrm_size = plot_ep - plot_sp
        # adj_ratio = (chrm_size - sw_size) / chrm_size
        # adj_POS = plot_sp + (ch.POS - plot_sp) * adj_ratio

        # Calculate, sSV/totalSV, G-statistic, and Δ(allele frequency) of a sliding window
        while sw_end <= plot_ep:
            # sw_df: sSVs in a sliding window; sw_dfT: all SVs in a sliding window
            sw_df = ch[(ch.POS>=sw_str) & (ch.POS<=sw_end)]
            # a = sw_df[fb_ld].mean()
            # b = sw_df[sb_ld].mean()
            # c = sw_df['sSV'].sum()
            # d = sw_df['sSV'].mean()
            # e = sw_df['G_S'].sum()
            # f = sw_df['Delta_AF'].abs().sum()
            # g = sw_df['Delta_AF'].sum()

            row_in_sw_df = len(sw_df.index)         # number of SVs in a sliding window
            x.append(sw_str)                        # Append the start point of the sliding window to x
            y.append(sw_df['sSV'].sum().item())     # Append number of sSVs of the sliding window to y
            y_total.append(row_in_sw_df)            # Append number of totalSVs of the sliding window to y_total

            # 'try/exception' cannot catch the EmptyDataError; seems it only works when reading a .csv/.tsv file, not an empty subset of a existing dataframe. DivisionByZero generates 'nan' for a series or an array.
            if row_in_sw_df >= min_SVs:
                y_ratio.append(sw_df['sSV'].mean().item())    # Append the sSV/totalSV ratio of the sliding window to y_ratio
                y5.append(sw_df['G_S'].sum().item()/row_in_sw_df)
                if num_ipfiles == 1 and parent1 != 'ref':
                    y7.append(sw_df['Delta_AF'].abs().sum().item()/row_in_sw_df)
                else:
                    y7.append(sw_df['Delta_AF'].sum().item()/row_in_sw_df)

                if test == True and 'GS_Thrshld' in hdr:
                    y6.append(sw_df['GS_Thrshld'].sum().item()/row_in_sw_df)
                    y8.append(sw_df['DAF_CI_LB'].sum().item()/row_in_sw_df)
                    y9.append(sw_df['DAF_CI_UB'].sum().item()/row_in_sw_df)
                    # useful info of the sliding window
                    row_contents = [chrm, sw_str, sw_df[fb_ld].mean().item(), sw_df[sb_ld].mean().item(), sw_df['sSV'].sum().item(), row_in_sw_df, y_ratio[-1], y5[-1], y6[-1], y7[-1], y8[-1], y9[-1]]
                else:
                    row_contents = [chrm, sw_str, sw_df[fb_ld].mean().item(), sw_df[sb_ld].mean().item(), sw_df['sSV'].sum().item(), row_in_sw_df, y_ratio[-1], y5[-1], y7[-1]]
            else:
                wm_list.append(['No SV in the sliding window', i, sw_str])
                zeroSV(y_ratio)
                zeroSV(y5)
                zeroSV(y7)

                # The mean of an empty column is 'nan', not zero, and int(nan) generates a ValueError
                if test == True and 'GS_Thrshld' in hdr:
                    zeroSV(y6)
                    zeroSV(y8)
                    zeroSV(y9)
                    row_contents = [chrm, sw_str, 0, 0, sw_df['sSV'].sum().item(), row_in_sw_df, y_ratio[-1], y5[-1], y6[-1], y7[-1], y8[-1], y9[-1]]
                else:
                    row_contents = [chrm, sw_str, 0, 0, sw_df['sSV'].sum().item(), row_in_sw_df, y_ratio[-1], y5[-1], y7[-1]]

            sw_rows.append(row_contents)

            if i not in sw_dict:
                sw_dict[i] = [row_contents]
            else:
                sw_dict[i].append(row_contents)

            sw_str += icrs
            sw_end += icrs

        # Replace the 'empty' values at the beginning of the lists with nearest non-empty value
        for yl in [y_ratio, y5, y6, y7, y8, y9]:
            if 'empty' in yl:
                replace_zero(yl)

        # Find the first non-empty value in the sw_dict
        pIndex = 0
        while sw_dict[i][pIndex][6] == 'empty':
            pIndex += 1

        # Replace the empty value(s) with the above non-empty value in the sw_dict
        j = 0
        while j < pIndex:
            sw_dict[i][j][6] = sw_dict[i][pIndex][6]
            sw_dict[i][j][7] = sw_dict[i][pIndex][7]
            sw_dict[i][j][8] = sw_dict[i][pIndex][8]
            
            if test == True:
                sw_dict[i][j][9] = sw_dict[i][pIndex][9]
                sw_dict[i][j][10] = sw_dict[i][pIndex][10]
                sw_dict[i][j][11] = sw_dict[i][pIndex][11]

            j += 1

        # Smoothing data of a chromosome or a selected region
        sg_y_ratio = savgol_filter(y_ratio, smth_wl, poly_order)
        sg_y5 = savgol_filter(y5, smth_wl, poly_order)
        sg_y7 = savgol_filter(y7, smth_wl, poly_order)

        sg_y_ratio_list.extend(sg_y_ratio)
        sg_gs_list.extend(sg_y5)
        sg_daf_list.extend(sg_y7)

        if test == True and 'GS_Thrshld' in hdr:
            sg_y6 = savgol_filter(y6, smth_wl, poly_order)
            sg_y8 = savgol_filter(y8, smth_wl, poly_order)
            sg_y9 = savgol_filter(y9, smth_wl, poly_order)
            sg_t_gs_list.extend(sg_y6)
            sg_nt_daf_list.extend(sg_y8)
            sg_pt_daf_list.extend(sg_y9)

        # Handle the plot with a single column (chromosome)
        if len(selected_chrms) == 1:
            # Set up x-ticks
            axs[0].set_xticks(np.arange(plot_sp, x[-1], xt_pro[0]))
            ticks = axs[0].get_xticks()*xt_pro[1]
            axs[0].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0].set_ylabel('Number of SVs')
                axs[1].set_ylabel(r'sSV/totalSV')
                axs[2].set_ylabel('$G$-statistic')
                axs[3].set_ylabel('\u0394$AF$')

            # Plot sSVs and total SVs against their genomic positions
            axs[0].plot(x, y, c=curve_color)
            axs[0].plot(x, y_total, c=ttl_sv_color)
            if chrm.isdigit() == True:
                axs[0].set_title('Chr'+chrm)
            elif '.' in chrm:
                axs[0].set_title('Chr'+chrm_prnt_dict[chrm])
            else:
                axs[0].set_title(chrm)

            # Plot sSV/totalSV, G-statistic, and Δ(allele frequency) against their genomic positions
            if smoothing == True:
                # sSVs/totalSVs via Fisher's exact test
                axs[1].plot(x, sg_y_ratio, c=curve_color)
                # axs[1].scatter(adj_POS, ch.FE_P, marker=',', s=0.4, c=bg_color, zorder=-1)
                # axs[1].scatter(adj_POS, 1-ch.FE_P, marker=',', s=0.4, c=bg_color, zorder=-1)

                # G-statistic
                axs[2].plot(x, sg_y5, c=curve_color)
                # axs[2].scatter(adj_POS, ch.G_S, marker=',', s=0.4, c=bg_color, zorder=-1)

                # Δ(allele frequency)
                axs[3].plot(x, sg_y7, c=curve_color)
                #axs[3].scatter(adj_POS, ch.Delta_AF, marker=',', s=0.4, c=bg_color, zorder=-1)

                if test == True and 'GS_Thrshld' in hdr:
                    # G-statistic threshold at the SV level
                    axs[2].plot(x, sg_y6, c=sv_threshold_color)
                    # Δ(allele frequency) confidence interval at the SV level

                    if num_ipfiles == 1 and parent1 != 'ref':
                        axs[3].plot(x, sg_y9, c=sv_threshold_color)
                    elif num_ipfiles == 2 or parent1 == 'ref':
                        axs[3].plot(x, sg_y8, c=sv_threshold_color)
                        axs[3].plot(x, sg_y9, c=sv_threshold_color)
            else:
                # sSVs/totalSVs via Fisher's exact test
                axs[1].plot(x, y_ratio, c=curve_color)
                # axs[1].scatter(adj_POS, ch.FE_P, marker=',', s=0.4, c=bg_color, zorder=-1)
                # axs[1].scatter(adj_POS, 1-ch.FE_P, marker=',', s=0.4, c=bg_color, zorder=-1)

                # G-statistic
                axs[2].plot(x, y5, c=curve_color)
                # axs[2].scatter(adj_POS, ch.G_S, marker=',', s=0.4, c=bg_color, zorder=-1)

                # Δ(allele frequency)
                axs[3].plot(x, y7, c=curve_color)
                # axs[3].scatter(adj_POS, ch.Delta_AF, marker=',', s=0.4, c=bg_color, zorder=-1)

                if test == True and 'GS_Thrshld' in hdr:
                    # G-statistic threshold at the SV level
                    axs[2].plot(x, y6, c=sv_threshold_color)
                    # Δ(allele frequency) confidence interval at the SV level
                    if num_ipfiles == 1 and parent1 != 'ref':
                        axs[3].plot(x, sg_y9, c=sv_threshold_color)
                    elif num_ipfiles == 2 or parent1 == 'ref':
                        axs[3].plot(x, sg_y8, c=sv_threshold_color)
                        axs[3].plot(x, sg_y9, c=sv_threshold_color)

            # Add the 99.5 percentile line as the threshold, x[-1] is the startpoint of the last sliding window of a chromosome
            axs[1].plot([plot_sp, x[-1]], [thrshld_fe, thrshld_fe], c=sw_threshold_color)
            # axs[2].plot([plot_sp, x[-1]], [thrshld_gs, thrshld_gs], c=sw_threshold_color)
            axs[2].plot([plot_sp, x[-1]], [sv_ci_gs, sv_ci_gs], c=sw_threshold_color)
            # axs[3].plot([plot_sp, x[-1]], [thrshld_af, thrshld_af], c=sw_threshold_color)
            # axs[3].plot([plot_sp, x[-1]], [thrshld_af*(-1), thrshld_af*(-1)], c=sw_threshold_color)
            if num_ipfiles == 1 and parent1 != 'ref':
                # axs[3].plot([plot_sp, x[-1]], [thrshld_af_abs, thrshld_af_abs], c=sw_threshold_color)
                axs[3].plot([plot_sp, x[-1]], [sv_ci_af_abs, sv_ci_af_abs], c=sw_threshold_color)
            elif num_ipfiles == 2 or parent1 == 'ref':
                # axs[3].plot([plot_sp, x[-1]], [thrshld_af_n, thrshld_af_n], c=sw_threshold_color)
                # axs[3].plot([plot_sp, x[-1]], [thrshld_af_p, thrshld_af_p], c=sw_threshold_color)
                axs[3].plot([plot_sp, x[-1]], [sv_ci_af_lb, sv_ci_af_lb], c=sw_threshold_color)
                axs[3].plot([plot_sp, x[-1]], [sv_ci_af_ub, sv_ci_af_ub], c=sw_threshold_color)
            # axs[1].plot(x, sm_thresholds_sw, c='m')

        # Handle the plot with multiple columns (chromosomes)
        else:
            # Set up x-ticks
            axs[0,i-1].set_xticks(np.arange(plot_sp, x[-1], xt_pro[0]))
            ticks = axs[0,i-1].get_xticks()*xt_pro[1]
            axs[0,i-1].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0,i-1].set_ylabel('Number of SVs')
                axs[1,i-1].set_ylabel(r'sSV/totalSV')
                axs[2,i-1].set_ylabel('$G$-statistic')
                axs[3,i-1].set_ylabel('\u0394$AF$')

            # Plot sSV and totalSVs against their genomic positions
            axs[0,i-1].plot(x, y, c=curve_color)
            axs[0,i-1].plot(x, y_total, c=ttl_sv_color)
            if chrm.isdigit() == True:
                axs[0,i-1].set_title('Chr'+chrm)
            elif '.' in chrm:
                axs[0,i-1].set_title('Chr'+chrm_prnt_dict[chrm])
            else:
                axs[0,i-1].set_title(chrm)

            # Plot sSV/totalSV, G-statistic, and Δ(allele frequency) against their genomic positions
            if smoothing == True:
                # sSVs/totalSVs via Fisher's exact test
                axs[1,i-1].plot(x, sg_y_ratio, c=curve_color)
                # axs[1,i-1].scatter(adj_POS, ch.FE_P, marker=',', s=0.4, c=bg_color, zorder=-1)
                # axs[1,i-1].scatter(adj_POS, 1-ch.FE_P, marker=',', s=0.4, c=bg_color, zorder=-1)

                # G-statistic
                axs[2,i-1].plot(x, sg_y5, c=curve_color)
                # axs[2,i-1].scatter(adj_POS, ch.G_S, marker=',', s=0.4, c=bg_color, zorder=-1)

                # Δ(allele frequency)
                axs[3,i-1].plot(x, sg_y7, c=curve_color)
                # axs[3,i-1].scatter(adj_POS, ch.Delta_AF, marker=',', s=0.4, c=bg_color, zorder=-1)

                if test == True and 'GS_Thrshld' in hdr:
                    # G-statistic threshold at the SV level
                    axs[2,i-1].plot(x, sg_y6, c=sv_threshold_color)
                    # Δ(allele frequency) confidence interval at the SV level
                    if num_ipfiles == 1 and parent1 != 'ref':
                        axs[3,i-1].plot(x, sg_y9, c=sv_threshold_color)
                    elif num_ipfiles == 2 or parent1 == 'ref':
                        axs[3, i-1].plot(x, sg_y8, c=sv_threshold_color)
                        axs[3, i-1].plot(x, sg_y9, c=sv_threshold_color)
            else:
                # sSVs/totalSVs via Fisher's exact test
                axs[1,i-1].plot(x, y_ratio, c=curve_color)
                # axs[1,i-1].scatter(adj_POS, ch.FE_P, s=0.4, marker=',', c=bg_color, zorder=-1)
                # axs[1,i-1].scatter(adj_POS, 1-ch.FE_P, s=0.4, marker=',', c=bg_color, zorder=-1)
                
                # G-statistic
                axs[2,i-1].plot(x, y5, c=curve_color)
                # axs[2,i-1].scatter(adj_POS, ch.G_S, marker=',', s=0.4, c=bg_color, zorder=-1)

                # Δ(allele frequency)
                axs[3,i-1].plot(x, y7, c=curve_color)
                # axs[3,i-1].scatter(adj_POS, ch.Delta_AF, marker=',', s=0.4, c=bg_color, zorder=-1)

                if test == True and 'GS_Thrshld' in hdr:
                    # G-statistic threshold at the SV level
                    axs[2,i-1].plot(x, y6, c=sv_threshold_color)
                    # Δ(allele frequency) confidence interval at the SV level
                    if num_ipfiles == 1 and parent1 != 'ref':
                        axs[3,i-1].plot(x, sg_y9, c=sv_threshold_color)
                    elif num_ipfiles == 2 or parent1 == 'ref':
                        axs[3, i-1].plot(x, sg_y8, c=sv_threshold_color)
                        axs[3, i-1].plot(x, sg_y9, c=sv_threshold_color)
            # Add the 99.5 percentile line as the threshold, x[-1] is the startpoint of the last sliding window of a chromosome
            axs[1,i-1].plot([plot_sp, x[-1]], [thrshld_fe, thrshld_fe], c=sw_threshold_color)
            # axs[2,i-1].plot([plot_sp, x[-1]], [thrshld_gs, thrshld_gs], c=sw_threshold_color)
            axs[2,i-1].plot([plot_sp, x[-1]], [sv_ci_gs, sv_ci_gs], c=sw_threshold_color)
            # axs[3,i-1].plot([plot_sp, x[-1]], [thrshld_af, thrshld_af], c=sw_threshold_color)
            # axs[3,i-1].plot([plot_sp, x[-1]], [thrshld_af*(-1), thrshld_af*(-1)], c=sw_threshold_color)
            if num_ipfiles == 1 and parent1 != 'ref':
                # axs[3,i-1].plot([plot_sp, x[-1]], [thrshld_af_abs, thrshld_af_abs], c=sw_threshold_color)
                axs[3,i-1].plot([plot_sp, x[-1]], [sv_ci_af_abs, sv_ci_af_abs], c=sw_threshold_color)
            elif num_ipfiles == 2 or parent1 == 'ref':
                # axs[3,i-1].plot([plot_sp, x[-1]], [thrshld_af_n, thrshld_af_n], c=sw_threshold_color)
                # axs[3,i-1].plot([plot_sp, x[-1]], [thrshld_af_p, thrshld_af_p], c=sw_threshold_color)
                axs[3,i-1].plot([plot_sp, x[-1]], [sv_ci_af_lb, sv_ci_af_lb], c=sw_threshold_color)
                axs[3,i-1].plot([plot_sp, x[-1]], [sv_ci_af_ub, sv_ci_af_ub], c=sw_threshold_color)
            # axs[1, i-1].plot(x, sm_thresholds_sw, c='m')

        ratio_peak_list.append(max(y_ratio))

        # Identify genomic regions above the threshold
        zigzag = []     # List of peaks above the threshold
        # Handle the case in which an QTL is at the very beginning of the chromosome
        if sw_dict[i][0][6] >= thrshld_fe:
            # Add the CHROM ID and the startpoint to the genomic region above the threshold
            sv_region.append(sw_dict[i][0][:2])
            if sw_dict[i][0][6] >= sw_dict[i][1][6]:
                # Add the first peak above the threshold to zigzag
                zigzag.append(sw_dict[i][0][1:])
            num_of_sws = 1

        m = 0       # The index of the sliding window of a chromosome
        while m < len(sw_dict[i]) - 1:
            if sw_dict[i][m][6] < thrshld_fe and sw_dict[i][m+1][6] >= thrshld_fe:
                # Add the CHROM ID and the startpoint to the genomic region above the threshold
                sv_region.append(sw_dict[i][m+1][:2])
                num_of_sws = 1
            elif sw_dict[i][m][6] >= thrshld_fe:
                # A sliding window is considered as a peak if its sSV/totalSV is greater than or equal to the threshold and greater than those of the flanking sliding windows
                if m >= 1 and max(sw_dict[i][m-1][6], sw_dict[i][m+1][6]) <= sw_dict[i][m][6]:
                    # Add a peak above the threshold to zigzag
                    zigzag.append(sw_dict[i][m][1:])
                if sw_dict[i][m+1][6] > thrshld_fe:
                    num_of_sws += 1
                elif sw_dict[i][m+1][6] < thrshld_fe:
                    # Add the endpoint, all peaks, the number of sliding widows
                    sv_region[-1].extend([sw_dict[i][m][1], zigzag, num_of_sws])
                    zigzag = []
            m += 1

        # Handle the case in which a QTL is at the end of the chromosome
        if sw_dict[i][-1][6] >= thrshld_fe:
            sv_region[-1].extend([sw_dict[i][-1][1], zigzag, num_of_sws])

        i += 1

    n = 0
    for eml in sv_region:
        converted_peaks = []
        for peak in eml[3]:
            # Remove 'np.int' or 'np.float' from the list elements
            converted_peak = [x.item() if hasattr(x, 'item') else x for x in peak]
            converted_peaks.append(converted_peak)
        sv_region[n][3] = converted_peaks

        n += 1

    header_results = ['CHROM','Start','End','Peaks', 'NumOfSWs']
    pd.DataFrame(sv_region, columns=header_results).to_csv(sv_region_file, index=False)

    if test == True and 'GS_Thrshld' in hdr:
        sw_data_frame = pd.DataFrame(sw_rows, columns=['CHROM', 'POS', fb_id+'.AvgLD', sb_id+'.AvgLD', 'sSV', 'totalSV', r'sSV/totalSV', 'GS', 'GS_Thrshld', 'Delta_AF', 'DAF_CI_LB', 'DAF_CI_UB'])
        sw_data_frame['smthedGSThrshld'] = sg_t_gs_list
        sw_data_frame['smthedDAFNCI'] = sg_nt_daf_list
        sw_data_frame['smthedDAFPCI'] = sg_pt_daf_list
    else:
        sw_data_frame = pd.DataFrame(sw_rows, columns=['CHROM', 'POS', fb_id+'.AvgLD', sb_id+'.AvgLD', 'sSV', 'totalSV', r'sSV/totalSV', 'GS', 'Delta_AF'])

    sw_data_frame['smthedRatio'] = sg_y_ratio_list
    sw_data_frame['smthedGS'] = sg_gs_list
    sw_data_frame['smthedDAF'] = sg_daf_list
    sw_data_frame.to_csv(sw_file, index=False)

    misc.append(['List of the peaks of the chromosomes', ratio_peak_list])

    wrnLog = os.path.join(results, 'wrnLog.csv')
    with open(wrnLog, 'w', newline='') as out_f1:
        xie1 = csv.writer(out_f1)
        xie1.writerow(['Type', 'Chr', 'Position', 'Warning Message'])
        xie1.writerows(wm_list)

    num_sv_on_chr_file = os.path.join(results, 'num_sv_on_chr_file.csv')
    with open(num_sv_on_chr_file, 'w', newline='') as out_f2:
        xie2 = csv.writer(out_f2)
        xie2.writerow(['Chromosome', 'Num of sSVs', 'Num of totalSVs', r'sSV/totalSV'])
        xie2.writerows(num_sv_on_chr)

    # Handle the plot with a single column (chromosome)
    if len(selected_chrms) == 1:
        fig.align_ylabels(axs[:])
    # Handle the plot with multiple columns (chromosomes)
    else:
        fig.align_ylabels(axs[:, 0])

    # fig.tight_layout(pad=0.15, rect=[0, 0.035, 1, 1])     0.0264,0.054,0.076,0.002
    fig.subplots_adjust(top=t_margin, bottom=b_margin, left=l_margin, right=r_margin, hspace=h_gap, wspace=v_gap)
    fig.suptitle(f'Genomic position ({xt_pro[2]})', y=g_pos, ha='center', va='bottom')
    fig.text(0.001, 0.995, 'a', weight='bold', ha='left', va='top')
    fig.text(0.001, 0.680, 'b', weight='bold', ha='left', va='bottom')
    fig.text(0.001, 0.465, 'c', weight='bold', ha='left', va='bottom')
    fig.text(0.001, 0.244, 'd', weight='bold', ha='left', va='bottom')

    fig.savefig(os.path.join(results, 'PyBSASeq.pdf'))
    fig.savefig(os.path.join(results, 'PyBSASeq.eps'))
    fig.savefig(os.path.join(results, 'PyBSASeq.svg'))
    fig.savefig(os.path.join(results, 'PyBSASeq.png'), dpi=600)
    print(f'Plotting completed, time elapsed: {(time.time()-t0)/60} minutes.')


def bsaseq_plot_sw(df):
    print('Plotting via the sliding window algorithm')
    with open(thrshld_file, 'r') as du:
        threshold_list = du.readline().strip().split(' ')
        temp = [float(i) for i in threshold_list]
        thrshld_fe = temp[0]
        thrshld_gs = temp[1]
        thrshld_af_n= temp[2]
        thrshld_af_p = temp[3]
        thrshld_af_abs = temp[4]
        sv_ci_af_lb = temp[5]
        sv_ci_af_ub = temp[6]
        sv_ci_gs = temp[7]
        sv_ci_af_abs = temp[8]

    tl = [i for i in selected_chrms]
    for x in tl:
        y = df[df.CHROM==x]
        # print(x, len(y.index))
        if len(y.index) < min_sv:
            selected_chrms.remove(x)

    # print(selected_chrms)
    if selected_chrms == []:
        print('Your dataset contains too few SVs. No analysis is performed')
        sys.exit()

    # Plot layout setup
    height_ratio = [1,0.8,0.8,0.8]
    xt_pro = xticks_property(chrmSzL)
    fig, axs = plt.subplots(nrows=len(height_ratio), ncols=len(selected_chrms), figsize=(20, 12.8), sharex='col', sharey='row', gridspec_kw={'width_ratios': chrmSzL, 'height_ratios': height_ratio})

    hdr = df.columns.values.tolist()
    # Analyze each chromosome separately
    i = 1
    for chrm in selected_chrms:
        plot_sp = chrmSzD[i-1][0]               # The start point of the selected genomic region in the plot
        plot_ep = chrmSzD[i-1][1]

        ch = df[(df.CHROM==chrm) & (df.POS>=plot_sp) & (df.POS<=plot_ep)]
        # ch = df[df.CHROM==chrm]

        # Handle the plot with a single column (chromosome)
        if len(selected_chrms) == 1:
            # Set up x-ticks
            axs[0].set_xticks(np.arange(plot_sp, plot_ep, xt_pro[0]))
            ticks = axs[0].get_xticks()*xt_pro[1]
            axs[0].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0].set_ylabel('Number of SVs')
                axs[1].set_ylabel(r'sSV/totalSV')
                axs[2].set_ylabel('$G$-statistic')
                axs[3].set_ylabel('\u0394$AF$')

            # Plot sSVs and total SVs against their genomic positions
            axs[0].plot(ch.POS, ch.totalSV, c=ttl_sv_color)
            axs[0].plot(ch.POS, ch.sSV, c=curve_color)

            if chrm.isdigit() == True:
                axs[0].set_title('Chr'+chrm)
            elif '.' in chrm:
                axs[0].set_title('Chr'+chrm_prnt_dict[chrm])
            else:
                axs[0].set_title(chrm)

            # Plot sSV/totalSV, G-statistic, and Δ(allele frequency) against their genomic positions
            if smoothing == True:
                # sSVs/totalSVs via Fisher's exact test
                axs[1].plot(ch.POS, ch.smthedRatio, c=curve_color)
                # G-statistic
                axs[2].plot(ch.POS, ch.smthedGS, c=curve_color)
                # Δ(allele frequency)
                axs[3].plot(ch.POS, ch.smthedDAF, c=curve_color)
                if test == True and 'smthedGSThrshld' in hdr:
                    # G-statistic threshold at the SV level
                    axs[2].plot(ch.POS, ch.smthedGSThrshld, c=sv_threshold_color)
                    # Δ(allele frequency) confidence intervals at the SV level
                    if num_ipfiles == 1 and parent1 != 'ref':
                        axs[3].plot(ch.POS, ch.smthedDAFPCI, c=sv_threshold_color)
                    elif num_ipfiles == 2 or parent1 == 'ref':
                        axs[3].plot(ch.POS, ch.smthedDAFPCI, c=sv_threshold_color)
                        axs[3].plot(ch.POS, ch.smthedDAFNCI, c=sv_threshold_color)
            else:
                 # sSVs/totalSVs via Fisher's exact test
                axs[1].plot(ch.POS, ch[r'sSV/totalSV'], c=curve_color)
                # G-statistic
                axs[2].plot(ch.POS, ch.GS, c=curve_color,)
                # Δ(allele frequency)
                axs[3].plot(ch.POS, ch.Delta_AF, c=curve_color,)
                if test == True and 'smthedGSThrshld' in hdr:
                    # G-statistic threshold at the SV level
                    axs[2].plot(ch.POS, ch.GS_Thrshld, c=sv_threshold_color)
                    # Δ(allele frequency) confidence intervals at the SV level
                    if num_ipfiles == 1 and parent1 != 'ref':
                        axs[3].plot(ch.POS, ch.DAF_CI_UB, c=sv_threshold_color)
                    elif num_ipfiles == 2 or parent1 == 'ref': 
                        axs[3].plot(ch.POS, ch.DAF_CI_LB, c=sv_threshold_color)
                        axs[3].plot(ch.POS, ch.DAF_CI_UB, c=sv_threshold_color)

            # Add the 99.5 percentile line as the threshold, plot_ep is the startpoint of the last sliding window of a chromosome
            axs[1].plot([plot_sp, plot_ep], [thrshld_fe, thrshld_fe], c=sw_threshold_color)
            # axs[2].plot([plot_sp, plot_ep], [thrshld_gs, thrshld_gs], c=sw_threshold_color)
            axs[2].plot([plot_sp, plot_ep], [sv_ci_gs, sv_ci_gs], c=sw_threshold_color)
            # axs[3].plot([plot_sp, plot_ep], [thrshld_af, thrshld_af], c=sw_threshold_color)
            # axs[3].plot([plot_sp, plot_ep], [thrshld_af*(-1), thrshld_af*(-1)], c=sw_threshold_color)
            if num_ipfiles == 1 and parent1 != 'ref':
                # axs[3].plot([plot_sp, plot_ep], [thrshld_af_abs, thrshld_af_abs], c=sw_threshold_color)
                axs[3].plot([plot_sp, plot_ep], [sv_ci_af_abs, sv_ci_af_abs], c=sw_threshold_color)
            elif num_ipfiles == 2 or parent1 == 'ref':
                # axs[3].plot([plot_sp, plot_ep], [thrshld_af_n, thrshld_af_n], c=sw_threshold_color)
                # axs[3].plot([plot_sp, plot_ep], [thrshld_af_p, thrshld_af_p], c=sw_threshold_color)
                axs[3].plot([plot_sp, plot_ep], [sv_ci_af_lb, sv_ci_af_lb], c=sw_threshold_color)
                axs[3].plot([plot_sp, plot_ep], [sv_ci_af_ub, sv_ci_af_ub], c=sw_threshold_color)
            # axs[1].plot(x, sm_thresholds_sw, c='m')

        # Handle the plot with multiple columns (chromosomes)
        else:
            # Set up x-ticks
            axs[0,i-1].set_xticks(np.arange(plot_sp, plot_ep, xt_pro[0]))
            ticks = axs[0,i-1].get_xticks()*xt_pro[1]
            axs[0,i-1].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0,i-1].set_ylabel('Number of SVs')
                axs[1,i-1].set_ylabel(r'sSV/totalSV')
                axs[2,i-1].set_ylabel('$G$-statistic')
                axs[3,i-1].set_ylabel('\u0394$AF$')

            # Plot sSV and totalSVs against their genomic positions
            axs[0,i-1].plot(ch.POS, ch.totalSV, c=ttl_sv_color)
            axs[0,i-1].plot(ch.POS, ch.sSV, c=curve_color)
            if chrm.isdigit() == True:
                axs[0,i-1].set_title('Chr'+chrm)
            elif '.' in chrm:
                axs[0,i-1].set_title('Chr'+chrm_prnt_dict[chrm])
            else:
                axs[0,i-1].set_title(chrm)

            # Plot sSV/totalSV, G-statistic, and Δ(allele frequency) against their genomic positions
            if smoothing == True:
                # sSVs/totalSVs via Fisher's exact test
                axs[1,i-1].plot(ch.POS, ch.smthedRatio, c=curve_color)
                # G-statistic
                axs[2,i-1].plot(ch.POS, ch.smthedGS, c=curve_color)
                # Δ(allele frequency)
                axs[3,i-1].plot(ch.POS, ch.smthedDAF, c=curve_color)
                if test == True and 'smthedGSThrshld' in hdr:
                    # G-statistic threshold at the SV level
                    axs[2,i-1].plot(ch.POS, ch.smthedGSThrshld, c=sv_threshold_color)
                    # Δ(allele frequency) confidence intervals at the SV level
                    if num_ipfiles == 1 and parent1 != 'ref':
                        axs[3,i-1].plot(ch.POS, ch.smthedDAFPCI, c=sv_threshold_color)
                    elif num_ipfiles == 2 or parent1 == 'ref':
                        axs[3,i-1].plot(ch.POS, ch.smthedDAFPCI, c=sv_threshold_color)
                        axs[3,i-1].plot(ch.POS, ch.smthedDAFNCI, c=sv_threshold_color)
            else:
                # sSVs/totalSVs via Fisher's exact test
                axs[1,i-1].plot(ch.POS, ch[r'sSV/totalSV'], c=curve_color)
                # G-statistic
                axs[2,i-1].plot(ch.POS, ch.GS, c=curve_color,)
                # Δ(allele frequency)
                axs[3,i-1].plot(ch.POS, ch.Delta_AF, c=curve_color,)
                if test == True and 'smthedGSThrshld' in hdr:
                    # G-statistic threshold at the SV level
                    axs[2,i-1].plot(ch.POS, ch.GS_Thrshld, c=sv_threshold_color)
                    # Δ(allele frequency) confidence intervals at the SV level
                    if num_ipfiles == 1 and parent1 != 'ref':
                        axs[3,i-1].plot(ch.POS, ch.DAF_CI_UB, c=sv_threshold_color)
                    elif num_ipfiles == 2 or parent1 == 'ref':
                        axs[3,i-1].plot(ch.POS, ch.DAF_CI_LB, c=sv_threshold_color)
                        axs[3,i-1].plot(ch.POS, ch.DAF_CI_UB, c=sv_threshold_color)

            # Add the 99.5 percentile line as the threshold, plot_ep is the startpoint of the last sliding window of a chromosome
            axs[1,i-1].plot([plot_sp, plot_ep], [thrshld_fe, thrshld_fe], c=sw_threshold_color)
            # axs[2,i-1].plot([plot_sp, plot_ep], [thrshld_gs, thrshld_gs], c=sw_threshold_color)
            axs[2,i-1].plot([plot_sp, plot_ep], [sv_ci_gs, sv_ci_gs], c=sw_threshold_color)
            # axs[3,i-1].plot([plot_sp, plot_ep], [thrshld_af, thrshld_af], c=sw_threshold_color)
            # axs[3,i-1].plot([plot_sp, plot_ep], [thrshld_af*(-1), thrshld_af*(-1)], c=sw_threshold_color)
            if num_ipfiles == 1 and parent1 != 'ref':
                # axs[3,i-1].plot([plot_sp, plot_ep], [thrshld_af_abs, thrshld_af_abs], c=sw_threshold_color)
                axs[3,i-1].plot([plot_sp, plot_ep], [sv_ci_af_abs, sv_ci_af_abs], c=sw_threshold_color)
            elif num_ipfiles == 2 or parent1 == 'ref':
                # axs[3,i-1].plot([plot_sp, plot_ep], [thrshld_af_n, thrshld_af_n], c=sw_threshold_color)
                # axs[3,i-1].plot([plot_sp, plot_ep], [thrshld_af_p, thrshld_af_p], c=sw_threshold_color)
                axs[3,i-1].plot([plot_sp, plot_ep], [sv_ci_af_lb, sv_ci_af_lb], c=sw_threshold_color)
                axs[3,i-1].plot([plot_sp, plot_ep], [sv_ci_af_ub, sv_ci_af_ub], c=sw_threshold_color)
            # axs[1, i-1].plot(x, sm_thresholds_sw, c='m')

        i += 1     
        
    # Handle the plot with a single column (chromosome)
    if len(selected_chrms) == 1:
        fig.align_ylabels(axs[:])
    # Handle the plot with multiple columns (chromosomes)
    else:
        fig.align_ylabels(axs[:, 0])

    # fig.tight_layout(pad=0.15, rect=[0, 0.035, 1, 1])     0.0264,0.054,0.076,0.002
    fig.subplots_adjust(top=t_margin, bottom=b_margin, left=l_margin, right=r_margin, hspace=h_gap, wspace=v_gap)
    fig.suptitle(f'Genomic position ({xt_pro[2]})', y=g_pos, ha='center', va='bottom')
    fig.text(0.001, 0.995, 'a', weight='bold', ha='left', va='top')
    fig.text(0.001, 0.680, 'b', weight='bold', ha='left', va='bottom')
    fig.text(0.001, 0.465, 'c', weight='bold', ha='left', va='bottom')
    fig.text(0.001, 0.244, 'd', weight='bold', ha='left', va='bottom')

    fig.savefig(os.path.join(results, 'PyBSASeq.pdf'))
    fig.savefig(os.path.join(results, 'PyBSASeq.eps'))
    fig.savefig(os.path.join(results, 'PyBSASeq.svg'))
    fig.savefig(os.path.join(results, 'PyBSASeq.png'), dpi=600)
    print(f'Plotting completed, time elapsed: {(time.time()-t0)/60} minutes.')


def peak(l):
    # print(type(l))
    ratio_list = []
    for sub_l in l:
        # print(sub_l)
        ratio_list.append(sub_l[5])

    return l[ratio_list.index(max(ratio_list))]


def pk_list(l):
    peak_list = []
    for sub_l in l:
        chromosome = sub_l[0]
        if sub_l[3] != [] and sub_l[4] > 10:
            tempL = peak(sub_l[3])
            if not isinstance(chromosome, str):         # isinstance(chromosome, int)
                peak_list.append([str(chromosome), tempL[0], tempL[5]])
            else:
               peak_list.append([chromosome, tempL[0], tempL[5]]) 

    return peak_list


def accurate_threshold_sw(l, df):
    peaks = []
    for sub_l in l:
        peak_sw = df[(df.CHROM == sub_l[0]) & (df.POS >= sub_l[1]) & (df.POS <= sub_l[1]+sw_size-1)]
        __, pvalue_tt = ttest_rel(peak_sw[fb_af], peak_sw[sb_af])
        sSV_peak_sw = peak_sw[peak_sw.FE_P<alpha]

        sSV, totalSV = len(sSV_peak_sw.index), len(peak_sw.index)
        try:
            ratio = sSV / totalSV
            thresholds = thresholds_sw(peak_sw)

            # peaks.append([sub_l[0], sub_l[1], int(peak_sw[fb_ld].mean()), int(peak_sw[sb_ld].mean()), sSV, totalSV, ratio, thresholds_sw(peak_sw)])
            peaks.append([sub_l[0], sub_l[1], peak_sw[fb_ld].mean(), peak_sw[sb_ld].mean(), sSV, totalSV, ratio, thresholds[0], peak_sw.G_S.mean(), thresholds[1], peak_sw.Delta_AF.mean(), thresholds[2][0], thresholds[2][1], thresholds[3], pvalue_tt])
        except ZeroDivisionError:
            print('This sliding window does not contain any SV')

    header_results = ['CHROM','POS', fb_id+'.AvgLD', sb_id+'.AvgLD', 'sSV', 'totalSV', r'sSV/totalSV', 'Threshold_sSV', 'GS', 'Threshold_GS', 'DAF', 'DAF_CI_LB', 'DAF_CI_UB', 'Threshold_DAF', 'pvalue_tt']
    peak_df = pd.DataFrame(peaks, columns=header_results)
    peak_df['Significance_sSV'] = np.where(peak_df[r'sSV/totalSV']>peak_df['Threshold_sSV'], 1, 0)
    peak_df['Significance_GS'] = np.where(peak_df.GS>peak_df.Threshold_GS, 1, 0)
    if num_ipfiles == 1 and parent1 != 'ref':
        peak_df['Significance_AF'] = np.where((peak_df.DAF>peak_df.Threshold_DAF), 1, 0)
    elif num_ipfiles == 2 or parent1 == 'ref':
        peak_df['Significance_AF'] = np.where((peak_df.DAF>peak_df.DAF_CI_UB) | (peak_df.DAF<peak_df.DAF_CI_LB), 1, 0)
    peak_df['Significance_TT'] = np.where(peak_df['pvalue_tt']<alpha, 1, 0)

    peak_df.to_csv(os.path.join(results, args.output), index=False)


def accurate_threshold_gw(l, df):
    peaks = []
    for sub_l in l:
        peak_sw = df[(df.CHROM == sub_l[0]) & (df.POS >= sub_l[1]) & (df.POS <= sub_l[1]+sw_size-1)]
        sSV_peak_sw = peak_sw[peak_sw.FE_P<alpha]
        sSV, totalSV = len(sSV_peak_sw.index), len(peak_sw.index)
        ratio = sSV / totalSV

        # peaks.append([sub_l[0], sub_l[1], int(peak_sw[fb_ld].mean()), int(peak_sw[sb_ld].mean()), sSV, totalSV, ratio, thrshld_fe])
        peaks.append([sub_l[0], sub_l[1], peak_sw[fb_ld].mean(), peak_sw[sb_ld].mean(), sSV, totalSV, ratio, thrshld_fe])

    header_results = ['CHROM','POS', fb_id+'.AvgLD', sb_id+'.AvgLD', 'sSV', 'totalSV', r'sSV/totalSV', 'Threshold']
    peak_df = pd.DataFrame(peaks, columns=header_results)
    peak_df.to_csv(os.path.join(results, args.output), index=False)


def impute_di(df):
    df['DI'] = di(df)
    df.to_csv(os.path.join(dgns_path, 'bulk_df_di_before.csv'))

    df.loc[df.DI==0, 'Delta_AF'] = np.nan
    df.loc[df.DI==0, 'G_S'] = np.nan
    df.loc[df.DI==0, 'FE_P'] = np.nan

    df['Delta_AF'].ffill(inplace=True)
    df['G_S'].ffill(inplace=True)
    df['FE_P'].ffill(inplace=True)

    # Calculate allele frequency
    # df[fb_af] = df[fb_ad_alt] / df[fb_ld]
    # df[sb_af] = df[sb_ad_alt] / df[sb_ld]

    # df.loc[df.DI==0, fb_af] = np.nan
    # df.loc[df.DI==0, sb_af] = np.nan
    # 
    # df[fb_af].ffill(inplace=True)
    # df[sb_af].ffill(inplace=True)
    # for chrm in selected_chrms:
    #     df.loc[df['CHROM']==chrm][fb_af].interpolate('nearest', inplace=True)
    #     df.loc[df['CHROM']==chrm][sb_af].interpolate('nearest', inplace=True)

    # d_len = len(df.loc[df.DI==0].index)
    # d_fb_ld_mean = df.loc[df.DI==0][fb_ld].mean()
    # d_fb_ld_std = df.loc[df.DI==0][fb_ld].std()
    # d_sb_ld_mean = df.loc[df.DI==0][sb_ld].mean()
    # d_sb_ld_std = df.loc[df.DI==0][sb_ld].std()
    # 
    # df.loc[df.DI==0, fb_ld] = np.rint(np.random.normal(d_fb_ld_mean, d_fb_ld_std, d_len))
    # df.loc[df.DI==0, sb_ld] = np.rint(np.random.normal(d_sb_ld_mean, d_sb_ld_std, d_len))
    # 
    # df = read_filter(df)
    # 
    # df.loc[df.DI==0, fb_ad_alt] = (df.loc[df.DI==0][fb_af] * df.loc[df.DI==0][fb_ld]).round().astype(int)
    # df.loc[df.DI==0, fb_ad_ref] = df.loc[df.DI==0][fb_ld] - df.loc[df.DI==0][fb_ad_alt]
    # df.loc[df.DI==0, sb_ad_alt] = (df.loc[df.DI==0][sb_af] * df.loc[df.DI==0][sb_ld]).round().astype(int)
    # df.loc[df.DI==0, sb_ad_ref] = df.loc[df.DI==0][sb_ld] - df.loc[df.DI==0][sb_ad_alt]

    # bsa_di = df[df.DI==1]
    # bsa_d = df[df.DI==0]

    # d_len = len(bsa_d.index)

    # d_fb_ld_mean = bsa_d[fb_ld].mean()
    # d_fb_ld_std = bsa_d[fb_ld].std()
    # d_sb_ld_mean = bsa_d[sb_ld].mean()
    # d_sb_ld_std = bsa_d[sb_ld].std()
            
    # bsa_d[fb_ld] = np.rint(np.random.normal(d_fb_ld_mean, d_fb_ld_std, d_len))
    # bsa_d[sb_ld] = np.rint(np.random.normal(d_sb_ld_mean, d_sb_ld_std, d_len))

    '''
    Got the warning below using the above code:
    SettingWithCopyWarning: A value is trying to be set on a copy of a slice from a DataFrame. Try using .loc[row_indexer,col_indexer] = value instead
    '''

    # bsa_d = read_filter(bsa_d)

    # bsa_d[fb_ld] = np.random.normal(d_fb_ld_mean, d_fb_ld_std, d_len).round().astype(int)
    # bsa_d[sb_ld] = np.random.normal(d_sb_ld_mean, d_sb_ld_std, d_len).round().astype(int)
    # bsa_d[fb_ad_alt] = (bsa_d[fb_af] * bsa_d[fb_ld]).round().astype(int)
    # bsa_d[fb_ad_ref] = bsa_d[fb_ld] - bsa_d[fb_ad_alt]
    # bsa_d[sb_ad_alt] = (bsa_d[sb_af] * bsa_d[sb_ld]).round().astype(int)
    # bsa_d[sb_ad_ref] = bsa_d[sb_ld] - bsa_d[sb_ad_alt]
            
    # df = pd.concat([bsa_di, bsa_d]).sort_index()
    df.to_csv(os.path.join(dgns_path, 'bulk_df_di_after.csv'))

    return df


def calculate_statistics(df):
    # Calculate simulated ALT reads for each SV under null hypothesis
    df[sm_fb_ad_alt] = np.random.binomial(df[fb_ld], fb_freq)
    df[sm_fb_ad_ref] = df[fb_ld] - df[sm_fb_ad_alt]
    df[sm_sb_ad_alt] = np.random.binomial(df[sb_ld], sb_freq)
    df[sm_sb_ad_ref] = df[sb_ld] - df[sm_sb_ad_alt]

    # Calculate allele frequency
    df[fb_af] = df[fb_ad_alt] / df[fb_ld]
    df[sb_af] = df[sb_ad_alt] / df[sb_ld]
    df['Delta_AF'] = df[sb_af] - df[fb_af]

    if num_ipfiles == 1 and parent1 != 'ref':
        df['Delta_AF'] = df['Delta_AF'].abs()

    # Calculate G-statistic
    df['G_S'] = g_statistic_array(df[fb_ad_ref], df[fb_ad_alt], df[sb_ad_ref], df[sb_ad_alt])

    fb_af_ci, sb_af_ci = fb_id + '.AF_CI', sb_id + '.AF_CI'
    # fb_af_ci_ub, sb_af_ci_ub = fb_af_ci+'_ub', sb_af_ci+'_ub'
    # fb_af_ci_lb, sb_af_ci_lb = fb_af_ci+'_lb', sb_af_ci+'_lb'

    if 'fisher' in sys.modules:
        # Create new columns for Fisher's exact test P-values and simulated P-values
        print('Perform Fisher\'s exact test on each SV.')
        fb_ad_alt_arr = df[fb_ad_alt].to_numpy(dtype=np.uint32)
        fb_ad_ref_arr = df[fb_ad_ref].to_numpy(dtype=np.uint32)
        sb_ad_alt_arr = df[sb_ad_alt].to_numpy(dtype=np.uint32)
        sb_ad_ref_arr = df[sb_ad_ref].to_numpy(dtype=np.uint32)

        __, __, df['FE_P'] = pvalue_npy(fb_ad_alt_arr, fb_ad_ref_arr, sb_ad_alt_arr, sb_ad_ref_arr)
        # df['FE_OR'] = (fb_ad_alt_arr * sb_ad_ref_arr) / (fb_ad_ref_arr * sb_ad_alt_arr)

        sm_fb_ad_alt_arr = df[sm_fb_ad_alt].to_numpy(dtype=np.uint32)
        sm_fb_ad_ref_arr = df[sm_fb_ad_ref].to_numpy(dtype=np.uint32)
        sm_sb_ad_alt_arr = df[sm_sb_ad_alt].to_numpy(dtype=np.uint32)
        sm_sb_ad_ref_arr = df[sm_sb_ad_ref].to_numpy(dtype=np.uint32)

        __, __, df['sm_FE_P'] = pvalue_npy(sm_fb_ad_alt_arr, sm_fb_ad_ref_arr, sm_sb_ad_alt_arr, sm_sb_ad_ref_arr)
        # df['sm_FE_OR'] = (sm_fb_ad_alt_arr * sm_sb_ad_ref_arr) / (sm_fb_ad_ref_arr * sm_sb_ad_alt_arr)

        print(f'Fisher\'s exact test completed, time elapsed: {(time.time()-t0)/60} minutes.')

        if test == True:
            print('Calculate the \u0394(allele frequency) and G-statistic thresholds of each SV.')
            df['STAT'] = df.apply(sv_ci, axis=1)

            # Create new columns containing allele frequency, Δ(allele frequency) confidence intervals, and G-statistic thresholds
            # df[[fb_af_ci, sb_af_ci, 'DAF_CI', 'GS_CI']] = pd.DataFrame(df.STAT.values.tolist(), index=df.index)
            df[['DAF_CI', 'GS_Thrshld', 'DAF_abs_Thrshld']] = pd.DataFrame(df.STAT.values.tolist(), index=df.index)

            print(f'Calculating thresholds of \u0394(allele frequency) and G-statistic completed, time elapsed: {(time.time()-t0)/60} minutes.')

    elif 'scipy' in sys.modules:
        print('The module \'Fisher\' (https://github.com/brentp/fishers_exact_test) is not installed on your computer, \'fisher_exact\' from \'scipy.stats\' will be used instead, but the former is much faster in processing large datasets.')
        print('Perform Fisher\'s exact test on each SV')
        df['STAT'] = df.apply(sv_ci_scipy, axis=1)
        df.to_csv('output.csv', columns=['CHROM','STAT'],index=False)

        if test == True:
            # Create new columns for Fisher's exact test results, allele frequency, Δ(allele frequency) confidence intervals, and G-statistic thresholds
            # df[['fisher_exact', 'sm_FE', fb_af_ci, sb_af_ci, 'DAF_CI', 'GS_CI']] = pd.DataFrame(df.STAT.values.tolist(), index=df.index)
            # df[['fisher_exact', 'sm_FE', 'DAF_CI', 'GS_CI']] = pd.DataFrame(df.STAT.values.tolist(), index=df.index)
            df[['FE_P', 'sm_FE_P', 'DAF_CI', 'GS_Thrshld', 'DAF_abs_Thrshld']] = pd.DataFrame(df.STAT.values.tolist(), index=df.index)
        else:
            # df[['fisher_exact', 'sm_FE']] = pd.DataFrame(df.STAT.values.tolist(), index=df.index)
            df[['FE_P', 'sm_FE_P']] = pd.DataFrame(df.STAT.values.tolist(), index=df.index)
        
        # Create new columns for Fisher's exact test P-values or simulated P-values
        # df['FE_P'] = df['fisher_exact'].apply(lambda x: x[1]).astype(float)
        # df['sm_FE_P'] = df['sm_FE'].apply(lambda x: x[1]).astype(float)

        print(f'Fisher\'s exact test and calculating thresholds of \u0394(allele frequency) and G-statistic completed, time elapsed: {(time.time()-t0)/60} minutes.')

    if test == True:
        # Create new columns for 99% Δ(allele frequency) confidence intervals, and 99.5 percentile G-statistic thresholds
        # df[fb_af_ci_lb] = df[fb_af_ci].apply(lambda x: x[0]).astype(float)
        # df[fb_af_ci_ub] = df[fb_af_ci].apply(lambda x: x[1]).astype(float)
        # df[sb_af_ci_lb] = df[sb_af_ci].apply(lambda x: x[0]).astype(float)            # df[sb_af_ci_ub] = df[sb_af_ci].apply(lambda x: x[1]).astype(float)
        df['DAF_CI_LB'] = df['DAF_CI'].apply(lambda x: x[0]).astype(float)
        df['DAF_CI_UB'] = df['DAF_CI'].apply(lambda x: x[1]).astype(float)
        # df['GS_Thrshld'] = df['GS_CI'].apply(lambda x: x[1]).astype(float)

        # falseSVs = df[((df[fb_af] < df[fb_af_ci_lb]) & (df[sb_af] < df[sb_af_ci_lb])) | ((df[fb_af] > df[fb_af_ci_ub]) & (df[sb_af] > df[sb_af_ci_ub]))]

        # df = df.drop(index=falseSVs.index)

        # falseSVs.to_csv(os.path.join(filtering_path, 'falseSVs.csv'), index=None)
        # df.to_csv(os.path.join(filtering_path, 'SVs.csv'), index=None)

    # To make the mirror image of the Δ(allele frequency) curve for comparison if necessary
    if mirror_index == -1:
        df['Delta_AF'] = df['Delta_AF'] * -1
        if test == True:
            df['DAF_CI_UB'] = df['DAF_CI_UB'] * -1
            df['DAF_CI_LB'] = df['DAF_CI_LB'] * -1

    # The above calculation may generate 'NA' value(s) for some SVs. Remove SVs with such 'NA' value(s)
    df.dropna(inplace=True)
    misc.append(['Number of SVs after drop of SVs with calculation-generated NA value', len(df.index)])

    misc.append(['Dataframe filtered with genotype quality scores', len(df.index)])

    # Make the saved file smaller via eliminating some fields that can be easily calculated.
    # Reorganize the columns.
    if test == True:
        reorder_columns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', fb_gt, fb_ad, fb_ad_ref, fb_ad_alt, fb_ld, sm_fb_ad_alt, fb_af, fb_gq, sb_gt, sb_ad, sb_ad_ref, sb_ad_alt, sb_ld, sm_sb_ad_alt, sb_af, sb_gq, 'Delta_AF', 'FE_P', 'sm_FE_P', 'G_S', 'DAF_CI_LB', 'DAF_CI_UB', 'GS_Thrshld', 'DAF_abs_Thrshld', 'DAF_CI']
        # reorder_columns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', fb_gt, fb_ad, fb_ad_ref, fb_ad_alt, fb_ld, sm_fb_ad_alt, fb_af, fb_af_ci_lb, fb_af_ci_ub, fb_af_ci, fb_gq, sb_gt, sb_ad, sb_ad_ref, sb_ad_alt, sb_ld, sm_sb_ad_alt, sb_af, sb_af_ci_lb, sb_af_ci_ub, sb_af_ci, sb_gq, 'Delta_AF', 'FE_P', 'sm_FE_P', 'G_S', 'DAF_CI_LB', 'DAF_CI_UB', 'GS_Thrshld', 'DAF_CI', 'GS_CI', 'STAT']
    else:
        reorder_columns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', fb_gt, fb_ad, fb_ad_ref, fb_ad_alt, fb_ld, sm_fb_ad_alt, fb_af, fb_gq, sb_gt, sb_ad, sb_ad_ref, sb_ad_alt, sb_ld, sm_sb_ad_alt, sb_af, sb_gq, 'Delta_AF', 'FE_P', 'sm_FE_P', 'G_S']
    # Remove unnecessary columns and reorganize the columns
    # reorder_columns = ['CHROM', 'POS', 'REF', 'ALT', fb_gt, fb_ad_ref, fb_ad_alt, fb_ld, sm_fb_ad_alt, fb_gq, sb_gt, sb_ad_ref, sb_ad_alt, sb_ld, sm_sb_ad_alt, sb_gq, 'FE_P', 'sm_FE_P']

    df = df[reorder_columns]
    df_fep = df[df.FE_P>0.8]
    print('Save the original and the calculated values to a file.')
    df.to_csv(oi_file, index=None)
    # df_fep.to_csv(fep_file, index=None)
    print(f'File writing completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return df


def gw_thresholds(df):
    if 'fisher' in sys.modules:
        thresholds = thresholds_gw(df)
    elif 'scipy' in sys.modules:
        thresholds = thresholds_gw_approximate(df)

    return thresholds


def peak_verification(l, df):
    if 'fisher' in sys.modules:
        print('Verify potential significant peaks/valleys.')
        accurate_threshold_sw(l, df)
        print(f'Peak verification completed, time elapsed: {(time.time()-t0)/60} minutes.')
    elif 'scipy' in sys.modules:
        accurate_threshold_gw(l, df)
        print('Please install the module \'Fisher\' if more precise thresholds of the QTL loci are desired.')


def seg_statistics(df):
    # num_cols = [fb_ad_ref, fb_ad_alt, sb_ad_ref, sb_ad_alt, fb_ld, sb_ld]
    # df['Seg'] = df.POS // seg_size
    group_df = df.groupby(['CHROM', 'Seg'])[num_cols].mean().reset_index()
    
    group_df['Delta_AF'] = group_df[sb_ad_alt] / group_df[sb_ld] - group_df[fb_ad_alt] / group_df[fb_ld]
    group_df['G_S'] = g_statistic_array(group_df[fb_ad_ref], group_df[fb_ad_alt], group_df[sb_ad_ref], group_df[sb_ad_alt])
    __, __, group_df['FE_P'] = pvalue_npy(group_df[fb_ad_ref].values.astype(np.uint32), group_df[fb_ad_alt].values.astype(np.uint32), group_df[sb_ad_ref].values.astype(np.uint32), group_df[sb_ad_alt].values.astype(np.uint32))

    # Got a message 'KeyError: 'FE_P'' because of the line 'bulk_df['sSV'] = np.where(bulk_df['FE_P'] < alpha, 1, 0)' if without the line below
    df = df.drop(['Delta_AF', 'G_S', 'FE_P'], axis=1)
    df = df.merge(group_df[['CHROM', 'Seg', 'Delta_AF', 'G_S', 'FE_P']], on=['CHROM', 'Seg'], how='left')

    return df


def average_seg_statistics(df):
    steps = int(seg_size / ncrmntl_stp)
    daf_arr, gs_arr, fep_arr = np.zeros((steps, len(df.index))), np.zeros((steps, len(df.index))), np.zeros((steps, len(df.index)))
    
    for i in range(steps):
        df['Seg'] = (df.POS + ncrmntl_stp*i) // seg_size
        df = seg_statistics(df)
        daf_arr[i], gs_arr[i], fep_arr[i] = df['Delta_AF'], df['G_S'], df['FE_P']

    # print(daf_arr.shape, gs_arr.shape, fep_arr.shape)
    # print(np.std(daf_arr, axis=0), np.std(gs_arr, axis=0), np.std(fep_arr, axis=0))

    df = df.drop(['Delta_AF', 'G_S', 'FE_P'], axis=1)
    df['Delta_AF'] = np.mean(daf_arr, axis=0)
    df['G_S'] = np.mean(gs_arr, axis=0)
    df['FE_P'] = np.mean(fep_arr, axis=0)
    
    return df


t0 = time.time()

# Font settings for plotting
plt.rc('font', family='Arial', size=22)     # controls default text sizes
plt.rc('axes', titlesize=22)                # fontsize of the axes title
plt.rc('axes', labelsize=22)                # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)               # fontsize of the tick labels
plt.rc('ytick', labelsize=20)               # fontsize of the tick labels
plt.rc('legend', fontsize=20)               # legend fontsize
plt.rc('figure', titlesize=22)              # fontsize of the figure title
# plt.tick_params(labelsize=20)

ap = argparse.ArgumentParser()
ap.add_argument('-i', '--input', required=False, help='file names of the GATK4-generated tsv files: parental, bulk', type=lambda s: [t for t in s.split(',')], default='parents.tsv,bulks.tsv')
ap.add_argument('-o', '--output', required=False, help='file name of the output csv file', default='BSASeq.csv')
ap.add_argument('-b', '--bulksizes', required=False, help='bulk sizes: first_bulk,second_bulk', type=lambda s: [int(t) for t in s.split(',')], default='430,385')
ap.add_argument('-p', '--popstruct', required=False, choices=['F2','RIL','BC'], help='population structure', default='F2')
ap.add_argument('-v', '--pvalues', required=False, help='cutoff p-values: real_data,simulation', type=lambda s: [float(t) for t in s.split(',')], default='0.01,0.05')
ap.add_argument('-r', '--replication', type=int, required=False, help='the number of replications for threshold calculation', default=10000)
ap.add_argument('-s', '--slidingwindow', required=False, help='size,incremental_step', type=lambda s: [int(t) for t in s.split(',')], default='2000000,10000')
ap.add_argument('-g', '--gaps', required=False, help='gaps between subplots: horizontal,vertical', type=lambda s: [float(t) for t in s.split(',')], default='0.028,0.056,0.0264,0.054,0.076,0.002,0.002')
ap.add_argument('-m', '--smoothing', required=False, help='smoothing parameters: window_len,polyorder', type=lambda s: [int(t) for t in s.split(',')], default='51,3')
ap.add_argument('-t','--smooth', type=bool, required=False, help='smooth the plot', default=True)
ap.add_argument('--chromosome_order', required=False, help='manually set chromosome order', type=lambda s: [t for t in s.split(',')], default=False)
ap.add_argument('-d', '--read_length', type=int, required=False, help='make deta point independent of each other', default=0)
ap.add_argument('-e', '--region', required=False, help='interested region(s): chrm,start,end', type=lambda s: [int(t) for t in s.split(',')], default='-1')
ap.add_argument('-c', '--misc', required=False, help='cut-off GQ value, minimum SVs in a sliding window, extremely high read, and mirror index of Δ(allele frequency)', type=lambda s: [int(t) for t in s.split(',')], default='20,1,1')
ap.add_argument('-a','--adjust_gap', type=bool, required=False, help='adjust gaps between subplot', default=False)
ap.add_argument('--parent', required=False, type=str, help='The genome sequence of a parent is used as the reference for sv calling', default='ref')
ap.add_argument('--seg', required=False, help='Divide a chromosome into segments of the same size', type=lambda s: [int(t) for t in s.split(',')], default='600,0')

args = ap.parse_args()
chrm_order = args.chromosome_order
input_files = args.input
pop_struct = args.popstruct
rep = args.replication
sr_length = args.read_length
fb_size, sb_size = args.bulksizes[0], args.bulksizes[1]
alpha, sm_alpha = args.pvalues[0], args.pvalues[1]
sw_size, incremental_step = args.slidingwindow[0], args.slidingwindow[1]
h_gap, v_gap, t_margin, b_margin, l_margin, r_margin, g_pos = args.gaps[0], args.gaps[1], 1-args.gaps[2], args.gaps[3], args.gaps[4], 1-args.gaps[5], args.gaps[6]
smoothing = args.smooth
smth_wl, poly_order = args.smoothing[0], args.smoothing[1]
region = args.region
gq_value, min_SVs, mirror_index = args.misc[0], args.misc[1], args.misc[2]
gap_adjust = args.adjust_gap
parent1 = args.parent
seg_size, ncrmntl_stp = args.seg[0], args.seg[1]
print(seg_size, ncrmntl_stp)

num_ipfiles = len(input_files)
misc = []
min_frag_size = sw_size + smth_wl * incremental_step    # Minimum chromosome size allowed
# additional_peaks = ''
percentile_ci = [alpha*100/2, 100-alpha*100/2]          # [0.5, 99.5, 2.5, 97.5, 5.0, 95.0]
percentile_th = 100-alpha*100/2
print(percentile_ci, percentile_th)
curve_color, sv_threshold_color, sw_threshold_color = 'k', 'magenta', 'r'
ttl_sv_color, bg_color = 'b', 'gray'
test = False
min_sv = 500    # Ploting will not be performed if the selected chromosome/region has less SVs than min_sv

# Obtain the frequencies of the ALT allele in both bulks
fb_freq = sm_allelefreq(pop_struct, fb_size, rep)
sb_freq = sm_allelefreq(pop_struct, sb_size, rep)

# Set paths for input/output files
path = os.getcwd()
current_dt = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
oi_file = os.path.join(path, 'Results', 'sv_fagz.csv')
oi_file_di = os.path.join(path, 'Results', 'sv_fagz_di.csv')
bulk_homo_svs_file = os.path.join(path, 'Results', 'bulk_homo_svs.csv')
ssv_file = os.path.join(path, 'Results', 'ssv.csv')
sw_file = os.path.join(path, 'Results', 'sliding_windows.csv')
sv_region_file = os.path.join(path, 'Results', 'sv_region.csv')
thrshld_file = os.path.join(path, 'Results', 'threshold.txt')
# status_file = os.path.join(path, 'Results', 'COMPLETE.txt')
results = os.path.join(path, 'Results', current_dt)
dgns_path = os.path.join(path, results, 'temp')
ed_file = os.path.join(path, results, 'highly_ed.csv')

if ncrmntl_stp > 0:
    selected_oi_file = oi_file_di
else:
    selected_oi_file = oi_file

if not os.path.exists(results):
    os.makedirs(results)

selected_chrms = None
if gap_adjust==True:
    if os.path.isfile(thrshld_file) and os.path.isfile(sw_file):
        sw_dataframe = pd.read_csv(sw_file, dtype={'CHROM':str})

        temp = select_chrms(sw_dataframe, region)
        selected_chrms = temp[0]
        chrmSzD = temp[1]
        chrmSzL = temp[2]

        bsaseq_plot_sw(sw_dataframe)
    else:
        print('You need to run the script without the --gap-adjust option to generate the required files first')
        sys.exit()
else:
    if os.path.isfile(selected_oi_file):
        # additional_peaks = input('Do you want to have additional peaks identified (yes or no)?\n')
        # print('\n')
        bulk_df = pd.read_csv(selected_oi_file, dtype={'CHROM':str})

        bulk_list = bulk_names(bulk_df)
        fb_id, sb_id = bulk_list[0], bulk_list[1]
        fb_ad, sb_ad = fb_id + '.AD', sb_id + '.AD'
        fb_ad_ref, fb_ad_alt = fb_ad + '_REF', fb_ad + '_ALT'
        sb_ad_ref, sb_ad_alt = sb_ad + '_REF', sb_ad + '_ALT'
        fb_af, sb_af = fb_id + '.AF', sb_id + '.AF'
        fb_ld, sb_ld = fb_id + '.LD', sb_id + '.LD'
        sm_fb_ad_ref, sm_fb_ad_alt = 'sm_' + fb_ad_ref, 'sm_' + fb_ad_alt
        sm_sb_ad_ref, sm_sb_ad_alt = 'sm_' + sb_ad_ref, 'sm_' + sb_ad_alt

        bulk_df[sm_fb_ad_ref] = bulk_df[fb_ld] - bulk_df[sm_fb_ad_alt]
        bulk_df[sm_sb_ad_ref] = bulk_df[sb_ld] - bulk_df[sm_sb_ad_alt]

        temp = select_chrms(bulk_df, region)
        selected_chrms = temp[0]
        chrmSzD = temp[1]
        chrmSzL = temp[2]

        if os.path.isfile(thrshld_file) and os.path.isfile(sw_file) and os.path.isfile(sv_region_file):
            sw_dataframe = pd.read_csv(sw_file, dtype={'CHROM':str})
            bsaseq_plot_sw(sw_dataframe)

            sr_df = pd.read_csv(sv_region_file, dtype={'CHROM':str})
            sr_df['Peaks'] = sr_df['Peaks'].apply(ast.literal_eval)
            sv_region = sr_df.values.tolist()
        else:
            temp = select_chrms(bulk_df, region)
            selected_chrms = temp[0]
            chrmSzD = temp[1]
            chrmSzL = temp[2]

            sv_per_sw = int(len(bulk_df.index) * sw_size / sum(chrmSzL))

            # Calculate or retrieve the thresholds.
            thrshlds = gw_thresholds(bulk_df)
            thrshld_fe = thrshlds[0]
            thrshld_gs = thrshlds[1]
            thrshld_af_n = thrshlds[2][0]
            thrshld_af_p = thrshlds[2][1]
            thrshld_af_abs = thrshlds[3]
            sv_cis = sv_ci_gw(bulk_df)
            sv_ci_af_lb = sv_cis[0]
            sv_ci_af_ub = sv_cis[1]
            sv_ci_gs = sv_cis[2]
            sv_ci_af_abs = sv_cis[3]

            with open(thrshld_file, 'w') as xie:
                xie.write(' '.join([str(thrshld_fe), str(thrshld_gs), str(thrshld_af_n), str(thrshld_af_p), str(thrshld_af_abs), str(sv_ci_af_lb), str(sv_ci_af_ub), str(sv_ci_gs), str(sv_ci_af_abs)]))

            # Identify likely trait-associated SVs. This column is needed for plotting
            bulk_df = bulk_df.copy()
            bulk_df['sSV'] = np.where(bulk_df['FE_P'] < alpha, 1, 0)

            bsaseq_plot(bulk_df)
    else:
        m = 0
        for file in input_files:
            in_file = os.path.join(path, file)

            if '/' in file:
                sample = file.split('/')[-1]
            else:
                sample = file

            sample_name = sample.split('.')[0].lower()
            filtering_path = os.path.join(path, 'Results', 'FilteredSVs', sample_name)

            if in_file.endswith('.tsv'):
                separator = '\t'
            elif in_file.endswith('.csv'):
                separator = ','
            else:
                print('The input file should be in either the csv or the tsv format and should contain the corresponding file extension (.csv ot .tsv)')
                sys.exit()

            if not os.path.exists(filtering_path):
                os.makedirs(filtering_path)

            # Generte an SV dataframe from the GATK4-generated tsv file
            raw_sv_df = pd.read_csv(in_file, delimiter=separator, encoding='utf-8', dtype={'CHROM':str})

            if selected_chrms == None:         # Create chromosome list only once
                temp = select_chrms(raw_sv_df, region)
                selected_chrms = temp[0]
                chrmSzD = temp[1]
                chrmSzL = temp[2]

            if selected_chrms == []:
                print('No valid chromosomal names were entered.')
                sys.exit()

            # Create a numeric ID for each chromosome, which can be used to sort the dataframe numerically by chromosome number
            chrm_dict = {}
            for i in range(1, len(selected_chrms)+1):
                chrm_dict[selected_chrms[i-1]] = i

            raw_sv_df['ChrmSortID'] = raw_sv_df['CHROM']
            raw_sv_df['ChrmSortID'] = raw_sv_df['ChrmSortID'].map(chrm_dict)
            # raw_sv_df.replace({'ChrmSortID': chrm_dict}, inplace=True)
            # raw_sv_df['ChrmSortID'].replace(chrm_dict, inplace=True)

            # Obtain the bulk IDs from the header
            bulk_list = bulk_names(raw_sv_df)
            misc.extend([['Bulk ID - '+sample_name, bulk_list], ['Number of SVs in the entire dataframe - '+sample_name, len(raw_sv_df.index)]])

            fb_id, sb_id = bulk_list[0], bulk_list[1]
            fb_gt, sb_gt = fb_id + '.GT', sb_id + '.GT'
            fb_ad, sb_ad = fb_id + '.AD', sb_id + '.AD'
            fb_gq, sb_gq = fb_id + '.GQ', sb_id + '.GQ'

            # Check if any required field is missing in the input file
            required_fields, missing_fields = ['CHROM', 'POS', 'REF', 'ALT', fb_gt, sb_gt, fb_ad, fb_gq, sb_ad, sb_gq], []

            for elmt in required_fields:
                if elmt not in header:
                    missing_fields.append(elmt)

            if missing_fields !=[]:
                if len(missing_fields) == 1:
                    print('The following required field is missing: ', missing_fields)
                else:
                    print('The following required fields are missing: ', missing_fields)

                print('Please remake the input file to include the missing field(s).')
                sys.exit()

            # These variables will be used later, so the order in the sample list is very important
            fb_af, sb_af = fb_id + '.AF', sb_id + '.AF'
            fb_ad_ref, fb_ad_alt = fb_ad + '_REF', fb_ad + '_ALT'
            sb_ad_ref, sb_ad_alt = sb_ad + '_REF', sb_ad + '_ALT'
            fb_gt_ref, fb_gt_alt = fb_gt + '_REF', fb_gt + '_ALT'
            sb_gt_ref, sb_gt_alt = sb_gt + '_REF', sb_gt + '_ALT'
            fb_ld, sb_ld = fb_id + '.LD', sb_id + '.LD'
            sm_fb_ad_ref, sm_fb_ad_alt = 'sm_' + fb_ad_ref, 'sm_' + fb_ad_alt
            sm_sb_ad_ref, sm_sb_ad_alt = 'sm_' + sb_ad_ref, 'sm_' + sb_ad_alt

            sv_df = sv_filtering(raw_sv_df)

            # Create a unique ID for each SV (row) with its 'CHROM_POS'
            sv_df['ID'] = sv_df.CHROM.astype(str) + '_' + sv_df.POS.astype(str)

            if num_ipfiles == 2:
                if m == 0:
                    homo_svs = sv_df[((sv_df[fb_ad_ref]==0) & (sv_df[sb_ad_alt]==0)) | ((sv_df[fb_ad_alt]==0) & (sv_df[sb_ad_ref]==0))]
                    hetero_svs = sv_df.drop(index=homo_svs.index)

                    # sv_df.to_csv(os.path.join(filtering_path, 'Parents.csv'), index=None)
                    homo_svs.to_csv(os.path.join(filtering_path, 'homo_svs.csv'), index=None)
                    hetero_svs.to_csv(os.path.join(filtering_path, 'hetero_svs.csv'), index=None)

                    # Save the parent-specific variable names for later use
                    p_fb_gt, p_sb_gt = fb_gt, sb_gt
                    p_fb_id, p_sb_id = fb_id, sb_id
                    p_filtering_path = filtering_path
                else:
                    # Find the common SVs between the homoSVs and bulkSVs
                    bulk_df = sv_df[sv_df.ID.isin(homo_svs.ID)]
                    np_svs = sv_df.drop(index=bulk_df.index)
                    bulk_df.to_csv(os.path.join(filtering_path, 'bulk_df_before.csv'), index=None)
                    np_svs.to_csv(os.path.join(filtering_path, 'np_svs.csv'), index=None)

                    ref_svs = homo_svs[homo_svs.ID.isin(bulk_df.ID)]
                    ref_svs.to_csv(os.path.join(p_filtering_path, 'ref_svs.csv'), index=None)

                    if bulk_df.ID.tolist() != ref_svs.ID.tolist():
                        print('The ID sets of bulk_df and ref_svs are not the same.')
                        sys.exit()

                    bulk_df = bulk_df.copy()
                    # Add the genotypes of the parents to the bulk SV dataset
                    bulk_df['p_REF'] = ref_svs[p_fb_gt].str.split('/|\\|',expand=True)[0].tolist()
                    bulk_df['p_ALT'] = ref_svs[p_sb_gt].str.split('/|\\|',expand=True)[0].tolist()

                    # Swap the AD values and genotypes of the REF/ALT alleles in each bulk if the genotype of parent1 is different from the REF base
                    bulk_df[fb_gt] = np.where(bulk_df.p_REF==bulk_df.REF, bulk_df[fb_gt], bulk_df[fb_gt_alt] + '/' + bulk_df[fb_gt_ref])
                    bulk_df[sb_gt] = np.where(bulk_df.p_REF==bulk_df.REF, bulk_df[sb_gt], bulk_df[sb_gt_alt] + '/' + bulk_df[sb_gt_ref])
                    bulk_df[fb_ad] = np.where(bulk_df.p_REF==bulk_df.REF, bulk_df[fb_ad], bulk_df[fb_ad_alt].astype(str) + ',' + bulk_df[fb_ad_ref].astype(str))
                    bulk_df[sb_ad] = np.where(bulk_df.p_REF==bulk_df.REF, bulk_df[sb_ad], bulk_df[sb_ad_alt].astype(str) + ',' + bulk_df[sb_ad_ref].astype(str))
                    bulk_df[[fb_ad_ref, fb_ad_alt]] = bulk_df[fb_ad].str.split(',', expand=True).astype(int)
                    bulk_df[[sb_ad_ref, sb_ad_alt]] = bulk_df[sb_ad].str.split(',', expand=True).astype(int)
                    bulk_df['ad_Swap'] = np.where(bulk_df.p_REF==bulk_df.REF, 1, -1)
                    bulk_df.to_csv(os.path.join(filtering_path, 'bulk_df_after.csv'), index=None)
            else:
                bulk_df = sv_df.copy()

            m += 1

        misc.append(['Header', header])
        misc.extend([['Chromosome ID', selected_chrms]])
        misc.append(['Chromosome sizes', chrmSzL])

        bulk_df = sv_filtering_final(bulk_df)
        bulk_df.to_csv(os.path.join(filtering_path, 'bulk_svs.csv'), index=False)

        bulk_df = calculate_statistics(bulk_df)

        num_cols = [fb_ad_ref, fb_ad_alt, fb_ld, sb_ad_ref, sb_ad_alt, sb_ld]

        if sr_length > 0:
            bulk_df = bulk_df.copy()
            bulk_df['DI'] = di(bulk_df)
            bulk_df = bulk_df[bulk_df.DI==1]
            # bulk_df = impute_di(bulk_df)
            
        elif ncrmntl_stp > 0 and 'fisher' in sys.modules:
            bulk_df = average_seg_statistics(bulk_df)
            bulk_df.to_csv(oi_file_di, index=None)

        # Identify likely trait-associated SVs. This column is needed for plotting
        bulk_df = bulk_df.copy()
        bulk_df['sSV'] = np.where(bulk_df['FE_P'] < alpha, 1, 0)

        essential_fields = ['CHROM', 'POS', 'REF', 'ALT', fb_ad_ref, fb_ad_alt, fb_ld, sb_ad_ref, sb_ad_alt, sb_ld]
        ssv = bulk_df[bulk_df['FE_P'] < alpha]
        ssv.to_csv(ssv_file, columns=essential_fields, index=None)

        # Calculate the average number of SVs in a sliding window
        sv_per_sw = int(len(bulk_df.index) * sw_size / sum(chrmSzL))

        misc.append(['Average SVs per sliding window', sv_per_sw])
        misc.append([f'Average locus depth in bulk {fb_id}', bulk_df[fb_ld].mean()])
        misc.append([f'Average locus depth in bulk {sb_id}', bulk_df[sb_ld].mean()])

        # Calculate or retrieve the thresholds.
        thrshlds = gw_thresholds(bulk_df)
        thrshld_fe = thrshlds[0]
        thrshld_gs = thrshlds[1]
        thrshld_af_n = thrshlds[2][0]
        thrshld_af_p = thrshlds[2][1]
        thrshld_af_abs = thrshlds[3]
        sv_cis = sv_ci_gw(bulk_df)
        sv_ci_af_lb = sv_cis[0]
        sv_ci_af_ub = sv_cis[1]
        sv_ci_gs = sv_cis[2]
        sv_ci_af_abs = sv_cis[3]

        with open(thrshld_file, 'w') as xie:
            xie.write(' '.join([str(thrshld_fe), str(thrshld_gs), str(thrshld_af_n), str(thrshld_af_p), str(thrshld_af_abs), str(sv_ci_af_lb), str(sv_ci_af_ub), str(sv_ci_gs), str(sv_ci_af_abs)]))

        # Perform plotting
        bsaseq_plot(bulk_df)

    peaklst = pk_list(sv_region)

    peak_verification(peaklst, bulk_df)

    misc.append(['Running time', [(time.time()-t0)/60]])

    with open(os.path.join(results, 'misc_info.csv'), 'w', newline='') as outf:
        xie = csv.writer(outf)
        xie.writerows(misc)

print('\nIf two or more peaks/valleys and all the values in between are beyond the confidence intervals/thresholds, only the highest peak or the lowerest valley will be identified as the peak/valley of this region. The positions of the other peaks/valleys can be identified and their significance can be verified by rerunning the script using the region option. An example is provied below:')
print('python PyBSASeq.py -i parents.csv,bulks.csv --region 1,1000001,4000000,1,6000001,9000000,10,20000001,40000000')
print('For the region argument, the first digit is the chromosome ID, while the second and the third digits are the startpoint and the endpoint of the chromosome region of interest, respectively. You can add as many regions as you need.\n')
# Usage for changing gap sizes:
# $ python /media/data/Code/Python/PyBSASeq/PyBSASeq.py -i /media/data/Code/Python/PyBSASeq/Data/sub_parents.csv,/media/data/Code/Python/PyBSASeq/Data/sub_bulks.csv -a True -g '0.028,0.009,0.0264,0.054,0.076,0.002,0.002'
