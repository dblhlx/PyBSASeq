"""
Created on Fri Oct  5 08:22:16 2018
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
from scipy.signal import savgol_filter


def smAlleleFreq(popStruc, sizeOfBulk, rep):
    '''
    This function is for the calculation of the allele frequency in the bulks via simulation
    under the null hypothesis.

    An AA, Aa, and aa individual carries 0%, 50%, and 100% of the alt (a) allele, respectively.
    If A/a is not associated with the trait (null hypothesis), the AA:Aa:aa ratios are 0.25:0.5:0.25, 0.5:0:0.5, 
    and 0.5:0.5:0, respectively, in a F2 population, in a RIL population, and in a back crossed population.
    '''
    freqL = []
    pop = [0.0, 0.5, 1.0]

    if popStruc == 'F2':
        prob = [0.25, 0.5, 0.25]
    elif popStruc == 'RIL':
        prob = [0.5, 0.0, 0.5]
    elif popStruc == 'BC':
        prob = [0.5, 0.5, 0.0]

    for __ in range(rep):
        altFreq = np.random.choice(pop, sizeOfBulk, p=prob).mean()
        freqL.append(altFreq)

    return sum(freqL)/len(freqL)


def chrmFiltering(df, chromosomeList):
    # Many reference genomes contain unmapped fragments that tend to be small and are not informative for SNP-trait association, filtering them out makes the chromosome list more readable.
    chrmSizeD, smallChrmL = {}, []
    for chrmID in chromosomeList:
        chrmSize = df[df.CHROM==chrmID]['POS'].max()

        if chrmSize <= minFragSize:
            smallChrmL.append(chrmID)
        else:
            chrmSizeD[chrmID] = [1, chrmSize]

    for chrmID in smallChrmL:
        chromosomeList.remove(chrmID)

    return [chrmSizeD, chromosomeList, smallChrmL]


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


def snpFiltering(df):
    print('Perform SNP filtering')
    global misc

    df = df.copy()

    # Identify SNPs not informative and remove these SNPs from the dataframe
    df_Ignored = df[~df.CHROM.isin(chrmIDL)]
    df_Ignored.to_csv(os.path.join(filteringPath, 'ignored.csv'), index=None)
    df = df.drop(index=df_Ignored.index)

    # Identify SNPs with an 'NA' value(s) and remove these SNPs from the dataframe
    df_NA = df[df.isnull().any(axis=1)]
    df_NA.to_csv(os.path.join(filteringPath, 'na.csv'), index=None)
    df.dropna(inplace=True)
    misc.append(['Number of SNPs after NA drop', len(df.index)])

    # Filter out SNPs with a low genotype quality score
    df_lowq = df[(df[fb_GQ]<gqValue) | (df[sb_GQ]<gqValue)]
    df_lowq.to_csv(os.path.join(filteringPath, 'lowqualitySNP.csv'), index=None)
    df = df.drop(index=df_lowq.index)

    # Filter out 1 bp deletions
    inDel_1 = df[df.ALT.str.contains(r'\*')]
    df = df.drop(index=inDel_1.index)

    # Identify SNPs with more than one ALT allele
    df_mALT = df[df.ALT.str.contains(',')]

    # Identify SNPs with a single ALT allele
    df_1ALT = df.drop(index=df_mALT.index)

    # Identify one-ALT SNPs with zero REF read in both bulks
    df_1ALT_Fake = df_1ALT[(df_1ALT[fb_AD].str.startswith('0')) & \
        (df_1ALT[sb_AD].str.startswith('0'))]
    df_1ALT_Fake.to_csv(os.path.join(filteringPath, '1altFake.csv'), index=None)

    # Remove one-ALT SNPs with zero REF read in both bulks
    df_1ALT_Real = df_1ALT.drop(index=df_1ALT_Fake.index)
    df_1ALT_Real.to_csv(os.path.join(filteringPath, '1altReal.csv'), index=None)

    # Using 'str.count' make the code below simpler and easier to understand
    df_3ALT = df_mALT[df_mALT.ALT.str.count(',') >= 2]
    df_3ALT.to_csv(os.path.join(filteringPath, '3ALT.csv'), index=None)
    df_2ALT = df_mALT.drop(index=df_3ALT.index)

    # A two-ALT SNP is a real SNP if the REF read is zero in both bulks
    df_2ALT_Real = df_2ALT[(df_2ALT[fb_AD].str.startswith('0')) & \
        (df_2ALT[sb_AD].str.startswith('0'))]

    # The two-ALT SNP may be caused by allele heterozygosity if the REF read is not zero
    # Repetitive sequences in the genome or sequencing artifacts are other possibilities
    df_2ALT_Het = df_2ALT.drop(index=df_2ALT_Real.index)
    df_2ALT_Het.to_csv(os.path.join(filteringPath, '2ALT_Het.csv'), index=None)

    # Making a copy to suppress the warning message. Updating the AD values of these SNPs is required
    df_2ALT_Real = df_2ALT_Real.copy()
    df_2ALT_Real.to_csv(os.path.join(filteringPath, '2altReal_Before.csv'), index=None)

    # Update the AD values of the above SNPs by removing the REF read that is zero (i.e. remove '0,' from '0,x,y')
    if not df_2ALT_Real.empty:
        df_2ALT_Real[fb_AD] = df_2ALT_Real[fb_AD].str.slice(start=2)
        df_2ALT_Real[sb_AD] = df_2ALT_Real[sb_AD].str.slice(start=2)
        df_2ALT_Real[['REF', 'ALT']] = df_2ALT_Real['ALT'].str.split(',', expand=True)
        df_2ALT_Real.to_csv(os.path.join(filteringPath, '2altReal_After.csv'), index=None)

    # Concatenate 1ALT_Real and 2ALT_Real
    snp = pd.concat([df_1ALT_Real, df_2ALT_Real])

    # Identify inDels, REF/Alt allele with more than 1 base
    inDel_2 = snp[(snp.REF.str.len()>1) | (snp.ALT.str.len()>1)]

    # Remove InDels from the SNP dataframe
    snp = snp.drop(index=inDel_2.index)

    df_InDel = pd.concat([inDel_1, inDel_2])
    df_InDel.to_csv(os.path.join(filteringPath, 'InDel.csv'), index=None)

    snp.sort_values(['ChrmSortID', 'POS'], inplace=True)

    # Obtain REF reads, ALT reads, and locus reads of each SNP
    snp[[fb_AD_REF, fb_AD_ALT]] = snp[fb_AD].str.split(',', expand=True).astype(int)
    snp[fb_LD] = snp[fb_AD_REF] + snp[fb_AD_ALT]
    snp[[sb_AD_REF, sb_AD_ALT]] = snp[sb_AD].str.split(',', expand=True).astype(int)
    snp[sb_LD] = snp[sb_AD_REF] + snp[sb_AD_ALT]

    # Filter out the SNPs with zero locus reads in either bulk
    snp_0LD = snp[(snp[fb_LD] == 0) | (snp[sb_LD] == 0)]
    snp_0LD.to_csv(os.path.join(filteringPath, '0ld.csv'), index=None)
    snp = snp.drop(index=snp_0LD.index)

    # Filter out SNPs in which GT and AD are not consistent
    snp[[fb_GT_REF, fb_GT_ALT]] = snp[fb_GT].str.split('/|\|', expand=True)
    snp[[sb_GT_REF, sb_GT_ALT]] = snp[sb_GT].str.split('/|\|', expand=True)
    gt_ad = snp[(snp[fb_GT_REF]==snp[fb_GT_ALT]) & (snp[sb_GT_REF]==snp[sb_GT_ALT]) & (snp[fb_GT_REF]==snp[sb_GT_REF])]
    snp = snp.drop(index=gt_ad.index)
    gt_ad.to_csv(os.path.join(filteringPath, 'gt_ad.csv'), index=None)

    # Filter out SNPs in which bulk GT and REF/ALT are not consistent
    gt = snp[~(((snp[fb_GT_REF]==snp.REF) | (snp[fb_GT_REF]==snp.ALT)) & \
               ((snp[fb_GT_ALT]==snp.REF) | (snp[fb_GT_ALT]==snp.ALT)) & \
               ((snp[sb_GT_REF]==snp.REF) | (snp[sb_GT_REF]==snp.ALT)) & \
               ((snp[sb_GT_ALT]==snp.REF) | (snp[sb_GT_ALT]==snp.ALT)))]
    gt.to_csv(os.path.join(filteringPath, 'gt.csv'), index=None)
    snp = snp.drop(index=gt.index)

    # A SNP with very high LD is likely from the repetitive genomic sequence
    # Remove SNPs that could be from the repetitive elements
    fb_repLD, sb_repLD = hiReadFold * snp[fb_LD].mean(), hiReadFold * snp[sb_LD].mean() 
    snpRep = snp[(snp[fb_LD]>fb_repLD) | (snp[sb_LD]>sb_repLD)]
    snpRep.to_csv(os.path.join(filteringPath, 'repetitiveSeq_WoP.csv'), index=None)
    snp = snp.drop(index=snpRep.index)

    print(f'SNP filtering completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return snp


def gStatistic_Array(o1, o3, o2, o4):
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


def statisticsFE(row):
    '''
    Perform Fisher's exact test for each SNP in the dataset using its actual AD values in both bulks, 
    and perform Fisher's exact test and estimate the thresholds of Δ(allele frequency) and G-statistic 
    value using its simulated AD values.
    This function will be used if the module 'Fisher' is not installed. It is slow for large dataset.
    '''
    try:
        fe = fisher_exact([[row[fb_AD_REF], row[fb_AD_ALT]], [row[sb_AD_REF], row[sb_AD_ALT]]])
    except TypeError:
        fe = ('NA', 'NA')

    # Perform Fisher's exact test for each SNP using the simulated REF/ALT reads
    try:
        sm_FE = fisher_exact([[row[sm_fb_AD_REF], row[sm_fb_AD_ALT]], [row[sm_sb_AD_REF], row[sm_sb_AD_ALT]]])
    except TypeError:
        sm_FE = ('NA', 'NA')

    # Create an array with 10000 (rep) simulated ALT reads - first bulk
    sm_yb_Alt_Array = np.random.binomial(row[fb_LD], fb_Freq, rep)
    yb_LD_Array = np.full(rep, row[fb_LD])
    # Create an array with 10000 (rep) simulated ALT reads - second bulk
    sm_eb_Alt_Array = np.random.binomial(row[sb_LD], sb_Freq, rep)
    eb_LD_Array = np.full(rep, row[sb_LD])

    # Create simulated allele frequency and Δ(allele frequency) arrays of a SNP
    sm_yb_AF_Array = sm_yb_Alt_Array/yb_LD_Array
    sm_eb_AF_Array = sm_eb_Alt_Array/eb_LD_Array
    sm_DAF_Array = sm_eb_AF_Array - sm_yb_AF_Array

    # Create a G-statistic array
    sm_GS_Array = gStatistic_Array(sm_yb_Alt_Array, yb_LD_Array-sm_yb_Alt_Array, sm_eb_Alt_Array, eb_LD_Array-sm_eb_Alt_Array)

    # Obtain the percentile of the above arrays for threshold estimation
    ci_yb_AF = np.percentile(sm_yb_AF_Array, percentile_list)
    ci_eb_AF = np.percentile(sm_eb_AF_Array, percentile_list)
    ci_DAF = np.percentile(sm_DAF_Array, percentile_list)
    ci_GS = np.percentile(sm_GS_Array, percentile_list)

    return [fe, sm_FE, ci_yb_AF, ci_eb_AF, ci_DAF, ci_GS]


def statistics(row):
    '''
    Estimate the thresholds of Δ(allele frequency) and G-statistic value using its simulated AD values 
    for each SNP in the dataset.
    This function will be used if the module 'Fisher' is available.
    '''
    # Create an array with 10000 (rep) simulated ALT reads of a SNP - first bulk
    sm_yb_Alt_Array = np.random.binomial(row[fb_LD], fb_Freq, rep)
    yb_LD_Array = np.full(rep, row[fb_LD])
    # Create an array with 10000 (rep) simulated ALT reads of a SNP - second bulk
    sm_eb_Alt_Array = np.random.binomial(row[sb_LD], sb_Freq, rep)
    eb_LD_Array = np.full(rep, row[sb_LD])

    # Create simulated allele frequency and Δ(allele frequency) arrays of a SNP
    sm_yb_AF_Array = sm_yb_Alt_Array/yb_LD_Array
    sm_eb_AF_Array = sm_eb_Alt_Array/eb_LD_Array
    sm_DAF_Array = sm_eb_AF_Array - sm_yb_AF_Array

    # Create a G-statistic array of a SNP
    sm_GS_Array = gStatistic_Array(sm_yb_Alt_Array, yb_LD_Array-sm_yb_Alt_Array, sm_eb_Alt_Array, eb_LD_Array-sm_eb_Alt_Array)

    # Obtain the percentile of the above arrays
    ci_yb_AF = np.percentile(sm_yb_AF_Array, percentile_list)
    ci_eb_AF = np.percentile(sm_eb_AF_Array, percentile_list)
    ci_DAF = np.percentile(sm_DAF_Array, percentile_list)
    ci_GS = np.percentile(sm_GS_Array, percentile_list)

    return [ci_yb_AF, ci_eb_AF, ci_DAF, ci_GS]


def smThresholds_proximal(DF):
    # Calculate the sSNP/totalSNP threshold via re-sampling if 'Fisher' is not installed
    print('Calculate the threshold of sSNPs/totalSNPs.')
    ratioLi = []        # List containing simulated sSNP/totalSNP ratios
    for __ in range(rep):
        # Sample the number of SNPs equal to the average number of SNPs per sliding window
        sm_SNP_SMPL = DF.sample(snpPerSW, replace=True)
        # Identify sSNPs in the above SNP subset based on the simulated AD values
        sm_sSNP_SMPL = sm_SNP_SMPL[sm_SNP_SMPL['sm_FE_P']<smAlpha]
        ratioLi.append(len(sm_sSNP_SMPL.index)/snpPerSW)

    misc.append(['Genome-wide sSNP/totalSNP ratio threshold', np.percentile(ratioLi, percentile_list)])
    print(f'Threshold calculation completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return np.percentile(ratioLi, percentile_list)


def smThresholds_gw(DF):
    # Calculate the genome-wide sSNP/totalSNP threshold if 'Fisher' is available
    print('Calculate the threshold of sSNPs/totalSNPs.')
    gw_ratioLi = []     # List containing simulated sSNP/totalSNP ratios
    for __ in range(rep):
        # Sample the number of SNPs equal to the average number of SNPs per sliding window
        sm_SNP_SMPL = DF.sample(snpPerSW, replace=True)
        # Convert the LD columns of the above SNP subset to arrays
        gw_sm_fb_LD_Arr = sm_SNP_SMPL[fb_LD].to_numpy().astype(np.uint)
        gw_sm_sb_LD_Arr = sm_SNP_SMPL[sb_LD].to_numpy().astype(np.uint)

        # Create arrays containing simulated AD values for the above SNP subset
        gw_sm_fb_AD_ALT_Arr = np.random.binomial(sm_SNP_SMPL[fb_LD], fb_Freq).astype(np.uint)
        gw_sm_fb_AD_REF_Arr = gw_sm_fb_LD_Arr - gw_sm_fb_AD_ALT_Arr
        gw_sm_sb_AD_ALT_Arr = np.random.binomial(sm_SNP_SMPL[sb_LD], sb_Freq).astype(np.uint)
        gw_sm_sb_AD_REF_Arr = gw_sm_sb_LD_Arr - gw_sm_sb_AD_ALT_Arr

        # Calculate the P-value of each SNP using the simulated AD values
        __, __, gw_sm_FE_P_Arr = pvalue_npy(gw_sm_fb_AD_ALT_Arr, gw_sm_fb_AD_REF_Arr, gw_sm_sb_AD_ALT_Arr, gw_sm_sb_AD_REF_Arr)
        # gw_sm_FE_OR_Arr = (gw_sm_fb_AD_ALT_Arr * gw_sm_sb_AD_REF_Arr) / (gw_sm_fb_AD_REF_Arr * gw_sm_sb_AD_ALT_Arr)

        # Based on the above P-value array, create an array in which the sSNP has value 1 while the non-sSNP has value 0, 
        # and its mean is the simulated sSNP/totalSNP ratio of the above SNP subset
        sSNP_Arr = np.where(gw_sm_FE_P_Arr<smAlpha, 1, 0)
        gw_ratioLi.append(np.mean(sSNP_Arr))

    misc.append(['Genome-wide sSNP/totalSNP ratio threshold', np.percentile(gw_ratioLi, percentile_list)])
    print(f'Threshold calculation completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return np.percentile(gw_ratioLi, percentile_list)


def smThresholds_sw(DF):
    # Calculate the sliding window-specific threshold
    sw_ratioLi = []     # List containing simulated sSNP/totalSNP ratios

    # Convert the LD columns of the sliding window to arrays
    sw_fb_LD_Arr = DF[fb_LD].to_numpy().astype(np.uint)
    sw_sb_LD_Arr = DF[sb_LD].to_numpy().astype(np.uint)

    for __ in range(rep):
        # Create arrays containing simulated AD values for the sliding window
        sw_sm_fb_AD_ALT_Arr = np.random.binomial(DF[fb_LD], fb_Freq).astype(np.uint)
        sw_sm_fb_AD_REF_Arr = sw_fb_LD_Arr - sw_sm_fb_AD_ALT_Arr
        sw_sm_sb_AD_ALT_Arr = np.random.binomial(DF[sb_LD], sb_Freq).astype(np.uint)
        sw_sm_sb_AD_REF_Arr = sw_sb_LD_Arr - sw_sm_sb_AD_ALT_Arr

        # Calculate the P-value of each SNP using the simulated AD values
        __, __, sw_sm_FE_P_Arr = pvalue_npy(sw_sm_fb_AD_ALT_Arr, sw_sm_fb_AD_REF_Arr, sw_sm_sb_AD_ALT_Arr, sw_sm_sb_AD_REF_Arr)
        # sw_sm_FE_OR_Arr = (sw_sm_fb_AD_ALT_Arr * sw_sm_sb_AD_REF_Arr) / (sw_sm_fb_AD_REF_Arr * sw_sm_sb_AD_ALT_Arr)

        # Based on the above P-value array, create an array in which the sSNP has value 1 while the non-sSNP has value 0, 
        # and its mean is the simulated sSNP/totalSNP ratio of the sliding window
        sSNP_Arr = np.where(sw_sm_FE_P_Arr<smAlpha, 1, 0)

        sw_ratioLi.append(np.mean(sSNP_Arr))

    return np.percentile(sw_ratioLi, percentile_list)


def zeroSNP(li):
    # Replace 'divide by zero' with the nearnest value. 'li' is a list
    if li != []:
        # Assign the previous value to the empty sliding window if the list is not empty
        li.append(li[-1])
    else:
        # Assign 'empty' to the first sliding window that is empty
        li.append('empty')  


def replaceZero(li):
    # Replace the 'empty' placeholders at the begining of the list with the nearnest non-empty value
    # Find the first non-empty value
    i = 0
    while li[i]=='empty':
        i += 1

    # Replace the 'empty' value(s) in the list with the above non-empty value
    j = 0
    while j < i:
        li[j] = li[i]
        j += 1


def bsaseqPlot(chrmIDL, datafr):
    '''
    wmL: list of warning messages
    swDict: a dictionary with the chromosome ID as its keys; the value of each key is a list containing
            the chromosome ID, the sSNP/totalSNP ratio in each sliding window, and the startpoint of the 
            sliding window
    misc: miscellaneous information
    snpRegion: genomic region above the threshold
    swDataFrame: DataFrame containing a sliding window
    '''
    print('Prepare SNP data for plotting via the sliding window algorithm')
    global misc
    global snpRegion, swDataFrame
    sg_yRatio_List = []     # Smoothed sSNP/totalSNP ratios of a chromosome or a selected genomic region
    swRows = []
    wmL, swDict, snpRegion = [], {}, []

    # Analyze each chromosome separately
    numOfSNPOnChr = []      # List containing the sSNP/totalSNP, # of sSNPs, and # of totalSNPs of a chromosome
    ratioPeakL = []         # List containing the peak of each chromosome 
    i = 1
    for chrmID in chrmIDL:
        ch = datafr[datafr.CHROM==chrmID]
        numOfSNPOnChr.append([chrmID, ch['sSNP'].sum(), len(ch.index), ch['sSNP'].sum()/len(ch.index)])

        regStart = chrmSzD[chrmID][0]   # The startpoint of the selected chromosomal region
        regEnd = chrmSzD[chrmID][1]     # The endpoint of the selected chromosomal region
        swStr = regStart                # The begining of the window
        swEnd = regStart+swSize         # The end of the window 
        icrs = incrementalStep          # The incremental step
        plotSP = regStart               # The start point of the selected genomic region in the plot
        x = []              # Sliding window position on the chromosome
        y = []              # Number of sSNPs in the sliding window
        yT = []             # Number of total SNPs in the sliding window
        yRatio = []         # sSNP/totalSNP ratio of the sliding window
        # fb_ZRatio = []      # sSNPs/totalSNPs in the first bulk via z-test
        # sb_ZRatio = []      # sSNPs/totalSNPs in the second bulk via z-test
        y5 = []             # G-statistic
        y6 = []             # The threshold of the G-statistic
        y7 = []             # y7: Δ(allele frequency)
        y8, y9 = [], []     # The confidence interval of the Δ(allele frequency)

        # Calculate, sSNP/totalSNP, G-statistic, and Δ(allele frequency) of a sliding window
        while swEnd <= regEnd:
            # swDF: Dataframe containing the sliding window
            swDF = ch[(ch.POS>=swStr) & (ch.POS<=swEnd)]

            rowInSwDF = len(swDF.index)       # number of SNPs in a sliding window
            x.append(swStr)                   # Append the starpoint of the sliding window to x
            y.append(swDF['sSNP'].sum())      # Append number of sSNPs of the sliding window to y
            yT.append(rowInSwDF)              # Append number of totalSNPs of the sliding window to yT

            # 'try/exception' cannot catch the EmptyDataError; seems it only works when reading a .csv/.tsv file, not an empty subset of a existing dataframe. DivisionByZero generates 'nan' for a series or an array.
            if rowInSwDF >= minSNPs:
                yRatio.append(swDF['sSNP'].mean())    # Append the sSNP/totalSNP ratio of the sliding window to yRatio
                # fb_ZRatio.append(swDF['fb_Z'].mean())
                # sb_ZRatio.append(swDF['sb_Z'].mean())
                y5.append(swDF['G_S'].sum()/rowInSwDF)
                y6.append(swDF['GS_CI0995'].sum()/rowInSwDF)
                y7.append(swDF['Delta_AF'].sum()/rowInSwDF)
                y8.append(swDF['DAF_CI0005'].sum()/rowInSwDF)
                y9.append(swDF['DAF_CI0995'].sum()/rowInSwDF)

                # useful info of the sliding window
                rowContents = [chrmID, swStr, int(swDF[fb_LD].mean()), int(swDF[sb_LD].mean()), swDF['sSNP'].sum(), rowInSwDF, yRatio[-1], y5[-1], y6[-1], y7[-1], y8[-1], y9[-1]]
            else:
                wmL.append(['No SNP in the sliding window', i, swStr])
                zeroSNP(yRatio)
                # zeroSNP(fb_ZRatio)
                # zeroSNP(sb_ZRatio)
                zeroSNP(y5)
                zeroSNP(y6)
                zeroSNP(y7)
                zeroSNP(y8)
                zeroSNP(y9)

                # The mean of an empty column is 'nan', not zero, and int(nan) generates a ValueError
                rowContents = [chrmID, swStr, 0, 0, swDF['sSNP'].sum(), rowInSwDF, yRatio[-1], y5[-1], y6[-1], y7[-1], y8[-1], y9[-1]]

            swRows.append(rowContents)

            if i not in swDict:
                swDict[i] = [rowContents]
            else:
                swDict[i].append(rowContents)

            swStr += icrs
            swEnd += icrs

        # Replace the 'empty' values at the begining of the lists with nearest non-empty value
        for yl in [yRatio, y5, y6, y7, y8, y9]:
            if 'empty' in yl:
                replaceZero(yl)

        # Find the first non-empty value in the swDict
        pIndex = 0
        while swDict[i][pIndex][6] == 'empty':
            pIndex += 1

        # Replace the empty value(s) with the above non-empty value in the swDict
        j = 0
        while j < pIndex:
            swDict[i][j][6] = swDict[i][pIndex][6]
            swDict[i][j][7] = swDict[i][pIndex][7]
            swDict[i][j][8] = swDict[i][pIndex][8]
            swDict[i][j][9] = swDict[i][pIndex][9]
            swDict[i][j][10] = swDict[i][pIndex][10]
            j += 1

        # Smoothing data of a chromosome or a selected region
        sg_yRatio = savgol_filter(yRatio, smthWL, polyOrder)
        sg_y5 = savgol_filter(y5, smthWL, polyOrder)
        sg_y6 = savgol_filter(y6, smthWL, polyOrder)
        sg_y7 = savgol_filter(y7, smthWL, polyOrder)
        sg_y8 = savgol_filter(y8, smthWL, polyOrder)
        sg_y9 = savgol_filter(y9, smthWL, polyOrder)
        # sg_fb_ZRatio = savgol_filter(fb_ZRatio, smthWL, polyOrder)
        # sg_sb_ZRatio = savgol_filter(sb_ZRatio, smthWL, polyOrder)

        sg_yRatio_List.extend(sg_yRatio)

        # Handle the plot with a single column (chromosome)
        if len(chrmIDL) == 1:
            # Set up x-ticks
            axs[0].set_xticks(np.arange(0, max(x), xt_pro[0]))
            ticks = axs[0].get_xticks()*xt_pro[1]
            axs[0].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0].set_ylabel('Number of SNPs')
                axs[1].set_ylabel(r'sSNP/totalSNP')
                axs[2].set_ylabel('G-statistic')
                axs[3].set_ylabel('\u0394AF')

            # Plot sSNPs and total SNPs against their genomic positions
            axs[0].plot(x, y, c='k')
            axs[0].plot(x, yT, c='b')
            if chrmID.isdigit() == True:
                axs[0].set_title('Chr'+chrmID)
            else:
                axs[0].set_title(chrmID)

            # Plot sSNP/totalSNP, G-statistic, and Δ(allele frequency) against their genomic positions
            if smoothing == True:
                # sSNPs/totalSNPs via Fisher's exact test
                axs[1].plot(x, sg_yRatio, c='k')

                # G-statistic
                axs[2].plot(x, sg_y5, c='k')
                axs[2].plot(x, sg_y6, c='r')

                # Δ(allele frequency)
                axs[3].plot(x, sg_y7, c='k')
                axs[3].plot(x, sg_y8, c='r')
                axs[3].plot(x, sg_y9, c='r')

                # sSNPs/totalSNPs plot via z-test
                # axs[1].plot(x, sg_fb_ZRatio, c='c')
                # axs[1].plot(x, sg_sb_ZRatio, c='g')

            else:
                # sSNPs/totalSNPs via Fisher's exact test
                axs[1].plot(x, yRatio, c='k')

                # G-statistic
                axs[2].plot(x, y5, c='k')
                axs[2].plot(x, y6, c='r')

                # Δ(allele frequency)
                axs[3].plot(x, y7, c='k')
                axs[3].plot(x, y8, c='r')
                axs[3].plot(x, y9, c='r')

                # sSNPs/totalSNPs plot via z-test
                # axs[1].plot(x, fb_ZRatio, c='c')
                # axs[1].plot(x, sb_ZRatio, c='g')

            # Add the 99.5 percentile line as the threshold, x[-1] is the startpoint of the last sliding window of a chromosome
            axs[1].plot([plotSP, x[-1]], [thrshld, thrshld], c='r')
            # axs[1].plot(x, smThresholds_sw, c='m')

        # Handle the plot with multiple columns (chromosomes)
        else:
            # Set up x-ticks
            axs[0,i-1].set_xticks(np.arange(0, max(x), xt_pro[0]))
            ticks = axs[0,i-1].get_xticks()*xt_pro[1]
            axs[0,i-1].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0,i-1].set_ylabel('Number of SNPs')
                axs[1,i-1].set_ylabel(r'sSNP/totalSNP')
                axs[2,i-1].set_ylabel('G-statistic')
                axs[3,i-1].set_ylabel('\u0394AF')

            # Plot sSNP and totalSNPs against their genomic positions
            axs[0,i-1].plot(x, y, c='k')
            axs[0,i-1].plot(x, yT, c='b')
            if chrmID.isdigit() == True:
                axs[0,i-1].set_title('Chr'+chrmID)
            else:
                axs[0,i-1].set_title(chrmID)

            # Plot sSNP/totalSNP, G-statistic, and Δ(allele frequency) against their genomic positions
            if smoothing == True:
                # sSNPs/totalSNPs via Fisher's exact test
                axs[1,i-1].plot(x, sg_yRatio, c='k')

                # G-statistic
                axs[2,i-1].plot(x, sg_y5, c='k')
                axs[2,i-1].plot(x, sg_y6, c='r')

                # Δ(allele frequency)
                axs[3,i-1].plot(x, sg_y7, c='k')
                axs[3,i-1].plot(x, sg_y8, c='r')
                axs[3,i-1].plot(x, sg_y9, c='r')

                # sSNPs/totalSNPs plot via z-test
                # axs[1,i-1].plot(x, sg_fb_ZRatio, c='c')
                # axs[1,i-1].plot(x, sg_sb_ZRatio, c='g')
            else:
                # sSNPs/totalSNPs via Fisher's exact test
                axs[1,i-1].plot(x, yRatio, c='k')
                
                # G-statistic
                axs[2,i-1].plot(x, y5, c='k')
                axs[2,i-1].plot(x, y6, c='r')

                # Δ(allele frequency)
                axs[3,i-1].plot(x, y7, c='k')
                axs[3,i-1].plot(x, y8, c='r')
                axs[3,i-1].plot(x, y9, c='r')

                # sSNPs/totalSNPs plot via z-test
                # axs[1,i-1].plot(x, fb_ZRatio, c='c')
                # axs[1,i-1].plot(x, sb_ZRatio, c='g')

            # Add the 99.5 percentile line as the threshold, x[-1] is the startpoint of the last sliding window of a chromosome
            axs[1,i-1].plot([plotSP, x[-1]], [thrshld, thrshld], c='r')
            # axs[1, i-1].plot(x, smThresholds_sw, c='m')

        ratioPeakL.append(max(yRatio))

        # Identify genomic regions above the threshold
        zigzag = []     # List of peaks above the threshold
        # Handle the case in which an QTL is at the very begining of the chromosome
        if swDict[i][0][6] >= thrshld:
            # Add the CHROM ID and the startpoint to the genomic region above the threshold
            snpRegion.append(swDict[i][0][:2])
            if swDict[i][0][6] >= swDict[i][1][6]:
                # Add the first peak above the threshold to zigzag
                zigzag.append(swDict[i][0][1:])
            numOfSWs = 1

        m = 0       # The index of the sliding window of a chromosome
        while m < len(swDict[i]) - 1:
            if swDict[i][m][6] < thrshld and swDict[i][m+1][6] >= thrshld:
                # Add the CHROM ID and the startpoint to the genomic region above the threshold
                snpRegion.append(swDict[i][m+1][:2])
                numOfSWs = 1
            elif swDict[i][m][6] >= thrshld:
                # A sliding window is considered as a peak if its sSNP/totalSNP is greater than or equal to the threshold and greater than those of the flanking sliding windows
                if m >= 1 and max(swDict[i][m-1][6], swDict[i][m+1][6]) <= swDict[i][m][6]:
                    # Add a peak above the threshold to zigzag
                    zigzag.append(swDict[i][m][1:])
                if swDict[i][m+1][6] > thrshld:
                    numOfSWs += 1
                elif swDict[i][m+1][6] < thrshld:
                    # Add the endpoint, all peaks, the number of sliding widows
                    snpRegion[-1].extend([swDict[i][m][1], zigzag, numOfSWs])
                    zigzag = []
            m += 1
        # Handle the case in which an QTL is nearby the end of the chromosome
        if swDict[i][-1][6] >= thrshld:
            snpRegion[-1].extend([swDict[i][-1][1], zigzag, numOfSWs])

        i += 1

    headerResults = ['CHROM','QTLStart','QTLEnd','Peaks', 'NumOfSWs']
    pd.DataFrame(snpRegion, columns=headerResults).to_csv(os.path.join(results, 'snpRegion.csv'), index=False)

    swDataFrame = pd.DataFrame(swRows, columns=['CHROM', 'sw_Str', fbID+'.AvgLD', sbID+'.AvgLD', 'sSNP', 'totalSNP', r'sSNP/totalSNP', 'GS', 'GS_CI0995', 'Delta_AF', 'DAF_CI0005', 'DAF_CI0995'])
    swDataFrame['smthedRatio'] = sg_yRatio_List

    swDataFrame.to_csv(os.path.join(results, 'slidingWindows_WoP.csv'), index=False)

    misc.append(['List of the peaks of the chromosomes', ratioPeakL])

    wrnLog = os.path.join(results, 'wrnLog.csv')
    with open(wrnLog, 'w', newline='') as outF1:
        xie1 = csv.writer(outF1)
        xie1.writerow(['Type', 'Chr', 'Position', 'Warning Message'])
        xie1.writerows(wmL)

    numOfSNPOnChrFile = os.path.join(results, 'numOfSNPOnChrFile_WoP.csv')
    with open(numOfSNPOnChrFile, 'w', newline='') as outF2:
        xie2 = csv.writer(outF2)
        xie2.writerow(['Chromosome', 'Num of sSNPs', 'Num of totalSNPs', r'sSNP/totalSNP'])
        xie2.writerows(numOfSNPOnChr)

    print(f'Plotting completed, time elapsed: {(time.time()-t0)/60} minutes.')


def peak(l):
    ratioList = []
    for subL in l:
        ratioList.append(subL[5])

    return l[ratioList.index(max(ratioList))]


def pkList(l):
    # l: genomic regions above the threshold identified in the plot function
    peakList = []
    for subL in l:
        if subL[3] != [] and subL[4] > 10:
            tempL = peak(subL[3])
            peakList.append([subL[0], tempL[0]])

    return peakList


def accurateThreshold_sw(l, df):
    peaks = []
    for subL in l:
        peakSW = df[(df.CHROM == subL[0]) & (df.POS >= subL[1]) & (df.POS <= subL[1]+swSize-1)]
        sSNP_PeakSW = peakSW[peakSW.FE_P<alpha]

        sSNP, totalSNP = len(sSNP_PeakSW.index), len(peakSW.index)
        ratio = sSNP / totalSNP

        peaks.append([subL[0], subL[1], int(peakSW[fb_LD].mean()), int(peakSW[sb_LD].mean()), sSNP, totalSNP, ratio, smThresholds_sw(peakSW)[1]])

    headerResults = ['CHROM','sw_Str', fbID+'.AvgLD', sbID+'.AvgLD', 'sSNP', 'totalSNP', r'sSNP/totalSNP', 'Threshold']
    pd.DataFrame(peaks, columns=headerResults).to_csv(os.path.join(results, args.output), index=False)


def accurateThreshold_gw(l, df):
    peaks = []
    for subL in l:
        peakSW = df[(df.CHROM == subL[0]) & (df.POS >= subL[1]) & (df.POS <= subL[1]+swSize-1)]
        sSNP_PeakSW = peakSW[peakSW.FE_P<alpha]

        sSNP, totalSNP = len(sSNP_PeakSW.index), len(peakSW.index)
        ratio = sSNP / totalSNP

        peaks.append([subL[0], subL[1], int(peakSW[fb_LD].mean()), int(peakSW[sb_LD].mean()), sSNP, totalSNP, ratio, thrshld])

    headerResults = ['CHROM','sw_Str', fbID+'.AvgLD', sbID+'.AvgLD', 'sSNP', 'totalSNP', r'sSNP/totalSNP', 'Threshold']
    pd.DataFrame(peaks, columns=headerResults).to_csv(os.path.join(results, args.output), index=False)


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

# Construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument('-i', '--input', required=False, help='file name of the GATK4-generated tsv file', default='snp100SE.tsv')
ap.add_argument('-o', '--output', required=False, help='file name of the output csv file', default='BSASeq_WoP.csv')
ap.add_argument('-b', '--bulksizes', required=False, help='bulk sizes: first_bulk,second_bulk', type=lambda s: [int(t) for t in s.split(',')], default='430,385')
ap.add_argument('-p', '--popstrct', required=False, choices=['F2','RIL','BC'], help='population structure', default='F2')
ap.add_argument('-v', '--pvalues', required=False, help='cutoff p-values: real_data,simulation', type=lambda s: [float(t) for t in s.split(',')], default='0.01,0.1')
ap.add_argument('-r', '--replication', type=int, required=False, help='the number of replications for threshold calculation', default=10000)
ap.add_argument('-s', '--slidingwindow', required=False, help='size,incremental_step', type=lambda s: [int(t) for t in s.split(',')], default='2000000,10000')
ap.add_argument('-g', '--gaps', required=False, help='gaps between subplots: horizontal,vertical', type=lambda s: [float(t) for t in s.split(',')], default='0.028,0.056')
ap.add_argument('-m', '--smoothing', required=False, help='smoothing parameters: window_len,polyorder', type=lambda s: [int(t) for t in s.split(',')], default='51,3')
ap.add_argument('--smooth', type=bool, required=False, help='smooth the plot', default=False)
ap.add_argument('-e', '--region', required=False, help='interested region(s): chrm,start,end', type=lambda s: [int(t) for t in s.split(',')], default='-1')
ap.add_argument('-c', '--misc', required=False, help='cut-off GQ value, minimum SNPs in a sliding window, extremely high read, and mirror index of Δ(allele frequency)', type=lambda s: [int(t) for t in s.split(',')], default='20,5,6,-1')

args = ap.parse_args()

popStr = args.popstrct
rep = args.replication
fb_Size, sb_Size = args.bulksizes[0], args.bulksizes[1]
alpha, smAlpha = args.pvalues[0], args.pvalues[1]
swSize, incrementalStep = args.slidingwindow[0], args.slidingwindow[1]
hGap, wGap = args.gaps[0], args.gaps[1]
smoothing = args.smooth
smthWL, polyOrder = args.smoothing[0], args.smoothing[1]
region = args.region
gqValue, minSNPs, hiReadFold, mrIndex = args.misc[0], args.misc[1], args.misc[2], args.misc[3]
minFragSize = swSize + smthWL * incrementalStep     # Minimum chromosome size allowed
additionalPeaks = ''
percentile_list = [0.5, 99.5, 2.5, 97.5, 5.0, 95.0]

# Obtain the ALT allele frequencies in both bulks under the null hypothesis
fb_Freq = smAlleleFreq(popStr, fb_Size, rep)
sb_Freq = smAlleleFreq(popStr, sb_Size, rep)

# Set paths for input/output files
path = os.getcwd()
inFile, oiFile = os.path.join(path, args.input), os.path.join(path, 'Results', 'snp_fagz_WoP.csv')
currentDT = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
results = os.path.join(path, 'Results', currentDT)
filteringPath = os.path.join(path, 'Results', 'FilteredSNPs')

# Create folders if not exist
if not os.path.exists(results):
    os.makedirs(results)

if not os.path.exists(filteringPath):
    os.makedirs(filteringPath)

if inFile.endswith('.tsv'):
    separator = '\t'
elif inFile.endswith('.csv'):
    separator = ','

# Generte a SNP dataframe from the GATK4-generated tsv file
snpRawDF = pd.read_csv(inFile, delimiter=separator, encoding='utf-8', dtype={'CHROM':str})

# Create a chromosome list, which can be very long because of the unmapped fragments
chrmRawList = snpRawDF['CHROM'].unique().tolist()

# Filter out chromosomes and unmapped fragments smaller than the sliding window
# Make the chromosome list more readable and meaningful
chrmCheck = chrmFiltering(snpRawDF, chrmRawList)
chrmSizeDict = chrmCheck[0]
chrmList = chrmCheck[1]
smallFrags = chrmCheck[2]

# Print the chromosome list on the screen to let the user to select desired chromosome(s)
rightInput = ''
if region[0] == -1:
    if smallFrags != []:
        print('\nThe chromosomes/fragments below are filtered out because of their small sizes:')
        print(smallFrags, '\n')

    print(f'The chromosomes below are greater than {minFragSize} bp and are suitable for BSA-Seq analysis:')
    print(chrmList,'\n')
    print('Although a subset of the above chromosomes can be selected for analysis in next step, it is strongly recommended to have all the chromosomes included when run the script the first time.\n')

    while rightInput.lower() != 'yes':
        inputString = input('Enter the names of the desired chromosomes IN ORDER and separate each name with a comma or press the ENTER/RETURN key if the above list is what you want:\n')
        if inputString == '':
            chrmIDL = chrmList
        else:
            chrmIDL = [x.strip() for x in inputString.split(',')]

        print('Invalid chromosome names, if any, will be removed.\n')

        # Filter out possible invalid chromosome name(s)
        invalidChrmL = []
        for ch in chrmIDL:
            if ch not in chrmList:
                invalidChrmL.append(ch)

        if invalidChrmL != []:
            print('These chromosome names are invalid:', invalidChrmL,'\n')

        for ch in invalidChrmL:
            chrmIDL.remove(ch)

        print('Sorted chromosome list:')
        print(chrmIDL,'\n')

        rightInput = input('Are the chromosome names in the above list in the right order (yes or no)?\n')
        print('\n')

    # Create dictionary/list containing the sizes of all the chromosomes
    chrmSzD, chrmSzL = {}, []
    for ch in chrmIDL:
        chrmSzD[ch] = chrmSizeDict[ch]
        chrmSzL.append(chrmSizeDict[ch][1])

# Handle the cases in which a chromosomal region is entered
else:
    chrmSzD, chrmIDL, chrmSzL = {}, [], []

    l, n = len(region), len(region) % 3
    if n != 0:
        region = region[0:l-n]
        print('Each chromosome name should be followed by the starting and ending points of the interested region.')

    i = 0
    while i < len(region):
        if region[i] == 1000:
            chID = 'X'
        elif region[i] == 1001:
            chID = 'Y'
        elif region[i] == 1002:
            chID = 'Z'
        elif region[i] == 1003:
            chID = 'W'
        elif region[i] == 1004:
            chID = 'U'
        elif region[i] == 1005:
            chID = 'V'
        else:
            chID = str(region[i])

        if chID not in chrmList:
            print(chID, 'is not a valid chromosome name.')
        else:
            if region[i+1] < 1:
                region[i+1] = 1
            if region[i+2] > chrmSizeDict[chID][1]:
                region[i+2] = chrmSizeDict[chID][1]

            if region[i+2] - region[i+1] > minFragSize:
                chrmSzD[chID] = [region[i+1], region[i+2]]
                chrmIDL.append(chID)
                chrmSzL.append(region[i+2]-region[i+1])
            else:
                print(f'The size of the interested region on chromosome {chID} should be greater than', minFragSize, 'bp.')

        i += 3

# Exit if chrmIDL is empty
if chrmIDL == []:
    print('No valid chromosome names were entered.')
    sys.exit()

xt_pro = xticks_property(chrmSzL)

# Create a numeric ID for each chromosome, which can be used to sort the dataframe numerically by chromosome
chrmDict = {}
for i in range(1, len(chrmIDL)+1):
    chrmDict[chrmIDL[i-1]] = i

# Replace the string IDs with numeric IDs 
snpRawDF['ChrmSortID'] = snpRawDF['CHROM']
snpRawDF['ChrmSortID'].replace(chrmDict, inplace=True)

# Obtain the bulk IDs from the header of the DataFrame
header = snpRawDF.columns.values.tolist()
bulks, misc, missingFlds = [], [], False
# Extract the bulk IDs from the AD fields
try:
    for ftrName in header:
        if ftrName.endswith('.AD'):
            bulks.append(ftrName.split('.')[0])

    fbID, sbID = bulks[0], bulks[1]
    fb_GT, sb_GT = fbID+'.GT', sbID+'.GT'
    fb_AD, sb_AD = fbID+'.AD', sbID+'.AD'
    fb_GQ, sb_GQ = fbID+'.GQ', sbID+'.GQ'
except (NameError, IndexError):
    print('The allele depth (AD) field is missing. Please include the AD field in the input file.')
    sys.exit()

# Check if any required field is missing in the input file
requiredFields = ['CHROM', 'POS', 'REF', 'ALT', fb_GT, fb_AD, fb_GQ, sb_GT, sb_AD, sb_GQ]
missingFields = []
for elmt in requiredFields:
    if elmt not in header:
        missingFields.append(elmt)
if missingFields !=[]:
    if len(missingFields) == 1:
        print('The following required field is missing: ', missingFields)
    else:
        print('The following required fields are missing: ', missingFields)

    print('Please remake the input file to include the missing field(s).')
    sys.exit()

misc.append(['Header', header])
misc.extend([['Bulk ID', bulks], ['Number of SNPs in the entire dataframe', len(snpRawDF.index)]])
misc.extend([['Chromosome ID', chrmIDL]])
misc.append(['Chromosome sizes', chrmSzL])

fb_AF, sb_AF = fbID+'.AF', sbID+'.AF'
fb_AD_REF, fb_AD_ALT = fb_AD + '_REF', fb_AD + '_ALT'
sb_AD_REF, sb_AD_ALT = sb_AD + '_REF', sb_AD + '_ALT'
fb_GT_REF, fb_GT_ALT = fb_GT + '_REF', fb_GT + '_ALT'
sb_GT_REF, sb_GT_ALT = sb_GT + '_REF', sb_GT + '_ALT'
fb_LD, sb_LD = fbID+'.LD', sbID+'.LD'
sm_fb_AD_REF, sm_fb_AD_ALT = 'sm_'+fb_AD_REF, 'sm_'+fb_AD_ALT
sm_sb_AD_REF, sm_sb_AD_ALT = 'sm_'+sb_AD_REF, 'sm_'+sb_AD_ALT
fb_AF_CI, sb_AF_CI = fbID+'.AF_CI', sbID+'.AF_CI'
fb_AF_CI0995, sb_AF_CI0995 = fb_AF_CI+'0995', sb_AF_CI+'0995'
fb_AF_CI0005, sb_AF_CI0005 = fb_AF_CI+'0005', sb_AF_CI+'0005'

if os.path.isfile(os.path.join(path, 'Results', 'COMPLETE.txt')) == False:
    bsaSNPs = snpFiltering(snpRawDF)

    if len(bsaSNPs.index) == 0:
        print('The dataset is empty')
        sys.exit()

    bsaSNPs.to_csv(os.path.join(filteringPath, 'bsaSNPs_WoP.csv'), index=None)

    # # Remove low LD SNPs to meet the z-test sample size requirement
    # bsaSNPs = bsaSNPs[(bsaSNPs[fb_LD]>=30) & (bsaSNPs[sb_LD]>=30)]

    # Calculate z scores
    bsaSNPs['fb_ZScore'] = (bsaSNPs[fb_AD_ALT] / bsaSNPs[fb_LD] - fb_Freq) * np.sqrt(bsaSNPs[fb_LD] / (fb_Freq * (1 - fb_Freq)))
    bsaSNPs['sb_ZScore'] = (bsaSNPs[sb_AD_ALT] / bsaSNPs[sb_LD] - sb_Freq) * np.sqrt(bsaSNPs[sb_LD] / (sb_Freq * (1 - fb_Freq)))

    # Calculate simulated ALT reads for each SNP under null hypothesis
    bsaSNPs[sm_fb_AD_ALT] = np.random.binomial(bsaSNPs[fb_LD], fb_Freq)
    bsaSNPs[sm_fb_AD_REF] = bsaSNPs[fb_LD] - bsaSNPs[sm_fb_AD_ALT]
    bsaSNPs[sm_sb_AD_ALT] = np.random.binomial(bsaSNPs[sb_LD], sb_Freq)
    bsaSNPs[sm_sb_AD_REF] = bsaSNPs[sb_LD] - bsaSNPs[sm_sb_AD_ALT]

    # Calculate allele frequency and Δ(allele frequency)
    bsaSNPs[fb_AF] = bsaSNPs[fb_AD_ALT]/bsaSNPs[fb_LD]
    bsaSNPs[sb_AF] = bsaSNPs[sb_AD_ALT]/bsaSNPs[sb_LD]
    bsaSNPs['Delta_AF'] = bsaSNPs[sb_AF] - bsaSNPs[fb_AF]

    # Calculate G-statistic
    bsaSNPs['G_S'] = gStatistic_Array(bsaSNPs[fb_AD_REF], bsaSNPs[fb_AD_ALT], bsaSNPs[sb_AD_REF], bsaSNPs[sb_AD_ALT])

    try:
        from fisher import pvalue_npy
        # Convert the AD columns to arrays
        print('Perform Fisher\'s exact test.')
        fb_AD_ALT_Arr = bsaSNPs[fb_AD_ALT].to_numpy(dtype=np.uint)
        fb_AD_REF_Arr = bsaSNPs[fb_AD_REF].to_numpy(dtype=np.uint)
        sb_AD_ALT_Arr = bsaSNPs[sb_AD_ALT].to_numpy(dtype=np.uint)
        sb_AD_REF_Arr = bsaSNPs[sb_AD_REF].to_numpy(dtype=np.uint)
        # Calculate the P-values via Fisher's exact test
        __, __, bsaSNPs['FE_P'] = pvalue_npy(fb_AD_ALT_Arr, fb_AD_REF_Arr, sb_AD_ALT_Arr, sb_AD_REF_Arr)
        # bsaSNPs['FE_OR'] = (fb_AD_ALT_Arr * sb_AD_REF_Arr) / (fb_AD_REF_Arr * sb_AD_ALT_Arr)
        # Convert the simulated AD columns to arrays
        sm_fb_AD_ALT_Arr = bsaSNPs[sm_fb_AD_ALT].to_numpy(dtype=np.uint)
        sm_fb_AD_REF_Arr = bsaSNPs[sm_fb_AD_REF].to_numpy(dtype=np.uint)
        sm_sb_AD_ALT_Arr = bsaSNPs[sm_sb_AD_ALT].to_numpy(dtype=np.uint)
        sm_sb_AD_REF_Arr = bsaSNPs[sm_sb_AD_REF].to_numpy(dtype=np.uint)
        # Calculate the P-values using the simulated ADs
        __, __, bsaSNPs['sm_FE_P'] = pvalue_npy(sm_fb_AD_ALT_Arr, sm_fb_AD_REF_Arr, sm_sb_AD_ALT_Arr, sm_sb_AD_REF_Arr)
        # bsaSNPs['sm_FE_OR'] = (sm_fb_AD_ALT_Arr * sm_sb_AD_REF_Arr) / (sm_fb_AD_REF_Arr * sm_sb_AD_ALT_Arr)

        print(f'Fisher\'s exact test completed, time elapsed: {(time.time()-t0)/60} minutes.')
        print('Calculate thresholds of \u0394(allele frequency) and G-statistic.')
        # Obtain the threshold of G-statistic and Δ(allele frequency)
        bsaSNPs['STAT'] = bsaSNPs.apply(statistics, axis=1)

        # Create new columns for Fisher's exact test results, allele frequency, Δ(allele frequency) confidence intervals, and G-statistic thresholds
        bsaSNPs[[fb_AF_CI, sb_AF_CI, 'DAF_CI', 'GS_CI']] = pd.DataFrame(bsaSNPs.STAT.values.tolist(), index=bsaSNPs.index)

        print(f'Calculating thresholds of \u0394(allele frequency) and G-statistic completed, time elapsed: {(time.time()-t0)/60} minutes.')

    except ImportError:
        from scipy.stats import fisher_exact
        print('Looks like the module \'Fisher\' (https://github.com/brentp/fishers_exact_test) is not installed on your computer. Running the script would take much shorter time with this module installed.')
        # Obtain the FE test P-values, the simulate P-values, and the thresholds of G-statistic and Δ(allele frequency)
        bsaSNPs['STAT'] = bsaSNPs.apply(statisticsFE, axis=1)

        # Create new columns for Fisher's exact test results, allele frequency, Δ(allele frequency) confidence intervals, and G-statistic thresholds
        bsaSNPs[['fisher_exact', 'sm_FE', fb_AF_CI, sb_AF_CI, 'DAF_CI', 'GS_CI']] = pd.DataFrame(bsaSNPs.STAT.values.tolist(), index=bsaSNPs.index)

        # Create new columns for Fisher's exact test P-values or simulated P-values
        bsaSNPs['FE_P'] = bsaSNPs['fisher_exact'].apply(lambda x: x[1]).astype(float)
        bsaSNPs['sm_FE_P'] = bsaSNPs['sm_FE'].apply(lambda x: x[1]).astype(float)

        print(f'Fisher\'s exact test and calculating thresholds of \u0394(allele frequency) and G-statistic completed, time elapsed: {(time.time()-t0)/60} minutes.')

    # Create new columns for 99% Δ(allele frequency) confidence intervals, and 99.5 percentile G-statistic thresholds
    bsaSNPs[fb_AF_CI0005] = bsaSNPs[fb_AF_CI].apply(lambda x: x[0]).astype(float)
    bsaSNPs[fb_AF_CI0995] = bsaSNPs[fb_AF_CI].apply(lambda x: x[1]).astype(float)
    bsaSNPs[sb_AF_CI0005] = bsaSNPs[sb_AF_CI].apply(lambda x: x[0]).astype(float)
    bsaSNPs[sb_AF_CI0995] = bsaSNPs[sb_AF_CI].apply(lambda x: x[1]).astype(float)
    bsaSNPs['DAF_CI0005'] = bsaSNPs['DAF_CI'].apply(lambda x: x[0]).astype(float)
    bsaSNPs['DAF_CI0995'] = bsaSNPs['DAF_CI'].apply(lambda x: x[1]).astype(float)
    bsaSNPs['GS_CI0995'] = bsaSNPs['GS_CI'].apply(lambda x: x[1]).astype(float)

    # To make the mirror image of the Δ(allele frequency) curve for comparison if necessary
    if mrIndex == -1:
        bsaSNPs['Delta_AF'] = bsaSNPs['Delta_AF'] * -1
        bsaSNPs['DAF_CI0995'] = bsaSNPs['DAF_CI0995'] * -1
        bsaSNPs['DAF_CI0005'] = bsaSNPs['DAF_CI0005'] * -1

    # Reorganize the columns
    reorderColumns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', fb_GT, fb_AD, fb_AD_REF, fb_AD_ALT, fb_LD, 'fb_ZScore', sm_fb_AD_ALT, fb_AF, fb_AF_CI0005, fb_AF_CI0995, fb_AF_CI, fb_GQ, sb_GT, sb_AD, sb_AD_REF, sb_AD_ALT, sb_LD, 'sb_ZScore', sm_sb_AD_ALT, sb_AF, sb_AF_CI0005, sb_AF_CI0995, sb_AF_CI, sb_GQ, 'Delta_AF', 'FE_P', 'sm_FE_P', 'G_S', 'DAF_CI0005', 'DAF_CI0995', 'GS_CI0995', 'DAF_CI', 'GS_CI', 'STAT']
    # Remove unnecessary columns and reorganize the columns
    # reorderColumns = ['CHROM', 'POS', 'REF', 'ALT', fb_GT, fb_AD_REF, fb_AD_ALT, fb_LD, 'fb_ZScore', sm_fb_AD_ALT, fb_GQ, sb_GT, sb_AD_REF, sb_AD_ALT, sb_LD, 'sb_ZScore', sm_sb_AD_ALT, sb_GQ, 'FE_P', 'sm_FE_P']
    bsaSNPs = bsaSNPs[reorderColumns]
    bsaSNPs.to_csv(oiFile, index=None)

    with open(os.path.join(path, 'Results', 'COMPLETE.txt'), 'w') as xie:
        xie.write('Statistical calculation is completed!')
else:
    additionalPeaks = input('Do you want to have additional peaks identified (yes or no)?\n')
    print('\n')
    ttlSNPs = pd.read_csv(oiFile, dtype={'CHROM':str})
    bsaSNPs = ttlSNPs[ttlSNPs.CHROM.isin(chrmIDL)]

# falseSNPs = bsaSNPs[((bsaSNPs[fb_AF] < bsaSNPs[fb_AF_CI0005]) & (bsaSNPs[sb_AF] < bsaSNPs[sb_AF_CI0005])) | \
#     ((bsaSNPs[fb_AF] > bsaSNPs[fb_AF_CI0995]) & (bsaSNPs[sb_AF] > bsaSNPs[sb_AF_CI0995]))]

# bsaSNPs = bsaSNPs.drop(index=falseSNPs.index)

# falseSNPs.to_csv(os.path.join(filteringPath, 'falseSNPs.csv'), index=None)
# bsaSNPs.to_csv(os.path.join(filteringPath, 'SNPs.csv'), index=None)

# The above calculation may generate 'NA' value(s) for some SNPs. Remove SNPs with such 'NA' value(s)
bsaSNPs.dropna(inplace=True)
misc.append(['Number of SNPs after drop of SNPs with calculation-generated NA value', len(bsaSNPs.index)])

# A SNP with its absolute z score value greater than 2.575 is a significant SNP
# bsaSNPs['fb_Z'] = np.where(bsaSNPs['fb_ZScore'].abs() > 2.575, 1, 0)
# bsaSNPs['sb_Z'] = np.where(bsaSNPs['sb_ZScore'].abs() > 2.575, 1, 0)

misc.append(['Dataframe filtered with genotype quality scores', len(bsaSNPs.index)])

# Calculate the average number of SNPs in a sliding window
snpPerSW = int(len(bsaSNPs.index) * swSize / sum(chrmSzL))

misc.append(['Average SNPs per sliding window', snpPerSW])
misc.append([f'Average locus depth in bulk {fbID}', bsaSNPs[fb_LD].mean()])
misc.append([f'Average locus depth in bulk {sbID}', bsaSNPs[sb_LD].mean()])

# Calculate or retrieve the threshold. The thresholds are normally in the range from 0.12 to 0.12666668
if os.path.isfile(os.path.join(path, 'Results', 'threshold.txt')) == False:
    if 'fisher' in sys.modules:
        thrshld = smThresholds_gw(bsaSNPs)[1]
    else:
        thrshld = smThresholds_proximal(bsaSNPs)[1]

    with open(os.path.join(path, 'Results', 'threshold.txt'), 'w') as xie:
        xie.write(str(thrshld))
else:
    with open(os.path.join(path, 'Results', 'threshold.txt'), 'r') as du:
        thrshld = float(du.readline().strip())

# Identify likely trait-associated SNPs
bsaSNPs['sSNP'] = np.where(bsaSNPs['FE_P'] < alpha, 1, 0)

# Plot layout setup
heightRatio = [1,0.8,0.8,0.8]
fig, axs = plt.subplots(nrows=len(heightRatio), ncols=len(chrmIDL), figsize=(20, 12.8), sharex='col', sharey='row', 
        gridspec_kw={'width_ratios': chrmSzL, 'height_ratios': heightRatio})

# Perform plotting
bsaseqPlot(chrmIDL, bsaSNPs)

if 'fisher' not in sys.modules:
    try:
        from fisher import pvalue_npy
    except ImportError:
        print('The module \'Fisher\' is not installed on your computer.')

# Handle the plot with a single column (chromosome)
if len(chrmIDL) == 1:
    fig.align_ylabels(axs[:])
# Handle the plot with multiple columns (chromosomes)
else:
    fig.align_ylabels(axs[:, 0])

# fig.tight_layout(pad=0.15, rect=[0, 0.035, 1, 1])
fig.subplots_adjust(top=0.9736, bottom=0.054, left=0.076, right=0.998, hspace=hGap, wspace=wGap)
fig.suptitle(f'Genomic position ({xt_pro[2]})', y=0.002, ha='center', va='bottom')
fig.text(0.001, 0.995, 'a', weight='bold', ha='left', va='top')
fig.text(0.001, 0.680, 'b', weight='bold', ha='left', va='bottom')
fig.text(0.001, 0.465, 'c', weight='bold', ha='left', va='bottom')
fig.text(0.001, 0.244, 'd', weight='bold', ha='left', va='bottom')

fig.savefig(os.path.join(results, 'PyBSASeq_WoP.pdf'))
fig.savefig(os.path.join(results, 'PyBSASeq_WoP.png'), dpi=600)

# Create a list of (CHROM, peak position)
peaklst = pkList(snpRegion)

if additionalPeaks.lower() == 'yes':
    with open(os.path.join(path, 'additionalPeaks.txt'), 'r') as inF:
        for line in inF:
            if not line.startswith('#'):
                a = line.rstrip().split()

                if a[0] in chrmIDL:
                    additionalSW = swDataFrame[(swDataFrame.CHROM==a[0]) & (swDataFrame.sw_Str >= int(a[1])) & (swDataFrame.sw_Str <= int(a[2]))]
                    peakSW = additionalSW[additionalSW[r'sSNP/totalSNP'] == additionalSW[r'sSNP/totalSNP'].max()]

                for row in peakSW.itertuples():
                    peaklst.append([row.CHROM, row.sw_Str])

peaklst = sorted(peaklst, key = lambda x: (int(x[0]), int(x[1])))

try:
    accurateThreshold_sw(peaklst, bsaSNPs)
    print(f'Peak verification completed, time elapsed: {(time.time()-t0)/60} minutes.')
except NameError:
    accurateThreshold_gw(peaklst, bsaSNPs)
    print('Please install the module \'Fisher\' if more precise thresholds of the QTL loci are desired.')

misc.append(['Running time', [(time.time()-t0)/60]])

with open(os.path.join(results, 'misc_info.csv'), 'w', newline='') as outF:
    xie = csv.writer(outF)
    xie.writerows(misc)

print('\nIf two or more peaks and all the values in between are greater than the threshold, these peaks would be recognized as a single peak. You can rerun the script if additional peaks are desired to be identified, a file \'additionalPeaks.txt\' (template provided) containing the chromosome ID and the range info (start and end) of the interested regions needs to be created in the working directory.\n')