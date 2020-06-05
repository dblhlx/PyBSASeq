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
    An AA/Aa/aa individual carries 0%, 50%, and 100% of the alt (a) allele, respectively.
    The AA:Aa:aa ratios are 0.25:0.5:0.25, 0.5:0:0.5, and 0.5:0.5:0, repectively, in a F2 
    population, in a RIL population, and in a back crossed population if A/a does not affect 
    the trait in the population (null hypothesis)
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


def snpFiltering(df):
    print('Perform SNP filtering')
    global misc

    df = df.copy()

    # Identify SNPs not interested and remove these SNPs from the dataframe
    df_Ignored = df[~df.CHROM.isin(chrmIDL)]
    df_Ignored.to_csv(os.path.join(filteringPath, 'ignored.csv'), index=None)
    df = df.drop(index=df_Ignored.index)

    # Identify SNPs with an 'NA' value(s) and remove these SNPs from the dataframe
    df_NA = df[df.isnull().any(axis=1)]
    df_NA.to_csv(os.path.join(filteringPath, 'na.csv'), index=None)
    df.dropna(inplace=True)
    misc.append(['Number of SNPs after NA drop', len(df.index)])

    # Filter out SNPs with a low genotype quality score
    df_lowq = df[(df[fb_GQ]<20) | (df[sb_GQ]<20)]
    df_lowq.to_csv(os.path.join(filteringPath, 'lowqualitySNP.csv'), index=None)
    df = df.drop(index=df_lowq.index)

    # Fiter out InDels
    inDel_1 = df[df['ALT'].str.contains(r'\*')]
    df = df.drop(index=inDel_1.index)

    # Identify SNPs with a single ALT allele
    df_1ALT = df[~df['ALT'].str.contains(',')]

    # Identify SNPs with more than one ALT allele
    df_mALT = df.drop(index=df_1ALT.index)

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

    # The two-ALT SNP may be cuased by allele heterozygosity if the REF read is not zero
    # Repetitive sequences in the genome or sequencing artifacts are other possibilities
    df_2ALT_Het = df_2ALT.drop(index=df_2ALT_Real.index)
    df_2ALT_Het.to_csv(os.path.join(filteringPath, '2ALT_Het_Before.csv'), index=None)

    # Making a copy to suppress the warning message. Updating the AD values of these SNPs is required
    df_2ALT_Real = df_2ALT_Real.copy()
    df_2ALT_Real.to_csv(os.path.join(filteringPath, '2altReal_Before.csv'), index=None)

    # Update the AD values of the above SNPs by removing the REF read that is zero (i.e. remove '0,' from '0,x,y')
    df_2ALT_Real[fb_AD] = df_2ALT_Real[fb_AD].str.slice(start=2)
    df_2ALT_Real[sb_AD] = df_2ALT_Real[sb_AD].str.slice(start=2)
    df_2ALT_Real.to_csv(os.path.join(filteringPath, '2altReal_After.csv'), index=None)

    df_2ALT_Het.to_csv(os.path.join(filteringPath, '2ALT_Het_After.csv'), index=None)

    # Concatenate 1ALT_Real and 2ALT_Real
    snp = pd.concat([df_1ALT_Real, df_2ALT_Real])

    # In case the input file contains Indels
    try:
        # Identify inDels, REF/Alt allele with more than 1 base
        inDel_2 = snp[(snp['REF'].str.len()>1) | (snp['ALT'].str.split(',', expand=True)[0].str.len()>1) | \
        (snp['ALT'].str.split(',', expand=True)[1].str.len()>1)]

    except KeyError:
        print('All the ALT loci contain only one allele.')
        # Identify inDels, REF/Alt allele with more than 1 base
        inDel_2 = snp[(snp['REF'].str.len()>1) | (snp['ALT'].str.len()>1)]

    df_InDel = pd.concat([inDel_1, inDel_2])
    df_InDel.to_csv(os.path.join(filteringPath, 'InDel.csv'), index=None)

    # Remove InDels from the SNP dataframe (inDel_1 is already removed)
    snp = snp.drop(index=inDel_2.index)

    snp.sort_values(['ChrmSortID', 'POS'], inplace=True)

    # Obtain REF reads, ALT reads, and locus reads of each SNP
    snp[[fb_AD_REF, fb_AD_ALT]] = snp[fb_AD].str.split(',', expand=True).astype(int)
    snp[fb_LD] = snp[fb_AD_REF] + snp[fb_AD_ALT]
    snp[[sb_AD_REF, sb_AD_ALT]] = snp[sb_AD].str.split(',', expand=True).astype(int)
    snp[sb_LD] = snp[sb_AD_REF] + snp[sb_AD_ALT]

    # A SNP with its LD 6 times higher than the average LD will be eliminated from the dataset 
    fb_repLD, sb_repLD = 6 * snp[fb_LD].mean(), 6 * snp[sb_LD].mean()

    # Remove SNPs that could be from the repetitive elements 
    snpRep = snp[(snp[fb_LD]>fb_repLD) | (snp[sb_LD]>sb_repLD)]
    snpRep.to_csv(os.path.join(filteringPath, 'repetitiveSeq.csv'), index=None)

    snp = snp[(snp[fb_LD]<=fb_repLD) & (snp[sb_LD]<=sb_repLD)]

    # Filter out the SNPs with zero locus reads in either bulk
    snp_0LD = snp[(snp[fb_LD] == 0) | (snp[sb_LD] == 0)]
    snp_0LD.to_csv(os.path.join(filteringPath, '0ld.csv'), index=None)
    snp = snp.drop(index=snp_0LD.index)

    print(f'SNP filtering completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return snp


def gStatistic_Array(o1, o3, o2, o4):
    # Calculate G-statistc using numpy array input
    np.seterr(all='ignore')

    e1 = np.where(o1+o2+o3+o4!=0, (o1+o2)*(o1+o3)/(o1+o2+o3+o4), 0)
    e2 = np.where(o1+o2+o3+o4!=0, (o1+o2)*(o2+o4)/(o1+o2+o3+o4), 0)
    e3 = np.where(o1+o2+o3+o4!=0, (o3+o4)*(o1+o3)/(o1+o2+o3+o4), 0)
    e4 = np.where(o1+o2+o3+o4!=0, (o3+o4)*(o2+o4)/(o1+o2+o3+o4), 0)

    llr1 = np.where(o1/e1>0, 2*o1*np.log(o1/e1), 0.0)
    llr2 = np.where(o2/e2>0, 2*o2*np.log(o2/e2), 0.0)
    llr3 = np.where(o3/e3>0, 2*o3*np.log(o3/e3), 0.0)
    llr4 = np.where(o4/e4>0, 2*o4*np.log(o4/e4), 0.0)

    return np.where(e1*e2*e3*e4==0, 0.0, llr1+llr2+llr3+llr4)


def statisticsFE(row):
    # Perform Fisher's exact test for each SNP using the actual REF/ALT reads
    try:
        fe = fisher_exact([[row[fb_AD_REF], row[fb_AD_ALT]], [row[sb_AD_REF], row[sb_AD_ALT]]])
    except TypeError:
        fe = ('NA', 'NA')

    # Perform Fisher's exact test for each SNP using the simulated REF/ALT reads
    try:
        sm_FE = fisher_exact([[row[sm_fb_AD_REF], row[sm_fb_AD_ALT]], [row[sm_sb_AD_REF], row[sm_sb_AD_ALT]]])
    except TypeError:
        sm_FE = ('NA', 'NA')

    # Create an array with 10000 (rep) simulated ALT reads of a SNP - first bulk
    sm_yb_Alt_Array = np.random.binomial(row[fb_LD], fb_Freq, rep)
    yb_LD_Array = np.full(rep, row[fb_LD])
    # Create an array with 10000 (rep) simulated ALT reads of a SNP - second bulk
    sm_eb_Alt_Array = np.random.binomial(row[sb_LD], sb_Freq, rep)
    eb_LD_Array = np.full(rep, row[sb_LD])

    # Create arraies of SNP indices and Δ(allele frequency) of a SNP
    sm_yb_AF_Array = sm_yb_Alt_Array/yb_LD_Array
    sm_eb_AF_Array = sm_eb_Alt_Array/eb_LD_Array
    sm_DAF_Array = sm_eb_AF_Array - sm_yb_AF_Array

    # Create a G-statistic array of a SNP
    sm_GS_Array = gStatistic_Array(sm_yb_Alt_Array, yb_LD_Array-sm_yb_Alt_Array, sm_eb_Alt_Array, eb_LD_Array-sm_eb_Alt_Array)

    # Obtain the percentile of the above arraies
    ci_yb_AF = np.percentile(sm_yb_AF_Array, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])
    ci_eb_AF = np.percentile(sm_eb_AF_Array, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])
    ci_DAF = np.percentile(sm_DAF_Array, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])
    ci_GS = np.percentile(sm_GS_Array, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])

    return [fe, sm_FE, ci_yb_AF, ci_eb_AF, ci_DAF, ci_GS]


def statistics(row):
    # Create an array with 10000 (rep) simulated ALT reads of a SNP - first bulk
    sm_yb_Alt_Array = np.random.binomial(row[fb_LD], fb_Freq, rep)
    yb_LD_Array = np.full(rep, row[fb_LD])
    # Create an array with 10000 (rep) simulated ALT reads of a SNP - second bulk
    sm_eb_Alt_Array = np.random.binomial(row[sb_LD], sb_Freq, rep)
    eb_LD_Array = np.full(rep, row[sb_LD])

    # Create arraies of SNP indices and Δ(allele frequency) of a SNP
    sm_yb_AF_Array = sm_yb_Alt_Array/yb_LD_Array
    sm_eb_AF_Array = sm_eb_Alt_Array/eb_LD_Array
    sm_DAF_Array = sm_eb_AF_Array - sm_yb_AF_Array

    # Create a G-statistic array of a SNP
    sm_GS_Array = gStatistic_Array(sm_yb_Alt_Array, yb_LD_Array-sm_yb_Alt_Array, sm_eb_Alt_Array, eb_LD_Array-sm_eb_Alt_Array)

    # Obtain the percentile of the above arraies
    ci_yb_AF = np.percentile(sm_yb_AF_Array, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])
    ci_eb_AF = np.percentile(sm_eb_AF_Array, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])
    ci_DAF = np.percentile(sm_DAF_Array, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])
    ci_GS = np.percentile(sm_GS_Array, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])

    return [ci_yb_AF, ci_eb_AF, ci_DAF, ci_GS]


# Using this function for the calculation of the threshold if 'fisher' is not available
def smThresholds_proximal(DF):
    print('Calculate the threshold of sSNPs/totalSNPs.')
    ratioLi = []
    for __ in range(rep):
        sm_SNP_SMPL = DF.sample(snpPerSW, replace=True)
        sm_sSNP_SMPL = sm_SNP_SMPL[sm_SNP_SMPL['sm_FE_P']<smAlpha]

        ratioLi.append(len(sm_sSNP_SMPL.index)/snpPerSW)

    misc.append(['Genome-wide sSNP/totalSNP ratio threshold', np.percentile(ratioLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])])
    print(f'Threshold calculation completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return np.percentile(ratioLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])


def smThresholds_gw(DF):
    # For the calculation of the genome-wide threshold
    print('Calculate the threshold of sSNPs/totalSNPs.')
    gw_ratioLi = []
    for __ in range(rep):
        sm_SNP_SMPL = DF.sample(snpPerSW, replace=True)

        gw_sm_fb_AD_ALT_Arr = np.random.binomial(sm_SNP_SMPL[fb_LD], fb_Freq).astype(np.uint)
        gw_sm_fb_AD_REF_Arr = sm_SNP_SMPL[fb_LD].to_numpy().astype(np.uint) - gw_sm_fb_AD_ALT_Arr
        gw_sm_sb_AD_ALT_Arr = np.random.binomial(sm_SNP_SMPL[sb_LD], sb_Freq).astype(np.uint)
        gw_sm_sb_AD_REF_Arr = sm_SNP_SMPL[sb_LD].to_numpy().astype(np.uint) - gw_sm_sb_AD_ALT_Arr

        __, __, gw_sm_FE_P_Arr = pvalue_npy(gw_sm_fb_AD_ALT_Arr, gw_sm_fb_AD_REF_Arr, gw_sm_sb_AD_ALT_Arr, gw_sm_sb_AD_REF_Arr)
        # gw_sm_FE_OR_Arr = (gw_sm_fb_AD_ALT_Arr * gw_sm_sb_AD_REF_Arr) / (gw_sm_fb_AD_REF_Arr * gw_sm_sb_AD_ALT_Arr)

        sSNP_Arr = np.where(gw_sm_FE_P_Arr<smAlpha, 1, 0)

        gw_ratioLi.append(np.mean(sSNP_Arr))

    misc.append(['Genome-wide sSNP/totalSNP ratio threshold', np.percentile(gw_ratioLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])])
    print(f'Threshold calculation completed, time elapsed: {(time.time()-t0)/60} minutes.')

    return np.percentile(gw_ratioLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])


def smThresholds_sw(DF):
    # For the calculation of the sliding window-specific threshold
    sw_ratioLi = []

    sw_fb_LD_Arr = DF[fb_LD].to_numpy().astype(np.uint)
    sw_sb_LD_Arr = DF[sb_LD].to_numpy().astype(np.uint)

    for __ in range(rep):
        # Create new columns for Fisher's exact test simulated P-values
        sw_sm_fb_AD_ALT_Arr = np.random.binomial(DF[fb_LD], fb_Freq).astype(np.uint)
        sw_sm_fb_AD_REF_Arr = sw_fb_LD_Arr - sw_sm_fb_AD_ALT_Arr
        sw_sm_sb_AD_ALT_Arr = np.random.binomial(DF[sb_LD], sb_Freq).astype(np.uint)
        sw_sm_sb_AD_REF_Arr = sw_sb_LD_Arr - sw_sm_sb_AD_ALT_Arr

        __, __, sw_sm_FE_P_Arr = pvalue_npy(sw_sm_fb_AD_ALT_Arr, sw_sm_fb_AD_REF_Arr, sw_sm_sb_AD_ALT_Arr, sw_sm_sb_AD_REF_Arr)
        # sw_sm_FE_OR_Arr = (sw_sm_fb_AD_ALT_Arr * sw_sm_sb_AD_REF_Arr) / (sw_sm_fb_AD_REF_Arr * sw_sm_sb_AD_ALT_Arr)

        sSNP_Arr = np.where(sw_sm_FE_P_Arr<smAlpha, 1, 0)

        sw_ratioLi.append(np.mean(sSNP_Arr))

    return np.percentile(sw_ratioLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])


def zeroSNP(li):
    # Replace 'divide by zero' with the nearnest value. Use 'empty' as a placeholder if it is the first element of the list
    if li != []:
        li.append(li[-1])   # Assign the previous value to the empty sliding window if the list is not empty
    else:
        li.append('empty')  # Assign 'empty' to the first sliding windows that is empty


def replaceZero(li):
    # Replace the 'empty' placeholders at the begining of the list with the nearnest non-empty value
    i = 0
    while li[i]=='empty':
        i += 1

    j = 0
    while j < i:
        li[j] = li[i]
        j += 1


def bsaseqPlot(chrmIDL, datafr):
    '''
    wmL: list of warning messages
    swDict: a dictionary with the chromosome ID as its keys; the value of each key is a list containing
            the chromosome ID, the sSNP/totalSNP ratio in each sliding window, and the midpoint of
            the sliding window
    '''
    print('Prepare SNP data for plotting via the sliding window algorithm')
    global misc
    global snpRegion, swDataFrame
    sg_yRatio_List = []
    swRows = []
    wmL, swDict, snpRegion = [], {}, []

    # Analyze each chromsome separately
    numOfSNPOnChr, ratioPeakL = [], []
    i = 1
    for chrmID in chrmIDL:
        ch = datafr[datafr.CHROM==chrmID]
        numOfSNPOnChr.append([chrmID, ch['sSNP'].sum(), len(ch.index), ch['sSNP'].sum()/len(ch.index)])

        regStart = chrmSzD[chrmID][0]
        regEnd = chrmSzD[chrmID][1]

        # Sliding window. swStr: the begining of the window; swEnd: the end of the window; icrs: incremental step
        swStr, swEnd, icrs = regStart, regStart+swSize, incrementalStep
        plotSP = regStart
        # x and y are lists, each sliding window represents a single data point
        x, y, yT, yRatio = [], [], [], []
        # fb_ZRatio, sb_ZRatio = [], []
        y5, y6, y7, y8, y9 = [], [], [], [], []
        while swEnd <= regEnd:
            # A single sliding window - a dataframe
            # swDF: sSNPs in a sliding window; swDFT: all SNPs in a sliding window
            swDF = ch[(ch.POS>=swStr) & (ch.POS<=swEnd)]

            rowInSwDF = len(swDF.index)       # number of SNPs in a sliding window
            x.append(swStr)                   # Append the starpoint of a sliding window to x
            y.append(swDF['sSNP'].sum())      # Append number of sSNPs in a sliding window to y
            yT.append(rowInSwDF)              # Append number of totalSNPs in a sliding window to yT

            # len(swDL) or len(swDLT) would be zero if no SNP in a sliding window
            # 'try/exception' cannot catch the EmptyDataError; seems it only works when reading a .csv/.tsv file, not an empty subset of a existing dataframe. DivisionByZero generates 'nan' for a series or an array.
            if rowInSwDF >= 5:
                yRatio.append(swDF['sSNP'].mean())    # Append the ratio of sSNP/totalSNP in a sliding window to yRatio
                # fb_ZRatio.append(swDF['fb_Z'].mean())
                # sb_ZRatio.append(swDF['sb_Z'].mean())
                y5.append(swDF['G_S'].sum()/rowInSwDF)
                y6.append(swDF['GS_CI0995'].sum()/rowInSwDF)
                y7.append(swDF['Delta.AF'].sum()/rowInSwDF)
                y8.append(swDF['DAF_CI0005'].sum()/rowInSwDF)
                y9.append(swDF['DAF_CI0995'].sum()/rowInSwDF)

                rowContents = [chrmID, swStr, int(swDF[fb_LD].mean()), int(swDF[sb_LD].mean()), swDF['sSNP'].sum(), rowInSwDF, yRatio[-1]]
            else:
                wmL.append(['No SNP', i, swStr])
                zeroSNP(yRatio)
                # zeroSNP(fb_ZRatio)
                # zeroSNP(sb_ZRatio)
                zeroSNP(y5)
                # print(chrmID, y5[0])
                zeroSNP(y6)
                zeroSNP(y7)
                zeroSNP(y8)
                zeroSNP(y9)

                # The mean of an empty column is 'nan', not zero, and int(nan) generates a ValueError
                rowContents = [chrmID, swStr, 0, 0, swDF['sSNP'].sum(), rowInSwDF, yRatio[-1]]

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

        pIndex = 0
        while swDict[i][pIndex][6] == 'empty':
            pIndex += 1

        j = 0
        while j < pIndex:
            swDict[i][j][6] = swDict[i][pIndex][6]
            j += 1

        # Data smoothing
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
            axs[0].set_xticks(np.arange(0, max(x), 10000000))
            ticks = axs[0].get_xticks()*1e-7
            axs[0].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0].set_ylabel('Number of SNPs')
                axs[1].set_ylabel(r'sSNP/totalSNP')
                axs[2].set_ylabel('G-statistic')
                axs[3].set_ylabel('\u0394(allele frequency)')

            # SNP plot
            axs[0].plot(x, y, c='k')
            axs[0].plot(x, yT, c='b')
            axs[0].set_title('Chr'+chrmID)

            # sSNP/totalSNP plot
            if smoothing == True:
                # sSNPs/totalSNPs plot via Fisher's exact test
                axs[1].plot(x, sg_yRatio, c='k')

                # G-statistic plot
                axs[2].plot(x, sg_y5, c='k')
                axs[2].plot(x, sg_y6, c='r')

                # Δ(allele frequency) plot
                axs[3].plot(x, sg_y7, c='k')
                axs[3].plot(x, sg_y8, c='r')
                axs[3].plot(x, sg_y9, c='r')

                # sSNPs/totalSNPs plot via z-test
                # axs[1].plot(x, sg_fb_ZRatio, c='c')
                # axs[1].plot(x, sg_sb_ZRatio, c='g')

            else:
                # sSNPs/totalSNPs plot via Fisher's exact test
                axs[1].plot(x, yRatio, c='k')

                # G-statistic plot
                axs[2].plot(x, y5, c='k')
                axs[2].plot(x, y6, c='r')

                # Δ(allele frequency) plot
                axs[3].plot(x, y7, c='k')
                axs[3].plot(x, y8, c='r')
                axs[3].plot(x, y9, c='r')

                # sSNPs/totalSNPs plot via z-test
                # axs[1].plot(x, fb_ZRatio, c='c')
                # axs[1].plot(x, sb_ZRatio, c='g')

            # Add the 99.5 percentile line as threshold, x[-1] is the midpoint of the last sliding window of a chromosome
            axs[1].plot([plotSP, x[-1]], [thrshld, thrshld], c='r')
            # axs[1].plot(x, smThresholds_sw, c='m')

        # Handle the plot with multiple columns (chromosomes)
        else:
            # Set up x-ticks
            axs[0,i-1].set_xticks(np.arange(0, max(x), 10000000))
            ticks = axs[0,i-1].get_xticks()*1e-7
            axs[0,i-1].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0,i-1].set_ylabel('Number of SNPs')
                axs[1,i-1].set_ylabel(r'sSNP/totalSNP')
                axs[2,i-1].set_ylabel('G-statistic')
                axs[3,i-1].set_ylabel('\u0394(allele frequency)')

            # SNP plot
            axs[0,i-1].plot(x, y, c='k')
            axs[0,i-1].plot(x, yT, c='b')
            axs[0,i-1].set_title('Chr'+chrmID)

            # sSNP/totalSNP plot
            if smoothing == True:
                # sSNPs/totalSNPs plot
                axs[1,i-1].plot(x, sg_yRatio, c='k')

                # G-statistic plot
                axs[2,i-1].plot(x, sg_y5, c='k')
                axs[2,i-1].plot(x, sg_y6, c='r')

                # Δ(allele frequency) plot
                axs[3,i-1].plot(x, sg_y7, c='k')
                axs[3,i-1].plot(x, sg_y8, c='r')
                axs[3,i-1].plot(x, sg_y9, c='r')

                # sSNPs/totalSNPs plot via z-test
                # axs[1,i-1].plot(x, sg_fb_ZRatio, c='c')
                # axs[1,i-1].plot(x, sg_sb_ZRatio, c='g')
            else:
                # sSNPs/totalSNPs plot via Fisher's exact test
                axs[1,i-1].plot(x, yRatio, c='k')
                
                # G-statistic plot
                axs[2,i-1].plot(x, y5, c='k')
                axs[2,i-1].plot(x, y6, c='r')

                # Δ(allele frequency) plot
                axs[3,i-1].plot(x, y7, c='k')
                axs[3,i-1].plot(x, y8, c='r')
                axs[3,i-1].plot(x, y9, c='r')

                # sSNPs/totalSNPs plot via z-test
                # axs[1,i-1].plot(x, fb_ZRatio, c='c')
                # axs[1,i-1].plot(x, sb_ZRatio, c='g')

            # Add the 99.5 percentile line as threshold, x[-1] is the midpoint of the last sliding window of a chromosome
            axs[1,i-1].plot([plotSP, x[-1]], [thrshld, thrshld], c='r')
            # axs[1, i-1].plot(x, smThresholds_sw, c='m')

        ratioPeakL.append(max(yRatio))

        # Identify genomic regions related to the trait
        m, peaks = 0, []
        # Handle the case in which an QTL is at the very begining of the chromosome
        if swDict[i][0][6] >= thrshld:
            snpRegion.append(swDict[i][0][:2])
            if swDict[i][0][6] >= swDict[i][1][6]:
                peaks.append(swDict[i][0][1:])
            numOfSWs = 1

        while m < len(swDict[i]) - 1:
            if swDict[i][m][6] < thrshld and swDict[i][m+1][6] >= thrshld:
                snpRegion.append(swDict[i][m+1][:2])
                numOfSWs = 1
            elif swDict[i][m][6] >= thrshld:
                # A sliding window is considered as a peak if its sSNP/totalSNP is greater than or equal to the threshold and greater than those of the flanking sliding windows
                if m >= 1 and max(swDict[i][m-1][6], swDict[i][m+1][6]) <= swDict[i][m][6]:
                    peaks.append(swDict[i][m][1:])
                if swDict[i][m+1][6] > thrshld:
                    numOfSWs += 1
                elif swDict[i][m+1][6] < thrshld:
                    snpRegion[-1].extend([swDict[i][m][1], peaks, numOfSWs])
                    peaks = []
            m += 1
        # Handle the case in which an QTL is nearby the end of the chromosome
        if swDict[i][-1][6] >= thrshld:
            snpRegion[-1].extend([swDict[i][-1][1], peaks, numOfSWs])

        i += 1

    headerResults = ['CHROM','QTLStart','QTLEnd','Peaks', 'NumOfSWs']
    pd.DataFrame(snpRegion, columns=headerResults).to_csv(os.path.join(results, 'snpRegion.csv'), index=False)

    swDataFrame = pd.DataFrame(swRows, columns=['CHROM', 'sw_Str', fbID+'.AvgLD', sbID+'.AvgLD', 'sSNP', 'toatalSNP', r'sSNP/totalSNP'])
    swDataFrame['smthedRatio'] = sg_yRatio_List

    swDataFrame.to_csv(os.path.join(results, 'slidingWindows.csv'), index=False)

    misc.append(['List of the peaks of the chromosomes', ratioPeakL])

    wrnLog = os.path.join(results, 'wrnLog.csv')
    with open(wrnLog, 'w', newline='') as outF1:
        xie1 = csv.writer(outF1)
        xie1.writerow(['Type', 'Chr', 'Position', 'Warning Message'])
        xie1.writerows(wmL)

    numOfSNPOnChrFile = os.path.join(results, 'numOfSNPOnChrFile.csv')
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
ap.add_argument('-o', '--output', required=False, help='file name of the output csv file', default='BSASeq.csv')
ap.add_argument('-b', '--bulksizes', required=False, help='bulk sizes: first_bulk,second_bulk', type=lambda s: [int(t) for t in s.split(',')], default='430,385')
ap.add_argument('-p', '--popstrct', required=False, choices=['F2','RIL','BC'], help='population structure', default='F2')
ap.add_argument('-v', '--pvalues', required=False, help='cutoff p-values: real_data,simulation', type=lambda s: [float(t) for t in s.split(',')], default='0.01,0.1')
ap.add_argument('-r', '--replication', type=int, required=False, help='the number of replications for threshold calculation', default=10000)
ap.add_argument('-s', '--slidingwindow', required=False, help='size,incremental_step', type=lambda s: [int(t) for t in s.split(',')], default='2000000,10000')
ap.add_argument('-g', '--gaps', required=False, help='gaps between subplots: horizontal,vertical', type=lambda s: [float(t) for t in s.split(',')], default='0.028,0.092')
ap.add_argument('-m', '--smoothing', required=False, help='smoothing parameters: window_len,polyorder', type=lambda s: [int(t) for t in s.split(',')], default='51,3')
ap.add_argument('--smooth', type=bool, required=False, help='smooth the plot', default=False)
ap.add_argument('-e', '--region', required=False, help='interested region(s): chrm,start,end', type=lambda s: [int(t) for t in s.split(',')], default='-1')

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
minFragSize = swSize + smthWL * incrementalStep
additionalPeaks = ''

# Obtain the frequencies of the ALT allele in both bulks
fb_Freq = smAlleleFreq(popStr, fb_Size, rep)
sb_Freq = smAlleleFreq(popStr, sb_Size, rep)

path = os.getcwd()
inFile, oiFile = os.path.join(path, args.input), os.path.join(path, 'Results', 'snp_fagz.csv')
currentDT = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
results = os.path.join(path, 'Results', currentDT)
filteringPath = os.path.join(path, 'Results', 'FilteredSNPs')

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
        inputString = input('Enter the names of the desired chromosomes in order and separate each name with a comma:\n')
        chrmIDL = [x.strip() for x in inputString.split(',')]

        print('Invalid chromsome names, if any, will be removed.\n')

        # Filter out possible invalid chromosome name(s)
        invalidChrmL = []
        for ch in chrmIDL:
            if ch not in chrmList:
                invalidChrmL.append(ch)

        if invalidChrmL != []:
            print('These chrmosome names are invalid:', invalidChrmL,'\n')

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

# Handle the cases in which interested chromosomal region is entered
else:
    chrmSzD, chrmIDL, chrmSzL = {}, [], []

    l, n = len(region), len(region) % 3
    if n != 0:
        region = region[0:l-n]
        print('Each chromosome name should be followsed by the starting and ending points of the interested region.')

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

if chrmIDL == []:
    print('No valid chromosome names were entered.')
    sys.exit()

# Create a numeric ID for each chromosome, which can be used to sort the dataframe numerically by chromosome
chrmDict = {}

for i in range(1, len(chrmIDL)+1):
    chrmDict[chrmIDL[i-1]] = i

snpRawDF['ChrmSortID'] = snpRawDF['CHROM']
snpRawDF['ChrmSortID'].replace(chrmDict, inplace=True)

header = snpRawDF.columns.values.tolist()

# Obtain the bulk IDs from the header
bulks, misc, missingFlds = [], [], False

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
fb_LD, sb_LD = fbID+'.LD', sbID+'.LD'
sm_fb_AD_REF, sm_fb_AD_ALT = 'sm_'+fb_AD_REF, 'sm_'+fb_AD_ALT
sm_sb_AD_REF, sm_sb_AD_ALT = 'sm_'+sb_AD_REF, 'sm_'+sb_AD_ALT

fb_AF_CI, sb_AF_CI = fbID+'.AF_CI', sbID+'.AF_CI'
fb_AF_CI0995, sb_AF_CI0995 = fb_AF_CI+'0995', sb_AF_CI+'0995'
fb_AF_CI0005, sb_AF_CI0005 = fb_AF_CI+'0005', sb_AF_CI+'0005'

if os.path.isfile(os.path.join(path, 'Results', 'COMPLETE.txt')) == False:
    bsaSNPs = snpFiltering(snpRawDF)

    # Remove low LD SNPs to meet the z-test sample size requirment
    # bsaSNPs = bsaSNPs[(bsaSNPs[fb_LD]>=30) & (bsaSNPs[sb_LD]>=30)]

    # Calculate z scores
    bsaSNPs['fb_ZScore'] = (bsaSNPs[fb_AD_ALT] / bsaSNPs[fb_LD] - fb_Freq) * np.sqrt(bsaSNPs[fb_LD] / (fb_Freq * (1 - fb_Freq)))
    bsaSNPs['sb_ZScore'] = (bsaSNPs[sb_AD_ALT] / bsaSNPs[sb_LD] - sb_Freq) * np.sqrt(bsaSNPs[sb_LD] / (sb_Freq * (1 - fb_Freq)))

    # Calculate simulated ALT reads for each SNP under null hypothesis
    bsaSNPs[sm_fb_AD_ALT] = np.random.binomial(bsaSNPs[fb_LD], fb_Freq)
    bsaSNPs[sm_fb_AD_REF] = bsaSNPs[fb_LD] - bsaSNPs[sm_fb_AD_ALT]
    bsaSNPs[sm_sb_AD_ALT] = np.random.binomial(bsaSNPs[sb_LD], sb_Freq)
    bsaSNPs[sm_sb_AD_REF] = bsaSNPs[sb_LD] - bsaSNPs[sm_sb_AD_ALT]

    # Calculate allele frequency
    bsaSNPs[fb_AF] = bsaSNPs[fb_AD_ALT]/bsaSNPs[fb_LD]
    bsaSNPs[sb_AF] = bsaSNPs[sb_AD_ALT]/bsaSNPs[sb_LD]
    bsaSNPs['Delta.AF'] = bsaSNPs[sb_AF] - bsaSNPs[fb_AF]

    # Calculate G-statistic
    bsaSNPs['G_S'] = gStatistic_Array(bsaSNPs[fb_AD_REF], bsaSNPs[fb_AD_ALT], bsaSNPs[sb_AD_REF], bsaSNPs[sb_AD_ALT])

    try:
        from fisher import pvalue_npy
        # Create new columns for Fisher's exact test P-values and simulated P-values
        print('Perform Fisher\'s exact test.')
        fb_AD_ALT_Arr = bsaSNPs[fb_AD_ALT].to_numpy(dtype=np.uint)
        fb_AD_REF_Arr = bsaSNPs[fb_AD_REF].to_numpy(dtype=np.uint)
        sb_AD_ALT_Arr = bsaSNPs[sb_AD_ALT].to_numpy(dtype=np.uint)
        sb_AD_REF_Arr = bsaSNPs[sb_AD_REF].to_numpy(dtype=np.uint)

        __, __, bsaSNPs['FE_P'] = pvalue_npy(fb_AD_ALT_Arr, fb_AD_REF_Arr, sb_AD_ALT_Arr, sb_AD_REF_Arr)
        # bsaSNPs['FE_OR'] = (fb_AD_ALT_Arr * sb_AD_REF_Arr) / (fb_AD_REF_Arr * sb_AD_ALT_Arr)

        sm_fb_AD_ALT_Arr = bsaSNPs[sm_fb_AD_ALT].to_numpy(dtype=np.uint)
        sm_fb_AD_REF_Arr = bsaSNPs[sm_fb_AD_REF].to_numpy(dtype=np.uint)
        sm_sb_AD_ALT_Arr = bsaSNPs[sm_sb_AD_ALT].to_numpy(dtype=np.uint)
        sm_sb_AD_REF_Arr = bsaSNPs[sm_sb_AD_REF].to_numpy(dtype=np.uint)

        __, __, bsaSNPs['sm_FE_P'] = pvalue_npy(sm_fb_AD_ALT_Arr, sm_fb_AD_REF_Arr, sm_sb_AD_ALT_Arr, sm_sb_AD_REF_Arr)
        # bsaSNPs['sm_FE_OR'] = (sm_fb_AD_ALT_Arr * sm_sb_AD_REF_Arr) / (sm_fb_AD_REF_Arr * sm_sb_AD_ALT_Arr)

        print(f'Fisher\'s exact test completed, time elapsed: {(time.time()-t0)/60} minutes.')

        print('Calculate thresholds of \u0394(allele frequency) and G-statistic.')
        bsaSNPs['STAT'] = bsaSNPs.apply(statistics, axis=1)

        # Create new columns for Fisher's exact test results, allele frequency, Δ(allele frequency) confidence intervals, and G-statistic thresholds
        bsaSNPs[[fb_AF_CI, sb_AF_CI, 'DAF_CI', 'GS_CI']] = pd.DataFrame(bsaSNPs.STAT.values.tolist(), index=bsaSNPs.index)

        print(f'Calculating thresholds of \u0394(allele frequency) and G-statistic completed, time elapsed: {(time.time()-t0)/60} minutes.')

    except ImportError:
        from scipy.stats import fisher_exact
        print('Perform Fisher\'s exact test. This step can take a few hours; the more SNPs in the dataset or the higher the sequencing depth, the longer will it take.')
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

    # Reorgnaize the columns
    reorderColumns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', fb_GT, fb_AD, fb_AD_REF, fb_AD_ALT, fb_LD, 'fb_ZScore', sm_fb_AD_ALT, fb_AF, fb_AF_CI0005, fb_AF_CI0995, fb_AF_CI, fb_GQ, sb_GT, sb_AD, sb_AD_REF, sb_AD_ALT, sb_LD, 'sb_ZScore', sm_sb_AD_ALT, sb_AF, sb_AF_CI0005, sb_AF_CI0995, sb_AF_CI, sb_GQ, 'Delta.AF', 'FE_P', 'sm_FE_P', 'G_S', 'DAF_CI0005', 'DAF_CI0995', 'GS_CI0995', 'DAF_CI', 'GS_CI', 'STAT']

    # Remove unnecessary columns and reorgnaize the columns
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

# Calculate or retrieve the threshold. The threshoslds are normally in the range from 0.12 to 0.12666668
if os.path.isfile(os.path.join(path, 'Results', 'threshold.txt')) == False:
    if 'fisher' in sys.modules:
        thrshld = smThresholds_gw(bsaSNPs)[1]
    else:
        thrshld = smThresholds_proximal(bsaSNPs)[1]

    with open(os.path.join(path, 'Results', 'threshold.txt'), 'w') as xie:
        xie.write(str(thrshld))
else:
    with open(os.path.join(path, 'threshold.txt'), 'r') as du:
        thrshld = float(du.readline().strip())

# Identify likely trait-associated SNPs
bsaSNPs['sSNP'] = np.where(bsaSNPs['FE_P'] < alpha, 1, 0)

# Plot layout setup
heightRatio = [1,0.8,0.8,0.8]
fig, axs = plt.subplots(nrows=len(heightRatio), ncols=len(chrmIDL), figsize=(20, 18), sharex='col', sharey='row', 
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
fig.subplots_adjust(top=0.98, bottom=0.043, left=0.064, right=0.995, hspace=hGap, wspace=wGap)
fig.suptitle('Genomic position (\u00D710 Mb)', y=0.002, ha='center', va='bottom')
fig.text(0.001, 0.995, 'A', weight='bold', ha='left', va='top')
fig.text(0.001, 0.688, 'B', weight='bold', ha='left', va='bottom')
fig.text(0.001, 0.465, 'C', weight='bold', ha='left', va='bottom')
fig.text(0.001, 0.244, 'D', weight='bold', ha='left', va='bottom')

fig.savefig(os.path.join(results, 'PyBSASeq.pdf'))
# fig.savefig(os.path.join(results, 'PyBSASeq.png'), dpi=600)

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