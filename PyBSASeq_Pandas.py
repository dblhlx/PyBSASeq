"""
Created on Fri Oct  5 08:22:16 2018
@author: Jianbo Zhang
"""
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import os
import datetime
import argparse
import time
import csv
import sys


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
    # Many reference genomes contain unmapped fragments that tend to be small and are not informative to SNP-trait association, filtering them out makes the chromosome list more readable 
    # Additionally, users may enter wrong chromosome names that could lead to unexpected behaviors 
    chrmSizeL, smallOrWrongChrmL = [], []
    for chrmID in chromosomeList:
        # Handle the case in which a wrong chromosome name was entered by the user
        if chrmID in chrmRawList:
            chrmSize = df[df.CHROM==chrmID]['POS'].max()
        else:
            smallOrWrongChrmL.append(chrmID)
            continue

        if chrmSize <= swSize:
            smallOrWrongChrmL.append(chrmID)
        else:
            chrmSizeL.append(chrmSize)

    for chrmID in smallOrWrongChrmL:
        chromosomeList.remove(chrmID)

    return [chrmSizeL, chromosomeList]


def snpFiltering(df):
    print('Perform SNP filtering')
    global snpDF, misc

    df = df.copy()

    # Identify unmapped SNPs
    df_Unmapped = df[~df.CHROM.isin(chrmIDL)]
    df_Unmapped.to_csv(os.path.join(filteringPath, 'unmapped.csv'), index=None)

    # Remove unmapped SNPs
    df = df[df.CHROM.isin(chrmIDL)]

    # Identify SNPs with 'NA' value(s)
    df_NA = df[df.isnull().any(axis=1)]
    df_NA.to_csv(os.path.join(filteringPath, 'na.csv'), index=None)

    # Remove SNPs with 'NA' value(s)
    df.dropna(inplace=True)
    misc.append(['Number of SNPs after NA drop', len(df.index)])

    # Identify SNPs with a single ALT allele
    df_1ALT = df[~df['ALT'].str.contains(',')]

    # Identify one-ALT SNPs with zero REF read in both bulks
    df_1ALT_Fake = df_1ALT[(df_1ALT[fb_AD].str.split(',',expand=True)[0]=='0') & \
        (df_1ALT[sb_AD].str.split(',',expand=True)[0]=='0')]
    df_1ALT_Fake.to_csv(os.path.join(filteringPath, '1altFake.csv'), index=None)

    # Remove one-ALT SNPs with zero REF read in both bulks
    df_1ALT_Real = df_1ALT[~((df_1ALT[fb_AD].str.split(',',expand=True)[0]=='0') & \
        (df_1ALT[sb_AD].str.split(',',expand=True)[0]=='0'))]
    df_1ALT_Real.to_csv(os.path.join(filteringPath, '1altReal.csv'), index=None)

    # Identify SNPs with more than one ALT allele
    df_mALT = df[df['ALT'].str.contains(',')]

    # Identify SNPs with two ALT alleles. Split ALT with ',' once, and the second element should not contain a ','
    df_2ALT = df_mALT[~df_mALT['ALT'].str.split(',', n=1, expand=True)[1].str.contains(',')]

    # A two-ALT SNP is a real SNP if the REF read is zero in both bulks
    # Making a copy to suppress the warning message. Updating the AD values of these SNPs is required
    df_2ALT_Real = df_2ALT[(df_2ALT[fb_AD].str.split(',',expand=True)[0]=='0') & \
        (df_2ALT[sb_AD].str.split(',',expand=True)[0]=='0')].copy()

    # Update the AD values of the above SNPs by removing the REF read which is zero
    df_2ALT_Real[fb_AD] = df_2ALT_Real[fb_AD].str.slice(start=2)
    df_2ALT_Real[sb_AD] = df_2ALT_Real[sb_AD].str.slice(start=2)
    df_2ALT_Real.to_csv(os.path.join(filteringPath, '2altReal.csv'), index=None)

    # The two-ALT SNP may be cuased by allele heterozygosity if the REF read in not zero
    # Repetitive sequences in the genome or sequencing artifacts are other possibilities
    df_2ALT_Het = df_2ALT[~((df_2ALT[fb_AD].str.split(',',expand=True)[0]=='0') & \
        (df_2ALT[sb_AD].str.split(',',expand=True)[0]=='0'))]

    # Identify SNPs with three or more ALT alleles
    df_3omALT = df_mALT[df_mALT['ALT'].str.split(',', n=1, expand=True)[1].str.contains(',')]

    snpDF_Het = pd.concat([df_3omALT, df_2ALT_Het])
    snpDF_Het.to_csv(os.path.join(filteringPath, 'heterozygousLoci.csv'), index=None)

    # Concatenate 1ALT_Real and 2ALT_Real
    snpDF = pd.concat([df_1ALT_Real, df_2ALT_Real])

    # In case the input file contains Indels
    try:
        # Identify inDels, REF/Alt allele with more than 1 base
        df_InDel = snpDF[(snpDF['REF'].str.len()>1) | (snpDF['ALT'].str.split(',', expand=True)[0].str.len()>1) | \
        (snpDF['ALT'].str.split(',', expand=True)[1].str.len()>1)]

        df_InDel.to_csv(os.path.join(filteringPath, 'InDel.csv'), index=None)

        # Create the SNP dataframe
        snpDF = snpDF[~((snpDF['REF'].str.len()>1) | (snpDF['ALT'].str.split(',', expand=True)[0].str.len()>1) | \
        (snpDF['ALT'].str.split(',', expand=True)[1].str.len()>1))]

    except KeyError:
        print('All the ALT loci contain only one allele')
        # Identify inDels, REF/Alt allele with more than 1 base
        df_InDel = snpDF[(snpDF['REF'].str.len()>1) | (snpDF['ALT'].str.len()>1)]

        df_InDel.to_csv(os.path.join(filteringPath, 'InDel.csv'), index=None)

        # Create the SNP dataframe
        snpDF = snpDF[~((snpDF['REF'].str.len()>1) | (snpDF['ALT'].str.len()>1))]


    snpDF.sort_values(['ChrmSortID', 'POS'], inplace=True)

    print(f'SNP filtering completed, time elapsed: {(time.time()-t0)/3600} hours')


def statistics(row):
    # Perform Fisher's exact test for each SNP using the actual REF/ALT reads
    try:
        fe = fisher_exact([[row[fb_AD_REF], row[fb_AD_ALT]], [row[sb_AD_REF], row[sb_AD_ALT]]])
    except TypeError:
        fe = 'NA'

    # Perform Fisher's exact test for each SNP using the simulated REF/ALT reads
    try:
        sm_FE = fisher_exact([[row[sm_fb_AD_REF], row[sm_fb_AD_ALT]], [row[sm_sb_AD_REF], row[sm_sb_AD_ALT]]])
    except TypeError:
        sm_FE = 'NA'

    return [fe, sm_FE]


def smThresholds(DF):
    print('Calculate the threshold of ltaSNPs/totalSNPs.')
    ratioLi = []
    for __ in range(rep):
        sm_SNP_SMPL = DF.sample(snpPerSW, replace=True)
        sm_ltaSNP_SMPL = sm_SNP_SMPL[sm_SNP_SMPL['sm_FE_P']<smAlpha]

        ratioLi.append(len(sm_ltaSNP_SMPL.index)/snpPerSW)

    misc.append(['Genome-wide ltaSNP/totalSNP ratio threshold', np.percentile(ratioLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])])
    print(f'Threshold calculation completed, time elapsed: {(time.time()-t0)/3600} hours')

    return np.percentile(ratioLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])


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


def bsaseqPlot(chrmIDL, datafr, datafrT):
    '''
    wmL: list of warning messages
    points: a dictionary with the chromosome ID as its keys; the value of each key is a list containing
            the chromosome ID, the ltaSNP/totalSNP ratio in each sliding window, and the midpoint of
            the sliding window
    '''
    print('Prepare SNP data for plotting via the sliding window algorithm')
    global misc
    wmL, points = [], {}

    # Analyze each chromsome separately
    numOfSNPOnChr, ratioPeakL, snpRegion = [], [], []
    i = 1
    for chrmID in chrmIDL:
        ch = datafr[datafr.CHROM==chrmID]
        chT = datafrT[datafrT.CHROM==chrmID]
        numOfSNPOnChr.append([chrmID, len(ch.index), len(chT.index), len(ch.index)/len(chT.index)])

        # Sliding window. swStr: the begining of the window; swEnd: the end of the window; icrs: incremental step
        swStr, swEnd, icrs = 1, swSize, 10000
        plotSP = swEnd/2
        # x and y are lists, each sliding window represents a single data point
        x, y, yT, yRatio = [], [], [], []
        while swEnd <= chT['POS'].max():
            # A single sliding window - a dataframe
            # swDF: ltaSNPs in a sliding window; swDFT: all SNPs in a sliding window
            swDF = ch[(ch.POS>=swStr) & (ch.POS<=swEnd)]
            swDFT = chT[(chT.POS>=swStr) & (chT.POS<=swEnd)]

            rowInSwDF = len(swDF.index)                 # number of ltaSNPs in a sliding window
            rowInSwDFT = len(swDFT.index)               # number of total SNPs in a sliding window

            x.append((swStr+swEnd)/2)                   # Append the midpoint of a sliding window to x
            y.append(rowInSwDF)                         # Append number of ltaSNPs in a sliding window to y
            yT.append(rowInSwDFT)                       # Append number of totalSNPs in a sliding window to yT

            # len(swDL) or len(swDLT) would be zero if no SNP in a sliding window
            try:
                yRatio.append(float(rowInSwDF/rowInSwDFT))  # Append the ratio of ltaSNP/totalSNP in a sliding window to yRatio
            except ZeroDivisionError as e:
                wmL.append(['No SNP', i, int((swStr+swEnd)/2), e])
                zeroSNP(yRatio)

            if i not in points:
                points[i] = [[chrmID, yRatio[-1], int((swStr+swEnd)/2), rowInSwDFT]]
            else:
                points[i].append([chrmID, yRatio[-1], int((swStr+swEnd)/2), rowInSwDFT])

            swStr += icrs
            swEnd += icrs

        # Replace the 'empty' values at the begining of the lists with nearest non-empty value
        if 'empty' in yRatio:
            replaceZero(yRatio)

        pIndex = 0
        while points[i][pIndex][1] == 'empty':
            pIndex += 1

        j = 0
        while j < pIndex:
            points[i][j][1] = points[i][pIndex][1]
            j += 1

        # Handle the plot with a single column (chromosome)
        if len(chrmIDL) == 1:
            # SNP plot
            axs[0].plot(x, y, c='k')
            axs[0].plot(x, yT, c='b')
            axs[0].set_title('Chr'+chrmID)

            # ltaSNP/totalSNP plot
            axs[1].plot(x, yRatio, c='k')
            # sg_yRatio = savgol_filter(yRatio, 51, 3)
            # axs[1].plot(x, sg_yRatio, c='r')

            # Add the 99.5 percentile line as threshold, x[-1] is the midpoint of the last sliding window of a chromosome
            axs[1].plot([plotSP, x[-1]], [thrshld, thrshld], c='r')

            # Set up x-ticks
            axs[1].set_xticks(np.arange(0, max(x), 10000000))
            ticks = axs[1].get_xticks()*1e-7
            axs[1].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0].set_ylabel('Number of SNPs')
                axs[1].set_ylabel(r'ltaSNP/totalSNP')

        # Handle the plot with multiple columns (chromosomes)
        else:
            # SNP plot
            axs[0,i-1].plot(x, y, c='k')
            axs[0,i-1].plot(x, yT, c='b')
            axs[0,i-1].set_title('Chr'+chrmID)

            # ltaSNP/totalSNP plot
            axs[1,i-1].plot(x, yRatio, c='k')
            # sg_yRatio = savgol_filter(yRatio, 51, 3)
            # axs[1,i-1].plot(x, sg_yRatio, c='r')

            # Add the 99.5 percentile line as threshold, x[-1] is the midpoint of the last sliding window of a chromosome
            axs[1,i-1].plot([plotSP, x[-1]], [thrshld, thrshld], c='r')

            # Set up x-ticks
            axs[1,i-1].set_xticks(np.arange(0, max(x), 10000000))
            ticks = axs[1,i-1].get_xticks()*1e-7
            axs[1,i-1].set_xticklabels(ticks.astype(int))

            # Add ylabels to the first column of the subplots
            if i==1:
                axs[0,i-1].set_ylabel('Number of SNPs')
                axs[1,i-1].set_ylabel(r'ltaSNP/totalSNP')

        ratioPeakL.append(max(yRatio))

        # Identify genomic regions related to the trait
        m, peaks = 0, []
        # Handle the case in which an QTL is at the very begining of the chromosome
        if points[i][0][1] >= thrshld:
            snpRegion.append(points[i][0][:3])
            if points[i][0][1] > points[i][1][1]:
                peaks.append(points[i][0][1:])
            numOfSWs = 1

        while m < len(points[i]) - 1:
            if points[i][m][1] < thrshld and points[i][m+1][1] >= thrshld:
                snpRegion.append(points[i][m+1][:3])
                numOfSWs = 1
            elif points[i][m][1] >= thrshld:
                # A sliding window is considered as a peak if its ltaSNP/totalSNP is greater than or equal to the threshold and greater than those of the flanking sliding windows
                if m >= 1 and max(points[i][m-1][1], points[i][m+1][1]) < points[i][m][1]:
                    peaks.append(points[i][m][1:])
                if points[i][m+1][1] > thrshld:
                    numOfSWs += 1
                elif points[i][m+1][1] < thrshld:
                    snpRegion[-1].extend([points[i][m][2], peaks, numOfSWs])
                    peaks = []
            m += 1
        # Handle the case in which an QTL is nearby the end of the chromosome
        if points[i][-1][1] >= thrshld:
            snpRegion[-1].extend([points[i][-1][2], peaks, numOfSWs])

        i += 1

    misc.append(['List of the peaks of the chromosomes', ratioPeakL])

    headerResults = ['CHROM',r'ltaSNP/totalSNP','QTLStart','QTLEnd','Peaks', 'NumOfSWs']
    pd.DataFrame(snpRegion, columns=headerResults).to_csv(os.path.join(results, args['output']), index=False)

    wrnLog = os.path.join(results, 'wrnLog.csv')
    with open(wrnLog, 'w', newline='') as outF1:
        xie1 = csv.writer(outF1)
        xie1.writerow(['Type', 'Chr', 'Position', 'Warning Message'])
        xie1.writerows(wmL)

    numOfSNPOnChrFile = os.path.join(results, 'numOfSNPOnChrFile.csv')
    with open(numOfSNPOnChrFile, 'w', newline='') as outF2:
        xie2 = csv.writer(outF2)
        xie2.writerow(['Chromosome', 'Num of ltaSNPs', 'Num of totalSNPs', r'ltaSNP/totalSNP'])
        xie2.writerows(numOfSNPOnChr)

    print(f'Plotting completed, time elapsed: {(time.time()-t0)/3600} hours')


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
ap.add_argument('-i', '--input', required=False, help='file name of the GATK4-generated tsv file', default='smallTestFile.tsv')
ap.add_argument('-o', '--output', required=False, help='file name of the output csv file', default='BSASeq.csv')
ap.add_argument('-f', '--fbsize', type=int, required=False, help='number of individuals in the first bulk', default=430)
ap.add_argument('-s', '--sbsize', type=int, required=False, help='number of individuals in the second bulk', default=385)
ap.add_argument('-p', '--popstrct', required=False, choices=['F2','RIL','BC'], help='population structure', default='F2')
ap.add_argument('--alpha', type=float, required=False, help='p-value for fisher\'s exact test', default=0.01)
ap.add_argument('--smalpha', type=float, required=False, help='p-value for calculating threshold', default=0.1)
ap.add_argument('-r', '--replication', type=int, required=False, help='the number of replications for threshold calculation', default=10000)
ap.add_argument('--swsize', type=int, required=False, help='sliding windows size', default=2000000)
ap.add_argument('--hgap', type=float, required=False, help='distance between rows of subplots', default=0.028)
ap.add_argument('--wgap', type=float, required=False, help='distance between columns of subplots', default=0.092)
args = vars(ap.parse_args())

popStr = args['popstrct']
rep = args['replication']
fb_Size, sb_Size = args['fbsize'], args['sbsize']
alpha, smAlpha = args['alpha'], args['smalpha']
swSize = args['swsize']
hGap, wGap = args['hgap'], args['wgap']

fb_Freq = smAlleleFreq(popStr, fb_Size, rep)
sb_Freq = smAlleleFreq(popStr, sb_Size, rep)

path = os.getcwd()
inFile, oiFile = os.path.join(path, args['input']), os.path.join(path, 'snp100SE_fe.csv')
currentDT = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
results = os.path.join(path, 'Results', currentDT)
filteringPath = os.path.join(path, 'FilteredSNPs')

if not os.path.exists(results):
    os.makedirs(results)

if not os.path.exists(filteringPath):
    os.makedirs(filteringPath)

# Generte a SNP dataframe from the GATK4-generated tsv file
snpRawDF = pd.read_csv(inFile, delimiter='\t', encoding='utf-8', dtype={'CHROM':str})

# Create a chromosome list, which can be very long because of the unmapped fragments
chrmRawList = list(set(snpRawDF['CHROM'].tolist()))

# Filter out chromosomes and unmapped fragments smaller than the sliding window
# Make the chromosome list more readable and meaningful
chrmList = chrmFiltering(snpRawDF, chrmRawList)[1]

print(chrmList)
print('\n')
print(f'Above is the chromosome list, from which you can select the chromosome name(s) for the next step. Chromosomes or unmapped fragments smaller than the sliding window ({swSize} bp) are filtered out. Adjust the sliding window size with option \'--swsize\' to include desired small chromosomes.\n')

# Print the chromosome list on the screen to let the user to select desired chromosome(s)
rightInput = ''
while rightInput.lower() != 'yes':
    inputString = input('Enter chromosome names in order and separate each name with a comma:\n')
    chrmIDL = [x.strip() for x in inputString.split(',')]
    print('Sorted chromosome list:')
    print(chrmIDL)
    rightInput = input('Are the chromosome names in the above list in the right order (yes or no)?\n')

# Filter out possible wrong chromosome name(s) and chromosomes smaller than the sliding window
# Create a list containing the sizes of all the chromosomes
chrmCheck = chrmFiltering(snpRawDF, chrmIDL)
chrmSzL, chrmIDL = chrmCheck[0], chrmCheck[1]

if chrmIDL == []:
    print('No valid chromosome name was entered')
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
    fb_AD, sb_AD = fbID+'.AD', sbID+'.AD'
    fb_GQ, sb_GQ = fbID+'.GQ', sbID+'.GQ'
except (NameError, IndexError):
    print('The allele depth (AD) field is missing. Please include the AD field in the input file.')
    sys.exit()

# Check if any required field is missing in the input file
requiredFields = ['CHROM', 'POS', 'REF', 'ALT', fb_AD, fb_GQ, sb_AD, sb_GQ]
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

fb_SI, sb_SI = fbID+'.SI', sbID+'.SI'
fb_AD_REF, fb_AD_ALT = fb_AD + '_REF', fb_AD + '_ALT'
sb_AD_REF, sb_AD_ALT = sb_AD + '_REF', sb_AD + '_ALT'
fb_LD, sb_LD = fbID+'.LD', sbID+'.LD'
sm_fb_AD_REF, sm_fb_AD_ALT = 'sm_'+fb_AD_REF, 'sm_'+fb_AD_ALT
sm_sb_AD_REF, sm_sb_AD_ALT = 'sm_'+sb_AD_REF, 'sm_'+sb_AD_ALT

if os.path.isfile(os.path.join(path, 'COMPLETE.txt')) == False:
    snpFiltering(snpRawDF)

    # Obtain REF reads, ALT reads, and locus reads of each SNP
    snpDF[[fb_AD_REF, fb_AD_ALT]] = snpDF[fb_AD].str.split(',', expand=True).astype(int)
    snpDF[fb_LD] = snpDF[fb_AD_REF] + snpDF[fb_AD_ALT]

    snpDF[[sb_AD_REF, sb_AD_ALT]] = snpDF[sb_AD].str.split(',', expand=True).astype(int)
    snpDF[sb_LD] = snpDF[sb_AD_REF] + snpDF[sb_AD_ALT]

    # Filter out the SNPs with zero locus reads in either bulk
    snpDF_0LD = snpDF[~((snpDF[fb_LD]>0) & (snpDF[sb_LD]>0))]
    snpDF = snpDF[(snpDF[fb_LD]>0) & (snpDF[sb_LD]>0)]
    snpDF_0LD.to_csv(os.path.join(filteringPath, '0ld.csv'), index=None)

    # Calculate simulated ALT reads for each SNP under null hypothesis
    snpDF[sm_fb_AD_ALT] = np.random.binomial(snpDF[fb_LD], fb_Freq)
    snpDF[sm_fb_AD_REF] = snpDF[fb_LD] - snpDF[sm_fb_AD_ALT]
    snpDF[sm_sb_AD_ALT] = np.random.binomial(snpDF[sb_LD], sb_Freq)
    snpDF[sm_sb_AD_REF] = snpDF[sb_LD] - snpDF[sm_sb_AD_ALT]

    print('Perform Fisher\'s exact test. This step can take a few hours; the more SNPs in the dataset or the higher the sequencing depth, the longer will it take.')
    snpDF['STAT'] = snpDF.apply(statistics, axis=1)
    print(f'Fisher\'s exact test completed, time elapsed: {(time.time()-t0)/3600} hours')

    # Create new columns for Fisher's exact test results
    snpDF[['fisher_exact', 'sm_FE']] = pd.DataFrame(snpDF.STAT.values.tolist(), index=snpDF.index)

    # Create new columns for Fisher's exact test P-values or simulated P-values
    snpDF['FE_P'] = snpDF['fisher_exact'].apply(lambda x: x[1]).astype(float)
    snpDF['sm_FE_P'] = snpDF['sm_FE'].apply(lambda x: x[1]).astype(float)

    # Remove unnecessary columns and reorgnaize the columns
    reorderColumns = ['CHROM', 'POS', 'REF', 'ALT', fb_AD_REF, fb_AD_ALT, fb_LD, sm_fb_AD_ALT, fb_GQ, sb_AD_REF, sb_AD_ALT, sb_LD, sm_sb_AD_ALT,  sb_GQ, 'FE_P', 'sm_FE_P']

    snpDF = snpDF[reorderColumns]

    snpDF.to_csv(oiFile, index=None)

    with open(os.path.join(path, 'COMPLETE.txt'), 'w') as xie:
        xie.write('Statistical calculation is completed!')
else:
    snpDF = pd.read_csv(oiFile, dtype={'CHROM':str})

# The above calculation may generate 'NA' value(s) for some SNPs. Remove SNPs with such 'NA' value(s)
snpDF.dropna(inplace=True)
misc.append(['Number of SNPs after drop of SNPs with calculation-generated NA value', len(snpDF.index)])

# Filter out SNPs with a low genotype quality score
snpDF = snpDF[(snpDF[fb_GQ]>=20) & (snpDF[sb_GQ]>=20)]

misc.append(['Dataframe filtered with genotype quality scores', len(snpDF.index)])

# Calculate the average number of SNPs in a sliding window
snpPerSW = int(len(snpDF.index) * swSize/sum(chrmSzL))

misc.append(['Average SNPs per sliding window', snpPerSW])
misc.append([f'Average locus depth in bulk {fbID}', snpDF[fb_LD].mean()])
misc.append([f'Average locus depth in bulk {sbID}', snpDF[sb_LD].mean()])

# Calculate or retrieve the threshold. The threshoslds are normally in the range from 0.12 to 0.12666668
if os.path.isfile(os.path.join(path, 'threshold.txt')) == False:
    thrshld = smThresholds(snpDF)[1]
    with open(os.path.join(path, 'threshold.txt'), 'w') as xie:
        xie.write(str(thrshld))
else:
    with open(os.path.join(path, 'threshold.txt'), 'r') as du:
        thrshld = float(du.readline().strip())

# Identify likely trait-associated SNPs
fe = snpDF[snpDF['FE_P']<alpha]

# Plot layout setup
heightRatio = [1,0.8]
fig, axs = plt.subplots(nrows=len(heightRatio), ncols=len(chrmIDL), figsize=(20, 10), sharex='col', sharey='row', 
        gridspec_kw={'width_ratios': chrmSzL, 'height_ratios': heightRatio})

# Perform plotting
bsaseqPlot(chrmIDL, fe, snpDF)

# Handle the plot with a single column (chromosome)
if len(chrmIDL) == 1:
    fig.align_ylabels(axs[:])
# Handle the plot with multiple columns (chromosomes)
else:
    fig.align_ylabels(axs[:, 0])

# fig.tight_layout(pad=0.15, rect=[0, 0.035, 1, 1])
fig.subplots_adjust(top=0.96, bottom=0.073, left=0.064, right=0.995, hspace=hGap, wspace=wGap)
fig.suptitle('Genomic position (\u00D710 Mb)', y=0.002, ha='center', va='bottom')
fig.text(0.001, 0.995, 'A', weight='bold', ha='left', va='top')
fig.text(0.001, 0.435, 'B', weight='bold', ha='left', va='bottom')

fig.savefig(os.path.join(results, 'PyBSASeq.pdf'))
# fig.savefig(os.path.join(results, 'PyBSASeq.png'), dpi=600)

misc.append(['Running time', [(time.time()-t0)/3600]])

with open(os.path.join(results, 'misc_info.csv'), 'w', newline='') as outF:
    xie = csv.writer(outF)
    xie.writerows(misc)