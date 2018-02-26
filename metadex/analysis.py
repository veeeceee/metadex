import glob
import math
import numpy as np
import pandas as pd
import sys
import os
import errno
import itertools
from skbio.stats import subsample_counts
from functools import reduce
from Bio import Entrez
import time
import seaborn as sns
import re
#import requests
from scipy import stats
############################################################################################
'''
DETERMING THE RELATIONSHIP BETWEEN CHANGE IN DIVERSITY AND CHANGE IN FUNCTION











'''
############################################################################################

pd.options.mode.chained_assignment = None
###SEARCH BY GENE###
#takes dataframe
def search_by_gene(query, dfCounts):
    """
    Summary: gets the subset of a dataframe by the query of interest
    Args:
        query (str): search term of interest, use [space] for OR & | [pipe] for AND
        dfCounts (pandas.DataFrame): counts dataframe
    Returns:
        geneSorted (pandas.DataFrame): dataframe for genes containing query
    """

    queryList = query.replace(" ", "|").replace("+", " ")
    geneSorted = dfCounts[dfCounts['gene function'].str.contains(queryList)]
    return geneSorted

def merge_dfs2(ldf, rdf):
    """
    Summary: merges two counts (with hierarchy)
    Args:
        ldf (pandas.DataFrame): left dataframe
        rdf (pandas.DataFrame): right dataframe
    Returns:
        pd.merge(ldf, rdf) (pandas.DataFrame): merged dataframe
    """

    print('yep, again...')
    return ldf.merge(rdf, how='outer', on=['gene function', 'organism', 'phylum', 'class', 'order', 'family'])

def merge_all_counts(studyName):
    """
    Summary: merges all counts w/ hierarchy
    Args:
        studyName (str): study name/directory name (ideally one and the same)
    Returns:
        mergedCounts (pandas.DataFrame): merged dataframe
    """

    annotatedCountsList = glob.glob(studyName + '/norm_taxonomy/*norm_taxonomy.csv')
    annotatedDFsList = [pd.read_csv(count) for count in annotatedCountsList]
    [annDF.pop('Unnamed: 0') for annDF in annotatedDFsList]
    mergedCounts  = reduce(merge_dfs2, annotatedDFsList)
    mergedCounts.fillna(0)
    mergedCounts.to_csv(str(studyName) + '_allcounts.csv')
    return mergedCounts

def split_hierarchically(geneSorted):
    """
    Summary: breaks out the Entrez lineage column into separate columns for each taxonomic level
    Args:
        geneSorted (pandas.DataFrame): annotated counts with column for hierarchy
    Returns:
        taxList (list): list of counts dataframes, each at a different taxonomic level
    """

    taxList = []
    taxList.append(geneSorted.groupby('phylum').sum())
    taxList.append(geneSorted.groupby('class').sum())
    taxList.append(geneSorted.groupby('order').sum())
    taxList.append(geneSorted.groupby('family').sum())
    return taxList

def subset_by_gene(studyName, geneName, testMerge):
    """
    Summary: subsets counts dataframe by query (gene name/gene family name)
    Args:
        studyName (str): name of study directory
        geneName (str): query of interest
        testMerge (pandas.DataFrame): merged annotated counts
    Returns:
        geneDF (pandas.DataFrame): dataframe of counts only for gene/query of interest
    """

    os.chdir(studyName)
    geneDF = search_by_gene(geneName, testMerge)
    geneDF.to_csv(str(geneName)+'_'+'allcounts.csv')
    os.chdir('..')
    return geneDF

##############################################################################################
####PLOT RENYI PROFILE###
#input is a data frame of values with rows representing taxa (in the individual columns) and columns representing the individual samples, as labelled with [group]_[samplenumber]

###MAKE COUNTS PROPORTIONAL (FOR RENYI PROFILE ETC)###
def make_counts_proportional(counts):
    """
    Summary: converts numerical counts to fractional representation
    Args:
        counts (pandas.Series): counts
    Returns:
        counts (pandas.Series): counts represented proportionally
    """

    counts = counts/counts.sum()
    return counts
###CALCULATE RENYI ENTROPY###
def renyientropy(px, alpha):
    """
    Summary: calculates renyi entropy of array px at alpha
    Args:
        px (array-like object): distribution  of counts
        alpha (): alpha value at which renyi is to be calculated 
    Returns:
        renyi (): renyi entropy of px at alpha 
    """

    if alpha < 0:
        raise ValueError("alpha must be a non-negative real number")
    elif alpha == 0:
        renyi = np.log2(len(px[px > 0])+.00000000000001)
    elif alpha == 1:
        renyi = -np.sum(np.nan_to_num(px*np.log2(px+.00000000000001)))
    elif alpha == np.inf:
        renyi = -np.log2(px.max())
    else:
        renyi = 1/(1-alpha) * np.log2((px**alpha).sum())
    return renyi
#####CALCULATE RENYI PROFILE#####
def calculate_renyi_profile(lvlCounts):
    """
    Summary: calculates renyi profile of a counts distribution
    Args:
        lvlCounts (array-like object): counts distribution
    Returns:
        profile_values (list): renyi profile of lvlCounts
    """

    profile_values = []
    if lvlCounts.sum() == 0:
        profile_values = [0,0,0,0,0,0,0,0,0,0,0]
    else:
        lvlCounts = lvlCounts/lvlCounts.sum()
        counts_distribution = lvlCounts

        profile_values.append(renyientropy(counts_distribution, 0))
        profile_values.append(renyientropy(counts_distribution, .25))
        profile_values.append(renyientropy(counts_distribution, .5))
        profile_values.append(renyientropy(counts_distribution, 1))
        profile_values.append(renyientropy(counts_distribution, 2))
        profile_values.append(renyientropy(counts_distribution, 4))
        profile_values.append(renyientropy(counts_distribution, 8))
        profile_values.append(renyientropy(counts_distribution, 16))
        profile_values.append(renyientropy(counts_distribution, 32))
        profile_values.append(renyientropy(counts_distribution, 64))
        profile_values.append(renyientropy(counts_distribution, np.inf))
    return profile_values


###RENYI PROFILE VIA SEABORN###(10 AUG 2017)###
def renyi_profile_of_level(levelName,countsDF):
    """
    Summary: calculate renyi profile at a taxonomic level
    Args:
        levelName (str): phlyum, class, order, or family
        countsDF (pandas.DataFrame): dataframe of counts at taxonomic level
    Returns:
    """

    alphasList = []
    renyiProfileList = []
    bigGroupList = []
    #write similarly to make_renyiprofile_list
    for i in range(len(countsDF.columns.tolist())):
        renyiProfileList.insert(i, calculate_renyi_profile(countsDF.ix[:, i]))
        alphasList.insert(i, ['0','0.25','0.5','1','2','2^2','2^3','2^4','2^5','2^6', 'Inf'])
        bigGroupList.insert(i , [countsDF.columns[i].split('_')[0]] * 11)

    flatAlphaList = list(itertools.chain.from_iterable(alphasList))
    flatGroupList = list(itertools.chain.from_iterable(bigGroupList))
    flatRenyiList = list(itertools.chain.from_iterable(renyiProfileList))
    levelList = [levelName] * len(flatAlphaList)
    countsDF.columns.str.split('_').str[0]

    renyiDF = pd.DataFrame({'level': levelList,'group':flatGroupList,'alpha':flatAlphaList,'H_alpha':flatRenyiList})
    renyiDF.to_csv(levelName +'_renyi_profile.csv')

    return renyiDF

def plot_renyi(geneDF, level):
    """
    Summary: plots all four renyi profiles for a gene of interest: phylum, class, order, and family
    Args:
        geneList (list): list of dataframes, one in each taxonomic level
    Returns:
    """

    renyiDF = renyi_profile_of_level(level, geneDF)
    sns.set()
    sns.set_context("paper", font_scale=2.5)
    sns.factorplot(x="alpha", y="H_alpha", col="level",  hue="group", data= renyiDF,
                   ci="sd", palette="muted", size=15, aspect=.85)
    return



def get_diversity_measures(geneName, searchString, studyName, level):
    """
    Summary: creates renyi plot, donut, and outputs code for treemap all as visualisations of diversity for gene of interest
    Args:
        geneName (str): actual gene name
        searchString (str): query for gene of interest
        studyName (str): name of study directory
    Returns:
        Renyi plots
        Donut plots
        code to generate treemaps via Brunel
    """

    print('merging pre-processed counts:')
    mergedCounts = merge_all_counts(studyName)
    print('subsetting merged counts table by gene '+ geneName + '...')
    geneDF = subset_by_gene(studyName,searchString, mergedCounts)
    avgGeneDF = output_treemap(geneDF, geneName, studyName)
    geneLvlDF = geneDF.groupby(level).sum()
    geneLvlDF.to_csv(geneName + level  +'_counts.csv')

    print('''
    quantitative measure of diversity: Renyi profile by taxonomic level

                  |
                  |   ^                                    ^
                  |   |                                    |
    Diversity     | richness                            evenness
    [H_alpha(X)]  |   |                                    |
                  |   v                                    v
                  |______________________________________________
                   0                   alpha                 inf

    The Renyi Profile is a generalised version of the Shannon Entropy, where
    the Renyi entropy of order alpha:

                    log(sum(i=1->n, (p_i)^alpha))
    H_alpha(X) =    -----------------------------
                             1 - alpha

    unless alpha = 1, in which case H_alpha(X) is the Shannon Entropy:

    H_1(X) = - sum(i=1->n, p_i*log(p_i))

    If there are multiple communities, they can be ranked using the renyi profile
    (see above- vertical axis allows for partial or complete diversity ordering)

    Is the change in functional abundance across environments due to a change in diversity?
    This provides a quantitative answer.

    ''')
    print(geneName)
    plot_renyi(geneLvlDF, level)
    return avgGeneDF



'''PLOTTING RENYI AND TREEMAP'''


#Renyi profile: using seaborn, reading CSV
#Donut: using Bokeh, reading CSV
#Treemap: using Brunel, reading CSV

#Renyi CSV
#input:
#output: dataframe or csv with columns for each group's renyi profile, and rows for each value in the renyi profilealpha
'''
Renyi will be plotted via seaborn
Group  alpha  H_alpha(X)
------|------|----------
      |  0   |
      |  ... |
      |  inf |

'''
#Brunel Treemap CSV
'''

%brunel data('.csv') treemap x(phylum) y(class) color(phylum) size(counts) sort(counts) label(class) tooltip(#all)
%brunel data('.csv') treemap x(class) y(order) color(class) size(counts) sort(counts) label(order) tooltip(#all)
%brunel data('.csv') treemap x(order) y(family) color(order) size(counts) sort(counts) label(family) tooltip(#all)
'''
'''







'''
#brunel command to plot: data('') x(alpha)
def average_by_groups(grpsList, taxListLevel):
    grpsMeanList = []
    for i in range(len(grpsList.tolist())):
        grpsMeanList.insert(i, taxListLevel.filter(like=grpsList.tolist()[i]).mean(1))
    return grpsMeanList

def output_treemap(geneSorted, geneName, studyName):
    """
    Summary: generates code to make treemaps
    Args:
        geneSorted (pandas.DataFrame): gene counts DF with hierarchy
        geneName (str): actual gene name
        studyName (str): name of study directory
    Returns:
        None, prints out code that can be copied into Jupyter as-is to generate treemaps via Brunel
    """

    cols = pd.Series([c for c in geneSorted.columns if '_'  in  c.lower()])
    grpsList = cols.str.split("_").str[0].unique()

    geneAvg = average_by_groups(grpsList, geneSorted)
    for i in range(len(geneAvg)):
        geneSorted[str(grpsList[i])] = geneAvg[i]
    cols = pd.Series([c for c in geneSorted.columns if '_' not in  c.lower()])
    geneSortedMean = geneSorted[cols]
    geneMeanList = []
    grpsList = grpsList.tolist()
    for k in range(len(geneAvg)):
        geneSortedMean[grpsList[k]] = geneAvg[k]
    donutList = []
    fileName = (studyName + '_' + geneName.replace(' ', '') + '_avg.csv')
    geneSortedMean.to_csv(fileName)
    print('''
    INSTRUCTIONS FOR GENERATING DIVERSITY TREEMAPS
    First copy the import statements and the read_csv into a cell.
    For each figure you want to generate, copy each line starting with "%brunel"  into its own  iPython/Jupyter cell:

    ''')
    variableName = geneName.replace(' ', '')
    print('''
    #copy this into a cell
    import brunel
    import pandas as pd
    {0} = pd.read_csv('{1}')
    '''.format(variableName , fileName))
    for e in range(len(grpsList)):
        print('''
        %brunel data('{0}') treemap x(phylum) y(class) color(phylum) size({1}) sort({1}) label(organism, {1}) tooltip(#all) :: width=2000, height=2000

        %brunel data('{0}') treemap x(class) y(order) color(class) size({1}) sort({1}) label(organism, {1}) tooltip(#all) :: width=2000, height=2000

        %brunel data('{0}') treemap x(order) y(family) color(order) size({1}) sort({1}) label(organism, {1}) tooltip(#all) :: width=2000, height=2000
        '''.format(variableName, grpsList[e], fileName))
        #geneMeanList[k] grpsList[k]
    return geneSortedMean

def search_by_taxa(query, level, dfCounts):
    """
    Summary: gets the subset of a dataframe by the query of interest
    Args:
        query (str): search term of interest, use [space] for OR & | [pipe] for AND
        dfCounts (pandas.DataFrame): counts dataframe
    Returns:
        geneSorted (pandas.DataFrame): dataframe for genes containing query
    """

    queryList = query.replace(" ", "|").replace("+", " ")
    geneSorted = dfCounts[dfCounts[level].str.contains(queryList)==True]
    return geneSorted

def subset_by_taxa(studyName, level, query, testMerge):
    """
    Summary: subsets counts dataframe by query (gene name/gene family name)
    Args:
        studyName (str): name of study directory
        level (str): taxonomic level of interest
        testMerge (pandas.DataFrame): merged annotated counts
    Returns:
        geneDF (pandas.DataFrame): dataframe of counts only for gene/query of interest
    """

    os.chdir(studyName)
    geneDF = search_by_taxa(query, level, testMerge)
    geneDF.to_csv(str(query)+'_'+'allcounts.csv')
    os.chdir('..')
    return geneDF


