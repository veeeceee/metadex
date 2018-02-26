import glob
import math
import numpy as np
import pandas as pd
import sys
import os
import errno
from skbio.stats import subsample_counts
from functools import reduce
from Bio import Entrez
import time
import re
from scipy import stats
import pkg_resources
pd.options.mode.chained_assignment = None
############################################################################################

def convert_annotation_to_counts(tableName, groupName, sampleID):
    """
    Summary: converts annotation files gotten from API to non-normalised counts data for both function and organism

    Args:
        tableName (str): sample's annotation table as gotten from MG-RAST API
        groupName (str): sample's group
        sampleID  (str): sample's ID


    Returns:
        counts (pandas.DataFrame): pandas DataFrame with columns  'gene function', 'organism', and 'count'
    """
    print('reading ' + str(tableName))
    rawData = pd.read_table(tableName)
    rawData['organism'] = rawData['semicolon separated list of annotations'].str.extract("org\w+\=\[([^\]]+)\]", expand = False)
    rawData['gene function'] = rawData['semicolon separated list of annotations'].str.extract("func\w+\=\[([^\]]+)\]", expand = False)
    print('filtering by max percentage identity...')
    rawData = keep_max_pctID_only(rawData)
    counts = get_and_save_counts(rawData, groupName, sampleID)
    return counts

def get_group_counts_API(studyName):
    """
    Summary: converts all annotations in a study, with group awareness, to raw counts data


    Args:
        studyName (str): name of the directory under which annotations are stored

    Returns:
        grpsCountsList (list): list containing each counts dataFrame for each sample. This is also output to CSV format in a folder 'counts'

    """
    os.chdir(str(studyName))
    print('loading files in ' + studyName)
    metagenomeList = glob.glob('*.tsv')
    filenames = [metagenome.split('.')[0] for metagenome in metagenomeList]
    groupName = [name.split('_')[0] for name in filenames]
    sampleID = [name.split('_')[1] for name in filenames]
    print('converting annotations to counts...')
    grpCountsList = [convert_annotation_to_counts(metagenomeList[i], groupName[i], sampleID[i]) for i in range(len(metagenomeList))]
    os.chdir('..')
    return grpCountsList

############################################################################################
'''
TRANSFORMING ANNOTATIONS INTO NORMALISED COUNTS DATA











'''
############################################################################################

####CREATE PATH IN DIRECTORY##
def create_path(path):
    """
    Summary: creates a path, with error checking

    Args:
        path (str): name of folder/directory to be created

    Returns:
        None, creates folder
    """

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

####READ CSV FILE INTO PANDAS DATAFRAME####
def make_df(filename):
    """
    Summary:

    Args:

    Returns:
    """

    #df = pd.read_csv(filename, header=None)
    #return df

    df = pd.read_csv(filename) #index_col=[0,1]
    return df


###KEEP ANNOTATIONS WITH MAX PERCENT IDENTITY###
def keep_max_pctID_only(table):
    """
    Summary: filters annotations by maximum percentage identity

    Args:
        table (pandas.DataFrame): annotations table, as DataFrame

    Returns:
        maxPctID (pandas.DataFrame): filtered annotations table, as DataFrame
    """
    mask = table.groupby('query sequence id').agg('idxmax')
    maxPctID = table.loc[mask['percentage identity']].reset_index()
    return maxPctID

###MERGE FUNCTION AND ORGANISM ANNOTATION TABLES###
def merge_annotation_tables(table1, table2):
    """
    Summary:

    Args:

    Returns:
    """
    mergedTable = pd.merge(table1, table2, on = 'query sequence id').sort_values(['semicolon separated list of annotations_x', 'semicolon separated list of annotations_y'], ascending = [0,0])
    mergedTable.rename(columns={'semicolon separated list of annotations_x':'gene function','semicolon separated list of annotations_y':'organism'})
    groupByQuery = mergedTable.groupby('query sequence id')
    return mergedTable

###GET RAW COUNTS FOR MERGED ANNOTATION TABLES###
def get_and_save_counts(mergedTable, groupName, sampleID):
    """
    Summary:

    Args:

    Returns:
    """
    print('sorting values...')
    counts1 = mergedTable.sort_values(['gene function', 'organism']).groupby(['gene function', 'organism']).size()
    newFileName = '%(group)s_%(sample)s_counts.csv' % \
            {"group": groupName, "sample": sampleID}
    print('writing ' + newFileName  +  ' to file')
    counts1.to_csv(newFileName)
    print('written ' + str(newFileName))
    #counts = pd.read_csv(newFileName)
    counts1.columns = ['gene function', 'organism', groupName + '_' + sampleID]
    return counts1

####GO FROM TWO TABLES (FXN AND ORG) TO COUNTS####
def counts_within_metagenome(fxnTable, orgTable, groupName, sampleID):
    """
    Summary:

    Args:

    Returns:
    """
    mergedTable = merge_annotation_tables(fxnTable, orgTable)
    counts = get_and_save_counts(mergedTable, groupName, sampleID)
    return counts


#####CONVERTING TAB ANNOTATIONS FROM MG-RAST INTO RAW COUNTS DATA#####
def get_counts_within_metagenome(fxnTable, orgTable, groupName, sampleID):
    """
    Summary:

    Args:

    Returns:
    """
    fxnData = pd.read_csv(fxnTable, delimiter = '\t', usecols=['query sequence id', 'percentage identity', 'bit score', 'semicolon separated list of annotations'])
    orgData = pd.read_csv(orgTable, delimiter = '\t', usecols=['query sequence id', 'percentage identity', 'bit score', 'semicolon separated list of annotations'])

    fxnDataMax = keep_max_pctID_only(fxnData)
    orgDataMax = keep_max_pctID_only(orgData)

    counts = counts_within_metagenome(fxnDataMax, orgDataMax, groupName, sampleID)

    return counts

###CONVERT TAB ANNOTATIONS TO COUNTS WITHIN A GROUP####
def get_counts_within_group(studyName, groupName):
    """
    Summary:

    Args:

    Returns:
    """
    fxnDataList = glob.glob(os.path.join(studyName + '/' + groupName, '*function*.tab'))
    orgDataList = glob.glob(os.path.join(studyName + '/' + groupName, '*organism*.tab'))
    grpCountsList = [get_counts_within_metagenome(fxnDataList[i], orgDataList[i], groupName, str(i+1)) for i in range(len(fxnDataList))]
    return grpCountsList

#####GET ALL THE COUNTS (FROM TAB ANNOTATIONS) OF ALL THE GROUPS IN A STUDY#####
def get_all_group_counts(studyName):
    """
    Summary:

    Args:

    Returns:
    """
    grpsList = [poo.split('/')[1] for poo in glob.glob(os.path.join(studyName, '*'))]
    fullCountsList = [get_counts_within_group(studyName, grp) for grp in grpsList]
    return fullCountsList
########################################NORMALISATION###############################################
'''Normalisation by subsampling and recodification'''
###FIND SAMPLING DEPTH FOR SUBSAMPLING###
def find_sampling_depth(rawCounts):
    """
    Summary: determines median counts of all samples and sets as sampling depth

    Args:
        rawCounts (list): list of all raw counts DataFrames, each corresponding to a sample

    Returns:
        samplingDepth (int): median number of counts across all samples
    """
    sums = []
    print('finding sampling depth...')
    for i in range(len(rawCounts)):
        rawCounts[i] = rawCounts[i].set_index(['gene function', 'organism'])
        sums.insert(i, int(rawCounts[i].sum()))
    lenSeries = pd.Series(sums)
    samplingDepth = int(lenSeries.median())
    return samplingDepth

####WRITE RAREFIED AND RECODED COUNTS TO FILE####
def rarefy_and_recode(filenames, rawCounts, samplingDepth):
    """
    Summary: subsamples all samples to the median (average of 100 times),  

    Args:
        filenames ()
        rawCounts ()
        samplingDepth ()

    Returns:
    """
    for i in range(len(rawCounts)):
        subsampleList = []
        if int(rawCounts[i].sum()) < samplingDepth:
            meanSubsample = rawCounts[i]
        else:
            for j in range(100):
                sample = subsample_counts(rawCounts[i].transpose().values[0], samplingDepth)
                subsampleList.insert(j, sample)
            print("completed 100 subsamples for sample number " + str(i))
            meanSubsample = pd.Series(subsampleList).mean()
            #recodification: setting all values less than 1.01 to zero
            meanSubsample[meanSubsample < 1.01] = 0
        sampleName = filenames[i].split('.')[0]
        rawCounts[i][sampleName] = meanSubsample
        newFileName = sampleName + "_norm.csv"
        create_path('normalised_counts')
        rawCounts[i].to_csv(os.path.join('normalised_counts', newFileName))
        print("written " + newFileName + " to file.")
    return

#####NORMALISE BY R&R (TO MEDIAN)#####
def normalise_rnr(studyName):
    """
    Summary: normalises all counts in study by rarefaction and recodification

    Args:
        studyName (str): name of directory/study containing raw non-normalised counts

    Returns: None, outputs normalised counts into studyName/normalised_counts
    """
    os.chdir(str(studyName))
    print('loading files in ' + studyName)
    metagenomeList = glob.glob('*counts.csv')
    flatCountsList = [pd.read_csv(metagenome, names=['gene function', 'organism', 'counts']) for metagenome in metagenomeList]
    #flatCountsList = sum(fullCountsList, [])
    samplingDepth = find_sampling_depth(flatCountsList)
    filenames = glob.glob(os.path.join('counts/*.csv'))
    rarefy_and_recode(metagenomeList, flatCountsList, samplingDepth)
    os.chdir('..')
    return
############################################################################################
'''
ANNOTATING COUNTS DATA

Communicates with Entrez E-utilities to annotate hierarchies









'''
############################################################################################
''' Insert the proper version of the get_taxid function '''
def get_taxid(species):
    """
    Summary: gets Entrex taxid for given lineage name 

    Args:
        species (str): species, genus, etc

    Returns:
        (str): taxid for given name
    """
    species = species.split(';')[0].replace(" ", "+").replace("''", " ").replace(":", " ").replace("(", " ").replace(")", " ")#.split('+')
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml", usehistory = "n")
    record = Entrez.read(search)
    if (record['IdList'] == []):
        name_fixed = record['ErrorList']['PhraseNotFound'][0]#.split("+")
        #redo = get_taxid(name_fixed[0] + " " + name_fixed[1])  ''' make this function non-recursive'''
        search2 = Entrez.esearch(term = name_fixed, db = "taxonomy", retmode = "xml")
        record = Entrez.read(search2)
        if (record['IdList'] == []):
            thirdTry = species.split("+")[0] + " " + species.split("+")[1]
            search3 = Entrez.esearch(term = thirdTry, db = "taxonomy", retmode = "xml")
            record = Entrez.read(search3)
        else:
            record = record
    else:
        record = record

    return record['IdList'][0]

def get_tax_data(taxid):
    """
    Summary: access lineage in Entrez database for given taxid

    Args:
        taxid (str): the taxid for which lineage information is desired

    Returns:
        JSON entry for Entrez taxid

    """
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)


def get_tax_hierarchy(species):
    """
    Summary: return taxonomic hierarchy of species

    Args:
        species (str): species, genus, etc for which lineage is desired

    Returns:
        lineage (dict): dictionary with {taxonomic level : name}
    """
    try:
        taxid = get_taxid(species)
        data = get_tax_data(taxid)
    except:
        time.sleep(50)
        taxid = get_taxid(species)
        data = get_tax_data(taxid)
    lineage = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in ['phylum', 'class', 'order', 'family', 'genus']}
    return lineage





#HIERARCHY ANNOTATION#
#input is a series of species names







def get_taxonomic_hierarchy(organismSeries, userEmail):
    """
    Summary: gets all the taxonomic hierarchies/lineages for a series of organism names

    Args: 
        organismSeries (pandas.Series): Series containing  organisms
        userEmail (string): email of tool developer

    Returns:
    """

    Entrez.email = userEmail
    Entrez.tool = "Vibhu-species2taxon"
    if not Entrez.email:
        print("you must add your email address")
        sys.exit(2)

    lineage_list = []

    print('parsing taxonomic data...') #starting the parser

    for species in organismSeries:
        print ('\t' + species.split(";")[0]) #print progress
        lineage = get_tax_hierarchy(species)
        lineage_list.append(str(lineage))
        time.sleep(5)

    return lineage_list


##ANNOTATE ENTIRE METAGENOMES DATA FRAME WITH HIERARCHY##
#input is a dataframe in which the numerical columns are each samples
def annotate_with_hierarchy(dfCounts):
    """
    Summary: annotates normalised counts with lineage for each organism

    Args:
        dfCounts (pandas.DataFrame): counts dataframe

    Returns:
        dfCounts (pandas.DataFrame): counts dataframe with columns for taxonomic hierarchy
    """
    io = pkg_resources.resource_stream("metadex", "lineage_dict.csv")
    dfCounts['lineage'] = dfCounts['organism'].str.split(";").str[0]

    speciesList = dfCounts['lineage'].unique().tolist()
    speciesSeries = pd.Series(speciesList)
    lineageDF = pd.read_csv(io)
    dictSpeciesList = lineageDF['species'].tolist()
    diffList = list(set(speciesList)-set(dictSpeciesList))
    initlineageDict = dict(zip(lineageDF['species'].tolist(), lineageDF['lineage'].tolist()))



    lineageList = get_taxonomic_hierarchy(pd.Series(diffList), 'vibhuc@me.com')

    difflineageDict = dict(zip(diffList, lineageList))
    lineageDict = initlineageDict.copy()
    lineageDict.update(difflineageDict)
    lineageDF = pd.DataFrame({'name':lineageDict.keys(), 'lineage':lineageDict.values()})
    lineageDF.to_csv('updated_lineage_dict.csv')
    dfCounts['lineage'].replace(lineageDict, inplace = True)

    dfCounts['phylum'] = dfCounts['lineage'].str.extract("'p\w+\'\: \'([A-z]\w+)\'", expand = False)
    dfCounts['class'] = dfCounts['lineage'].str.extract("'c\w+\'\: \'([A-z]\w+)\'", expand = False)
    dfCounts['order'] = dfCounts['lineage'].str.extract("'o\w+\'\: \'([A-z]\w+)\'", expand = False)
    dfCounts['family'] = dfCounts['lineage'].str.extract("'f\w+\'\: \'([A-z]\w+)\'", expand = False)
    dfCounts.pop('lineage')
    return dfCounts

def annotate_all_counts(studyName):
    """
    Summary: annotates all counts within a study with taxonomic hierarchy

    Args:
        studyName (str): name of directory within which counts are located

    Returns:
        None, outputs annotated files to norm_taxonomy folder
    """

    normFilesList = glob.glob(os.path.join(studyName + '/normalised_counts', '*norm.csv'))
    normFileNames = [normFile.strip().split('/')[-1].split('.')[0] for normFile in normFilesList]
    normDFsList = [pd.read_csv(normCount) for normCount in normFilesList]
    [normDF.pop('counts') for normDF in normDFsList]
    annotatedList = [annotate_with_hierarchy(normDF) for normDF in normDFsList]
    os.chdir(studyName)
    os.mkdir('norm_taxonomy')
    os.chdir('norm_taxonomy')
    [normDFsList[i].to_csv(normFileNames[i] + '_taxonomy.csv') for i in range(len(normDFsList))]
    return annotatedList



##############################PRE-PROCESSING FOR GENE-ONLY COUNTS MATRIX##############################
###RUN BORUTA TO DETERMINE WHICH GENES HAVE CHANGES MOST###
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from boruta import BorutaPy

def find_genes_of_interest(studyName, groupsList, geneCounts, lvl1pct=70, lvl2pct=70, lvl3pct=60, fxnpct=40):
    """
    Summary: uses Boruta machine learning method to roughly determine potential genes of interest. requires tab-separated  matrix from MG-RAST analysis page

    Args:
        studyName (str): directory (study name)
        geneCountsName (str): filename for tab separated matrix
        lvl1pct (int): threshold for Boruta on level 1
        lvl2pct (int): threshold for Boruta on level 2
        lvl3pct (int): threshold for Boruta on level 3
        fxnpct (int): threshold for Boruta on gene name


    Returns: None, outputs files with tentative genes/gene families of interest

    """

    #geneCounts = pd.read_table(geneCountsName, header=0)#, header=0)#header=0
    numGeneCounts = geneCounts.select_dtypes(include=[np.number])
    Y = numGeneCounts.transpose().index.str.split('_').str[0].values
    samplingDepth = numGeneCounts.sum().median()
    os.chdir(studyName)
    for i in range(len(numGeneCounts.columns)):
        subsampleList = []
        if int(numGeneCounts[numGeneCounts.columns[i]].sum()) < samplingDepth:
            meanSubsample = numGeneCounts[numGeneCounts.columns[i]]
        else:
            for j in range(100):
                sample = subsample_counts(numGeneCounts[numGeneCounts.columns[i]].transpose().values, int(samplingDepth))
                subsampleList.insert(j, sample)
            print("completed 100 subsamples for sample number " + str(i))
            meanSubsample = pd.Series(subsampleList).mean()
            #recodification: setting all values less than 1.01 to zero
            meanSubsample[meanSubsample < 1.01] = 0
        meanSubsample = 100 *  meanSubsample/meanSubsample.sum()
        numGeneCounts[numGeneCounts.columns[i]] = meanSubsample
    numGeneCounts['level1'] = geneCounts['level1']
    numGeneCounts['level2'] = geneCounts['level2']
    numGeneCounts['level3'] = geneCounts['level3']
    numGeneCounts['function'] = geneCounts['function']
    countsLvl1 = numGeneCounts.groupby('level1').sum()
    countsLvl2 = numGeneCounts.groupby('level2').sum()
    countsLvl3 = numGeneCounts.groupby('level3').sum()
    countsLvl4 = numGeneCounts.groupby('function').sum()
    levelList = [countsLvl1, countsLvl2, countsLvl3, countsLvl4]
    countsLvl1.to_csv(studyName + 'genes_lvl1.csv')
    countsLvl2.to_csv(studyName + 'genes_lvl2.csv')
    countsLvl3.to_csv(studyName + 'genes_lvl3.csv')
    countsLvl4.to_csv(studyName + 'genes_function.csv')
    groupsDict = dict(enumerate(pd.Series(groupsList).unique()))
    dictGroups = {y:x for x,y in groupsDict.items()}
    rf = RandomForestClassifier(n_jobs = -1, class_weight = 'balanced', max_depth=3)

    X = countsLvl1.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(lvl1pct))
    feat_selector.fit(X, Y)
    if len(countsLvl1[feat_selector.support_]) > 0:
        countsLvl1[feat_selector.support_].to_csv('level1_tentative.csv')
    countsLvl1[feat_selector.support_weak_].to_csv('level1_tentative_weak.csv')

    X = countsLvl2.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(lvl2pct), max_iter=300)
    feat_selector.fit(X, Y)
    if len(countsLvl2[feat_selector.support_]) > 0:
        countsLvl2[feat_selector.support_].to_csv('level2_tentative.csv')
    countsLvl2[feat_selector.support_weak_].to_csv('level2_tentative_weak.csv')

    X = countsLvl3.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(lvl3pct), max_iter=500)
    feat_selector.fit(X, Y)
    if len(countsLvl3[feat_selector.support_]) > 0:
        countsLvl3[feat_selector.support_].to_csv('level3_tentative.csv')
    countsLvl3[feat_selector.support_weak_].to_csv('level3_tentative_weak.csv')

    X = countsLvl4.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(fxnpct), max_iter=700)
    feat_selector.fit(X, Y)
    if len(countsLvl4[feat_selector.support_]) > 0:
        countsLvl4[feat_selector.support_].to_csv('level4_tentative_.csv')
    countsLvl4[feat_selector.support_weak_].to_csv('level4_tentative_weak.csv')
    os.chdir('..')

def find_subsystems_of_interest(studyName, groupsList, geneCounts, level, percentage):
    """
    Summary: uses Boruta machine learning method to roughly determine potential genes of interest. requires tab-separated  matrix from MG-RAST analysis page

    Args:
        studyName (str): directory (study name)
        groupsList (list): list of group names
        level (str): subsystems level at which to run Boruta
        percentage (int): threshold for Boruta feature selection


    Returns: None, outputs files with tentative genes/gene families of interest

    """


    numGeneCounts = geneCounts.select_dtypes(include=[np.number])
    Y = numGeneCounts.transpose().index.str.split('_').str[0].values
    samplingDepth = numGeneCounts.sum().median()
    os.chdir(studyName)
    for i in range(len(numGeneCounts.columns)):
        subsampleList = []
        if int(numGeneCounts[numGeneCounts.columns[i]].sum()) < samplingDepth:
            meanSubsample = numGeneCounts[numGeneCounts.columns[i]]
        else:
            for j in range(100):
                sample = subsample_counts(numGeneCounts[numGeneCounts.columns[i]].transpose().values, int(samplingDepth))
                subsampleList.insert(j, sample)
            print("completed 100 subsamples for sample number " + str(i))
            meanSubsample = pd.Series(subsampleList).mean()
            #recodification: setting all values less than 1.01 to zero
            meanSubsample[meanSubsample < 1.01] = 0
        meanSubsample = 100 *  meanSubsample/meanSubsample.sum()
        numGeneCounts[numGeneCounts.columns[i]] = meanSubsample
    numGeneCounts['level1'] = geneCounts['level1']
    numGeneCounts['level2'] = geneCounts['level2']
    numGeneCounts['level3'] = geneCounts['level3']
    numGeneCounts['function'] = geneCounts['function']
    countsLvl = numGeneCounts.groupby(level).sum()
    groupsDict = dict(enumerate(pd.Series(groupsList).unique()))
    dictGroups = {y:x for x,y in groupsDict.items()}
    rf = RandomForestClassifier(n_jobs = -1, class_weight = 'balanced', max_depth=3)
    X = countsLvl.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(percentage))
    feat_selector.fit(X, Y)
    if len(countsLvl[feat_selector.support_]) > 0:
        countsLvl[feat_selector.support_].to_csv(str(level)+'_tentative.csv')
    countsLvl[feat_selector.support_weak_].to_csv(str(level)+'_tentative_weak.csv')
    os.chdir('..')








def merge_dfs1(ldf, rdf):
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

def merge_all_counts2(studyName):
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
    mergedCounts  = reduce(merge_dfs1, annotatedDFsList)
    mergedCounts.fillna(0)
    mergedCounts.to_csv(str(studyName) + '_allcounts.csv')
    return mergedCounts

def find_taxa_of_interest(studyName, groupsList, level, percentage):
    """
    Summary: uses Boruta machine learning method to roughly determine potential lineages of interest. requires tab-separated  matrix from MG-RAST analysis page

    Args:
        studyName (str): directory (study name)
        groupsList (list): list of group names
        level (str): taxonomic level at which to run Boruta
        percentage (int): threshold for Boruta feature selection

    Returns: None, outputs files with tentative genes/gene families of interest

    """


    mergedCounts = merge_all_counts2(studyName)
    mergedCounts = mergedCounts.fillna(0)
    os.chdir(studyName)
    countsLvl = mergedCounts.groupby(level).sum()
    countsLvl.to_csv(studyName + '_' + level + '.csv')
    groupsDict = dict(enumerate(pd.Series(groupsList).unique()))
    dictGroups = {y:x for x,y in groupsDict.items()}

    Y = pd.Series(groupsList).values
    rf = RandomForestClassifier(n_jobs = -1, class_weight = 'balanced', max_depth=3)

    X = countsLvl.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(percentage))
    feat_selector.fit(X, Y)
    if len(countsLvl[feat_selector.support_]) > 0:
        countsLvl[feat_selector.support_].to_csv(level + '_tentative.csv')
    countsLvl[feat_selector.support_weak_].to_csv(level + '_tentative_weak.csv')
    os.chdir('..')

def find_organisms_of_interest(studyName, groupsList,  lvl1pct=70, lvl2pct=70, lvl3pct=60, lvl4pct=40):
    """
    Summary: uses Boruta machine learning method to roughly determine potential genes of interest. requires tab-separated  matrix from MG-RAST analysis page

    Args:
        studyName (str): directory (study name)
        geneCountsName (str): filename for tab separated matrix
        lvl1pct (int): threshold for Boruta on level 1
        lvl2pct (int): threshold for Boruta on level 2
        lvl3pct (int): threshold for Boruta on level 3
        fxnpct (int): threshold for Boruta on gene name


    Returns: None, outputs files with tentative genes/gene families of interest

    """

    #geneCounts = pd.read_table(geneCountsName, header=0)#, header=0)#header=0
    mergedCounts = merge_all_counts2(studyName)
    mergedCounts = mergedCounts.fillna(0)
    os.chdir(studyName)
    countsLvl1 = mergedCounts.groupby('phylum').sum()
    countsLvl2 = mergedCounts.groupby('class').sum()
    countsLvl3 = mergedCounts.groupby('order').sum()
    countsLvl4 = mergedCounts.groupby('family').sum()
    levelList = [countsLvl1, countsLvl2, countsLvl3, countsLvl4]
    countsLvl1.to_csv(studyName + '_phylum.csv')
    countsLvl2.to_csv(studyName + '_class.csv')
    countsLvl3.to_csv(studyName + '_order.csv')
    countsLvl4.to_csv(studyName + '_family.csv')
    groupsDict = dict(enumerate(pd.Series(groupsList).unique()))
    dictGroups = {y:x for x,y in groupsDict.items()}
    Y = pd.Series(groupsList).values
    rf = RandomForestClassifier(n_jobs = -1, class_weight = 'balanced', max_depth=3)

    X = countsLvl1.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(lvl1pct))
    feat_selector.fit(X, Y)
    if len(countsLvl1[feat_selector.support_]) > 0:
        countsLvl1[feat_selector.support_].to_csv('phylum_1tentative.csv')
    countsLvl1[feat_selector.support_weak_].to_csv('phylum_1tentative_weak.csv')

    X = countsLvl2.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(lvl2pct), max_iter=300)
    feat_selector.fit(X, Y)
    if len(countsLvl2[feat_selector.support_]) > 0:
        countsLvl2[feat_selector.support_].to_csv('class_1tentative.csv')
    countsLvl2[feat_selector.support_weak_].to_csv('class_1tentative_weak.csv')

    X = countsLvl3.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(lvl3pct), max_iter=500)
    feat_selector.fit(X, Y)
    if len(countsLvl3[feat_selector.support_]) > 0:
        countsLvl3[feat_selector.support_].to_csv('order_1tentative.csv')
    countsLvl3[feat_selector.support_weak_].to_csv('order_1tentative_weak.csv')

    X = countsLvl4.transpose().values
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, perc=int(lvl4pct), max_iter=700)
    feat_selector.fit(X, Y)
    if len(countsLvl4[feat_selector.support_]) > 0:
        countsLvl4[feat_selector.support_].to_csv('family_1tentative_.csv')
    countsLvl4[feat_selector.support_weak_].to_csv('family_1tentative_weak.csv')
    os.chdir('..')

