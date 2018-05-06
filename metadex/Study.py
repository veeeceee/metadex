
from .preprocessing import *
from .analysis import *
from pkg_resources import resource_string
import numbers

class Study:

    def __init__(self, name=None, counts=None, groups=None, iDs=None, geneMatrix=None, interestingLvl1=None, interestingLvl2=None, interestingLvl3=None, interestingFxn=None, interestingPhylum=None, interestingClass=None, interestingOrder=None, interestingFamily=None):
        self.name = name

    def read_annotations(self):
        self.counts = get_group_counts_API(self.name)
        filenames = glob.glob(self.name +'/*.csv')
        self.groups = [filename.split('/')[1].split('_')[0] for filename in filenames]
        self.iDs = [filename.split('/')[1].split('_')[1] for filename in filenames]



    def normalise_counts(self):
        normalise_rnr(self.name)
        normalisedList = glob.glob(os.path.join('normalised_counts/*.csv'))
        normalisedDFList = [pd.read_csv(norm) for norm in normalisedList]
        self.counts = normalisedDFList

    def load_counts(self, folder):
        filenames = glob.glob(folder +'/*.csv')
        self.counts = [pd.read_csv(filename) for filename in filenames]
        self.groups = [filename.split('/')[-1].split('_')[0] for filename in filenames]
        self.iDs = [filename.split('/')[-1].split('_')[1] for filename in filenames]

    def add_lineage_info(self):
        self.counts = annotate_all_counts(self.name)


    def load_gene_matrix(self, filename):
        self.geneMatrix = pd.read_table(filename, header=0)

    def determine_genes_of_interest(self, lvl1pct=70, lvl2pct=70, lvl3pct=60, fxnpct=40):
        find_genes_of_interest(self.name, self.groups, self.geneMatrix, lvl1pct, lvl2pct, lvl3pct, fxnpct)
        lvl1List = [pd.read_csv(i) for i in glob.glob(self.name + '/level1*tentative*.csv')]
        lvl1DF = pd.concat(lvl1List, ignore_index=True)
        lvl1DF.columns = lvl1DF.columns.str.split('_').str[0]
        self.interestingLvl1 = lvl1DF.groupby(by=lvl1DF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        lvl2List = [pd.read_csv(i) for i in glob.glob(self.name + '/level2*tentative*.csv')]
        lvl2DF = pd.concat(lvl2List, ignore_index=True)
        lvl2DF.columns = lvl2DF.columns.str.split('_').str[0]
        self.interestingLvl2 = lvl2DF.groupby(by=lvl2DF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        lvl3List = [pd.read_csv(i) for i in glob.glob(self.name + '/level3*tentative*.csv')]
        lvl3DF = pd.concat(lvl3List, ignore_index=True)
        lvl3DF.columns = lvl3DF.columns.str.split('_').str[0]
        self.interestingLvl3 = lvl3DF.groupby(by=lvl3DF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        fxnList = [pd.read_csv(i) for i in glob.glob(self.name + '/level4*tentative*.csv')]
        fxnDF = pd.concat(fxnList, ignore_index=True)
        fxnDF.columns = fxnDF.columns.str.split('_').str[0]
        self.interestingFxn = fxnDF.groupby(by=fxnDF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])

    def determine_subsystems_of_interest(self, level, percentage):
        find_subsystems_of_interest(self.name, self.groups, self.geneMatrix, level, percentage)
        lvlList = [pd.read_csv(i) for i in glob.glob(self.name + '/' + level + '*tentative*.csv')]
        lvlDF = pd.concat(lvlList, ignore_index=True)
        lvlDF.columns = lvlDF.columns.str.split('_').str[0]
        if level == 'level1':
            self.interestingLvl1 = lvlDF.groupby(by=lvlDF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        elif level == 'level2':
            self.interestingLvl2 = lvlDF.groupby(by=lvlDF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        elif level == 'level3':
            self.interestingLvl3 = lvlDF.groupby(by=lvlDF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        elif level == 'function':
            self.interestingFxn = lvlDF.groupby(by=lvlDF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])

    def determine_taxa_of_interest(self, level, percentage):
        find_taxa_of_interest(self.name, self.groups, level, percentage)
        lvlList = [pd.read_csv(i) for i in glob.glob(self.name + '/' + level + '*tentative*.csv')]
        lvlDF = pd.concat(lvlList, ignore_index=True)
        lvlDF.columns = lvlDF.columns.str.split('_').str[0]
        if level == 'phylum':
            self.interestingPhylum = lvlDF.groupby(by=lvlDF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        elif level == 'class':
            self.interestingClass = lvlDF.groupby(by=lvlDF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        elif level == 'order':
            self.interestingOrder = lvlDF.groupby(by=lvlDF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        elif level == 'family':
            self.interestingFamily = lvlDF.groupby(by=lvlDF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])

    def determine_organisms_of_interest(self, lvl1pct=70, lvl2pct=70, lvl3pct=60, lvl4pct=40):
        find_organisms_of_interest(self.name, self.groups, lvl1pct, lvl2pct, lvl3pct, lvl4pct)
        lvl1List = [pd.read_csv(i) for i in glob.glob(self.name + '/phylum*tentative*.csv')]
        lvl1DF = pd.concat(lvl1List, ignore_index=True)
        lvl1DF.columns = lvl1DF.columns.str.split('_').str[0]
        self.interestingPhylum = lvl1DF.groupby(by=lvl1DF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        lvl2List = [pd.read_csv(i) for i in glob.glob(self.name + '/class*tentative*.csv')]
        lvl2DF = pd.concat(lvl2List, ignore_index=True)
        lvl2DF.columns = lvl2DF.columns.str.split('_').str[0]
        self.interestingClass = lvl2DF.groupby(by=lvl2DF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        lvl3List = [pd.read_csv(i) for i in glob.glob(self.name + '/order*tentative*.csv')]
        lvl3DF = pd.concat(lvl3List, ignore_index=True)
        lvl3DF.columns = lvl3DF.columns.str.split('_').str[0]
        self.interestingOrder = lvl3DF.groupby(by=lvl3DF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])
        lvl4List = [pd.read_csv(i) for i in glob.glob(self.name + '/family*tentative*.csv')]
        lvl4DF = pd.concat(lvl4List, ignore_index=True)
        lvl4DF.columns = lvl4DF.columns.str.split('_').str[0]
        self.interestingFamily = lvl4DF.groupby(by=lvl4DF.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])

    def examine_subsystem(self, queryLevel, query, outputLevel):
        return self.geneMatrix[self.geneMatrix[queryLevel] == query].groupby(outputLevel).sum()

    def focus_on_gene(self, searchString):
        mergedCounts = merge_all_counts(self.name)
        mergedCounts = mergedCounts.fillna(0)
        geneDF = subset_by_gene(self.name, searchString, mergedCounts)
        return geneDF

    def visualise_diversity_for_gene(self, geneName, searchString, level='phylum'):
        get_diversity_measures(geneName, searchString, self.name, level)

    def focus_on_organisms(self, level, query):
        mergedCounts = merge_all_counts(self.name)
        mergedCounts = mergedCounts.fillna(0)
        orgDF = subset_by_taxa(self.name, level, query, mergedCounts)
        return orgDF

    def get_all_fxns_for_organisms(self, level, query):
        return self.focus_on_organisms(level, query).groupby([level, 'gene function'], as_index=False).sum().groupby('gene function', as_index = False).sum()

    def focus_on_fxn_in_organisms(self, level, orgQuery, geneQuery):
        allFxnTable = self.get_all_fxns_for_organisms(level, orgQuery)
        return allFxnTable[allFxnTable['gene function'].str.contains(geneQuery)]
