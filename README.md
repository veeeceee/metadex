MetaDEx: The Metagenomic Data Explorer
=======

Comparative metagenomics is based on determining features of interest from differentially abundant counts. Differential abundance does not imply importance to one's study, and can be the result of differences in diversity. A key question of interest in comparative metagenomics is to identify the role of environment on the functionome. Essentially, which gene functions are undergoing selection pressures in a given environment?

MetaDEx is a simple Python package for determining environmentally relevant functions and subsystems for comparative metagenomic data. The package receives annotations from MG-RAST and continues with a straightforward workflow for determining genes and subsystems of interest.


Get Study Information from MG-RAST
----------------------------------

MG-RAST assigns each metagenome a metagenome ID (mgm...). For each
metagenome you would like to include in your study, you will need:
- metagenome ID 
- corresponding assigned group 
- corresponding assigned sample ID (within the group can be 1, 2, 3; or something more specfic)

Download Gene Counts Matrix from MG-RAST
----------------------------------------

1.  Go to **Analyze** section of MG-RAST
2.  Select all metagenomes in project
3.  Select matrix
4.  Export to CSV (tab-separated CSV)

Getting Annotations from MG-RAST
--------------------------------

Based on the information above, we can download the full counts (the
main information) directly via API:

Here is an example for a public data set from the MG-RAST data repository:

    metadex.get_all_async('lagoon study', {'mgm4739968.3':'nador_lagoon', 'mgm4739969.3':'oualidia_lagoon', 'mgm4739970.3': 'oualidia_lagoon', 'mgm4739971.3':'nador_lagoon'}, 'RefSeq', evalue=5, identity=60, length=15) 

To see how metadex handles this data set, you can find the full example in the accompanying Jupyter Notebook.


This function asynchronously calls the MG-RAST API to get all metagenomic annotations for a study from the MG-RAST servers and stores each sample as a CSV in the study directory.

### Input
Args:
- study (str): the name of the directory where all the counts will be located
- metagenomeGroupDict (dict): dictionary of the format {MGRAST id: group}
- source (str): annotation source, e.g. RefSeq, GenBank, SEED. See MG-RAST API docs
- evalue (int): e-value cut-off. See MG-RAST API docs
- identity (int): percentage identity cut-off for annotations. See MG-RAST API docs
- length (int): length parameter for MG-RAST API call. See MG-RAST API docs


Output
------
Returns:
        None, saves all annotations in folder named study


Reading Annotations from Study into Metadex
-------------------------------------------

    metadex.Study.read_annotations(studyName)

### Input

studyName should be the name of the directory containing individual
counts files, each with the formatting of \*\[group\]\_\[ID\]*.tsv
where* represents wildcard. If you used the get annotations function,
'folder' should be the same as the name of your study

### Output

Load Gene Counts Matrix into Study
----------------------------------

Load the Gene counts matrix, because it will help with finding genes of
interest:

    metadex.Study.load_gene_matrix(filename)

Remember the gene counts matrix you downloaded from MG-RAST's 'Analyze'
interface? Load the separated text file using its location.

Loading Study Counts into Metadex
---------------------------------

If you have counts files, you can load them into a Metadex study:

    metadex.Study.load_counts(folder) 

### Input

folder should be the name of the directory containing individual counts
files, each with the formatting of \*\[group\]\_\[ID\]*.csv where*
represents wildcard.

Can accept any folder name, so can load counts data that has been
annotated and/or normalised.

### Output


Creates a Study, which contains:

Study.counts =&gt; list of counts DataFrames Study.groups =&gt; list of
corresponding groups Study.iDs =&gt; list of corresponding ID within the
groups

Normalising Counts
------------------

Counts need to be normalised. Metadex takes care of this too:

    metadex.Study.normalise_counts()

### Input

Study.counts

### Output

normalises Study.counts via rarefaction and recodification Study.counts
will be replaced with the rarefied and recodified version

Annotating with Lineage Information
-----------------------------------

The counts data has species info, but this may obscure trends that occur
along higher taxonomic levels. Thus, metadex updates the counts data
with correspondin information for phylum, class, order, and family:

    metadex.Study.add_lineage_info()

This step requires persistent communication with the Entrez servers, and
will take a while. Thankfully it outputs each query to the console.
\#\#\# Input Study.counts

### Output

Study.counts will have additional columns denoting phylum, class, order,
and family for each annotation

Determining Subsystems of Interest
----------------------------------

While MG-RAST allows you to visualise genes that have changed using
their 'Analyze' interface, Metadex leverages the Boruta feature
selection method to determine differences between groups and samples and
each level of functional hierarchy.

Boruta uses a specific threshold for feature selection (see [here]
(http://danielhomola.com/2015/05/08/borutapy-an-all-relevant-feature-selection-method/)) and each can be manipulated as per your needs.

    metadex.Study.determine_subsystems_of_interest(level, percentage)

At each functional level, Boruta will provide selected and weakly
selected features. These should help guide the user as to potential
genes of interest that merit further examination.

Focusing on Gene (User Query)
-----------------------------

Once the user has a potential gene or gene family of interest in mind,
they can zoom in to that subset of the data.

    metadex.Study.focus_on_gene(searchString) 

Visualising Diversity for Gene
------------------------------

Understanding how a gene seen in an environment is distributed is a key
insight to understanding the link between gene function and environment.
Metadex provides ways to depict both the quantitative and qualitative
facts of this relationship within one' study:

    metadex.Study.visualise_diversity_for_gene(geneName, searchString)
