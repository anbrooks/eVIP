
# What is eVIP?

Expression-based variant impact phenotyping (eVIP) is an approach to predict functional impacts of mutations by comparing gene expression signatures induced by wild-type vs. mutated ORFs.



# How to cite eVIP?

The method was orginally described in :

High-throughput Phenotyping of Lung Cancer Somatic Mutations
Berger AH*, Brooks AN*, Wu X*, Shrestha Y, Chouinard C, Piccioni F, Bagul M, Kamburov A, Imielinski M, Hogstrom L, Zhu C, Yang X, Pantel S, Sakai R, Watson J, Kaplan N, Campbell JD, Singh S, Root DE, Narayan R, Natoli T, Lahr DL, Tirosh I, Tamayo P, Getz G, Wong B, Doench J, Subramanian A, Golub TR, Meyerson M, Boehm JS.
Cancer Cell, Aug 8;30(2):214-28

# Requirements

1) cmap/io/gct.py (https://github.com/cmap/l1ktools): Classes to interact with .gct and .gctx files.

    Tip: These tools also have many requirements. The Anaconda Python install contains most of these requirements. The module blessings can be installed by running:
     
    conda install -c https://conda.anaconda.org/mutirri blessings
    
2) rpy2: http://rpy2.bitbucket.org/

# Overview of eVIP pipeline

1) **eVIP_corr.py** : Uses a gene expression table (should be normalized and filtered, can use filterGeneExpressionTable.py) and creates a spearman correlation matrix and calculates zscores.  If using  -zscore_gct, not --infile, the spearman correlation  matrix will be calculated from the inputted zscores. 
  
2) **eVIP_compare.py** : Uses the spearman correlation rank matrix to perform calculatins.

3) **eVIP_predict.py** : Uses the  calculations in the  algorithm to generate the predictions.
  
4+5) **eVIP_viz.py** and **VIP_sparkler.py** : visualization 

The options for these scripts are included in run_eVIP.py 

For eVIP pathway analysis, use **-eVIPP** option in run_eVIP.py. This  runs eVIP  to get functional predictions for pathways in a variant. 



# Quickstart

To run the eVIP pipeline, use run_eVIP.py. We recommend the input data (--infile) has the low expressed genes filtered out and is log2 transformed.  This can be done by using the filterGeneExpression.py

Usage: 
    
    run_eVIP.py [-h] [--infile INFILE] [-zscore_gct ZSCORE_GCT] 
        -out_directory OUT_DIRECTORY -sig_info SIG_INFO -c C -r R -num_reps NUM_REPS 
        [-ie_filter IE_FILTER] [-ie_col IE_COL] [-i I] [-allele_col ALLELE_COL] [-conn_null CONN_NULL] 
        [-conn_thresh CONN_THRESH][-mut_wt_rep_thresh MUT_WT_REP_THRESH][-disting_thresh DISTING_THRESH] 
        [-mut_wt_rep_rank_diff MUT_WT_REP_RANK_DIFF] [-use_c_pval] [-cell_id CELL_ID] [-plate_id PLATE_ID] 
        [-ref_allele_mode] [-x_thresh X_THRESH] [-y_thresh Y_THRESH] [-annotate] [-by_gene_color BY_GENE_COLOR] 
        [-pdf PDF] [-xmin XMIN] [-xmax XMAX] [-ymin YMIN] [-ymax YMAX] [-viz_ymin VIZ_YMIN][-viz_ymax VIZ_YMAX] 
        [-corr_val CORR_VAL] [-eVIPP][-JSON JSON] [-min_genes MIN_GENES] [-viz_off][-sparkler_off]
                    
Example:

python run_eVIP.py --infile /path/to/data/data_geneExp_filtered.txt -out_directory /path/to/data/eVIP_out -sig_info /path/to/data/sig.info -c /path/to/data/controls.grp -r /path/to/data/comparisons.tsv -num_reps 3 -i 100 -allele_col allele -x_thresh 1.3 -y_thresh 1.3 -by_gene_color /path/to/data/gene_label.tsv -ymin -2 -ymax 2 -corr_val "spearman"

Arguments:

      -h, --help            show this help message and exit
      --infile INFILE       Input txt file (filtered and log transformed data).
      -zscore_gct ZSCORE_GCT
                            Zscore input gct file (use instead of --infile)
      -out_directory OUT_DIRECTORY
                            Path to directory for eVIP output files
      -sig_info SIG_INFO    sig info file with gene information and distil
                            information
      -c C                  .grp file containing allele names of control
                            perturbations. If this file is given, a null will be
                            calculated from these
      -r R                  File explicitly indicating which comparisons to do.
                            Assumes the file has a header and it is ignored. The
                            first column is the reference allele and second column
                            is test allele. If this file is not given, then the
                            reference alleles are assumed to be WT and inferred
                            from the allele names.
      -num_reps NUM_REPS    Number of replicates expected for each allele. DEF=3
      -ie_filter IE_FILTER  Threshold for infection efficiency. Any wildtype or
                            mutant alleles having an ie below this threshold, will
                            be removed
      -ie_col IE_COL        Name of the column in the sig_info file with infection
                            efficiency information. DEF=x_ie_a549
      -i I                  Number of iterations to run. DEF=1000
      -allele_col ALLELE_COL
                            Column name in sig_info file that indicates the allele
                            names.DEF=x_mutation_status
      -conn_null CONN_NULL  Optional file containing connectivity null values from
                            a previous run. Should end in _conn_null.txt
      -conn_thresh CONN_THRESH
                            P-value threshold for connectivity vs null.
                            DEFAULT=0.05
      -mut_wt_rep_thresh MUT_WT_REP_THRESH
                            P-value threshold for comparison of WT and mut
                            robustness. DEFAULT=0.05
      -disting_thresh DISTING_THRESH
                            P-value threshold that tests if mut and wt reps are
                            indistinguishable from each other. DEFAULT=0.05
      -mut_wt_rep_rank_diff MUT_WT_REP_RANK_DIFF
                            The minimum difference in median rankpoint WT and mut
                            to consider a difference. DEF=0
      -use_c_pval           Will use corrected p-value instead of raw p-val
      -cell_id CELL_ID      Optional: Will only look at signatures from this cell
                            line. Helps to filter sig_info file.
      -plate_id PLATE_ID    Optional: Will only look at signatures from this plate
      -ref_allele_mode      Sparkler+Viz: Instead of organizing plots by gene,
                            will use the wt column to determine what are the
                            reference alleles.
      -x_thresh X_THRESH    Sparkler: Threshold of significance
      -y_thresh Y_THRESH    Sparkler: Threshold of impact direction
      -annotate             Sparkler: Will add allele labels to points.
      -by_gene_color BY_GENE_COLOR
                            Sparkler: File containing labels and colors for gene-
                            centric plot.
      -pdf PDF              Sparkler + Viz: Will print plots in pdf format instead
                            of png.
      -xmin XMIN            Sparkler: Min value of x-axis. DEF=0
      -xmax XMAX            Sparkler: Max value of x-axis. DEF=4
      -ymin YMIN            Sparkler: Min value of y-axis. DEF=-3
      -ymax YMAX            Sparkler: Min value of y-axis. DEF=3
      -viz_ymin VIZ_YMIN    Viz: Minimum y-value of rep value. DEF=-100
      -viz_ymax VIZ_YMAX    Viz: Maximum y-value of rep value. DEF=100
      -corr_val CORR_VAL    Viz: String used to label the correlation value. DEF=
                            'row median rankpoints'
      -eVIPP                Use this option when doing pathway analysis, must also
                            have JSON file
      -JSON JSON            JSON file created by create_pathway_JSON.py. Contains
                            dictionary of pathways and the associated ids
      -min_genes MIN_GENES  Minimum amount of pathway genes found in data to run
                            eVIPP on. DEF = 5
      -viz_off              Will not perform eVIP viz step
      -sparkler_off         Will not perform eVIP sparkler step

# Required Files

Required files: 
1. data (zscores or spearman correlation matrix input (.gct format))
2. sig_info
3. reference test file

For eVIPP, a JSON  of pathway to gene information is also required. This can be created by using create_pathway_JSON.py to convert MsigDB gmt files (http://software.broadinstitute.org/gsea/downloads.jsp) to JSON format.  

# Sample Files

Following are examples of for an expirement looking to characterize a variant called "GENE_MUT". 
There are 4 replicates of the "GENE_MUT"(the variant we are testing), WT gene ("GENE_WT"), and GFP (the control).

### Sig Info (-sig_info)

An example sig_info file of an experiment for a prediction of "GENE_MUT"
4 replicates, a GFP control


```python
distil_id	sig_id	pert_mfc_desc	cell_id	allele	293_ie
GENE_WT_4|GENE_WT_3|GENE_WT_2|GENE_WT_1	GENE_WT	GENE	293	GENE_WT	1.0
GENE_MUT_4|GENE_MUT_3|GENE_MUT_2|GENE_MUT_1	GENE_MUT	GENE	293	GENE_MUT	1.0
GFP_4|GFP_3|GFP_2|GFP_1	GFP	GFP	293	GFP	1.0
```

###Reference File (-r)

An example reference file (.tsv)


```python
wt	mutant
GENE	GENE_MUT
```

###Controls File (-c)

In this example experiment, the file would just contain GFP. (.grp)


```python
GFP
```

### Gene Labels (-by_gene_color)

This file is used  to group genes into categories when creating sparkler plots. 
The  gene labels are:

    TSG(tumor suppressor gene)
    ONC(oncogene)
    TSG_noTP53
    ONC-NEG (oncogene-negative)


```python
gene	label
GENE	TSG
```

#Tips

The scripts can also be run individually. For each python script, you can view all options by simply typing: python <script.py> -h 

If plots made by eVIP_viz.py are blank, adjust the min and max.

# Coming soon...

Tutorial and test files coming soon.
