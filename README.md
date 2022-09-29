# Tutorial_for_Aboretum-based_DE_geneset_analysis
A tutorial which includes Arboretum &amp; findTransitionGenesets with an example dataset.

- [0. About this tutorial](#step-0-about-this-tutorial)
- [1. Prepare input tree](#step-1-prepare-input-tree)
- [2. Prepare input order](#step-2-prepare-input-order-file)
- [3. Prepare input OGID](#step-3-prepare-input-OGID-file)
- [4. Prepare input value](#step-4-prepare-input-value-files)
- [5. Run initialization clustering / prepare input config](#step-5-run-initialization-clutering-and-prepare-config-file)
- [6. Run Arboretum](#step-6-run-arboretum)
- [7. Run findTransitioningGenesets](#step-7-run-findTransitionGenesets)
- [8. Result visualization](#step-8-visualize-the-results)

---

### \[Step 0\] About this tutorial

This tutorial explains how to run the ***Arboretum*** and ***findTransitionGeneset***. Arboretum is a multi-task clustering framework which uses hierarchical relationship of samples simultaneously while grouping the genes into a finite number of expression states. In this tutorial, we will do the **clustering of genes** based on the **pseudo-bulk expression** of each dataset, i.e. average values of gene expressions within each dataset (cell cluster, in single cell setting). Once we have defined these **gene expression states** for all cell clusters, we can identify genes with interesting patterns of expression on the hierarchy and group the genes based on the similarity of their patterns to identify the **differential expressing (DE) genesets**. 

`input_files` folder contains tutorial dataset which consists of following files :
| Dataset | Description| File name | 
| :---    | :---  | :--- |
| c1_matrix | 155 cells x 20840 genes, with row & column headers | input_files/c1_matrix.txt |
| c2_matrix | 166 cells x 20840 genes, with row & column headers | input_files/c2_matrix.txt |
| c3_matrix | 126 cells x 20840 genes, with row & column headers | input_files/c3_matrix.txt |
| c4_matrix | 265 cells x 20840 genes, with row & column headers | input_files/c4_matrix.txt |
| c5_matrix | 155 cells x 20840 genes, with row & column headers | input_files/c5_matrix.txt |
| Gene ID | Gene names of 20840 genes | input_files/allgenenames.txt |
| Tree | [newick format](https://evolution.genetics.washington.edu/phylip/newicktree.html) tree text | input_files/newick_tree.txt 

- **Note**: the codes and scripts of this tutorial were fully tested and checked as working with the format above. So please firslty prepare your data as like above if you want to work with your own data.

`code` folder contains the main programs of this tutorial, which could be also downloaded from its own site:
| Program | Description | Original site |
| :---    | :---  | :--- |
| learnMoE | Initialization GMM clustering | https://github.com/Roy-lab/learnMoE |
| Arboretum | Multi-task gene clustering | https://github.com/Roy-lab/Arboretum2.0 |
| findTransitionGenesets | Hierarchical clustering of patterns identified from Arboretum | https://github.com/Roy-lab/clade-specific_gene_sets |
| reformatSpeciesTree | Format newick tree to Arboretum input | - |

- **Note**: All programs are written in C++ and the compling could be fulfilled by "Makefile" in each directory with the "make" command:


`script` folder contains several shell and PERL scripts which is used for the specific formatting of the data while running the wrapper scripts.

The key steps of clustering are performed by the following wrapper scripts:
| Wrapper script name | Running program | Output |
| :---    | :---  | :--- |
| step1_run_GMM_and_prep_config.sh | learnMoE | arb_input/merged_gmm_k#/, arb_input/config_k#.txt |
| step2_run_arboretum.sh | Arboretum | arb_output/ |
| step3_run_findTransitionGenesets.sh | findTransitionGenesets | transition_genesets_output/ |

---

### \[Step 1\] Prepare input tree

Arboretum requires a relationship tree structure between samples for the modeling. The format of the input tree is a text file consist of 3 columns (tab delimited). Each row is explaining the relationship of a parent node and a child node. One parental node always have 2 children nodes (left and right). The left/right children nodes will be dealt as equivalent, i.e. there is no order between children nodes. An ancestral node could be a child of another superordinate ancestral node. Following example is showing a tree consists of 3 nodes (c1, c2, c3) where c1 and c2 are more close. "Anc" means ancestral node.
```
c1 (TAB) left (TAB) Anc2
c2 (TAB) right (TAB) Anc2
Anc2 (TAB) left (TAB) Anc1
c3 (TAB) right (TAB) Anc1
```

To deal with a tree easier, especially for a complex tree with many leaves, we can use a program `reformatSpeciesTree` in the `code` folder. This program takes newick-style text file as a input and prints out the tree file for the Arboretum as the output. 

**Running command:**
```
- USAGE: reformatSpeciesTree [newick_tree.txt]
./code/reformatSpeciesTree input_files/newick_tree.txt > arb_input/tree.txt
```
`input_files/newick_tree.txt` shows the newick tree of relationship among 5 sampples in this tutorial. The output `arb_input/tree.txt` is a formatted tree for the Arboretum.

---

### \[Step 2\] Prepare input order file

This text files is a simple list of sample names WITHOUT Ancestral nodes ("Anc#") with format like below:
```
c2
c4
c3
c1
c5
```

**Running command:**
```
cut -f1 arb_input/input_tree.txt |grep -v Anc > arb_input/orders.txt
```
This script used the tree file prepared in the previous step. `arb_input/orders.txt` is the order file of our 5 samples example dataset, an input file of the Arboretum.

---

### \[Step 3\] Prepare input OGID file

OGID stands for "OrthoGroup ID" and it is a list of orthogroups(OGs) with a profiled list of corresponding gene IDs per species (wihch is a legacy of the original Arboretum). Regardless of the its naming, this file is used for the mapping of genes' identities from different datasets.
The format of this file is 
```
OGID (TAB) dataset1_GeneID,dataset2_GeneID,dataset3_GeneID,... 
```
(tab delimited between OGID and gene IDs, comma delimited among genes)

Note that the order of dataset-specific geneID should be same to the order file from [Step 2](#step-2-prepare-input-order-file). Also, from now on, the gene IDs of each dataset is named as "dataset_GeneID" (e.g. c1_AAA).

**Running command:**
```
- USAGE: generating_OGID.sh [orders.txt] [allgenenames.txt]
sh script/generating_OGID.sh arb_input/orders.txt input_files/allgenenames.txt > arb_input/OGID.txt
```
This script generates OGID based on the order of genes in "allgenenames.txt" file and shapes the gene IDs into "dataset#_geneID" as well as ordering it matched to the order file. The output `arb_input/OGID.txt` is an input file of the Arboretum.

---

### \[Step 4\] Prepare input value files

We are using pseudo-bulk expression vector of genes per samples, which is the average values of gene expressions within each dataset matrix so that each dataset has 1 vector values. The format of value file is a tab-delimited text file of 2 columns, the 1st column as "dataset_geneID" and the 2nd column is the mean expression values of all cells within the dataset. 
```
c1_A1BG	(TAB) 0.673579
c1_A1BG-AS1 (TAB) 0.246118
c1_A1CF (TAB) 0.000000
c1_A2M  (TAB) 0.899409
...
```

**Running command:**
```
- USAGE: generating_meanvals.pl [genes-by-cells_matrix.txt] ["dataset-name"]
perl script/generating_meanvals.pl input_files/c1_matrix.txt c1
mv c1_meanval.txt arb_input/
perl script/generating_meanvals.pl input_files/c2_matrix.txt c2
mv c2_meanval.txt arb_input/
... (Do tihs for all datset matrices.)
```
"arb_input/c#_meanval.txt" files are the values used as the source values for the Arboretum clustering.
- **Note**: The script expects the matrix text file as a tab-delimited [genes x cells] matrix WITH row and column headers. If the user's dataset has given as [cells x genes], one can transpose it by using "script/transpose_matrix.pl" like:
```
- USAGE: transpose_matrix.pl [(tab-delimited) matrix.txt]
perl script/transpose_matrix.pl c1_matrix.txt > c1_matrix_transposed.txt
```

---

### \[Step 5\] Run initialization clutering and prepare config file


**Running command:**
```
sh step1_run_GMM_and_prep_config.sh 3 arb_input/orders.txt 5
```


---

### \[Step 6\] Run Arboretum


**Running command:**
```
sh step2_run_arboretum.sh 3 arb_input/orders.txt arb_input/OGID.txt arb_input/input_tree.txt arb_input/config.txt c1 arb_output/
```

---

### \[Step 7\] Run findTransitionGenesets


**Running command:**
```
sh step3_run_findTransitionGenesets.sh arb_output/ arb_input/orders.txt arb_input/OGID.txt c1 transition_genesets_output

```

---

### \[Step 8\] Visualize the results


**Running command:**
```
sh script/draw_heatmap.sh transition_genesets_output/ordered_clusterset_means.txt 3 2
sh script/draw_heatmap.sh transition_genesets_output/clusterset101.txt 3 2
```
