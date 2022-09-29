# Tutorial_for_Aboretum-based_DE_geneset_analysis
A tutorial which includes Arboretum &amp; findTransitionGenesets with an example dataset.

- [0. About this tutorial](#step-0-about-this-tutorial)
- [1. Prepare input tree](#step-1-prepare-input-tree)
- [2. Prepare input order](#step-2-prepare-input-order-file)
- [3. Prepare input value](#step-3-prepare-input-value-files)
- [4. Initialization clustering / prepare input config](#step-4-do-initialization-clutering-and-prepare-config-file)
- [5. Run Arboretum](#step-5-run-arboretum)
- [6. Run findTransitioningGenesets](#step-6-run-findTransitionGenesets)
- [7. Result visualization](#step-7-visualize-the-results)

---

### \[Step 0\] About this tutorial

This tutorial explains how to run the ***Arboretum*** and ***findTransitionGeneset*** based on the result from Arboretum. Arboretum was originally developed for the gene exprssion analysis of multiple species, but the multi-task framework could be also used in other types of comparative analysis including comparative network anslsys or analysis of different cell groups.<br>
In this tutorial, we are assuming that we wanted to identify differential expressing geneset which exhibits complex patterns across multiple samples. Before going into the procedure, Please check the following initial files are correctly existing.

- [ ] `input_files` folder contains 5 different data matrices (`c1_datamat.txt`, ... , `c5_datamat.txt`) where the cells are different but the genes are common acorss all samples. Also, there is a file for tree relationship written as newick format (`newick.txt`).
- [ ] `code` folder contains the key programs of this procedure. Programs written in C++ are needed to be compiled via "make".
- [ ] `script` folder contains several shell and perl scripts which help out doing the formatting of the data.

---

### \[Step 1\] Prepare input tree

Arboretum requires a relationship tree structure between samples for the modeling. The format of the input tree is a text file consist of 3 columns (tab delimited). Each row is explaining the relationship of a parent node and a child node. One parental node always have 2 children nodes (left and right). The left/right children nodes will be dealt as equivalent, i.e. there is no order between children nodes. An ancestral node could be a child of another superordinate ancestral node. Following example is showing a tree consists of 3 nodes (c1, c2, c3) where c1 and c2 are more close. "Anc" means ancestral node.
```
c1 (TAB) left (TAB) Anc2
c2 (TAB) right (TAB) Anc2
Anc2 (TAB) left (TAB) Anc1
c3 (TAB) right (TAB) Anc1
```
To deal with a tree easier, especially for a complex tree with many leaves, we can use a program `reformatSpeciesTree` in the `code` folder. This program takes newick-style text file as a input and prints out the tree file for the Arboretum as the output. Usage of the program is like below:
```
./code/reformatSpeciesTree input_files/newick_tree.txt > arb_input/tree.txt
```
`input_files/newick_tree.txt` exhibits the newick tree of relationship among 5 sampples in this tutorial. The final input format of tree for Arboretum is `arb_input/tree.txt`.

---

### \[Step 2\] Prepare input order file

This text files is a simple list of sample names WITHOUT Ancestral nodes ("AncX"). This could be generated based on the tree file above by following command:
```
cut -f1 arb_input/input_tree.txt |grep -v Anc > arb_input/orders.txt
```
`arb_input/orders.txt` shows how the order file looks.

---

### \[Step 3\] Prepare input value files

```
perl script/generating_meanvals.pl input_files/c1_matrix.txt c1
mv c1_meanval.txt arb_input/
perl script/generating_meanvals.pl input_files/c2_matrix.txt c2
mv c2_meanval.txt arb_input/
...
```

---

### \[Step 4\] Do initialization clutering and prepare config file

```
sh step1_run_GMM_and_prep_config.sh 3 arb_input/orders.txt 5
```


---

### \[Step 5\] Run Arboretum

```
sh step2_run_arboretum.sh 3 arb_input/orders.txt arb_input/OGID.txt arb_input/input_tree.txt arb_input/config.txt c1 arb_output/
```

---

### \[Step 6\] Run findTransitionGeneset

```
sh step3_run_findTransitionGenesets.sh arb_output/ arb_input/orders.txt arb_input/OGID.txt c1 transition_genesets_output

```

---

### \[Step 7\] Visualize the results

```
sh script/draw_heatmap.sh transition_genesets_output/ordered_clusterset_means.txt 3 2
sh script/draw_heatmap.sh transition_genesets_output/clusterset101.txt 3 2
```
