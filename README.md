# Tutorial_for_Aboretum-based_DE_geneset_analysis
A tutorial which includes Arboretum &amp; findTransitionGenesets with an example dataset.

- [0. About this tutorial](#step-0-about-this-tutorial)
- [1. Prepare input tree](#step-1-prepare-input-tree)
- [2. Prepare input order](#step-2-prepare-input-order-file)
- [3. Prepare input value](#step-3-prepare-input-value-files)
- [4. Run TMF](#step-4-run-tmf)
- [5. Generate output clusters](#step-5-generate-output-clusters)

### \[Step 0\] About this tutorial

This tutorial explains how to run the ***Arboretum*** and ***findTransitionGeneset*** based on the result from Arboretum. Arboretum was originally developed for the gene exprssion analysis of multiple species, but the multi-task framework could be also used in other types of comparative analysis including comparative network anslsys or analysis of different cell groups.<br>
In this tutorial, we are assuming that we wanted to identify differential expressing geneset which exhibits complex patterns across multiple samples. Before going into the procedure, Please check the following initial files are correctly existing.

- [ ] `input_files` folder contains 5 different data matrices (`c1_datamat.txt`, ... , `c5_datamat.txt`) where the cells are different but the genes are common acorss all samples. Also, there is a file for tree relationship written as newick format (`newick.txt`).
- [ ] `code` folder contains the key programs of this procedure. Programs written in C++ are needed to be compiled via "make".
- [ ] `script` folder contains several shell and perl scripts which help out doing the formatting of the data.



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



### \[Step 2\] Prepare input order file

This text files is a simple list of sample names WITHOUT Ancestral nodes ("AncX"). This could be generated based on the tree file above by following command:
```
cut -f1 arb_input/input_tree.txt |grep -v Anc > arb_input/orders.txt
```
`arb_input/orders.txt` shows how the order file looks.



### \[Step 3\] Prepare input value files



### \[Step 4\] Run TMF


### \[Step 5\] Generate output clusters


See [Parameters](#parameters) for more details.
```
./run_tmf input/tree/toy_tree.txt 120 2 -o output/ -a 10 -l 200
```
- `input/toy_tree.txt` specifies the tree file, which contains file locations to individual task matrices (paths are relative to location of run_tmf executable location). 
- `120` is the number of features/columns in each task matrix, which has to be be the same across all tasks. 
- `2` = k, the smaller dimensions of U and V. 
-	[Optional] `-o output/` will put all output files to output/ directory. Check out the example output directory in the repo. By default output will be saved to current directory.
-	[Optional] `-a 10` will set the alpha (strength of regularization to parent node) to be 10. Default is alpha = 10.
- [Optional] `-l 200` will set lambda (strength of sparsity constraint) to be 200. By default there is no sparsity constraint, i.e., lambda = 0.

#### Input tree file format
See example in input/toy_tree.txt.
```
1 3 A input/toy/A.txt 95
2 3 B input/toy/B.txt 80
3 -1 root N/A N/A
```
- Column 1: **node ID**; start from 1 and move up.
- Column 2: **parent node ID**; right now the implementation will only work correctly if you ID all children nodes before a parent node (so start from the lowest depth of tree, move to next level, till you hit the root, which should be the last node ID.)
- Column 3: **node alias**, used as prefix to U and V output files.
- Column 4: **location of input matrix file for leaf nodes**, relative to where the `run_tmf` executable is. Set as N/A for non-leaf nodes.
- Column 5: **number of rows/data points in each input matrix**. Set as N/A for non-leaf nodes.
