# Phylommand
Phylommand (**Phylo**genetics c**ommand** line) is a software package for **creating**, **manipulating**, and/or **getting statistics** from **trees** or working with **pairwise alignments**.
Phylommand is a commandline toolkit designed to be easily integrated in workflow pipelines.

# Table of Contents
- [The programs](#treebender)
  - [Treebender](#treebender)
  - [Treeator](#treeator)
  - [Contree](#contree)
  - [Pairalign](#pairalign)
  - [Rudisvg](#rudisvg)
- [Examples]
- [Download]
- [Compile/Build]
- [Command reference](#treebender-reference)
  - [Treebender reference](#treebender-reference)
  - [Treeator reference](#treeator-reference)
  - [Contree reference](#contree-reference)
  - [Pairalign reference](#pairalign-reference)
  - [Rudisvg reference](#rudisvg-reference)

## Treebender
Is the major program in phylommand to manipulate and get statistics from individual trees. You may get the branch length, get the sum of the branch lengths, set all branch lengths to a given value, multiply branch lengths in whole tree, a defined clade, or after a certain distance from the root, set short branches to zero, get depth of tree, distance to the root of all tips, get patristic distances, get the number of nodes with more than given support, get the log of the product of the support values (if support is probabilities = clade credibility), get statistics on internal node values (including values given in treeannotator/figtree format), clear internal nodes, get the numbers of the branches (used by some other functions), get the number of tips, get the tip names, change tip names, drop tips from the tree, test if a group is monophyletic, root on midpoint or using outgroup, get an outgroup that maximize the proportion of given tips in the outgroup compared to the tips that are not in the given set (or root on this group directly), ladderize tree, do nni branch swapping, generate random topologies, cluster tips based on tree,  split trees, get a subset of trees from a file, get matrix representation of a (set of) tree(s), and convert tree format between newick and nexus. 

## Treeator
Treeator is the main program to construct and evaluate trees in phylommand. It can construct neighbor-joining trees from distance matrices as produced by pairalign (see below). It can also calculate the parsimony score and do parsimony ancestral state reconstruction for a tree and datamatrix, as well as do parsimony stepwise addition to construct a tree. In addition it is possible to calculate likelihood scores and get the normalized likelihood for each trait at each internal node for one character and a given model, parameters, and tree. If compiled with NLOPT (see [installation instructions](Phylommand installation)), it is possible to optimize all or selected parameters of the substitution model. Treeator can read fasta, and sequential nexus and phylip formats but also require the alphabet for the character traits. The DNA, protein, and binary (0 1 -) alphabets are hardcoded into treeator, but other alphabets need to be given in a file (look in the example_files for an example of the DNA alphabet).
 
## Contree
Contree mainly compares trees, but it can also estimate decisiveness ([Sanderson et al. 2010](http://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-10-155). It is possible to compare all trees in a file against each other, or compare one set of trees against another. Topologies can be compared to calculate support values, conflicts between splits with more than given support can be determined, tips not present in the other tree can be identified, and [Robinson-Foulds metric](https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric) can be calculated. 

## Pairalign
Pairalign calculates statistics, and performs pairwise alignments from DNA sequences. It can calculate the proportion different sites (or similarity) and Jukes-Cantor distance, and the difference between the two. And it can output these sequential or as a matrix. It can calculate approximate MAD scores ([Smith et al. 2009](http://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-9-37)) to determine at which taxonomic level multiple sequences can be "easily" aligned, and perform clustering based on sequence distance.

## Rudisvg
Rudisvg is a rudimentary svg viewer that can be used to view trees output as svg by treebender (except for colors and different fonts). It is not meant to be a complete svg viewer but only a way to get a quick look at a tree.You should be able to move the tree around using the arrow keys, and zoom in and out with + and -.

<img src="https://github.com/mr-y/phylommand/wiki/rudisvg_tree.png" width="300" alt="Tree displayed by rudisvg">

# Examples

These examples are based on the EXAMPLE file distributed with phylommand. The examples assumes a
basic knowledge in how to work on the command line, for example that > will pipe
the output from stdout (the screen) to a file which is the standard way of
getting the output to a file in phylommand, e.g.:

    treebender --output svg tree_file.tree > tree_file.svg
 
Many cases also assume that you are working in a bash shell or similar.

To make and view a neighbour joining tree based on a mafft alignment:

    mafft alignment_file.fst | pairalign -A -j -n -m | treeator -n | treebender \
    --output svg | rudisvg

To create monophyletic OTUs based on the branch lengths in a tree (c.f. virtual
taxa sensu [Öpick et al 2009](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-8137.2009.02920.x/abstract)):

    treebender --cluster branch_length:0.03 tree_file.tree 

Get the parsimony score of all nearest neighbour interchange trees from:

    treebender --nni all tree_file.tree | treeator -f alignment_file.fst -p

Get alignable groups based on MAD score ([Smith et al. 2009](http://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-9-37)):

    pairalign --group alignment_groups alignment_file_with_taxon_string.fst

The MAD score does not care about gaps, this means that if the pairwise
alignments are over aligned (inserting more gaps than there should be to match
the same type of bases to each other), alignment_groups may over estimate what
can be aligned. Sometimes a multiple sequence alignment may be better in this
respect, and programs like, for example, muscle and mafft may even be faster
than pairalign in aligning. So it is possible to do something like:

    mafft alignment_file_with_taxon_string.fst | pairalign --group \
    alignment_groups -A

Delete old support values, and draw new one on a topology:

    treebender --clear_internal_node_labels tree_file.tree | contree -d \
    tree_file.trees.nex -s

Get conflict with supported above 70 based on an interval of 100 trees, compared
the support in a given tree, output in HTML format:

    treebender --interval 533-632 tree_file.trees.nex | contree -d - -a \
    tree_file.nex | contree -c 70 -d tree_file.tree --html

Get the maximum clade credibility tree in a set of trees, removing burn-in (5000
trees), and with help of awk:

    treebender --interval 5001 tree_file.trees | contree -a | treebender \
    --clade_credibility | awk -v max=0 '{if($1>max){print NR " " $1; max=$1}}'

The first value is which tree after burn-in, and the second the log sum of
support values. The last printed values should be for the maximum clade
credibility (MCC) tree. The second value is not the clade credibility since we
do not divide the support by the number of trees.

If we want to change the names of many or all tips in a tree it may be convenient 
to use a file with the name changes. We can get a nice start for such a file
getting the tip names from treebender:

    treebender -t '|,\n' tree_file.tree > name_changes.txt

We can then then edit name_changes.txt, putting our new names between the '|' and
',' (you need to add | for the last name though; depending on your system and
text editor you may want to use other line break than \n), and then run:

    treebender -c file:name_changes.txt tree_file.tree

Do not forget to add new name or delete the row if you do not want to change that
name, otherwise these tip names will be deleted. If you like the output in svg
format just add --output svg. Start files for the other arguments that take
taxa as input can be gotten similarly:

    treebender -t ,\\n tree_file.tree > taxa.txt

To be able to pipe all information for the --group analyses it is possible to
give the taxonomy in the sequence names. However, pairalign does not output the
taxonomy as part of the sequence name. This may be irritating sometimes when you
make a pairwise alignment to use for a later --group analysis. To get the
taxonomy from the first fasta file into a separate file that can be used for the
analysis you can do:

    grep '>' alignment_file_with_taxon_string.fst | sed 's/^>//' | awk \
    'BEGIN { FS="|" } {print $2 "|" $1}' > taxonomy_file.txt

If you want to summarize the rates in a MCC tree generated in with the BEAST
package you can do:

    treebender --internal_node_stats 1:rate_median tree_file.nex

It is also possible to do:

    treebender --internal_node_stats 1:rate tree_file.trees.nex

if you want to look at the rates over all trees in the tree distribution
produced by BEAST. This does however quickly become incomprehensible, so you
will probably need a script to do something meaningful with the output. Tracer
is probably also an better alternative, in most cases.

To reconstruct the ancestral state of the first trait in a DNA data matrix you
can do:

    treeator -a dna -t tree_file.tree --get_state_at_nodes \
    -l alignment_file.fst

If you want to make it a one rate model you can do:

    treeator -a dna -t tree_file.tree --get_state_at_nodes -l \
    -m 0,0,0,0,0,0,0,0,0,0,0,0 alignment_file.fst

or as an symmetric model:

    treeator -a dna -t tree_file.tree --get_state_at_nodes -l \
    -m 0,1,2,0,3,4,1,3,5,2,4,5 alignment_file.fst

Notice that these are not time reversible models. You may change the alphabet
given to -a if you are working with other type of characters.


# Treebender
Treebender is a command line program for manipulating trees. The program take a tree in newick or nexus format as indata through standard in or from a file and process it according to given options.

This manual is originally based on the help function (-h) of treebender but with successive changes.

## Usage
    treebender [arguments] < file.tree
    treebender [arguments] file.tree

For the second alternative you need to be careful so treebender does not interpret the filename as an extra argument to a switch. If this happen treebender will expect input from standard in and it will appear as nothing is happening. This can be avoided by giving the filename after the switch --file/-f. When taxa should be given as extra arguments they can be given in a text following the format for the argument. Newline and carriage returns will be ignored. The file name should be given behind the word file and colon, e.g.:

    treebender -d file:file_name.txt

## Arguments
### --branch_lengths / -a
Print branch lengths, the separator can be given as first argument after the switch, e.g.:

    treebender -a '\n'

Default is ','. If the switch r is given as second argument after a colon (:), e.g.:

    treebender -a '\n:r'

The value of the root branch will be printed as well, if n (default) is given it will not.  A separator between output from different trees can be given after another colon.
### --change_names / -c
Change the name of tips. Tip names to be changed should be given pairwise with present name first and new name second, separated by '|'. Separate pairs should be separated by ',' e.g.:

    treebender -c 'taxon1|new1,taxon2|new2'

The quotation marks are required to cancel the special meaning of |. If several tips have the same name all will be changed. Changes later in the list will be effected by changes made earlier in the list, e.g.

    treebender -c 'taxon1|new1,new1|new2'

will change the name of taxon1 to new2.

### --clade_credibility
Give the log of the product of the support of all clades.

### --clear_internal_node_labels
Delete the internal node labels

### --cluster
Get clusters based on method, e.g.:

    treebender --cluster branch_length

#### Available methods:
##### branch_length 
Separate clusters by single link clustering based on phylogenetic distance. Cut off should be given after colon, e.g.:

    treebender --cluster branch_length:0.03

##### long_branch
Returns taxa in clades on branches longer than cut off. Cut off should be given after colon (:).

##### tip_name
Cluster taxa based on taxon annotation. Should be followed after a colon by a single character that separates different parts of the tip name (default ' ') and after another colon (:) a number giving which position in the name should be used for clustering, (default 1), e.g.:

    treebender --cluster tip_name:\|:5

### --depth / -D
Get the longest distance from the root to any of the tips.

### --distances_to_root / -z
Output the number of nodes and branch length distance to the root for each tip. The separator between tip name and value can be specified, separated by colon, e.g.:

    treebender -p ",:|"

Default is newline and tab. A separator between output from different trees can be given after another colon

### --drop_tips / -d
Drop the given tips from the tree, e.g.

    treebender -d taxon1,taxon2,taxon3

### --get_tip_names / -t
Get the names of the tips in the tree, a separator can be specified, e.g.:

    treebender -t \\n

will give each name on separate rows (',' is the default separator).  A separator between output from different trees can be given after another colon

### --get_branch_numbers
Assign branch numbers as node labels.

### --get_relaxed_outgroup
Get the taxa in the clade that include the largest fraction of the difference between number of taxa included in the given group and number not included in the group divided by the total number in the group. Taxa given as comma separated string (see --drop_tips).

### --file / -f
Give file name, e.g.:

    treebender -f file.tree.

### --format
Give format of input, e.g.:

    treebender --format nexus

If no format is given and the input is a file treebender will try to guess the format, if it is through standard in it will assume newick format.

### --help / -h
Print this help.

### --internal_node_stats
Print stats about the values on the internal nodes. Counts nodes with value above given value, e.g.:

    treebender --internal_node_stats 1.96

Default: 1.0. If extra stats are given in FigTree/treeanotator format the parameter to summarize can be given behind colon, e.g.:

    treebender --internal_node_stats 1.96:rate

,or

    treebender --internal_node_stats :rate

### --interval
Only print the trees in the interval. Interval given as first tree to print - last tree to print, e.g.:

    treebender --interval 10-100

, or just the first tree to print, e.g.:

    treebender --interval 1000

### --inverse / -i
Inverse the string of taxa, e.g. drop all tips but the given. E.g:

    treebender -d taxon1,taxon2,taxon3 -i

### --is_monophyletic
Test if the given taxa form a monophyletic group, e.g.

    treebender --is_monophyletic taxon1,taxon2

### --ladderize / -l
Laddrize the tree. If followed by l - left ladderize, if followed by r - right ladderize (default), e.g.

    treebender -l r

### --matrix_representation
Present a fasta-formated matrix with splits of trees coded as characters. Intended for matrix representation parsimony.

### --mid_point_root / -m
Root the tree at the mid point.

### --multiply_branch_lengths / -u
Multiply each branch in the tree with the given value, e.g.

    treebender -u 3.5

Default 1.0.

### --multiply_branch_lengths_clade / -V
Multiply branches in clades defined by the most recent common ancestor of comma separated taxa. Separate clade with colon. E.g.

    treebender -V 3:Taxon_1,Taxon_2:Taxon_3,Taxon_4

### --multiply_branch_lengths_until / -U
Multiply branches in tree up until given distance (cut off) from root with the given value (separated by colon), e.g. 

    treebender -U 2:40 

Default 1.0:0.0.

### --n_supported
Get the number of nodes with higher support than given value. Should be followed by value, e.g.

    --n_supported 70.0

### --nni
Perform nearest neighbor interchange. If a integer is given as extra argument the interchange will be done on that branch (use --get_branch_numbers to get branch numbers). If 0 or no extra argument is given a branch will be selected randomly. If 'all' is given NNI will be performed for all branches, e.g.

    treebender --nni 4

, or

    treebender --nni all

### --no_branch_length / -0
Do not print branch lengths. If there are no branch lengths in input tree the default is to print zero length branches in the out tree. This argument override this and print no branch lengths.

### --null_short_branches
Set branches with shorter than given value to 0

### --number_of_taxa / -n
Get the number of tips/taxa in the tree.

### --outgroup_root / -o
Root using most recent common ancestor of given taxa, e.g.

    treebender -o taxa1,taxa2

### --output
Give tree format for output, nexus (nex or x for short), newick (new or w for short), or svg e.g.:

    treebender --output x

Default w. For svg extra graphical commands can be given after a colon (:). Each command should be separated by &, and commands and arguments should be separated by colon. Possible commands are: 'width' set width of figure; 'height' set hight of figure; 'offset' set offset between tips and tip label; 'strokewidth' set the width of the branches; 'fontsize' sets the size of the font used; 'font' sets which font to use; and 'tipcolor' sets the color of the tip labels given in parenthesis directly behind the color. 'width' and 'height' are mandatory if you want to set any other parameter than tip color. E.g.

    treebender --output 'svg:width:300&height:400&tipcolor:red(taxon1,taxon2,taxon3)green(taxon4)'

### --patristic_distances / -p
Get the total patristic distance to all other taxa in the tree for each taxon, the separator between different taxa, and the separator between taxon name and value can be specified (separated by colon) e.g.:

    treebender -p ",: | "

Default is new line and space. A separator between output from different trees can be given after another colon.

### --random_tree / -r
Get a random topology (no branch lengths) with given number of taxa, e.g. -r 20 (default 0). Number of random trees can be given behind a colon (:), e.g.:

    treebender -r 20:100

### --read_figtree_annotations
Will read annotations in FigTree/treeanotator format (e.g. [&rate=1.0,height=3.0]).

### --relaxed_outgroup_root
Will root on the group defined as for --get_relaxed_outgroup.

### --set_branch_lengths / -b
Set all branches in the tree to the given value, e.g.:

    treebender -b 0.5 (default 1.0)

### --split
Splits tree based on the longest branch (longest_branch/l) or the mid point (mid_point/m) until a stop criterion set by --split_stop is reached. Which derived tree to split in each iteration can be set after :. Either the tree with the longest branch (l; default for longest_branch split) or the tree with most tips (n; default for mid_point split).

### --split_stop
Sets criterion for when to stop splitting trees, either at a maximum number of trees (max_tree_number/t) or when all trees have fewer than a certain number of tips (max_tree_size/s). The number should be given togather with the specific criterion after :.

### --rooted / -R
Sets if the tree should be considered as rooted or not (only matters when splitting trees).

### --sum_branch_length / -s
Get the sum of the branch lengths in the tree (including root branch if length for this is given).

### --verbose / -v
Get additional output.


