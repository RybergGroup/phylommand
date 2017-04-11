# Phylommand
Phylommand (**Phylo**genetics c**ommand** line) is a software package for **creating**, **manipulating**, and/or **getting statistics** from **trees** or working with **pairwise alignments**.
Phylommand is a commandline toolkit designed to be easily integrated in workflow pipelines.

## Table of Contents
- [The programs](#treebender)
  - [Treebender](#treebender)
  - [Treeator](#treeator)
  - [Contree](#contree)
  - [Pairalign](#pairalign)
  - [Rudisvg](#rudisvg)
- [Download](#download)
- [Compile/Build](#compilebuild)
- [Examples](#examples)
- [Miscellaneous](#miscellaneous)


# The programs
All programs accept input from standard in, or a file. This means that input can
be piped from one (or another) program to the other. The program that deals with
trees accept trees in newick or nexus format; treeator also accept character
data matrices in fasta, sequential phylip, and sequential nexus format, and space
separated left triangular distance matrices; for decisiveness estimates contree
require a text file with the genes for each taxon on a separate row. Output is
printed to stdout (the screen) and can be piped to a file using '>' in most
command line environments.

## Treebender
Is the major program in phylommand to manipulate and get statistics from individual trees. You may get the branch length, get the sum of the branch lengths, set all branch lengths to a given value, multiply branch lengths in whole tree, a defined clade, or after a certain distance from the root, set short branches to zero, get depth of tree, distance to the root of all tips, get patristic distances, get the number of nodes with more than given support, get the log of the product of the support values (if support is probabilities = clade credibility), get statistics on internal node values (including values given in treeannotator/figtree format), clear internal nodes, get the numbers of the branches (used by some other functions), get the number of tips, get the tip names, change tip names, drop tips from the tree, test if a group is monophyletic, root on midpoint or using outgroup, get an outgroup that maximize the proportion of given tips in the outgroup compared to the tips that are not in the given set (or root on this group directly), ladderize tree, do nni branch swapping, generate random topologies, cluster tips based on tree,  split trees, get a subset of trees from a file, get matrix representation of a (set of) tree(s), and convert tree format between newick and nexus. 

[Complete command line reference for Treebender in Wiki](https://github.com/mr-y/phylommand/wiki/treebender)

## Treeator
Treeator is the main program to construct and evaluate trees in phylommand. It can construct neighbor-joining trees from distance matrices as produced by pairalign (see below). It can also calculate the parsimony score and do parsimony ancestral state reconstruction for a tree and datamatrix, as well as do parsimony stepwise addition to construct a tree. In addition it is possible to calculate likelihood scores and get the normalized likelihood for each trait at each internal node for one character and a given model, parameters, and tree. If compiled with NLOPT (see [installation instructions](Phylommand installation)), it is possible to optimize all or selected parameters of the substitution model. Treeator can read fasta, and sequential nexus and phylip formats but also require the alphabet for the character traits. The DNA, protein, and binary (0 1 -) alphabets are hardcoded into treeator, but other alphabets need to be given in a file (look in the example_files for an example of the DNA alphabet).
 
[Complete command line reference for Treeator in Wiki](https://github.com/mr-y/phylommand/wiki/treeator)
 
## Contree
Contree mainly compares trees, but it can also estimate decisiveness ([Sanderson et al. 2010](http://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-10-155). It is possible to compare all trees in a file against each other, or compare one set of trees against another. Topologies can be compared to calculate support values, conflicts between splits with more than given support can be determined, tips not present in the other tree can be identified, and [Robinson-Foulds metric](https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric) can be calculated. 

[Complete command line reference for Contree in Wiki](https://github.com/mr-y/phylommand/wiki/contree)

## Pairalign
Pairalign calculates statistics, and performs pairwise alignments from DNA sequences. It can calculate the proportion different sites (or similarity) and Jukes-Cantor distance, and the difference between the two. And it can output these sequential or as a matrix. It can calculate approximate MAD scores ([Smith et al. 2009](http://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-9-37)) to determine at which taxonomic level multiple sequences can be "easily" aligned, and perform clustering based on sequence distance.

[Complete command line reference for Pairalign in Wiki](https://github.com/mr-y/phylommand/wiki/pairalign)

## Rudisvg
Rudisvg is a rudimentary svg viewer that can be used to view trees output as svg by treebender (except for colors and different fonts). It is not meant to be a complete svg viewer but only a way to get a quick look at a tree.You should be able to move the tree around using the arrow keys, and zoom in and out with + and -.

<img src="https://github.com/mr-y/phylommand/wiki/rudisvg_tree.png" width="300" alt="Tree displayed by rudisvg">

# Download
There are pre-compiled binaries for the following OS:
- [Windows](https://github.com/mr-y/phylommand/tree/master/bin)
- [OsX](https://github.com/mr-y/phylommand/tree/master/bin)

(For Linux and if you have problems with the pre-compiled files, see [Compile/Build](#compilebuild) below)

To be able to execute the phylommand programs from any folder, download and move the programs 
to a folder in your PATH, e.g.

```bash
# Check what dirs are in your PATH
echo $PATH

# Usually /usr/local/bin/ is there and is a suitable dir
# otherwise choose another, then:

sudo mv treebender /usr/local/bin/
sudo mv treeator /usr/local/bin/
sudo mv pairalign /usr/local/bin/
sudo mv contree /usr/local/bin/
```
**Note:** You might need to install XQuarts(https://www.xquartz.org/) on some osX distributions for rudisvg to work.

on most UNIX like systems (e.g. LINUX and OS X), or put the folder with the programs into your PATH.

# Compile/Build
Phylommand is written in C++ and can be compiled using make and the gnu compilers. The basic version with the four core programs can be compiled executing make on the command line in the folder with the phylommand code (you get the files either by cloning the github repo or downloading it as a zip-file):

    make

Phylommand has been successfully compiled on linux (Ubuntu) and MacOS X, and also on Windows using MinGW, but then with the WIN=YES option:

    make WIN=YES

On Bash on Ubuntu on Windows (available for Windows 10 anniversary update), phylommand has been installed just as on linux (i.e. without the WIN=YES option).
There are three major additional options for phylommand.

PTHREADS=YES enables the use of multi-threading when working with pairwise alignments (not compatible with WIN=YES).

NLOPT=YES require that the library [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt) is installed. On Ubuntu 16.04 this is available as a package an can be installed from on the command line by:

     sudo apt-get install libnlopt-dev

This option is needed to optimize model parameters using Maximum Likelihood, if you do not want to do that using phylommand you do not need this option.

RUDISVG=YES will compile a rudimentary SVG viewer (rudisvg) as well. This option require the x11 developmental libraries. If rudisvg is installed on bash on Ubuntu on Windows you will need a X server as well, for example [Xming](https://sourceforge.net/projects/xming/), and set the DISPLAY variable:

    export DISPLAY=localhost:0.0

There is no need to do this on OS X or ordinary linux distributions (like Ubuntu).


# Examples
These examples are based on the [example_files/](example_files/) distributed with phylommand. The examples assumes a
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

#Miscellaneous
Phylommand work on fully resolved (bifurcating) trees and will arbitrarily
resolve any polytomies. This needs to be especially considered when defining
or queering monophyletic groups.

Phylommand does not handle : ; , ( ) [ ] and white space in tip names very well.
So avoid these characters, even if the names are surrounded by ' or ". It handles
white space better in nexus format than in newick format, using the translate
NEXUS command.

Treeator uses the NLOPT LN_NELDERMEAD algorithm for optimization (Nelder JA &
Mead R. 1965. A simplex method for function minimization. The Computer J. 7:308
-313).

Treeator only do likelihood on one character at the time. If given a multi
character matrix it will calculate the likelihood of the first character. This
may change in future editions, so do not count on this behaviour. 

Treeator calculate likelihoods even if there are negative branches in the tree.
So look out for warnings about negative branches.

The neighbour-joining implementation will produce the optimal solution, even if
this include negative branches.

When colon (:) is used to separate options in extra arguments, back slash (\) can
be used as an escape character or to help denote newline (\n/\r) and tab (\t).
The special meaning of the back slash can be escaped by a back slash.

It is possible to run rudisvg on ubuntu on windows using the Xming X11 server. You
then need to set the DISPLAY parameter in your bash window by:

    export DISPLAY=:0

This option is not officially supported by Windows and may not work on all
machines.

The MAD score calculated in pairalign is only approximate, the precision needs
to be set at compile time. To change the precision change the value of
precision at line 33 in align_group.h (default 10000 i.e. 0.001). If clustering
the calculation of MAD only include pairs that have a lower similarity than the
cut off. If all your sequences are clustered together, the alignment group will
be empty (this is also the case if no taxonomy was read).

The calculation of decisiveness is based on randomly generated topologies. And is
thus only an estimation. To increase the precision increase the number of
iterations.

Empty trees are printed as -1;. You may get this, for example, if you do nni on
branch 1, and the sister lineage only have one tip (i.e. is a terminal branch).

If treebender gets a tip name that it does not recognize it usually just ignore
it. So if a taxon to recognize the most recent common ancestor is misspelt it is
just ignored (show must go on).

When drawing printing trees in SVG (and HTML) and there is no branch length, all
branches will be printed with the length 0. To change this you can use the
--set_branch_lengths option in treebender, for example setting all branch lengths
to 1.



