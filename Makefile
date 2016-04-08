# Makefile to compile phylcommand
# using gcc and g++

PP=g++
CC=gcc
#CFLAGS
PHYLOMMAND = treebender clustertree alignmentgroups pairalign treespliter superstat conftree treeator
TREE = tree.cpp
ALIGN = align_group.cpp
CLUST = clustertree.cpp
SEQPAIR = seqpair.cpp
DECISIVE = decisiveness.cpp
SQLITE = sqlite3.c
TREEB = treebender.cpp
CLUSTTREE = clustertree_main.cpp
CLUSTFLAG = -DDATABASE
ALIGNMENT = alignmentgroups.cpp
ALIGNFLAGS = -DPTHREAD
PAIRALIGN = pairalign.cpp
TREESPLIT = treespliter.cpp
SUPER = superstat.cpp
CONFTREE = conftree.cpp
TREEATOR = treeator.cpp
STRINGTREE = string_tree.cpp
NJTREE = nj_tree.cpp
#SUPPORTFUNCTIONS = support_functions.cpp
FILE_PARSER = file_parser.cpp

OTREE = tree.o treebender.o string_tree.o file_parser.o # support_functions.o
OCLUSTTREE = clustertree.o tree.o sqlite3.o clustertree_main.o string_tree.o # support_functions.o
OALIGNMENT = seqpair.o align_group.o sqlite3.o alignmentgroups.o
OPAIRALIGN = seqpair.o pairalign.o
OSPLIT = tree.o treespliter.o string_tree.o # support_functions.o
OSUPER = superstat.o tree.o decisiveness.o string_tree.o # support_functions.o
OCONFTREE = conftree.o tree.o string_tree.o # support_functions.o
OTREEATOR = treeator.o tree.o string_tree.o nj_tree.o # support_functions.o
SQLOFLAGS = -ldl -lpthread

# treeator.cpp tree.cpp string_tree.cpp nj_tree.cpp

all: $(PHYLOMMAND)

treebender: $(OTREE)
	$(PP) -o treebender $(OTREE)

clustertree: $(OCLUSTTREE)
	$(PP) -o clustertree $(OCLUSTTREE) $(SQLOFLAGS)

alignmentgroups: $(OALIGNMENT)
	$(PP) -o alignmentgroups $(OALIGNMENT) $(SQLOFLAGS)

pairalign: $(OPAIRALIGN)
	$(PP) -o pairalign $(OPAIRALIGN)

treespliter: $(OSPLIT)
	$(PP) -o treespliter $(OSPLIT)

superstat: $(OSUPER)
	$(PP) -o superstat $(OSUPER)

conftree: $(OCONFTREE)
	$(PP) -o conftree $(OCONFTREE)

treeator: $(OTREEATOR)
	$(PP) -o treeator $(OTREEATOR)

treebender.o: $(TREEB)
	$(PP) -c $(TREEB)

clustertree_main.o: $(CLUSTTREE)
	$(PP) $(CLUSTFLAG) -c $(CLUSTTREE)

alignmentgroups.o: $(ALIGNMENT)
#	$(PP) $(ALIGNFLAGS) -c $(ALIGNMENT)
	$(PP) -c $(ALIGNMENT)
pairalign.o: $(PAIRALIGN)
	$(PP) -c $(PAIRALIGN)

treespliter.o: $(TREESPLIT)
	$(PP) -c $(TREESPLIT)

superstat.o: $(SUPER)
	$(PP) -c $(SUPER)

conftree.o: $(CONFTREE)
	$(PP) -c $(CONFTREE)

treeator.o: $(TREEATOR)
	$(PP) -c $(TREEATOR)

tree.o: $(TREE)
	$(PP) -c $(TREE)

align_group.o: $(ALIGN)
	$(PP) -c $(ALIGN)

clustertree.o: $(CLUST)
	$(PP) $(CLUSTFLAG) -c $(CLUST)

seqpair.o: $(SEQPAIR)
	$(PP) -c $(SEQPAIR)

decisiveness.o: $(DECISIVE)
	$(PP) -c $(DECISIVE)

sqlite3.o: $(SQLITE)
	$(CC) -c $(SQLITE)

string_tree.o: $(STRINGTREE)
	$(PP) -c $(STRINGTREE)

nj_tree.o: $(NJTREE)
	$(PP) -c $(NJTREE)

#support_functions.o: $(SUPPORTFUNCTIONS)
#	$(PP) -c $(SUPPORTFUNCTIONS)

file_parser.o: $(FILE_PARSER)
	$(PP) -c $(FILE_PARSER)
