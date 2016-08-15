# Makefile to compile phylcommand
# using gcc and g++
# If compiling on windows using MinGW add WIN=YES to make command

PP=g++ -std=c++11
CC=gcc

# Options
NLOPT= YES # if NLOPT set to NO set TREEATORFLAGS = -DNOLNLOPT
WIN = NO # add WIN=YES if compiling on windows using MinGW
PTHREADS = NO # set to yes to compile using pthreads, PTHREADS=-DPTREADS
DATABASE = YES # set to compile with database

PHYLOMMAND = treebender clustertree alignmentgroups pairalign contree treeator # anconstruction treesplitter superstat 
#CFLAGS
EXTRAS= # add EXTRAS=-DDEBUG for debug mode
TREEATORFLAGS = -DNLOPT
TREEATORLINKFLAGS = -lnlopt -lm
ifeq ($(NLOPT),NO)
    TREEATORFLAGS=
    TREEATORLINKFLAGS=
endif
SQLOFLAGS = -ldl -lpthread
ifeq ($(WIN),YES)
    SQLOFLAGS=
endif
TREE = tree.cpp
ALIGN = align_group.cpp
CLUST = clustertree.cpp
SEQPAIR = seqpair.cpp
DECISIVE = decisiveness.cpp
SQLITE = sqlite3.c
SQLITEO = sqlite3.o
TREEB = treebender.cpp
CLUSTTREE = clustertree_main.cpp
DATABASEFLAG = -DDATABASE
ifeq ($(DATABASE),NO)
    DATABASEFLAG =
    SQLITEO =
    SQLITEFLAGS =
endif
ALIGNMENT = alignmentgroups.cpp
ALIGNFLAGS =
ifeq ($(PTHREADS),YES)
    ALIGNFLAGS = -DPTHREAD
endif
PAIRALIGN = pairalign.cpp
# TREESPLIT = treesplitter.cpp
# SUPER = superstat.cpp
CONTREE = contree.cpp
TREEATOR = treeator.cpp
#ANCON = anconstruction.cpp
STRINGTREE = string_tree.cpp
NJTREE = nj_tree.cpp
SIMPLEML = simpleML.cpp
MARTH = marth/marth.cpp
#SUPPORTFUNCTIONS = support_functions.cpp
FILE_PARSER = file_parser.cpp
MATRIXPARS = matrix_parser.cpp
SEQDB = seqdatabase.cpp
INDEXEDFST = indexedfasta.cpp

OTREE = tree.o treebender.o string_tree.o file_parser.o matrix_parser.o # support_functions.o
OCLUSTTREE = clustertree.o tree.o clustertree_main.o string_tree.o matrix_parser.o file_parser.o $(SQLITEO) # support_functions.o
OALIGNMENT = seqpair.o align_group.o seqdatabase.o indexedfasta.o alignmentgroups.o $(SQLITEO)
OPAIRALIGN = seqpair.o pairalign.o
#OSPLIT = tree.o treesplitter.o string_tree.o matrix_parser.o # support_functions.o
#OSUPER = superstat.o tree.o decisiveness.o string_tree.o matrix_parser.o # support_functions.o
OCONTREE = contree.o tree.o decisiveness.o string_tree.o matrix_parser.o file_parser.o # support_functions.o
OTREEATOR = treeator.o tree.o string_tree.o nj_tree.o simpleML.o marth.o matrix_parser.o file_parser.o $(TREEATORLINKFLAGS) # support_functions.o

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

#treesplitter: $(OSPLIT)
#	$(PP) -o treesplitter $(OSPLIT)

#superstat: $(OSUPER)
#	$(PP) -o superstat $(OSUPER)

contree: $(OCONTREE)
	$(PP) -o contree $(OCONTREE)

treeator: $(OTREEATOR)
	$(PP) -o treeator $(OTREEATOR) $(TREEATORFLAGS)

treebender.o: $(TREEB)
	$(PP) -c $(TREEB) $(EXTRAS)

clustertree_main.o: $(CLUSTTREE)
	$(PP) $(DATABASEFLAG) -c $(CLUSTTREE) $(EXTRAS)

alignmentgroups.o: $(ALIGNMENT)
	$(PP) $(DATABASEFLAG) $(ALIGNFLAGS) -c $(ALIGNMENT) $(EXTRAS)
#	$(PP) $(DATABASEFLAG) -c $(ALIGNMENT) $(EXTRAS)

pairalign.o: $(PAIRALIGN)
	$(PP) -c $(PAIRALIGN) $(EXTRAS)

#treesplitter.o: $(TREESPLIT)
#	$(PP) -c $(TREESPLIT) $(EXTRAS)

#superstat.o: $(SUPER)
#	$(PP) -c $(SUPER) $(EXTRAS)

contree.o: $(CONTREE)
	$(PP) -c $(CONTREE) $(EXTRAS)

treeator.o: $(TREEATOR)
	$(PP) -c $(TREEATOR) $(EXTRAS)

#anconstruction.o: $(ANCON)
#	$(PP) -c $(ANCON) $(EXTRAS)

tree.o: $(TREE)
	$(PP) -c $(TREE) $(EXTRAS)

align_group.o: $(ALIGN)
	$(PP) -c $(ALIGN) $(EXTRAS)

clustertree.o: $(CLUST)
	$(PP) $(DATABASEFLAG) -c $(CLUST) $(EXTRAS)

seqpair.o: $(SEQPAIR)
	$(PP) -c $(SEQPAIR) $(EXTRAS)

decisiveness.o: $(DECISIVE)
	$(PP) -c $(DECISIVE) $(EXTRAS)

ifeq ($(DATABASE), YES)
sqlite3.o: $(SQLITE)
	$(CC) -c $(SQLITE)
endif

string_tree.o: $(STRINGTREE)
	$(PP) -c $(STRINGTREE) $(EXTRAS)

nj_tree.o: $(NJTREE)
	$(PP) -c $(NJTREE) $(EXTRAS)

simpleML.o: $(SIMPLEML)
	$(PP) -c $(SIMPLEML) $(EXTRAS)

marth.o: $(MARTH)
	$(PP) -c $(MARTH) $(EXTRAS)

#support_functions.o: $(SUPPORTFUNCTIONS)
#	$(PP) -c $(SUPPORTFUNCTIONS)

file_parser.o: $(FILE_PARSER)
	$(PP) -c $(FILE_PARSER) $(EXTRAS)

matrix_parser.o: $(MATRIXPARS)
	$(PP) -c $(MATRIXPARS) $(EXTRAS)

seqdatabase.o: $(SEQDB)
	$(PP) $(DATABASEFLAG) -c $(SEQDB) $(EXTRAS)

indexedfasta.o: $(INDEXEDFST)
	$(PP) -c $(INDEXEDFST) $(EXTRAS)
