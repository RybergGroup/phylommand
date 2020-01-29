# Makefile to compile phylommand
# using gcc and g++
# If compiling on windows using MinGW add WIN=YES to make command

PP=g++ -std=c++11
CC=gcc

# Options
NLOPT= NO # if NLOPT set to NO set TREEATORFLAGS = -DNOLNLOPT
WIN = NO # add WIN=YES if compiling on windows using MinGW
PTHREADS = NO # set to YES to compile using pthreads, PTHREADS=-DPTREADS
DATABASE = NO # set to YES compile with database
RUDISVG = NO # set to  YES to compile rudisvg

ifeq ($(RUDISVG),YES)
PHYLOMMAND = treebender pairalign contree treeator rudisvg
else
PHYLOMMAND = treebender pairalign contree treeator
endif
#CFLAGS
EXTRAS= # add EXTRAS=-DDEBUG for debug mode
TREEATORFLAGS=
TREEATORLINKFLAGS=
RUDISVGLINKFLAGS = -I/usr/X11R6/include -I/usr/X11R6/include/X11 -L/usr/X11R6/lib -L/usr/X11R6/lib/X11 -lX11 #-lX11
ifeq ($(NLOPT),YES)
    TREEATORFLAGS = -DNLOPT
    TREEATORLINKFLAGS = -lnlopt -lm
endif
TREE = tree.cpp
ALIGN = align_group.cpp
CLUST = clustertree.cpp
SEQPAIR = seqpair.cpp
DECISIVE = decisiveness.cpp
SQLITE = sqlite3.c
SQLITEO =
TREEB = treebender.cpp
DATABASEFLAG =
SQLITEFLAGS =
PAIRALIGNFLAGS =
ifeq ($(DATABASE),YES)
    SQLITEO = sqlite3.o
    DATABASEFLAG = -DDATABASE
    SQLOFLAGS = -ldl -lpthread
    PAIRALIGNFLAGS = -ldl -lpthread
endif
ALIGNFLAGS =
ifeq ($(PTHREADS),YES)
    ALIGNFLAGS = -DPTHREAD
    PAIRALIGNFLAGS = -ldl -lpthread
endif
ifeq ($(WIN),YES)
    SQLOFLAGS=
endif
PAIRALIGN = pairalign.cpp
CONTREE = contree.cpp
TREEATOR = treeator.cpp
RUDISVGCPP = rudisvg.cpp
PUGIXML = pugixml/pugixml.cpp
STRINGTREE = string_tree.cpp
NJTREE = nj_tree.cpp
SIMPLEML = simpleML.cpp
MARTH = marth/marth.cpp
FILE_PARSER = file_parser.cpp
ARGV_PARSER = argv_parser.cpp
MATRIXPARS = matrix_parser.cpp
SEQDB = seqdatabase.cpp
INDEXEDFST = indexedfasta.cpp
SUB_MODEL = sub_model.cpp
RATE_MODEL = rate_model.cpp
OTREE = tree.o treebender.o string_tree.o file_parser.o argv_parser.o matrix_parser.o clustertree.o rate_model.o $(SQLITEO)
OPAIRALIGN = seqpair.o pairalign.o align_group.o seqdatabase.o argv_parser.o indexedfasta.o $(SQLITEO)
OCONTREE = contree.o tree.o decisiveness.o string_tree.o matrix_parser.o file_parser.o rate_model.o
OTREEATOR = treeator.o tree.o string_tree.o nj_tree.o simpleML.o marth.o matrix_parser.o file_parser.o argv_parser.o sub_model.o rate_model.o
ifeq ($(RUDISVG),YES)
ORUDISVG = rudisvg.o pugixml.o
endif

all: $(PHYLOMMAND)

treebender: $(OTREE)
	$(PP) -o treebender $(OTREE) $(SQLOFLAGS)

pairalign: $(OPAIRALIGN)
	$(PP) -o pairalign $(OPAIRALIGN) $(PAIRALIGNFLAGS)

contree: $(OCONTREE)
	$(PP) -o contree $(OCONTREE)

treeator: $(OTREEATOR)
	$(PP) -o treeator $(OTREEATOR) $(TREEATORLINKFLAGS)

ifeq ($(RUDISVG),YES)
rudisvg: $(ORUDISVG)
	$(PP) -o rudisvg $(ORUDISVG) $(RUDISVGLINKFLAGS)
endif

treebender.o: $(TREEB)
	$(PP) $(DATABASEFLAG) -c $(TREEB) $(EXTRAS)

pairalign.o: $(PAIRALIGN)
	$(PP) $(ALIGNFLAGS) $(DATABASEFLAG) -c $(PAIRALIGN) $(EXTRAS)

contree.o: $(CONTREE)
	$(PP) -c $(CONTREE) $(EXTRAS)

treeator.o: $(TREEATOR)
	$(PP) -c $(TREEATOR) $(TREEATORFLAGS) $(EXTRAS)

ifeq ($(RUDISVG),YES)
rudisvg.o: $(RUDISVGCPP)
	$(PP) -c $(RUDISVGCPP) $(EXTRAS)
endif

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

sub_model.o: $(SUB_MODEL)
	$(PP) -c $(SUB_MODEL) $(EXTRAS)

rate_model.o: $(RATE_MODEL)
	$(PP) -c $(RATE_MODEL) $(EXTRAS)

ifeq ($(DATABASE),YES)
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

file_parser.o: $(FILE_PARSER)
	$(PP) -c $(FILE_PARSER) $(EXTRAS)

argv_parser.o: $(ARGV_PARSER)
	$(PP) -c $(ARGV_PARSER) $(EXTRAS)

matrix_parser.o: $(MATRIXPARS)
	$(PP) -c $(MATRIXPARS) $(EXTRAS)

seqdatabase.o: $(SEQDB)
	$(PP) $(DATABASEFLAG) -c $(SEQDB) $(EXTRAS)

indexedfasta.o: $(INDEXEDFST)
	$(PP) -c $(INDEXEDFST) $(EXTRAS)

ifeq ($(RUDISVG),YES)
pugixml.o: $(PUGIXML)
	$(PP) -c $(PUGIXML) $(EXTRAS)
endif

clean:
	$(RM) $(PHYLOMMAND) *.o
