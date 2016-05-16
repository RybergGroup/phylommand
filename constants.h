#ifndef PHYLOMMAND_CONSTANTS
#define PHYLOMMAND_CONSTANTS
static const unsigned int SIZE = 24;
static const char VERSION[] = "0.2.8";


namespace nexus_command {
    const char END = 'e';
    const char FIXAGE = 'f';
    const char MRCA = 'm';
    const char NON = 'n';
    const char TRANSLATE = 'r';
    const char TAXSET = 's';
    const char TREE = 't';
}

namespace nexus_block {
    const char TAXA = 'A';
    const char DATA = 'D';
    const char ERROR = 'E';
    const char NON = 'N';
    const char PHYLOMMAND = 'P';
    const char TREES = 'T';
    const char UNKNOWN = 'U';
}

#endif //PHYLOMMAND_CONSTANTS
