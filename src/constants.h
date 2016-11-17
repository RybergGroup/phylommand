/********************************************************************
Copyright (C) 2016 Martin Ryberg

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

contact: martin.ryberg@ebc.uu.se
*********************************************************************/

#ifndef PHYLOMMAND_CONSTANTS
#define PHYLOMMAND_CONSTANTS
#include <map>
#include <bitset>

static const unsigned int SIZE = 24;
static const char VERSION[] = "1.0";
static const char YEAR[] = "2016";

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
