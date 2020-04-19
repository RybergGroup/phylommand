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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <ctime>
#include <vector>
#include "pugixml/pugixml.hpp"
#include "drawing.h"
#include "app_display.h"

using namespace std;

/*void sleep(unsigned int mseconds) {
    clock_t goal = mseconds + clock();
    while (goal > clock());
}*/

/*class app_display {
    public:
	app_display() {
	    // Prepair window
	    #ifdef DEBUG
	    cerr << "Opening display" << endl;
	    #endif // DEBUG
	    dis = XOpenDisplay(NULL);
	    screen = DefaultScreen(dis);
	    black_int = BlackPixel(dis, screen); 
	    white_int = WhitePixel(dis, screen);
	    colormap = DefaultColormap(dis, 0);
	    scale = 1.0;
	    font = BadName;
	}
	~app_display() {
	    XUnloadFont(dis,font);
	    XFreeGC(dis, gc);
	    XDestroyWindow(dis, win);
	    XCloseDisplay(dis);
	}

	void start_win (const long colour) {
	    win = XCreateSimpleWindow(dis,RootWindow(dis,0),1,1,500,500, 0, colour, colour);

	    XSetStandardProperties(dis,win, "rudisvg", "i!", None, NULL, 0, NULL);

	    XSelectInput(dis, win, ExposureMask | KeyPressMask | ButtonPressMask);

	    // To close window without error
	    Atom delWindow = XInternAtom( dis, "WM_DELETE_WINDOW", 0 );
	    XSetWMProtocols(dis , win, &delWindow, 1);
	    #ifdef DEBUG
	    cerr << "Maping window" << endl;
	    #endif // DEBUG
	    XMapWindow(dis, win);
	    XFlush(dis);
	    //sleep(5000);
	    #ifdef DEBUG
	    cerr << "Setting windows options" << endl;
	    #endif // DEBUG
	    gc = XCreateGC(dis, win, 0, 0);
	}
	void start_win( ) { start_win(white_int); }
	void start_white_win( ) { start_win(white_int); }
	void start_black_win( ) { start_win(black_int); }
	void set_default_col( const char hex_col[8] ) {
	    set_FG_col(hex_to_int_col(hex_col));
	}
	void set_default_FG ( ) {
	    set_FG_col( default_col);
	}
	void set_FG_col ( unsigned int col ) {
	    XSetForeground(dis, gc, col);
	    present_col = col;
	}
	unsigned long hex_to_int_col(const char hex_col[8]) {
	    XColor FGcol;
	    XParseColor(dis, colormap, hex_col, &FGcol);
	    XAllocColor(dis, colormap, &FGcol);
	    return FGcol.pixel;
	}
	bool is_default_col( unsigned long query ) { return query==default_col; }
	void draw_line ( const int x1, const int y1, const int x2, const int y2 ) {
	    XDrawLine(dis, win, gc, (x1-x_setoff)*scale, (y1-y_setoff)*scale, (x2-x_setoff)*scale, (y2-y_setoff)*scale);
	}
	void draw_rectangle ( const int x, const int y, const int width, const int height) {
	    XDrawRectangle(dis, win, gc, (x-x_setoff)*scale, (y-y_setoff)*scale, width*scale, height*scale);
	}
	void draw_arc (const int x, const int y, const unsigned int width, const unsigned int height, const int angle1, const int angle2) {
	    XDrawArc(dis, win, gc, (x-x_setoff)*scale, (y-y_setoff)*scale, width*scale, height*scale, angle1, angle2);
	}
	void draw_string ( string& text, const int x, const int y, const unsigned int text_size, const unsigned long col) {
	    unsigned int font_size(0);
	    unsigned int size(0);
	    if (size == 0 || text_size != font_size) {
		font_size = text_size;
		if (font_size >0) size = font_size;
		else size = 10;
		size *= scale;
		bool increase(false);
		bool end(false);
		int Nfonts(0);
		while (font == BadName || font == BadAlloc || !Nfonts) {
		    //if (font != BadName) XUnloadFont(dis,font);
		    stringstream fontname;
		    #ifdef DEBUG
		    cerr << font << ' ' << BadName << ' ' << BadAlloc << endl;
		    #endif //DEBUG
		    fontname << "*x" << size;
		    #ifdef DEBUG
		    cerr << fontname.str()  << endl;
		    #endif //DEBUG
		    //char** possible_fonts = XListFonts(dis, fontname.str().c_str(),5,&Nfonts);
		    XListFonts(dis, fontname.str().c_str(),5,&Nfonts);
		    #ifdef DEBUG
		    cerr << "N: " << Nfonts  << endl;
		    #endif //DEBUG
		    if (Nfonts)
			font = XLoadFont(dis, fontname.str().c_str());    
		    else {
			if (((size < 5 && !increase) || (size > 30 && increase)) && end) {
			    fontname << "*x*";
			    font = XLoadFont(dis, fontname.str().c_str());
			    break;
			}
			else if (size > 30 && increase) { --size; increase=false; end = true; }
			else if (size < 5 && !increase)  { ++size; increase=true; end = true; }
			else if (increase) ++size;
			else --size;
		    }
		}
		#ifdef DEBUG
		cerr << "Final font: " << font << ' ' << BadName << endl;
		#endif //DEBUG
		#ifdef DEBUG
		int FontSet = XSetFont(dis,gc,font);
		cerr << FontSet << endl;
		#else
		XSetFont(dis,gc,font);
		#endif //DEBUG
	    }
	///////////
	    if (col && col != present_col) {
		set_FG_col(col);
	    }
	    else if (present_col != default_col) {
		set_default_FG();
	    }
	    XDrawString(dis, win, gc, (x-x_setoff)*scale, (y-y_setoff)*scale, text.c_str(),text.length());
	}
	int next_event () {	
	    XNextEvent(dis, &report);
	    return report.type;
	}
	void clear_window() {
	    XClearWindow(dis, win);
	}
	void add_to_x_setoff( int value ) { x_setoff += value; }
	void add_to_y_setoff( int value ) { y_setoff += value; }
	void change_scale( float value ) { scale *= value; }
    private:
	// Graphical variables
	Display* dis;
	int screen;
	Window win;
	XEvent report;
	// colours
	GC gc;
	unsigned long black_int;
	unsigned long white_int;
	unsigned long default_col;
	unsigned long present_col;
	Colormap colormap;
	char hex_col[8];
	// dislocation
	int x_setoff;
	int y_setoff;
	float scale;
	// Text
	Font font;
};*/

const char DEFAULT_COL[] = "#000000";

void set_svg_colour_hex( char* hex, const char* col );

void display_drawing( app_display& window, drawing& fig ); 
bool pars_svg(istream& input, drawing& fig, app_display& window);
void help();
int main (int argc, char *argv[]) {
    // Variables for general behaviour
    bool quiet(true);
    string file_name;
    ifstream infile;
    istream* input = &cin;
    // Pars arguments
    for (int i=1; i < argc; ++i) {
	if ( !strcmp(argv[i],"-f") || !strcmp(argv[i],"--file") ) {
	    if (i+1 < argc && argv[i+1][0] != '-') {
		file_name = argv[++i];
	    }
	    else {
		cerr << "-f/--file needs to be followed by a file name." << endl;
		return 1;
	    }
	}
	else if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"--verbose")) quiet = false;
	else if ( !strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") ) {
	    help();
	    return 0;
	}
	else if ( i == argc-1 && argv[i][0] != '-' && file_name.empty() ) {
	    file_name = argv[i];
	}
	else if (i < argc) {
	    std::cerr << "The program was called with the following command:" << endl;
	    for (int j=0; j<argc; ++j) std::cerr << argv[j] << ' ';
	    std::cerr << endl;
	    std::cerr << "Argument " << argv[i] << " not recognized. For available arguments give -h or --help." << endl;
	    return 0;
	}

    }
    // Open infile if one given
    if (!file_name.empty()) {
	infile.open(file_name.c_str(), std::ifstream::in);
	if (infile.good())
	    input = &infile;
	else {
	    cerr << "Could not open file: " << file_name << endl;
	    return 1;
	}
    }
    // declaring app window
    app_display window;

    // Pars SVG
    drawing figure;
    if (pars_svg(*input, figure, window) && !quiet) std::cerr << "Parsed figure." << std::endl;
    else if (!quiet) {
        std::cerr << "Unable to pars figure." << std::endl;
	return 1;
    }

    // Open display
    window.start_win();
    window.set_default_col("#000000");
    window.set_default_FG();

    #ifdef DEBUG
    cerr << "Entering loop" << endl;
    #endif // DEBUG
    if (!quiet) {
	cerr << "Zoom in with the + or p key and out with the - or m key, move right and left" << endl;
	cerr << "with arrow keys (Swedish keyboard)." << endl;
    }
    while(1) {
	int happening = window.next_event();
	if (happening == Expose) {
                display_drawing( window, figure );
	}
	else if (happening == KeyPress) {
		#ifdef DEBUG
		    std::cout << "Key presed: " << XLookupKeysym(&report.xkey,0) << std::endl;
		#endif //DEBUG
                if (window.get_key_event_symbol()==XK_space || window.get_key_event_symbol()==65307 || window.get_key_event_symbol()==113) { // space + esc + q?
                    break;
                }
		else if (window.get_key_event_symbol() == 65364) { // down
		    window.add_to_y_setoff( -10);
		    window.clear_window();
		    display_drawing( window, figure );
		}
		else if (window.get_key_event_symbol() == 65362) { // up
		    window.add_to_y_setoff( 10 );
		    window.clear_window();
		    display_drawing( window, figure );
		}
		else if (window.get_key_event_symbol() == 65363) { // right
		    window.add_to_x_setoff(-10);
		    window.clear_window();
		    display_drawing( window, figure );
		}
		else if (window.get_key_event_symbol() == 65361) { // left
		    window.add_to_x_setoff(10);
		    window.clear_window();
		    display_drawing( window, figure );
		}
		else if (window.get_key_event_symbol() == 43 || window.get_key_event_symbol() == 112) { // +? or p
		    window.change_scale( 1.1 );
		    window.clear_window();
		    display_drawing( window, figure );
		}
		else if (window.get_key_event_symbol() == 45 || window.get_key_event_symbol() == 109) { //-? or m
		    window.change_scale( 0.9 );
		    window.clear_window();
		    display_drawing( window, figure );
		}
	        /*else {
		    std::cout << window.get_key_event_symbol() << std::endl;
		}*/	
	}
	else if ( happening == ClientMessage)
	    break;
    }
    return(0);
}

void set_svg_colour_hex( char* hex, const char* col ) {
    if (col == 0 || col[0] == '\0') strcpy(hex,"#000000");
    if (col[0] == '#') {
	unsigned int c(0);
	while (col[c]) ++c;
	if (c == 7) strcpy(hex,col);
    }
    unsigned int i(0);
    if (strcmp(col, "black") == 0) strcpy(hex,"#000000");
    else if (strcmp(col, "navy") == 0) strcpy(hex,"#000080");
    else if (strcmp(col, "darkblue") == 0) strcpy(hex,"#00008B");
    else if (strcmp(col, "mediumblue") == 0) strcpy(hex,"#0000CD");
    else if (strcmp(col, "blue") == 0) strcpy(hex,"#0000FF");
    else if (strcmp(col, "darkgreen") == 0) strcpy(hex,"#006400");
    else if (strcmp(col, "green") == 0) strcpy(hex,"#008000");
    else if (strcmp(col, "teal") == 0) strcpy(hex,"#008080");
    else if (strcmp(col, "darkcyan") == 0) strcpy(hex,"#008B8B");
    else if (strcmp(col, "deepskyblue") == 0) strcpy(hex,"#00BFFF");
    else if (strcmp(col, "darkturquoise") == 0) strcpy(hex,"#00CED1");
    else if (strcmp(col, "mediumspringgreen") == 0) strcpy(hex,"#00FA9A");
    else if (strcmp(col, "lime") == 0) strcpy(hex,"#00FF00");
    else if (strcmp(col, "springgreen") == 0) strcpy(hex,"#00FF7F");
    else if (strcmp(col, "cyan") == 0) strcpy(hex,"#00FFFF");
    else if (strcmp(col, "aqua") == 0) strcpy(hex,"#00FFFF");
    else if (strcmp(col, "midnightblue") == 0) strcpy(hex,"#191970");
    else if (strcmp(col, "dodgerblue") == 0) strcpy(hex,"#1E90FF");
    else if (strcmp(col, "lightseagreen") == 0) strcpy(hex,"#20B2AA");
    else if (strcmp(col, "forestgreen") == 0) strcpy(hex,"#228B22");
    else if (strcmp(col, "seagreen") == 0) strcpy(hex,"#2E8B57");
    else if (strcmp(col, "darkslategray") == 0) strcpy(hex,"#2F4F4F");
    else if (strcmp(col, "darkslategrey") == 0) strcpy(hex,"#2F4F4F");
    else if (strcmp(col, "limegreen") == 0) strcpy(hex,"#32CD32");
    else if (strcmp(col, "mediumseagreen") == 0) strcpy(hex,"#3CB371");
    else if (strcmp(col, "turquoise") == 0) strcpy(hex,"#40E0D0");
    else if (strcmp(col, "royalblue") == 0) strcpy(hex,"#4169E1");
    else if (strcmp(col, "steelblue") == 0) strcpy(hex,"#4682B4");
    else if (strcmp(col, "darkslateblue") == 0) strcpy(hex,"#483D8B");
    else if (strcmp(col, "mediumturquoise") == 0) strcpy(hex,"#48D1CC");
    else if (strcmp(col, "indigo") == 0) strcpy(hex,"#4B0082");
    else if (strcmp(col, "darkolivegreen") == 0) strcpy(hex,"#556B2F");
    else if (strcmp(col, "cadetblue") == 0) strcpy(hex,"#5F9EA0");
    else if (strcmp(col, "cornflowerblue") == 0) strcpy(hex,"#6495ED");
    else if (strcmp(col, "mediumaquamarine") == 0) strcpy(hex,"#66CDAA");
    else if (strcmp(col, "dimgrey") == 0) strcpy(hex,"#696969");
    else if (strcmp(col, "dimgray") == 0) strcpy(hex,"#696969");
    else if (strcmp(col, "slateblue") == 0) strcpy(hex,"#6A5ACD");
    else if (strcmp(col, "olivedrab") == 0) strcpy(hex,"#6B8E23");
    else if (strcmp(col, "slategrey") == 0) strcpy(hex,"#708090");
    else if (strcmp(col, "slategray") == 0) strcpy(hex,"#708090");
    else if (strcmp(col, "lightslategray") == 0) strcpy(hex,"#778899");
    else if (strcmp(col, "lightslategrey") == 0) strcpy(hex,"#778899");
    else if (strcmp(col, "mediumslateblue") == 0) strcpy(hex,"#7B68EE");
    else if (strcmp(col, "lawngreen") == 0) strcpy(hex,"#7CFC00");
    else if (strcmp(col, "chartreuse") == 0) strcpy(hex,"#7FFF00");
    else if (strcmp(col, "aquamarine") == 0) strcpy(hex,"#7FFFD4");
    else if (strcmp(col, "maroon") == 0) strcpy(hex,"#800000");
    else if (strcmp(col, "purple") == 0) strcpy(hex,"#800080");
    else if (strcmp(col, "olive") == 0) strcpy(hex,"#808000");
    else if (strcmp(col, "gray") == 0) strcpy(hex,"#808080");
    else if (strcmp(col, "grey") == 0) strcpy(hex,"#808080");
    else if (strcmp(col, "skyblue") == 0) strcpy(hex,"#87CEEB");
    else if (strcmp(col, "lightskyblue") == 0) strcpy(hex,"#87CEFA");
    else if (strcmp(col, "blueviolet") == 0) strcpy(hex,"#8A2BE2");
    else if (strcmp(col, "darkred") == 0) strcpy(hex,"#8B0000");
    else if (strcmp(col, "darkmagenta") == 0) strcpy(hex,"#8B008B");
    else if (strcmp(col, "saddlebrown") == 0) strcpy(hex,"#8B4513");
    else if (strcmp(col, "darkseagreen") == 0) strcpy(hex,"#8FBC8F");
    else if (strcmp(col, "lightgreen") == 0) strcpy(hex,"#90EE90");
    else if (strcmp(col, "mediumpurple") == 0) strcpy(hex,"#9370DB");
    else if (strcmp(col, "darkviolet") == 0) strcpy(hex,"#9400D3");
    else if (strcmp(col, "palegreen") == 0) strcpy(hex,"#98FB98");
    else if (strcmp(col, "darkorchid") == 0) strcpy(hex,"#9932CC");
    else if (strcmp(col, "yellowgreen") == 0) strcpy(hex,"#9ACD32");
    else if (strcmp(col, "sienna") == 0) strcpy(hex,"#A0522D");
    else if (strcmp(col, "brown") == 0) strcpy(hex,"#A52A2A");
    else if (strcmp(col, "darkgray") == 0) strcpy(hex,"#A9A9A9");
    else if (strcmp(col, "darkgrey") == 0) strcpy(hex,"#A9A9A9");
    else if (strcmp(col, "lightblue") == 0) strcpy(hex,"#ADD8E6");
    else if (strcmp(col, "greenyellow") == 0) strcpy(hex,"#ADFF2F");
    else if (strcmp(col, "paleturquoise") == 0) strcpy(hex,"#AFEEEE");
    else if (strcmp(col, "lightsteelblue") == 0) strcpy(hex,"#B0C4DE");
    else if (strcmp(col, "powderblue") == 0) strcpy(hex,"#B0E0E6");
    else if (strcmp(col, "firebrick") == 0) strcpy(hex,"#B22222");
    else if (strcmp(col, "darkgoldenrod") == 0) strcpy(hex,"#B8860B");
    else if (strcmp(col, "mediumorchid") == 0) strcpy(hex,"#BA55D3");
    else if (strcmp(col, "rosybrown") == 0) strcpy(hex,"#BC8F8F");
    else if (strcmp(col, "darkkhaki") == 0) strcpy(hex,"#BDB76B");
    else if (strcmp(col, "silver") == 0) strcpy(hex,"#C0C0C0");
    else if (strcmp(col, "mediumvioletred") == 0) strcpy(hex,"#C71585");
    else if (strcmp(col, "indianred") == 0) strcpy(hex,"#CD5C5C");
    else if (strcmp(col, "peru") == 0) strcpy(hex,"#CD853F");
    else if (strcmp(col, "chocolate") == 0) strcpy(hex,"#D2691E");
    else if (strcmp(col, "tan") == 0) strcpy(hex,"#D2B48C");
    else if (strcmp(col, "lightgray") == 0) strcpy(hex,"#D3D3D3");
    else if (strcmp(col, "lightgrey") == 0) strcpy(hex,"#D3D3D3");
    else if (strcmp(col, "thistle") == 0) strcpy(hex,"#D8BFD8");
    else if (strcmp(col, "orchid") == 0) strcpy(hex,"#DA70D6");
    else if (strcmp(col, "goldenrod") == 0) strcpy(hex,"#DAA520");
    else if (strcmp(col, "palevioletred") == 0) strcpy(hex,"#DB7093");
    else if (strcmp(col, "crimson") == 0) strcpy(hex,"#DC143C");
    else if (strcmp(col, "gainsboro") == 0) strcpy(hex,"#DCDCDC");
    else if (strcmp(col, "plum") == 0) strcpy(hex,"#DDA0DD");
    else if (strcmp(col, "burlywood") == 0) strcpy(hex,"#DEB887");
    else if (strcmp(col, "lightcyan") == 0) strcpy(hex,"#E0FFFF");
    else if (strcmp(col, "lavender") == 0) strcpy(hex,"#E6E6FA");
    else if (strcmp(col, "darksalmon") == 0) strcpy(hex,"#E9967A");
    else if (strcmp(col, "violet") == 0) strcpy(hex,"#EE82EE");
    else if (strcmp(col, "palegoldenrod") == 0) strcpy(hex,"#EEE8AA");
    else if (strcmp(col, "lightcoral") == 0) strcpy(hex,"#F08080");
    else if (strcmp(col, "khaki") == 0) strcpy(hex,"#F0E68C");
    else if (strcmp(col, "aliceblue") == 0) strcpy(hex,"#F0F8FF");
    else if (strcmp(col, "honeydew") == 0) strcpy(hex,"#F0FFF0");
    else if (strcmp(col, "azure") == 0) strcpy(hex,"#F0FFFF");
    else if (strcmp(col, "sandybrown") == 0) strcpy(hex,"#F4A460");
    else if (strcmp(col, "wheat") == 0) strcpy(hex,"#F5DEB3");
    else if (strcmp(col, "beige") == 0) strcpy(hex,"#F5F5DC");
    else if (strcmp(col, "whitesmoke") == 0) strcpy(hex,"#F5F5F5");
    else if (strcmp(col, "mintcream") == 0) strcpy(hex,"#F5FFFA");
    else if (strcmp(col, "ghostwhite") == 0) strcpy(hex,"#F8F8FF");
    else if (strcmp(col, "salmon") == 0) strcpy(hex,"#FA8072");
    else if (strcmp(col, "antiquewhite") == 0) strcpy(hex,"#FAEBD7");
    else if (strcmp(col, "linen") == 0) strcpy(hex,"#FAF0E6");
    else if (strcmp(col, "lightgoldenrodyellow") == 0) strcpy(hex,"#FAFAD2");
    else if (strcmp(col, "oldlace") == 0) strcpy(hex,"#FDF5E6");
    else if (strcmp(col, "red") == 0) strcpy(hex,"#FF0000");
    else if (strcmp(col, "fuchsia") == 0) strcpy(hex,"#FF00FF");
    else if (strcmp(col, "magenta") == 0) strcpy(hex,"#FF00FF");
    else if (strcmp(col, "deeppink") == 0) strcpy(hex,"#FF1493");
    else if (strcmp(col, "orangered") == 0) strcpy(hex,"#FF4500");
    else if (strcmp(col, "hotpink") == 0) strcpy(hex,"#FF69B4");
    else if (strcmp(col, "coral") == 0) strcpy(hex,"#FF7F50");
    else if (strcmp(col, "darkorange") == 0) strcpy(hex,"#FF8C00");
    else if (strcmp(col, "lightsalmon") == 0) strcpy(hex,"#FFA07A");
    else if (strcmp(col, "orange") == 0) strcpy(hex,"#FFA500");
    else if (strcmp(col, "lightpink") == 0) strcpy(hex,"#FFB6C1");
    else if (strcmp(col, "pink") == 0) strcpy(hex,"#FFC0CB");
    else if (strcmp(col, "gold") == 0) strcpy(hex,"#FFD700");
    else if (strcmp(col, "peachpuff") == 0) strcpy(hex,"#FFDAB9");
    else if (strcmp(col, "navajowhite") == 0) strcpy(hex,"#FFDEAD");
    else if (strcmp(col, "moccasin") == 0) strcpy(hex,"#FFE4B5");
    else if (strcmp(col, "mistyrose") == 0) strcpy(hex,"#FFE4E1");
    else if (strcmp(col, "blanchedalmond") == 0) strcpy(hex,"#FFEBCD");
    else if (strcmp(col, "papayawhip") == 0) strcpy(hex,"#FFEFD5");
    else if (strcmp(col, "lavenderblush") == 0) strcpy(hex,"#FFF0F5");
    else if (strcmp(col, "seashell") == 0) strcpy(hex,"#FFF5EE");
    else if (strcmp(col, "cornsilk") == 0) strcpy(hex,"#FFF8DC");
    else if (strcmp(col, "lemonchiffon") == 0) strcpy(hex,"#FFFACD");
    else if (strcmp(col, "floralwhite") == 0) strcpy(hex,"#FFFAF0");
    else if (strcmp(col, "snow") == 0) strcpy(hex,"#FFFAFA");
    else if (strcmp(col, "yellow") == 0) strcpy(hex,"#FFFF00");
    else if (strcmp(col, "lightyellow") == 0) strcpy(hex,"#FFFFE0");
    else if (strcmp(col, "ivory") == 0) strcpy(hex,"#FFFFF0");
    else if (strcmp(col, "white") == 0) strcpy(hex,"#FFFFFF");
    else strcpy(hex,"#000000");
}

bool pars_svg(istream& input, drawing& fig, app_display& window){
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load(input);

    if(result) {
	unsigned int n_objects(0);
        pugi::xml_node objects = doc.child("svg");
	for (pugi::xml_node_iterator object = objects.begin(); object != objects.end(); ++object) {
            if (!strcmp("line", object->name())) {
	        int x1 = atoi(object->attribute("x1").value());
	        int y1 = atoi(object->attribute("y1").value());
	        int x2 = atoi(object->attribute("x2").value());
	        int y2 = atoi(object->attribute("y2").value());
		fig.add_line(x1,y1,x2,y2);
		++n_objects;
	    }
	    else if (!strcmp("polyline", object->name())) {
                const std::string points = object->attribute("points").value();
		unsigned int length = points.length();
                std::string value;
		int x1=0;
		int y1=0;
		int x2=0;
		int y2=0;
		bool start = true;
		for (unsigned int i=0; i<=length; ++i) {
                    if ((points[i] == ' ' || i == length) && !value.empty()) {
			if (!start) {
                            y2 = atoi(value.c_str());
                            fig.add_line(x1,y1,x2,y2);
			    x1 = x2;
			    y1 = y2;
			}
			else {
			    y1 = atoi(value.c_str());
			    start = false;
			}
			value.clear();
		    }
		    else if (points[i] == ',') {
			if (value.empty()) {
                            std::cerr << "Polyline x value parsing error." << std::endl;
			    return false;
			}
                        if (start) {
			    x1 = atoi(value.c_str());
			}
			else x2 = atoi(value.c_str());
			value.clear();
		    }
		    else value+=points[i];
		}
		++n_objects;
	    }
	    else if (!strcmp("text", object->name())) {
	        int x = atoi(object->attribute("x").value());
	        int y = atoi(object->attribute("y").value());
		/// Get colour in hex
		char hex_col[8];
		string style(object->attribute("style").value());
		string col(object->attribute("fill").value());
		if (!style.empty()) {
		    string key("");
		    string value("");
		    char mode('k'); 
		    for (unsigned int i=0; i < style.size(); ++i) {
			if (style[i] == ':') mode = 'v';
			else if (style[i] == ';') {
			    mode = 'k';
			    if (col.empty() && key.compare("fill") == 0 && !value.empty()) {
				col = value;
			    }
			    key = "";
			    value = "";
			}
			else if (mode == 'k') key += style[i];
			else if (mode == 'v') value += style[i];
		    }
		}
		if (!col.empty()) {
		    for (string::size_type i=0; i < col.size(); ++i) col[i] = tolower(col[i]);
		    set_svg_colour_hex(hex_col,col.c_str());
		}
		else strcpy(hex_col,DEFAULT_COL);
		/////
		unsigned int font_size = atoi(object->attribute("font-size").value());
		string font = object->attribute("font").value();
		string text_string = object->child_value();
		fig.add_text(x,y,text_string,font,font_size,window.hex_to_int_col(hex_col));
		++n_objects;
	    }
	    else std::cerr << "Unknown object: " << object->name() << std::endl;
	    if (!(n_objects % 1000)) std::cout << "Added " << n_objects << " objects." << std::endl;
	}
    }
    else return false;
    return true;
}

void display_drawing( app_display& window, drawing& fig ) {
    fig.init();
    drawing::line line;
    //unsigned long present_col = default_col;
    while (fig.more_lines()) {
        line = fig.get_line();
        window.draw_line(line.start.x, line.start.y, line.end.x, line.end.y);
    }
    drawing::rectangle rectangle;
    while (fig.more_rectangles()){
        rectangle = fig.get_rectangle();
	window.draw_rectangle(rectangle.upper_left.x, rectangle.upper_left.y, rectangle.width, rectangle.height);
    }
    drawing::ellipse ellipse;
    while (fig.more_ellipses()) {
        ellipse = fig.get_ellipse();
	window.draw_arc(ellipse.upper_left.x, ellipse.upper_left.y, ellipse.width, ellipse.height, ellipse.start_angle, ellipse.end_angle);
    }
    drawing::texttest text;
    // Font test
    //unsigned int font_size(0);
    //unsigned int size(0);
    while (fig.more_texts()) {
        text = fig.get_text();
	window.draw_string(text.characters, text.start.x, text.start.y, text.font_size, text.x11col);
    }
}

void help() {
    cout << "Rudisvg 0.1.0 is a rudimentary svg viewer. It can read and draw lines," << endl;
    cout << "rectangles, ellipsoids, and text. It only display in one color. Zoom in with the" << endl;
    cerr << "+ or p key and out with the - or m key, move right and left with arrow keys" << endl;
    cout << "(Swedish keyboard)." << endl;
    cout << "(c) Martin Ryberg 2016." << endl << endl;
    cout << "Usage:" << endl << "rudisvg < file.svg" << endl << "rudisvg file.svg" << endl << endl;
    cout << "Arguments:" << endl;
    cout << "--file / -f [file]    give svg file name." << endl;
    cout << "--help / -h           print this help." << endl;
    cout << "--verbose / -v        give additional output (to command line)." << endl;

}
