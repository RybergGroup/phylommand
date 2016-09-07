#include <stdio.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <ctime>
#include <vector>
#include "pugixml/pugixml.hpp"
#include "drawing.h"

using namespace std;

void sleep(unsigned int mseconds) {
    clock_t goal = mseconds + clock();
    while (goal > clock());
}

void close_x(Display* dis, Window &win, GC &gc) {
    XFreeGC(dis, gc);
    XDestroyWindow(dis, win);
    XCloseDisplay(dis);
}

void display_drawing(Display* dis, Window& win, GC& gc, drawing& fig, int x_setoff, int y_setoff, float scale); 
bool pars_svg(istream& input, drawing& fig);

int main (int argc, char *argv[]) {
    // Variables for general behaviour
    bool quiet(false);

    // Pars SVG
    drawing figure;
    if (pars_svg(std::cin, figure) && !quiet) std::cerr << "Parsed figure." << std::endl;
    else if (!quiet) {
        std::cerr << "Unable to pars figure." << std::endl;
	return 1;
    }

    // Graphical variables
    Display* dis;
    int screen;
    Window win;
    XEvent report;
    GC green_gc;

    // Prepair window
    #ifdef DEBUG
    cerr << "Opening display" << endl;
    #endif // DEBUG
    dis = XOpenDisplay(NULL);
    screen = DefaultScreen(dis);

    // Prepair colors
    unsigned long black;//, white; 
    black = BlackPixel(dis, screen);
    //white = WhitePixel(dis, screen);
    XColor green_col;
    Colormap colormap;
    char green[] = "#00FF00";

    //variables to move around
    int x_setoff(0);
    int y_setoff(0);
    float scale(1.0);

    // Start window
    #ifdef DEBUG
    cerr << "Creating window" << endl;
    #endif // DEBUG
    win = XCreateSimpleWindow(dis,RootWindow(dis,0),1,1,500,500, 0, black, black);

    XSetStandardProperties(dis,win, "Incomplete SVG", "i!", None, NULL, 0, NULL);

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

    colormap = DefaultColormap(dis, 0);
    green_gc = XCreateGC(dis, win, 0, 0);


    XParseColor(dis, colormap, green, &green_col);
    XAllocColor(dis, colormap, &green_col);
    XSetForeground(dis, green_gc, green_col.pixel);

    #ifdef DEBUG
    cerr << "Entering loop" << endl;
    #endif // DEBUG
    if (!quiet) {
	cerr << "Zoom in and out with + and - key, move right and left with arrow keys (Swedish" << endl;
	cerr << "keyboard)." << endl;
    }
    while(1) {
        XNextEvent(dis, &report);
	if (report.type == Expose) {
                display_drawing(dis, win, green_gc, figure,x_setoff,y_setoff,scale);
	}
	else if (report.type == KeyPress) {
		#ifdef DEBUG
		    std::cout << "Key presed: " << XLookupKeysym(&report.xkey,0) << std::endl;
		#endif //DEBUG
                if (XLookupKeysym(&report.xkey,0)==XK_space || XLookupKeysym(&report.xkey,0)==65307 || XLookupKeysym(&report.xkey,0)==113) { // space + esc + q?
                    break;
                }
		else if (XLookupKeysym(&report.xkey,0) == 65364) { // neråt
		    y_setoff -= 10;
		    XClearWindow(dis, win);
		    display_drawing(dis, win, green_gc, figure,x_setoff,y_setoff,scale);
		}
		else if (XLookupKeysym(&report.xkey,0) == 65362) { // uppåt
		    y_setoff += 10;
		    XClearWindow(dis, win);
		    display_drawing(dis, win, green_gc, figure,x_setoff,y_setoff,scale);
		}
		else if (XLookupKeysym(&report.xkey,0) == 65363) { // höger
		    x_setoff -= 10;
		    XClearWindow(dis, win);
		    display_drawing(dis, win, green_gc, figure,x_setoff,y_setoff,scale);
		}
		else if (XLookupKeysym(&report.xkey,0) == 65361) { // vänster
		    x_setoff += 10;
		    XClearWindow(dis, win);
		    display_drawing(dis, win, green_gc, figure,x_setoff,y_setoff,scale);
		}
		else if (XLookupKeysym(&report.xkey,0) == 45) { // +?
		    scale *= 1.1;
		    XClearWindow(dis, win);
		    display_drawing(dis, win, green_gc, figure,x_setoff,y_setoff,scale);
		}
		else if (XLookupKeysym(&report.xkey,0) == 47) { //-?
		    scale *= 0.9;
		    XClearWindow(dis, win);
		    display_drawing(dis, win, green_gc, figure,x_setoff,y_setoff,scale);
		}
	        /*else {
		    std::cout << XLookupKeysym(&report.xkey,0) << std::endl;
		}*/	
	}
	else if (report.type == ClientMessage)
	    break;
    }
    close_x(dis, win, green_gc);
    return(0);
}

bool pars_svg(istream& input, drawing& fig){
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
		//std::cout << "Added line: " << x1 << ',' << y1 << ' ' << x2 << ',' << y2 << std::endl;
	    }
	    else if (!strcmp("polyline", object->name())) {
                const std::string points = object->attribute("points").value();
                //std::cout << points << endl;
		unsigned int length = points.length();
                std::string value;
		int x1=0;
		int y1=0;
		int x2=0;
		int y2=0;
		bool start = true;
		for (unsigned int i=0; i<=length; ++i) {
                    if ((points[i] == ' ' || i == length) && !value.empty()) {
			//std::cout << "x1,y1: " << x1 << ',' << y1 << "  x2, value: " << x2 << ',' << value << std::endl;
			if (!start) {
                            y2 = atoi(value.c_str());
                            fig.add_line(x1,y1,x2,y2);
			    //std::cout << "Added line: " << x1 << ',' << y1 << ' ' << x2 << ',' << y2 << std::endl;
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
			    //std::cout << "Read first x: " << x1 << std::endl;
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
		std::string text_string = object->child_value();
		//std::cout << "Splunch: " << object->child_value() << std::endl;
		fig.add_text(x,y,text_string);
		//std::cout << "Added text string: " << text_string << std::endl;
		++n_objects;
	    }
	    else std::cerr << "Unknown object: " << object->name() << std::endl;
	    if (!(n_objects % 1000)) std::cout << "Added " << n_objects << " objects." << std::endl;
	}
    }
    else return false;
    return true;
}

void display_drawing(Display* dis, Window& win, GC& gc, drawing& fig, int x_setoff, int y_setoff, float scale) {
    fig.init();
    drawing::line line;
    while (fig.more_lines()) {
        line = fig.get_line();
        XDrawLine(dis, win, gc, (line.start.x-x_setoff)*scale, (line.start.y-y_setoff)*scale, (line.end.x-x_setoff)*scale, (line.end.y-y_setoff)*scale);
    }
    drawing::rectangle rectangle;
    while (fig.more_rectangles()){
        rectangle = fig.get_rectangle();
	XDrawRectangle(dis, win, gc, (rectangle.upper_left.x-x_setoff)*scale, (rectangle.upper_left.y-y_setoff)*scale, rectangle.width*scale, rectangle.height*scale);
    }
    drawing::ellipse ellipse;
    while (fig.more_ellipses()) {
        ellipse = fig.get_ellipse();
	XDrawArc(dis, win, gc, (ellipse.upper_left.x-x_setoff)*scale, (ellipse.upper_left.y-y_setoff)*scale, ellipse.width*scale, ellipse.height*scale, ellipse.start_angle, ellipse.end_angle);
    }
    drawing::texttest text;
    // Font test
    Font font(BadName);
    unsigned int size = 10*scale;
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
    ///////////
    while (fig.more_texts()) {
        text = fig.get_text();
	XDrawString(dis, win, gc, (text.start.x-x_setoff)*scale, (text.start.y-y_setoff)*scale, text.characters.c_str(),text.characters.length());
    }
    XUnloadFont(dis,font);
}

