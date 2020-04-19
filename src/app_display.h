#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
using namespace std;

class app_display {
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
	int get_event_type() { return report.type; }
	KeySym get_key_event_symbol ( ) { return XLookupKeysym(&report.xkey,0); }
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
};

