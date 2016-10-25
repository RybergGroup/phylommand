#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>

using namespace std;

class drawing {
    public:
        struct coordinat {
            int x;
            int y;
        };
        struct line {
            coordinat start;
            coordinat end;
        };
        struct ellipse {
            coordinat upper_left;
            unsigned int width;
            unsigned int height;
            int start_angle;
            int end_angle;
        };
        struct rectangle {
            coordinat upper_left;
            unsigned int width;
            unsigned int height;
        };
	struct texttest {
            coordinat start;
	    string characters;
	    unsigned int font_size;
	    string font;
	};
        void init() {
	    line_it = lines.begin();
	    rectangle_it = rectangles.begin();
	    ellipse_it = ellipses.begin();
            text_it = texts.begin();
	};
        bool more_lines() {
            if (line_it == lines.end()) return false;
            else return true;
        };
	bool more_rectangles() {
            if (rectangle_it == rectangles.end()) return false;
	    else return true;
	};
	bool more_ellipses() {
            if (ellipse_it == ellipses.end()) return false;
	    else return true;
	};
	bool more_texts() {
            if(text_it == texts.end()) return false;
	    else return true;
	};
        line get_line(){
	    ++line_it;
	    return *(line_it-1);
	};
	rectangle get_rectangle(){
            ++rectangle_it;
	    return *(rectangle_it-1);
	};
	ellipse get_ellipse(){
            ++ellipse_it;
	    return *(ellipse_it-1);
	};
	texttest get_text(){
            ++text_it;
	    return *(text_it-1);
	};
        void add_line(const int start_x, const int start_y, const int end_x, const int end_y) {
            line add;
            add.start.x = start_x;
            add.start.y = start_y; //comment
            add.end.x = end_x;
            add.end.y = end_y;
            lines.push_back(add);
	};
	void add_rectangle(const int x, const int y, const unsigned int width, const unsigned int height){
	    rectangle add;
	    add.upper_left.x = x;
	    add.upper_left.y = y;
	    add.width = width;
	    add.height = height;
	    rectangles.push_back(add);
	};
	void add_ellipse(const int x, const int y, const unsigned int width, const unsigned int height, const int start_angle, const int end_angle){
	    ellipse add;
	    add.upper_left.x = x;
	    add.upper_left.y = y;
	    add.width = width;
	    add.height = height;
	    add.start_angle = start_angle;
	    add.end_angle = end_angle;
	    ellipses.push_back(add);
	};
	void add_text(const int x, const int y, const string& characters, const string& font, const unsigned int font_size) {
	    texttest add;
	    add.start.x = x;
	    add.start.y = y;
	    add.font = font;
	    add.font_size = font_size;
            add.characters = characters;
            texts.push_back(add);
	};
    private:
        vector<line> lines;
	vector<rectangle> rectangles;
	vector<ellipse> ellipses;
	vector<texttest> texts;
        vector<line>::iterator line_it;
	vector<rectangle>::iterator rectangle_it;
	vector<ellipse>::iterator ellipse_it;
	vector<texttest>::iterator text_it;
};
