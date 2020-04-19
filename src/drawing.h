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
	    unsigned long x11col;
        };
        struct ellipse {
            coordinat upper_left;
            unsigned int width;
            unsigned int height;
            int start_angle;
            int end_angle;
	    unsigned long x11col;
        };
        struct rectangle {
            coordinat upper_left;
            unsigned int width;
            unsigned int height;
	    unsigned long x11col;
        };
	struct texttest {
            coordinat start;
	    string characters;
	    unsigned int font_size;
	    string font;
	    unsigned long x11col;
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
	void add_text(const int x, const int y, const string& characters, const string& font, const unsigned int font_size, unsigned long col) {
	    texttest add;
	    add.start.x = x;
	    add.start.y = y;
	    add.font = font;
	    add.font_size = font_size;
            add.characters = characters;
	    add.x11col=col;
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
