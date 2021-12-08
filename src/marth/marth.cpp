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

#include "marth.h"

bool marth::square_matrix::set_value (unsigned int x, unsigned int y, double value) {
    if (x>=dimentions || y>=dimentions)
        return false;
    else {
        values[x][y] = value;
        return true;
    };
}

double marth::square_matrix::get_value (unsigned int x, unsigned int y) {
    if (x>=dimentions || y>=dimentions)
        return 0.0;
    else
        return values [x][y];
}

bool marth::square_matrix::copy_matrix( const square_matrix* other){
    this->delete_matrix();
    this->dimentions = other->dimentions;
    this->create_matrix();
    for (unsigned int i=0; i < dimentions; ++i){
        for (unsigned int j=0; j < dimentions; ++j){
            this->values[i][j] = other->values[i][j];
        }
    }
    return true;
}

bool marth::square_matrix::multiply (const square_matrix* other) {
    if (this->dimentions != other->dimentions) return false;
    square_matrix prod(this->dimentions);
    for (unsigned int i=0; i<dimentions; ++i) {
        for (unsigned int j=0; j<dimentions; ++j) {
	    prod.values[i][j] = 0;
	    for (unsigned int k=0; k<dimentions; ++k)
                prod.values[i][j] += values[i][k]*other->values[k][j];
        }
    }
    return this->copy_matrix(&prod);
}

void marth::square_matrix::multiply ( const double k ) {
    for (unsigned int i=0; i < dimentions; ++i) {
        for (unsigned int j=0; j < dimentions; ++j) {
	    values[i][j] *= k;
	}
    }
}

bool marth::square_matrix::add ( const square_matrix* other) {
    if (this->dimentions != other->dimentions) return false;
    for (unsigned int i=0; i < dimentions; ++i) {
        for (unsigned int j=0; j < dimentions; ++j) {
	    values[i][j] += other->values[i][j];
        }
    }
    return true;
}

void marth::square_matrix::add ( const double k ) {
    for (unsigned int i=0; i < dimentions; ++i) {
        for (unsigned int j=0; j < dimentions; ++j) {
            values[i][j] += k;
        }
    }
}

void marth::square_matrix::identity () {
    for (unsigned int i=0; i < dimentions; ++i) {
        for (unsigned int j=0; j < dimentions; ++j) {
            if (i==j) values[i][j] = 1;
	    else values[i][j] = 0;
	}
    }
}

void marth::square_matrix::diagonal ( double* vector ) {
    for (unsigned int i=0; i < dimentions; ++i) {
        for (unsigned int j=0; j < dimentions; ++j) {
            if (i==j) vector[i] = values[i][j];
        }
    }
}

void marth::square_matrix::power (const unsigned int k) {
    square_matrix other;
    other.copy_matrix(this);
    for (unsigned int i=1; i < k; ++i)
        multiply(&other);
}

void marth::square_matrix::exponential ( square_matrix* exponential, const double r ) {
    exponential->delete_matrix();
    exponential->dimentions = this->dimentions;
    exponential->create_matrix();
    exponential->identity();
    square_matrix Q(this->dimentions);
    for (unsigned int i=1; i < 21; ++i) {
	for (unsigned int x=0; x < Q.dimentions; ++x)
	    for (unsigned int y=0; y < Q.dimentions; ++y)
		Q.values[x][y] = values[x][y] * r;
	Q.power(i);
	Q.multiply(1.00/factorial(i));
	exponential->add(&Q);
    }
}

void marth::square_matrix::exponential ( square_matrix* exponential, const double r, unsigned int k ) {
    exponential->delete_matrix();
    exponential->dimentions = this->dimentions;
    exponential->create_matrix();
    exponential->identity();
    square_matrix Q(this->dimentions);
    for (unsigned int x=0; x < Q.dimentions; ++x)
	for (unsigned int y=0; y < Q.dimentions; ++y)
	    Q.values[x][y] = (values[x][y] * r)/pow(2,k);
    exponential->add(&Q);
    for (unsigned int i=0; i < k; ++i) exponential->power(2);
}

void marth::square_matrix::print(std::ostream& out_stream) {
    for (unsigned int x=0; x < dimentions; ++x){
        for (unsigned int y=0; y < dimentions; ++y){
            out_stream << values[x][y];
            if (y < dimentions-1) out_stream << ' ';
        }
        if (x < dimentions-1) out_stream << std::endl;
    }
}

/*** General functions ***/

unsigned long int marth::factorial (const unsigned int k) {
    if (k == 0) return 1;
    unsigned long int prod = k;
    for (unsigned int i = k-1; i > 0; --i)
	prod *= i;
    return prod;
}

/*** Probability density distributions ***/
double inline marth::gammad(double x, int k, double b) { //Takes value for which to calculate probability (x), shape (k), and inv. scale (b).
   return pow(x,k-1)*(exp(-x/(1/b))/(pow(1/b,k)*factorial(k-1)));
}

double inline marth::uniformd (double x, double a, double b) {
   if (a > b) return 0;
   if (x < a) return 0;
   if (x > b) return 0;
   else return 1/(b-a);
}

double inline marth::normald (double x, double mean, double sd) {
   if (sd==0 && x==mean) return 1;
   else if (sd==0 && x!=mean) return 0;
   else return (1/sqrt(2*3.14159265*pow(sd,2)))*exp(-(pow(x-mean,2)/(2*pow(sd,2))));
}

/***/

/* Random values */
void inline marth::init_rand () { srand(time(NULL)); }; // initiates random number generatior from clock
void inline marth::init_rand ( unsigned int seed ) { srand(seed); }; // initiates random number generator with given number

double marth::normalr() {
   double s=2;
   double x;
   double y;
   while (!(s<1)) {
      x = (rand()/((double)RAND_MAX/2))-1;
      y = (rand()/((double)RAND_MAX/2))-1;
      s = (pow(x,2)+pow(y,2));
   }
   return x*sqrt(-2*log(s)/s);
}

double inline marth::uniformr(double a=0.0, double b=1.0) {
    if (a>b) return 0.0;
    return ((rand()/double(RAND_MAX))*(b-a))+a;
}

double inline marth::uniformr () { return uniformr( 0.0, 1.0 ); }; // draw random value between 0 and 1 from a uniform distribution
double marth::expr ( double mean ) { // draw random value from exponential distribution
    if ( mean < 0.0 )  return 0.0;
    return -log(rand()/double(RAND_MAX))/mean;
};

