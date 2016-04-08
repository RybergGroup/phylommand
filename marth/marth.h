#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

namespace marth {
    // Classes
    class square_matrix { // class to store square matrices with values in double
	public:
        square_matrix () {
            dimentions = 4;
            create_matrix();
        };
        square_matrix (const unsigned int dim) {
            dimentions = dim;
            create_matrix(); 
        };
        ~square_matrix () {
            delete_matrix();
        };
        bool set_value (const unsigned int x, const unsigned int y, const double value);
        double get_value (const unsigned int x, const unsigned int y);
	int get_dimentions () { return dimentions; };
        bool multiply (const square_matrix* other);
	void multiply ( const double k );
	bool add ( const square_matrix* other);
	void add ( const double k );
        void power(const unsigned int k);
	void exponential ( square_matrix* exponential, const double r ); // is not correct for long branches
	void exponential ( square_matrix* exponential, const double r, unsigned int k ); // recommended
	void identity ();
	void diagonal ( double* vector ); // require pointer to vector that is as long as dimentions
        void print();
        bool copy_matrix( const square_matrix* other);
	bool reset( const unsigned int dim ) {
	    if (delete_matrix()) {
		dimentions = dim;
		return create_matrix();
	    }
	    else return false;
	}
	private:
        double** values;
        unsigned int dimentions;
        bool create_matrix() {
            values = new double* [dimentions];
            values[0] = new double [dimentions*dimentions];
            for (unsigned int i=1; i<dimentions; ++i)
                values[i]=values[0]+i*dimentions;
	    return true;
        };
        bool delete_matrix () {
            delete [] values[0];
            delete [] values;
	    return true;
        };
    };

    /*** General functions ***/
    unsigned long int factorial (const unsigned int k); // does not work for numbers over 20

    // Probability density distributions
    double gammad(double x, int k, double b); // returns probability density for x given shape (k) and scale (b) parameters
    double uniformd(double x, double a, double b); // Returns the probability density for x given a uniform distribution between a and b
    double normald (double x, double mean, double sd); // Returns the probability density for x for given mean and standart deviation

    // Random number distributions
    void init_rand (); // initiates random number generatior from clock
    void init_rand ( unsigned int seed ); // initiates random number generator with given number
    double normalr (); // draw random value from standard normal distribution
    double uniformr (double a, double b); // draw random value from uniform distribution betwwen values a and b. If b < a p=0.0
    double uniformr (); // draw random value between 0 and 1 from a uniform distribution
    double expr ( double mean ); // draw random value from exponential distribution
};
