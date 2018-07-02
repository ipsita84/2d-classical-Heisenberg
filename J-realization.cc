//  g++ -std=c++11 -Wall -O3 J-realization.cc -o testo

#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen(std::time(0)); // time(0) changes seed every time you run
//not to change seed: use 
//boost::random::mt19937 gen;
using namespace std;

typedef
boost::multi_array < double, 2 > array_2d;
// typedef keyword allows you to create an alias fo a data type



int main()
{
	array_2d J(boost::extents[3][3]);
	double j_isotropic =- 20;
	double sym1 = -0 ;
	double sym2 = -0 ;
	double sym3 = -0 ;
	double antisym1 = -200 ;
	double antisym2 = -0 ;
	double antisym3 = -0 ;

	ofstream gout("J2.dat");	// Opens a file for output


	J[0][0] = j_isotropic ;
    J[1][1] = j_isotropic ;
    J[2][2] = j_isotropic ;


  	J[0][1] = sym1 +antisym1 ;
    J[1][0] = sym1 - antisym1 ;
    
  	J[0][2] = sym2 +antisym2 ;   
    J[2][0] = sym2 - antisym2 ;
   
  	J[2][1] = sym3 + antisym3 ;
    J[1][2] =  sym3 - antisym3;

for (unsigned int i = 0; i < 3; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			gout << J[i][j] << endl;
		}
}

	gout.close();
	return 0;
}



