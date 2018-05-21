// g++ -std=c++11 -Wall -O3 autocorrelation.cc -o test
//Reference: https://ocw.mit.edu/courses/chemistry/5-74-introductory-quantum-mechanics-ii-spring-2004/lecture-notes/11.pdf or https://arxiv.org/pdf/0905.1629.pdf

#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h> 
#include <array>

using namespace std;

const unsigned int N_mc = 1e6;

int main(int argc, char const * argv[])
{      ifstream gin("test.dat");
       std::array <double, N_mc> energy;
       for (unsigned int i=1; i <N_mc+1; ++i )
       {
			gin>> energy[i];

        }
	gin.close();
       
        ofstream fout("autocorr.dat");	// Opens a file for output
	
	double en_avg(0);
	double en_square_avg(0);
        for (unsigned int i=0; i < N_mc; ++i )
        {
			en_square_avg = en_square_avg + energy[i]*energy[i];
			en_avg = en_avg + energy[i];

         }
         en_avg = en_avg / N_mc ;
         en_square_avg = en_square_avg / N_mc ;
         cout << "avg energy" << en_avg <<endl ;
         cout << "avg energy square" << en_square_avg  << endl;
         
        double autosum(0); 
        for (unsigned int i=1; i < N_mc-1; ++i )
        {
        
               double sum(0);
               for (unsigned int j=0; j < N_mc -i; ++j )
               {
			sum += (energy[j]- en_avg)*(energy[j+i] - en_avg) ;

                } 
                fout<< i << '\t' << sum/ (N_mc - i) <<endl ;
                autosum +=  sum/(N_mc - i);

         }
         
        double autocorr_time =  autosum / (en_square_avg - en_avg * en_avg);
        cout << autocorr_time;

       
        fout << autocorr_time << endl;
        fout.close();

	return 0;
}
