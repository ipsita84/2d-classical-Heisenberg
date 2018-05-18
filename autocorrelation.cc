// g++ -std=c++11 -Wall -O3 autocorrelation.cc -o test

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
	
	double en_avg(0);
	double en_square_avg(0);
        for (unsigned int i=0; i < N_mc; ++i )
        {
			en_square_avg = en_square_avg + energy[i]*energy[i];
			en_avg = en_avg + energy[i];

         }
         en_avg = en_avg / N_mc ;
         en_square_avg = en_square_avg / N_mc ;
         
        double autosum(0); 
        for (unsigned int i=0; i < N_mc; ++i )
        {
        
               double sum(0);
               for (unsigned int j=1+i;j < N_mc -i; ++j )
               {
			sum += energy[j]*energy[j+i] - en_square_avg ;

                }
                autosum +=  sum/(N_mc - i);

         }
         
        double autocorr_time =  autosum / (en_square_avg - en_avg * en_avg);
        cout << autocorr_time;

       	ofstream fout("0.dat", std::fstream::app);	// Opens a file for output
        fout << autocorr_time << endl;
        fout.close();

	return 0;
}
