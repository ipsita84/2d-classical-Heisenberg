// g++ -std=c++11 -Wall -O3 energy.cc -o testo

//warming up system for first N_mc updates
//averaging energy for the next N_mc updates
//Metropolis algorithm employed
//Parameters that can be changed for different runs:
//J, axis1, axis2, N_mc

#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/lexical_cast.hpp>
#include <math.h> 
#include <array>

const double pi = acos(-1.0);

// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using namespace std;

typedef boost::multi_array < double, 3 > array_3d;
typedef boost::multi_array < double, 2 > array_2d;
// typedef keyword allows you to create an alias fo a data type


const unsigned int axis1 = 8, axis2 = axis1;
// above assigns length along each dimension of the 2d configuration
const unsigned int sys_size = axis1 * axis2;

//No.of Monte Carlo updates we want
const unsigned int N_mc = 1e5;

const double beta=1;


const std::array <double, 3> h={0,0,1};

//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_3d sitespin, array_2d J);
double nn_energy(array_3d sitespin,  array_2d J, unsigned int row, unsigned int col);

int main(int argc, char const * argv[])
{


        array_2d J(boost::extents[3][3]);
        
	//Read the random signed bonds for a particular stored realization
	ifstream gin("J.dat");
	
        for (unsigned int comp1=0; comp1<3; ++comp1)
        {
            for (unsigned int comp2=0; comp2<3; ++comp2)
            {  

			gin>>J[comp1][comp2];
			
		    
	     }
         }

	gin.close();




//	cout << "Enter beta" << endl;
//	cin >> beta;

	ofstream fout("Energy.dat");	// Opens a file for output
        ofstream gout("test.dat");	// Opens a file for output


		// Create a 3d array that is comp*axis1 * axis2
		array_3d sitespin(boost::extents[3][axis1][axis2]);

		// stores the spin configuration of the system
		//initial state chosen by random no. generator above
		for (unsigned int i = 0; i < axis1; ++i)
		{	for (unsigned int j = 0; j < axis2; ++j)
                        {       double theta = roll_coin(0,pi);
                                double phi = roll_coin(0,2*pi);
				sitespin[0][i][j] = sin(theta)*cos(phi);
				sitespin[1][i][j] = sin(theta)*sin(phi);
				sitespin[2][i][j] = cos(theta);
                         }
                }  

		double energy = energy_tot(sitespin, J);
		double en_sum(0);


		for (unsigned int i = 1; i <=1e5+N_mc; ++i)
		{
			for (unsigned int j = 1; j <= sys_size; ++j)
			{
				//Now choose a random spin site with site no.=label
				unsigned int label, row, col ;
				label = roll_coin(1, sys_size);

				if (label % axis2 == 0)
				{
					row = (label / axis2) - 1;
					col = axis2 -1 ;
				}
				else
				{
					col = label % axis2 - 1;
					row = (label-col-1)/axis2;
				}


                                double s0 = sitespin[0][row][col];
                                double s1 = sitespin[1][row][col];
                                double s2 = sitespin[2][row][col];
				double energy_old =energy ;
				double energy_minus_rnd_site =energy_old - nn_energy(sitespin,J, row, col);
		
				
				
				double r1 = 0.5*random_real(0, 1)/beta;
				double r2 = 0.5*random_real(0, 1)/beta;
				double r3 = 0.5*random_real(0, 1)/beta;
 
                                 
                                double tot = sqrt( pow( s0+ r1, 2)+pow( s1+ r1, 2)+pow(s2 + r1, 2) );
                                sitespin[0][row][col] = (s0+r1)/tot;
                                sitespin[1][row][col] = (s1+r2)/tot;
                                sitespin[2][row][col]= (s2+r3)/tot;
                                double energy_new = energy_minus_rnd_site +  nn_energy(sitespin,J, row, col);
                                double energy_diff = energy_new - energy_old;
				double acc_ratio = exp(-1.0 * energy_diff* beta);
 
				//Generate a random no. r such that 0 < r < 1
			        double r =  random_real(0, 1) ;	
				//Spin flipped if r <= acceptance ratio

				if (r <= acc_ratio)
				{
					energy = energy_new ;
				}

				if (r > acc_ratio)
				{
					sitespin[0][row][col] = s0;
                                        sitespin[1][row][col] = s1;
                                        sitespin[2][row][col] = s2;
					energy = energy_old ;

				}
			}

			if (i > 1e5){ en_sum += energy;
                                      gout << i<< '\t'  << energy  << endl;}
		}

		fout << beta
		     << '\t' << en_sum / N_mc << endl;

        gout.close();
	fout.close();
	return 0;
}

//function to generate random integer
// between 2 integers a & b, including a & b
int roll_coin(int a, int b)
{
	boost::random::uniform_int_distribution <> dist(a, b);
	return dist(gen);
}

//function to generate random real no.
// between 2 integers a & b, including a & excluding b

double random_real(int a, int b)
{
	boost::random::uniform_real_distribution <> dist(a, b);
	// uniform_real_distribution: continuous uniform distribution
	//on some range [min, max) of real number
	return dist(gen);
}

//function to calculate total energy
//for a given spin configuration
//with periodic boundary conditions

double energy_tot(array_3d sitespin, array_2d J)
{
	double energy = 0;
	
	for (unsigned comp  = 0; comp < 3; ++comp)

	 {      for (unsigned int i = 0; i < axis1 ; ++i)
	       {
		     for (unsigned int j = 0; j < axis2 ; ++j)
		     {
			energy += h[comp]*sitespin[comp][i][j];
		     }
	        }
             
         }

	for (unsigned comp1  = 0; comp1 < 3; ++comp1)
        {
	    for (unsigned comp2  = 0; comp2 < 3; ++comp2)
            {

	       for (unsigned int i = 0; i < axis1 - 1; ++i)
	       {
		     for (unsigned int j = 0; j < axis2 - 1; ++j)
		     {
			energy += J[comp1][comp2]*sitespin[comp1][i][j]*sitespin[comp2][i+1][j];
			energy += J[comp1][comp2]*sitespin[comp1][i][j]*sitespin[comp2][i][j+1];
		     }
	        }
             }
         }

	//periodic boundary conditions
	for (unsigned comp1  = 0; comp1 < 3; ++comp1)
        {
	    for (unsigned comp2  = 0; comp2 < 3; ++comp2)
            {
	       for (unsigned int j = 0; j < axis2; ++j)
		energy += J[comp1][comp2]*sitespin[comp1][axis1-1][j] * sitespin[comp2][0][j];
            }
         }


	for (unsigned comp1  = 0; comp1 < 3; ++comp1)
        {
	    for (unsigned comp2  = 0; comp2 < 3; ++comp2)
            {

	       for (unsigned int i = 0; i < axis1; ++i)
		    energy += J[comp1][comp2]*sitespin[comp1][i][axis2-1] * sitespin[comp2][i][0];
            }
         }

	return energy;
}

//Calculating interaction energy change for spin
//at random site->(row,col) with its nearest neighbours
double nn_energy(array_3d sitespin,  array_2d J, unsigned int row, unsigned int col)
{
	double nn_en = 0;
	for (unsigned comp1  = 0; comp1 < 3; ++comp1)
        {
	    for (unsigned comp2  = 0; comp2 < 3; ++comp2)
            {
	      if (row > 0 && row < axis1 - 1)
	      {
		nn_en += J[comp1][comp2] * sitespin[comp2][row][col] * sitespin[comp1][row-1][col];
		nn_en += J[comp1][comp2]*sitespin[comp1][row][col] * sitespin[comp2][row+1][col];
	      }

	      if (col > 0 && col < axis2 - 1)
	      {
		nn_en += J[comp1][comp2]*sitespin[comp2][row][col] * sitespin[comp1][row][col-1];
		nn_en += J[comp1][comp2]*sitespin[comp1][row][col] * sitespin[comp2][row][col+1];
	      }

	      if (row == 0)
	      {
		nn_en += J[comp1][comp2]*sitespin[comp2][0][col] * sitespin[comp1][axis1-1][col];
		nn_en += J[comp1][comp2]* sitespin[comp1][0][col] * sitespin[comp2][1][col];
	      }

	      if (row == axis1-1)
	      {
		nn_en += J[comp1][comp2]* sitespin[comp2][axis1 - 1][col] * sitespin[comp1][axis1-2][col];
	 	nn_en += J[comp1][comp2]* sitespin[comp1][axis1-1][col] * sitespin[comp2][0][col];
	      }

	      if (col == 0)
	      {
		nn_en += J[comp1][comp2]* sitespin[comp2][row][0] * sitespin[comp1][row][axis2-1];
		nn_en += J[comp1][comp2]* sitespin[comp1][row][0] * sitespin[comp2][row][1];
	      }

	      if (col == axis2-1)
	      {
		nn_en += J[comp1][comp2]* sitespin[comp2][row][axis2-1] * sitespin[comp1][row][axis2-2];
		nn_en += J[comp1][comp2]* sitespin[comp1][row][axis2-1] * sitespin[comp2][row][0];
	      }
           }
        }
	return nn_en;
}

