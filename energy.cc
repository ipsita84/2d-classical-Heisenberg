// g++ -std=c++11 -Wall -O3 energy.cc -o testo
// vim: set ai et cin ts=4 sw=4 tw=80:

//warming up system for first N_mc updates
//averaging energy for the next N_mc updates
//Metropolis algorithm employed
//Parameters that can be changed for different runs:
//J, axis1, axis2, N_mc

#include <iostream>
#include <iomanip>
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


const unsigned int axis1 = 10, axis2 = axis1;
// above assigns length along each dimension of the 2d configuration
const unsigned int sys_size = axis1 * axis2;

//No.of Monte Carlo updates we want
const unsigned int N_mc = 1e5;

const double beta=1;

//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_3d sitespin, array_2d J, std::array <double, 3> h);
double nn_energy(array_3d sitespin,  array_2d J, std::array <double, 3> h,
       unsigned int row, unsigned int col);

int main(int argc, char const * argv[])
{


    array_2d J(boost::extents[3][3]);
    double mx=0, my=0, mz=0;
    std::array <double, 3> h = {0,0,0};
    //std::array <double, N_mc> energy_array =  {0}, mx_array =  {0}, my_array =  {0}, mz_array =  {0};

    //Read the random signed bonds for a particular stored realization
    ifstream gin("J0.dat");
    ofstream f1out("mag0-test.dat",std::fstream::app);	// Opens a file for output
    ofstream fout("Energy0-test.dat", std::fstream::app);


    for (unsigned int comp1=0; comp1<3; ++comp1)
    {
        for (unsigned int comp2=0; comp2<3; ++comp2)
        {

            gin>>J[comp1][comp2];


        }
    }
    gin.close();

    //ofstream gout("test.dat");	// Opens a file for output

    array_3d sitespin(boost::extents[3][axis1][axis2]);

    // stores the spin configuration of the system
    //initial state chosen by random no. generator above
    for (unsigned int i = 0; i < axis1; ++i)
    {
        for (unsigned int j = 0; j < axis2; ++j)
        {
            double theta = roll_coin(0,pi);
            double phi = roll_coin(0,2*pi);
            sitespin[0][i][j] = sin(theta)*cos(phi);
            sitespin[1][i][j] = sin(theta)*sin(phi);
            sitespin[2][i][j] = cos(theta);
            //double s= pow(sitespin[0][i][j],2)+pow(sitespin[1][i][j],2) +pow(sitespin[2][i][j],2);
            //cout << s << endl ;
            double checksum= pow(sitespin[0][i][j],2)
                               +pow(sitespin[1][i][j],2)
                               +pow(sitespin[2][i][j],2);
             if (checksum > 1) {printf (" initial %f error \n", checksum);}

        }
    }
    double energy(0);
    double en_sum;
    unsigned int moves_accepted(0);

    for (unsigned int hsteps=0; hsteps<1; ++hsteps)
    {
        h[0] = 0 + hsteps*0.5;
        energy = energy_tot(sitespin, J, h);
        en_sum =0;
        std::array <double, N_mc> energy_array =  {0};
        std::array <double, N_mc> mx_array={0}, my_array ={0},mz_array ={0};
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
                double energy_minus_rnd_site =energy_old - nn_energy(sitespin,J,h, row, col);
                double r0 = 0.5*random_real(-1, 1)/beta;
                double r1 = 0.5*random_real(-1, 1)/beta;
                double r2 = 0.5*random_real(-1, 1)/beta;
                double tot = pow( s0+r0, 2)+pow( s1+ r1, 2)+pow(s2 + r2, 2);
                //printf ("tot %f \n",tot);

                sitespin[0][row][col] = (s0+r0)/sqrt(tot);
                sitespin[1][row][col] = (s1+r1)/sqrt(tot);
                sitespin[2][row][col]= (s2+r2)/sqrt(tot);

              double checksum= pow(sitespin[0][row][col],2)
                               +pow(sitespin[1][row][col],2)
                               +pow(sitespin[2][row][col],2);
             if (checksum > 1.00001) {printf ("%f error \n", checksum);}

                double energy_new = energy_minus_rnd_site +  nn_energy(sitespin,J, h,row, col);
                double energy_diff = energy_new - energy_old;
                double acc_ratio = exp(-1.0 * energy_diff* beta/200);
                double r =  random_real(0, 1) ;	//Generate a random no. r such that 0 < r < 1
                //Spin flipped if r <= acceptance ratio
                if (r <= acc_ratio)
                {
                    energy = energy_new ;
                    moves_accepted = moves_accepted +1;// cout << "energy changed" <<endl;
                }
                if (r > acc_ratio)
                {
                    sitespin[0][row][col] = s0;
                    sitespin[1][row][col] = s1;
                    sitespin[2][row][col] = s2;
                    energy = energy_old ;

                }
                //gout << i+j << '\t'  << energy << '\t' << moves_accepted << endl;
              //double checksum= pow(sitespin[0][row][col],2)
                               //+pow(sitespin[1][row][col],2)
                              // +pow(sitespin[2][row][col],2);
             //if (checksum > 1) {printf ("%f error \n", checksum);}
            }

            if (i > 1e5)
            {
                en_sum += energy;
                energy_array[i-N_mc -1] = energy;
                //double rat  = 1.0* moves_accepted/(i*sys_size);
                //gout<< energy << endl;
                for (unsigned int l = 0; l < axis1; ++l)
                {
                    for (unsigned int j = 0; j < axis2; ++j)
                    {
                        mx += sitespin[0][l][j] ;
                        mx_array[i-N_mc -1] += sitespin[0][l][j] ;
                        my += sitespin[1][l][j] ;
                        my_array[i-N_mc -1] += sitespin[1][l][j] ;
                        mz += sitespin[2][l][j] ;
                        mz_array[i-N_mc -1] += sitespin[2][l][j] ;
                    }
                }

            }
        }


        double sigma_en = 0, sigma_mx = 0, sigma_my = 0, sigma_mz = 0;
        for (unsigned i=0; i< N_mc; i++)
        {
            sigma_en += (energy_array[i] - en_sum/ N_mc) 
                        * (energy_array[i] - en_sum/ N_mc);
            sigma_mx += (mx_array[i] - mx/ N_mc) * (mx_array[i] - mx/ N_mc) ;
            sigma_my += (my_array[i] - my/ N_mc) * (my_array[i] - my/ N_mc) ;
            sigma_mz += (mz_array[i] - my/ N_mc) * (mz_array[i] - mz/ N_mc) ;
        }

        fout.setf( ios_base::fixed, ios_base::floatfield );
        fout.precision(2);
        fout << setw(6) << h[0];
        fout.precision(7);
        fout << setw(15)
             << en_sum / N_mc << setw(15)
             << sqrt(sigma_en) / N_mc << endl;
        // printing energy to file "Energy.dat"

        f1out.setf( ios_base::fixed, ios_base::floatfield );
        f1out.precision(2);
        f1out << setw(6) << h[0];
        f1out.precision(7);
        f1out << setw(15) << mx/(sys_size*N_mc)
              << setw(15) << sqrt(sigma_mx)/(sys_size*N_mc)
              << setw(15) << my/(sys_size*N_mc)
              << setw(15) << sqrt(sigma_my)/(sys_size*N_mc)
              << setw(15) << mz/(sys_size*N_mc)
              << setw(15) << sqrt(sigma_mz)/(sys_size*N_mc)  << endl;
        // printing magnetization to file "mag.dat"
        mx=0;
        my=0;
        mz=0;
    }

    //gout.close();
    fout.close();
    f1out.close();
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

double energy_tot(array_3d sitespin, array_2d J, std::array <double, 3> h)
{
    double energy = 0;

    for (unsigned comp  = 0; comp < 3; ++comp)

    {
        for (unsigned int i = 0; i < axis1 ; ++i)
        {
            for (unsigned int j = 0; j < axis2 ; ++j)
            {
                energy += -h[comp]*sitespin[comp][i][j];
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

//Calculating interaction & on-site energy change for spin
//at random site->(row,col) with its nearest neighbours
double nn_energy(array_3d sitespin,  array_2d J, std::array <double, 3> h,
                 unsigned int row, unsigned int col)
{
    double nn_en = 0;

    // Calculate the on-site energy
    for (unsigned comp  = 0; comp < 3; ++comp)

    {
        nn_en += -h[comp]*sitespin[comp][row][col];

    }

    for (unsigned comp1  = 0; comp1 < 3; ++comp1)
    {
        for (unsigned comp2  = 0; comp2 < 3; ++comp2)
        {
            if (row > 0 && row < axis1 - 1)
            {
                nn_en += J[comp1][comp2] * sitespin[comp2][row][col] * sitespin[comp1][row
                         -1][col];
                nn_en += J[comp1][comp2]*sitespin[comp1][row][col] * sitespin[comp2][row
                         +1][col];
            }

            if (col > 0 && col < axis2 - 1)
            {
                nn_en += J[comp1][comp2]*sitespin[comp2][row][col] * sitespin[comp1][row][col
                         -1];
                nn_en += J[comp1][comp2]*sitespin[comp1][row][col] * sitespin[comp2][row][col
                         +1];
            }

            if (row == 0)
            {
                nn_en += J[comp1][comp2]*sitespin[comp2][0][col] * sitespin[comp1][axis1
                         -1][col];
                nn_en += J[comp1][comp2]* sitespin[comp1][0][col] * sitespin[comp2][1][col];
            }

            if (row == axis1-1)
            {
                nn_en += J[comp1][comp2]* sitespin[comp2][axis1 - 1][col] *
                         sitespin[comp1][axis1-2][col];
                nn_en += J[comp1][comp2]* sitespin[comp1][axis1-1][col] *
                         sitespin[comp2][0][col];
            }

            if (col == 0)
            {
                nn_en += J[comp1][comp2]* sitespin[comp2][row][0] * sitespin[comp1][row][axis2
                         -1];
                nn_en += J[comp1][comp2]* sitespin[comp1][row][0] * sitespin[comp2][row][1];
            }

            if (col == axis2-1)
            {
                nn_en += J[comp1][comp2]* sitespin[comp2][row][axis2-1] *
                         sitespin[comp1][row][axis2-2];
                nn_en += J[comp1][comp2]* sitespin[comp1][row][axis2-1] *
                         sitespin[comp2][row][0];
            }
        }
    }
    return nn_en;
}

