//1D harmonic oscillator with k=1, m=1
//the simulation is used Euler algorithm [euler_algo()]
//and Velocity-Verlet algorithm [VV_algo()]
//the data will been stored in csv files 
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;




void euler_algo(double t); //Euler algorithm

void VV_algo(double t); //Velocity-Verlet algorithm

int main()
{
    using namespace std;
    
    // euler_algo(0.1);
    euler_algo(0.01);
    // VV_algo(0.1);
    // VV_algo(0.01);

    return 0;
}

void euler_algo(double t) //t is the steps length
{
    long steps = 10000/t;
    double * x = new double [steps]; //x is the position of oscillator
    double * v = new double [steps]; //v is the velocity of oscillator

    ofstream outfile; //store the x v data to the files
    string filename = "Euler_algo_"+to_string(t) + ".csv";
    outfile.open(filename);
    
    outfile << "x,v"<<endl;

    x[0] = 0.0;
    v[0] = 1.0;

    outfile << x[0] << "," << v[0] << endl;

    for (uint i=0; i < steps-1; i++)
    {
        x[i+1] = x[i] + v[i]*t;
        v[i+1] = v[i] - x[i]*t; // F(i) = - x(i)
        outfile << x[i+1] << "," << v[i+1] << endl;
        
        // x[0] = x[1];
        // v[0] = v[1];
    }

    delete [] x;
    delete [] v;

    outfile.close();
}

void VV_algo(double t) //t is the steps length
{
    int steps = 10000/t;
    double x[2]; //x is the position of oscillator
    double v[2]; //v is the velocity of oscillator

    ofstream outfile; //store the x v data to the files
    string filename = "VV_algo_"+to_string(t)+".csv";
    outfile.open(filename);
    
    outfile << "x,v" << endl;

    x[0] = 0.0;
    v[0] = 1.0;

    outfile << x[0] << "," << v[0] << endl;

    for (uint i=0; i < steps-1; i++)
    {
        x[1] = x[0] + v[0]*t - 0.5*x[0]*pow(t,2);// F(i) = - x(i)
        v[1] = v[0] - 0.5*(x[0]+x[1])*t; //F(i+1) = -x(i+1)
        
        outfile << x[1] << "," << v[1] << endl;

        x[0] = x[1];
        v[0] = v[1];
    }

    outfile.close();
}