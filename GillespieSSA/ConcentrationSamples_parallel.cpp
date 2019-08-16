#include <fstream>
#include <math.h>
#include <random>
#include <mpi.h>

const int Ns = 1e5;     // Number of SSA runs
const double W = 90.0;  // System size
const double ts = 1e3;  // Switching time

//------ RNG ------//
std::random_device rd;
int seed = rd();
std::mt19937 gen(seed);
std::uniform_real_distribution<double> dis(0.0, 1.0);

//------ PARAMETERS ------//
const double a0 = 0.5;      
const double a = 3.0;
const double K = 9.0;
const double n = 3.0;

//------ INPUT BOUNDS ------//
const double umin = 0.12;
const double umax = 0.19;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // the processor number
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); // number of processors

    //------ INPUT PARAMETER ------//
    double u = umin + (umax-umin)*rank/(nprocs-1); // input on the processor rank

    //------ OUTPUT FILE ------//
    std::string OUTFILE = std::string("samples/W=") + std::to_string(W) + std::string("_u=") + std::to_string(u) + std::string(".txt");

    //------ AVERAGING OVER SSA RUNS (ENSAMBLE AVERAGE) ------//
    double xhist[Ns] = {0};
    for(int s=0; s<Ns; s++){
        // initialize
        double x = 1.5*W; // arbitrary IC
        double t = 0;
        // run until uncorrelation with the IC (2 switching times)
        while(t < 2*ts){
            //------ PROPENSITY ------//
            double r = W*(a0 + a*pow((u+x/W),n)/(K + pow((u+x/W),n)));
            //------ TIME UPDATE ------//
            t = t + log(1/dis(gen))/(r+x);
            //------ STATE UPDATE ------//
            if(dis(gen) < r/(r+x)){x = x+1;}
            else x = x-1;
        }
        // measure end value
        std::ofstream fout;
        fout.open(OUTFILE, std::ios_base::app); // appends to file instead of deleting
        fout << x/W << '\n';
        fout.close();
    }

    MPI_Finalize();
}
