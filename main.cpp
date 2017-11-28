#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <string>
#include "bond.h"
#include "site.h"
#include "lattice.h"
#include "montecarlo.h"
#include "gaussiandeviate.h"

using namespace std;
using std::ofstream; using std::string;

void one_run(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, long int latticeseed, long int MCseed, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void one_run(int L1, int L2, int L3, int eqsteps, int mcsteps_inbin, int no_of_bins, long int latticeseed, long int MCseed, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void run_for_several_betas(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, long int latticeseed, long int MCseed, double betamin, double betamax, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool dobootstrap, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void run_for_betasgiven(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, long int latticeseed, long int MCseed, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool dobootstrap, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void run_for_betasgiven2(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, long int latticeseed, long int MCseed, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep,  bool dobootstrap, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void run_for_betasgiven_diffdirs(int L1, int L2, int L3, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, long int latticeseed, long int MCseed, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep,  bool dobootstrap, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);

// Lattice functions
// fcc or flexible
void extract_yline(int L, string ylinefilenamePrefix); // Extracting the lattice sites along (0,y,0) in the fcc
void extract_yline(int L1, int L2, int L3, string ylinefilenamePrefix);
void extract_xline(int L, string ylinefilenamePrefix); // fcc only
void extract_zline(int L, string ylinefilenamePrefix); // fcc only
void extract_qyline(int L, string ylinefilenamePrefix); // Extracting the lattice sites along (0,K,0) in the fcc
void extract_qyline(int L1, int L2, int L3, string latticefilenamePrefix); // Extension
void extract_qxline(int L, string ylinefilenamePrefix); // fcc only
void extract_qzline(int L, string ylinefilenamePrefix); // fcc only
void extract_qdline(int L, string ylinefilenamePrefix); // fcc only
void extract_yline_shifted(int L, double xshift, double zshift, string ylinefilenamePrefix); // Extracting the lattice sites along (0,y,0) in the fcc
void diagline(int L, string latticefilenamePrefix);
void lattice_coordinates_straightforward(int L, char type_lattice, string latticefilenamePrefix);
void lattice_coordinates_straightforward(int L1, int L2, int L3, char type_lattice, string latticefilenamePrefix);
void lattice_coordinates_xyz_lines(int L, string latticefilenamePrefix);  // Only for fcc
void lattice_coordinates_xyz_lines(int L1, int L2, int L3, string latticefilenamePrefix);
void reciprocallattice_coordinates(int L, char type_lattice, string latticefilenamePrefix);
void reciprocallattice_coordinates(int L1, int L2, int L3, char type_lattice, string latticefilenamePrefix);
void reciprocallattice_coordinates_xyzline(int L, char type_lattice, string latticefilenamePrefix); // Only for fcc
// cubic
void cubic_extract_xyzlines(int L, string latticefilenamePrefix);
// quadratic
void quadratic_extract_xylines(int L, string latticefilenamePrefix);
void printyneighbours_fcc(int L, string latticefilenamePrefix);
void printnearestneighbours_fcc(int L, string latticefilenamePrefix);

// Test functions
void test_betagenerator(int beta_n, int betamin, int betamax);
void test_fftw(int L, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_fftw_againstsims(int L, int eqsteps, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, char type_lattice, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_fftw_againstsims_av(int L, int eqsteps, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, char type_lattice, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_fcc_extended(int L, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_fcc_extended_diffdims(int L1, int L2, int L3, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void test_dm();
void test_fcc_extended_yopen();
void test_fcc_extended_yopen_throughMC(int L, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void checkneighbours(int L, char type_lattice, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
void testnestnearestneighbour_chain(int L);
double dm_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double Dx, double Dy, double Dz);
double J_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double J);
void printrandomnumbers();

double converttemp(double intemp);


int main()
{
    bool DEBUG = true;
    if(DEBUG)    cout << "In main" << endl;

    // Input parameters
    int L = 8; // The program is going to be slow if we run for many particles on a 3D lattice

    int L1 = 6;
    int L2 = 8;
    int L3 = 6;

    // bools to determine system type
    bool isotropic    = true;
    bool sianisotropy = true;  // This one does not change its energy unless Dix, Diy and Diz are not all equal.
    bool magfield     = false;
    bool dm           = false;  // dmdiffdirs = false
    bool nextnearest  = true;

    // Bool to determine periodicity
    bool periodic     = true; // To determine whether we have periodic boundary conditions or not
                              // Only implemented for the chain so far

    // Selecting the lattice type
    // F: face-centered cubic (fcc); E: fcc with different directions; Y: fcc with open BCs in the y-dir
    // C: cubic; D: cubic with different directions
    // Q:quadratic; R: quadratic with different directions
    // O: chain;
    char type_lattice = 'Y';
    // If periodic is false, that means we get a grid with open boundary conditions. Currently,
    // that is only implemented for the chain.

    // Setting the strengths of the Hamiltonian terms
    // Magnetic field terms
    double hx = 1;    double hy = 2;    double hz = 7;
    // Single-ion anisotropy terms
    double Dix = 0.34;    double Diy = 1.82;    double Diz = 0; // Jensen's parameters
    //double Dix = 0.34;    double Diy = 1.92;    double Diz = 0; // Li's parameters
    //double Dix = 0.413;    double Diy = 1.42;    double Diz = 0; // Toft-Pettersen's LiNiPO4 parameters
    // DM terms
    //double Dx = 1.82;     double Dy = 0;    double Dz = 0;
    double Dx = 0.0;     double Dy = 0.32;    double Dz = 0;       // Toft-Pettersen's LiNiPO4 parameters

    // Heisenberg terms
    // The nearest neighbour coupling for O, F, C, Q
    double J = 2.28; // Nearest neighbour coupling
    //double J = 1;

    // These are the nearest neighbour couplings for: D, R
    // These are the next-nearest neighbour couplings for: O, F, Y
    //double Jx  = 0.67;    double Jy  = 0;    double Jz  = 0; // Simplified couplings, chain
    double Jx  = 0;    double Jy  = 0.67;    double Jz  = -0.05; // Jensen's parameters
    //double Jx  = 0;    double Jy  = 0.59;    double Jz  = -0.11; // Li's parameters
    //double Jx  = 0;    double Jy  = 0.67;    double Jz  = -0.06; // Toft-Pettersen's LiNiPO4 parameters


    // These are the nearest neighbour couplings for E, Y
    double Jxy = 0.30;    double Jxz = -0.11;    double Jyz = 1.04; // Jensen's parameters
    //double Jxy = 0.26;    double Jxz = -0.16;    double Jyz = 0.94; // Li's parameters
    //double Jxy = 0.32;    double Jxz = -0.11;    double Jyz = 1.00; // Toft-Pettersen's LiNiPO4 parameters

    vector<double> sitestrengthsin = vector<double>(6);
    sitestrengthsin[0] = hx;    sitestrengthsin[1] = hy;    sitestrengthsin[2] = hz;
    sitestrengthsin[3] = Dix;
    sitestrengthsin[4] = Diy;   sitestrengthsin[5] = Diz;

    vector<double> heisenbergin = vector<double>(7);
    heisenbergin[0] = J;
    heisenbergin[1] = Jx;       heisenbergin[2] = Jy;       heisenbergin[3] = Jz;
    heisenbergin[4] = Jxy;      heisenbergin[5] = Jxz;      heisenbergin[6] = Jyz;

    vector<double> dm_in = vector<double>(3);
    dm_in[0] = Dx;              dm_in[1] = Dy;              dm_in[2] = Dz;

    // Bools to determine printing
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = true; // This is set to false if we run for several betas
    bool dobootstrap = true; // Often set in the function

    // Run parameters
    int eqsteps = 10000;  // Number of steps in the equilibration procedure
    int mcsteps_inbin = 10000;  // MCsteps per bin.
    int no_of_bins = 100;     // The number of bins.

    // Could also convert from T to beta // That is probably easier
    double T = 30;
    double beta = converttemp(T); // beta = 11.6045221/T;

    // Input parameters specifically for run_for_several_betas
    int beta_n = 25;
    double betamin = 0.78;
    double betamax = 0.8;
    int betanset = 5;
    vector<double> betas = vector<double>(betanset);
    betas[0] = 0.5; betas[1] = 1.0; betas[2] = 2.0; betas[3] = 10.0; betas[4] = 50.0;

    // By default, the run_for_several_betas-functions do not calculate the correlation function
    //--------------------------------Running for several betas--------------------------------------//
    // A few examples on how to run
    long int latticeseed, MCseed;
    string filenamePrefix;
    /*// Running with several different random seeds
    L = 6;
    beta_n = 5;
    //betamin = 1.1; betamax = 2.0;
    mcsteps_inbin = 100000;
    eqsteps = 100000;
    betamin = 0.808;    betamax = 0.818;
    //mcsteps_inbin = 10000;
    
    latticeseed = 71; MCseed = 31;
    // Jensen
    filenamePrefix = "fcc6x6x6yopen_beta0p808to0p818_Nbeta5_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_seed31_latticeseed71_slowcool";

    // Li
    //filenamePrefix = "fcc6x6x6yopen_beta0p01to1p4_Nbeta100_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc10000_bins100_seed31_latticeseed71_slowcool";

    // Toft-Pettersen
    //filenamePrefix = "fcc6x6x6yopen_beta0p01to1p4_Nbeta100_nnJyz1p0_nnJxy0p32_nnJxzm0p11_nnnJy0p67_nnnJzm0p06_sianDx0p413_Dy1p42_DMDy0p32_eq10000_mc10000_bins100_seed31_latticeseed71_slowcool";

    run_for_several_betas(L, eqsteps, mcsteps_inbin, no_of_bins, beta_n, latticeseed, MCseed, betamin, betamax, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);


    latticeseed = 47; MCseed = 83;
    // Jensen
    filenamePrefix = "fcc6x6x6yopen_beta0p808to0p818_Nbeta5_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_seed83_latticeseed47_slowcool";

    // Li
    //filenamePrefix = "fcc6x6x6yopen_beta0p01to1p4_Nbeta100_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc10000_bins100_seed83_latticeseed47_slowcool";

    // Toft-Pettersen
    //filenamePrefix = "fcc6x6x6yopen_beta0p01to1p4_Nbeta100_nnJyz1p0_nnJxy0p32_nnJxzm0p11_nnnJy0p67_nnnJzm0p06_sianDx0p413_Dy1p42_DMDy0p32_eq10000_mc10000_bins100_seed83_latticeseed47_slowcool";   

    //run_for_several_betas(L, eqsteps, mcsteps_inbin, no_of_bins, beta_n, latticeseed, MCseed, betamin, betamax, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);


    latticeseed = 23; MCseed = 59;
    // Jensen
    filenamePrefix = "fcc6x6x6yopen_beta0p808to0p818_Nbeta5_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_seed59_latticeseed23_slowcool";
    
    // Li
    //filenamePrefix = "fcc6x6x6yopen_beta0p01to1p4_Nbeta100_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc10000_bins100_seed59_latticeseed23_slowcool";

    // Toft-Pettersen
    //filenamePrefix = "fcc6x6x6yopen_beta0p01to1p4_Nbeta100_nnJyz1p0_nnJxy0p32_nnJxzm0p11_nnnJy0p67_nnnJzm0p06_sianDx0p413_Dy1p42_DMDy0p32_eq10000_mc10000_bins100_seed59_latticeseed23_slowcool";

    //run_for_several_betas(L, eqsteps, mcsteps_inbin, no_of_bins, beta_n, latticeseed, MCseed, betamin, betamax, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    */

    // We can also set a few values of beta manually and run for them
    //run_for_betasgiven(L, eqsteps, mcsteps_inbin, no_of_bins, beta_n, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep,  dobootstrap, type_lattice, filenamePrefix, betas, sitestrengthsin, heisenbergin, dm_in);

    //run_for_betasgiven2(L, eqsteps, mcsteps_inbin, no_of_bins, betanset, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep, bool dobootstrap, type_lattice, filenamePrefix, betas, sitestrengthsin, heisenbergin, dm_in);
    //run_for_betasgiven_diffdirs(L1, L2, L3, eqsteps, mcsteps_inbin, no_of_bins, betanset, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, periodic, printeveryMCstep,  dobootstrap, type_lattice, filenamePrefix, betas, sitestrengthsin, heisenbergin, dm_in);

    //----------------------------------Running for one beta-----------------------------------------//
    // 
    L1 = 20;
    L2 = 20;
    L3 = 10;
    L = 6;
    eqsteps = 10000;
    mcsteps_inbin = 10000;
    no_of_bins = 100;
    latticeseed = 71;
    MCseed = 31;
    T = 5.0;
    beta = 11.6045221/T;
    type_lattice= 'Y';

    //sianisotropy = false;
    calculatespincorrelationfunction = true; // Can comment this out just to make it faster
    beta = converttemp(11);
    filenamePrefix = "chain100p_Jnn2p28_Jnnn0p67_0p5xsian_11K_eq100000_mc100000_bins100_seed31_latticeseed71_slowcool_spinconfigtrue_with_spcorrf_221117";
    //one_run(L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, beta, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    
    beta = converttemp(10);
    filenamePrefix = "littletest";
    //filenamePrefix = "fcc20x20x20yopen_Jensenlike_nnnJy0p91_T10K_eq100000_mc100000_bins100_seed31_latticeseed71_slowcool_spinconfigtrue_with_spcorrf_lastminute_II_right_fixed";
    one_run(L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, beta, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    

    
    //-------------------------------------Test functions--------------------------------------------//
    //L = 2;
    //test_fftw(L, sitestrengthsin, heisenbergin, dm_in);
    //test_fftw_againstsims(L, eqsteps, beta, isotropic, sianisotropy, magfield, dm, periodic, type_lattice, sitestrengthsin, heisenbergin, dm_in);
    //test_fftw_againstsims_av(L, eqsteps, beta, isotropic, sianisotropy, magfield, dm, periodic, type_lattice, sitestrengthsin, heisenbergin, dm_in);
    //test_fcc_extended(L, isotropic, sianisotropy, magfield, dm, periodic, sitestrengthsin, heisenbergin, dm_in);
    //test_fcc_extended_diffdims(L1, L2, L3, isotropic, sianisotropy, magfield, dm, periodic, sitestrengthsin, heisenbergin, dm_in);
    //test_dm();
    //L = 4;
    //test_fcc_extended_yopen();
    //test_fcc_extended_yopen_throughMC(L, sitestrengthsin, heisenbergin, dm_in);
    //checkneighbours(L, type_lattice, periodic, sitestrengthsin, heisenbergin, dm_in);
    //testnestnearestneighbour_chain(L);

    //---------------------------Finding coordinates, lines, etc------------------------------------//
    // We only need to find the line (0,K,0) once for each fcc L1xL2xL3 Lattice
    //L = 12;
    L1 = 20;
    L2 = 20;
    L3 = 20;
    type_lattice = 'E';
    double xshift = 2;
    double zshift = 1;
    string latticefilenamePrefix = "fcc20x20x20";
    string latticefilenamePrefixd = "L6";
    //string latticefilenamePrefix = "test";
    // Line finding functions
    // fcc only:
    //extract_yline(L, latticefilenamePrefix);
    //extract_yline(L1, L2, L3, latticefilenamePrefix);
    //extract_xline(L, latticefilenamePrefix);
    //extract_zline(L, latticefilenamePrefix);
    extract_qyline(L, latticefilenamePrefix);
    //extract_qyline(L1, L2, L3, latticefilenamePrefix);
    //extract_qxline(L, latticefilenamePrefix);
    //extract_qzline(L, latticefilenamePrefix);
    //extract_qdline(L, latticefilenamePrefix);
    //extract_yline_shifted(L, xshift, zshift, latticefilenamePrefix);
    //diagline(L, latticefilenamePrefixd);
    //printyneighbours_fcc(L, latticefilenamePrefix);
    //printnearestneighbours_fcc(L, latticefilenamePrefix);
    // Lattice coordinates
    //lattice_coordinates_straightforward(L, type_lattice, latticefilenamePrefix);
    //lattice_coordinates_straightforward(L1, L2, L3, type_lattice, latticefilenamePrefix);
    //reciprocallattice_coordinates(L, type_lattice, latticefilenamePrefix);
    //reciprocallattice_coordinates(L1, L2, L3, type_lattice, latticefilenamePrefix);
    //lattice_coordinates_xyz_lines(L, latticefilenamePrefix);
    //lattice_coordinates_xyz_lines(L1, L2, L3, latticefilenamePrefix);
    //reciprocallattice_coordinates_xyzline(L, type_lattice, latticefilenamePrefix);
    //isitonplane(filenamePrefix);
    //isitonline(filenamePrefix);
    string indexfilenamePrefix = "fcc6x6x6_index126";
    int maxindex = 126;
    //findlinethroughmax(L, maxindex, indexfilenamePrefix);
    //cubic_extract_xyzlines(L, latticefilenamePrefix);
    //quadratic_extract_xylines(L, latticefilenamePrefix);

    //printrandomnumbers();

}

void one_run(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, long int latticeseed, long int MCseed, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
     bool dobootstrap = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.debugmode(true);
    // Run Metropolis algorithm with beta
    mymc.runmetropolis(beta);
    mymc.endsims();
    mymc.test_couplings_strengths(); // Just to test it. Can remove later
    cout << "beta = " << beta << endl;

    double T = 11.6045221/beta;
    cout << "T = " << T << endl;
    //cout << "L = " << L << endl;
}

void one_run(int L1, int L2, int L3, int eqsteps, int mcsteps_inbin, int no_of_bins,long int latticeseed, long int MCseed, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool calculatespincorrelationfunction, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
     bool dobootstrap = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L1, L2, L3, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.debugmode(true);
    // Run Metropolis algorithm
    mymc.runmetropolis(beta);
    mymc.endsims();
    mymc.test_couplings_strengths(); // Just to test it. Can remove later
    cout << "beta = " << beta << endl;
    //cout << "L = " << L << endl;
}

void run_for_several_betas(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, long int latticeseed, long int MCseed, double betamin, double betamax, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool dobootstrap, char type_lattice, string filenamePrefix, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    bool calculatespincorrelationfunction = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    //mymc.debugmode(true); // Debugging options
    //mymc.majordebugtrue();

    vector<double> betas = vector<double>(beta_n);
    double deltabeta = (betamax-betamin)/(beta_n-1.0);

    for(int i=0; i<beta_n; i++)
    {
        betas[i] = betamin + deltabeta*i;
        mymc.runmetropolis(betas[i]);
        mymc.reset_energy();
        cout << "Temp " << i << " done" << endl;
    }
    mymc.endsims();
    mymc.test_couplings_strengths(); // Just to test it. Can remove later
}

void run_for_betasgiven(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, long int latticeseed, long int MCseed, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool dobootstrap, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    bool calculatespincorrelationfunction = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    //mymc.debugmode(true);
    //mymc.majordebugtrue();

    for(int i=0; i<beta_n; i++)
    {
        mymc.runmetropolis(betas[i]);
        mymc.reset_energy();
        cout << betas[i] << " done" << endl;
    }
    mymc.endsims();
    mymc.test_couplings_strengths(); // Just to test it. Can remove later
}

void run_for_betasgiven2(int L, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, long int latticeseed, long int MCseed, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool dobootstrap, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    bool calculatespincorrelationfunction = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    //mymc.debugmode(true);
    //mymc.majordebugtrue();

    for(int i=0; i<beta_n; i++)
    {
        mymc.runmetropolis(betas[i]);
        mymc.reset_energy();
    }
    mymc.endsims();
    mymc.test_couplings_strengths(); // Just to test it. Can remove later
}

void run_for_betasgiven_diffdirs(int L1, int L2, int L3, int eqsteps, int mcsteps_inbin, int no_of_bins, int beta_n, long int latticeseed, long int MCseed, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, bool printeveryMCstep, bool dobootstrap, char type_lattice, string filenamePrefix, vector<double> betas, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    bool calculatespincorrelationfunction = false;
    // Initializing Monte Carlo
    MonteCarlo mymc(L1, L2, L3, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    //mymc.debugmode(true);
    //mymc.majordebugtrue();

    for(int i=0; i<beta_n; i++)
    {
        mymc.runmetropolis(betas[i]);
        mymc.reset_energy();
    }
    mymc.endsims();
}

void extract_yline(int L, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (0,L,0)
    long int latticeseed = 59;
    Lattice mylattice = Lattice(L, latticeseed, false, false, false, false); // We are only interested in the sites
    vector<int> ysites = mylattice.fccyline();

    ofstream ylineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_yline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    ylineFile.open(filename);
    delete filename;

    int N = ysites.size();

    for(int i=0; i<N; i++)
    {
        ylineFile << ysites[i]  << endl;
    }
    ylineFile.close();
}

void extract_yline(int L1, int L2, int L3, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (0,L,0)
    long int latticeseed = 59;
    Lattice mylattice = Lattice(L1, L2, L3, latticeseed, false, false, false, false); // We are only interested in the sites
    vector<int> ysites = mylattice.fccyline();

    ofstream ylineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_yline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    ylineFile.open(filename);
    delete filename;

    int N = ysites.size();

    for(int i=0; i<N; i++)
    {
        ylineFile << ysites[i]  << endl;
    }
    ylineFile.close();
}

void extract_xline(int L, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (0,L,0)
    long int latticeseed = 59;
    Lattice mylattice = Lattice(L, latticeseed, false, false, false, false); // We are only interested in the sites
    vector<int> xsites = mylattice.fccxline();

    ofstream xlineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_xline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xlineFile.open(filename);
    delete filename;

    int N = xsites.size();

    for(int i=0; i<N; i++)
    {
        xlineFile << xsites[i]  << endl;
    }
    xlineFile.close();
}

void extract_zline(int L, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (0,L,0)
    long int latticeseed = 59;
    Lattice mylattice = Lattice(L, latticeseed, false, false, false, false); // We are only interested in the sites
    vector<int> zsites = mylattice.fcczline();

    ofstream zlineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_zline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    zlineFile.open(filename);
    delete filename;

    int N = zsites.size();

    for(int i=0; i<N; i++)
    {
        zlineFile << zsites[i]  << endl;
    }
    zlineFile.close();
}

void extract_yline_shifted(int L, double xshift, double zshift, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (0,L,0)
    long int latticeseed = 59;
    Lattice mylattice = Lattice(L, latticeseed, false, false, false, false); // We are only interested in the sites
    mylattice.fcc_helical_initialize_extended();
    vector<int> ysites = mylattice.fccyline_shifted(xshift, zshift);

    ofstream ylineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_yline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    ylineFile.open(filename);
    delete filename;

    int N = ysites.size();

    for(int i=0; i<N; i++)
    {
        ylineFile << ysites[i]  << endl;
    }
    ylineFile.close();
}

void extract_qyline(int L, string ylinefilenamePrefix)
{
    cout << "In extract_qyline" << endl;
    // Getting the Lattice sites along (0,K,0)
    long int latticeseed = 59;
    Lattice mylattice = Lattice(L, latticeseed, false, false, false, false); // We are only interested in the sites
    mylattice.fcc_helical_initialize_extended();
    vector<int> ysites = mylattice.fccqyline();

    ofstream ylineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_qyline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    ylineFile.open(filename);
    delete filename;

    int N = ysites.size();

    for(int i=0; i<N; i++)
    {
        ylineFile << ysites[i]  << endl;
    }
    ylineFile.close();
}

void extract_qyline(int L1, int L2, int L3, string latticefilenamePrefix)
{
    cout << "In extract_qyline" << endl;
    // Getting the Lattice sites along (0,K,0)
    long int latticeseed = 59;
    Lattice mylattice = Lattice(L1, L2, L3, latticeseed, false, false, false, false); // We are only interested in the sites
    mylattice.fcc_helical_initialize_extended();
    vector<int> ysites = mylattice.fccqyline();

    ofstream ylineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_qyline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    ylineFile.open(filename);
    delete filename;

    int N = ysites.size();

    for(int i=0; i<N; i++)
    {
        ylineFile << ysites[i]  << endl;
    }
    ylineFile.close();
}

void extract_qxline(int L, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (qx,0,0)
    long int latticeseed = 59;
    Lattice mylattice = Lattice(L, latticeseed, false, false, false, false); // We are only interested in the sites
    mylattice.fcc_helical_initialize_extended();
    vector<int> xsites = mylattice.fccqxline();

    ofstream xlineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_qxline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xlineFile.open(filename);
    delete filename;

    int N = xsites.size();

    for(int i=0; i<N; i++)
    {
        xlineFile << xsites[i]  << endl;
    }
    xlineFile.close();
}

void extract_qzline(int L, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along (0,0,qz)
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, false, false, false, false); // We are only interested in the sites
    mylattice.fcc_helical_initialize_extended();
    vector<int> zsites = mylattice.fccqzline();

    ofstream zlineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_qzline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    zlineFile.open(filename);
    delete filename;

    int N = zsites.size();

    for(int i=0; i<N; i++)
    {
        zlineFile << zsites[i]  << endl;
    }
    zlineFile.close();
}

void extract_qdline(int L, string ylinefilenamePrefix)
{
    // Getting the Lattice sites along a diagonal line
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, false, false, false, false); // We are only interested in the sites
    mylattice.fcc_helical_initialize_extended();
    vector<int> dsites = mylattice.fccqdline();

    ofstream dlineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_qdline.txt", ylinefilenamePrefix.c_str() );   // Create filename with prefix and ending
    dlineFile.open(filename);
    delete filename;

    int N = dsites.size();

    for(int i=0; i<N; i++)
    {
        dlineFile << dsites[i]  << endl;
    }
    dlineFile.close();
}

void diagline(int L, string latticefilenamePrefix)
{
    long int seed = 59;
    // Getting the Lattice sites along the diagonal in the quadratic lattice
    Lattice qlattice = Lattice(L, seed, false, false, false, false); // We are only interested in the sites
    qlattice.quadratic_helical_initialize();
    vector<int> dqsites = qlattice.diagline_quadr();

    // Getting the Lattice sites along the diagonal in the simple cubic lattice
    Lattice clattice = Lattice(L, seed, false, false, false, false); // We are only interested in the sites
    clattice.cubic_helical_initialize();
    vector<int> dcsites = clattice.diagline_cubic();

    // Getting the Lattice sites along the diagonal in the face centered cubic lattice
    Lattice flattice = Lattice(L, seed, false, false, false, false); // We are only interested in the sites
    flattice.fcc_helical_initialize_extended();
    vector<int> dfsites = flattice.diagline_cubic();

    ofstream dqlineFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_quadr_diagline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    dqlineFile.open(filename);
    delete filename;

    ofstream dclineFile;
    char *filename2 = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename2, "%s_cubic_diagline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    dclineFile.open(filename2);
    delete filename2;

    ofstream dflineFile;
    char *filename3 = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename3, "%s_fcc_diagline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    dflineFile.open(filename3);
    delete filename3;

    int N = dqsites.size();

    for(int i=0; i<N; i++)        dqlineFile << dqsites[i]  << endl;
    for(int i=0; i<N; i++)        dclineFile << dcsites[i]  << endl;
    for(int i=0; i<N; i++)        dflineFile << dfsites[i]  << endl;

    dqlineFile.close();
    dclineFile.close();
    dflineFile.close();
}

void cubic_extract_xyzlines(int L, string latticefilenamePrefix)
{
    // Getting the Lattice sites on the form of (x,0,0), (0,y,0), (0,0,z)
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, false, false, false, false); // We are only interested in the sites
    mylattice.cubic_helical_initialize();
    vector<int> xsites = mylattice.cubicxline();
    vector<int> ysites = mylattice.cubicyline();
    vector<int> zsites = mylattice.cubiczline();

    ofstream xlineFile;
    char *filenamex = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamex, "%s_xline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xlineFile.open(filenamex);
    delete filenamex;

    ofstream ylineFile;
    char *filenamey = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamey, "%s_yline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    ylineFile.open(filenamey);
    delete filenamey;

    ofstream zlineFile;
    char *filenamez = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamez, "%s_zline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    zlineFile.open(filenamez);
    delete filenamez;

    cout << "All the files are created" << endl;

    int Nx = xsites.size(); // Just in case we aren't considering a quadratic lattice

    cout << "And the dim size x extracted" << endl;

    for(int i=0; i<Nx; i++)
    {
        cout << "In loop over x-dim. i = " << i << endl;
        xlineFile << xsites[i] << endl;
    }

    int Ny = ysites.size();

    for(int i=0; i<Ny; i++)
    {
        cout << "In loop over y-dim. i = " << i << endl;
        ylineFile << ysites[i] << endl;
    }

    int Nz = zsites.size();

    for(int i=0; i<Nz; i++)
    {
        cout << "In loop over z-dim. i = " << i << endl;
        zlineFile << zsites[i] << endl;
    }
    xlineFile.close();
    ylineFile.close();
    zlineFile.close();
}

void quadratic_extract_xylines(int L, string latticefilenamePrefix)
{
    // Getting the Lattice sites on the form of (x,0,0), (0,y,0)
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, false, false, false, false); // We are only interested in the sites
    mylattice.quadratic_helical_initialize();
    vector<int> xsites = mylattice.quadrxline();
    vector<int> ysites = mylattice.quadryline();

    ofstream xlineFile;
    char *filenamex = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamex, "%s_xline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xlineFile.open(filenamex);
    delete filenamex;

    ofstream ylineFile;
    char *filenamey = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamey, "%s_yline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    ylineFile.open(filenamey);
    delete filenamey;

    cout << "All the files are created" << endl;

    // The sizes of the lines
    int Nx = xsites.size(); // Just in case we aren't considering a quadratic lattice
    int Ny = ysites.size();

    // Printing indices to file
    for(int i=0; i<Nx; i++)        xlineFile << xsites[i] << endl;
    for(int i=0; i<Ny; i++)        ylineFile << ysites[i] << endl;

    xlineFile.close();
    ylineFile.close();
}


void test_betagenerator(int beta_n,  int betamin, int betamax)
{
    vector<double> betas = vector<double>(beta_n);
    double deltabeta = (betamax-betamin)/(beta_n-1.0);

    cout << "betamin = " << betamin << "; betamax = " << betamax << "; number of values: " << beta_n << endl;
    cout << "deltabeta = " << deltabeta << endl;

    for(int i=0; i<beta_n; i++)
    {
        betas[i] = betamin + deltabeta*i;
        cout << "i = " << i << "; betas[i] = " << betas[i] << endl;
    }
}

void lattice_coordinates_straightforward(int L, char type_lattice, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, isotropic, sianisotropy, magfield, dm);

    // Type of Lattice
    if(type_lattice=='F')           mylattice.fcc_helical_initialize();          // F for fcc
    else if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended(); // E for extended
    else if(type_lattice=='Y')      mylattice.fcc_helical_initialize_extended_yopen(); // Y for y-dir
    else if(type_lattice=='C')      mylattice.cubic_helical_initialize();        // C for cubic
    else if(type_lattice=='Q')      mylattice.quadratic_helical_initialize();    // Q for quadratic
    else if(type_lattice=='O')      mylattice.chain_periodic_initialize();       // O for one-dimensional

    // Printing information about the q-vectors straight away
    ofstream xyzFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_xyz.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFile.open(filename);
    delete filename;

    xyzFile << "xyz-positions are listed units of grid lengths" << endl;

    cout << "File for printing xyz-coord to file is initiated" << endl;

    double x;
    double y;
    double z;

    int N  = mylattice.N;

    if(mylattice.dim==1) // Not really neccessary. Not double tested, either
    {
        for(int i=0; i<N; i++)    xyzFile << "Running index: " << i << "; Index : i: " << i << "; Position : x = " << i << endl;
    }
    if(mylattice.dim==2)
    {

        cout << "For quadratic, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];
            vector<double> pos = mylattice.sitepositions[i];
            x = ns[0];
            y = ns[1];
            xyzFile << "Running index: " << i << "; Index : i = " << x << ";   j = " << y << ";   Position : x = ";
            x = pos[0];
            y = pos[1];
            xyzFile << x << ";   y = " << y << endl;
        }
        cout << "Done printing to xyzFile" << endl;

    }

    else if(mylattice.dim==3)
    {
        cout << "For cubic or fcc, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns     = mylattice.sitecoordinates[i];
            vector<double> pos = mylattice.sitepositions[i];
            //cout << "ns retrieved" << endl;
            // These must be changed if we change into possibly setting L1, L2, L3 different
            // Don't really need b1, b2, b3, could just use a1, a2, a3 multiplied by 2*M_PI...
            // Could be more general, but we don't need to play around with our lattices that much...
            x = ns[0];
            y = ns[1];
            z = ns[2];
            // Print to file. Site number, qx, qy, qz.
            xyzFile << "Running index: " << i << ";  Index :   i = " << x << ";   j = " << y <<  ";   k = " << z <<"; Position :   x = ";
            x = pos[0];
            y = pos[1];
            z = pos[2];
            xyzFile << x << ";   y = " << y << ";   z = " << z << endl;
        }
        //cout << "Done printing to qFile" << endl;
    }
    xyzFile.close();
}

void lattice_coordinates_straightforward(int L1, int L2, int L3, char type_lattice, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;
    long int seed = 59;
    Lattice mylattice = Lattice(L1, L2, L3, seed, isotropic, sianisotropy, magfield, dm);

    // Type of Lattice
    if(type_lattice=='F')           mylattice.fcc_helical_initialize();          // F for fcc
    else if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended(); // E for extended
    else if(type_lattice=='Y')      mylattice.fcc_helical_initialize_extended_yopen(); // Y for y-dir
    else if(type_lattice=='C')      mylattice.cubic_helical_initialize();        // C for cubic
    else if(type_lattice=='Q')      mylattice.quadratic_helical_initialize();    // Q for quadratic
    else if(type_lattice=='O')      mylattice.chain_periodic_initialize();       // O for one-dimensional

    // Printing information about the q-vectors straight away
    ofstream xyzFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_xyz.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFile.open(filename);
    delete filename;

    xyzFile << "xyz-positions are listed units of grid lengths" << endl;

    cout << "File for printing xyz-coord to file is initiated" << endl;

    double x;
    double y;
    double z;

    int N  = mylattice.N;

    if(mylattice.dim==1) // Not really neccessary. Not double tested, either
    {
        for(int i=0; i<N; i++)    xyzFile << "Running index: " << i << "; Index : i: " << i << "; Position : x = " << i << endl;
    }
    if(mylattice.dim==2)
    {

        cout << "For quadratic, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];
            vector<double> pos = mylattice.sitepositions[i];
            x = ns[0];
            y = ns[1];
            xyzFile << "Running index: " << i << "; Index : i = " << x << ";   j = " << y << ";   Position : x = ";
            x = pos[0];
            y = pos[1];
            xyzFile << x << ";   y = " << y << endl;
        }
        cout << "Done printing to xyzFile" << endl;

    }

    else if(mylattice.dim==3)
    {
        cout << "For cubic or fcc, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns     = mylattice.sitecoordinates[i];
            vector<double> pos = mylattice.sitepositions[i];
            //cout << "ns retrieved" << endl;
            // These must be changed if we change into possibly setting L1, L2, L3 different
            // Don't really need b1, b2, b3, could just use a1, a2, a3 multiplied by 2*M_PI...
            // Could be more general, but we don't need to play around with our lattices that much...
            x = ns[0];
            y = ns[1];
            z = ns[2];
            // Print to file. Site number, qx, qy, qz.
            xyzFile << "Running index: " << i << ";  Index :   i = " << x << ";   j = " << y <<  ";   k = " << z <<"; Position :   x = ";
            x = pos[0];
            y = pos[1];
            z = pos[2];
            xyzFile << x << ";   y = " << y << ";   z = " << z << endl;
        }
        //cout << "Done printing to qFile" << endl;
    }
    xyzFile.close();
}

void lattice_coordinates_xyz_lines(int L, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, isotropic, sianisotropy, magfield, dm);
    mylattice.fcc_helical_initialize_extended();

    ofstream xyzFilex;
    char *filenamex = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamex, "%s_xyzxline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFilex.open(filenamex);
    delete filenamex;

    ofstream xyzFiley;
    char *filenamey = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamey, "%s_xyzyline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFiley.open(filenamey);
    delete filenamey;

    ofstream xyzFilez;
    char *filenamez = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamez, "%s_xyzzline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFilez.open(filenamez);
    delete filenamez;

    xyzFilex << "xyz-positions are listed units of grid lengths" << endl;
    xyzFiley << "xyz-positions are listed units of grid lengths" << endl;
    xyzFilez << "xyz-positions are listed units of grid lengths" << endl;

    cout << "File for printing xyz-coord a la x-, y- and z-line to file is initiated" << endl;

    int N  = mylattice.N;

    int i,j,k;
    for(int n=0; n<N; n++)
    {
        vector<int> ns = mylattice.sitecoordinates[n];
        i = ns[0];
        j = ns[1];
        k = ns[2];

        vector<double> posx = mylattice.giveposition_fcc_lines(i,j,k,'x');
        vector<double> posy = mylattice.giveposition_fcc_lines(i,j,k,'y');
        vector<double> posz = mylattice.giveposition_fcc_lines(i,j,k,'z');

        // Print to file. Site number, qx, qy, qz.
        double xlx = posx[0];    double ylx = posy[0];    double zlx = posz[0];
        double xly = posx[1];    double yly = posy[1];    double zly = posz[1];
        double xlz = posx[2];    double ylz = posy[2];    double zlz = posz[2];

        xyzFilex << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        xyzFilex << "; Position: x = " << xlx << "; y = " << xly << "; z = " << xlz << endl;
        xyzFiley << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        xyzFiley << "; Position: x = " << ylx << "; y = " << yly << "; z = " << ylz << endl;
        xyzFilez << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        xyzFilez << "; Position: x = " << zlx << "; y = " << zly << "; z = " << zlz << endl;
    }
    xyzFilex.close();
    xyzFiley.close();
    xyzFilez.close();
    cout << "Done!" << endl;
}

void lattice_coordinates_xyz_lines(int L1, int L2, int L3, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;
    long int seed = 59;
    Lattice mylattice = Lattice(L1, L2, L3, seed, isotropic, sianisotropy, magfield, dm);
    mylattice.fcc_helical_initialize_extended();

    ofstream xyzFilex;
    char *filenamex = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamex, "%s_xyzxline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFilex.open(filenamex);
    delete filenamex;

    ofstream xyzFiley;
    char *filenamey = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamey, "%s_xyzyline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFiley.open(filenamey);
    delete filenamey;

    ofstream xyzFilez;
    char *filenamez = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamez, "%s_xyzzline.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyzFilez.open(filenamez);
    delete filenamez;

    xyzFilex << "xyz-positions are listed units of grid lengths" << endl;
    xyzFiley << "xyz-positions are listed units of grid lengths" << endl;
    xyzFilez << "xyz-positions are listed units of grid lengths" << endl;

    cout << "File for printing xyz-coord a la x-, y- and z-line to file is initiated" << endl;

    int N  = mylattice.N;

    int i,j,k;
    for(int n=0; n<N; n++)
    {
        vector<int> ns = mylattice.sitecoordinates[n];
        i = ns[0];
        j = ns[1];
        k = ns[2];

        vector<double> posx = mylattice.giveposition_fcc_lines(i,j,k,'x');
        vector<double> posy = mylattice.giveposition_fcc_lines(i,j,k,'y');
        vector<double> posz = mylattice.giveposition_fcc_lines(i,j,k,'z');

        // Print to file. Site number, qx, qy, qz.
        double xlx = posx[0];    double ylx = posy[0];    double zlx = posz[0];
        double xly = posx[1];    double yly = posy[1];    double zly = posz[1];
        double xlz = posx[2];    double ylz = posy[2];    double zlz = posz[2];

        xyzFilex << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        xyzFilex << "; Position: x = " << xlx << "; y = " << xly << "; z = " << xlz << endl;
        xyzFiley << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        xyzFiley << "; Position: x = " << ylx << "; y = " << yly << "; z = " << ylz << endl;
        xyzFilez << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        xyzFilez << "; Position: x = " << zlx << "; y = " << zly << "; z = " << zlz << endl;
    }
    xyzFilex.close();
    xyzFiley.close();
    xyzFilez.close();
    cout << "Done!" << endl;
}

void reciprocallattice_coordinates(int L, char type_lattice, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, isotropic, sianisotropy, magfield, dm);

    if(type_lattice=='F')           mylattice.fcc_helical_initialize();          // F for fcc
    else if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended(); // E for extended
    else if(type_lattice=='C')      mylattice.cubic_helical_initialize();        // C for cubic
    else if(type_lattice=='Q')      mylattice.quadratic_helical_initialize();    // Q for quadratic
    else if(type_lattice=='O')      mylattice.chain_periodic_initialize();       // O for one-dimensional

    // Printing information about the q-vectors straight away
    ofstream qFile;
    char *filenameq = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenameq, "%s_qs.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    qFile.open(filenameq);
    delete filenameq;

    cout << "File for printing q-vector to file is initiated" << endl;

    double qx;
    double qy;
    double qz;

    int L1 = mylattice.L1;
    int L2 = mylattice.L2;
    int L3 = mylattice.L3;
    int N  = mylattice.N;

    if(mylattice.dim==1) // Not really neccessary. Not double tested, either
    {
        for(int i=0; i<N; i++)    qFile << "Running index = " << i << ";  qx = " << 2*M_PI/L1*i << endl;
    }
    if(mylattice.dim==2)
    {

        cout << "For quadratic, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];
            //cout << "ns retrieved" << endl;
            qx = ns[0]*mylattice.b1[0]/L1; // Not so sure about this, just need it to compile
            qy = ns[1]*mylattice.b2[1]/L2; // Double check
            //cout << "qvec set" << endl;
            qFile << "Running index = " << i << "; Indices: i = " << ns[0] << "; j =  " << ns[1] << "; Positions in q-space:  qx = " << qx << ";  qy =  " << qy << " " << endl;

        }
        cout << "Done printing to qFile" << endl;

    }

    else if(mylattice.dim==3)
    {
        cout << "For cubic or fcc, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];

            qx = ns[0]*mylattice.b1[0]/L1 + ns[1]*mylattice.b2[0]/L2 + ns[2]*mylattice.b3[0]/L3;
            qy = ns[0]*mylattice.b1[1]/L1 + ns[1]*mylattice.b2[1]/L2 + ns[2]*mylattice.b3[1]/L3;
            qz = ns[0]*mylattice.b1[2]/L1 + ns[1]*mylattice.b2[2]/L2 + ns[2]*mylattice.b3[2]/L3;
            // Print to file. Site number, qx, qy, qz.
            qFile << "Running index = " << i << "; Indices: i = " << ns[0] << ";  j = " << ns[1] << "; k = " << ns[2] << "; Positions in q-space:   qx = " << qx << ";   qy  = " << qy << ";    qz = " << qz << endl;
        }
        //cout << "Done printing to qFile" << endl;
    }
    qFile.close();
}

void reciprocallattice_coordinates(int L1, int L2, int L3, char type_lattice, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;
    long int seed = 59;
    Lattice mylattice = Lattice(L1, L2, L3, seed, isotropic, sianisotropy, magfield, dm);

    if(type_lattice=='F')           mylattice.fcc_helical_initialize();                // F for fcc
    else if(type_lattice=='E')      mylattice.fcc_helical_initialize_extended();       // E for extended
    else if(type_lattice=='Y')      mylattice.fcc_helical_initialize_extended_yopen(); // Y for open in y-dir
    else if(type_lattice=='C')      mylattice.cubic_helical_initialize();              // C for cubic
    else if(type_lattice=='D')      mylattice.cubic_helical_initialize_extended();     // D for cubic extended
    else if(type_lattice=='Q')      mylattice.quadratic_helical_initialize();          // Q for quadratic
    else if(type_lattice=='Q')      mylattice.quadratic_helical_initialize_extended(); // R for quadr extended
    else if(type_lattice=='O')      mylattice.chain_periodic_initialize();             // O for one-dim

    // Printing information about the q-vectors straight away
    ofstream qFile;
    char *filenameq = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenameq, "%s_qs.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    qFile.open(filenameq);
    delete filenameq;

    cout << "File for printing q-vector to file is initiated" << endl;

    double qx;
    double qy;
    double qz;

    // L1, L2 and L3 are already given
    int N  = mylattice.N;

    if(mylattice.dim==1) // Not really neccessary. Not double tested, either
    {
        for(int i=0; i<N; i++)    qFile << "Running index = " << i << ";  qx = " << 2*M_PI/L1*i << endl;
    }
    if(mylattice.dim==2)
    {

        cout << "For quadratic, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];
            //cout << "ns retrieved" << endl;
            qx = ns[0]*mylattice.b1[0]/L1; // Not so sure about this, just need it to compile
            qy = ns[1]*mylattice.b2[1]/L2; // Double check
            //cout << "qvec set" << endl;
            qFile << "Running index = " << i << "; Indices: i = " << ns[0] << "; j =  " << ns[1] << "; Positions in q-space:  qx = " << qx << ";  qy =  " << qy << " " << endl;
        }
        cout << "Done printing to qFile" << endl;
    }
    else if(mylattice.dim==3)
    {
        cout << "For cubic or fcc, in if-test" << endl;
        for(int i=0; i<N; i++) // Possibly only up to N/2.
        {
            vector<int> ns = mylattice.sitecoordinates[i];

            qx = ns[0]*mylattice.b1[0]/L1 + ns[1]*mylattice.b2[0]/L2 + ns[2]*mylattice.b3[0]/L3;
            qy = ns[0]*mylattice.b1[1]/L1 + ns[1]*mylattice.b2[1]/L2 + ns[2]*mylattice.b3[1]/L3;
            qz = ns[0]*mylattice.b1[2]/L1 + ns[1]*mylattice.b2[2]/L2 + ns[2]*mylattice.b3[2]/L3;
            // Print to file. Site number, qx, qy, qz.
            qFile << "Running index = " << i << "; Indices: i = " << ns[0] << ";  j = " << ns[1] << "; k = " << ns[2] << "; Positions in q-space:   qx = " << qx << ";   qy  = " << qy << ";    qz = " << qz << endl;
        }
        //cout << "Done printing to qFile" << endl;
    }
    qFile.close();
}

void reciprocallattice_coordinates_xyzline(int L, char type_lattice, string latticefilenamePrefix)
{
    bool isotropic = false; bool sianisotropy = false; bool magfield = false; bool dm = false;
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, isotropic, sianisotropy, magfield, dm);
    mylattice.fcc_helical_initialize_extended(); // We default this as well.

    ofstream qFilex;
    char *filenamex = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamex, "%s_xlineq_verbose.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    qFilex.open(filenamex);
    delete filenamex;

    ofstream qFiley;
    char *filenamey = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamey, "%s_ylineq_verbose.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    qFiley.open(filenamey);
    delete filenamey;

    ofstream qFilez;
    char *filenamez = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamez, "%s_zlineq_verbose.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    qFilez.open(filenamez);
    delete filenamez;

    qFilex << "q-vectors are listed units of grid lengths" << endl;
    qFiley << "q-vectors are listed units of grid lengths" << endl;
    qFilez << "q-vectors are listed units of grid lengths" << endl;

    cout << "Files for printing q-vector a la x-, y- and z-line to file are initiated" << endl;

    int N  = mylattice.N;

    int i, j, k;
    for(int n=0; n<N; n++)
    {
        vector<int> ns = mylattice.sitecoordinates[n];
        i = ns[0];
        j = ns[1];
        k = ns[2];

        vector<double> qvecx = mylattice.giveqvector_fcc_lines(i,j,k,'x');
        vector<double> qvecy = mylattice.giveqvector_fcc_lines(i,j,k,'y');
        vector<double> qvecz = mylattice.giveqvector_fcc_lines(i,j,k,'z');

        // Print to file. Site number, qx, qy, qz.
        double qxlx = qvecx[0];    double qylx = qvecy[0];    double qzlx = qvecz[0];
        double qxly = qvecx[1];    double qyly = qvecy[1];    double qzly = qvecz[1];
        double qxlz = qvecx[2];    double qylz = qvecy[2];    double qzlz = qvecz[2];

        qFilex << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        qFilex << "; Vector: qx = " << qxlx << "; qy = " << qxly << "; qz = " << qxlz << endl;
        qFiley << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        qFiley << "; Vector: qx = " << qylx << "; qy = " << qyly << "; qz = " << qylz << endl;
        qFilez << "Running index: " << n << "; Indices: i = " << i << "; j = " << j << "; k = " << k;
        qFilez << "; Vector: qx = " << qzlx << "; qy = " << qzly << "; qz = " << qzlz << endl;
    }
    qFilex.close();
    qFiley.close();
    qFilez.close();

    cout << "DONE!" << endl;
}


void test_fftw(int L, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    // All these parameters are irrelevant for this simple test
    int eqsteps = 1000; int mcsteps_inbin = 1000; int no_of_bins= 100;
    long int latticeseed = 59; long int MCseed = 79;
    bool isotropic = true; bool sianisotropy = false; bool magfield = false; bool dm = false;
    bool nextnearest = false;
    bool periodic = true;

    // These we want to turn off
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = false;
    bool dobootstrap = false;


    // Set these
    char type_lattice = 'O';
    string filenamePrefix = "discard"; // Want to know that this file is unimportant

    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.testFFTW(); // Test FFTW
    mymc.endsims();  // Close the file
}

void test_fftw_againstsims(int L, int eqsteps, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, char type_lattice, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    // Input irrelevant for this test
    int mcsteps_inbin = 20000;
    int no_of_bins    = 1;
    // Not strictly irrelevant, but not crucial. Don't need to test for all seeds
    long int latticeseed = 59; long int MCseed = 79;

    // These we want to turn off
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = false; // We don't need this for the testing function
    bool dobootstrap = false;

    // So it's easy to throw it away afterwards
    string filenamePrefix = "discard"; // Want to know that this file is unimportant

    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.compareFFTW_withmanual(beta); // Test FFTW
    mymc.endsims();  // Close the file
}

void test_fftw_againstsims_av(int L, int eqsteps, double beta, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool nextnearest, bool periodic, char type_lattice, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    // Input irrelevant for this test
    int mcsteps_inbin = 20;
    int no_of_bins    = 1;

    // Random number generator seeds
    long int latticeseed = 59; long int MCseed = 79;

    // These we want to turn off
    bool printeveryMCstep = false;
    bool calculatespincorrelationfunction = false; // We don't need this for the testing function
    bool dobootstrap = false;

    // So it's easy to throw it away afterwards
    string filenamePrefix = "discard"; // Want to know that this file is unimportant

    // Initializing Monte Carlo
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, dobootstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);
    mymc.compareFFTW_withmanual_av(beta); // Test FFTW
    mymc.endsims();  // Close the file
}

void test_fcc_extended(int L, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    long int seed = 59;
    Lattice mylattice(L, seed, isotropic, sianisotropy, magfield, dm);
    mylattice.setstrengths(sitestrengthsin, heisenbergin, dm_in);
    mylattice.fcc_helical_initialize_extended();

    // Test different sites and their neighbours
    // Neighbours OK
    /*
    int n = 3;
    int no_of_neighbours = mylattice.no_of_neighbours;
    vector<Bond> neighbours =  mylattice.sites[n].bonds;
    vector<double> siteposn = mylattice.sitepositions[n];
    vector<int> sitecoordn = mylattice.sitecoordinates[n];
    cout << "Position , site " << n << " : [" <<  siteposn[0] << "," << siteposn[1] << "," << siteposn[2] << "]" << endl;
    cout << "Coordinates , site " << n << " : [" <<  sitecoordn[0] << "," << sitecoordn[1] << "," << sitecoordn[2] << "]" << endl;

    for(int i=0;i<no_of_neighbours;i++)
    {
        int l = neighbours[i].siteindex2;
        vector<double> siteposl = mylattice.sitepositions[l];
        vector<int> sitecoordl = mylattice.sitecoordinates[l];
        cout << "Position , neighbour " << i << " (site " << l << ") : [" <<  siteposl[0] << "," << siteposl[1] << "," << siteposl[2] << "]" << endl;
        cout << "Coordinates , neighbour " << i << "(site " << l << ") [" <<  sitecoordl[0] << "," << sitecoordl[1] << "," << sitecoordl[2] << "]" << endl;
    }
    */

    // Testing the implementation of J and the sian terms.
    /*
    int no_of_neighbours = mylattice.no_of_neighbours;
    // n = 5
    cout << "For n = 5" << endl;
    vector<Bond> site5bonds = mylattice.sites[5].bonds;
    cout << "Single-ion anisotropy terms: Dix = " << mylattice.sites[5].Dix << "; Diy = " << mylattice.sites[5].Diy << "; Diz = " << mylattice.sites[5].Diz << endl;
    for(int i=0; i<no_of_neighbours; i++)        cout << "i = " << i << "; J = " << site5bonds[i].J  << endl;

    int n;
    if(L<3)    n = 7;
    else       n = 21;
    cout << "Now for " << n << endl;
    vector<Bond> sitebonds = mylattice.sites[n].bonds;
    cout << "Single-ion anisotropy terms: Dix = " << mylattice.sites[n].Dix << "; Diy = " << mylattice.sites[n].Diy << "; Diz = " << mylattice.sites[n].Diz << endl;
    for(int i=0; i<no_of_neighbours; i++)        cout << "i = " << i << "; J = " << sitebonds[i].J  << endl;
    */

    // Testing the position of each point
    int N = mylattice.N;
    for(int i=0; i<N;i++)
    {
        vector<double> pos = mylattice.sitepositions[i];
        vector<int>    crd = mylattice.sitecoordinates[i];
        cout << "Site " << i << endl;
        cout << "Coordinates: [" << crd[0] << " , " << crd[1] << " , " << crd[2] << "]" << endl;
        cout << "Position:    [" << pos[0] << " , " << pos[1] << " , " << pos[2] << "]" << endl << endl;
    }

    // Finding the neighbours of each point:
    int no_of_neighbours = mylattice.no_of_neighbours;
    for(int i=0; i<N;i++)
    {
        cout << "Site " << i << endl;
        for(int j=0; j<no_of_neighbours; j++)
        {
            cout << "Neighbour " << j << ": " << mylattice.sites[i].bonds[j].siteindex2 << endl;
        }
    }

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<no_of_neighbours; j++)
        {
            int n1 = mylattice.sites[i].bonds[j].siteindex1;
            int n2 = mylattice.sites[i].bonds[j].siteindex2;
            cout << "Spin " << i <<  ", bond no. " << j << endl;
            cout << "From site " << n1 << " to site " << n2 << endl;
            cout << "The Heisenberg interaction is: " << mylattice.sites[i].bonds[j].J << endl << endl;

        }
    }
}

void test_fcc_extended_diffdims(int L1, int L2, int L3, bool isotropic, bool sianisotropy, bool magfield, bool dm, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    long int seed = 59;
    Lattice mylattice(L1, L2, L3, seed, isotropic, sianisotropy, magfield, dm);
    mylattice.setstrengths(sitestrengthsin, heisenbergin, dm_in);
    mylattice.fcc_helical_initialize_extended();

    // Test different sites and their neighbours
    // Neighbours OK
    /*
    int n = 3;
    int no_of_neighbours = mylattice.no_of_neighbours;
    vector<Bond> neighbours =  mylattice.sites[n].bonds;
    vector<double> siteposn = mylattice.sitepositions[n];
    vector<int> sitecoordn = mylattice.sitecoordinates[n];
    cout << "Position , site " << n << " : [" <<  siteposn[0] << "," << siteposn[1] << "," << siteposn[2] << "]" << endl;
    cout << "Coordinates , site " << n << " : [" <<  sitecoordn[0] << "," << sitecoordn[1] << "," << sitecoordn[2] << "]" << endl;

    for(int i=0;i<no_of_neighbours;i++)
    {
        int l = neighbours[i].siteindex2;
        vector<double> siteposl = mylattice.sitepositions[l];
        vector<int> sitecoordl = mylattice.sitecoordinates[l];
        cout << "Position , neighbour " << i << " (site " << l << ") : [" <<  siteposl[0] << "," << siteposl[1] << "," << siteposl[2] << "]" << endl;
        cout << "Coordinates , neighbour " << i << "(site " << l << ") [" <<  sitecoordl[0] << "," << sitecoordl[1] << "," << sitecoordl[2] << "]" << endl;
    }
    */

    // Testing the implementation of J and the sian terms.
    /*
    int no_of_neighbours = mylattice.no_of_neighbours;
    // n = 5
    cout << "For n = 5" << endl;
    vector<Bond> site5bonds = mylattice.sites[5].bonds;
    cout << "Single-ion anisotropy terms: Dix = " << mylattice.sites[5].Dix << "; Diy = " << mylattice.sites[5].Diy << "; Diz = " << mylattice.sites[5].Diz << endl;
    for(int i=0; i<no_of_neighbours; i++)        cout << "i = " << i << "; J = " << site5bonds[i].J  << endl;

    int n;
    if(L<3)    n = 7;
    else       n = 21;
    cout << "Now for " << n << endl;
    vector<Bond> sitebonds = mylattice.sites[n].bonds;
    cout << "Single-ion anisotropy terms: Dix = " << mylattice.sites[n].Dix << "; Diy = " << mylattice.sites[n].Diy << "; Diz = " << mylattice.sites[n].Diz << endl;
    for(int i=0; i<no_of_neighbours; i++)        cout << "i = " << i << "; J = " << sitebonds[i].J  << endl;
    */

    // Testing the position of each point
    int N = mylattice.N;
    for(int i=0; i<N;i++)
    {
        vector<double> pos = mylattice.sitepositions[i];
        vector<int>    crd = mylattice.sitecoordinates[i];
        cout << "Site " << i << endl;
        cout << "Coordinates: [" << crd[0] << " , " << crd[1] << " , " << crd[2] << "]" << endl;
        cout << "Position:    [" << pos[0] << " , " << pos[1] << " , " << pos[2] << "]" << endl << endl;
    }

    // Finding the neighbours of each point:
    int no_of_neighbours = mylattice.no_of_neighbours;
    for(int i=0; i<N;i++)
    {
        cout << "Site " << i << endl;
        for(int j=0; j<no_of_neighbours; j++)
        {
            cout << "Neighbour " << j << ": " << mylattice.sites[i].bonds[j].siteindex2 << endl;
        }
    }
}

void test_dm()
{

    double x1, y1, z1, x2, y2, z2;
    double Dx, Dy, Dz, J;

    // Input variables:

    // Spin (unnormalized)
    x1 = 1;
    y1 = 2;
    z1 = 3;

    x2 = 3;
    y2 = 5;
    z2 = 1;

    // Direction of DM
    Dx = 1.0;
    Dy = 1.0;
    Dz = 1.0;

    // Heisenberg coupling
    J = 1.0;

    double den1 = sqrt(x1*x1+y1*y1+z1*z1);
    double den2 = sqrt(x2*x2+y2*y2+z2*z2);

    double spin1x = x1/den1;
    double spin1y = y1/den1;
    double spin1z = z1/den1;

    double spin2x = x2/den2;
    double spin2y = y2/den2;
    double spin2z = z2/den2;

    double onetotwo = dm_oneneighbour_fortest(spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, Dx, Dy, Dz);
    double twotoone = dm_oneneighbour_fortest(spin2x, spin2y, spin2z, spin1x, spin1y, spin1z, Dx, Dy, Dz);
    double orthogonalitytest = dm_oneneighbour_fortest(spin1x, spin1y, spin1z, spin1x, spin1y, spin1z, Dx, Dy, Dz);

    cout << "D*(S1 x S2): " << std::setprecision(std::numeric_limits<double>::digits10 + 1) <<  onetotwo << endl;
    cout << "D*(S2 x S1): " << twotoone << endl;
    cout << "Orthogonality test: D*(S1 x S1): " << orthogonalitytest << endl;
    cout << "D*(S1 x S2)+D*(S2 x S1): " << onetotwo+twotoone << endl;

    double Jterm = J_oneneighbour_fortest(spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, J);

    cout << "J*S1*S2 + D*(S1 x S2): " << Jterm+onetotwo << endl;

    // Energy difference, spin flipped
    //Renaming spins to fit the formula
    double sx = spin1x;
    double sy = spin1y;
    double sz = spin1z;

    double sxk = spin2x;
    double syk = spin2y;
    double szk = spin2z;

    // Changing a spin arbitrarily
    double xt, yt, zt, dent;
    xt = 1;    yt = 3;    zt = 1;
    // Normalizing
    dent = sqrt(xt*xt+yt*yt+zt*zt);
    xt = xt/dent; yt = yt/dent; zt = zt/dent;
    cout << "xt: " << xt << "; yt = " << yt << "; zt = " << zt << endl;
    double sx_t = xt; double sy_t= yt; double sz_t = zt;

    double energy_new = dm_oneneighbour_fortest(xt, yt, zt, spin2x, spin2y, spin2z, Dx, Dy, Dz);

    cout << "Energy_new = " << energy_new << endl;
    double energy_diff = -(Dx*((sy-sy_t)*szk-syk*(sz-sz_t))+Dy*((sz-sz_t)*sxk-szk*(sx-sx_t))+Dz*((sx-sx_t)*syk-(sy-sy_t)*sxk));
    cout << "energy difference, dm, according to class MonteCarlo: " << energy_diff << endl;
    cout << "energy difference dm, according to our energies before and after: " << energy_new-onetotwo << endl;
    cout << "Difference between the two approaches: " << energy_diff-(energy_new-onetotwo) << endl;

    cout << "energy_new, extrapolated from MonteCarlo: " << energy_diff+onetotwo << endl;

}

double dm_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double Dx, double Dy, double Dz)
{
    return Dx*(y1*z2-y2*z1) + Dy*(z1*x2-z2*x1) +Dz*(x1*y2-y1*x2);
}

double J_oneneighbour_fortest(double x1, double y1, double z1, double x2, double y2, double z2, double J)
{
    return J*(x1*x2+y1*y2+z1*z2);
}



void test_fcc_extended_yopen()
{
    long int seed = 59;
    Lattice mylattice(4, seed, false, false, false, false);
    mylattice.fcc_helical_initialize_extended_yopen();

    int N = mylattice.N;
    int no_of_neighbours;
    int neighbour;
    vector<int> y_edge_sites;
    int no_y_edge_sites = 0;
    for(int n=0; n<N; n++)
    {
        /*
        no_of_neighbours = mylattice.sites[n].no_of_neighbours_site;
        if(no_of_neighbours==8)
        {
            y_edge_sites.push_back(n);
            no_y_edge_sites++;
            cout << "Site " << n << ": No of neighbours = " << no_of_neighbours << endl;
        }
        */
        for(int i=0; i<no_of_neighbours; i++)
        {
            neighbour = mylattice.sites[n].bonds[i].siteindex2;
            //int J = mylattice.sites[n].bonds[i].J;
            string dir = mylattice.sites[n].bonds[i].direction;
            //cout << "Neighbour " << i << ": " << neighbour << "; " << dir << " = " << J << endl;
            bool increasing = mylattice.sites[n].bonds[i].increasing;
            if(no_of_neighbours==8)
            {
                if(increasing) cout << "Neighbour no. " << i << ", index " << neighbour << "; dir = " << dir << " increasing?: yes" << endl;
                else           cout << "Neighbour no. " << i << ", index " << neighbour << "; dir = " << dir << " increasing?: no" << endl;
            }

        }
        if(no_of_neighbours==8)    cout << endl;
    }

    int yesite;
    double yesx, yesy, yesz;
    vector<double> posvec;
    cout << "Sites at the edge of the y-range:" << endl;
    for(int i=0; i<no_y_edge_sites; i++)
    {
        yesite = y_edge_sites[i];
        posvec = mylattice.sitepositions[yesite];
        yesx = posvec[0]; yesy = posvec[1]; yesz = posvec[2];
        cout << "Site index: " << yesite << " Coordinates: [ " << yesx << " , " << yesy << " , " << yesz << " ]" << endl;
    }

    /*
    cout << "For site 16:" << endl;
    posvec = mylattice.sitepositions[16];
    yesx = posvec[0]; yesy = posvec[1]; yesz = posvec[2];
    cout << "Site index: " << 16 << " Coordinates: [ " << yesx << " , " << yesy << " , " << yesz << " ]" << endl;


    //vector<int> endpointswegotfromtest;
    int no_eptest;
    int ypos;
    cout << "Endpoints we got from the if-test:" << endl;
    for(int i=0; i<N; i++)
    {
        posvec = mylattice.sitepositions[i];
        ypos = posvec[1];
        if(ypos==0 || ypos==(mylattice.L2-1))
        {
            //endpointswegotfromtest.push_back(i);
            no_eptest++;
            cout << i << endl;
        }
    }
    //vector<double> posvec = mylattice.sitepositions[k];
    //int ypos =  posvec[1];
    //if(ypos==0 || ypos==(mylattice.L2-1)) cout << "Endpoint, no of neighbours:" << nneighbours << endl;
    */

    int l;
    int eureka = 13;
    cout << "Neighbours of " << eureka << endl;
    int nnbeureka = mylattice.sites[eureka].no_of_neighbours_site;
    for(int i=0; i<nnbeureka; i++)
    {
        l = mylattice.sites[eureka].bonds[i].siteindex2;
        cout << "Neighbour " << i << ": Site " << l << endl;
    }
    cout << "Site " << eureka << " has " << nnbeureka << " neighbours" << endl;
}

void test_fcc_extended_yopen_throughMC(int L, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    int eqsteps = 1; int mcsteps_inbin = 1; int no_of_bins = 1;
    bool isotropic = true; // This is the one we are testing now
    long int latticeseed = 59; long int MCseed = 79;
    bool sianisotropy = false; bool magfield = false; bool dm = false; bool nextnearest = false;
    bool periodic = true; bool printeveryMCstep = false; bool calculatespincorrelationfunction = false;
    bool doboolstrap = false;
    char type_lattice = 'Y';
    string filenamePrefix = "tezter";
    MonteCarlo mymc(L, L, L, eqsteps, mcsteps_inbin, no_of_bins, latticeseed, MCseed, isotropic, sianisotropy, magfield, dm, nextnearest, periodic, printeveryMCstep, calculatespincorrelationfunction, doboolstrap, type_lattice, filenamePrefix, sitestrengthsin, heisenbergin, dm_in);

    mymc.testyopenfcc();
}

void checkneighbours(int L, char type_lattice, bool periodic, vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in)
{
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, false, false, false, false); // We only look at the neighbours
    mylattice.setstrengths(sitestrengthsin, heisenbergin, dm_in);

    // We only look at periodic functions here
    if(periodic)
    {
        if(type_lattice=='E')         mylattice.fcc_helical_initialize_extended();
        else if(type_lattice=='Y')    mylattice.fcc_helical_initialize_extended_yopen();
        else if(type_lattice=='D')    mylattice.cubic_helical_initialize_extended();
        else if(type_lattice=='R')    mylattice.quadratic_helical_initialize_extended();
        else if(type_lattice=='F')    mylattice.fcc_helical_initialize();
        else if(type_lattice=='C')    mylattice.cubic_helical_initialize();
        else if(type_lattice=='Q')    mylattice.quadratic_helical_initialize();
        else if(type_lattice=='O')    mylattice.chain_periodic_initialize();
    }
    else
    {
        if(type_lattice=='O')    mylattice.chain_open_initialize();
        else                     cout << "Error! For open BCs, only the chain is implemented!!!" << endl;
    }

    if(type_lattice=='Y')         periodic = false;

    //cout << "mylattice.dim = " << mylattice.dim << endl;

    int neighbour;
    int no_of_neighbours;
    int N                = mylattice.N;
    bool inspectall = true;
    if(periodic) no_of_neighbours = mylattice.no_of_neighbours;
    Bond nbond;
    if(inspectall)
    {
        //N = 52;
        cout << "In inspectall" << endl;
        for(int n=0; n<N; n++)
        {
            cout << "The neighbours of spin " << n << endl;
            if(!periodic)    no_of_neighbours = mylattice.sites[n].no_of_neighbours_site;
            cout << "No of neighbours = " << no_of_neighbours << endl;
            for(int i=0; i<no_of_neighbours; i++)
            {
                nbond = mylattice.sites[n].bonds[i];
                neighbour = nbond.siteindex2;
                cout << "Neighbour " << i << ": " << neighbour << "; " << nbond.direction << " = " << nbond.J << endl;
            }
            cout << endl;
        }
    }
    else
    {
        int n = 1; // Some spin in the system
        cout << "n = " << n << endl << "Neighbours: ";
        for(int i=0; i<no_of_neighbours; i++)
        {
            neighbour = mylattice.sites[n].bonds[i].siteindex2;
            cout << "Neighbour " << i << ": " << neighbour << endl;
        }
    }

}

void printyneighbours_fcc(int L, string latticefilenamePrefix)
{

    ofstream yneighbourFile;
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_yneighbours_eachpoint.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    yneighbourFile.open(filename);
    delete filename;

    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, false, false, false, false); // We only look at the neighbours
    mylattice.fcc_helical_initialize_extended(); // So that we have next nearest neighbours

    int N = mylattice.N;
    int yneigh;
    for(int i=0; i<N; i++)
    {
        yneigh = mylattice.sites[i].nextnearesty[1].siteindex2; // This is for periodic BCs
        yneighbourFile << i << " " << yneigh << endl;
    }
    yneighbourFile.close();
}


void printnearestneighbours_fcc(int L, string latticefilenamePrefix)
{

    ofstream xyneighbourFile;
    char *filenamexy = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamexy, "%s_xyneighbours_eachpoint.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xyneighbourFile.open(filenamexy);
    delete filenamexy;

    ofstream xzneighbourFile;
    char *filenamexz = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamexz, "%s_xzneighbours_eachpoint.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    xzneighbourFile.open(filenamexz);
    delete filenamexz;

    ofstream yzneighbourFile;
    char *filenameyz = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenameyz, "%s_yzneighbours_eachpoint.txt", latticefilenamePrefix.c_str() );   // Create filename with prefix and ending
    yzneighbourFile.open(filenameyz);
    delete filenameyz;

    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, false, false, false, false); // We only look at the neighbours
    mylattice.fcc_helical_initialize_extended(); // So that we have next nearest neighbours

    int N = mylattice.N;
    int neigh;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            neigh = mylattice.sites[i].bonds[j].siteindex2; // This is for periodic BCs
            // The way it is ordered:
            if(j==4 || j==5 || j==6 || j==7)              xyneighbourFile << i << " " << neigh << endl;
            if(j==0 || j==1 || j==10 || j==11)            xzneighbourFile << i << " " << neigh << endl;
            if(j==2 || j==3 || j==8 || j==9)              yzneighbourFile << i << " " << neigh << endl;
        }

    }
    xyneighbourFile.close();
    xzneighbourFile.close();
    yzneighbourFile.close();
}

void testnestnearestneighbour_chain(int L)
{
    long int seed = 59;
    Lattice mylattice = Lattice(L, seed, false, false, false, false);
    mylattice.chain_periodic_initialize();

    for(int i=0; i<L; i++)
    {
        int nbm = mylattice.sites[i].nextnearesty[0].siteindex2;
        int nbp = mylattice.sites[i].nextnearesty[1].siteindex2;

        cout << "Spin " << i << ": next-nearest neighbours: " << nbm << ", " << nbp << endl;
    }
}

void printrandomnumbers()
{
    ofstream rannumFile;
    string filenamePrefix = "randomnumbergenerator";
    char *filename = new char[1000];                                // File name can have max 1000 characters
    sprintf(filename, "%s_test.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    rannumFile.open(filename);
    delete filename;

    long int seed = 59;
    for(int i=0; i<1e5; i++)    rannumFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << ran2(&seed) << endl;
    rannumFile.close();
}

double converttemp(double intemp)
{
    return 11.6045221/intemp;
}
