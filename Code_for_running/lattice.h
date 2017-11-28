#ifndef LATTICE_H
#define LATTICE_H
#include <iostream>
#include <math.h>
#include <cmath>
#include "bond.h"
#include "site.h"
#include "gaussiandeviate.h"
//#include <armadillo>  // in this case, enter LIBS += -larmadillo -llapack -lblas in .pro-file

using namespace std; // Now, I can remove all stds.

class Lattice
{
public:

    double hx, hy, hz, Dix, Diy, Diz;     // Site strengths
    double J, Jx, Jy, Jz, Jxy, Jxz, Jyz;  // Bond strengths, Heisenberg
    double Dx, Dy, Dz;                    // Bond strengths, DM

    bool isotropic,  dm;         // Bools for n.n. terms
    bool sianisotropy, magfield; // Bools for site terms
    bool notperiodic;            // Bool to keep track of neighbours
    bool dimequal;               // Bool for specifying whether we have dimensions of equal length, i.e. LxLxL.
    bool systemstrengthsgiven;   // Bool for indicating whether system strengths are given. Actions will be taken if not
    bool extended;               // Bool for whether we use an extended class or not (different attributes)

    int dim, L, N, no_of_neighbours;
    int L1, L2, L3;              // For when we have unequal dimensions
    long int seed;

    // Typedefs
    typedef vector<int> vecint;
    typedef vector<vecint> intmatrix;

    vector<int> dimlengths;
    vector<int> yline;
    vector<double> a1, a2, a3; // primitive vectors. Mainly double because of the fcc
    vector<double> b1, b2, b3; // Primitive vectors of the reciprocal lattice

    //std::vector<Bond> bonds;
    std::vector<Site> sites;
    //std::vector<std::vector<double> > sitepositions_chain; // In case we introduce a grid length
    std::vector<std::vector<double> > sitepositions;
    std::vector<std::vector<int> >    sitecoordinates;
    std::vector<std::vector<int> >    siteneighbours;

    std::vector<bool>   yhalfsite_vec;

    // Initialization
    Lattice();
    Lattice(int L, long int seed, bool isotropic, bool sianisotropy, bool magfield, bool dm);
    Lattice(int L1, int L2, int L3, long int seed, bool isotropic, bool sianisotropy, bool magfield, bool dm);

    // Debugging
    //void setmajordebug();

    int findneighbour(int n, int toi, int toj, int tok);
    int findneighbour2D(int n, int toi, int toj);

    //Lattice grid functions
    //void chain_2p_periodic_initialize();
    void chain_periodic_initialize();
    void chain_open_initialize();
    void quadratic_helical_initialize();
    void quadratic_helical_initialize_extended();
    void cubic_helical_initialize();
    void cubic_helical_initialize_extended();
    void fcc_helical_initialize();
    void fcc_helical_initialize_extended();
    void fcc_helical_initialize_extended_yopen();

    // Feed interaction functions
    void setstrengths(vector<double> sitestrengthsin, vector<double> heisenbergin, vector<double> dm_in);
    void givestrengths_automatic();

    std::vector<double> giveposition_fcc_lines(int i, int j, int k, char letter);
    std::vector<double> giveqvector_fcc_lines(int i, int j, int k, char letter);
    std::vector<double> giveqvector_fcc(int i, int j, int k);
    std::vector<int> fccyline(); // Line of points in the (0,y,0)-direction
    std::vector<int> fccxline(); // Line of points in the (x,0,0)-direction
    std::vector<int> fcczline(); // Line of points in the (0,0,z)-direction
    std::vector<int> fccqyline(); // Line of points in the (0,y,0)-direction
    std::vector<int> fccqxline(); // Line of points in the (x,0,0)-direction
    std::vector<int> fccqzline(); // Line of points in the (0,0,z)-direction
    std::vector<int> fccqdline();
    std::vector<int> cubicyline(); // Line of points in the (0,y,0)-direction
    std::vector<int> cubicxline(); // Line of points in the (x,0,0)-direction
    std::vector<int> cubiczline(); // Line of points in the (0,0,z)-direction
    std::vector<int> quadrxline(); // Line of points in the (x,0)-direction
    std::vector<int> quadryline(); // Line of points in the (0,y)-direction
    std::vector<int> fccyline_shifted(double xshift, double zshift); //Line (0,y,0)+(a,b,c)
    std::vector<int> diagline_cubic();
    std::vector<int> diagline_quadr();
    std::vector<int> diagline_fcc();

};

#endif // LATTICE_H
