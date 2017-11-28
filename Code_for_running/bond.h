#ifndef BOND_H
#define BOND_H
#include<iostream>
#include <vector>
#include <string>
//#include <site.h> // Or do I?

using namespace std;

class Bond
{
public:
    // Bond quantities
    double J;

    // Bond type
    bool increasing; // For the sign of DM
    string direction;

    // The sites which the bond connects
    int siteindex1, siteindex2; // Only need the latter, I think, but not probably not a big deal

    Bond();
    Bond(int siteindex1, int siteindex2, double J, bool increasing);  // Initializing nearest neighbour bond
    Bond(int siteindex1, int siteindex2, double J, bool increasing, string direction);
};

#endif // BOND_H
