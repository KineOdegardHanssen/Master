#ifndef SITE_H
#define SITE_H
#include <iostream>
#include <vector>
#include <bond.h>


class Site
{
public:

    int index;
    double spinx, spiny, spinz;

    std::vector<Bond> bonds;
    std::vector<Bond> nextnearesty;
    std::vector<Bond> nextnearestz;

    // For open boundary conditions
    int no_of_neighbours_site;
    int no_of_nneighbours_site;

    // Initializers
    // For periodic boundary conditions
    Site(int n, double spinx, double spiny, double spinz, std::vector<Bond> bonds);
    // For open boundary conditions
    Site(int n, int no_of_neighbours_site, double spinx, double spiny, double spinz, std::vector<Bond> bonds);
    // For including next nearest neighbour interactions
    Site(int n, double spinx, double spiny, double spinz, std::vector<Bond> bonds, std::vector<Bond> nextnearesty, std::vector<Bond> nextnearestz);
    Site(int n, int no_of_neighbours_site, int no_of_nneighbours_site, double spinx, double spiny, double spinz, std::vector<Bond> bonds, std::vector<Bond> nextnearesty);
    Site(int n, int no_of_neighbours_site, int no_of_nneighbours_site, double spinx, double spiny, double spinz, std::vector<Bond> bonds, std::vector<Bond> nextnearesty, std::vector<Bond> nextnearestz);

};

#endif // SITE_H
