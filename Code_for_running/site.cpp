#include "site.h"

Site::Site(int n, double spinx, double spiny, double spinz, std::vector<Bond> bonds)
{
    index = n;
    this->spinx = spinx;
    this->spiny = spiny;
    this->spinz = spinz;
    this->bonds = bonds;
}

Site::Site(int n, int no_of_neighbours_site, double spinx, double spiny, double spinz, std::vector<Bond> bonds)
{   // I guess I don't need this anymore...
    index = n;
    this->no_of_neighbours_site = no_of_neighbours_site;
    this->spinx = spinx;
    this->spiny = spiny;
    this->spinz = spinz;
    this->bonds = bonds;
}


Site::Site(int n, double spinx, double spiny, double spinz, std::vector<Bond> bonds, std::vector<Bond> nextnearesty, std::vector<Bond> nextnearestz)
{
    index = n;
    this->spinx = spinx;   // Have spinx in a new class, States?
    this->spiny = spiny;
    this->spinz = spinz;
    this->bonds = bonds;
    this->nextnearesty = nextnearesty;
    this->nextnearestz = nextnearestz;
}

Site::Site(int n, int no_of_neighbours_site, int no_of_nneighbours_site, double spinx, double spiny, double spinz, std::vector<Bond> bonds, std::vector<Bond> nextnearesty)
{
    index = n;
    this->spinx = spinx;   // Have spinx in a new class, States?
    this->spiny = spiny;
    this->spinz = spinz;
    this->bonds = bonds;
    this->nextnearesty = nextnearesty;
    this->no_of_neighbours_site = no_of_neighbours_site;
    this->no_of_nneighbours_site = no_of_nneighbours_site;
}

Site::Site(int n, int no_of_neighbours_site, int no_of_nneighbours_site, double spinx, double spiny, double spinz, std::vector<Bond> bonds, std::vector<Bond> nextnearesty, std::vector<Bond> nextnearestz)
{
    index = n;
    this->spinx = spinx;   // Have spinx in a new class, States?
    this->spiny = spiny;
    this->spinz = spinz;
    this->bonds = bonds;
    this->nextnearesty = nextnearesty;
    this->nextnearestz = nextnearestz;
    this->no_of_neighbours_site = no_of_neighbours_site;
    this->no_of_nneighbours_site = no_of_nneighbours_site;
}
