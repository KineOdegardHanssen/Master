#include "bond.h"

Bond::Bond()
{
}

// Standard initializer
// Works for both nearest and next nearest neighbour Bond
Bond::Bond(int siteindex1, int siteindex2, double J, bool increasing)
{
    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;
    this->J = J;
    this->increasing = increasing; // Is this neccessary here?
}

// Including a string for tests
Bond::Bond(int siteindex1, int siteindex2, double J, bool increasing, string direction)
{
    this->J = J;
    this->siteindex1 = siteindex1;
    this->siteindex2 = siteindex2;
    this->increasing = increasing;
    this->direction = direction;
}
