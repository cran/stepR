#include <R.h>
#include <Rinternals.h>

/***************
* class Jump
* describes a jump in a block signal
* Thomas Hotz, 2007-2008
***************/

class Jump {
  public:
    Jump(int n, int ri, double im); // constructor
    Jump(); // default constructor initialises to "before" the data
    
    int number; // the (ordinal) number of the jump
    int rightIndex; // the index of the right end of the block described
    double improve; // the improvement that the inclusion of this jump brought
};
