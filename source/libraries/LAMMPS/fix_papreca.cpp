#include "fix_papreca.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "error.h"
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPAPRECA::FixPAPRECA(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  if (narg != 3) error->all(FLERR, "Illegal papreca command. Use this fix as: fix papreca all papreca");
  if ( strcmp( arg[0] , "papreca" ) != 0 ) error->all(FLERR, "Illegal papreca command. Group id MUST be papreca. Only use this fix as: fix papreca all papreca");
  if ( strcmp( arg[1] , "all" ) != 0 ) error->all(FLERR, "Illegal papreca command. This fix has to be applied to the all group. Only use this fix as: fix papreca all papreca");
}

/* ---------------------------------------------------------------------- */

FixPAPRECA::~FixPAPRECA() {
  // Destructor if needed
}

/* ---------------------------------------------------------------------- */

int FixPAPRECA::setmask() {
  int mask = 0;
  return mask;
}


/* ---------------------------------------------------------------------- */

void FixPAPRECA::init() {

  // Request a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL )->set_id(1);
  
  // Request a half neighbor list
  neighbor->add_request( this )->set_id(2);
  
}

void FixPAPRECA::init_list(int id, NeighList *ptr)
{

 if( id == 1 ){
   nlist_full = ptr;
 }else if( id == 2 ){
   nlist_half = ptr;
 }else{
   error->all(FLERR, "Error in papreca.cpp in function init_list: ptr could not be assigned to neiblist");
 }
  
}

