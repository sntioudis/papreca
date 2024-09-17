#ifdef FIX_CLASS
FixStyle(papreca, FixPAPRECA)
#else

#ifndef LMP_FIX_PAPRECA_H
#define LMP_FIX_PAPRECA_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPAPRECA : public Fix {
 public:
  FixPAPRECA(class LAMMPS *, int, char **);
  ~FixPAPRECA();
  int setmask() override;
  void init() override;
  void init_list(int id, NeighList *ptr) override;
 private:
  class NeighList *nlist_half;
  class NeighList *nlist_full;
};

}

#endif
#endif

