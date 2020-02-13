#include "cctk.h"
#include "params.h"

namespace TPID {

#if ID_PARS == 0
  double par_m_plus =  0.5538461538461539;
  double par_m_minus =  0.4461538461538462;
#elif ID_PARS == 1
  double par_m_plus = 1.5;
  double par_m_minus = 1.0;
#elif ID_PARS == 2
  double par_m_plus = 0.4824;
  double par_m_minus = 0.4824;
#endif

}
