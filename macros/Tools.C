#include "Tools.h"


double GetBeta(double mass, double momentum) {
  double bg = momentum/mass;
  double beta = sqrt(bg*bg/(1+bg*bg));
  return beta;
}
