#pragma once

#include "GBInclinationBase.h"
#include "GBMisorientationHelper.h"

class GBCombinedAnisotropyMaterial : public GBInclinationBase, protected GBMisorientationHelper
{
public:
  static InputParameters validParams();

  GBCombinedAnisotropyMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  enum GBMode
  {
    COS = 0,
    INC = 1,
    MISO = 2,
    FULL = 3
  };

  Real computeCosineOnlyGBE(Real theta_inc, Real polar_inc) const;

  Real computeInclinationGBE(Real theta_inc, Real polar_inc) const;

  Real computeMisorientationGBE(const MisorientationData & miso) const;

  Real computeFullGBE(Real theta_inc, Real polar_inc, const MisorientationData & miso) const;

protected:
  const int _gb_mode;

  MaterialProperty<std::vector<Real>> & _gamma_ij;
};
