#include "GBCombinedAnisotropyMaterial.h"

registerMooseObject("PhaseFieldApp", GBCombinedAnisotropyMaterial);

InputParameters
GBCombinedAnisotropyMaterial::validParams()
{
  InputParameters params = GBInclinationBase::validParams();

  params += GBMisorientationHelper::validParams();

  MooseEnum gb_mode("cos=0 inc=1 miso=2 full=3", "full");

  params.addParam<MooseEnum>("gb_mode",
                             gb_mode,
                             "Which grain-boundary property model to use. "
                             "'cos' and 'inc' use only inclination data. "
                             "'miso' and 'full' also use misorientation data.");

  params.addClassDescription(
      "Parent material combining inclination and optional misorientation data "
      "to compute pairwise grain-boundary properties.");

  return params;
}

GBCombinedAnisotropyMaterial::GBCombinedAnisotropyMaterial(const InputParameters & parameters)
  : GBInclinationBase(parameters),
    GBMisorientationHelper(parameters),
    _gb_mode(getParam<MooseEnum>("gb_mode")),
    _gamma_ij(declareProperty<std::vector<Real>>("gamma_ij"))
{
  if ((_gb_mode == MISO || _gb_mode == FULL) && !misorientationEnabled())
    mooseError(name(),
               ": gb_mode='miso' or gb_mode='full' requires "
               "enable_misorientation=true.");
}

void
GBCombinedAnisotropyMaterial::computeQpProperties()
{
  // Builds theta_ij, polar_ij, ij_i, ij_j, no_ij_pairs, and derivatives.
  GBInclinationBase::computeQpProperties();

  const auto & theta = _theta_ij[_qp];
  const auto & polar = _polar_ij[_qp];
  const auto & ij_i = _ij_i[_qp];
  const auto & ij_j = _ij_j[_qp];

  auto & gamma = _gamma_ij[_qp];

  gamma.assign(theta.size(), 0.0);

  if (_no_ij_pairs[_qp])
    return;

  for (std::size_t k = 0; k < theta.size(); ++k)
  {
    if (theta[k] == -1.0)
      continue;

    const std::size_t i = ij_i[k];
    const std::size_t j = ij_j[k];

    const Real theta_inc = theta[k];
    const Real polar_inc = polar[k];

    switch (_gb_mode)
    {
      case COS:
      {
        gamma[k] = computeCosineOnlyGBE(theta_inc, polar_inc);
        break;
      }

      case INC:
      {
        gamma[k] = computeInclinationGBE(theta_inc, polar_inc);
        break;
      }

      case MISO:
      {
        const auto & miso = getMisorientationData(i, j);
        gamma[k] = computeMisorientationGBE(miso);
        break;
      }

      case FULL:
      {
        const auto & miso = getMisorientationData(i, j);
        gamma[k] = computeFullGBE(theta_inc, polar_inc, miso);
        break;
      }

      default:
        mooseError(name(), ": unknown gb_mode = ", _gb_mode);
    }
  }
}

Real
GBCombinedAnisotropyMaterial::computeCosineOnlyGBE(const Real theta_inc, const Real polar_inc) const
{
  // Fill in later.
  // Uses:
  //   theta_inc = inclination azimuth angle
  //   polar_inc = inclination polar angle

  return 1.0;
}

Real
GBCombinedAnisotropyMaterial::computeInclinationGBE(const Real theta_inc,
                                                    const Real polar_inc) const
{
  // Fill in later.
  // More general inclination-only function.

  return 1.0;
}

Real
GBCombinedAnisotropyMaterial::computeMisorientationGBE(const MisorientationData & miso) const
{
  // Fill in later.
  // Available:
  //   miso.theta    = misorientation angle
  //   miso.polar_ax = misorientation-axis polar angle
  //   miso.azim_ax  = misorientation-axis azimuth angle
  //   miso.q        = quaternion
  //   miso.qnorm    = vector-part norm before axis normalization

  return 1.0;
}

Real
GBCombinedAnisotropyMaterial::computeFullGBE(const Real theta_inc,
                                             const Real polar_inc,
                                             const MisorientationData & miso) const
{
  // Fill in later.
  // Combines:
  //   inclination angle / polar angle
  //   misorientation angle / axis

  return 1.0;
}
