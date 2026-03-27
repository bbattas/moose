#pragma once

#include "Material.h"
#include "GrainTracker.h"
#include "EulerAngleTxtFileReader.h"
// #include "EBSDReader.h"

/**
 * Visualize the location of grain boundaries in a polycrystalline simulation.
 */
class GBMisorientationTxt : public Material
{
public:
  static InputParameters validParams();

  GBMisorientationTxt(const InputParameters & parameters);

  // Returned by getMisorientationFromQuaternion()
  struct MisorientationResult
  {
    double theta;    // |misorientation| in radians
    double polar_ax; // polar
    double azim_ax;  // azimuth
    Eigen::Quaternion<Real> q;
    Real qnorm; // norm of minimizing quaternion (sanity check)
  };

protected:
  /// Necessary override. This is where the property values are set.
  virtual void computeQpProperties() override;

  /// Function to obtain the total line number in misorientation angle file
  virtual unsigned int getTotalLineNum() const;
  /// Function to obtain line number for a given grain pair
  virtual unsigned int getLineNum(unsigned int grain_i, unsigned int grain_j);
  // /// Function to get the GB type for triple junctions
  // virtual Real getTripleJunctionType();
  /// Function to convert symmetry matrix to quaternion form
  void rotationSymmetryToQuaternion(const double O[3][3], Eigen::Quaternion<Real> & q);
  /// Function to define the symmetry operator
  void defineSymmetryOperator();
  /// Function to return the misorientation of two quaternions
  MisorientationResult getMisorientationFromQuaternion(const Eigen::Quaternion<Real> & qi,
                                                       const Eigen::Quaternion<Real> & qj);
  /// Get the Misorientation angle
  void getMisorientationAngles();

  /// Grain tracker object
  const GrainTracker & _grain_tracker;

  /// EBSD reader user object
  // const EBSDReader & _ebsd_reader;
  const EulerAngleProvider & _euler;

  /// Parameters to calculate the Misorientation angle file
  std::vector<Real> _misorientation_angles;
  std::vector<Real> _miso_ax_polar;
  std::vector<Real> _miso_ax_azimuth;
  std::vector<Eigen::Quaternion<Real>> _qmin;
  std::vector<Real> _other;

  /// parameters to store the EBSD id and corresponding value on GB
  std::vector<unsigned int> _gb_pairs;
  std::vector<Real> _gb_op_pairs;
  // std::vector<unsigned int> _gb_i_vals;

  /// order parameters
  const unsigned int _op_num;
  const std::vector<const VariableValue *> _vals;

  /// The parameters to calculate the misorientation
  std::vector<Eigen::Quaternion<Real>> _sym_quat;
  int _o_sym = 24;
  std::vector<EulerAngles> _euler_angle;
  std::vector<Eigen::Quaternion<Real>> _quat_angle;

  // Misorientation function weight (mostly for later when vs inc)
  const Real _m_weight;

  // Inputs
  // Other Free Energy terms
  const MaterialProperty<Real> & _gbe_iso;
  const MaterialProperty<Real> & _kappa;
  const MaterialProperty<Real> & _mu;

  // Output Finals
  MaterialProperty<Real> & _int_width;
  MaterialProperty<Real> & _gamma;

  // TESTING
  MaterialProperty<Real> & _gtnum;
  MaterialProperty<Real> & _other_out;
  //
  MaterialProperty<Real> & _eul_a;
  MaterialProperty<Real> & _eul_b;
  MaterialProperty<Real> & _eul_c;
  //
  MaterialProperty<Real> & _quat_a;
  MaterialProperty<Real> & _quat_b;
  MaterialProperty<Real> & _quat_c;
  MaterialProperty<Real> & _quat_d;
  MaterialProperty<Real> & _quat_mag;
  //
  MaterialProperty<Real> & _miso;
  MaterialProperty<Real> & _miso_polar;
  MaterialProperty<Real> & _miso_azim;
  //
  MaterialProperty<Real> & _miso_ang_en;
  MaterialProperty<Real> & _miso_ax_en;
  //
  MaterialProperty<Real> & _twist;
  MaterialProperty<Real> & _tilt;
  //
  MaterialProperty<Real> & _f_mis;

  // L
  const MaterialProperty<Real> & _L0;
  MaterialProperty<Real> & _L;
};
