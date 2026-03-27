#pragma once

#include "EulerAngleProvider.h"
#include <vector>

// Forward declaration

/**
 * Read a set of Euler angles from a txt file
 */
class EulerAngleTxtFileReader : public EulerAngleProvider
{
public:
  static InputParameters validParams();

  EulerAngleTxtFileReader(const InputParameters & parameters);

  virtual const EulerAngles & getEulerAngles(unsigned int) const;
  virtual unsigned int getGrainNum() const;

  virtual void initialize() {}
  virtual void execute() {}
  virtual void finalize() {}

protected:
  void readFile();

  FileName _file_name;
  std::vector<EulerAngles> _angles;
};
