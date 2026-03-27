#include "EulerAngleTxtFileReader.h"

#include <fstream>

registerMooseObject("PhaseFieldApp", EulerAngleTxtFileReader);

InputParameters
EulerAngleTxtFileReader::validParams()
{
  InputParameters params = EulerAngleProvider::validParams();
  params.addClassDescription("Read Euler angle data from a file and provide it to other objects.");
  params.addRequiredParam<FileName>("file_name", "Euler angle data file name");
  return params;
}

EulerAngleTxtFileReader::EulerAngleTxtFileReader(const InputParameters & params)
  : EulerAngleProvider(params), _file_name(getParam<FileName>("file_name"))
{
  readFile();
}

unsigned int
EulerAngleTxtFileReader::getGrainNum() const
{
  return _angles.size();
}

const EulerAngles &
EulerAngleTxtFileReader::getEulerAngles(unsigned int i) const
{
  mooseAssert(i < getGrainNum(), "Requesting Euler angles for an invalid grain id");
  return _angles[i];
}

void
EulerAngleTxtFileReader::readFile()
{
  // Read in Euler angles from _file_name
  std::ifstream inFile(_file_name.c_str());
  if (!inFile)
    mooseError("Can't open ", _file_name);

  // Skip first 1 lines
  for (unsigned int i = 0; i < 1; ++i)
    inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // Loop over grains
  EulerAngles a;
  while (inFile >> a.phi1 >> a.Phi >> a.phi2)
    _angles.push_back(EulerAngles(a));
}
