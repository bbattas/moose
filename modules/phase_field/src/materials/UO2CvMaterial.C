#include "UO2CvMaterial.h"
#include "DelimitedFileReader.h"

registerMooseObject("PhaseFieldApp", UO2CvMaterial);

InputParameters
UO2CvMaterial::validParams()
{
  InputParameters params = DerivativeParsedMaterialHelper::validParams();
  params.addClassDescription(
      "Calculates equilibrium_vacancy_concentration for UO2+/-X "
      "to use with GrandPotentialSinteringMaterial.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addRequiredCoupledVar("T", "Temperature");
  params.addCoupledVar("OU", "Variable of O/U ratio in UO2+/-x");
  params.addRequiredCoupledVar("c", "Void phase order parameter");
  // params.addRequiredParam
  params.addRequiredParam<std::vector<FileName>>("gb_se_csv",
                                    "List of .csv file names for the segregation energies.");
  params.addParam<Real>("gb_conc_max", 0.5, "Maximum value for grain boundary"
                                            "vacancy concentration");
  params.addParam<Real>("bulk_conc_max", 0.25, "Maximum value for grain boundary"
                                            "vacancy concentration");
  // params.addParam<FileName>("gb2_se_csv",
  //                                   "File name for the 2nd set of segregation energies");
  // params.addParam<FileName>("gb3_se_csv",
  //                                   "File name for the 3rd set of segregation energies");
  return params;
}

UO2CvMaterial::UO2CvMaterial(const InputParameters & parameters)
  : DerivativeParsedMaterialHelper(parameters),
    _T(coupledValue("T")),
    _OU(coupledValue("OU")),
    _phi(coupledValue("c")),
    _kb(8.617343e-5), // Boltzmann constant in eV/K
    _op_num(coupledComponents("v")),
    _vals_name(_op_num),
    _cv_name(getParam<std::string>("f_name")),
    _cv_out(declareProperty<Real>(_cv_name)),
    _dcv(_op_num),
    _d2cv(_op_num),
    _gb_csv(getParam<std::vector<FileName>>("gb_se_csv")),
    _gb_max(getParam<Real>("gb_conc_max")),
    _b_max(getParam<Real>("bulk_conc_max"))
    // _gb2_csv(getParam<FileName>("gb2_se_csv")),
    // _gb3_csv(getParam<FileName>("gb3_se_csv"))
{
  if (_op_num == 0)
    mooseError("Model requires op_num > 0");

  _vals.resize(_op_num);//
  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _vals_name[i] = getVar("v", i)->name();
    _vals[i] = &coupledValue("v", i);
    _dcv[i] = &declarePropertyDerivative<Real>(_cv_name, _vals_name[i]);
    _d2cv[i].resize(_op_num);

    for (unsigned int j = 0; j <= i; ++j)
    {
      _d2cv[j][i] = &declarePropertyDerivative<Real>(_cv_name, _vals_name[j], _vals_name[i]);
    }
  }


  // Read the data from CSV file
  // MAKE SURE THE PATH IS CORRECT
  // Sigma 11 GB Data (trimmed)
  // MooseUtils::DelimitedFileReader sigma11_reader("/Users/bbattas/projects/marmot/csv_include/sigma11_se.csv");
  // sigma11_reader.read();
  // _sigma11_data = sigma11_reader.getData(1);
  // Sigma 9 GB Data (trimmed)
  // MooseUtils::DelimitedFileReader sigma9_reader("/Users/bbattas/projects/marmot/csv_include/sigma9_se.csv");
  // sigma9_reader.read();
  // _sigma9_data = sigma9_reader.getData(1);

  // MooseUtils::DelimitedFileReader reader(_gb_csv[0]);
  // reader.read();
  // // auto tempData = reader.getData(1);
  // _gb_data.push_back(reader.getData(1)); //= reader.getData(1);
  for (unsigned int i = 0; i < _gb_csv.size(); ++i)
  {
    MooseUtils::DelimitedFileReader reader(_gb_csv[i]);
    reader.read();
    _gb_data.push_back(reader.getData(1));
  }


}

Real lowT_line(Real plusx, Real T){
  if (plusx > 0.0) {
    Real a = -9.92395296e-01;
    Real b = 1.61688991e+00;
    Real c = 8.38156005e-04;
    Real d = 5.78661025e+00;
    Real f = 1.54217248e-02;
    // a,b,c,d,f = [-9.92395296e-01, 1.61688991e+00, 8.38156005e-04, 5.78661025e+00, 1.54217248e-02]
    Real intercept = d + (a - d)/(std::pow((1 + std::pow((plusx/c),b)),f));
    Real slope = -3.06712649e-05 * std::pow(log10(plusx),2) - 6.83078465e-04 * log10(plusx) - 1.87719009e-04;
    return slope*T + intercept;
  }
  else if (plusx == 0.0) {
    Real intercept = 2.7177198288000017;
    Real slope = -3.126169800000156e-05;
    return slope*T + intercept;
  }
  else { //if (plusx < 0.0)
    mooseError("lowT_line cannot use O/U less than 2.0");
  }
}

Real bulkC(Real OU, Real T){
  Real _kb = 8.617343e-5;
  Real plusx = OU - 2.0;
  if (plusx > 0) {
    return exp(-lowT_line(0,T)/(_kb*T)) + exp(-lowT_line(plusx,T)/(_kb*T));
  }
  else if (plusx == 0.0) {
    return exp(-lowT_line(0,T)/(_kb*T));
  }
  else if (plusx < 0){
    Real top = log10(exp(-lowT_line(0,T)/(_kb*T)) + exp(-lowT_line(abs(plusx),T)/(_kb*T)));
    Real mid = log10(exp(-lowT_line(0,T)/(_kb*T)));
    return pow(10,(mid - (top - mid)));
  }
  else {
    mooseError("bulkC stoichiometry invalid?");
  }
}

Real Efv_effective(Real OU, Real T){
  Real _kb = 8.617343e-5;
  return - _kb * T * log(bulkC(OU,T));
  }

Real gbC(Real OU, Real T, std::vector<Real> gb_se_data){
  Real _kb = 8.617343e-5;
  Real Efv_bulk = Efv_effective(OU, T);
  std::vector<Real> bulk_vector(gb_se_data.size());
  std::fill(bulk_vector.begin(), bulk_vector.end(), Efv_bulk);
  std::transform(gb_se_data.begin( ), gb_se_data.end( ), bulk_vector.begin( ), gb_se_data.begin( ),std::plus<double>( ));
  std::vector<Real> gb_conc_vector = gb_se_data;
  std::transform(gb_conc_vector.begin(), gb_conc_vector.end(), gb_conc_vector.begin(), [&](double &i){
        return exp(- i / (_kb * T));
    });
  std::transform(gb_conc_vector.begin(), gb_conc_vector.end(), gb_se_data.begin(), gb_conc_vector.begin(),
    [](auto const& a, auto const& b) {
      if (b <= 0.0)
        return 1.0;
      return a;
    }
  );
  return std::accumulate(gb_conc_vector.begin(), gb_conc_vector.end(), 0.0) / gb_conc_vector.size();
}



void
UO2CvMaterial::computeQpProperties()
{

  Real c_GB = 0.0;
  for (unsigned int i = 0; i < _gb_csv.size(); ++i)
  {
    c_GB += gbC(_OU[_qp], _T[_qp], _gb_data[i]);
  }
  c_GB = c_GB / _gb_csv.size();
  c_GB = std::min(_gb_max,c_GB);
  Real c_B = std::min(bulkC(_OU[_qp], _T[_qp]), _b_max);

  Real bounds = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
  {
    bounds += (*_vals[i])[_qp] * (*_vals[i])[_qp];
  }
  bounds += _phi[_qp] * _phi[_qp];

  // Solid phase equilibrium vacancy concentration calculation
  _cv_out[_qp] = c_B + 4.0 * (c_GB - c_B) * (1.0 - bounds) * (1.0 - bounds);
  // Derivatives wrt OPs
  for (unsigned int i = 0; i < _op_num; ++i)
  {
    (*_dcv[i])[_qp] = - 16 * (*_vals[i])[_qp] * (c_GB - c_B) * (1.0 - bounds);
    for (unsigned int j = i; j < _op_num; ++j)
    {
      (*_d2cv[i][j])[_qp] = 32 * (*_vals[i])[_qp] * (*_vals[j])[_qp] * (c_GB - c_B);
    }
  }
  // _cv_out[_qp] = std::accumulate(_gb_data[1].begin(), _gb_data[1].end(), 0.0) / _gb_data[1].size();

}
