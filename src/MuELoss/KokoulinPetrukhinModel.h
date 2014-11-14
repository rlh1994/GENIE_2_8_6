//____________________________________________________________________________
/*!

\class    genie::mueloss::KokoulinPetrukhinModel

\brief    Kokoulin-Petrukhin model for the energy loss of muons due to direct
          e+e- pair production.
          Concrete implementation of the MuELossI interface.

\ref      W.Lohmann, R.Kopp and R.Voss,
          Energy Loss of Muons in the Energy Range 1-10000 GeV, CERN 85-03

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  December 10, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KOKOULIN_PETRUKHIN_MODEL_H_
#define _KOKOULIN_PETRUKHIN_MODEL_H_

#include "MuELoss/MuELossI.h"
#include "Numerical/GSFunc.h"

namespace genie {

class IntegratorI;

namespace mueloss {

class KokoulinPetrukhinModel : public MuELossI
{
public:

  KokoulinPetrukhinModel();
  KokoulinPetrukhinModel(string config);
  virtual ~KokoulinPetrukhinModel();

  //! implement the MuELossI interface
  double        dE_dx   (double E, MuELMaterial_t material) const;
  MuELProcess_t Process (void) const { return eMupPairProduction; }

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
  const IntegratorI * fIntegrator;
};

} // mueloss namespace
} // genie   namespace

//____________________________________________________________________________
/*!
\class    genie::mueloss::KokoulinPetrukhinIntegrand

\brief    Auxiliary scalar function for the internal integration in Kokulin
          Petrukhin model

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  December 10, 2003
*/
//____________________________________________________________________________

namespace genie   {
namespace mueloss {

class KokoulinPetrukhinIntegrand : public GSFunc
{
public:
  KokoulinPetrukhinIntegrand(double E, double Z);
  ~KokoulinPetrukhinIntegrand();

  double operator () (const vector<double> & x);

private:
  double fE;
  double fZ;
};

} // mueloss namespace
} // genie   namespace

#endif  // _KOKOULIN_PETRUKHIN_MODEL_
