//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 05, 2008 - CA
   In vrs 2.3.1, most of the driver implementation code was factored out 
   to the new GAtmoFlux base class so as to be shared by the functionally 
   identical GBartolAtmoFlux driver
 @ Feb 05, 2008 - Chris Backhouse (Oxford)
   Fixed a bug in bin definitions (the TH2 constructor takes an array of 
   bin lower edges, but the bin centres were being passed in instead).
 @ Feb 23, 2010 - CA
   Build bin arrays at ctor. Re-structuring and clean-up.

*/
//____________________________________________________________________________

#include <cassert>
#include <fstream>

#include <TH2D.h>
#include <TMath.h>

#include "FluxDrivers/GHondaAtmoFlux.h"
#include "Messenger/Messenger.h"

using std::ifstream;
using std::ios;
using namespace genie;
using namespace genie::flux;

//____________________________________________________________________________
GHondaAtmoFlux::GHondaAtmoFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
       << "Instantiating the Fluka-3D atmospheric neutrino flux driver";

  this->SetBinSizes();
  this->Initialize();
}
//___________________________________________________________________________
GHondaAtmoFlux::~GHondaAtmoFlux()
{

}
//___________________________________________________________________________
void GHondaAtmoFlux::SetBinSizes(void)
{
// Generate the correct cos(theta) and energy bin sizes
// The flux is given in 40 bins of cos(zenith angle) from -1.0 to 1.0
// (bin width = 0.05) and 61 equally log-spaced energy bins (20 bins 
// per decade), with Emin = 0.100 GeV.
//

  fCosThetaBins  = new double [kGHondaNumCosThetaBins  + 1];
  fEnergyBins    = new double [kGHondaNumLogEvBins     + 1];

  double dcostheta = 
      (kGHondaCosThetaMax - kGHondaCosThetaMin) /
      (double) kGHondaNumCosThetaBins;

  double logEmax = TMath::Log10(1.);
  double logEmin = TMath::Log10(kGHondaEvMin);
  double dlogE = 
      (logEmax - logEmin) / 
      (double) kGHondaNumLogEvBinsPerDecade;

  for(unsigned int i=0; i<= kGHondaNumCosThetaBins; i++) {
     fCosThetaBins[i] = kGHondaCosThetaMin + i * dcostheta;
     if(i != kGHondaNumCosThetaBins) {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: CosTheta bin " << i+1 
         << ": lower edge = " << fCosThetaBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: CosTheta bin " << kGHondaNumCosThetaBins 
         << ": upper edge = " << fCosThetaBins[kGHondaNumCosThetaBins];
     }
  }

  for(unsigned int i=0; i<= kGHondaNumLogEvBins; i++) {
     fEnergyBins[i] = TMath::Power(10., logEmin + i*dlogE);
     if(i != kGHondaNumLogEvBins) {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: Energy bin " << i+1 
         << ": lower edge = " << fEnergyBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: Energy bin " << kGHondaNumLogEvBins 
         << ": upper edge = " << fEnergyBins[kGHondaNumLogEvBins];
     }
  }

  for(unsigned int i=0; i< kGHondaNumLogEvBins; i++) {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: Energy bin " << i+1
         << ": bin centre = " << (fEnergyBins[i] + fEnergyBins[i+1])/2.;
  }

  fNumCosThetaBins = kGHondaNumCosThetaBins;
  fNumEnergyBins   = kGHondaNumLogEvBins;
}
//____________________________________________________________________________
bool GFlukaAtmo3DFlux::FillFluxHisto2D(TH2D * histo, string filename)
{
  LOG("Flux", pNOTICE) << "Loading: " << filename;

  if(!histo) {
     LOG("Flux", pERROR) << "Null flux histogram!";
     return false;
  }
  ifstream flux_stream(filename.c_str(), ios::in);
  if(!flux_stream) {
     LOG("Flux", pERROR) << "Could not open file: " << filename;
     return false;
  }

  int    ibin;
  double energy, costheta, flux;
  char   j1, j2;

  double scale = 1.0; // 1.0 [m^2], OR 1.0e-4 [cm^2]

  while ( !flux_stream.eof() ) {
    flux = 0.0;
    flux_stream >> energy >> j1 >> costheta >> j2 >> flux;
    if( flux>0.0 ){
      LOG("Flux", pINFO)
        << "Flux[Ev = " << energy 
        << ", cos8 = " << costheta << "] = " << flux;
      // note: reversing the Fluka sign convention for zenith angle
      ibin = histo->FindBin( (Axis_t)energy, (Axis_t)(-costheta) );   
      histo->SetBinContent( ibin, (Stat_t)(scale*flux) );
    }
  }
  return true;
}
//___________________________________________________________________________
