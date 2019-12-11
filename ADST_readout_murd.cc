//Libraries necesaries for script

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

#include <SDEvent.h>

//#define AUGER_RADIO_ENABLED
#include <RdEvent.h>
#include <RdRecShower.h>
#include "ShowerRRecDataQuantities.h"
#include "RdRecShowerParameterStorageMap.h"

//#define AUGER_AMIGA_ENABLED
//#include "/path/to/your/offline/Modules/General/RecDataWriterNG/MD2ADST.h"
//#include "/cr/data01/holt-e/offline_141017/v3r3/Modules/General/RecDataWriterNG/MD2ADST.h"

//#include "/cr/aera02/user/sarmiento/offline/trunk/Modules/General/RecDataWriterNG/MD2ADST.h"

#include <MDEventADST.h>
#include <MdSimCounter.h>

#include <RecEventFile.h>
#include <RecEvent.h>
#include <GenShower.h>
#include <Detector.h>
#include <Shower.h>
#include <DetectorGeometry.h>

#include <boost/algorithm/string.hpp>
#include "TF1.h"
#include "TMath.h"
#include "Math/DistFunc.h"

#include <sys/types.h>
#include <dirent.h>
#include <getopt.h>
#include <vector>

using namespace std;
//using namespace otoa;

//Read file ADST from path
void HelpMsg(char* av[]) {
   cerr << "Usage : "<< endl;
   cerr << "------------------------------------------------------------------------"<< endl;
   cerr << av[0] << " PATH_TO_ADST_FILES/YYYY/MM/ADST_YYYY_*.root OUTPUTNAMETAG"<< endl;
}


///main .....................................................................................

int main(int argc, char *argv[]) {
  map<int,int> m_Self_MC;

  //---------------------------------------------------------------
  //-------- read the adst files from the ascii input file --------
  //---------------------------------------------------------------
  if( argc < 2 ){ HelpMsg(argv); abort(); }

  cout << argv[1] << endl;
  RecEventFile *file = new RecEventFile(argv[1]);

  //Create a new RecEvent
  RecEvent *theevent = new RecEvent();
  file->SetBuffers(&theevent); //link the buffer

  ///--------------READ-DATA------------------------------------
  while(file->ReadNextEvent()==RecEventFile::eSuccess) {

  ///-----------------MC------------------------------------------------------

  const GenShower& MCShower = theevent->GetGenShower();
  //const MCEvent& MCEvent = theevent->GetMCEvent();
  double MCXmax = MCShower.GetXmaxInterpolated(); //Don't use GetXmax()!!!
  const TVector3 MCCore = MCShower.GetCoreUTMCS();
  const TVector3 MCAxisSiteCS = MCShower.GetAxisSiteCS();
  const TVector3 MCCoreSiteCS = MCShower.GetCoreSiteCS();
  double MCEnergy = MCShower.GetEnergy();
  double MCZenith = MCShower.GetZenith();
  double MCAzimuth = MCShower.GetAzimuth();

  string MCRunNr = theevent->GetEventId();
  MCRunNr = MCRunNr.substr(17,6);
  cout << "runnr " << MCRunNr << endl;

  ///----------------------SD-------------------------------------------------
  //    event ids
  const string EventID(theevent->GetAugerId());
  //UInt_t Date = theevent->GetYYMMDD();
  const SDEvent& SDEvent = theevent->GetSDEvent();
  const SdRecShower& sdrecshower = SDEvent.GetSdRecShower();
  double RecSuccessSD = SDEvent.HasLDF();
  double SdEnergy = sdrecshower.GetEnergy();
  const TVector3& SDCore = sdrecshower.GetCoreUTMCS();
  double SDZenith = TMath::RadToDeg()*sdrecshower.GetZenith();
  double SDAzimuth = TMath::RadToDeg()*sdrecshower.GetAzimuth();
  const TVector3& SdAxis = sdrecshower.GetPlaneFrontAxis();

  ///----------------------RADIO-------------------------------------------------
  // reconstruction success, zenith, azimuth, energy
  //  const RdEvent& RdEvent = theevent->GetRdEvent();
  const RdRecShower& rdrecshower = theevent->GetRdEvent().GetRdRecShower();
  double RecSuccessRd = rdrecshower.GetParameter(revt::eFitSuccess);
  double NrOfSignalStations = rdrecshower.GetParameter(revt::eNumberOfStationsWithPulseFound);
  double RecSuccessGeoRd = rdrecshower.GetParameter(revt::eFitSuccess);
  double GeoCeDXmax = rdrecshower.GetParameter(revt::eGeoCeLDFDxmax);
  double GeoCeDXmaxError = rdrecshower.GetParameterError(revt::eGeoCeLDFDxmax);
  double GeoCeErad = rdrecshower.GetParameter(revt::eGeoCeLDFErad);
  double GeoCeEradError = rdrecshower.GetParameterError(revt::eGeoCeLDFErad);
  double GeoCeLDFChi2 = rdrecshower.GetParameter(revt::eGeoCeLDFChi2);
  double GeoCeLDFNDF = rdrecshower.GetParameter(revt::eGeoCeLDFNDF);
  double pvalueGeoCe = ROOT::Math::chisquared_cdf_c(GeoCeLDFChi2,GeoCeLDFNDF, 0);

  //Define variables of Radio
  double zenithRd = -1.;
  double azimuthRd = -1.;
  double RdCoreX = -1.;
  double RdCoreY = -1.;
  double RdCoreZ = -1.;
  double radioEnergy = -1.;
  double radioEnergyError = -1.;
  double twoDLDFChi2 = -1;
  int twoDLDFNDF = -1;
  double RdEnergy = -1;
  double GeoAngle = -1.;
  double Alpha = -1.;
  double SE_RdLDF1 = -1.;
  double SE_RdLDF2 = -1.;
  double EE_RdLDF1 = -1.;
  double EE_RdLDF1error = -1.;
  double EE_RdLDF2 = -1.;
  double EE_RdLDF2error = -1.;
  double EE_RdLDF3 = -1.;
  double EE_RdLDF3error = -1.;
  double pvalue1dLDF = -1.;
  double EE_RdLDFNDF = -1.;
  /*double firtsAnt  = -1.;
  double firtsAmp = -1.;
  double antMaxAmp = -1.;
  double maxAmp = -1.;*/

  // fit results 2dldf fit
  double A = -1.;
  double Aerror = -1.;
  double sigma = -1.;
  double sigmaError = -1.;
  double core_x = -1.;
  double core_y = -1.;
  double chi2 = -1.;
  double ndf = -1.;
  double OnedLDFChi2 = -1;
  bool val2dLDFFit = 0;
  TFitResult fitresult = rdrecshower.Get2dLDFFitResult();
  val2dLDFFit = fitresult.IsValid();
  cout << "val2dLDFFit = " << val2dLDFFit << endl;

  TVector3 RdAxis = TVector3(0.,0.,0.);

  if ( RecSuccessRd==1 ){
    zenithRd = TMath::RadToDeg()*rdrecshower.GetZenith();
    azimuthRd = TMath::RadToDeg()*rdrecshower.GetAzimuth();
    // core position
    RdCoreX = rdrecshower.GetParameter(revt::eCoreX);
    RdCoreY = rdrecshower.GetParameter(revt::eCoreY);
    RdCoreZ = rdrecshower.GetParameter(revt::eCoreZ);

    const std::vector<TF1>& LDFs = rdrecshower.GetRdLDF();
    const TString Name=LDFs[0].TF1::GetExpFormula();
    SE_RdLDF1 = LDFs[0].TF1::GetParameter(0);
    SE_RdLDF2 = LDFs[0].TF1::GetParameter(1);
    EE_RdLDF1 = LDFs[1].TF1::GetParameter(0);
    EE_RdLDF1error = LDFs[1].TF1::GetParError(0);
    EE_RdLDF2 = LDFs[1].TF1::GetParameter(1);
    EE_RdLDF2error = LDFs[1].TF1::GetParError(1);
    EE_RdLDF3 = LDFs[1].TF1::GetParameter(2);
    EE_RdLDF3error = LDFs[1].TF1::GetParError(2);
    OnedLDFChi2 = LDFs[0].TF1::GetChisquare();
    EE_RdLDFNDF = LDFs[1].TF1::GetNDF();
    pvalue1dLDF = ROOT::Math::chisquared_cdf_c(OnedLDFChi2, EE_RdLDFNDF, 0);

    if ( rdrecshower.HasParameter( revt::eRadiationEnergy ) ){
      radioEnergy = rdrecshower.GetParameter(revt::eRadiationEnergy);
      radioEnergyError = rdrecshower.GetParameterError(revt::eRadiationEnergy);
      twoDLDFChi2 = rdrecshower.GetParameter(revt::e2DLDFfitSignalChi2);
      twoDLDFNDF = rdrecshower.GetParameter(revt::e2DLDFfitSignalNDF);
      cout << "radiation energy reconstructed to " << radioEnergy << " +/- " << radioEnergyError << " with Chi2 = " << twoDLDFChi2 << " and NDF = " << twoDLDFNDF << endl;
      RdEnergy = rdrecshower.GetParameter(revt::eEnergy);
      cout << "primary energy reconstructed from radio is " << RdEnergy << endl;

      A = fitresult.Parameter(0);
      Aerror = fitresult.ParError(0);
      sigma = fitresult.Parameter(3);
      sigmaError = fitresult.ParError(3);
      core_x = fitresult.Parameter(1);
      core_y = fitresult.Parameter(2);
      chi2 = fitresult.Chi2();
      ndf = fitresult.NFreeParameters();
      // chi2 of fitresult uses all stations of the fit. Thereby stations without a signal are considered as well, which do not provide information to the fit. This minimizes the chi2 artificially and its underestimated.
      
      cout << "fit result of the 2D LDF fit:\t\t A = " << A << " ±" << Aerror <<"\t\t sigma = " << sigma << "  " << sigmaError << "\t\t x = " << core_x << "\t\t y = " << core_y << "\t\t chi2 = " << chi2 << "\t\t ndf = " << ndf << "\t\t chi2/ndf = " << chi2/ndf << "\t\t DXmax = " << GeoCeDXmax <<endl;

    } else {
      cout << "radiation energy not reconstructed" << endl;
    }

/*      stringstream filenameAntennas("");
      filenameAntennas << "Antennas_" << argv[2] << ".txt";
      bool FileAntennasExists = false;
      if (ifstream(filenameAntennas.str().c_str())) FileAntennasExists = true;
      ofstream outfileAntennas;
      outfileAntennas.open(filenameAntennas.str().c_str(),ios::app);
      filenameAntennas.clear();
      if(!FileAntennasExists){
        outfileAntennas
               << "antenna id\tdistance MC core\tdistance MC axis\tmax amplitude\tintegrated signal\tlorentz angle\tMC energy"
           << endl;
  }

  const std::vector<RdRecStation> RdStations = theevent->GetRdEvent().GetRdStationVector();

  DetectorGeometry* detGeom = new DetectorGeometry();
  file->ReadDetectorGeometry(*detGeom);
  vector<double> maxAmplitudeMag(RdStations.size());
  vector<double> distAntennaAxis(RdStations.size());
  vector<double> signal(RdStations.size());

  for (unsigned int s=0; s < RdStations.size(); ++s ) {
    int e = RdStations.size();

    const RdRecStation& rdrecstation = RdStations[s];

    if ( !rdrecstation.HasPulse() || rdrecstation.IsSaturated() || rdrecstation.IsRejected() ) {
      //cout << "Antenna is rejected. Continue..." << endl;
      continue;
    }

    int sId = rdrecstation.GetId();
    //cout << "antenna id " << sId << endl;
    const TVector3 antennaPosition = detGeom->GetRdStationPosition(sId);

    //project the vector on shower direction
    double projVec = - MCAxisSiteCS[0] * antennaPosition[0]
                      - MCAxisSiteCS[1] * antennaPosition[1]
                      - MCAxisSiteCS[2] * antennaPosition[2];

    double distAntennaCore = (antennaPosition - MCCoreSiteCS).Mag();
    double distAntennaAxiss = detGeom->GetRdStationAxisDistance(sId, MCAxisSiteCS, MCCoreSiteCS);
    maxAmplitudeMag[s] = rdrecstation.GetParameter(revt::ePeakAmplitudeMag);
    double maxAmplitudeMagg = rdrecstation.GetParameter(revt::ePeakAmplitudeMag);
    //double integratedSignal = rdrecstation.GetParameter(revt::eIntegratedSignalMag);
    double integratedSignal = rdrecstation.GetParameter(revt::eSignal);
    signal[s] = rdrecstation.GetParameter(revt::eSignal);
    distAntennaAxis[s] = detGeom->GetRdStationAxisDistance(sId, MCAxisSiteCS, MCCoreSiteCS);

    outfileAntennas
           << sId
           << "\t"
           << distAntennaCore
           << "\t"
           << distAntennaAxiss
           << "\t"
           << maxAmplitudeMagg
           << "\t"
           << integratedSignal
           << "\t"
           << Alpha
           << "\t"
           << MCEnergy
           << endl;

    }

  vector<double> distancias;
  double maxAmp = 0;
  double antMaxAmp=0;
  for(unsigned int ii=0; ii < signal.size(); ii++){
    if (maxAmp < signal.at(ii)){
       maxAmp = signal.at(ii);
       antMaxAmp = distAntennaAxis.at(ii);
    }
  }

  double firtsAmp = std::numeric_limits<double>::infinity();
  double firtsAnt = std::numeric_limits<double>::infinity();
  for(unsigned int jj=0; jj < distAntennaAxis.size(); jj++){
    if ( (firtsAnt > distAntennaAxis.at(jj) ) && (signal.at(jj) > 0) ){
      firtsAnt = distAntennaAxis.at(jj);
      firtsAmp = signal.at(jj);
    }
  }

  cout << firtsAnt << " VVV " << firtsAmp << " "  << endl;
  cout << antMaxAmp << " KKK " << maxAmp << endl;
*/

  TVector3 MagFieldVec = rdrecshower.GetMagneticFieldVector();

  RdAxis = TVector3(rdrecshower.GetParameter(revt::eShowerAxisX), rdrecshower.GetParameter(revt::eShowerAxisY), rdrecshower.GetParameter(revt::eShowerAxisZ));
  double sineAlpha = RdAxis.Cross(MagFieldVec).Mag()/(RdAxis.Mag()*MagFieldVec.Mag());
  Alpha = TMath::RadToDeg()*std::asin(sineAlpha);
  }
  else{ cout << "Radio not successfully reconstructed!" << endl; }

  ///----------------------AMIGA-------------------------------------------------
  MDEvent& MdEvent = theevent->GetMDEvent();
  const MdRecShower& mdrecshower = theevent->GetMDEvent().GetMdRecShower();
  //double MDLDFSuccess = mdrecshower.GetMLDFRecStage();
  double MDLDFSuccess = mdrecshower.IsLdfReconstructed();

  // data from MdRecShower
  double MLDFChi2 = mdrecshower.GetMLDFChi2();
  double MLDFNdof = mdrecshower.GetMLDFNdof();
  double MLDFLikelihood = mdrecshower.GetMLDFLikelihood();
  //double MAlpha = mdrecshower.GetAlpha(); // fixed to 1
  //double MAlphaError = mdrecshower.GetAlphaError();
  double MBeta = mdrecshower.GetBeta();
  double MBetaError = mdrecshower.GetBetaError();
  double MBetaSystematics = mdrecshower.GetBetaSystematics();
  //double MGamma = mdrecshower.GetGamma(); // fixed to 1.85
  //double MGammaError = mdrecshower.GetGammaError();
  //double MR0 = mdrecshower.GetR0(); // fixed to 150
  //double MR0Error = mdrecshower.GetR0Error();
  double NMuRef = mdrecshower.GetNMuRef();
  double NMuRefError = mdrecshower.GetNMuRefError();
  //double NMuRefSystematics = mdrecshower.GetNMuRefSystematics(); // immer 0, nicht gefüllt
  double ReferenceDistance = mdrecshower.GetReferenceDistance();
  //const TVector3& Barycenter = mdrecshower.GetBarycenter();

  cout << "Md LDF reconstructed with reference distance of " << ReferenceDistance << ": beta = " << MBeta << " ± " << MBetaError << " ± " << MBetaSystematics << ", Nmu = " << NMuRef << " ± " << NMuRefError << endl;

  if (!MDLDFSuccess)
    cout << "MD Reconstruction not successful!" << endl;

///----------------------------------------------------------------------------------------------------
///write all parameters in one outputfile--------------------------------------------------------------
///----------------------------------------------------------------------------------------------------

///write out to AllEvents file----------------------------------

  stringstream filenameAllEvents("");
  filenameAllEvents << "AllEvents_" << argv[2] << ".txt";
  bool FileAllEventsExists = false;
  if (ifstream(filenameAllEvents.str().c_str())) FileAllEventsExists = true;
  ofstream outfileAllEvents;
  outfileAllEvents.open(filenameAllEvents.str().c_str(),ios::app);
  filenameAllEvents.clear();

  if(!FileAllEventsExists){
    outfileAllEvents
           << "# All events, no cuts applied after Offline reconstruction" << endl << "#" << endl << "#" << endl
           << "# "
           << "runnr\tMC energy\tMC zenith\tMC azimuth\tMC Xmax\tMC core x\tMC core y\tMC core z\t"
           << "Sd rec success\tSd energy\tSD zenith\tSD azimuth\tSD core x\tSD core y\tSD core z\t"
           << "Rd rec success\tStations with Pulse\tRd zenith\tRd azimuth\tRd core x\tRd core y\tRd core z\t"
	   << "rad energy\trad energy error\tlorentz angle\t2dLDF sigma\t2dLDF sigma error\t2dLDF Chi2\t2dLDf NDF\t"
	   << "Rd energy\tRd SE_LDF 1\tRd SE_LDF 2\t"
	   << "SE_RdLDF1\tSE_RdLDF2\tEE_RdLDF1\tEE_RdLDF1error\tEE_RdLDF2\tEE_RdLDF2error\tEE_RdLDF3\tOnedLDFChi2\tEE_RdLDFNDF\tpvalue1dLDF\t"
	   << "GeoCe DXmax\tGeoCe DXmax Error\tGeoCe Erad\tGeoCe Erad Error\tGeoCe LDFChi2\tGeoCe LDF ndf\tpvalue GeoCeLDF\t"
           << "Md rec success\tN_mu_ref\tN_mu_ref error\tMLDF Chi2\tMLDF NDF\tMLDF Likelihood\tM beta\tM beta error\tM beta syst"
           << endl;
  }

  outfileAllEvents
                   << MCRunNr
           << "\t" << MCEnergy
           << "\t" << MCZenith*(180./M_PI)
           << "\t" << MCAzimuth*(180./M_PI)
           << "\t" << MCXmax
           << "\t" << MCCore[0]
           << "\t" << MCCore[1]
           << "\t" << MCCore[2]

           << "\t" << RecSuccessSD
           << "\t" << SdEnergy
           << "\t" << SDZenith
           << "\t" << SDAzimuth
           << "\t" << SDCore[0]
           << "\t" << SDCore[1]
           << "\t" << SDCore[2]

           << "\t" << RecSuccessRd
	   << "\t" << NrOfSignalStations
           << "\t" << zenithRd
           << "\t" << azimuthRd
           << "\t" << RdCoreX
	   << "\t" << RdCoreY
	   << "\t" << RdCoreZ
	   << "\t" << radioEnergy
           << "\t" << radioEnergyError
           << "\t" << Alpha
           << "\t" << sigma
           << "\t" << sigmaError
           << "\t" << twoDLDFChi2
           << "\t" << twoDLDFNDF
           << "\t" << RdEnergy

           << "\t" << SE_RdLDF1
    	   << "\t" << SE_RdLDF2
	   << "\t" << EE_RdLDF1
           << "\t" << EE_RdLDF1error
	   << "\t" << EE_RdLDF2
	   << "\t" << EE_RdLDF2error
	   << "\t" << EE_RdLDF3
	   << "\t" << EE_RdLDF3error
           << "\t" << OnedLDFChi2
           << "\t" << EE_RdLDFNDF
           << "\t" << pvalue1dLDF
           /*<< "\t" << firtsAnt
           << "\t" << firtsAmp
           << "\t" << antMaxAmp
           << "\t" << maxAmp*/

           << "\t" << GeoCeDXmax
           << "\t" << GeoCeDXmaxError
           << "\t" << GeoCeErad
           << "\t" << GeoCeEradError
           << "\t" << GeoCeLDFChi2
           << "\t" << GeoCeLDFNDF
	   << "\t" << pvalueGeoCe

	   << "\t" << MDLDFSuccess
	   << "\t" << NMuRef
           << "\t" << NMuRefError
           << "\t" << MLDFChi2
           << "\t" << MLDFNdof
           << "\t" << MLDFLikelihood
           << "\t" << MBeta
           << "\t" << MBetaError
           << "\t" << MBetaSystematics

           << endl;

  outfileAllEvents.close();


///-----------------------------------CUTS----------------------------------
  
  stringstream cuts("");
  cuts << "# Cuts applied after Offline reconstruction:" << endl << "# ";

  ///radio------------------------------

  cout << "#### LDF1D parameteres ###" << endl;
  cout << "eps      eta     epsExt      a       b     chi2" << endl;
  cout << " " << SE_RdLDF1 << " " << SE_RdLDF2 << " " << EE_RdLDF1 << " " << EE_RdLDF1error << " " << EE_RdLDF2 << " " << EE_RdLDF3 << " " << OnedLDFChi2 << endl;
  cout << "###############" << endl;


  if ( !RecSuccessRd ) continue;
  cuts << "RecSuccessRd";
  if ( !rdrecshower.HasParameter( revt::eRadiationEnergy ) ) continue;
  cuts << ", radiation energy reconstructed";

  if (radioEnergy<=0.) {
    cout << "radio energy smaller than zero!" << endl << endl;
    continue;
  }
  cuts << ", E_Rd > 0";

  double diffAngleSdRd = TMath::RadToDeg()*std::asin(RdAxis.Cross(SdAxis).Mag()/(RdAxis.Mag()*SdAxis.Mag()));
  if( diffAngleSdRd > 5. ) {
    cout << "angle between Sd and Rd axis = " << diffAngleSdRd << endl << endl;
    continue;
  }
  cuts << ", diff reconstructed angle of Sd and Rd < 20";

  if( Alpha < 10. ) {
    cout << "alpha < 10." << endl << endl;
    continue;
  }
  cuts << ", alpha > 10";

  if(pvalue1dLDF < 0.05)
    continue;

  if( sigma > 270. ) {
    cout << "sigma higher than 270!" << endl << endl;
    continue;
  }
  cuts << ", sigma < 270";

  if( sigma < 60. ) {
    cout << "sigma smaller than 60!" << endl << endl;
    continue;
  }
  cuts << ", sigma > 60";

  if( sigmaError/sigma > 0.4 ) {
    cout << "sigmaError/sigma larger than 0.4" << endl << endl;
    continue;
  }
  cuts << ", relative uncertainty of sigma < 0.4";

/*
  if( EE_RdLDF1 < 0.00015 ) continue;
  cuts << ", epsilon smaller than 0.15E-03!";

    if( MCEnergy>5E18 ) {
    cout << "Sd energy < 10^17.5" << endl << endl;
    continue;
    }

  
  double chi2divbyndf = twoDLDFChi2 / twoDLDFNDF;
  double cutchi2 = 10.;
  if( chi2divbyndf > cutchi2 ) {
    cout << "Chi2 / NDF larger than " << cutchi2 << "!" << endl << endl;
    continue;
  }
  cuts << ", 2dLDF chi2/NDF < 10";
  */

  ///SD---------------------------------
  if( MCZenith*(180./M_PI)>60 ) {
    cout << " zenith angle > 60 degrees" << endl << endl;
    continue;
  }
  cuts << ", MC zenith < 60°";
  /*if( SdEnergy<17.5 ) {
    cout << "Sd energy < 10^17.5" << endl << endl;
    continue;
  }
  cuts << ", Sd energy > 10^17.5 eV";
  */


  ///MD---------------------------------
  if ( !MDLDFSuccess ) continue;
  cuts << ", MDLDFSuccess";

  //Cut on Unitary Cell
  double UniCellX[7]={450617.51, 449879.54, 449499.61, 449874.21, 450626.13, 450994.47, 450248.61};		//last element is center
  double UniCellY[7]={6113923.49, 6113943.05, 6114580.2, 6115223.69, 6115227.11, 6114576.65, 6114571.74};	//last element is center

  double m[6][6] = {};
  double b[6][6] = {};
  for ( int i=0; i<6; i++ ) {
    for ( int j=0; j<6; j++ ) {
      if ( i != j ) {
        m[i][j]=(UniCellY[i]-UniCellY[j])/(UniCellX[i]-UniCellX[j]);
        b[i][j]=UniCellY[i]-UniCellX[i]*m[i][j];
      }
      else if ( i==j ) {
        m[i][j]=0;
        b[i][j]=0;
      }
    }
  }

  bool InUniCell=false;

  if (UniCellX[1]<=SDCore[0] and SDCore[0]<=UniCellX[0] and UniCellY[1]<=SDCore[1] and SDCore[1]<=UniCellY[3]) InUniCell=true;
  else if (UniCellX[3]<=SDCore[0] and SDCore[0]<=UniCellX[1] and UniCellY[2]<=SDCore[1] and SDCore[1]<=UniCellY[3]) InUniCell=true;
  else if (UniCellX[0]<=SDCore[0] and SDCore[0]<=UniCellX[4] and UniCellY[5]<=SDCore[1] and SDCore[1]<=UniCellY[4]) InUniCell=true;

  else if (UniCellX[2]<=SDCore[0] and SDCore[0]<=UniCellX[3] and UniCellY[2]<=SDCore[1] and SDCore[1]<=(m[2][3]*SDCore[0]+b[2][3])) InUniCell=true;
  else if (UniCellX[2]<=SDCore[0] and SDCore[0]<=UniCellX[1] and UniCellY[2]>=SDCore[1] and SDCore[1]>=(m[2][1]*SDCore[0]+b[2][1])) InUniCell=true;
  else if (UniCellX[3]<=SDCore[0] and SDCore[0]<=UniCellX[4] and UniCellY[3]<=SDCore[1] and SDCore[1]<=(m[3][4]*SDCore[0]+b[3][4])) InUniCell=true;
  else if (UniCellX[4]<=SDCore[0] and SDCore[0]<=UniCellX[5] and UniCellY[5]<=SDCore[1] and SDCore[1]<=(m[4][5]*SDCore[0]+b[4][5])) InUniCell=true;
  else if (UniCellX[1]<=SDCore[0] and SDCore[0]<=UniCellX[0] and UniCellY[1]>=SDCore[1] and SDCore[1]>=(m[1][0]*SDCore[0]+b[1][0])) InUniCell=true;
  else if (UniCellX[0]<=SDCore[0] and SDCore[0]<=UniCellX[5] and UniCellY[5]>=SDCore[1] and SDCore[1]>=(m[0][5]*SDCore[0]+b[0][5])) InUniCell=true;

  /*if (InUniCell==false) {
    cout << "core not in unitary cell!" << endl << endl;
    //continue;
  }
  cuts << ", core inside AMIGA unitary cell";
  */

  cout << endl << "# cuts here!" << endl;

/// write to CutEvents file------------------

  stringstream filenameCutEvents("");
  filenameCutEvents << "CutEvents_" << argv[2] << ".txt";
  bool FileCutEventsExists = false;
  if (ifstream(filenameCutEvents.str().c_str())) FileCutEventsExists = true;
  ofstream outfileCutEvents;
  outfileCutEvents.open(filenameCutEvents.str().c_str(),ios::app);
  filenameCutEvents.clear();
  if(!FileCutEventsExists){
    outfileCutEvents
           << cuts.str()
           << "#" 
           << "runnr\tMC energy\tMC zenith\tMC azimuth\tMC Xmax\tMC core x\tMC core y\tMC core z\t"
           << "Sd rec success\tSd energy\tSD zenith\tSD azimuth\tSD core x\tSD core y\tSD core z\t"
           << "Rd rec success\tStations with Pulse\tRd zenith\tRd azimuth\tRd core x\tRd core y\tRd core z\t"
           << "rad energy\trad energy error\tlorentz angle\t2dLDF sigma\t2dLDF sigma error\t2dLDF Chi2\t2dLDf NDF\t"
           << "Rd energy\tRd SE_LDF 1\tRd SE_LDF 2\t"
           << "SE_RdLDF1\tSE_RdLDF2\tEE_RdLDF1\tEE_RdLDF1error\tEE_RdLDF2\tEE_RdLDF2error\tEE_RdLDF3\tOnedLDFChi2\tEE_RdLDFNDF\tpvalue1dLDF\t"
           << "GeoCe DXmax\tGeoCe DXmax Error\tGeoCe Erad\tGeoCe Erad Error\tGeoCe LDFChi2\tGeoCe LDF ndf\tpvalue 2dLDF\t"
           << "Md rec success\tN_mu_ref\tN_mu_ref error\tMLDF Chi2\tMLDF NDF\tMLDF Likelihood\tM beta\tM beta error\tM beta syst"
           << endl;
  }
      	outfileCutEvents
                   << MCRunNr
           << "\t" << MCEnergy
           << "\t" << MCZenith*(180./M_PI)
           << "\t" << MCAzimuth*(180./M_PI)
           << "\t" << MCXmax
           << "\t" << MCCore[0]
           << "\t" << MCCore[1]
           << "\t" << MCCore[2]

           << "\t" << RecSuccessSD
           << "\t" << SdEnergy
           << "\t" << SDZenith
           << "\t" << SDAzimuth
           << "\t" << SDCore[0]
           << "\t" << SDCore[1]
           << "\t" << SDCore[2]

           << "\t" << RecSuccessRd
           << "\t" << NrOfSignalStations
           << "\t" << zenithRd
           << "\t" << azimuthRd
           << "\t" << RdCoreX
           << "\t" << RdCoreY
           << "\t" << RdCoreZ
           << "\t" << radioEnergy
           << "\t" << radioEnergyError
           << "\t" << Alpha
           << "\t" << sigma
           << "\t" << sigmaError
           << "\t" << twoDLDFChi2
           << "\t" << twoDLDFNDF
           << "\t" << RdEnergy

           << "\t" << SE_RdLDF1
           << "\t" << SE_RdLDF2
           << "\t" << EE_RdLDF1
           << "\t" << EE_RdLDF1error
           << "\t" << EE_RdLDF2
           << "\t" << EE_RdLDF2error
           << "\t" << EE_RdLDF3
           << "\t" << EE_RdLDF3error
           << "\t" << OnedLDFChi2
           << "\t" << EE_RdLDFNDF
           << "\t" << pvalue1dLDF
           /*<< "\t" << firtsAnt
           << "\t" << firtsAmp
           << "\t" << antMaxAmp
           << "\t" << maxAmp*/

           << "\t" << GeoCeDXmax
           << "\t" << GeoCeDXmaxError
           << "\t" << GeoCeErad
           << "\t" << GeoCeEradError
           << "\t" << GeoCeLDFChi2
           << "\t" << GeoCeLDFNDF
	   << "\t" << pvalueGeoCe

           << "\t" << MDLDFSuccess
           << "\t" << NMuRef
           << "\t" << NMuRefError
           << "\t" << MLDFChi2
           << "\t" << MLDFNdof
           << "\t" << MLDFLikelihood
           << "\t" << MBeta
           << "\t" << MBetaError
           << "\t" << MBetaSystematics

           << endl;

      outfileCutEvents.close();
   }//while ReadNextEvent()

  exit(0);
}

