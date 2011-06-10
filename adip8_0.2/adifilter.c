#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include "fftw3.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TMath.h"
#include "radiation.h"
#include "paramanage.h"

using namespace std;
struct t_parameter parameter;
//reads in fft_% Tree from existing root file, filters, mixes, returns to time domain
int main(int argc,char *argv[])
{
  TApplication theApp("App",&argc,argv);
  TString cardname = theApp.Argv(1);

  //use Adipark parameter-reading routine to read control card
  char *parameter_filename = new char[255];
  strcpy(parameter_filename,cardname.Data());
  init_parameters();
  int parameter_file_found = 0;
  if (load_init_data(parameter_filename)==0)
  {

    char *parameter_filename_sub = new char[255];
    cout << "Try parameter file w/o suffix: ";
    strncat(parameter_filename_sub,parameter_filename,strlen(parameter_filename)-4);
    cout << parameter_filename_sub << endl;
    if (load_init_data(parameter_filename_sub) == 0)
    {
      cout << "Try parameter file w/o point: ";
      strcpy(parameter_filename_sub,"");
      strncat(parameter_filename_sub,parameter_filename,strlen(parameter_filename)-1);
      cout << parameter_filename_sub << endl;
      if (load_init_data(parameter_filename_sub) == 0)
            {
              cout << endl << "ERROR: missing some settings in INI-file!";
              cout << endl << endl << flush;
            }
      else
      {
        parameter_file_found = 1;
              cardname.Resize(strlen(parameter_filename)-1);
      }
    }
    else
    {
      parameter_file_found = 1;
      cardname.Resize(strlen(parameter_filename)-4);
    }
  }
  else  parameter_file_found = 1;
  
  double kT = (parameter.antenna_temp)*1.38e-8;//fJ, T=35K, amplifier Teff = 25 K
  Bool_t addNoise = kFALSE;
  if (kT > 0) addNoise = kTRUE;
  init_data(); //set default tl_data
  switch (parameter.rad_calc_mode)
  {
    case 2://parallel wires
           init_tl_data(kFALSE);
           break;
    case 3://square wg
           init_sq_wg_data(27.05e3*2*M_PI/c);
           //Warning!  Wave imp. for TE modes freq dep, not implemented properly
           break;
    case 4://circular wg
           init_circ_wg_data(27.05e3*2*M_PI/c);
           //Warning!  Wave imp. for TE modes freq dep, not implemented properly
           break;
    case 5://parallel plates/strips
           init_pp_data();
           break;
    case 6://coaxial cable
           init_coax_data();
           break;
    case 7://offset parallel wires
           init_tl_data(kTRUE);
           break;
  }

  TString rootname = cardname+TString(".root");
  TFile* tfin = new TFile(rootname,"read");
  
  TTree *runcard = (TTree*)tfin->Get("runcard");
  Int_t nEvents = runcard->GetEntries();
  Bool_t printDAT = kFALSE;
  Bool_t extendWF = kTRUE;

  //set root's global rand gen to Mersene Twister
  TRandom3 *r3 = new TRandom3(0);
  gRandom = r3;
  
  //output file, no extra points
  TFile* tfout = new TFile(cardname+"_filtered.root","recreate");
  //time-domain antenna data, post-filtering and mixing
  typedef struct {
    Long_t i;
    Double_t t, Ef, fW_t, nfW_t, sig, noise, vtot;
  } ANTINFO;
  static ANTINFO anti;//short 
  static ANTINFO wfi;//with added time points, like real waveform
  //fft output (after fileter) 
  typedef struct {
    Double_t Hz, outr, outi, fW_f, fJ_f, fJpHz;//from fftw
  } POWERCALC;
  static POWERCALC pc;//short
  static POWERCALC npc;//long time series with noise
  
 //file for Gray and Joe  with exteneded waveform
 TFile* tfex = new TFile(cardname+"_extended.root","recreate");
  TNtuple* filter = new TNtuple("filter","filter parameters","i:LO:SF:sigStart");
  float ntarray[4];
  
  for ( Int_t event = 1; event <=nEvents; event++ ) {
    //read in data (signal only)
    TTree *fft_in = (TTree*)tfin->Get(Form("fft_%d",event));
    //noise already added
    //TTree *fft_in = (TTree*)tfin->Get(Form("nfft_%d",event));
    Int_t nEntries = fft_in->GetEntries();
    cout <<  endl;
    cout <<  endl;
    cout << cardname << " open with " << nEntries << " freq. bins " << endl;

    //new trees
    tfout->cd();
    TTree* antTree= new TTree(Form("ant_%d",event),"filtered time domain results");
    antTree->Branch("anti", &anti, "i/L:t/D:Ef/D:fW_t/D:nfW_t/D:sig/D:noise/D:vtot/D");
    TTree* fftTree = new TTree(Form("fft_%d",event),"fft results after filter");
    fftTree->Branch("pc", &pc, "Hz/D:outr/D:outi/D:fW_f/D:fJ_f/D:fJpHz/D");

    tfex->cd();
    TTree* wfTree = new TTree(Form("wf_%d",event),"Extended time domain results with noise");
    wfTree->Branch("wfi", &wfi, "i/L:t/D:Ef/D:fW_t/D:nfW_t/D:sig/D:noise/D:vtot/D");
    TTree* nfftTree = new TTree(Form("nfft_%d",event),"freq domain results with extra samples");
    nfftTree->Branch("npc", &npc, "Hz/D:noutr/D:nouti/D:nfW_f/D:nfJ_f/D:nfJpHz/D");
 
  
    //initialize backward fft, take N/2+1 freq. points, create N time points 
    int N = 2*(nEntries-1); //max number of time points for this sim
    double *out_td;
    fftw_complex *in;
    fftw_plan p;
    out_td = (double*) fftw_malloc(sizeof(double) * N);
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
 
    //iniitalize other arrays
    int nF = 0;//nF=N/2+1
    Double_t freq = 0, fMax = 0;
    Double_t LO = 26.9e9;//Total Local oscillator, in Hz, minimum 24.5, max 26.4?
    Double_t SF = 6e8;//Sampling Freq, in Hz, max is 2.5e9, bandwidth is SF/2, goal is 0.2 GHz
    cout << "Reading in power spectrum: LO " << LO*1e-9 << " GHz, SF: " << SF*1e-6 << " MHz "<<endl;
    for(int i=0;i<nEntries;i++)
    {
      fft_in->GetEntry(i);
      freq = fft_in->FindLeaf("Hz")->GetValue();
      //this is mixing and filtering section
      if ( freq>=LO&&(freq-LO)<=SF/2 ) {
        in[nF][0] = fft_in->FindLeaf("outr")->GetValue();
        in[nF][1] = fft_in->FindLeaf("outi")->GetValue();
        nF++;//number of freq. points within filter
        fMax = freq-LO;
      }
    }
    cout << "intial freq points " << nEntries << " and time points " << N <<  endl; 
    cout << " filtered freq points: " << nF << " and time points 2*(nF-1) " << 2*(nF-1) << endl;
    int max = TMath::Min(N,2*(nF-1));
   //from freq. to time
    p = fftw_plan_dft_c2r_1d(max, in, out_td, FFTW_ESTIMATE);
    cout << "executing backward fft " << max << endl;
    fftw_execute(p);
    cout << "done backward fft" << endl;
  
    //write time domain to file
    Double_t delt = 1/2.0/fMax;//in seconds, =1/SF
    double sig[max];
    double noise[max];
    for (int j=0;j<max;j++)
    {
      out_td[j] = out_td[j]/2/(nEntries-1);//backward needs to be normalized
      anti.i = j;
      anti.t = (j+1)*delt;
      anti.Ef = out_td[j];//no noise
      sig[j] = sqrt(tl_data.Zc/tl_data.Zw)*1e-17*out_td[j];
      noise[j] = r3->Gaus(0, sqrt(kT*fMax*tl_data.Zc)*1e-19);
      out_td[j] = sig[j];
      if (addNoise) out_td[j] = sig[j]+noise[j];//units sqrt fW
      anti.sig = sig[j];
      anti.noise = noise[j];
      anti.vtot = sig[j]+noise[j];
      anti.fW_t = pow(sig[j],2)/tl_data.Zc*1e+38;//units of fW, no noise
      anti.nfW_t = pow(anti.vtot,2)/tl_data.Zc*1e+38;//units of fW, with noise
      tfout->cd();
      antTree->Fill();
    }
    //now back to freq domain for short sample
    cout << "initializing short forward fft " << max << endl;
    fftw_complex *out_s;
    fftw_plan p_forward;
    out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * max);
   //from time to freq. 
    p_forward = fftw_plan_dft_r2c_1d(max, out_td, out_s, FFTW_ESTIMATE);
    cout << "executing short forward fft " << max << endl;
    fftw_execute(p_forward);
    double *energyf, *powerf, *ESD;
    //short time seriew 
    energyf = (double*) fftw_malloc(sizeof(double) *( max/2+1));
    powerf = (double*) fftw_malloc(sizeof(double) * (max/2+1));
    ESD = (double*) fftw_malloc(sizeof(double) * (max/2+1));
    for (int j=0;j<=max/2;j++)
    {
      powerf[j] = 1e+38/tl_data.Zc/max/max*(pow(out_s[j][0],2)+pow(out_s[j][1],2));
      energyf[j] = 1e+38*delt/tl_data.Zc/max*(pow(out_s[j][0],2)+pow(out_s[j][1],2));
      ESD[j] = 1e+38*(pow(delt*out_s[j][0],2)+pow(delt*out_s[j][1],2))/tl_data.Zc;
    }
    //now add in negative frequencies 
    for (int j=max/2+1;j<max;j++)
    {
      powerf[max-j] += 1e+38/tl_data.Zc/max/max*(pow(out_s[max-j][0],2)+pow(out_s[max-j][1],2));
      energyf[max-j] += 1e+38*delt/tl_data.Zc/max*(pow(out_s[max-j][0],2)+pow(out_s[max-j][1],2));
      ESD[max-j] += 1e+38*(pow(delt*out_s[max-j][0],2)+pow(delt*out_s[max-j][1],2))/tl_data.Zc;
    }
    double delf = 1.0/max/delt;
    for (int j=0;j<max/2;j++)
    {
      pc.Hz = j*delf;
      pc.outr = out_s[j][0];
      pc.outi = out_s[j][1];
      pc.fW_f = powerf[j];
      pc.fJ_f = energyf[j];
      pc.fJpHz = ESD[j];
      tfout->cd();
      fftTree->Fill();
    }

    //add extra time points of noise
    int Ntot = max; //max number of time points for this sim
    double tlength = 0.0020; //duration of extended waveform in sec
    int start = 0;//start of signal in time bins
    if (extendWF ) { 
      Ntot = tlength/delt; //max number of time points for this sim
      start = int(r3->Rndm()*Ntot);
      cout << start << endl;
      cout << "Creating " << Ntot << " points, " << Ntot*delt << " s long with data at " << start*delt << endl;
    }
    ofstream outFile;
    if (printDAT) {
      outFile.open(cardname+"_VolVsT.dat");//fix
      outFile << "time [s]\t" << "voltage [sqrt(fW)]"  << endl;
      //outFile.setf(ios::fixed,ios::floatfield);
      outFile.precision(6);
      outFile.setf(ios::scientific,ios::floatfield);
    }
    //will be fft input
    double *wf;
    wf = (double*) fftw_malloc(sizeof(double) * Ntot);
    for (int j=0;j<Ntot;j++)
    {
      wfi.i = j;
      wfi.t = j*delt;
      wf[j]=0;
      if (addNoise) wf[j] = r3->Gaus(0, sqrt(kT*fMax*tl_data.Zc)*1e-19);
      wfi.sig = 0;
      wfi.Ef = 0;//no noise
      wfi.noise = wf[j];
      if (j>=start&&j<start+max) { 
        wf[j] = out_td[j-start];//units volts
        wfi.sig = sig[j-start];//units volts
        wfi.Ef = sig[j-start]*1e17/sqrt(tl_data.Zc/tl_data.Zw);//no noise
        wfi.noise = noise[j-start];
      }
      wfi.vtot = wf[j];
      wfi.fW_t = pow(wfi.sig,2)/tl_data.Zc*1e+38;//units of fW, no noise
      wfi.nfW_t = pow(wf[j],2)/tl_data.Zc*1e+38;//units of fW, includes noise
      if (printDAT) outFile << j*delt << "\t" << wf[j] << endl;
      tfex->cd();
      wfTree->Fill();
    }
    if (printDAT) outFile.close();
    //now back to freq domain
    fftw_complex *out_fd;
    fftw_plan p_flong;
    out_fd = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ntot);
   //from time to freq. 
    p_flong= fftw_plan_dft_r2c_1d(Ntot, wf, out_fd, FFTW_ESTIMATE);
    cout << "executing long forward fft " << Ntot << endl;
    fftw_execute(p_flong);
   
    cout << "done forward fft" << endl;
    //iniitalize other arrays
    double *nenergyf, *npowerf, *nESD;
    //long time seriew 
    nenergyf = (double*) fftw_malloc(sizeof(double) * (Ntot/2+1));
    npowerf = (double*) fftw_malloc(sizeof(double) * (Ntot/2+1));
    nESD = (double*) fftw_malloc(sizeof(double) * (Ntot/2+1));

    //calculate power and energy from long sample fft
    for (int j=0;j<=Ntot/2;j++)
    {
      //mean square amplitude, time-averaged power at bin f is fW_f
      //should have units of fW 
      npowerf[j] = 1e+38/tl_data.Zc/Ntot/Ntot*(pow(out_fd[j][0],2)+pow(out_fd[j][1],2));
      //energyf is time-integral square amplitude, total energy at bin f is fJ_f
      //should have units of fJ
      nenergyf[j] = 1e+38*delt/tl_data.Zc/Ntot*(pow(out_fd[j][0],2)+pow(out_fd[j][1],2));
      nESD[j] = 1e+38*(pow(delt*out_fd[j][0],2)+pow(delt*out_fd[j][1],2))/tl_data.Zc;
    }
    //now add in negative frequencies 
    for (int j=Ntot/2+1;j<Ntot;j++)
    {
      //mean square amplitude, time-averaged power at bin f is fW_f
      //should have units of fW 
      npowerf[Ntot-j] += 1e+38/tl_data.Zc/Ntot/Ntot*(pow(out_fd[Ntot-j][0],2)+pow(out_fd[Ntot-j][1],2));
      //time-integral square amplitude, total energy at bin f is fJ_f
      //should have units of fJ
      nenergyf[Ntot-j] += 1e+38*delt/tl_data.Zc/Ntot*(pow(out_fd[Ntot-j][0],2)+pow(out_fd[Ntot-j][1],2));
      nESD[Ntot-j] += 1e+38*(pow(delt*out_fd[Ntot-j][0],2)+pow(delt*out_fd[Ntot-j][1],2))/tl_data.Zc;
    }
    //now make histogram and fill tree
    delf = 1.0/Ntot/delt;
    for (int j=0;j<Ntot/2;j++)
    {
      npc.Hz = j*delf;
      npc.outr = out_fd[j][0];
      npc.outi = out_fd[j][1];
      npc.fW_f = npowerf[j];
      npc.fJ_f = nenergyf[j];
      npc.fJpHz = nESD[j];
      tfex->cd();
      nfftTree->Fill();
    }
    tfout->cd();
    fftTree->Write();
    antTree->Write();
    tfex->cd();
    ntarray[0] = event;
    ntarray[1] = LO;
    ntarray[2] = SF;
    ntarray[3] = start*delt;
    filter->Fill(ntarray);
    wfTree->Write();
    nfftTree->Write();

    fftw_destroy_plan(p);
    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_flong);
    fftw_free(in);
    fftw_free(out_td);
    fftw_free(out_fd);
    fftw_free(wf);
    fftw_free(out_s);
    fftw_free(ESD);
    fftw_free(powerf);
    fftw_free(energyf);
    fftw_free(nESD);
    fftw_free(npowerf);
    fftw_free(nenergyf);
  }
  filter->Write();
  tfout->Close();
  tfex->Close();
  tfin->Close();
  cout << "Ntuple filled " << endl;
  

}

