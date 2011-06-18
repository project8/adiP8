#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "fftw3.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TArc.h"
#include "TLine.h"
#include "TText.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TMarker.h"
#include "TText.h"
#include "TMath.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TRandom3.h"
#include "mag_pa_tool.h"
#include "paramanage.h"
#include "radiation.h"
#include "fft_fitter.h"

using namespace std;

struct t_parameter parameter;

//interpolated data at retarded time at the particle and calculated values
typedef struct {
  Long64_t i;
  Double_t t_ant, t_ret, x, y, z, ycen, zcen, rad;
  Double_t vx, vy, vz;
  Double_t Ef, fW_t, phase, dphdt, omega;
} INTERINFO;
int calculate_radiation(INTERINFO ii, double *in, double dir, double d_ant);

int main(int argc, char *argv[])
{
/*
 1- file reads in tracking information, loop over tracked points, 
 2- interpolates btw track times to get points evenly spaced in time at the antenna 
     and calculates parameters at the retarded time of the particle for forward and
     backward waves
 3- adds waves together if reflected from end of antenna 
 4- fourier transform signal at antenna
  values come in cm, us from adipark, but attempting to use SI units here.
*/

  TApplication theApp("App", &argc, argv);
  TROOT rsession("test", "test");
  //data from adipark tracking, no interpolation
  typedef struct {
    Double_t t, xcen, ycen, zcen, rad;    //units us, cm
    Double_t vx, vy, vz;                  //units m/s 
    Double_t ekin, eloss, b, phase, omega;//units eV, T, rad/us
  } TRACKINFO;
  static TRACKINFO ti;

  //interpolated data at retarded time at the particle and calculated values
  static INTERINFO fi;          //forward wave
  static INTERINFO bi;          //backward wave

  //time-domain antenna data similar to real data, with noise
  typedef struct {
    Long64_t i;
    Double_t t, Ef_for, Ef_bk, nEf;//units s,volts 
    Double_t fW_t, nfW_t;          //units fW
    Double_t sig, noise, vtot;     //units V
  } ANTINFO;
  static ANTINFO anti;

  //fourier transform results
  typedef struct {
    Double_t Hz;                //units Hz
    Double_t outr, outi;        //from fftw
    Double_t fW_f, fJ_f, fJpHz; //units fW, fJ
  } POWERCALC;
  static POWERCALC pc;          //without noise
  static POWERCALC npc;         //with noise


  TString cardname = theApp.Argv(1);
  char cdummy[1000];

  //use Adipark parameter-reading routine to read control card
  char *parameter_filename = new char[255];
  strcpy(parameter_filename, cardname.Data());
  init_parameters();
  int parameter_file_found = 0;
  if (load_init_data(parameter_filename) == 0) {

    char *parameter_filename_sub = new char[255];
    cout << "Try parameter file w/o suffix: ";
    strncat(parameter_filename_sub, parameter_filename, strlen(parameter_filename) - 4);
    cout << parameter_filename_sub << endl;
    if (load_init_data(parameter_filename_sub) == 0) {
      cout << "Try parameter file w/o point: ";
      strcpy(parameter_filename_sub, "");
      strncat(parameter_filename_sub, parameter_filename, strlen(parameter_filename) - 1);
      cout << parameter_filename_sub << endl;
      if (load_init_data(parameter_filename_sub) == 0) {
        cout << endl << "ERROR: missing some settings in INI-file!";
        cout << endl << endl << flush;
      } else {
        parameter_file_found = 1;
        cardname.Resize(strlen(parameter_filename) - 1);
      }
    } else {
      parameter_file_found = 1;
      cardname.Resize(strlen(parameter_filename) - 4);
    }
  } else {
    parameter_file_found = 1;
  }

  //open run file
  TString runname = cardname + TString(".run");
  TString ininame = cardname + TString(".ini");

  ifstream runfile;
  runfile.open(runname.Data());
  runfile.getline(cdummy, 1000);       //read and discard header line
  int repeat;
  double xi, yi, zi, ekin, thetai, phii, mass, charge, phasei;


  int i = 1;
  TString rootname = cardname + TString(".root");
  TFile *tfout = new TFile(rootname, "recreate");

  //set root's global rand gen to Mersene Twister
  TRandom3 *r3 = new TRandom3(0);
  gRandom = r3;
  double kT = (parameter.antenna_temp) * K_BOL;    //J, T=35K, amplifier Teff = 25 K

  TNtuple *runcard = new TNtuple("runcard", "run parameters", "i:repeat:xi:yi:zi:ekin:thetai:phii:mass:charge:phasei:nscatters:n_temp:imp:x_ant:t_del:atten:mean:duration:height:ta_pow:emean:eduration:eheight:status");
  float ntarray[25];
  double pars[7];

  //set up for fourier transform
  int N = parameter.fft_max_npts;
  Bool_t transform = kFALSE;
  if (parameter.fft_on == 1) {
    transform = kTRUE;
  }
  double tstep = parameter.fft_resample_tstep;   //in us
  double tmax = parameter.max_tof_in_sec * S2US;
  Bool_t transNoise = kFALSE;
  if (parameter.antenna_temp > 0) {
    transNoise = kTRUE;
  }
  //set up for radiation calculation
  double x_ant = parameter.antenna_pos;
  int dir = copysign(1.0, x_ant);
  Bool_t shift = kFALSE;
  if (parameter.rad_shift == 1) {
    shift = kTRUE;
  }
  double impedance = parameter.impedance;  //normalized load impedance
  double refCo = (impedance - 1) / (impedance + 1);
  if (std::isinf(impedance)) {
    refCo = 1;
  }
  init_data();                  //set tl geometry in radiation.c 
  switch (parameter.rad_calc_mode) {
    case 2:                      //parallel wires
      init_tl_data(kFALSE);
      break;
    case 3:                      //square wg
      init_sq_wg_data(OMEGA0 / Clight / M2CM);
      //Warning!  Wave imp. for TE modes freq dep, not implemented properly
      break;
    case 4:                      //circular wg
      init_circ_wg_data(OMEGA0 / Clight / M2CM);
      //Warning!  Wave imp. for TE modes freq dep, not implemented properly
      break;
    case 5:                      //parallel plates/strips
      init_pp_data();
      break;
    case 6:                      //coaxial cable
      init_coax_data();
      break;
    case 7:                      //offset parallel wires
      init_tl_data(kTRUE);
      break;
  }
  double EF2SIG = sqrt(tl_data.Zc/tl_data.Zw);
  double SIG2EF = sqrt(tl_data.Zw/tl_data.Zc);
  if (parameter.rad_atten == 0) {
    tl_data.att = 0;
  }
  //1- file reads in tracking information, loop over tracked points, 
  while (runfile >> repeat >> xi >> yi >> zi >> ekin >> thetai >> phii >> mass >> charge >> phasei) {
    cout << endl << endl << "interpolating for particle " << i << endl;
    //for each line of runfile we expect a track file 
    TString trackname = cardname + TString(".track");
    trackname += i;
    ifstream trackfile;
    trackfile.open(trackname.Data());

    //initialize Tree, fill with track information (doubles) 
    TTree *trackTree = new TTree(Form("track_%d", i), "adipark track results");
    trackTree->Branch("ti", &ti, "t/D:xcen/D:ycen/D:zcen/D:rad/D:vx/D:vy/D:vz/D:ekin/D:eloss/D:b/D:phase/D:omega/D");
    TTree *intFTree = new TTree(Form("interF_%d", i), "forward interpolated track results");
    intFTree->Branch("fi", &fi, "i/L:t_ant/D:t_ret/D:x/D:y/D:z/D:ycen/D:zcen/D:r/D:vx/D:vy/D:vz/D:Ef/D:fW_t/D:phase/D:dphdt/D:omega/D");
    TTree *intBTree = new TTree(Form("interB_%d", i), "backward interpolated track results");
    intBTree->Branch("bi", &bi, "i/L:t_ant/D:t_ret/D:x/D:y/D:z/D:ycen/D:zcen/D:r/D:vx/D:vy/D:vz/D:Ef/D:fW_t/D:phase/D:dphdt/D:omega/D");
    //initialize Tree, fill with "real data" (doubles) 
    TTree *wfTree = new TTree(Form("ant_%d", i), "antenna results with noise");
    wfTree->Branch("anti", &anti, "i/L:t/D:Ef_for/D:Ef_bk/D:nEf/D:fW_t/D:nfW_t/D:sig/D:noise/D:vtot/D");

    int status = 0;
    //initialize fft and power at antenna
    double *in_for, *in_bk, *in_n;
    fftw_complex *out;
    fftw_complex *out_n;
    fftw_plan p;
    in_for = (double *) fftw_malloc(sizeof(double) * N);
    in_bk = (double *) fftw_malloc(sizeof(double) * N);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    double *energyf, *powerf, *ESD, *nenergyf, *npowerf, *nESD;
    energyf = (double *) fftw_malloc(sizeof(double) * (N / 2 + 1));
    powerf = (double *) fftw_malloc(sizeof(double) * (N / 2 + 1));
    ESD = (double *) fftw_malloc(sizeof(double) * (N / 2 + 1));
    if (transNoise) {
      in_n = (double *) fftw_malloc(sizeof(double) * N);
      out_n = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
      nenergyf = (double *) fftw_malloc(sizeof(double) * (N / 2 + 1));
      npowerf = (double *) fftw_malloc(sizeof(double) * (N / 2 + 1));
      nESD = (double *) fftw_malloc(sizeof(double) * (N / 2 + 1));
    }

    //In this loop we will take steps of length tstep (generated by it loop) and read lines of the track-report file whenever we need to.

    //values of last tracking simulation point
    double last_simPhase;
    double last_simR;
    double last_simX[3];
    double last_simV[3];
    double last_simOm;
    double last_ek;
    double last_simT;
    double ep, ek, epar, eperp, egain, eloss, b, scatt, wq;
    //read in first line of file, has t = 0
    trackfile >> last_simX[0] >> last_simX[1] >> last_simX[2] >> last_simR >> ep >> last_ek >> epar >> eperp >> egain >> eloss >> b >> last_simPhase >> last_simT >> scatt >> wq >> last_simV[0] >> last_simV[1] >> last_simV[2] >> last_simOm;

    //values of central tracking simulation point
    double simPhase;
    double simR;
    double simX[3];
    double simV[3];
    double simOm;
    double simT;
    //read in second line of file
    trackfile >> simX[0] >> simX[1] >> simX[2] >> simR >> ep >> ek >> epar >> eperp >> egain >> eloss >> b >> simPhase >> simT >> scatt >> wq >> simV[0] >> simV[1] >> simV[2] >> simOm;

    //check for time steps smaller that tstep
    double minStep = tstep;
    double minIt = -1;
    if ((simT - last_simT) < minStep) {    //if real time step is smaller than interpolated time step
      minStep = simT - last_simT;
      minIt = 0;
    }
    //values of next tracking simulation point, will be initialized in loop
    double next_simPhase = 0;
    double next_simR = 0;
    double next_simX[3];
    double next_simV[3];
    double next_simOm;
    double next_ep, next_ek, next_epar, next_eperp, next_egain, next_eloss, next_b, next_scatt, next_wq;
    double next_simT;

    //interpolated values 
    Long64_t itf = 0;
    Long64_t itb = 0;
    Double_t ta0 = abs(x_ant / tl_data.vg);    //time to reach antenna in either direction
    //inter_T is time at antenna minus constant offset (ta0), evenly spaced, evaluate power here
    Double_t inter_T = tstep;
    //for_T is retarted time (at time of emission at the particle) for wave in +dir,evaluate pos/vel here
    Double_t for_T = inter_T;
    //for_T = t_ant + dir*x(t_r)/c, appx and num eval x(t_r) for initial value
    if (shift == 1) {
      for_T = inter_T + dir * (last_simX[0] + 1.e-4 * last_simV[0] * (inter_T - last_simT)) / tl_data.vg;
    }
    //bk_T is retarted time (emission tim at particle) for backward or reflected wave
    //if reflected, will add delay of n_del to antenna
    Int_t n_del = 2 * ta0 / tstep;    //number of tsteps in delay, about 80 
    Double_t t_del = 2 * ta0 - n_del * tstep;
    if (refCo == 0) {
      //if no reflection, no extra delay
      n_del = 0;
      t_del = 0;
    }
    Double_t bk_T = inter_T - t_del;
    //bk_T = t_ant -2*ta0 - dir*x(t_r)/c, appx and num eval x(t_r) for initial value
    if (shift == 1) {
      bk_T = inter_T - t_del - dir * (last_simX[0] + 1.e-4 * last_simV[0] * (inter_T - last_simT)) / tl_data.vg;
    }
    double inter_R, drdt;
    double inter_Phase, dphdt;
    double inter_Om, domdt;
    double inter_X[3], dxdt[3];
    double inter_V[3], dvdt[3];
    int iline = 0;
    int nscatters = 0;

    //loop over all lines in the track file starting with line 3
    while (trackfile >> next_simX[0] >> next_simX[1] >> next_simX[2] >> next_simR >> next_ep >> next_ek >> next_epar >> next_eperp >> next_egain >> next_eloss >> next_b >> next_simPhase >> next_simT >> next_scatt >> next_wq >> next_simV[0] >> next_simV[1] >> next_simV[2] >> next_simOm) {
      iline++;
      //check if real time step is smaller than interpolated time step
      if ((next_simT - simT) < minStep) {
        minStep = next_simT - simT;
        minIt = itf;
      }
      //now calculate derivates
      drdt = (next_simR - last_simR) / (next_simT - last_simT);
      dphdt = (next_simPhase - last_simPhase) / (next_simT - last_simT);
      domdt = (next_simOm - last_simOm) / (next_simT - last_simT);     //rad/us^2
      for (int j = 0; j < 3; j++) {
        dxdt[j] = (next_simX[j] - last_simX[j]) / (next_simT - last_simT);
      }
      for (int j = 0; j < 3; j++) {
        dvdt[j] = (next_simV[j] - last_simV[j]) / (next_simT - last_simT);
      }
      /*if (for_T < last_simT) {
        cout << "Problem!! for_T < last_simT" << " for " << itf;
        cout << " for_T= " << for_T << " last_simT " << last_simT << endl;
      }
      if (bk_T < last_simT) {
        cout << "Problem!! bk_T < last_simT" << " for " << itb;
        cout << " bk_T= " << bk_T << " last_simT " << last_simT << endl;
      }*/
      while ((for_T < simT)) {  //interpolate between simulated steps
        inter_T = (itf + 1) * tstep;
        //2- interpolates btw track times to get points evenly spaced in time at the antenna 
        //and calculate the power at the retarded time of the particle
        inter_R = last_simR + drdt * (for_T - last_simT);    //cm
        inter_Phase = last_simPhase + dphdt * (for_T - last_simT);    //radians 
        inter_Om = last_simOm + domdt * (for_T - last_simT);    //rad/usec, 
        //guiding center position in cm
        for (int j = 0; j < 3; j++) {
          inter_X[j] = last_simX[j] + dxdt[j] * (for_T - last_simT);
        }
        //guiding center velocity in m/s
        for (int j = 0; j < 3; j++) {
          inter_V[j] = last_simV[j] + dvdt[j] * (for_T - last_simT);
        }
        if (abs(inter_X[0]) > abs(x_ant)) {
          cout << "Warning! Passed Antenna! Exiting Interpolation Now!" << endl;     //passed antenna
          status = 1;
          break;                //passed antenna
        }
        //save track values at retarded time
        fi.i = itf;
        fi.t_ant = inter_T;
        fi.t_ret = for_T;
        fi.x = inter_X[0];      //in cm
        fi.y = inter_X[1] + inter_R * cos(inter_Phase);
        fi.z = inter_X[2] + inter_R * sin(inter_Phase);     //cm
        fi.ycen = inter_X[1];
        fi.zcen = inter_X[2];
        fi.rad = inter_R;
        fi.vx = inter_V[0] * 1.e-4;     //cm/us
        fi.vy = inter_V[1] * 1.e-4 - inter_R * inter_Om * sin(inter_Phase);
        fi.vz = inter_V[2] * 1.e-4 + inter_R * inter_Om * cos(inter_Phase);
        fi.phase = inter_Phase;
        fi.dphdt = dphdt;       //at inter_T
        fi.omega = inter_Om;    //at for_T
        status = calculate_radiation(fi, in_for, dir, x_ant);
        fi.Ef = in_for[itf];    //coeff of efield, C*Ohm/s
        fi.fW_t = 1.e+15 / tl_data.Zw * (pow(in_for[itf], 2));  //should have units of fW 
        intFTree->Fill();
        itf++;
        inter_T = (itf + 1) * tstep;
        for_T = inter_T;
        if (shift == 1) {
          for_T = inter_T + dir * inter_X[0] / tl_data.vg;
        }
        if (itf >= N) {
          break;
        }
        if (inter_T > tmax) {
          break;
        }
        if (status > 0) {
          break;
        }
        //if (zin[it] > 0 && yin[it] > 0 && it > 1000) break;//for tracking half orbit of offset wires
        //if (inter_X[2] < 0 && inter_X[1] > 0 && it > 1000) break;//for tracking half orbit of symmetric wires
      }
      while ((bk_T < simT)) {   //interpolate between simulated steps
        inter_T = (itb + 1) * tstep;
        //3- interpolates btw track times to get power reflected from end of antenna 
        //calculate interpolated points
        inter_R = last_simR + drdt * (bk_T - last_simT);
        inter_Phase = last_simPhase + dphdt * (bk_T - last_simT); //radians 
        inter_Om = last_simOm + domdt * (bk_T - last_simT); //rad/usec, 
        for (int j = 0; j < 3; j++) inter_X[j] = last_simX[j] + dxdt[j] * (bk_T - last_simT);
        //velocity in m/s
        for (int j = 0; j < 3; j++) {
          inter_V[j] = last_simV[j] + dvdt[j] * (bk_T - last_simT);
        }
        if (abs(inter_X[0]) > abs(x_ant)) {
          cout << "Warning! Passed Antenna! Exiting Interpolation Now!" << endl;      //passed antenna
          status = 1;
          break;                //passed antenna
        }
        //save track values at retarded time
        bi.i = itb;
        bi.t_ant = inter_T;
        bi.t_ret = bk_T;
        bi.x = inter_X[0];      //in cm
        bi.y = inter_X[1] + inter_R * cos(inter_Phase);
        bi.z = inter_X[2] + inter_R * sin(inter_Phase);     //cm
        bi.ycen = inter_X[1];
        bi.zcen = inter_X[2];
        bi.rad = inter_R;
        bi.vx = inter_V[0] * 1.e-4;    //cm/us
        bi.vy = inter_V[1] * 1.e-4 - inter_R * inter_Om * sin(inter_Phase);
        bi.vz = inter_V[2] * 1.e-4 + inter_R * inter_Om * cos(inter_Phase);
        bi.phase = inter_Phase;
        bi.dphdt = dphdt;       //at inter_T
        bi.omega = inter_Om;    //at for_T
        //bi.omega = inter_Om/(1+dir*inter_V[0]*1.e-4/c); 
        if (refCo == 0) {
          status = calculate_radiation(bi, in_bk, -dir, x_ant);
        } else {
          status = calculate_radiation(bi, in_bk, -dir, 3 * x_ant);
        }
        bi.Ef = in_bk[itb];     //coeff of efield,volts 
        bi.fW_t = 1.e+15 / tl_data.Zw * (pow(in_bk[itb], 2));    //should have units of fW 
        intBTree->Fill();
        itb++;
        inter_T = (itb + 1) * tstep;
        bk_T = inter_T - t_del;
        if (shift == 1) {
          bk_T = inter_T - t_del - dir * inter_X[0] / tl_data.vg;
        }
        if (itb >= N) {
          break;
        }
        if (inter_T > tmax) {
          break;
        }
        if (status > 0) {
          break;
        }
        //if (zin[it] > 0 && yin[it] > 0 && it > 1000) break;//for tracking half orbit of offset wires
        //if (inter_X[2] < 0 && inter_X[1] > 0 && it > 1000) break;//for tracking half orbit of symmetric wires
      }
      ti.t = simT;
      ti.xcen = simX[0];
      ti.ycen = simX[1];
      ti.zcen = simX[2];
      ti.rad = simR;
      ti.vx = simV[0];
      ti.vy = simV[1];
      ti.vz = simV[2];
      ti.ekin = ek;
      ti.eloss = eloss;
      ti.b = b;
      ti.phase = simPhase;
      ti.omega = simOm;
      trackTree->Fill();
      if (status > 0) {
        break;
      }
      if (itf >= N) {
        break;
      }
      if (itb >= N) {
        break;
      }
      if (inter_T > tmax) {
        break;
      }
      //if (simX[2] > 0 && simX[1] > 0 && it > 1000) break;//for tracking half orbit of offset wires
      //if (simX[2] < 0 && simX[1] > 0 && it > 1000) break;//for tracking half orbit of symmetric wires
      //reset last values 
      if (last_ek - ek > 1.0) {
        nscatters++;
      }
      last_simT = simT;
      last_simR = simR;
      last_simPhase = simPhase;
      last_simOm = simOm;
      for (int j = 0; j < 3; j++) {
        last_simX[j] = simX[j];
      }
      //for (int j=0; j<3; j++) last_simV[j] = dxdt[j];//calculat from position
      for (int j = 0; j < 3; j++) {
        last_simV[j] = simV[j];
      }
      last_ek = ek;
      //reset central values 
      simT = next_simT;
      simR = next_simR;
      simPhase = next_simPhase;
      simOm = next_simOm;
      for (int j = 0; j < 3; j++) {
        simX[j] = next_simX[j];
      }
      for (int j = 0; j < 3; j++) {
        simV[j] = next_simV[j];
      }
      ek = next_ek;
      b = next_b;
      eloss = next_eloss;
      ep = next_ep;
      epar = next_epar;
      eperp = next_eperp;
      egain = next_egain;
      scatt = next_scatt;
      wq = next_wq;
    }
    /*if (minStep < tstep) {      //if real time step is smaller than interpolated time step
      cout << "Warning!  simulated time step smaller than interpolated time step! ";
      cout << "minStep = " << minStep << " at it = " << minIt << endl;
    }*/
    intFTree->Write();
    intBTree->Write();
    trackTree->Write();
    trackfile.close();
    Long64_t it = TMath::Min(itf, itb);     //max time points
    Long64_t max = TMath::Min(Long64_t(N), it);  //max time points
    Long64_t maxf = floor(max / 2) + 1;        //max freq points
    inter_T = tstep;
    for (it = 0; it < max; it++) {
      //save tree at antenna with noise
      anti.i = it;
      anti.t = inter_T*US2S;//convert to seconds
      anti.Ef_for = in_for[it]; //coeff of efield, volts
      anti.Ef_bk = 0;
      if (n_del <= it) {
        anti.Ef_bk = refCo * in_bk[it - n_del];  //coeff of efield, volts
      }
      //use characteristic impedance to get voltage at antenna
      anti.sig = EF2SIG * in_for[it];     //in Volts
      if ((!(impedance == 1)) && it >= n_del) {
        anti.sig += EF2SIG * refCo * in_bk[it - n_del];     
        in_for[it] += refCo * in_bk[it - n_del];
      }
      anti.noise = r3->Gaus(0, TMath::Sqrt(kT / 2 / tstep / US2S * tl_data.Zc) );    //in Volts 
      anti.vtot = anti.sig + anti.noise;  //volts
      anti.nEf = SIG2EF * anti.vtot;    //Volts
      if (transNoise) {
        in_n[it] = anti.nEf;
      }
      anti.fW_t = pow(anti.sig, 2) / tl_data.Zc * 1.e+15;     //units of fW 
      anti.nfW_t = pow(anti.vtot, 2) / tl_data.Zc * 1.e+15;   //units of fW 
      wfTree->Fill();
      //finished with antenna data
      inter_T = (it + 2) * tstep;
    }
    wfTree->Write();
    //3- fourier transform signal at antenna, in_for[] is electric field amplitude
    if (transform) {
      cout << "done making time series at trackfile line " << iline << " t=" << inter_T << endl;
      p = fftw_plan_dft_r2c_1d(max, in_for, out, FFTW_ESTIMATE);
      cout << "executing fft for particle " << i << " with points " << max << endl;
      fftw_execute(p);
      cout << "done fft" << i << endl;
      /* cout << "max (time points) " << max << " max/2+1 (freq points) " << maxf << endl;
         cout << "zero freq: " << out[0][0] << " " << out[0][1] << " should be real" << endl; 
         cout << "nyquist freq -1: " << out[maxf-2][0] << " " << out[maxf-2][1] << endl; 
         cout << "nyquist freq: " << out[maxf-1][0] << " " << out[maxf-1][1] << " should be real if max is even" << endl;
         cout << "nyquist freq +1: " << out[maxf][0] << " " << out[maxf][1] << " should be zero" << endl; 
       */
      if (transNoise) {
        cout << "executing fft with noise " << i << " with points " << max << endl;
        p = fftw_plan_dft_r2c_1d(max, in_n, out_n, FFTW_ESTIMATE);
        fftw_execute(p);
      }

      TTree *fftTree = new TTree(Form("fft_%d", i), "fft results");
      fftTree->Branch("pc", &pc, "Hz/D:outr/D:outi/D:fW_f/D:fJ_f/D:fJpHz/D");
      TTree *nfftTree = new TTree(Form("nfft_%d", i), "fft results with noise");
      nfftTree->Branch("npc", &npc, "Hz/D:outr/D:outi/D:fW_f/D:fJ_f/D:fJpHz/D");
      Double_t delf = 1.0e6 / ti.t;    //now ti.t is max time, in Hz
      //calculate power for positive frequencies including 0 and nyquist
      for (int j = 0; j < maxf; j++) {
        //energyf is time-integral square amplitude, total energy at bin f is fJ_f
        //should have units of fJ
        energyf[j] = 1.e15 / tl_data.Zw * tstep * US2S / max * (pow(out[j][0], 2) + pow(out[j][1], 2));
        //powerf is mean square amplitude, time-averaged power at bin f is fW_f
        //should have units of fW 
        powerf[j] = 1.e15 / tl_data.Zw / max / max * (pow(out[j][0], 2) + pow(out[j][1], 2));
        //ESD is appr energy spectral density E(f), total energy at f is integral of ESD 
        //should have units of fJ per Hz
        ESD[j] = 1.e15 / tl_data.Zw * (pow(tstep * US2S* out[j][0], 2) + pow(tstep * US2S* out[j][1], 2));
        if (transNoise) {
          nenergyf[j] = 1.e15 / tl_data.Zw * tstep * US2S / max * (pow(out_n[j][0], 2) + pow(out_n[j][1], 2));
          npowerf[j] = 1.e15 / tl_data.Zw / max / max * (pow(out_n[j][0], 2) + pow(out_n[j][1], 2));
          nESD[j] = 1.e15 / tl_data.Zw * (pow(tstep * US2S * out_n[j][0], 2) + pow(tstep * US2S * out_n[j][1], 2));
        }
      }
      //now add in negative frequencies, essentially adding same power in twice b/c r2c FFTW
      for (int j = maxf; j < max; j++) {     //for k=1...maxf-1(excluding DC and nyquist)
        //time-integral square amplitude, total energy at bin f is fJ_f
        //should have units of fJ
        energyf[max - j] += 1.e15 / tl_data.Zw * tstep * US2S / max * (pow(out[max - j][0], 2) + pow(out[max - j][1], 2));
        //mean square amplitude, time-averaged power at bin f is fW_f
        //should have units of fW 
        powerf[max - j] += 1.e15 / tl_data.Zw / max / max * (pow(out[max - j][0], 2) + pow(out[max - j][1], 2));
        //appr energy spectral density E(f), total energy at f is integral of ESD 
        //should have units of fJ per Hz
        ESD[max - j] += 1.e15 / tl_data.Zw * (pow(tstep * US2S * out[max - j][0], 2) + pow(tstep * US2S * out[max - j][1], 2));
        if (transNoise) {
          nenergyf[max - j] += 1.e15 / tl_data.Zw * tstep * US2S / max * (pow(out_n[max - j][0], 2) + pow(out_n[max - j][1], 2));
          npowerf[max - j] += 1.e15 / tl_data.Zw / max / max * (pow(out_n[max - j][0], 2) + pow(out_n[max - j][1], 2));
          nESD[max - j] += 1.e15 / tl_data.Zw * (pow(tstep * US2S * out_n[max - j][0], 2) + pow(tstep * US2S * out_n[max - j][1], 2));
        }
      }

      //now make histogram and fill tree
      double binWidth = delf / 1.0e9;
      int binL = 22 / binWidth;
      int bins = 10 / binWidth;
      TH1F *h1 = new TH1F(Form("hPS_%d", i), "Power Spectrum", bins, (binL + 0.5) * binWidth, (binL + bins + 0.5) * binWidth);
      h1->SetXTitle("Freq [GHz]");
      h1->SetYTitle("Power per Freq bin");
      for (Int_t j = 0; j < maxf; j++) {
        pc.Hz = Double_t(j) * delf;
        pc.outr = out[j][0];
        pc.outi = out[j][1];
        pc.fJpHz = ESD[j];
        pc.fJ_f = energyf[j];
        pc.fW_f = powerf[j];
        h1->AddBinContent(h1->FindBin(Double_t(j) * delf / 1.0e9), powerf[j]);
        fftTree->Fill();
        if (transNoise) {
          npc.Hz = Double_t(j) * delf;
          npc.outr = out_n[j][0];
          npc.outi = out_n[j][1];
          npc.fJpHz = nESD[j];
          npc.fJ_f = nenergyf[j];
          npc.fW_f = npowerf[j];
          nfftTree->Fill();
        }
      }
      double f0 = h1->GetBinCenter(h1->GetMaximumBin())*1.e9;
      //cout << "Trees filled " << endl;
      //now make arrays and do fitting
      double *time_array = (double *) fftw_malloc(sizeof(double) * N);
      double *powert = (double *) fftw_malloc(sizeof(double) * N);
      double *freq_array = (double *) fftw_malloc(sizeof(double) * (N / 2 + 1));
      double *PSD = (double *) fftw_malloc(sizeof(double) * (N / 2 + 1));
      for (int k = 0; k < max; k++) {
        powert[k] = 1.e15 / tl_data.Zw * (pow(in_for[k], 2));    //should have units of fW 
        time_array[k] = (k + 1) * tstep * US2S;    //convert to s
      }
      for (int k = 0; k < maxf; k++) {
        freq_array[k] = k * delf;
        PSD[k] = ESD[k] / max / tstep / US2S;    //should have units of fW per Hz
      }

      cout << "fitting " << endl;
      fit_pow_to_cos(pars, time_array, powert, f0, max, i);
      fit_fft_to_sinc(pars, freq_array, PSD, f0, i);

      fftTree->Write();
      h1->Write();
      if (transNoise) {
        nfftTree->Write();
      }

      fftw_free(time_array);
      fftw_free(powert);
      fftw_free(freq_array);
      fftw_free(PSD);
    }
    ntarray[0] = i;
    ntarray[1] = repeat;
    ntarray[2] = xi;
    ntarray[3] = yi;
    ntarray[4] = zi;
    ntarray[5] = ekin;
    ntarray[6] = thetai;
    ntarray[7] = phii;
    ntarray[8] = mass;
    ntarray[9] = charge;
    ntarray[10] = phasei;
    ntarray[11] = nscatters;
    ntarray[12] = parameter.antenna_temp;
    ntarray[13] = impedance;
    ntarray[14] = x_ant;
    ntarray[15] = n_del * tstep;
    ntarray[16] = tl_data.att;
    if (transform) {
      ntarray[17] = pars[0];
      ntarray[18] = pars[1];
      ntarray[19] = pars[2];
      ntarray[20] = pars[3];
      ntarray[21] = pars[4];
      ntarray[22] = pars[5];
      ntarray[23] = pars[6];
      ntarray[24] = pars[7];
    }
    runcard->Fill(ntarray);
    cout << "done w/ particle " << i << endl;

    fftw_destroy_plan(p);
    fftw_free(in_for);
    fftw_free(in_bk);
    fftw_free(out);
    fftw_free(ESD);
    fftw_free(powerf);
    fftw_free(energyf);
    if (transNoise) {
      fftw_free(in_n);
      fftw_free(out_n);
      fftw_free(nESD);
      fftw_free(npowerf);
      fftw_free(nenergyf);
    }
    i++;
  }

  runcard->Write();
  tfout->Close();

}


int calculate_radiation(INTERINFO ii, double *in, double dir, double d_ant)
{
  /*
     <insert doc string here>
   */
  Double_t position[3], velocity[3], omega, phase, atten;
  position[0] = ii.x;
  position[1] = ii.y;
  position[2] = ii.z;
  velocity[0] = ii.vx;
  velocity[1] = ii.vy;
  velocity[2] = ii.vz;
  omega = ii.omega;
  phase = ii.phase;
  atten = tl_data.att;
  Long64_t it = ii.i;
  int status = 0;
  double efield[3];
  switch (parameter.rad_calc_mode) {
    case 0:
      in[it] = cos(phase) * exp(-atten * abs(d_ant - position[0]));
      break;
    case 1:
      in[it] = cos(phase) * exp(-atten * abs(d_ant - position[0]));
      break;
    case 2:                      //parallel wire transmission line
      status = get_tl_efield(position, efield);
      in[it] = coeff_of_t(efield, velocity, dir) * exp(-atten * abs(d_ant - position[0]));
      break;
    case 3:                      //square waveguide
      status = get_sq_wg_efield(position, efield);
      in[it] = coeff_of_t(efield, velocity, dir) * exp(-atten * abs(d_ant - position[0]));
      break;
    case 4:                      //circular waveguide
      status = get_circ_wg_efield(phase, position, efield);
      in[it] = coeff_of_t(efield, velocity, dir) * exp(-atten * abs(d_ant - position[0]));
      break;
    case 5:                      //parallel plates
      status = get_pp_efield(position, efield);
      in[it] = coeff_of_t(efield, velocity, dir) * exp(-atten * abs(d_ant - position[0]));
      break;
    case 6:                      //coaxial cables
      status = get_coax_efield(position, efield);
      in[it] = coeff_of_t(efield, velocity, dir) * exp(-atten * abs(d_ant - position[0]));
      break;
    case 7:                      //offset parallel wire transmission line
      status = get_tl_efield(position, efield);
      in[it] = coeff_of_t(efield, velocity, dir) * exp(-atten * abs(d_ant - position[0]));
      break;
  }
  return status;
}
