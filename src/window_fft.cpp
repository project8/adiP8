// Standard
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
// Root
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TLeaf.h"
#include "TObject.h"
// Custom
#include "window_fft.hpp"
#include "structs.h"

using namespace std;

// Methods..............................
// Construct & Destruct
window_fft::window_fft()
{
  fft_inited = 0;
}

window_fft::~window_fft()
{
  if (fft_inited == 1){
    cout << "fft parameters never freed, attempting now" << endl;
    cleanup_fft();
  }
}

/********General Use************************************/
/*******************************************************/
void window_fft::open_file()
{
  tfile = new TFile(filename, "update");
  if (tfile->IsZombie()) {
    cout << "File failed to open" << endl;
    exit(1);
  }
  intree = (TTree *) tfile->Get(intreename);
  fwdXform = new TTree("fwdXform","windowed transform of forward passes");
  bwdXform = new TTree("bwdXform","windowed transform of backward passes");
  int n;
  fwdXform->Branch("i", &n, "i/I");
  bwdXform->Branch("i", &n, "i/I");
  for (int i = 0; i < intree->GetEntries(); i++) {
    n=i;
    fwdXform->Fill();
    bwdXform->Fill();
  }
  fwdXform->Write("", TObject::kOverwrite);
  bwdXform->Write("", TObject::kOverwrite);
}

/*******************************************************/
void window_fft::close_file()
{
  tfile->Close();
}

/*******************************************************/
void window_fft::setup_fft()
{
  int N = intree->GetEntries();
  inXform = (double *) fftw_malloc(sizeof(double) * N);
  outXform = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
  p = fftw_plan_dft_r2c_1d(N, inXform, outXform, FFTW_ESTIMATE);
  fft_inited = 1;
}

/*******************************************************/
void window_fft::cleanup_fft()
{
  fftw_free(inXform);
  fftw_free(outXform);
  fftw_destroy_plan(p);
  fft_inited = 0;
}

/*******************************************************/
void window_fft::find_mc_passes()
{
  INTERINFO fi;
  intree->SetBranchAddress("fi", &fi);
  int wasin = 0;
  int i = 0;
  int fwdindex = 0;
  int bwdindex = 0;
  nfwdpasses = 0;
  nbwdpasses = 0;

  TNtuple *runcard = (TNtuple *) tfile->Get("runcard");
  runcard->GetEntry(0);
  double xend = runcard->FindLeaf("x_ant")->GetValue();

  while (i < intree->GetEntries()) {
    intree->GetEntry(i);
    if (((fi.x < xend) && (fi.x > -xend)) && (wasin == 0)) { //if just entered antenna record
      if (fi.vx > 0) {
        fwdpasses[0][fwdindex] = fi.i;
        nfwdpasses++;
      } else {
        bwdpasses[0][bwdindex] = fi.i;
        nbwdpasses++;
      }
    } else if (((fi.x > xend) || (fi.x < -xend)) && (wasin == 1)) { //if just left antenna record
      if (fi.vx > 0) {
        fwdpasses[1][fwdindex] = fi.i;
        fwdindex++;
      } else {
        bwdpasses[1][bwdindex] = fi.i;
        bwdindex++;
      }
    } else if ((i == intree->GetEntries()-1) && (wasin == 1)) { //if still inside but last entry record as exit
      if (fi.vx > 0) {
        fwdpasses[1][fwdindex] = fi.i;
        fwdindex++;
      } else {
        bwdpasses[1][bwdindex] = fi.i;
        bwdindex++;
      }
    }

    // Get ready for next pass
    i++;
    if ((fi.x > xend) || (fi.x < -xend)) {
      wasin = 0;
    } else {
      wasin = 1;
    }
  }
}

/*******************************************************/
void window_fft::dft_a_pass(int index, bool fwd)
{
  int start, stop;
  if (fwd) {
    start = fwdpasses[0][index];
    stop  = fwdpasses[1][index];
  } else {
    start = bwdpasses[0][index];
    stop  = bwdpasses[1][index];
  }
  bwdXform = new TTree("bwdXform","windowed transform of backward passes");
  int n;
  fwdXform->Branch("i", &n, "i/I");
  bwdXform->Branch("i", &n, "i/I");
  for (int i = 0; i < intree->GetEntries(); i++) {
    n=i;
    fwdXform->Fill();
    bwdXform->Fill();
  }
  fwdXform->Write("", TObject::kOverwrite);
  bwdXform->Write("", TObject::kOverwrite);
}

void window_fft::dft_a_pass(int start, int stop)
{
  INTERINFO fi;
  intree->SetBranchAddress("fi", &fi);
  for (int i = 0; i < intree->GetEntries(); i++) {
    if (i < start) {
      inXform[i] = 0;
    } else if (i > stop) {
      inXform[i] = 0;
    } else {
      intree->GetEntry(i);
      inXform[i] = fi.Ef;
    }
  }
  fftw_execute(p);
}

/*******************************************************/
void window_fft::write_pass(int pass_id, bool fwd)
{ 
  double in, real, cplx;
  if (fwd) {
    TBranch *inbranch = fwdXform->Branch(Form("pass_%d", pass_id), &in,   Form("pass_%d/D", pass_id));
    TBranch *rebranch = fwdXform->Branch(Form("real_%d", pass_id), &real, Form("real_%d/D", pass_id));
    TBranch *imbranch = fwdXform->Branch(Form("imag_%d", pass_id), &cplx, Form("imag_%d/D", pass_id));
    for (int i = 0; i < intree->GetEntries(); i++) {
      in = inXform[i];
      real = outXform[i][0];
      cplx = outXform[i][1];
      inbranch->Fill();
      rebranch->Fill();
      imbranch->Fill();
    }
    fwdXform->Write("", TObject::kOverwrite);
  } else {
    TBranch *inbranch = bwdXform->Branch(Form("pass_%d", pass_id), &in,   Form("pass_%d/D", pass_id));
    TBranch *rebranch = bwdXform->Branch(Form("real_%d", pass_id), &real, Form("real_%d/D", pass_id));
    TBranch *imbranch = bwdXform->Branch(Form("imag_%d", pass_id), &cplx, Form("imag_%d/D", pass_id));
    for (int i = 0; i < intree->GetEntries(); i++) {
      in = inXform[i];
      real = outXform[i][0];
      cplx = outXform[i][1];
      inbranch->Fill();
      rebranch->Fill();
      imbranch->Fill();
    }
    bwdXform->Write("", TObject::kOverwrite);
  }
}

/********Gets & Sets************************************/
/*******************************************************/
char *window_fft::get_filename()
{
  return filename;
};

void window_fft::set_filename(char *name)
{
  strcpy(filename,name);
};

/*******************************************************/
char *window_fft::get_intreename()
{
  return intreename;
};

void window_fft::set_intreename(char *name)
{
  strcpy(intreename,name);
}

/*******************************************************/
int window_fft::get_nfwdpasses()
{
  return nfwdpasses;
}

/*******************************************************/
int window_fft::get_nbwdpasses()
{
  return nbwdpasses;
}
