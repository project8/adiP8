#include "TFile.h"
#include "TTree.h"
#include "fftw3.h"

class window_fft
{
private:
  //IO (files, trees etc)
  char filename[255];
  TFile *tfile;
  char intreename[255];
  TTree *intree;
  TTree *fwdXform;
  TTree *bwdXform;

  int fwdpasses[2][10000];
  int nfwdpasses;
  int bwdpasses[2][10000];
  int nbwdpasses;

  //FFTW related
  int fft_inited;
  double *inXform;
  fftw_complex *outXform;
  fftw_plan p;
public:
  /**********Creator(s) and Destructor(s)******************/
  window_fft();
  ~window_fft();
  /**********General Use***********************************/
  void open_file();
  void close_file();
  void setup_fft();
  void cleanup_fft();
  void find_mc_passes();
  void dft_a_pass(int, bool);
  void dft_a_pass(int, int);
  void write_pass(int, bool);
  /**********Get & Set*************************************/
  char *get_filename();
  void  set_filename(char*);
  char *get_intreename();
  void  set_intreename(char*);
  int get_nfwdpasses();
  int get_nbwdpasses();
};
