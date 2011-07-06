#include "TFile.h"
#include "TTree.h"

class window_fft
{
private:
  char filename[255];
  char intreename[255];
  int fwdpasses[2][10000];
  int nfwdpasses;
  int bwdpasses[2][10000];
  int nbwdpasses;
  TFile *tfile;
  TTree *intree;
public:
  window_fft();
  ~window_fft();
  /**********General Use***********************************/
  void open_file();
  void close_file();
  void find_mc_passes();
  /**********Get & Set*************************************/
  char *get_filename();
  void  set_filename(char*);
  char *get_intreename();
  void  set_intreename(char*);
};
