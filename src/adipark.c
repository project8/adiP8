/***********************************************/
/*                  ADIPARK                    */
/*        Adiabatic Particle Tracking          */
/* July to September 2001, by Thomas Thümmler  */
/*   Thanks for Support From Chr. Weinheimer   */
/*             and Björn Flatt                 */
/*                                             */
/*Based on                                     */
/*               Transmission                  */
/*           Adiabatic Ray Tracing             */
/*       May,June 2000, by Lars Fickinger      */
/*                                             */
/*          Development going on...            */
/*                                             */
/***********************************************/
using namespace std;
#define current_version "Version 4.1 Beta 3 (valerius)"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cmath>
#include "paramanage.h"
#include "sim_core.h"
#include "sim_pilot.h"
#include "mag_pa_tool.h"
#include "el_pa_tool.h"
#include "sim_help.h"

struct t_parameter parameter;

int main(int argc, char *argv[])
{

  int choice;
  int init_status, parameter_file_found;
  char *parameter_filename = new char[255];

  cout << endl << endl;
  cout << "********************************************" << endl;
  cout << "**                                        **" << endl;
  cout << "**            A D I P A R K +             **" << endl;
  cout << "**                                        **" << endl;
  cout << "**      ADIabatic PARticle tracKing       **" << endl;
  cout << "**                                        **" << endl;
  cout << "**    written by Th. Thuemmler in 2001    **" << endl;
  cout << "**                                        **" << endl;
  cout << "**        Development Version 2004        **" << endl;
  cout << "**                                        **" << endl;
  cout << "********************************************" << endl;
  cout << "**     " << current_version << "      **" << endl;
  cout << "********************************************" << endl << endl;

  if (argc < 2) {
    cout << "Usage: adipark <parameter filename w/ or w/o suffix>" << endl;
  } else {
    strcpy(parameter_filename, argv[1]);
    init_parameters();
    parameter_file_found = 0;
    if (load_init_data(parameter_filename) == 0) {
      cout << "Try parameter file w/o suffix: ";
      strcpy(parameter_filename, "");
      strncat(parameter_filename, argv[1], strlen(argv[1]) - 4);
      cout << parameter_filename << endl;
      if (load_init_data(parameter_filename) == 0) {
        cout << "Try parameter file w/o point: ";
        strcpy(parameter_filename, "");
        strncat(parameter_filename, argv[1], strlen(argv[1]) - 1);
        cout << parameter_filename << endl;
        if (load_init_data(parameter_filename) == 0) {
          cout << endl << "ERROR: missing some settings in INI-file!";
          cout << endl << endl << flush;
        } else {
          parameter_file_found = 1;
        }
      } else {
        parameter_file_found = 1;
      }
    } else {
      parameter_file_found = 1;
    }

    if (parameter_file_found == 1) {
      strcpy(parameter.filename, parameter_filename);
      init_status = 0;

      if (parameter.enable_epot == 1) {
        alloc_electric_arrays();
        if (read_epot(parameter_filename) != 0) {
          cout << endl << "ERROR: reading of electric potential array failed!" << endl << endl;
          free_electric_arrays();
          init_status -= 10;
        } else {
          init_status += 1;
        }
      } else {
        init_status += 1;
      }

      if (parameter.use_mag_pa == 1) {
        alloc_mag_arrays();
        if (read_mag_pa(parameter_filename) != 0) {
          cout << endl << "ERROR: reading of magnetic potential array failed!" << endl << endl;
          free_mag_arrays();
          init_status -= 10;
        } else {
          init_status += 1;
        }
      } else {
        init_status += 1;
      }

      if (init_status > 0) {
        choice = parameter.adip_run_mode;
        switch (choice) {
          case 1:
            tracking_loop(parameter_filename);
            break;
          case 2:
            cout << "ATTENTION: transmission loop disabled!" << endl;
            // transmission_loop(parameter_filename);
            break;
          case 3:
            trapping_loop(parameter_filename);
            break;
          case 4:
            scan_field_gradient(parameter_filename);
            break;
          case 5:
            scan_Bfield(parameter_filename);
            break;
          case 6:
            mc_tracking(parameter_filename);
            break;
          default:
            cout << "ERROR: Command unknown!" << endl;
        }
        if (ENABLE_EPOT == 1) {
          free_electric_arrays();
        }
        if (USE_MAG_PA == 1) {
          free_mag_arrays();
        }
        cout << "Thank you for using ADIPARK+ (" << current_version << ")." << endl << endl;
      }
    }
  }
  return 0;
}
