/* $Id: magsource.c 235 2007-07-27 08:08:27Z s_voec01 $ */
/*
 * Magsource
 *
 * written by Sebastian Voecking
 *
 * Creates the source files from an input file with magfield3. All output
 * file names are prefixed by the input file name.
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "magfield3.h"
#include "paramanage.h"
struct t_parameter parameter;

int main(int argc, char **argv)
{
    const char *inputfile;
    char *basename;
    int len;
    int c;
    
    while ((c = getopt(argc, argv, "n:")) != -1) {
        switch (c) {
            case 'n':
            magfield3_set_nmaxmag(strtol(optarg, NULL, 0));
            break;
        }
    }

    if (argc - optind < 1) {
        fprintf(stderr, "Syntax: magsource [options] <input file>\n");
        return 1;
    }

    inputfile = argv[1];
    len = strlen(inputfile);

    if (len >= 4 && !strcmp(inputfile + len - 4, ".dat")) {
        basename = (char*) malloc(len - 3);
        strncpy(basename, inputfile, len - 4);
        basename[len - 4] = '\0';
    } else
        basename = strdup(inputfile);

    magfield3_set_prefix(basename);

    magfield3_input_coils(inputfile);
    magfield3_test_coils();
    magfield3_source();

    return 0;
}

//////////////////////////////////////////////////
