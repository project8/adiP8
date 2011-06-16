/* $Id: array.h 235 2007-07-27 08:08:27Z s_voec01 $ */

#ifndef ARRAY_H
#define ARRAY_H

#include <stdlib.h>

typedef struct {
  double *data;
  size_t size;
  void* next;
} Entry;

typedef struct {
  Entry* start;
  Entry* last;
  size_t size;
} Array;

Array* array_new();
void array_free(Array* array);
void array_append(Array* array, double x);
double* array_ptr(Array* array);

#endif

