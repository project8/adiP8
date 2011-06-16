/* $Id: array.c 235 2007-07-27 08:08:27Z s_voec01 $ */

#include "array.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define BLOCKSIZE 512

Entry* entry_new()
{
  Entry* entry = (Entry*) malloc(sizeof(Entry));
  entry->data = (double*) malloc(BLOCKSIZE*sizeof(double));
  entry->size = 0;
  entry->next = NULL;

  return entry;
}

Entry* entry_append(Entry* entry, double x)
{
  if(entry->size < BLOCKSIZE) {
    entry->data[entry->size] = x;
    entry->size++;
    return entry;
  }
  else {
    Entry *newentry = entry_new();
    entry->next = newentry;
    entry_append(newentry, x);
    return newentry;
  }
}

Array* array_new()
{
  Array *array = (Array*) malloc(sizeof(Array));
  array->start = NULL;
  array->last = NULL;
  array->size = 0;

  return array;
}

void array_free(Array* array)
{
  Entry *entry = NULL, *next = NULL;

  next = array->start;
  while(next) {
    entry = next;
    next = (Entry*) entry->next;
    free(entry->data);
    free(entry);
  }
  free(array);
}

void array_append(Array* array, double x)
{
  if(!array->last) {
    Entry *entry = entry_new();
    array->start = entry;
    array->last = entry;
  }

  Entry *old = array->last;
  array->last = entry_append(array->last, x);
  if(old != array->last)
    array->size++;
}

double* array_ptr(Array *array)
{
  double* ptr = (double*) malloc((array->size * BLOCKSIZE + array->last->size) 
      * sizeof(double));
  Entry *entry = array->start;
  int i = 0;
  
  while(entry) {
    memcpy(&ptr[i*BLOCKSIZE], entry->data, entry->size*sizeof(double));
    i++;
    entry = (Entry*) entry->next;
  }

  return ptr;
}

