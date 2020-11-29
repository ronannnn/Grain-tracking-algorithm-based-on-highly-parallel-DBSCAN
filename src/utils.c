/**
 * Author:   Ruonan Chen (ruonanch@usc.edu)
 * Date:     11/28/20
 * Filename: utils.c
 */

#include "stdlib.h"

static char *substring(const char *string, int position, int length) {
  char *pointer = calloc(length + 1, sizeof(char));
  for (int i = 0; i < length; i++) pointer[i] = string[i + position];
  return pointer;
}