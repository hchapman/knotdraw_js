#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "PMdef.h"
#include "PMconjugation.h"
#include "PMenlight.h"
#include "PMextract.h"
#include "PMplanmap.h"
#include "stats.h"
#include "PMemscripten.h"

#define PM_ALL 0
#define PM_PRIME 1
#define PM_CONN6 2
#define PM_BIQUART 3

int randomDiagram(int n_crossings, int n_components, int max_att, int dia_type, int seed, int32_t** vertData) {
  pmMap plmap;
  pmSize size;
  pmMethod meth;
  pmMemory mem;
  pm_edge *cur_e;
  pm_vertex *cur_v;
  int v_idx, e_idx, pos;

  //printf("Hello\n");

  size.e = 0;
  size.v = n_crossings;
  size.f = 0;

  size.minLoopComps = 0;
  size.maxLoopComps = 0;

  size.r = 0;
  size.g = 0;
  size.d = 0;
  size.t = 0;

  size.dgArr = NULL;

  switch (dia_type) {
  case(PM_ALL):
    size.m = PM_MAP_TYPE_QUART_2C;
    size.b = PM_BASIC_TYPE_QUART_2C;
    break;

  case(PM_PRIME):
    size.m = 5;
    size.b = 5;
    break;

  case(PM_CONN6):
    size.m = 6;
    size.b = 5;
    break;

  case(PM_BIQUART):
    size.m = 9;
    size.b = 9;
    break;
  }

  meth.core = 0;
  meth.pic = 0;
  meth.seed = seed;
  meth.verbose = 0;

  //printf("%d \n", meth.seed);

  pmInitRND(&meth);
  pmSetParameters(&size, &meth);
  pmMemoryInit(&size, &meth, &mem);
  pmExtendMemory(&size, &meth, &mem, 0);

  long numTry = 0;
  char done = 0;
  do {
    pmPlanMap(&size, &meth, &mem, &plmap);
    done = (n_components <= 0 ? 1 : pmStatGauss(&plmap) == n_components);
    numTry++;
  } while (!done && numTry < max_att);

  if (numTry >= max_att && !done) {
    pmFreeMap(&plmap);
    return 0;
  }

  *vertData = (int32_t*)malloc(sizeof(int32_t)*plmap.v*4);

  cur_v = plmap.root->from;
  v_idx = cur_v->label-1;
  cur_e = cur_v->root;
  pos = 0;
  while (cur_e != cur_v->root->prev) {
    e_idx = cur_e->label < 0 ? (-2*cur_e->label)-2 : (2*cur_e->label)-1;
    (*vertData)[4*v_idx+pos] = e_idx;
    cur_e = cur_e->next;
    pos += 1;
  }
  e_idx = cur_e->label < 0 ? (-2*cur_e->label)-2 : (2*cur_e->label)-1;
  (*vertData)[4*v_idx+pos] = e_idx;

  while (cur_v->next != NULL) {
    cur_v = cur_v->next;
    v_idx = cur_v->label-1;
    cur_e = cur_v->root;
    pos = 0;
    while (cur_e != cur_v->root->prev) {
      e_idx = cur_e->label < 0 ? (-2*cur_e->label)-2 : (2*cur_e->label)-1;
      (*vertData)[4*v_idx+pos] = e_idx;
      cur_e = cur_e->next;
      pos += 1;
    }
    e_idx = cur_e->label < 0 ? (-2*cur_e->label)-2 : (2*cur_e->label)-1;
    (*vertData)[4*v_idx+pos] = e_idx;
  }

  pmFreeMap(&plmap);

  return plmap.v;
}
