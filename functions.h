#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include "SistemasLineares.h"

/*teste se a alocacao de memoria funcioou*/
void testaMalloc(void *t);

/*aloca os coeficientes de uma matriz*/
real_t **alocaCoef(int n);

/*aloca os termos independentes*/
real_t *alocaVetor(int n);

/*encontra o maior termo de uma coluna*/
int encontraMaior(SistLinear_t *SL, int i);

/*troca a linha de um sistema linear*/
void trocaLinha(SistLinear_t *SL, int i, int maior);

/*calcula o vetor solucao*/
void retroS(SistLinear_t *SL, real_t *x);

/*libera um vetor solução*/
void liberaVetor(real_t *x);

SistLinear_t *copiaMatriz(SistLinear_t *SL);

void pulaLinha(int n);

real_t *residuo(SistLinear_t *SL, real_t *x);

void somaVetor(real_t *a, real_t *b, int n);

void copiaVetor(real_t *a, real_t *b, int n);

real_t maxDiff(real_t *a, real_t *b, int n);

void chamaRefinamento(SistLinear_t *SL, real_t *X, double *tTotal);

#endif