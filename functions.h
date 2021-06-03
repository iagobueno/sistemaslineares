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

/*da printf \n, n vezes*/
void pulaLinha(int n);

/*calcula o residuo*/
real_t *residuo(SistLinear_t *SL, real_t *x);

/*soma a + b em a*/
void somaVetor(real_t *a, real_t *b, int n);

/*copia a em b*/
void copiaVetor(real_t *a, real_t *b, int n);

/*calcula a diferenca maxima entre os elementos de um vetor a e b*/
real_t maxDiff(real_t *a, real_t *b, int n);

/*chama a funcao refinamento na main*/
void chamaRefinamento(SistLinear_t *SL, real_t *X, double *tTotal);

#endif