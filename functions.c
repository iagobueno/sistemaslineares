#include "SistemasLineares.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>

//testa se a alocacao dinamica deu certo
void testaMalloc(void *t){
    if(t == NULL){
        perror("alocacao de memoria invalida");
        exit(2);
    }
}

//aloca os coeficientes de uma matriz com a
//estratégia de vetor de ponteiros de linhas contíguas
real_t **alocaCoef(int n){
    real_t **c;    
    c = malloc( n * sizeof(real_t*));
    testaMalloc(c);
    c[0] = malloc( n * n * sizeof(real_t));
    testaMalloc(c[0]);

    int i, j;
    for(i = 0;i < n;i++)
        c[i] = c[0] + n * i;

    return c;
}

//aloca um vetor
real_t *alocaVetor(int n){
    real_t *vetor;
    
    vetor = malloc( n * sizeof(real_t));
    testaMalloc(vetor);

    return vetor;
}

//encontra o maior termo de uma coluna
int encontraMaior(SistLinear_t *SL, int i){
    int j, max = i;
    for(j = i+1; j < SL->n; j++){
               if(SL->A[j][i] > SL->A[max][i])
                    max = j;
    }
    return max;
}

//troca a linha de um sistema linear
void trocaLinha(SistLinear_t *SL, int i, int maior){
    int j = 0;
    real_t aux;

    //troca linhas da matriz
    for(j = 0; j < SL->n; j++){
        aux = SL->A[i][j];
        SL->A[i][j] = SL->A[maior][j];
        SL->A[maior][j] = aux;
    }

    //troca termos independentes
    aux = SL->b[i];
    SL->b[i] = SL->b[maior];
    SL->b[maior] = aux;
}

void retroS(SistLinear_t *SL, real_t *x){
    int i;
    for( i = SL->n - 1; i>= 0; i--){

        x[i] = SL->b[i];

        int j;
        for(j = i+1; j < SL->n; j++)
            x[i] -= SL->A[i][j] * x[j];
    
        x[i] /= SL->A[i][i];
    }   
}

void liberaVetor(real_t *x){
    free(x);
}