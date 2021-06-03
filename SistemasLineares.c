#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "functions.h"

/*!
\brief Essa função calcula a norma L2 do resíduo de um sistema linear 

\param SL Ponteiro para o sistema linear
\param x Solução do sistema linear
\param res Valor do resíduo

\return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res)
{
    int i;
    real_t r = 0.0;
    for(i = 0; i < (SL->n); i++){
        r += pow(res[i], 2.0);
    }
    return sqrt(r);

}

/*!
\brief Método da Eliminação de Gauss

\param SL Ponteiro para o sistema linear
\param x ponteiro para o vetor solução
\param tTotal tempo gasto pelo método

\return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss(SistLinear_t *SL, real_t *x, double *tTotal)
{
    //percorre as linhas da matriz
    int i;
    for(i = 0; i < (SL->n); i++){

        //pivoteamento parcial
        int maior = encontraMaior(SL, i);
        if(i != maior)
            trocaLinha(SL, i, maior);

        int k;
        for(k = i+1; k < SL->n; k++){
            
            //calcula M e zera o elemento
            double m = SL->A[k][i] / SL->A[i][i];
            SL->A[k][i] = 0;

            //percorre as colunas da matriz
            int j;
            for(j = i+1; j < SL->n; j++){
                SL->A[k][j] -= SL->A[i][j] * m;
            }

            //termo independente
            SL->b[k] -= SL->b[i] * m;
        }
    }
    retroS(SL, x);

    return 0;
}

/*!
\brief Método de Jacobi

\param SL Ponteiro para o sistema linear
\param x ponteiro para o vetor solução. Ao iniciar função contém
    valor inicial
\param tTotal tempo gasto pelo método

\return código de erro. Um nr positivo indica sucesso e o nr
de iterações realizadas. Um nr. negativo indica um erro:
-1 (não converge) -2 (sem solução)
*/
int gaussJacobi(SistLinear_t *SL, real_t *x, double *tTotal)
{
    //vetor x anterior
    real_t *y = alocaVetor(SL->n);

    int i, j, k;
    /*calcula solucao inicial*/
    for(i = 0; i < SL->n; i++)
        y[i] = 0;

    // printf("Y[0]:\n");
    // prnVetor(y, SL->n);

    real_t diffmax = SL->erro;
    for( k = 0; diffmax >= SL->erro && k < MAXIT; k++){
        
        /*calcula proximo vetor solucao*/
        for(i = 0; i < SL->n; i++){

            double num = 0;
            for(j = 0; j < SL->n; j++){
                if( i != j)
                    num += SL->A[i][j] * y[j];
            }

            x[i] = ( SL->b[i] - num  ) / SL->A[i][i];
        }
        //salva vetor atual no anterior para calcular o proximo
        //e calcula a diferenca
        real_t diff;
        for(i = 0; i < SL->n; i++){
            diff = fabs(x[i] - y[i]);
            if( diffmax < diff)
                diffmax = diff;
            y[i] = x[i];
        }

    }
    liberaVetor(y);

    return k;
}

/*!
\brief Método de Gauss-Seidel

\param SL Ponteiro para o sistema linear
\param x ponteiro para o vetor solução. Ao iniciar função contém
    valor inicial
\param tTotal tempo gasto pelo método

\return código de erro. Um nr positivo indica sucesso e o nr
de iterações realizadas. Um nr. negativo indica um erro:
-1 (não converge) -2 (sem solução)
*/
int gaussSeidel(SistLinear_t *SL, real_t *x, double *tTotal)
{
    //vetor x anterior
    real_t *y = alocaVetor(SL->n);

    int i, j, k = 1;
    /*calcula solucao inicial*/
    for(i = 0; i < SL->n; i++)
        x[i] = 0;

    // printf("Y[0]:\n");
    // prnVetor(y, SL->n);

    real_t diffmax = SL->erro;
    for( k = 0; diffmax >= SL->erro && k < MAXIT; k++){

        for(i = 0; i < SL->n; i++)
            y[i] = x[i];
        
        /*calcula proximo vetor solucao*/
        for(i = 0; i < SL->n; i++){

            double num = 0;
            for(j = 0; j < SL->n; j++){
                if( i != j)
                    num += SL->A[i][j] * x[j];
            }

            x[i] = ( SL->b[i] - num  ) / SL->A[i][i];
        }
        //salva vetor atual no anterior para calcular o proximo
        //e calcula a diferenca
        real_t diff;

        for(i = 0; i < SL->n; i++){
            diff = fabs(x[i] - y[i]);
            if( diffmax < diff)
                diffmax = diff;
            y[i] = x[i];
        }

    }

    liberaVetor(y);
    return k;
}

/*!
\brief Método de Refinamento

\param SL Ponteiro para o sistema linear
\param x ponteiro para o vetor solução. Ao iniciar função contém
    valor inicial para início do refinamento
\param tTotal tempo gasto pelo método

\return código de erro. Um nr positivo indica sucesso e o nr
de iterações realizadas. Um nr. negativo indica um erro:
-1 (não converge) -2 (sem solução)
*/
int refinamento(SistLinear_t *SL, real_t *x, double *tTotal)
{
    real_t *R, *w = alocaVetor(SL->n);
    R = residuo(SL, x);

    int i;
    for(i=0; i < MAXIT ;i++){
    
        // 2 - Calcular o resíduo e testar critério de parada(a)
        if(normaL2Residuo(SL, x, R) < 5.0)
            return i;

        copiaVetor(SL->b, R, SL->n);
        SistLinear_t *SL2;
        SL2 = copiaMatriz(SL);

        // 3 - Obter w resolvendo Aw = r
        eliminacaoGauss(SL2, w, tTotal);

        // 4 - Obter nova solução e testar critério de parada(b)
        somaVetor(x, w, SL->n);
        real_t diff = maxDiff(x , w, SL->n);
        if( diff < SL->erro )
            return i;

        liberaSistLinear(SL2);  
    }
    liberaVetor(w);

    return i;
}

/*!
\brief Alocaçao de memória 

\param n tamanho do SL

\return ponteiro para SL. NULL se houve erro de alocação
*/
SistLinear_t *alocaSistLinear(unsigned int n)
{
    SistLinear_t *SL = malloc ( sizeof(SistLinear_t));
    testaMalloc(SL);

    SL->n = n;

    SL->A = alocaCoef(n);
    SL->b = alocaVetor(n);

    return SL;
}

/*!
\brief Liberaçao de memória 

\param sistema linear SL
*/
void liberaSistLinear(SistLinear_t *SL)
{
    free(SL->b);
    SL->b = NULL;

    free(SL->A[0]);
    SL->A = NULL;

    free(SL->A);
    SL->A = NULL;

    free(SL);
    SL = NULL;
}

/*!
\brief Leitura de SL a partir de Entrada padrão (stdin).

\return sistema linear SL. NULL se houve erro (leitura ou alocação)
*/
SistLinear_t *lerSistLinear()
{
    int n, fim;
    real_t erro;

    fim = scanf("%d%f", &n, &erro);
    if(fim == EOF)
        exit(0);

    SistLinear_t *SL = alocaSistLinear(n);

    SL->erro = erro;

    int i, j;
    for(i = 0;i < n;i++){
        for(j = 0;j < n;j++)
            scanf("%f", SL->A[i]+j);
    }
    
    for(i = 0;i < n;i++)
        scanf("%f", SL->b+i);

    return SL;
}

// Exibe SL na saída padrão
void prnSistLinear(SistLinear_t *SL)
{
    int i, j;
    for(i = 0; i < SL->n; i++){
        for(j = 0; j < SL->n; j++)
            printf("%f ", SL->A[i][j]); 
        printf("\n");   
    }
}

// Exibe um vetor na saída padrão
void prnVetor(real_t *v, unsigned int n)
{
    int i;
    for(i = 0; i < n; i++){
        printf("%f ", v[i]);
    }
    printf("\n");
}
