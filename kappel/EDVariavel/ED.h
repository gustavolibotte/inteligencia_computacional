#ifndef ED_H_INCLUDED
#define ED_H_INCLUDED

#include <iomanip>

#define DIMENSOES 8          // Função que será otimizada

#define F 0.5;               // Controla a amplificação da variacao (xr2-xr3)
#define CR 0.9;              // Constante de crossover entre 0 e 1
#define NIT 300;             // número de iteracores

using namespace std;

double InfLim[Dimension] = {0.508, 0.0254, 0.0254, 2, 2, 2, 100, 200};
double SupLim[Dimension] = {1.27, 0.254, 0.762, 5, 5, 5, 101.999, 202.999};


// Sorteia um numero entre 0 e 1
double randDouble()
{
    double r = (double)rand()/(double)RAND_MAX;
    return r;
}

// inicializa população
double** inicializaPop(int tamPop, int dimensoes, double** regiaoDeBusca){
    double** matriz = new double*[tamPop];
    for (int i = 0; i < tamPop; i++){
        matriz[i] = new double[dimensoes];
    }
    for (int i = 0; i< tamPop; i++){
        for (int j = 0;  j < dimensoes; j++){
            matriz[i][j] = ((randDouble()*(regiaoDeBusca[j][1]-regiaoDeBusca[j][0])+regiaoDeBusca[j][0]));
        }
    }
    return matriz;
}

// inicializa regiao de busca
double** inicializaRegiao(int dimensoes){
    double** regiao = new double*[dimensoes];
    for (int i = 0; i < dimensoes; i++){
        regiao[i] = new double[2];
    }
    for (int i = 0; i< dimensoes; i++){
        regiao[i][0] = InfLim[i];
        regiao[i][1] = SupLim[i];
    }
    return regiao;
}

// insere regiao de busca
double** insereNaPop(int id, int dimensoes, double* individuo, double** pop){
    for (int i = 0; i < dimensoes; i++){
        pop[id][i] = individuo[i];
    }
    return pop;
}

// pega um individuo da populacao
double* pegaIndividuo(int id, int dimensoes, double** populacao){
    double* individuo = new double[dimensoes];
    for (int i = 0; i < dimensoes; i++){
        individuo[i] = populacao[id][i];
        if (individuo[i] < InfLim[i]){
            individuo[i] = InfLim[i];
        }
        if (individuo[i] > SupLim[i]){
            individuo[i] = SupLim[i];
        }
    }
    return individuo;
}

// soma individuos
double* somaIndividuos(int dimensoes, double* individuo1, double* individuo2){
    double* individuoSoma = new double[dimensoes];
    for (int i = 0; i < dimensoes; i++){
        individuoSoma[i] =  individuo1[i] + individuo2[i];
    }
    return individuoSoma;
}

// subtracao individuos
double* subtracaoIndividuos(int dimensoes, double* individuo1, double* individuo2){
    double* individuoSub = new double[dimensoes];
    for (int i = 0; i < dimensoes; i++){
        individuoSub[i] =  individuo1[i] - individuo2[i];
    }
    return individuoSub;
}

// multiplicacao K x individuo
double* multKindividuo(int dimensoes, double* individuo, double K){
    double* individuoK = new double[dimensoes];
    for (int i = 0; i < dimensoes; i++){
        individuoK[i] =  individuo[i]*K;
    }
    return individuoK;
}

double S1(double x1, double x2){
    return (1+(pow((x1+x2+1),2))*(19-14*x1+3*(pow(x1,2))-14*x2+6*x1*x2+3*(pow(x2,2))));
}
double S2(double x1, double x2){
    return 30+(pow((2*x1-3*x2),2))*(18-32*x1+12*(pow(x1,2))+48*x2-36*x1*x2+27*(pow(x2,2)));
}
//Defindo a Funcao Objetivo a ser usada no metodo.
double funcaoObjetivo(double* x){
//      HAMMER
//      minimoAnalitico = ?;
   return eval(x);
}

double* EvolucaoDiferencial(){
    int dimensoes = DIMENSOES;
    int tamPop = 300;
    double fatorF = F;
    double fatorCR = CR;
    int niter = NIT;

    double** regiaoDeBusca;
    regiaoDeBusca = inicializaRegiao(dimensoes);

    double** populacao;
    populacao = inicializaPop(tamPop, dimensoes, regiaoDeBusca);

    double* melhorX = pegaIndividuo(1, dimensoes, populacao);
    double melhorFit = funcaoObjetivo(melhorX);

    double** popMut = inicializaPop(tamPop, dimensoes, regiaoDeBusca);
    double** popRec = inicializaPop(tamPop, dimensoes, regiaoDeBusca);
    double** novaPop = inicializaPop(tamPop, dimensoes, regiaoDeBusca);
    double* vMut;
    double* v1;
    double* v2;
    double* v3;
    double* v4;
    double* FvAux;
    double* xAnterior;
    double fitU;

    xAnterior = melhorX;

    int i = 1;
    while((i <= niter)){
        printf("\n\n--------------------- Geracao %d --------------------\n\n", i);
        if (melhorX != xAnterior){
            xAnterior = melhorX;
        }

        // Mutacao
        for(int j = 0; j < tamPop; j++){
            int r1, r2, r3;
            r1 = (rand() % tamPop);
            r2 = (rand() % tamPop);
            r3 = j;

            v1 = pegaIndividuo(r1, dimensoes, populacao);
            v2 = pegaIndividuo(r2, dimensoes, populacao);
            v3 = pegaIndividuo(r3, dimensoes, populacao);

            v4 = subtracaoIndividuos(dimensoes, v1, v2);
            FvAux = multKindividuo(dimensoes, v4, fatorF);
            vMut = somaIndividuos(dimensoes, v3, FvAux);

            popMut = insereNaPop(j, dimensoes, vMut, popMut);
        }

        // Recombinacao
        int indiceAleatorio = (rand() % dimensoes);
        int individuoAleatorio = (rand() % tamPop);
        for(int j = 0; j < tamPop; j++){
            double nAleatorioAux = rand();
            v1 = pegaIndividuo(j, dimensoes, popMut);
            if ((nAleatorioAux <= fatorCR) | (j==individuoAleatorio)){
                v2 = pegaIndividuo(j, dimensoes, populacao);
                v1[indiceAleatorio] = v2[indiceAleatorio];
            }
            popRec = insereNaPop(j, dimensoes, v1, popRec);
        }

        // Selecao
        for(int j = 0; j < tamPop; j++){
            fitU = funcaoObjetivo(pegaIndividuo(j, dimensoes, popRec));
            if (fitU <= funcaoObjetivo(pegaIndividuo(j, dimensoes, populacao))){
                printf("\n-------------------------------- Melhorou!\n");
                if (fitU < 1.3){
                    FILE *fp;
                    fp=fopen("regioes.csv","a");
                    fprintf(fp, "%d;%7.10f;%7.10f;%7.10f;%7.10f;%7.10f;%7.10f;%d;%d;%7.10f\n",Num_of_Fitness_Evaluations,vect[0],vect[1],vect[2],vect[3],vect[4],vect[5],(int)floor(vect[6]),(int)floor(vect[7]),f);
                    fclose(fp);
                }
                novaPop = insereNaPop(j, dimensoes, pegaIndividuo(j, dimensoes, popRec), novaPop);
                if (fitU <= melhorFit){
                    melhorFit = fitU;
                    melhorX = pegaIndividuo(j, dimensoes, populacao);
                }
            } else {
                novaPop = insereNaPop(j, dimensoes, pegaIndividuo(j, dimensoes, populacao), novaPop);
            }
        }
        populacao = novaPop;
        i++;
    }

    //liberar memória
    for(int i = 0;i < tamPop; ++i)
        delete [] popMut[i];
    delete [] popMut;
    for(int i = 0;i < tamPop; ++i)
        delete [] popRec[i];
    delete [] popRec;
    for(int i = 0;i < tamPop; ++i)
        delete [] novaPop[i];
    delete [] novaPop;
    delete [] vMut;
    delete [] v1;
    delete [] v2;
    delete [] v3;
    delete [] v4;
    delete [] FvAux;
    delete [] xAnterior;

    return melhorX;
}

double* EvolucaoDiferencialDPS(){
    int dimensoes = DIMENSOES;
    double fatorF = F;
    double fatorCR = CR;
    int niter = NIT;

    int nf = 5;  // tamanho do conjunto F_set
    int ncr = 5; // tamanho do conjunto CR_set
    int nps = 5; // tamanho do conjunto PS_set

    int nComb = nf*ncr; // total de combinações possíveis

    double F_set[] = {0.3, 0.4, 0.5, 0.6, 0.7};
    double CR_set[] = {0.5, 0.6, 0.7, 0.8, 0.9};
    int PS_set[] = {100, 150, 200, 250, 300};
    int tamPop = PS_set[nps-1];

    int **combinacoes = new int*[nComb];
    for (int i=0;i<nComb;i++){
        combinacoes[i] = new int[2];
    }
    int aux = 0;
    for (int m = 0; m < nf; m++){
        for (int n = 0; n < ncr; n++){
            combinacoes[aux][0] = m;
            combinacoes[aux][1] = n;
            aux++;
        }
    }

    int *totalCombUsadas = new int[nComb];
    for (int i = 0; i<nComb; i++){
        totalCombUsadas[i] = 0;
    }
    int *totalCombSucess = new int[nComb];
    for (int i = 0; i<nComb; i++){
        totalCombSucess[i] = 0;
    }
    double *rankCombinacoes = new double[nComb];
    for (int i = 0; i<nComb; i++){
        rankCombinacoes[i] = 0;
    }

    int* comb_individual = new int[PS_set[nps-1]];

    double** regiaoDeBusca;
    regiaoDeBusca = inicializaRegiao(dimensoes);

    double** populacao;
    populacao = inicializaPop(tamPop, dimensoes, regiaoDeBusca);


    // somente para regioes
    double clone[8] = {0.8595113174, 0.1108324252, 0.6520854355, 2.2093628594, 2.3331228566, 4.4609036845, 100, 201};
    populacao = insereNaPop(0, dimensoes, clone , populacao);

    double* melhorX = pegaIndividuo(1, dimensoes, populacao);
    double melhorFit = funcaoObjetivo(melhorX);

    //int period = 0;
    int PS_period = 0;
    int eta = 3; // número de ciclos até reavaliar melhores combinações

    double** popMut = inicializaPop(tamPop, dimensoes, regiaoDeBusca);
    double** popRec = inicializaPop(tamPop, dimensoes, regiaoDeBusca);
    double** novaPop = inicializaPop(tamPop, dimensoes, regiaoDeBusca);
    double* vMut;
    double* v1;
    double* v2;
    double* v3;
    double* v4;
    double* FvAux;
    double* xAnterior;
    double fitU;

    xAnterior = melhorX;
    int melhorCombinacao=0;

    FILE *fp;
    fp=fopen("parametrosUsados.txt","a");
    fprintf(fp, "Geracao;F;CR\n");
    fclose(fp);


    int i = 1;
    while((PS_period <= niter)){
        printf("\n\n--------------------- Geracao %d --------------------\n\n", i);
        if (melhorX != xAnterior){
            xAnterior = melhorX;
        }

        // Sorteando combinações
        if (PS_period%eta==0){
            for(int j=0;j<tamPop;j++){
                 int idComb = rand() % nComb;
                 comb_individual[j] = idComb;
                 totalCombUsadas[idComb] = totalCombUsadas[idComb]+1;
            }
        }

        // Mutacao
        for(int j = 0; j < tamPop; j++){
            int r1, r2, r3;
            r1 = (rand() % tamPop);
            r2 = (rand() % tamPop);
            r3 = j;


            v1 = pegaIndividuo(r1, dimensoes, populacao);
            v2 = pegaIndividuo(r2, dimensoes, populacao);
            v3 = pegaIndividuo(r3, dimensoes, populacao);

            // Local onde F é definido
            if (PS_period%eta!=0){
                fatorF = F_set[combinacoes[melhorCombinacao][0]];
            }else{
                fatorF = F_set[combinacoes[comb_individual[j]][0]];
            }

            v4 = subtracaoIndividuos(dimensoes, v1, v2);
            FvAux = multKindividuo(dimensoes, v4, fatorF);
            vMut = somaIndividuos(dimensoes, v3, FvAux);

            popMut = insereNaPop(j, dimensoes, vMut, popMut);
        }

        // Recombinacao
        int indiceAleatorio = (rand() % dimensoes);
        int individuoAleatorio = (rand() % tamPop);
        for(int j = 0; j < tamPop; j++){
            double nAleatorioAux = rand();


            // Local onde CR é definido
            if (PS_period%eta!=0){
                fatorCR = CR_set[combinacoes[melhorCombinacao][1]];
            }else{
                fatorCR = CR_set[combinacoes[comb_individual[j]][1]];
            }

            v1 = pegaIndividuo(j, dimensoes, popMut);
            if ((nAleatorioAux <= fatorCR) | (j==individuoAleatorio)){
                v2 = pegaIndividuo(j, dimensoes, populacao);
                v1[indiceAleatorio] = v2[indiceAleatorio];
            }
            popRec = insereNaPop(j, dimensoes, v1, popRec);
        }

        // Selecao
        for(int j = 0; j < tamPop; j++){
            fitU = funcaoObjetivo(pegaIndividuo(j, dimensoes, popRec));
            if (fitU <= funcaoObjetivo(pegaIndividuo(j, dimensoes, populacao))){
                printf("\n\n-------------------------------- Melhorou!\n");

                if (fitU < 1.3){
                    FILE *fp;
                    fp=fopen("regioes.csv","a");
                    fprintf(fp, "%d;%7.10f;%7.10f;%7.10f;%7.10f;%7.10f;%7.10f;%d;%d;%7.10f\n",Num_of_Fitness_Evaluations,vect[0],vect[1],vect[2],vect[3],vect[4],vect[5],(int)floor(vect[6]),(int)floor(vect[7]),f);
                    fclose(fp);
                }

                // Salvando sucesso da combinação
                if (PS_period%eta==0){
                    totalCombSucess[comb_individual[j]] = totalCombSucess[comb_individual[j]]+1;
                }

                novaPop = insereNaPop(j, dimensoes, pegaIndividuo(j, dimensoes, popRec), novaPop);
                if (fitU <= melhorFit){
                    melhorFit = fitU;
                    melhorX = pegaIndividuo(j, dimensoes, populacao);
                }
            } else {
                novaPop = insereNaPop(j, dimensoes, pegaIndividuo(j, dimensoes, populacao), novaPop);
            }
        }

        // Calculando ranking de combinações
        if (PS_period%eta==0){
            for (int j = 0; j<nComb; j++){
                if (totalCombUsadas[j]!=0){
                    rankCombinacoes[j] = (double)totalCombSucess[j]/(double)totalCombUsadas[j];
                }
            }
        }

        // Descobringo a melhor combinação
        if (PS_period%eta==0){
            melhorCombinacao = 0;
            for (int j = 0; j<nComb; j++){
                if(rankCombinacoes[melhorCombinacao]<=rankCombinacoes[j]){
                    melhorCombinacao = j;
                }
            }
            FILE *fp;
            fp=fopen("parametrosUsados.txt","a");
            fprintf(fp, "%d;%f;%f\n",i, F_set[combinacoes[melhorCombinacao][0]], CR_set[combinacoes[melhorCombinacao][1]]);
            fclose(fp);
        }
        PS_period++;
        populacao = novaPop;
        i++;
    }

        //liberar memória
    for(int i = 0;i < tamPop; ++i)
        delete [] popMut[i];
    delete [] popMut;
    for(int i = 0;i < tamPop; ++i)
        delete [] popRec[i];
    delete [] popRec;
    for(int i = 0;i < tamPop; ++i)
        delete [] novaPop[i];
    delete [] novaPop;
    delete [] vMut;
    delete [] v1;
    delete [] v2;
    delete [] v3;
    delete [] v4;
    delete [] FvAux;
    delete [] xAnterior;

    return melhorX;
}


#endif // ED_H_INCLUDED
