/**************************************************************************/
/*                OTIMIZACOES EM CALCULOS DE REATORES                     */
/**************************************************************************/

#include "hammer.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "randomico.h"

#define Dimension 8

long int Num_of_Fitness_Evaluations;
int ix;

float Best_Fitness;
float Best_Fit_Coord[Dimension];
long int Best_Fit_Num_of_Fitness_Evaluations;

/* Lower and upper bounds of the function*/
float InfLim[Dimension] = {0.508, 0.025, 0.025, 2, 2, 2, 100, 200};
float SupLim[Dimension] = {1.27, 0.254, 0.762, 5, 5, 5, 101, 202};

float vect[Dimension];


int MFuel,MClad;
float Rf,Rc,Re,E1,E2,E3;
float f, fluxo, buckling, kef, ktest, dkef, dfluxo, dFP, dB2, FP, fluxo_alvo = 8e-5;

/*
*************************************************
   Finds a random between [Sup, Inf]
*************************************************
*/
float Randomico(float Sup, float Inf)
{
    float Temp,r;
    Temp = Sup - Inf;
    do
    {
        r = Inf + genrand_real2()*Temp;
    }
    while((r > Sup) || (r < Inf));
    return(r);
}

/*
*************************************************
   Finds an integer random between [Sup, Inf]
*************************************************
*/

int int_Randomico(int Sup, int Inf)
{
    float Temp;
    float r;
    int int_r;
    Temp = (float)(Sup - Inf);
    do
    {
        r = Inf + genrand_real2()*(Sup - Inf);
    }
    while((r > Sup) || (r < Inf));
    int_r = (int)r;
    if((r - (float)int_r) > 0.5000)
    {
        int_r = int_r + 1;
    }

    if(int_r < Inf)
    {
        int_r = Inf;
    }
    if(int_r > Sup)
    {
        int_r = Sup;
    }

    return(int_r);
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
float Densidade(int Material)
{
    switch (Material)
    {
    case 100 :
        return(18.71);
    case 101 :
        return(10.09);
    case 200 :
        return(2.71);
    case 201 :
        return(8.02);
    case 202 :
        return(6.50);
    }
    return -1;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void EscreveArquivoHammer(float Rfuel,float Rclad,float Req,float Enr1,
                          float Enr2,float Enr3,int MatFuel,int MatClad,
                          float Dfuel,float Dclad,float B2ini)
{
    FILE *fp;

    fp=fopen("HAMMER.DAT","wt");
    fprintf(fp,"2       11\n");
    fprintf(fp,"'HAMMER.OUT'\n");
    fprintf(fp,"'LITHELIB.BIN'\n");
    fprintf(fp,"'HELPLIB.BIN'\n");
    fprintf(fp,"11          2           1  11       11 11\n");
    fprintf(fp,"   0         2 3 11 CASO EXEMPLO\n");
    fprintf(fp,"   1        21      REGIAO 1\n");
    fprintf(fp,"   2  3 3    3 1    21125.E-04\n");
    fprintf(fp,"   3  1 1  %d 6    %10.8lf%10.4lf%10.6lf 11538.E-1\n",MatFuel,Rfuel,Dfuel,Enr1);
    fprintf(fp,"   4                92236.E+00                                  1.E-24\n");
    fprintf(fp,"   5                94239.E+00                                  1.E-24\n");
    fprintf(fp,"   6                94240.E+00                                  1.E-24\n");
    fprintf(fp,"   7                94241.E+00                                  1.E-24\n");
    fprintf(fp,"   8                54135.E+00                                  1.E-24\n");
    fprintf(fp,"   9                62149.E+00                                  1.E-24\n");
    fprintf(fp,"  10  2 2  %d      %10.8lf%10.4lf            332.E+00\n",MatClad,Rclad,Dclad);
    fprintf(fp,"  11  3 3  300 2  1 %10.8lf713151.E-6            327.E+00\n",Req);
    fprintf(fp,"  12                 5010.E+00                                125.E-07\n");
    fprintf(fp,"1 13                  304.E+00                               6632.E-07\n");
    fprintf(fp,"   1        22      REGIAO 2\n");
    fprintf(fp,"   2  3 3    3 1    12225.E-03\n");
    fprintf(fp,"1  3  1 1  %d 6    %10.8lf%10.4lf%10.6lf 11538.E-1\n",MatFuel,Rfuel,Dfuel,Enr2);
    fprintf(fp,"   1        23      REGIAO 3\n");
    fprintf(fp,"   2  3 3    3 1    21025.E-03\n");
    fprintf(fp,"1  3  1 1  %d 6    %10.8lf%10.4lf%10.6lf 11538.E-1\n",MatFuel,Rfuel,Dfuel,Enr3);
    fprintf(fp,"   1         1      CALCULO ESPACIAL\n");
    fprintf(fp,"   2  2 3    9          1.E+00\n");
    fprintf(fp,"   3  1     21  01  33859.E-03%10.6lf\n",B2ini);
    fprintf(fp,"   4  2     22  01  14961.E-03\n");
    fprintf(fp,"1  5  3     23  01  70866.E-04\n");
    fputc(13,fp);
    fclose(fp);
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void RodaHammer()
{
    system("HAMMER32");
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void LeSaidaHammer()
{

    FILE *fp;
    char palavra[130],strFluxo4[11],strFP[7],strkef[14];
    int it;
    fp=fopen("HAMMER.OUT","rt");
    while (!feof(fp))
    {
        fgets(palavra,130,fp);
        /* Procura identificador do kef */
        if (strstr(palavra, "KEFF="))
        {
            // Le kef
            for(it=0; it<10; it++)
            {
                strkef[it]=palavra[5+it];
            }
            strkef[10]='\0';
            kef=atof(strkef);
        }
        /* Procura identificador do Fluxo 4 e Fator de Pico*/
        if (strstr(palavra,"0REGION    FLUX 1      FLUX 2      FLUX 3      FLUX 4    REG. POWER  SPEC. POWER  PEAK/AVE POWER"))
        {
//			printf(palavra);
            while (!strstr(palavra,"REACTOR"))
            {
                fgets(palavra,132,fp);
//				printf(palavra);
            }
            // Le o fluxo 4
            for(it=0; it<11; it++)
            {
                strFluxo4[it]=palavra[44+it];
            }
            strFluxo4[11]='\0';
            fluxo=atof(strFluxo4);
            // Le o fator de pico
            for(it=0; it<7; it++)
            {
                strFP[it]=palavra[84+it];
            }
            strFP[7]='\0';
            FP=atof(strFP);
            if(FP<=1.0)  //...entao houve erro de leitura!
            {
                FP=10.0;
            }
        }
//      exit(0);
    }
    fclose(fp);
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
float Lekef()
{
    FILE *fp2;
    char palavra[130],strkef[14];
    int it;
    fp2=fopen("HAMMER.OUT","rt");
    while (!feof(fp2))
    {
        fgets(palavra,130,fp2);
        /* Procura identificador do kef */
        if (strstr(palavra, "KEFF="))
        {
            // Le kef
            for(it=0; it<10; it++)
            {
                strkef[it]=palavra[5+it];
            }
            strkef[10]='\0';
        }
    }
    fclose(fp2);
    return(atof(strkef));
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
/*void ExibeSaidaHammer() {
	clrscr();
	gotoxy(1,1);
	printf("                                          \n");
	printf(" G A H A:  Otimizacao de Projeto de Reatores\n");
	printf(" ===========================================\n");
	printf("                                          \n");
	printf(" Parametros da Celula sendo testados:\n");
	printf(" Rf=%6.4lf Rc=%6.4lf Re=%6.4lf E1=%6.4lf E2=%6.4lf E3=%6.4lf MF=%d MC=%d\n",Rf,Rc,Re,E1,E2,E3,MFuel,MClad);
	printf("                                          \n");
	printf(" *** Saida do HAMMER ***\n");
	printf(" f = %lf\n",f);
	printf(" Fluxo = %le\n",fluxo);
	printf(" FP = %le\n",FP);
	printf(" kef = %lf  ktest = %lf \n",kef,ktest);
	printf(" Rf = %5.3lf Rc = %5.3lf Re = %5.3lf   \n",Rf,Rc,Re);
	printf(" E1 = %5.3lf E2 = %5.3lf E3 = %5.3lf   \n",E1,E2,E3);
	printf(" MF = %d  MC = %d           \n",MFuel,MClad);
}*/

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
float SubModeracao(float ff)
{
    float dktest;
    float r = 1;

    EscreveArquivoHammer(Rf,Rc,1.03*Re,E1,E2,E3,MFuel,MClad,Densidade(MFuel),Densidade(MClad),buckling);
    RodaHammer();
    ktest=Lekef();
    dktest=abs((ktest-kef)/kef);

    if(ktest<kef)
    {
        ff=ff-r*(dktest/0.03);
    }
    return ff;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
float Fitness(float Ckef,float Cfluxo,float Cfp,float CB2)
{

    float fit;
    float rfp=1;

    fit=FP;
    if(dkef>Ckef) fit += rfp*dkef;
    if(dfluxo>Cfluxo) fit += rfp*dfluxo;
    fit = -fit;
    fit=SubModeracao(fit);

    return fit;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void IniciaCelula(float vrf,float vdc,float vre,float ve1,
                  float ve2,float ve3,float vmf,float vmc)
{
    Rf=vrf;
    Rc=vrf+vdc;
    Re=Rc+vre;
    E1=ve1;
    E2=ve2;
    E3=ve3;
    MFuel=(int)floor(vmf);
    MClad=(int)floor(vmc);
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
float eval(float vect[])
{
    IniciaCelula(vect[0],vect[1],vect[2],vect[3],vect[4],vect[5],vect[6],vect[7]);
    EscreveArquivoHammer(Rf,Rc,Re,E1,E2,E3,MFuel,MClad,Densidade(MFuel),
                         Densidade(MClad),3.73215);
    RodaHammer();
    LeSaidaHammer();
    dfluxo = fabs((fluxo_alvo-fluxo)/fluxo_alvo);
    dkef = fabs(1-kef);
    ktest=0;
    f=Fitness(0.01,0.01,0.01,0.01);
    f = -1.0*f;
    Num_of_Fitness_Evaluations++;

    if(Num_of_Fitness_Evaluations%40 == 0)
    {
        printf("\n\n\n\n\n# Eval = %d Best = %f\n\n\n\n\n", Num_of_Fitness_Evaluations, Best_Fitness);
    }

    if(f<Best_Fitness)
    {
        Best_Fitness = f;
        for(ix=0; ix<Dimension; ix++)
        {
            Best_Fit_Coord[ix] = vect[ix];
        }
        Best_Fit_Num_of_Fitness_Evaluations = Num_of_Fitness_Evaluations;
    }

    return(f);
}
