/**************************************************************************/
/*                OTIMIZACOES EM CALCULOS DE REATORES                     */
/**************************************************************************/

//#include <process.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <ostream>
#include <string>
#include <vector>

#pragma once

/* Function dimension*/
#define Dimension 8

long int Num_of_Fitness_Evaluations;
int ix;

double Best_Fitness;
double Best_Fit_Coord[Dimension];
long int Best_Fit_Num_of_Fitness_Evaluations;

/* Lower and upper bounds of the function*/
double InfLim[Dimension] = {0.2, 0.01, 0.01, 2, 2, 2, 100, 200};
double SupLim[Dimension] = {0.5, 0.1, 0.3, 5, 5, 5, 101, 202.999};

double vect[Dimension];


int 		MFuel,MClad;
double 	Rf,Rc,Re,E1,E2,E3;
double 	f, fluxo, buckling, kef, ktest, dkef, dfluxo, dFP, dB2, FP,
			fluxo_alvo = 8e-5;
double dktest;

double Densidade(int Material) {
	switch (Material){
		case 100 : return(18.71);
		case 101 : return(10.09);
		case 200 : return(2.71);
		case 201 : return(8.02);
		case 202 : return(6.50);
		/*case 100 : return(19.05);//18.71 U
		case 101 : return(10.97);//10.09 UO2
		case 200 : return(2.712);//aluminium
		case 201 : return(7.999);//ss 304
		case 202 : return(6.51);//zircaloy*/
	}
	return -1;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void EscreveArquivoHammer(double Rfuel,double Rclad,double Req,double Enr1,
								  double Enr2,double Enr3,int MatFuel,int MatClad,
								  double Dfuel,double Dclad,double B2ini) {
	FILE *fp;

	fopen_s(&fp,"HAMMER.DAT","wt");
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
void RodaHammer() {
//	gotoxy(1,19);
	system("HAMMER32");
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void LeSaidaHammer() {

	FILE *fp;
	char palavra[130],strFluxo4[11+1],strFP[7+1],strkef[14];
	int it;

	fopen_s(&fp, "HAMMER.OUT","rt");
	while (!feof(fp)){
		fgets(palavra,130,fp);

		/* Procura identificador do kef */
		if (strstr(palavra, "KEFF=")){
		// Le kef
			for(it=0;it<10;it++){
				strkef[it]=palavra[5+it];
			}
			strkef[10]='\0';
			kef=atof(strkef);
		}

		/* Procura identificador do Fluxo 4 e Fator de Pico*/
		if (strstr(palavra,"0REGION    FLUX 1      FLUX 2      FLUX 3      FLUX 4    REG. POWER  SPEC. POWER  PEAK/AVE POWER")){
//			printf(palavra);
			while (!strstr(palavra,"REACTOR")){
				fgets(palavra,132,fp);
//				printf(palavra);
			}
			// Le o fluxo 4
			for(it=0;it<11;it++){
				strFluxo4[it]=palavra[44+it];
			}
			strFluxo4[11]='\0';
			fluxo=atof(strFluxo4);
			// Le o fator de pico
			for(it=0;it<7;it++){
				strFP[it]=palavra[84+it];
			}
			strFP[7]='\0';
			FP=atof(strFP);
			if(FP<=1.0){ //...entao houve erro de leitura!
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
double Lekef() {
	FILE *fp2;
	char palavra[130],strkef[14];
	int it;

	fopen_s(&fp2, "HAMMER.OUT","rt");
	while (!feof(fp2)){
		fgets(palavra,130,fp2);
		/* Procura identificador do kef */
		if (strstr(palavra, "KEFF=")){
			// Le kef
			for(it=0;it<10;it++){
				strkef[it]=palavra[5+it];
			}
			strkef[10]='\0';
		}
	}
	fclose(fp2);
	return(atof(strkef));
}

double SubModeracao(double ff) {
		double r = 1;
		EscreveArquivoHammer(Rf,Rc,1.03*Re,E1,E2,E3,MFuel,MClad,Densidade(MFuel),Densidade(MClad),buckling);
		RodaHammer();
		ktest=Lekef();
		dktest=abs((ktest-kef)/kef);
//		gotoxy(1,14);
		if(ktest<kef){ // Se nao for submoderado
			ff=ff-r*(dktest/0.03); // penaliza a fitness com o a variacao de dk com dRe
			// obs: negativamente pois a fitness ja e negativa e o problema e
			// de maximizacao.
		}
		return ff;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
double Fitness(double Ckef,double Cfluxo,double Cfp,double CB2) {

	double fit;
	double rfp=1; //multiplicadores da penalizacao

	fit=FP;
	if(dkef>Ckef) fit += rfp*dkef;
	if(dfluxo>Cfluxo) fit += rfp*dfluxo;
 	//printf("Fit=%f\nFP=%f\ndkeff=%f\ndfluxo=%f",fit,FP,dkef,dfluxo);
	fit = -fit;
	fit=SubModeracao(fit);

	return fit;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void IniciaCelula(double vrf,double vdc,double vre,double ve1,
						double ve2,double ve3,double vmf,double vmc) {
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
double eval(std::vector<double> & vect){
	IniciaCelula(vect[0],vect[1],vect[2],vect[3],vect[4],vect[5],vect[6],vect[7]);
	EscreveArquivoHammer(Rf,Rc,Re,E1,E2,E3,MFuel,MClad,Densidade(MFuel),
   							Densidade(MClad),3.73215);
	RodaHammer();
	LeSaidaHammer();
	dfluxo = fabs((fluxo_alvo-fluxo)/fluxo_alvo);
	dkef = fabs(1-kef);
	//printf("\nFP=%f\nkeff=%f\nfluxo=%10.8f",FP,kef,fluxo);
	ktest=0; //se for testada a submoderacao, o valor de ktest sera !=0
	f=Fitness(0.01,0.01,0.01,0.01);
	f = -1.0*f;
//	ExibeSaidaHammer();

	Num_of_Fitness_Evaluations++;

	if(Num_of_Fitness_Evaluations%40 == 0){

		printf("\n\n\n\n\n# Eval = %d Best = %f\n\n\n\n\n", Num_of_Fitness_Evaluations, Best_Fitness);

	}

	if(f<Best_Fitness){
		
		Best_Fitness = f;

		for(ix=0; ix<Dimension; ix++) Best_Fit_Coord[ix] = vect[ix];

		Best_Fit_Num_of_Fitness_Evaluations = Num_of_Fitness_Evaluations;

	}

	return(f);
}

double eval2(double vect[]){
	IniciaCelula(vect[0],vect[1],vect[2],vect[3],vect[4],vect[5],vect[6],vect[7]);
	EscreveArquivoHammer(Rf,Rc,Re,E1,E2,E3,MFuel,MClad,Densidade(MFuel),
   							Densidade(MClad),3.73215);
	RodaHammer();
	LeSaidaHammer();
	dfluxo = fabs((fluxo_alvo-fluxo)/fluxo_alvo);
	dkef = fabs(1-kef);
	//printf("\nFP=%f\nkeff=%f\nfluxo=%10.8f",FP,kef,fluxo);
	ktest=0; //se for testada a submoderacao, o valor de ktest sera !=0
	f=Fitness(0.01,0.01,0.01,0.01);
	f = -1.0*f;
//	ExibeSaidaHammer();

	Num_of_Fitness_Evaluations++;

	if(Num_of_Fitness_Evaluations%40 == 0){

		printf("\n\n\n\n\n# Eval = %d Best = %f\n\n\n\n\n", Num_of_Fitness_Evaluations, Best_Fitness);

	}

	if(f<Best_Fitness){
		
		Best_Fitness = f;

		for(ix=0; ix<Dimension; ix++) Best_Fit_Coord[ix] = vect[ix];

		Best_Fit_Num_of_Fitness_Evaluations = Num_of_Fitness_Evaluations;

	}

	return(f);
}
