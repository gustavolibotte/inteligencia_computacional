#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <ostream>
#include <string>
#include <vector>
#include <Windows.h>

#pragma once;

const double MARGEM = 0.01;
const double SUB_MOD_DELTA = 0.03;
const double FLUXO_ALVO = 8e-5;

void _RodaHammer(std::wstring path);

class HammerInterface {
	void LeSaidaHammer();
	void EscreveArquivoHammer(double Rfuel,double Rclad,double Req,double Enr1,
								  double Enr2,double Enr3,int MatFuel,int MatClad,
								  double Dfuel,double Dclad,double B2ini);
	void IniciaCelula(double vrf,double vdc,double vre,double ve1,
						double ve2,double ve3,double vmf,double vmc);
	void HammerInterface::SubModeracao();
	void Lekef();

	std::wstring path;
	FILE *fpi;
	FILE *fpo;
public:
	int MFuel,MClad;
	double Rf,Rc,Re,E1,E2,E3;
	double f, fluxo, buckling, kef, ktest, dkef, dfluxo, dfpi, dB2, FP;
	double dktest;

	HammerInterface(std::wstring _path) : buckling(0.0) {
		if(!_path.empty() && _path.back() != L'\\') {
			path = _path + L'\\';
		} else {
			path = _path;
		}
	}

	void HammerInterface::eval(std::vector<double> vect);
	void HammerInterface::eval2(double vect[]);
};


double _Densidade(int Material) {
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

void HammerInterface::EscreveArquivoHammer(double Rfuel,double Rclad,double Req,double Enr1,
								  double Enr2,double Enr3,int MatFuel,int MatClad,
								  double Dfuel,double Dclad,double B2ini) {

	_wfopen_s(&fpi, (path + L"HAMMER.DAT").c_str(), L"wt");

	fprintf(fpi,"2       11\n");
	fprintf(fpi,"'HAMMER.OUT'\n");
	fprintf(fpi,"'LITHELIB.BIN'\n");
	fprintf(fpi,"'HELPLIB.BIN'\n");
	fprintf(fpi,"11          2           1  11       11 11\n");
	fprintf(fpi,"   0         2 3 11 CASO EXEMPLO\n");
	fprintf(fpi,"   1        21      REGIAO 1\n");
	fprintf(fpi,"   2  3 3    3 1    21125.E-04\n");
	fprintf(fpi,"   3  1 1  %d 6    %10.8lf%10.4lf%10.6lf 11538.E-1\n",MatFuel,Rfuel,Dfuel,Enr1);
	fprintf(fpi,"   4                92236.E+00                                  1.E-24\n");
	fprintf(fpi,"   5                94239.E+00                                  1.E-24\n");
	fprintf(fpi,"   6                94240.E+00                                  1.E-24\n");
	fprintf(fpi,"   7                94241.E+00                                  1.E-24\n");
	fprintf(fpi,"   8                54135.E+00                                  1.E-24\n");
	fprintf(fpi,"   9                62149.E+00                                  1.E-24\n");
	fprintf(fpi,"  10  2 2  %d      %10.8lf%10.4lf            332.E+00\n",MatClad,Rclad,Dclad);
	fprintf(fpi,"  11  3 3  300 2  1 %10.8lf713151.E-6            327.E+00\n",Req);
	fprintf(fpi,"  12                 5010.E+00                                125.E-07\n");
	fprintf(fpi,"1 13                  304.E+00                               6632.E-07\n");
	fprintf(fpi,"   1        22      REGIAO 2\n");
	fprintf(fpi,"   2  3 3    3 1    12225.E-03\n");
	fprintf(fpi,"1  3  1 1  %d 6    %10.8lf%10.4lf%10.6lf 11538.E-1\n",MatFuel,Rfuel,Dfuel,Enr2);
	fprintf(fpi,"   1        23      REGIAO 3\n");
	fprintf(fpi,"   2  3 3    3 1    21025.E-03\n");
	fprintf(fpi,"1  3  1 1  %d 6    %10.8lf%10.4lf%10.6lf 11538.E-1\n",MatFuel,Rfuel,Dfuel,Enr3);
	fprintf(fpi,"   1         1      CALCULO ESPACIAL\n");
	fprintf(fpi,"   2  2 3    9          1.E+00\n");
	fprintf(fpi,"   3  1     21  01  33859.E-03%10.6lf\n",B2ini);
	fprintf(fpi,"   4  2     22  01  14961.E-03\n");
	fprintf(fpi,"1  5  3     23  01  70866.E-04\n");
	fputc(13,fpi);
	fclose(fpi);
}

void HammerInterface::LeSaidaHammer() {

	char palavra[130],strFluxo4[11+1],strFP[7+1],strkef[14];
	int it;

	_wfopen_s(&fpo, (path + L"HAMMER.OUT").c_str(), L"rt");

	while (!feof(fpo)){
		fgets(palavra,130,fpo);

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
				fgets(palavra,132,fpo);
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
	}
	fclose(fpo);
}


void HammerInterface::IniciaCelula(double vrf,double vdc,double vre,double ve1,
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

void HammerInterface::SubModeracao() {
		double r = 1;
		EscreveArquivoHammer(Rf,Rc,(1.0 + SUB_MOD_DELTA)*Re,E1,E2,E3,MFuel,MClad,_Densidade(MFuel),_Densidade(MClad),buckling);
		_RodaHammer(path);
		Lekef();
		dktest=abs((ktest-kef)/kef);
}

void HammerInterface::Lekef() {

	char palavra[130],strkef[14];
	int it;

	_wfopen_s(&fpo, (path + L"HAMMER.OUT").c_str(), L"rt");
	while (!feof(fpo)){
		fgets(palavra,130,fpo);
		/* Procura identificador do kef */
		if (strstr(palavra, "KEFF=")){
			// Le kef
			for(it=0;it<10;it++){
				strkef[it]=palavra[5+it];
			}
			strkef[10]='\0';
		}
	}
	ktest = atof(strkef);
	fclose(fpo);
}

void HammerInterface::eval(std::vector<double> vect) {
	IniciaCelula(vect[0],vect[1],vect[2],vect[3],vect[4],vect[5],vect[6],vect[7]);
	EscreveArquivoHammer(Rf,Rc,Re,E1,E2,E3,MFuel,MClad,_Densidade(MFuel),
   							_Densidade(MClad),3.73215);
	_RodaHammer(path);
	LeSaidaHammer();
	dfluxo = fabs((FLUXO_ALVO - fluxo)/ FLUXO_ALVO);
	dkef = fabs(1-kef);
	ktest=0;
	SubModeracao();
	f = FP;
	f += dkef > MARGEM ? dkef : 0.0; 
	f += dfluxo > 0.01 ? dfluxo : 0.0;
	f += ktest < kef ? dktest / SUB_MOD_DELTA : 0.0;
}

void HammerInterface::eval2(double vect[]) {
	IniciaCelula(vect[0],vect[1],vect[2],vect[3],vect[4],vect[5],vect[6],vect[7]);
	EscreveArquivoHammer(Rf,Rc,Re,E1,E2,E3,MFuel,MClad,_Densidade(MFuel),
   							_Densidade(MClad),3.73215);
	_RodaHammer(path);
	LeSaidaHammer();
	dfluxo = fabs((FLUXO_ALVO - fluxo)/ FLUXO_ALVO);
	dkef = fabs(1-kef);
	ktest=0;
	SubModeracao();
	f = FP;
	f += dkef > MARGEM ? dkef : 0.0; 
	f += dfluxo > 0.01 ? dfluxo : 0.0;
	f += ktest < kef ? dktest / SUB_MOD_DELTA : 0.0;
}



void _RodaHammer(std::wstring path) {
	STARTUPINFO si;
    PROCESS_INFORMATION pi;

    ZeroMemory( &si, sizeof(si) );
    si.cb = sizeof(si);
    ZeroMemory( &pi, sizeof(pi) );


    // Start the child process. 
    if( !CreateProcess((path + L"HAMMER32.exe").c_str(),   // No module name (use command line)
        NULL,        // Command line
        NULL,           // Process handle not inheritable
        NULL,           // Thread handle not inheritable
        FALSE,          // Set handle inheritance to FALSE
        REALTIME_PRIORITY_CLASS,              // No creation flags
        NULL,           // Use parent's environment block
		path.c_str(),   // DO not Use parent's starting directory 
        &si,            // Pointer to STARTUPINFO structure
        &pi )           // Pointer to PROCESS_INFORMATION structure
    ) 
    {
        printf( "CreateProcess failed (%d).\n", GetLastError() );
        return;
    }

    // Wait until child process exits.
    WaitForSingleObject( pi.hProcess, INFINITE );

    // Close process and thread handles. 
    CloseHandle( pi.hProcess );
    CloseHandle( pi.hThread );
}