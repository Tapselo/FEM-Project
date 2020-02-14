// MES_PROJ.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "pch.h"
#include "Struktury.h"
#include "Struktury.cpp"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#define eps 1e-12

void print2D(double ** a, int n);
double* PdachUpdate(GlobalData gd, double** globalC, double* wektorT, double* globalP, double* Pdach);
void min_max(double* tab, int n);
double* gauss(int n, double ** tab, double *X);
void simulation(GlobalData &gd, double*t0, double**globalH, double** globalC, double *globalP);

int main()
{
	ifstream plik;
	//							WYBÓR TEST CASE
	//=====================================================================
			plik.open("data.txt");					//4x4	time:500 step:50
			//plik.open("data2.txt");				//31x31 time:100 step:1
	//=====================================================================
	double h, w,ro,cw,k,alfa,t_poczatkowe,t_otoczenia,nH, nW, time, step;
	element_uni elementUniwersalny;
	if (plik.good() == true) {
		//wczytywanie danych
		plik >> h;
		plik >> w;
		plik >> nH;
		plik >> nW;
		plik >> ro;
		plik >> cw;
		plik >> k;
		plik >> alfa;
		plik >> time;
		plik >> step;
		plik >> t_poczatkowe;
		plik >> t_otoczenia;
		//wypisuje wczytane dane
		cout << "Global Data:\nh:  " << h << "  w:" << w << "  nH:" << nH << "  nW:" << nW << "  ro:" << ro << "  cw:" << cw<<
									"  k:"<<k<<" alfa:"<<alfa<<" time of simulation:"<<time<<" time step:"<<step<<" t-start:"<<t_poczatkowe<<
									" t-otoczenia:"<<t_otoczenia<< endl << endl;
		GlobalData gd(h, w, nH, nW,ro,cw,k,alfa,time,step, t_poczatkowe, t_otoczenia);
		Grid g(gd.nN, gd.nE);
		//1 Zajecia (nody, siatka) - Generacja siatki
		g.fillGrid(&g, gd);		
		//2 Zajecia (element uniwersalny) - Obliczanie elementu uniwersalnego
		plik >> elementUniwersalny.n_pc;
		elementUniwersalny.makeElement_Uni();
		//ustalenie punktow calkowania i wag
		for (int i = 0; i < elementUniwersalny.n_pc*elementUniwersalny.n_pc; i++) {
			plik >> elementUniwersalny.pc[i][0];
			plik >> elementUniwersalny.wc[i][0];
			plik >> elementUniwersalny.pc[i][1];
			plik >> elementUniwersalny.wc[i][1];
		}
		elementUniwersalny.fillElement_Uni();
		//3 Zajecia (H C P)
		local *lokalTab=new local[gd.nE];
		double*** lokalneH =new double**[gd.nN];
		double*** lokalneC =new double**[gd.nN];
		double*** lokalneH_bc = new double**[gd.nN];
		double** lokalneP = new double*[gd.nE];
		for (int i = 0; i < gd.nE; i++) {
			lokalTab[i] = local(elementUniwersalny.n_pc);
			lokalTab[i].makeLocal(&elementUniwersalny, &gd);
			lokalneH[i] = lokalTab[i].calculateH(g, i);
			lokalneC[i] = lokalTab[i].calculateC();
			//Dodaje do H warunek brzegowy (Hbc)
			lokalneH_bc[i] = lokalTab[i].calculateHbc(g, i);
			for (int j = 0; j < elementUniwersalny.n_pc*elementUniwersalny.n_pc; j++) {
				for (int k = 0; k < elementUniwersalny.n_pc*elementUniwersalny.n_pc; k++) {
					lokalneH[i][j][k] += lokalneH_bc[i][j][k];
				}
			}
			lokalneP[i] = lokalTab[i].CalculateP(g, gd, i);
		}
		double **globalH = new double*[gd.nN];
		double **globalC = new double*[gd.nN];
		double *globalP = new double[gd.nN];
		for (int i = 0; i < gd.nN; i++) {
			globalH[i] = new double[gd.nN];
			globalC[i] = new double[gd.nN];
			globalP[i] = 0.;
			for (int j = 0; j < gd.nN; j++) {
				globalH[i][j] = 0.;
				globalC[i][j] = 0.;
			}
		}
		//Agregacja macierzy lokalnych
		for (int i = 0; i < gd.nE; i++) {
			for (int k = 0; k < 4; k++) {
				globalP[g.elements[i].ID[k]] += lokalneP[i][k];
				for (int z = 0; z < 4; z++) {
					globalH[g.elements[i].ID[k]][g.elements[i].ID[z]] += lokalneH[i][k][z];
					globalC[g.elements[i].ID[k]][g.elements[i].ID[z]] += lokalneC[i][k][z];
				}
			}
		}
		//Wektor temperatur poczatkowych(krok 0)
		double *t0 = new double[gd.nN];
		for (int i = 0; i < gd.nN; i++) {
			t0[i] = g.nodes[i].T;
		}
		// Hdach=H+C/step
		for (int i = 0; i < gd.nN; i++) {
			for (int j = 0; j < gd.nN; j++) {
				globalH[i][j] += (globalC[i][j] / gd.step);
			}
		}
		//Symulacja
		double czas = step;
		while (czas<=gd.time) {
			cout <<fixed<< "Czas: " << czas << "s.:" << endl;
			simulation(gd, t0, globalH, globalC, globalP);
			czas += gd.step;
		}
	}
	else {
		cout << "Input error" << endl;
	}
	plik.close();
	cout << endl;
	system("PAUSE");
	return 0;
}
void print2D(double ** a, int n) {
	for (int i = 0; i < n; i++)
	{
		cout << "\n|";
		for (int j = 0; j < n; j++) {
			cout << fixed;
			cout << a[j][i] << "   ";
		}
		cout << "|";
	}
	cout << endl;
}

double* PdachUpdate(GlobalData gd, double** globalC, double* t0, double* globalP, double* Pdach) {
	for (int i = 0; i < gd.nN; i++) {
		Pdach[i] = 0.;
		for (int j = 0; j < gd.nN; j++) {
			Pdach[i] += (globalC[i][j] / gd.step) * t0[j];
		}
		Pdach[i] += globalP[i];
	}
	return Pdach;
}

//Wyszukiwanie minimalnej i maksymalnej temperatury
void min_max(double* tab,int n) {
	double min = tab[0];
	double max = tab[0];
	for (int i = 0; i < n; i++) {
		if (min >= tab[i]) min = tab[i];
		if (max <= tab[i]) max = tab[i];
	}
	cout << "MIN_T=" << min << " MAX_T=" << max << endl;
}

//eliminacja gaussa+rozwiazanie ukladu
double* gauss(int n, double ** tab, double *X) {
	int i, j, k;
	double m, s;

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(tab[i][i]) < eps) return NULL;
			m = -tab[j][i] / tab[i][i];
			for (k = i + 1; k <= n; k++)
				tab[j][k] += m * tab[i][k];
		}
	}
	for (i = n - 1; i >= 0; i--)
	{
		s = tab[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= tab[i][j] * X[j];
		if (fabs(tab[i][i]) < eps) return NULL;
		X[i] = s / tab[i][i];
	}
	return X;
}

void simulation(GlobalData &gd, double*t0, double**globalH, double** globalC, double *globalP) {
	//Pdach=globalP+(C/dt)*t0
	double* Pdach = new double[gd.nN];
	Pdach = PdachUpdate(gd, globalC, t0, globalP, Pdach);

	//P dodaje jako nowa kolumne w H (powstaje HplusP do eliminacji gaussa)
	double **HplusP = new double*[gd.nN];
	for (int i = 0; i < gd.nN; i++) {
		HplusP[i] = new double[gd.nN + 1];
	}
	for (int i = 0; i < gd.nN; i++) {
		for (int j = 0; j < gd.nN; j++) {
			HplusP[i][j] = globalH[i][j];
		}
	}
	for (int i = 0; i < gd.nN; i++) {
		HplusP[i][gd.nN] = Pdach[i];
	}
	//Eliminacja gaussa + rozwiazanie ukladu rownan
	t0 = gauss(gd.nN, HplusP, t0);
	//Znalezienie min i max temp w wektorze wynikowym z gaussa
	min_max(t0, gd.nN);
}