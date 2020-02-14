#include "Struktury.h"
#include <iostream>
#include <iomanip>
using namespace std;

void Grid::fillGrid(Grid *g, GlobalData gd) {
	double dx = gd.w / (gd.nW -1);
	double dy = gd.h / (gd.nH -1);
	int n = 0; //dodatkowy licznik do zmiany nodow
	for (int i = 0; i < gd.nW; i++) { //tworzenie nodow
		for (int j = 0; j < gd.nH; j++)
		{	
			g->nodes[n].x = i * dx;
			g->nodes[n].y = j * dy;
			g->nodes[n].T = gd.t_start;
			if (i == 0 || j == 0 || i == (gd.nW - 1) || j == (gd.nH - 1)) { //ustalanie czy node lezy na granicy
				g->nodes[n].BC = true;
			}
			else {
				g->nodes[n].BC = false;
			}
			n++;
		}
	}
	int p = 0;
	for (int i = 0; i < gd.nE; i++) //tworzenie siatki elementow z nodow
	{
		if (i % (gd.nH - 1) == 0 && i != 0) p++; //przeskok do kolejnej kolumny lub pierwszy element w ogóle
		g->elements[i].nodes[0] = g->nodes[i+p];
		g->elements[i].nodes[1] = g->nodes[i+p + gd.nH];
		g->elements[i].nodes[2] = g->nodes[i+p + gd.nH + 1];
		g->elements[i].nodes[3] = g->nodes[i+p + 1];
		//uzupelnienie tablicy ID (potrzebne do agregacji)
		g->elements[i].ID[0] = i + p;
		g->elements[i].ID[1] = i + p +gd.nH;
		g->elements[i].ID[2] = i + p+gd.nH+1;
		g->elements[i].ID[3] = i + p +1;
	}
}

void element_uni::makeElement_Uni()
{
	int n = n_pc * n_pc;
	//tworze "1 wymiar"
	pc = new double*[n];
	wc = new double*[n];
	N = new double*[n];
	dN_dKsi = new double*[n];
	dN_dEta = new double*[n];
	Nbc = new double*[2*n];
	pc_bc = new double*[2*n];
	//dokladam "2 wymiar"
	for (int i = 0; i < n; i++) {
		pc[i] = new double[n_pc];
		wc[i] = new double[n_pc];
		N[i] = new double[n];
		dN_dKsi[i] = new double[n];
		dN_dEta[i] = new double[n];
	}
	for (int i = 0; i < 2 * n; i++) {
		Nbc[i] = new double[n];
		pc_bc[i] = new double[n_pc];
	}
}

void element_uni::fillElement_Uni() {
	pc_bc[0][0] = -1.*(1./sqrt(3)); pc_bc[1][0] = (1. / sqrt(3)); pc_bc[2][0] = 1.;					 pc_bc[3][0] = 1.;			   pc_bc[4][0] = (1. / sqrt(3));   pc_bc[5][0] = -1. * (1. / sqrt(3));   pc_bc[6][0] = -1.;			  pc_bc[7][0] = -1.;
	pc_bc[0][1] = -1.;			    pc_bc[1][1] = -1.;            pc_bc[2][1] = -1. * (1. / sqrt(3));pc_bc[3][1] = (1. / sqrt(3)); pc_bc[4][1] = 1.;			   pc_bc[5][1] = 1.;					 pc_bc[6][1] = (1. / sqrt(3)); pc_bc[7][1] = -1. * (1. / sqrt(3));

	for (int i = 0; i < n_pc*n_pc; i++) {
		N[0][i] = 0.25 * (1 - pc[i][0])*(1 - pc[i][1]);
		N[1][i] = 0.25 * (1 + pc[i][0])*(1 - pc[i][1]);
		N[2][i] = 0.25 * (1 + pc[i][0])*(1 + pc[i][1]);
		N[3][i] = 0.25 * (1 - pc[i][0])*(1 + pc[i][1]);

		dN_dKsi[0][i] = -1 * 0.25 * (1 - pc[i][1]);
		dN_dKsi[1][i] = 1 * 0.25 * (1 - pc[i][1]);
		dN_dKsi[2][i] = 1 * 0.25 * (1 + pc[i][1]);
		dN_dKsi[3][i] = -1 * 0.25 * (1 + pc[i][1]);

		dN_dEta[0][i] = -1 * 0.25 * (1 - pc[i][0]);
		dN_dEta[1][i] = -1 * 0.25 * (1 + pc[i][0]);
		dN_dEta[2][i] = 1 * 0.25 * (1 + pc[i][0]);
		dN_dEta[3][i] = 1 * 0.25 * (1 - pc[i][0]);
	}
	for (int i = 0; i < 2 * n_pc*n_pc; i++) {
		Nbc[i][0] = 0.25 * (1 - pc_bc[i][0])*(1 - pc_bc[i][1]);
		Nbc[i][1] = 0.25 * (1 + pc_bc[i][0])*(1 - pc_bc[i][1]);
		Nbc[i][2] = 0.25 * (1 + pc_bc[i][0])*(1 + pc_bc[i][1]);
		Nbc[i][3] = 0.25 * (1 - pc_bc[i][0])*(1 + pc_bc[i][1]);
	}
}

void local::makeLocal(element_uni* eu, GlobalData* gd)
{	
	int n = n_pc*n_pc;
	el_uni = eu;
	data = gd;
	//tablice
	dx_dEta = new double[n];
	dy_dEta = new double[n];
	dx_dKsi = new double[n];
	dy_dKsi = new double[n];
	//1 wymiar
	J = new double*[n];
	detJ = new double[n];
	dN_dx = new double*[n];
	dN_dy = new double*[n];
	H = new double*[n];
	Hpc = new double**[n];
	Nx_NxT = new double**[n];
	Ny_NyT = new double**[n];
	C = new double*[n];
	//2 wymiar
	for (int i = 0; i < n; i++) {
		Hpc[i] = new double*[n];
		Nx_NxT[i] = new double*[n];
		Ny_NyT[i] = new double*[n];
		J[i] = new double[n];
		dN_dx[i] = new double[n];
		dN_dy[i] = new double[n];
		H[i] = new double[n];
		C[i] = new double[n];
		dx_dEta[i] = 0.;
		dy_dEta[i] = 0.;
		dx_dKsi[i] = 0.;
		dy_dKsi[i] = 0.;
	}
	//3 wymiar
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Hpc[i][j] = new double[n];
			Nx_NxT[i][j] = new double[n];
			Ny_NyT[i][j] = new double[n];
			dN_dx[i][j] = 0.;
			dN_dy[i][j] = 0.;
			J[i][j] = 0.;
			H[i][j] = 0.;
			C[i][j] = 0.;
		}
	}
}
		
double** local::calculateH(Grid g, int a)
{
	//a=numer elementu w siatce g
	int n = n_pc * n_pc;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			dx_dKsi[i] += el_uni->dN_dKsi[j][i] * g.elements[a].nodes[j].x;
			dx_dEta[i] += el_uni->dN_dEta[j][i] * g.elements[a].nodes[j].x;
			dy_dKsi[i] += el_uni->dN_dKsi[j][i] * g.elements[a].nodes[j].y;
			dy_dEta[i] += el_uni->dN_dEta[j][i] * g.elements[a].nodes[j].y;
		}
	}

	//det jakobian dla kazdego punktu calkowania (i)
	for (int i = 0; i < n; i++) {
		detJ[i] = dx_dKsi[i] * dy_dEta[i] - dy_dKsi[i] * dx_dEta[i];
		//cout << detJ[i]<<"\t";
	}
	//Liczenie jakobianu przeksztalcenia
	// J tablica pomocnicza do mnozenia calosci (1/det * pierwsza macierz)
	for (int i = 0; i < n; i++) {
		J[i][0] = (1. / detJ[i])*dy_dEta[i];
		J[i][1] = (1. / detJ[i])*(-1)*dy_dKsi[i];
		J[i][2] = (1. / detJ[i])*(-1)*dx_dEta[i];
		J[i][3] = (1. / detJ[i])*dx_dKsi[i];
	}

	//J*druga macierz
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			dN_dx[i][j] += J[i][0] * el_uni->dN_dKsi[i][j] + J[i][1] * el_uni->dN_dEta[i][j];
			dN_dy[i][j] += J[i][2] * el_uni->dN_dKsi[i][j] + J[i][3] * el_uni->dN_dEta[i][j];
		}
	}

	for (int pc = 0; pc < n; pc++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				Nx_NxT[pc][i][j] = dN_dx[i][pc] * dN_dx[j][pc];
				Ny_NyT[pc][i][j] = dN_dy[i][pc] * dN_dy[j][pc];
			}
		}
	}

	//H w kazdym punkcie, czyli H = k * (n_x*n_xT + n_y*n_yT)
	for (int pc = 0; pc < n; pc++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				Hpc[pc][i][j] = data->k * (Nx_NxT[pc][i][j] + Ny_NyT[pc][i][j]);
			}
		}
	}

	//H elementu:
	for (int pc = 0; pc < n; pc++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				H[i][j] += Hpc[pc][i][j] * el_uni->wc[pc][0] * el_uni->wc[pc][1] * detJ[pc];
			}
		}
	}
	return H;
}

double** local::calculateC()
{
	for (int pc = 0; pc < n_pc*n_pc; pc++) {
		for (int i = 0; i < n_pc*n_pc; i++) {
			for (int j = 0; j < n_pc*n_pc; j++) {
				C[i][j] += data->cw * data->ro * el_uni->N[pc][i] * el_uni->N[pc][j] * el_uni->wc[pc][0] * el_uni->wc[pc][1] * detJ[pc];
			}
		}
	}
	//Wyswietlanie macierzy C
	//cout << "Macierz C:" << endl;
	//printMatrix2D(C);
	return C;
}

double** local::calculateHbc(Grid g, int a)
{
	bool PcBc[4];	//ustalam ktorych pc dotyczy bc zeby uproscic zapis
	for (int i = 0; i < 4; i++) { 
		if (g.elements[a].nodes[i].BC == true) {
			PcBc[i] = true;
		}
		else{
			PcBc[i] = false;
		}
	}
	//dlugosci bokow elementu (Tylko jesli elementy sa rowne tak jak u nas wszystkie takim samym kwadratem)
	double dx = g.elements[a].nodes[1].x - g.elements[a].nodes[0].x;
	double dy = g.elements[a].nodes[3].y - g.elements[a].nodes[0].y;
	//tablica do przechowywania warunkow brzegowych
	double**BC = new double*[n_pc*n_pc];
	for (int i = 0; i < n_pc*n_pc; i++) { //tworzenie i uzupelnianie tablic
		BC[i] = new double[n_pc*n_pc];
		for (int j = 0; j < n_pc*n_pc; j++) {
			BC[i][j] = 0.;
		}
	}
	if (PcBc[1] && PcBc[2]) {//prawa sciana
		//cout << "Prawa sciana" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				BC[i][j] += data->alfa * el_uni->Nbc[2][i] * el_uni->Nbc[2][j] * el_uni->wc[1][1] * (dy / 2.);
				BC[i][j] += data->alfa * el_uni->Nbc[3][i] * el_uni->Nbc[3][j] * el_uni->wc[2][1] * (dy / 2.);
			}
		}
	}
	if (PcBc[2] && PcBc[3]) {//gorna sciana
		//cout << "gorna sciana" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				BC[i][j] += data->alfa * el_uni->Nbc[4][i] * el_uni->Nbc[4][j] * el_uni->wc[2][0] * (dx / 2.);
				BC[i][j] += data->alfa * el_uni->Nbc[5][i] * el_uni->Nbc[5][j] * el_uni->wc[3][0] * (dx / 2.);
			}
		}
	}
	if (PcBc[3] && PcBc[0]) {//lewa sciana
		//cout << "lewa sciana" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				BC[i][j] += data->alfa * el_uni->Nbc[6][i] * el_uni->Nbc[6][j] * el_uni->wc[3][1] * (dy / 2.);
				BC[i][j] += data->alfa * el_uni->Nbc[7][i] * el_uni->Nbc[7][j] * el_uni->wc[0][1] * (dy / 2.);
			}
		}
	}
	if (PcBc[0] && PcBc[1]) {//dolna sciana
		//cout << "dolna sciana" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				BC[i][j] += data->alfa * el_uni->Nbc[0][i] * el_uni->Nbc[0][j] * el_uni->wc[0][0] * (dx / 2.);
				BC[i][j] += data->alfa * el_uni->Nbc[1][i] * el_uni->Nbc[1][j] * el_uni->wc[1][0] * (dx / 2.);
			}
		}
	}
	//printMatrix2D(BC);
	return BC;
}

double* local::CalculateP(Grid g, GlobalData gd, int a){
	//P liczymy analogicznie do Hbc, ró¿nica we wzorze
	bool PcBc[4];
	for (int i = 0; i < 4; i++) {
		if (g.elements[a].nodes[i].BC == true) {
			PcBc[i] = true;
		}
		else {
			PcBc[i] = false;
		}
	}
	double dx = g.elements[a].nodes[1].x - g.elements[a].nodes[0].x;
	double dy = g.elements[a].nodes[3].y - g.elements[a].nodes[0].y;
	double* P = new double[n_pc*n_pc];
	for (int i = 0; i < n_pc*n_pc; i++) {
		P[i] = 0.;
	}
	if (PcBc[1] && PcBc[2]) {//prawa sciana
		for (int i = 0; i < 4; i++) {
				P[i] += data->alfa * el_uni->Nbc[2][i] * data->t_otocz * el_uni->wc[1][1] * (dy / 2);
				P[i] += data->alfa * el_uni->Nbc[3][i] * data->t_otocz * el_uni->wc[2][1] * (dy / 2);
		}
	}
	if (PcBc[2] && PcBc[3]) {//gorna sciana
		for (int i = 0; i < 4; i++) {
				P[i] += data->alfa * el_uni->Nbc[4][i] * data->t_otocz * el_uni->wc[2][0] * (dx / 2);
				P[i] += data->alfa * el_uni->Nbc[5][i] * data->t_otocz * el_uni->wc[3][0] * (dx / 2);
		}
	}
	if (PcBc[3] && PcBc[0]) {//lewa sciana
		for (int i = 0; i < 4; i++) {
				P[i] += data->alfa * el_uni->Nbc[6][i] * data->t_otocz * el_uni->wc[3][1] * (dy / 2);
				P[i] += data->alfa * el_uni->Nbc[7][i] * data->t_otocz * el_uni->wc[0][1] * (dy / 2);
		}
	}
	if (PcBc[0] && PcBc[1]) {//dolna sciana
		for (int i = 0; i < 4; i++) {
				P[i] += data->alfa * el_uni->Nbc[0][i] * data->t_otocz * el_uni->wc[0][0] * (dx / 2);
				P[i] += data->alfa * el_uni->Nbc[1][i] * data->t_otocz * el_uni->wc[1][0] * (dx / 2);
		}
	}
	return P;
}

void local::printMatrix2D(double **m)
{
	int n = n_pc * n_pc;
	for (int i = 0; i < n; i++)
	{
		cout << "\n|";
		for (int j = 0; j < n; j++) {
			cout << fixed;
			cout << m[j][i] << "\t";
		}
		cout << "|";
	}
	cout << endl;
}
