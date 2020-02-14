#pragma once
#include <iostream>
using namespace std;


struct GlobalData { //dane wczytane z pliku
public:
	double h = 1, w = 1;
	int nH = 1, nW = 1, nN = 1, nE = 1;
	double k = 25;
	double ro = 1;
	double cw = 1;
	double alfa = 300;
	double time = 100;
	double step = 10;
	double t_start = 100;
	double t_otocz = 1200;
	GlobalData(double a, double b, int c, int d,double e, double f,double g,double x,double t,double s,double aa, double bb) {
		h = a;
		w = b;
		nH = c;
		nW = d;
		nN = nH * nW;
		nE = (nH - 1)*(nW - 1);
		ro = e;
		cw = f;
		k = g;
		alfa = x;
		time = t;
		step = s;
		t_start = aa;
		t_otocz = bb;
	}
};

struct Node {  //punkt na rogach elementu
	double x, y; //wspolrzedne
	double T; //wartosc (np. temperatura)
	bool BC; //czy lezy na granicy obiektu 
};

struct Element { //element z 4 rogami
	Node nodes[4]; //przechowuje 4 node
	int ID[4];	//przechowuje ID nodow z grid
};

struct Grid {  //siatka zlozona z elementow
public:
	int nN, nE;
	Node *nodes;
	Element *elements;

	Grid(int a, int b) {
		nN = a;
		nE = b;
		nodes = new Node[a];
		elements = new Element[b];
	}

	void fillGrid(Grid *g, GlobalData gd);
};

struct element_uni{ //przechowuje uniwersalne wartosci
public:
	//tablice:
	double **N, **dN_dKsi, **dN_dEta;
	//Punkty calkowania:
	double **pc;
	//Wagi calkowania:
	double **wc;
	//dodatkowo:
	int n_pc; //ilosc zmiennych w pkt calkowania (2=xy, 3=xyz)
	//punkty calkowania do warunkow brzegowych
	double **pc_bc;
	//wartosci funkcji ksztaltu w pc_bc
	double **Nbc;
	void makeElement_Uni();
	void fillElement_Uni();
};

struct local {
	int n_pc;
	Element *element;
	element_uni *el_uni;
	GlobalData *data;
	
	double *dx_dEta, *dy_dEta, *dx_dKsi, *dy_dKsi;

	double *detJ, **J;

	double **dN_dx, **dN_dy, ***Nx_NxT, ***Ny_NyT;

	double ***Hpc, **H;

	double **C;

	void makeLocal(element_uni*,GlobalData*);
	double** calculateH(Grid, int);
	double** calculateC();
	double** calculateHbc(Grid ,int);
	double* CalculateP(Grid, GlobalData, int);

	void printMatrix2D(double **m);

	//domyslne konstruktory
	local(int a) {
		n_pc = a;
	}
	local() {
		n_pc = 2;
	}
};