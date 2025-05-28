#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>     
#include <time.h>
#include<vector>
#include <fstream>
#include <unordered_map>
#include <filesystem>
#include <cstdlib>

using namespace std;

struct site4 {
    int position[4];
    int value_direction[4];
};

struct plaquette4 {
    int position[4];
    int value[6];
};



void lattice4_creation(vector<site4> &lattice, int n, int nt);
void lattice4Random_creation(vector<site4> &lattice, int n, int nt);
void plaquette4_creation(vector<plaquette4> &plaquettes, vector<site4> &lattice, int lattice_size, int n, int nt);
void plaquette4Random_creation(vector<plaquette4> &plaquettes, vector<site4> &lattice, int lattice_size, int n, int nt);
double sum4_all_plaquettes (vector<plaquette4> &plaquettes, int plaquettes_size);
void minimize4_energy(vector<plaquette4> &plaq, int plaquettes_size ,vector<site4> &lattice, int lattice_size, int n, int nt, double beta, int iterations, string fileE, string fileP, int &numberWilson, int Wilson_lenght);
void Wilson(int T, int R, vector<site4> &lattice, int n, int nt, double beta, string file);
void Polyakov(int T, vector<site4> &lattice, int lattice_size, int n, int nt, double beta, string file);
