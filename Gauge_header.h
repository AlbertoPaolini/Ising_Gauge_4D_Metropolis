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
    vector< int > plaquette_attached[4];
};

void lattice4_creation(vector<site4> &lattice, int n, int nt, vector< vector<int> > &map_nn);
void lattice4Random_creation(vector<site4> &lattice, int n, int nt, vector< vector<int> > &map_nn);
void plaquette4_creation(vector<int> &plaquettes, vector<site4> &lattice, int lattice_size, vector< vector<int> > &map_nn);
void test_map(vector<site4> &lattice, vector< vector<int> > &map_nn, int n, int nt);
double sum4_all_plaquettes (vector<int> &plaquettes, int plaquettes_size);
void minimize4_energy(vector< vector<int> > &map_nn, vector<int> &plaq, int plaquettes_size ,vector<site4> &lattice, int lattice_size, int n, int nt, double beta, int iterations, string fileE, string fileP, int &numberWilson, int Wilson_lenght);
void Wilson(int T, int R, vector<site4> &lattice, vector< vector<int> > &map_nn, double beta, int number);
void Polyakov(vector<site4> &lattice, int lattice_size, int n, int nt, double beta, string file);
