//
//                  ***GRADE_v1.00***
//  MyFunctions.hpp
//
//  Created by Farbod Mahmoudinobar and Cristiano L. Dias on 12/6/18.
//  Copyright © 2018 Farbod Mahmoudinobar and Cristiano L. Dias.
//
//  This file is part of GRADE.
//
//  GRADE is a free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  GRADE is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with GRADE.  If not, see <https://www.gnu.org/licenses/>.


#ifndef MyFunctions_hpp
#define MyFunctions_hpp

#include <stdio.h>

using namespace std;


void calc_Distance(int count_solvent, int count_solute, vector<vector<int>>& neigh_list, vector<vector<double>>& atom_Pos, double& boxX, double& boxY, double& boxZ, vector<int>& Nneigh, int Natoms, int topSolute, string time, double HBOND_DIST );

void ring_Finder(int count_solute, int Natoms, vector<int>& Nneigh, vector<vector<int>>& My_neigh, vector<vector<int>>& ring5_temp, vector<vector<int>>& ring6_temp, int firstSOL, int count_solvent, vector<vector<double>>& atom_Pos, double boxX, double boxY, double boxZ, double HBOND_DIST, double delta_p, double delta_h);

void find_shared_edges_ring5(int count_ring5, vector<vector<int>>& ring5, vector<vector<int>>& My_neigh_ring5, vector<unsigned long int>& N_ring5_neigh);

void find_shared_edges_ring6(int count_ring6, vector<vector<int>>& ring6, vector<vector<int>>& My_neigh_ring6, vector<unsigned long int>& N_ring6_neigh);

void find_shared_edges_ring6_ring5(int count_ring6, int count_ring5, vector<vector<int>>& ring5, vector<vector<int>>& ring6, vector<vector<int>>& My_neigh_ring6_ring5, vector<unsigned long int>& N_ring6_ring5_neigh );

bool  compare(vector<int>& ringA, vector<int>& ringB, int ringA_N, int ringB_N );

bool  compare_adjacant(vector<int>& ringA, vector<int>& ringB, int ringA_N, int ringB_N , vector<int>& base);

int remove_duplicates_map(vector<vector<int>>& cup512);

int remove_duplicates_map_rings(vector<vector<int>>& cups);

void coplanar_Points(vector<vector<int>>& ring , vector<vector<double>>& atoms, string time, vector<vector<int>>& ring_New, int THETA);

void cup_512_Finder(vector<vector<int>>& ring5,int count_ring5, vector<unsigned long int>& N_ring5_neigh, vector<vector<int>>& My_neigh_ring5, vector<vector<int>>& cup512 );

void cup_62512_Finder(vector<vector<int>>& ring6, int count_ring6, vector<vector<int>>& My_neigh_ring6_ring5, vector<unsigned long int> N_ring6_ring5_neigh, vector<vector<int>>& ring5, vector<unsigned long int> N_ring5_neigh, vector<vector<int>>& My_neigh_ring5, vector<vector<int>>& cup62512);

int cage_Finder(vector<vector<int>> cups, unsigned long count_cups, vector<vector<int>>& neighour_rings, vector<vector<int>>& cage, vector<vector<int>>& cage_rings, string time);

void print_vmd_cage_frings(vector<vector<int>> cups, int cage_count, vector<vector<int>> cage_rings, vector<vector<int>> ring5, vector<vector<int>> ring6, vector<vector<double>> atom_Pos, string time, string rawFilename , string box_size_xyz, vector<vector<double>> solutes, size_t & meth_counter, string solute1, int topSolute, string solute2, int count_solute2, int frameCounter);

double calc_F4(int count_solvent, int count_solute, vector<vector<int>>& My_neigh, vector<vector<double>>& atom_Pos, double& boxX, double& boxY, double& boxZ, vector<int>& Nneigh, int Natoms, int topSolute, string time, double HBOND_DIST );

int cage_Finder_64512(vector<vector<int>> cup62512, int count_62512_cups, vector<vector<int>>& cage_64512_rings);

void print_vmd_cage64512_frings(vector<vector<int>> cups, int cage_count, vector<vector<int>> cage_rings, vector<vector<int>> ring5, vector<vector<int>> ring6, vector<vector<double>> atom_Pos, string time, string rawFilename , string box_size_xyz, vector<vector<double>> solutes, size_t & meth_counter, string solute1, int topSolute, string solute2, int count_solute2, int frameCounter);

void print_usage();

#endif /* MyFunctions_hpp */
