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


#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "MyFunctions.hpp"

using namespace std;


    // New Function--------------------------------------------------------------------------------------

//using namespace std::chrono;

void calc_Distance(int count_solvent, int count_solute, vector<vector<int>>& My_neigh, vector<vector<double>>& atom_Pos, double& boxX, double& boxY, double& boxZ, vector<int>& Nneigh, int Natoms, int topSolute, string time, double HBOND_DIST ){
    
    string line, str1, str2, str3;  //str(i) and int(i) hold unused data.
    int neigh_counter=0;
    double dx, dy, dz, dist = 0.0;
    vector<double> temp_vect;
    vector<int> temp_vect2;
    
    
    
    //Calculate Pair Distances
    
    //cout << "solvent: " << count_solvent << ", solute: " << count_solute << "\n";
    
    temp_vect2={0};
    My_neigh.clear();
    Nneigh.clear();
    
    for (int i = 0 ; i < topSolute+1 ; i++){My_neigh.push_back(temp_vect2);}        //Fill My_neigh with zero's for solutes, i.e., make first topsolute lines of My_neigh file empty. This is needed for other calculations.
    
    temp_vect2.clear();
    
    for (int i = topSolute+1 ; i < topSolute + count_solvent ; i+=4)      // i is Oxygen index
    {
        neigh_counter = 0;
        temp_vect2.clear();
        temp_vect2.push_back(0);
        
        for (int j = topSolute+1 ; j < topSolute + count_solvent ; j+=4)       // j is 2nd Oxygen index
        {
            if(i != j)
            {
                dx = atom_Pos[j][0] - atom_Pos[i][0];
                dy = atom_Pos[j][1] - atom_Pos[i][1];
                dz = atom_Pos[j][2] - atom_Pos[i][2];
                
                
                if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                
                
                
                dist = sqrt(dx*dx + dy*dy + dz*dz);
                
                
                
                if (dist <= HBOND_DIST && dist > 0.18 )
                {
                    neigh_counter++;
                    
                    temp_vect2.push_back(j);
                }
            }
        }
        
        My_neigh.push_back(temp_vect2);
        Nneigh.push_back(neigh_counter);
        
        My_neigh.push_back({0});
        My_neigh.push_back({0});
        My_neigh.push_back({0});
        
        
    }
    
    // Writing output to a file. This is usefull for checking the results in early stages of writing the code. Mainly to make sure the atoms and rings are found correctly. After v1.12 these two files have been commented.
    

//    ofstream outFile;
//
//    int counter=0;
//
//    outFile.open("atom_pos.txt");
//    outFile << "t=\t" << time << "\n";
//
//
//    for (int i = 0  ; i < atom_Pos.size(); i++)
//    {
//        outFile << counter++ << ": ";
//        for (int j = 0 ; j < 3 ; j++)
//        {
//            outFile << atom_Pos[i][j] << ' ' ;
//        }
//        outFile << "\n";
//    }
//    outFile.close();
//
//
//    outFile.open("neigh_atoms.txt");
//    outFile << "t=\t" << time << "\n";
//
//
//    for (int i = 0 ; i < Nneigh.size() ; i++)
//    {
//        outFile << "Nneigh[" << topSolute+1 + (4*i) << "]= " << Nneigh[i] << "\n";
//        for (int j = 1 ; j < Nneigh[i]+1 ; j++)
//        {
//            if(Nneigh[i] != 0)
//            {
//                outFile << "My_neigh[" << topSolute+1 + (4*i) << "][" << j << "]= " << My_neigh[topSolute+1 + (4*i)][j] << "\n";
//            }
//        }
//        outFile << "\n";
//    }
//
//    outFile.close();
//

    
};



// New Function--------------------------------------------------------------------------------------
void ring_Finder(int count_solute, int Natoms, vector<int>& Nneigh, vector<vector<int>>& My_neigh, vector<vector<int>>& ring5_temp, vector<vector<int>>& ring6_temp, int topSolute, int count_solvent, vector<vector<double>>& atom_Pos, double boxX, double boxY, double boxZ, double HBOND_DIST, double delta_p, double delta_h){
    
    int q_2, q_3;
    int a_1, a_2, a_3, a_4, a_5, a_6;
    int k5=0, k6=0;
    vector<int> temp5, temp6 ;
    vector<vector<int>> m5, m6;
    vector<int> temp_vec5={0,0,0,0,0,0}, temp_vec6={0,0,0,0,0,0,0}, temp_vecn={0};
    vector<vector<int>> My_neigh_ring6, My_neigh_ring5, My_neigh_ring6_ring5;
    bool isPresent = false;
    double dist_a1_a4 = 0, dist_a1_a3 = 0 , dist_a2_a5 = 0, dist_a3_a6 = 0;
    double dx,dy,dz;
    //double delta1_p = 0.18;
    double Lrange5 = 1.6 * HBOND_DIST - delta_p;

    //double delta1_h = 0.3;
    double Lrange6 = 2 * HBOND_DIST - delta_h;
    
    
    My_neigh_ring5.clear();
    My_neigh_ring6.clear();
    My_neigh_ring6_ring5.clear();
    
    for( a_1 = topSolute+1; a_1 < (topSolute + count_solvent) ; a_1+=4)
    {
        for(int m_1 = 1 ; m_1 <= Nneigh[(a_1 -topSolute+1)/4] ; m_1++)
        {
            a_2 = My_neigh[a_1][m_1];
            if (a_2 != a_1)
            {
                for(int m_2 = 1 ; m_2 <= Nneigh[(a_2 -topSolute+1)/4] ; m_2++)
                {
                    a_3 = My_neigh[a_2][m_2];
                    
                    if(a_3 != a_1)
                    {
                        for(int m_3 = 1 ; m_3 <= Nneigh[(a_3 -topSolute+1)/4] ; m_3++)
                        {
                            a_4 = My_neigh[a_3][m_3];
                            if(a_4 != a_2 && a_4 != a_1)
                            {
                                for(int m_4 = 1 ; m_4 <= Nneigh[(a_4 -topSolute+1)/4] ; m_4++)
                                {
                                    a_5 = My_neigh[a_4][m_4];
                                    
                                    if(a_5 != a_3 && a_5 != a_2 && a_5 != a_1)
                                    {
                                        for(int m_5 = 1 ; m_5 <= Nneigh[(a_5 -topSolute+1)/4] ; m_5++)
                                        {
                                            
                                            for(q_2 = 1 ; q_2 <= Nneigh[(a_1 -topSolute+1)/4] ; q_2++)
                                            {
                                                if( a_5 == My_neigh[a_1][q_2-1])
                                                {
                                                    
                                                    //*************
                                                   /*This condition ensures that the ring is actualy a pentagon, not a collapsed ring.
                                                    This condition was added in v1.14. */
                                                    dx = atom_Pos[a_1][0] - atom_Pos[a_4][0];
                                                    dy = atom_Pos[a_1][1] - atom_Pos[a_4][1];
                                                    dz = atom_Pos[a_1][2] - atom_Pos[a_4][2];
                         
                                                    if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                                                    if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                                                    if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}

                                                    dist_a1_a4 = sqrt(dx*dx + dy*dy + dz*dz);
                                                    
                                                    //if (dist_a1_a4 > 0.37  && dist_a1_a4 < 0.53)
                                                    if( dist_a1_a4 > Lrange5  )
                                                    {
                                                        dx = atom_Pos[a_1][0] - atom_Pos[a_3][0];
                                                        dy = atom_Pos[a_1][1] - atom_Pos[a_3][1];
                                                        dz = atom_Pos[a_1][2] - atom_Pos[a_3][2];
                                                        
                                                        if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                                                        if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                                                        if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                                                        
                                                        dist_a1_a3 = sqrt(dx*dx + dy*dy + dz*dz);
                                                        //if (dist_a1_a3 > 0.37 && dist_a1_a3 < 0.53 )
                                                        if (dist_a1_a3 > Lrange5 )
                                                        {
                                                            dx = atom_Pos[a_2][0] - atom_Pos[a_5][0];
                                                            dy = atom_Pos[a_2][1] - atom_Pos[a_5][1];
                                                            dz = atom_Pos[a_2][2] - atom_Pos[a_5][2];
                                                            
                                                            if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                                                            if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                                                            if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                                                            
                                                            dist_a2_a5 = sqrt(dx*dx + dy*dy + dz*dz);
                                                            //if (dist_a2_a5 > 0.37  && dist_a2_a5 < 0.53)
                                                              if (dist_a2_a5 > Lrange5 )
                                                            {
                                                    //*************


                                                                temp5.clear();
                                                                
                                                                temp5.push_back(a_1);
                                                                temp5.push_back(a_2);
                                                                temp5.push_back(a_3);
                                                                temp5.push_back(a_4);
                                                                temp5.push_back(a_5);
                                                                
                                                                
                                                                
                                                                temp_vec5[0] = a_1;
                                                                temp_vec5[1] = a_2;
                                                                temp_vec5[2] = a_3;
                                                                temp_vec5[3] = a_4;
                                                                temp_vec5[4] = a_5;
                                                                temp_vec5[5] = a_1;
                                                                
                                                                sort(temp5.begin(), temp5.end());
                                                                
                                                                
                                                                if (k5==0) //put the first ring into matrix m5.
                                                                {
                                                                    m5.push_back(temp5);
                                                                    k5++;
                                                                    ring5_temp.push_back(temp_vec5);
                                                                }
                                                                
                                                                for (int i = 0 ; i < k5 ; i++)
                                                                {
                                                                    if(m5[i][0] == temp5[0] && m5[i][1] == temp5[1] && m5[i][2] == temp5[2] && m5[i][3] == temp5[3] && m5[i][4] == temp5[4])
                                                                    {
                                                                        isPresent=true ;
                                                                        break;
                                                                    }
                                                                    else isPresent=false;
                                                                }
                                                                
                                                                if (isPresent == false)
                                                                {
                                                                    m5.push_back(temp5);
                                                                    k5++;
                                                                    ring5_temp.push_back(temp_vec5);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }

                                            a_6 = My_neigh[a_5][m_5];
                                            
                                            if(a_6 != a_5)
                                            {
                                                for(q_3 = 1 ; q_3 <= Nneigh[(a_1 - topSolute+1)/4] ; q_3++)
                                                {
                                                    if( a_6!=a_1 && a_3!=a_2 &&  a_6!=a_2  && a_4!=a_3 && a_6!=a_3 && a_5!=a_4 && a_6!=a_4 && a_5!=a_6 && a_6==My_neigh[a_1][q_3-1] )
                                                    {
                                                        //*************
                                                        /*This condition ensures that the ring is actualy a hexagon, not a collapsed ring. It measures the 3 diagonal distances a1-a4, a2-a5 and a3-a6.                                                                 This condition was added in v1.14. */
                                                        dx = atom_Pos[a_1][0] - atom_Pos[a_4][0];
                                                        dy = atom_Pos[a_1][1] - atom_Pos[a_4][1];
                                                        dz = atom_Pos[a_1][2] - atom_Pos[a_4][2];
                                                        
                                                        if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                                                        if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                                                        if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                                                        
                                                        dist_a1_a4 = sqrt(dx*dx + dy*dy + dz*dz);
                                                        
                                                        //if (dist_a1_a4 > 0.4 && dist_a1_a4 < 0.75 )
                                                        if (dist_a1_a4 > Lrange6  )
                                                        {
                                                            dx = atom_Pos[a_2][0] - atom_Pos[a_5][0];
                                                            dy = atom_Pos[a_2][1] - atom_Pos[a_5][1];
                                                            dz = atom_Pos[a_2][2] - atom_Pos[a_5][2];
                                                            
                                                            if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                                                            if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                                                            if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                                                            
                                                            dist_a2_a5 = sqrt(dx*dx + dy*dy + dz*dz);
                                                            //if (dist_a2_a5 > 0.4 && dist_a2_a5 < 0.75 )
                                                            if (dist_a2_a5 > Lrange6 )
                                                            {
                                                                dx = atom_Pos[a_3][0] - atom_Pos[a_6][0];
                                                                dy = atom_Pos[a_3][1] - atom_Pos[a_6][1];
                                                                dz = atom_Pos[a_3][2] - atom_Pos[a_6][2];
                                                                
                                                                if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                                                                if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                                                                if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                                                                
                                                                dist_a3_a6 = sqrt(dx*dx + dy*dy + dz*dz);
                                                                //if (dist_a3_a6 > 0.4 && dist_a3_a6 < 0.75 )
                                                                if (dist_a3_a6 > Lrange6 )
                                                                {
                                                                    
                                                        //*************
                                                        
                                                        
                                                        temp6.clear();
                                                        temp6.push_back(a_1);
                                                        temp6.push_back(a_2);
                                                        temp6.push_back(a_3);
                                                        temp6.push_back(a_4);
                                                        temp6.push_back(a_5);
                                                        temp6.push_back(a_6);
                                                        
                                                        
                                                        temp_vec6[0]=a_1;
                                                        temp_vec6[1]=a_2;
                                                        temp_vec6[2]=a_3;
                                                        temp_vec6[3]=a_4;
                                                        temp_vec6[4]=a_5;
                                                        temp_vec6[5]=a_6;
                                                        temp_vec6[6]=a_1;
                                                        
                                                        sort(temp6.begin(), temp6.end());
                                                        
                                                        if (k6==0) //put the first ring into matrix m6.
                                                        {
                                                            m6.push_back(temp6);
                                                            ring6_temp.push_back(temp_vec6);
                                                        }
                                                        
                                                        
                                                        for (int i=0; i<=k6 ; i++)
                                                        {
                                                            
                                                            if(m6[i][0] == temp6[0] && m6[i][1] == temp6[1] && m6[i][2] == temp6[2] && m6[i][3] == temp6[3] && m6[i][4] == temp6[4] && m6[i][5] == temp6[5])
                                                            {
                                                                
                                                                isPresent=true ;
                                                                break;
                                                                
                                                            }
                                                            
                                                            else isPresent=false;
                                                            
                                                        }
                                                        
                                                        if (isPresent == false)
                                                        {
                                                            
                                                            m6.push_back(temp6);
                                                            
                                                            k6++;
                                                            
                                                            ring6_temp.push_back(temp_vec6);
                                                            
                                                        }
                                                    }
                                                }}}
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
   
    
//    ofstream out1, out2;
//    out1.open("ring5-temp.txt");
//    out2.open("ring6-temp.txt");
//
//    for (int i = 0 ; i < ring5_temp.size() ; i++)
//    {
//        out1 << i << "ring5 ";
//        out2 << i << "ring6 ";
//
//        for (int j = 0 ; j < ring5_temp[i].size(); j++)
//        {
//            out1 << ring5_temp[i][j] << " " ;
//            out2 << ring6_temp[i][j] << " " ;
//        }
//        out1 << "\n";
//        out2 << "\n";
//    }
//    out1.close();
//    out2.close();
    
}




// New Function--------------------------------------------------------------------------------------
void find_shared_edges_ring5(int count_ring5, vector<vector<int>>& ring5, vector<vector<int>>& My_neigh_ring5, vector<unsigned long int>& N_ring5_neigh){
    
    // Finding shared edges for 5-membered rings
    
    vector<int> temp_vecn;
    temp_vecn.clear();
    //int counter=0;
    int neigh_ring_counter=0;       //Number of ring-neighbours of each ring.
                                    //vector<int> N_ring5_neigh;      //Vector of ring-neighbours of all rings.
    vector<vector<int>> ring5_neigh;        //List of ring-neighbours of all rings.
    vector<int> temp_vec4={0,0,0,0};
    
    
    for (int i = 0 ; i < count_ring5 ; i++)
    {
        for (int j = 0 ; j < ring5[i].size() ; j++)
        {
            for (int k = 0 ; k < count_ring5 ; k++)
            {
                for (int l=0 ; l < ring5[k].size() ; l++)
                {
                    if ( i != k && j+1 < ring5[i].size() && ( ( ring5[i][j] == ring5[k][l] && ring5[i][j+1] == ring5[k][l+1]) || (ring5[i][j] == ring5[k][l+1] && ring5[i][j+1] == ring5[k][l]) ) )
                    {
                        neigh_ring_counter++;
                        temp_vec4[0] = i+1;
                        temp_vec4[1] = k+1;
                        temp_vec4[2] = ring5[i][j];
                        temp_vec4[3] = ring5[i][j+1];
                        
                        ring5_neigh.push_back(temp_vec4);
                        
                        temp_vecn.push_back(k);
                    }
                }
            }
        }
        N_ring5_neigh.push_back(neigh_ring_counter);
        neigh_ring_counter = 0;
        
        sort(temp_vecn.begin(), temp_vecn.end());
        temp_vecn.erase(unique(temp_vecn.begin(), temp_vecn.end()), temp_vecn.end());
        My_neigh_ring5.push_back(temp_vecn);
        temp_vecn.clear();
    }
    
    N_ring5_neigh.clear();
    for (int i = 0 ; i < count_ring5 ; i++)
    {
        N_ring5_neigh.push_back(My_neigh_ring5[i].size()) ;
    }
    
    
    
}

// New Function--------------------------------------------------------------------------------------

void find_shared_edges_ring6(int count_ring6, vector<vector<int>>& ring6, vector<vector<int>>& My_neigh_ring6, vector<unsigned long int>& N_ring6_neigh){
    
    vector<int> temp_vecn;
    temp_vecn.clear();
    int counter=0;
    int neigh_ring_counter=0;
    //vector<int> N_ring6_neigh;
    vector<vector<int>> ring6_neigh;
    vector<int> temp_vec4 = {0,0,0,0};
    
    
    for (int i = 0 ; i < count_ring6 ; i++)
    {
        for (int j = 0 ; j < ring6[i].size() ; j++)
        {
            for (int k = 0 ; k < count_ring6 ; k++)
            {
                for (int l=0 ; l < ring6[k].size() ; l++)
                {
                    if (i != k && j+1 < ring6[i].size() && ( ( ring6[i][j] == ring6[k][l] && ring6[i][j+1] == ring6[k][l+1]) || (ring6[i][j] == ring6[k][l+1] && ring6[i][j+1] == ring6[k][l]) ))
                    {
                        counter++;
                        neigh_ring_counter++;
                        temp_vec4[0] = i+1;
                        temp_vec4[1] = k+1;
                        temp_vec4[2] = ring6[i][j];
                        temp_vec4[3] = ring6[i][j+1];
                        
                        ring6_neigh.push_back(temp_vec4);
                        
                        temp_vecn.push_back(k);
                    }
                }
            }
        }
        N_ring6_neigh.push_back(neigh_ring_counter);
        neigh_ring_counter = 0;
        
        sort(temp_vecn.begin(), temp_vecn.end());
        temp_vecn.erase(unique(temp_vecn.begin(), temp_vecn.end()), temp_vecn.end());
        My_neigh_ring6.push_back(temp_vecn);
        temp_vecn.clear();
        
    }
    
    N_ring6_neigh.clear();
    for (int i = 0 ; i < count_ring6 ; i++)
    {
        N_ring6_neigh.push_back(My_neigh_ring6[i].size()) ;
    }
    
}

// New Function--------------------------------------------------------------------------------------

void find_shared_edges_ring6_ring5(int count_ring6, int count_ring5, vector<vector<int>>& ring5, vector<vector<int>>& ring6, vector<vector<int>>& My_neigh_ring6_ring5, vector<unsigned long int>& N_ring6_ring5_neigh ){
    
    vector<int> temp_vecn;
    temp_vecn.clear();
    int counter=0;
    int neigh_ring_counter=0;
    //vector<int> N_ring6_ring5_neigh;
    vector<vector<int>> ring6_ring5_neigh;
    vector<int> temp_vec4 = {0,0,0,0};
    
    
    for (int i = 0 ; i < count_ring6 ; i++)
    {
        for (int j = 0 ; j < ring6[i].size() ; j++)
        {
            for (int k = 0 ; k < count_ring5 ; k++)
            {
                for (int l = 0 ; l < ring5[k].size() ; l++)
                {
                    if ( i != k && j+1 < ring6[i].size() && ( (ring6[i][j] == ring5[k][l] && ring6[i][j+1] == ring5[k][l+1]) || (ring6[i][j] == ring5[k][l+1] && ring6[i][j+1] == ring5[k][l]) )  )
                    {
                        counter++;
                        neigh_ring_counter++;
                        temp_vec4[0] = i+1;
                        temp_vec4[1] = k+1;
                        temp_vec4[2] = ring6[i][j];
                        temp_vec4[3] = ring6[i][j+1];
                        
                        ring6_ring5_neigh.push_back(temp_vec4);
                        
                        temp_vecn.push_back(k);
                        
                    }
                }
            }
        }
        N_ring6_ring5_neigh.push_back(neigh_ring_counter);
        neigh_ring_counter = 0;
        
        sort(temp_vecn.begin(), temp_vecn.end());
        temp_vecn.erase(unique(temp_vecn.begin(), temp_vecn.end()), temp_vecn.end());
        My_neigh_ring6_ring5.push_back(temp_vecn);
        temp_vecn.clear();
        
    }
    
    N_ring6_ring5_neigh.clear();
    for (int i = 0 ; i < count_ring6 ; i++)
    {
        N_ring6_ring5_neigh.push_back(My_neigh_ring6_ring5[i].size()) ;
    }
    
}


// New Function--------------------------------------------------------------------------------------

// New Function--------------------------------------------------------------------------------------

//This function compares two rings for shared sides. Each side is made up of two consecutive atoms.
//Every two consecutive atoms from ringA are compared with every two consecutive atoms of ringB.
//For example, if ringA={1,3,5,7,9} and ringB={2,4,6,8,10}, {1,3} will be compared with both {2,4} and {4,2} since both cases should be considered.

bool  compare(vector<int>& ringA, vector<int>& ringB, int ringA_N, int ringB_N )
{
    int neigh_ring_counter=0;
    
    for (int i = 0 ; i < ringA_N  ; i++)
    {
        for (int j = 0 ; j < ringB_N  ; j++)
        {
            
            if  ( (ringA[i] == ringB[j] && ringA[i+1] == ringB[j+1]) || (ringA[i] == ringB[j+1] && ringA[i+1] == ringB[j]) )
            {
                neigh_ring_counter++;
                break;
            }
        }
    }
    if (neigh_ring_counter == 0)
    {
        return false;
    }
    else return true;
}       //close compare function

// New Function---------------------------------------------------------------------------------------

// This functions compares 3 rings. If two 5-rings have a common edge/side, then it compares the common side with the sides of base ring(6-ring).
// If the common side b/w two 5-rings are shared with the base ring, these two 5-rings are not neighbors and function returns "false", otherwise
// it returns "true".
bool  compare_adjacant(vector<int>& ringA, vector<int>& ringB, int ringA_N, int ringB_N , vector<int>& base)
{
    int neigh_ring_counter=0;
    
    for (int i = 0 ; i < ringA_N  ; i++)
    {
        for (int j = 0 ; j < ringB_N ; j++)
        {
            if  ( (ringA[i] == ringB[j] && ringA[i+1] == ringB[j+1]) || (ringA[i] == ringB[j+1] && ringA[i+1] == ringB[j] ) )       //compare sides of the two 5-rings. If they have a common side, compare the common side with base ring.
            {
                for (int k6 = 0 ; k6 < 6; k6++)     // Loop to go over all the sides of base ring.
                {
                    if( (ringA[i] == base[k6] && ringA[i+1] == base[k6+1])  || (ringA[i] == base[k6+1] && ringA[i+1] == base[k6]) )
                    {
                        return false;
                    }
                }
                neigh_ring_counter++;
                break;
            }
        }
    }
    if (neigh_ring_counter == 0)
    {
        return false;
    }
    else return true;
}       //close compare function

// New Function---------------------------------------------------------------------------------------


/* This Functions reads in a vector<vector<int>> which holds all the cages/cups/rings and removes the duplicates from them. Then it
 writes the new list without duplicates inside the original input file. This function returns an "int" which is the count number of
 the cups/rings without duplicates. Pay attention that the after this function is executed, the original input variable still holds the same number of rows/columns as before, but the first "Count" number of rows hold the results without any duplicate. So it is crucial to use the first "Count" rows in future functions. 
 */

int remove_duplicates_map(vector<vector<int>>& cups)
{
    unsigned long int cup_size;
    if(cups.size() != 0){cup_size = cups[0].size();}         //Number of rings in each cup, including base ring.
    else cup_size = 0;
    
    unsigned long int N_cups = cups.size();              //Number of cups before removing the duplicates.
    int count=0;                                                //Number of cups.
    vector<int> current;                                        //Cup which is being read and compared at each loop.
    map<vector<int>,int> mymap;                                 //A map which counts the frequencty of each cup( [cup,frequency] )
    int prev = 0 ;
    
    for (int i = 0 ; i < N_cups ; i++)      //Loop over all the cups.
    {
        for(int j = 0 ; j < cup_size ; j++){current.push_back(cups[i][j]);}      //Loop to put the read cup into "current".
        sort(current.begin(), current.end());
        int duplicate_count= 0;
        for(int j = 0 ; j < current.size()-1 ; j++ )
        {
            if(current[j] == current[j+1]) {duplicate_count++;}
        }
        
        if(duplicate_count == 0)
        {
            if(!mymap[current])                  //For each read cup, check if it exists in the map. If it does NOT, write it in the map.
                
            {
                
                mymap[current] = 1;
                count++;
                for (int k = 0 ; k < cup_size ; k++)
                {
                    cups[prev][k] = cups[i][k];   //If the cup does not exist in map, write it to first line of original cup variable.
                }
                prev++;
            }
            else {mymap[current]++;}        //If the cup already exists in the map, jump the writing step and add one to frequencty of that cup.
        }
        else duplicate_count = 0;
            
            current.clear();
    }
    
    return count;
}       // close remove_dupclicates_maps

// New Function---------------------------------------------------------------------------------------

//Start of remove duplicates for rings
int remove_duplicates_map_rings(vector<vector<int>>& cups)
{
    unsigned long int cup_size;
    if(cups.size() != 0){cup_size = cups[0].size();}         //Number of rings in each cup, including base ring.
    else cup_size = 0;
    
    unsigned long int N_cups = cups.size();              //Number of cups before removing the duplicates.
    int count=0;                                                //Number of cups.
    vector<int> current;                                        //Cup which is being read and compared at each loop.
    map<vector<int>,int> mymap;                                 //A map which counts the frequencty of each cup( [cup,frequency] )
    int prev = 0 ;
    
    for (int i = 0 ; i < N_cups ; i++)      //Loop over all the cups.
    {
        for(int j = 0 ; j < cup_size ; j++){current.push_back(cups[i][j]);}      //Loop to put the read cup into "current".
        sort(current.begin(), current.end());
        int duplicate_count= 0;
        for(int j = 1 ; j < current.size()-1 ; j++ )
        {
            if(current[j] == current[j+1]) {duplicate_count++;}
        }
        if(duplicate_count <= 1)
        {
            if(!mymap[current])                  //For each read cup, check if it exists in the map. If it does NOT, write it in the map.
                
            {
                
                mymap[current] = 1;
                count++;
                for (int k = 0 ; k < cup_size ; k++)
                {
                    cups[prev][k] = cups[i][k];   //If the cup does not exist in map, write it to first line of original cup variable.
                }
                prev++;
            }
            else {mymap[current]++;}        //If the cup already exists in the map, jump the writing step and add one to frequencty of that cup.
        }
        else duplicate_count = 0;
        
        current.clear();
    }
    
    return count;
}

//End of remove duplicates for rings

// New Function---------------------------------------------------------------------------------------

// New Function--------------------------------------------------------------------------------------

/*
 This function gets the ring list and atom positions as input and checks all the rings to be in one plane. It names the pentagon/hexagon in clockwise direction from A to E (pentagon) or F (hexagon) and creates vectors AB, AC, AD, AE and AF, calculates cross product of AB*AC, AD*AE, AF*AE and findes the angle between these vectors. If the angle is in an appropriate range, then the ring is in one plane. If the ring is in one plane, the ring will be written to the ring list, otherwise it will be removed.
 */

//To print the debugging output files, un-comment lines including "outFile".

//START OF NEW COPLANAR FUNCTION
void coplanar_Points_test(vector<vector<int>>& ring , vector<vector<double>>& atoms, string time, vector<vector<int>>& ring_New, int THETA)
{
#include <math.h>
    
    
    vector<double> A, B, C, D, E, F;    //Vertices of each ring and their positions(A=ai+aj+ak, etc).
    vector<double> AB, AC, AD, AE, AF, FD, DE;  //Edges of each ring as vectors(AB=B-A, etc).
    vector<double> ABC, ADE, AFE, ACD, ADF, FDE;       //Crossproduct results(ABC=AB*AC, etc).
    vector<int> temp_ring, temp_vec;
    
    unsigned long int N_rings = ring.size();
    unsigned long int N_vertices = ring[0].size();
    
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    
    
    double magnitude_ABC, magnitude_ADE, magnitude_AFE, magnitude_ACD, magnitude_ADF, magnitude_FDE;
    double theta1 = 0;
    double theta2 = 0;
    double theta3 = 0;
    double theta4 = 0;
    double theta5 = 0;
    
    const double PI = 3.14159265;
    
    int counter = 0;
    
    for (int r = 0 ; r < N_rings ; r++)
    {
        //Put atom positions in vectors A to F.
        A.push_back(atoms[ ring[r][0] ][0]);
        A.push_back(atoms[ ring[r][0] ][1]);
        A.push_back(atoms[ ring[r][0] ][2]);
        
        B.push_back(atoms[ ring[r][1] ][0]);
        B.push_back(atoms[ ring[r][1] ][1]);
        B.push_back(atoms[ ring[r][1] ][2]);
        
        C.push_back(atoms[ ring[r][2] ][0]);
        C.push_back(atoms[ ring[r][2] ][1]);
        C.push_back(atoms[ ring[r][2] ][2]);
        
        D.push_back(atoms[ ring[r][3] ][0]);
        D.push_back(atoms[ ring[r][3] ][1]);
        D.push_back(atoms[ ring[r][3] ][2]);
        
        E.push_back(atoms[ ring[r][4] ][0]);
        E.push_back(atoms[ ring[r][4] ][1]);
        E.push_back(atoms[ ring[r][4] ][2]);
        
        
        //Create vectors of edges.
        AB.push_back(B[0]-A[0]);
        AB.push_back(B[1]-A[1]);
        AB.push_back(B[2]-A[2]);
        
        AC.push_back(C[0]-A[0]);
        AC.push_back(C[1]-A[1]);
        AC.push_back(C[2]-A[2]);
        
        AD.push_back(D[0]-A[0]);
        AD.push_back(D[1]-A[1]);
        AD.push_back(D[2]-A[2]);
        
        AE.push_back(E[0]-A[0]);
        AE.push_back(E[1]-A[1]);
        AE.push_back(E[2]-A[2]);
        
        
        //Calculate cross products.
        
        ABC.push_back(AB[1]*AC[2] - AB[2]*AC[1]);
        ABC.push_back(AB[2]*AC[0] - AB[0]*AC[2]);
        ABC.push_back(AB[0]*AC[1] - AB[1]*AC[0]);
        
        ACD.push_back(AC[1]*AD[2] - AC[2]*AD[1]);
        ACD.push_back(AC[2]*AD[0] - AC[0]*AD[2]);
        ACD.push_back(AC[0]*AD[1] - AC[1]*AD[0]);
        
        ADE.push_back(AD[1]*AE[2] - AD[2]*AE[1]);
        ADE.push_back(AD[2]*AE[0] - AD[0]*AE[2]);
        ADE.push_back(AD[0]*AE[1] - AD[1]*AE[0]);
        
        //Calculate the angle between ABC and ADE from their dot product. call it theta1. if it -<theta1<+, then the ring is in one plane.
        
        for (int i = 0 ; i < 3 ; i++)       //Calculate theta1, angle between Abc and ACD.
        {
            sum1 += ABC.at(i) * ACD.at(i);
        }
        
        magnitude_ABC = sqrt(pow(ABC.at(0), 2) + pow(ABC.at(1), 2) + pow(ABC.at(2), 2));
        magnitude_ACD = sqrt(pow(ACD.at(0), 2) + pow(ACD.at(1), 2) + pow(ACD.at(2), 2));
        
        theta1 = acos( sum1 / (magnitude_ABC * magnitude_ACD) ) * 180.0 / PI;       //convert the result to degrees.
        
        //START OF NEW COMPARISON
        for (int i = 0 ; i < 3 ; i++)       //Calculate theta2, angle between ACD and ADE.
        {
            sum2 += ACD.at(i) * ADE.at(i);
        }
        
        magnitude_ACD = sqrt(pow(ACD.at(0), 2) + pow(ACD.at(1), 2) + pow(ACD.at(2), 2));
        magnitude_ADE = sqrt(pow(ADE.at(0), 2) + pow(ADE.at(1), 2) + pow(ADE.at(2), 2));
        
        theta2 = acos( sum2 / (magnitude_ADE * magnitude_ACD) ) * 180.0 / PI;       //convert the result to degrees.
                                                                                    //END OF NEW COMPARISON
        
        if(N_vertices ==7)      //Perform these only if rings with 6 members are provided.
        {
            
            F.push_back(atoms[ ring[r][5] ][0]);    //Put 6-th point in vector F only if rings with 6 members are provided as input.
            F.push_back(atoms[ ring[r][5] ][1]);
            F.push_back(atoms[ ring[r][5] ][2]);
            
            AF.push_back(F[0]-A[0]);
            AF.push_back(F[1]-A[1]);
            AF.push_back(F[2]-A[2]);
            
            FD.push_back(F[0]-D[0]);
            FD.push_back(F[1]-D[1]);
            FD.push_back(F[2]-D[2]);
            
            DE.push_back(D[0]-E[0]);
            DE.push_back(D[1]-E[1]);
            DE.push_back(D[2]-E[2]);
            
            AFE.push_back(AE[1]*AF[2] - AE[2]*AF[1]);
            AFE.push_back(AE[2]*AF[0] - AE[0]*AF[2]);
            AFE.push_back(AE[0]*AF[1] - AE[1]*AF[0]);
            
            for (int i = 0 ; i < 3 ; i++)       //Calculate theta3, angle between ADE and AFE from their dot product.
            {
                sum3 += ADE.at(i) * AFE.at(i);
            }
            
            magnitude_AFE = sqrt(pow(AFE.at(0), 2) + pow(AFE.at(1), 2) + pow(AFE.at(2), 2));
            
            theta3 = acos( sum3 / (magnitude_ADE * magnitude_AFE) ) * 180.0 / PI;       //Theta3 is the angle between ADE.AFE.
            
            
        }       //End of loop of 6 member rings calculations.
        
        
        if(N_vertices == 6)     //replacing original ring5 with planar rings.
        {
            
            if( theta1 < THETA  && theta2 < THETA)      //If a ring is planar, do these:
            {
                for (int i = 0 ; i < N_vertices ; i++)
                {
                    temp_ring.push_back(ring[r][i]);    //Put the ring atoms in temp_ring.
                    temp_vec.push_back(temp_ring[0]);
                    temp_ring.clear();
                }
                ring_New.push_back(temp_vec);
                temp_vec.clear();
                counter++;      //counter of planar rings.
            }
        }
        
        else if(N_vertices == 7 )       //replacing original ring6 with planar rings.
        {
            if ( theta1 < THETA && theta2 < THETA && theta3 < THETA )
            {
                for (int i = 0 ; i < N_vertices ; i++)
                {
                    temp_ring.push_back(ring[r][i]);    //Put the ring atoms in temp_ring.
                    temp_vec.push_back(temp_ring[0]);
                    temp_ring.clear();
                }
                ring_New.push_back(temp_vec);
                temp_vec.clear();
                counter++;      //counter of planar rings.
            }
        }
        //Zero all the variables used in each loop.
        
        sum1 = 0;
        sum2 = 0;
        sum3 = 0;
        sum4 = 0;
        sum5 = 0;
        theta1 = 0;
        theta2 = 0;
        theta3 = 0;
        theta4 = 0;
        theta5 = 0;
        A.clear();
        B.clear();
        C.clear();
        D.clear();
        E.clear();
        F.clear();
        AB.clear();
        AC.clear();
        AD.clear();
        AE.clear();
        AF.clear();
        FD.clear();
        DE.clear();
        ABC.clear();
        ACD.clear();
        ADE.clear();
        AFE.clear();
        ADF.clear();
        FDE.clear();
        magnitude_ABC = 0;
        magnitude_ADE = 0;
        magnitude_AFE = 0;
        magnitude_ACD = 0;
        magnitude_ADF = 0;
        magnitude_FDE = 0;
        
        
    }       // End of for loop for each ring.
    
}
//END OF NEW COPLANAR FUNCTION
// New Function-------------------------------------------------------------------------------------------

/*
 This Function finds cups (half of cage) of 5(12) cages. Sometimes cups and half_cups are used interchangeably("half" reminds me it is not a full cage and distinguishes them).
 It gets ring5 and all related info about ring5 and gives back cup512.
 */
void cup_512_Finder(vector<vector<int>>& ring5,int count_ring5, vector<unsigned long int>& N_ring5_neigh, vector<vector<int>>& My_neigh_ring5, vector<vector<int>>& cup512 ){
    
    // finding half of 5(12) cages
    
    vector<int> temp_vecn, temp_vec5={0,0,0,0,0,0,};
    temp_vecn.clear();
    int first_ring, second_ring, third_ring, fourth_ring, fifth_ring, sixth_ring;       //These refer to first, second, third, ... ring of cup. First ring is the base of the cup.
    vector<int> temp_ring1, temp_ring2, temp_ring3, temp_ring4, temp_ring5, temp_ring6, temp_ring7;
    vector<int> repeated_atoms;
    
    for (int index1 = 0 ; index1 < count_ring5 ; index1++)
    {
        first_ring = index1;
        if (N_ring5_neigh[first_ring] > 4)
        {            
            for (int index2 = 0 ; index2 < N_ring5_neigh[first_ring] ; index2++)
            {
                second_ring = My_neigh_ring5[first_ring][index2];
                if (N_ring5_neigh[second_ring] > 2)
                {
                    if(second_ring != first_ring)
                    {
                        for (int index3 = 0 ; index3 < N_ring5_neigh[second_ring] ; index3++)
                        {
                            if (My_neigh_ring5[second_ring][index3] != first_ring &&
                                My_neigh_ring5[second_ring][index3] != second_ring
                                )
                            {
                                if(compare(ring5[first_ring], ring5[My_neigh_ring5[second_ring][index3]], 5, 5))
                                {
                                    third_ring = My_neigh_ring5[second_ring][index3];
                                    
                                    if (N_ring5_neigh[third_ring] > 2)
                                    {
                                        for (int index4 = 0 ; index4 < N_ring5_neigh[third_ring]; index4++)
                                        {
                                            if (My_neigh_ring5[third_ring][index4] != first_ring &&
                                                My_neigh_ring5[third_ring][index4] != second_ring &&
                                                My_neigh_ring5[third_ring][index4] != third_ring &&
                                                compare(ring5[first_ring], ring5[My_neigh_ring5[third_ring][index4]], 5, 5) &&
                                                compare_adjacant(ring5[third_ring], ring5[My_neigh_ring5[third_ring][index4]], 5, 5, ring5[first_ring])
                                                )
                                            {
                                                fourth_ring = My_neigh_ring5[third_ring][index4];
                                                if (N_ring5_neigh[fourth_ring] > 2)
                                                {
                                                    for (int index5 = 0 ; index5 < N_ring5_neigh[fourth_ring]; index5++)
                                                    {
                                                        if (My_neigh_ring5[fourth_ring][index5] != first_ring &&
                                                            My_neigh_ring5[fourth_ring][index5] != second_ring &&
                                                            My_neigh_ring5[fourth_ring][index5] != third_ring  &&
                                                            My_neigh_ring5[fourth_ring][index5] != fourth_ring &&
                                                            compare(ring5[first_ring], ring5[My_neigh_ring5[fourth_ring][index5]], 5, 5) &&
                                                            compare_adjacant(ring5[fourth_ring], ring5[My_neigh_ring5[fourth_ring][index5]], 5, 5, ring5[first_ring])
                                                            )
                                                        {
                                                            fifth_ring = My_neigh_ring5[fourth_ring][index5];
                                                            
                                                            if (N_ring5_neigh[fifth_ring] != 0 && N_ring5_neigh[fifth_ring] > 2)
                                                            {
                                                                for (int index6 = 0 ; index6 < N_ring5_neigh[fifth_ring]; index6++)
                                                                {
                                                                    if (My_neigh_ring5[fifth_ring][index6] != first_ring  &&
                                                                        My_neigh_ring5[fifth_ring][index6] != second_ring &&
                                                                        My_neigh_ring5[fifth_ring][index6] != third_ring  &&
                                                                        My_neigh_ring5[fifth_ring][index6] != fourth_ring &&
                                                                        My_neigh_ring5[fifth_ring][index6] != fifth_ring  &&
                                                                        compare(ring5[second_ring], ring5[My_neigh_ring5[fifth_ring][index6]], 5, 5) &&
                                                                        compare_adjacant(ring5[fifth_ring], ring5[My_neigh_ring5[fifth_ring][index6]], 5, 5, ring5[first_ring])
                                                                        )
                                                                    {
                                                                        sixth_ring = My_neigh_ring5[fifth_ring][index6];
                                                                        
                                                                        if (compare(ring5[first_ring], ring5[sixth_ring], 5, 5))
                                                                        {
                                                                            
                                                                            temp_ring1.clear();
                                                                            temp_ring2.clear();
                                                                            temp_ring3.clear();
                                                                            temp_ring4.clear();
                                                                            temp_ring5.clear();
                                                                            temp_ring6.clear();
                                                                            
                                                                            
                                                                            for (int i_1 = 0; i_1 < 5; i_1++)
                                                                            {
                                                                                temp_ring1.push_back(ring5[first_ring ][i_1]);
                                                                                temp_ring2.push_back(ring5[second_ring][i_1]);
                                                                                temp_ring3.push_back(ring5[third_ring ][i_1]);
                                                                                temp_ring4.push_back(ring5[fourth_ring][i_1]);
                                                                                temp_ring5.push_back(ring5[fifth_ring ][i_1]);
                                                                                temp_ring6.push_back(ring5[sixth_ring ][i_1]);
                                                                                
                                                                            }
                                                                            
                                                                            repeated_atoms.clear();
                                                                            set_intersection(temp_ring1.begin(),temp_ring1.end(),temp_ring2.begin(),temp_ring2.end(),back_inserter(repeated_atoms));    //check ring1 and ring2 for overlaps
                                                                            
                                                                            if (repeated_atoms.size() <= 2)
                                                                            {
                                                                                repeated_atoms.clear();
                                                                                set_intersection(temp_ring2.begin(),temp_ring2.end(),temp_ring3.begin(),temp_ring3.end(),back_inserter(repeated_atoms));
                                                                                if (repeated_atoms.size() <= 2)
                                                                                {
                                                                                    repeated_atoms.clear();
                                                                                    set_intersection(temp_ring3.begin(),temp_ring3.end(),temp_ring4.begin(),temp_ring4.end(),back_inserter(repeated_atoms));
                                                                                    if (repeated_atoms.size() <= 2)
                                                                                    {
                                                                                        repeated_atoms.clear();
                                                                                        set_intersection(temp_ring4.begin(),temp_ring4.end(),temp_ring5.begin(),temp_ring5.end(),back_inserter(repeated_atoms));
                                                                                        if (repeated_atoms.size() <= 2)
                                                                                        {
                                                                                            repeated_atoms.clear();
                                                                                            set_intersection(temp_ring5.begin(),temp_ring5.end(),temp_ring6.begin(),temp_ring6.end(),back_inserter(repeated_atoms));
                                                                                            if (repeated_atoms.size() <= 2)
                                                                                            {
                                                                                                repeated_atoms.clear();
                                                                                                set_intersection(temp_ring6.begin(),temp_ring6.end(),temp_ring1.begin(),temp_ring1.end(),back_inserter(repeated_atoms));
                                                                                                
                                                                                                temp_vec5[0] = first_ring ;
                                                                                                temp_vec5[1] = second_ring;
                                                                                                temp_vec5[2] = third_ring ;
                                                                                                temp_vec5[3] = fourth_ring;
                                                                                                temp_vec5[4] = fifth_ring ;
                                                                                                temp_vec5[5] = sixth_ring ;
                                                                                                
                                                                                                cup512.push_back(temp_vec5);
                                                                                                
                                                                                                temp_ring1.clear();
                                                                                                temp_ring2.clear();
                                                                                                temp_ring3.clear();
                                                                                                temp_ring4.clear();
                                                                                                temp_ring5.clear();
                                                                                                temp_ring6.clear();
                                                                                                
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }       //end of finding half of 5(12) cages
    
    
}

// New Function--------------------------------------------------------------------------------------

// New Function--------------------------------------------------------------------------------------

void cup_62512_Finder(vector<vector<int>>& ring6, int count_ring6, vector<vector<int>>& My_neigh_ring6_ring5, vector<unsigned long int> N_ring6_ring5_neigh, vector<vector<int>>& ring5, vector<unsigned long int> N_ring5_neigh, vector<vector<int>>& My_neigh_ring5, vector<vector<int>>& cup62512){
    
    
    // finding half of 6(2)5(12) cages
    
    
    int first_ring, second_ring, third_ring, fourth_ring, fifth_ring, sixth_ring, seventh_ring;
    vector<int> temp_ring1, temp_ring2, temp_ring3, temp_ring4, temp_ring5, temp_ring6, temp_ring7;
    vector<int> repeated_atoms, temp_vec6={0,0,0,0,0,0,0};
    
    for (int index1 = 0 ; index1 < count_ring6 ; index1++)
    {
        first_ring = index1;
        if (N_ring6_ring5_neigh[first_ring] > 5)
        {
            for (int index2 = 0 ; index2 < N_ring6_ring5_neigh[first_ring] ; index2++)
            {
                second_ring = My_neigh_ring6_ring5[first_ring][index2];
                if (N_ring5_neigh[second_ring] > 3)
                {
                    if(second_ring != first_ring)
                    {
                        for (int index3 = 0 ; index3 < N_ring5_neigh[second_ring] ; index3++)
                        {
                            if (My_neigh_ring5[second_ring][index3] != first_ring &&
                                My_neigh_ring5[second_ring][index3] != second_ring
                                )
                            {
                                if(compare(ring6[first_ring], ring5[My_neigh_ring5[second_ring][index3]], 6, 5) &&
                                   compare_adjacant(ring5[second_ring], ring5[My_neigh_ring5[second_ring][index3]], 5, 5, ring6[first_ring]) )
                                {
                                    third_ring = My_neigh_ring5[second_ring][index3];
                                    if (N_ring5_neigh[third_ring] > 3)
                                    {
                                        for (int index4 = 0 ; index4 < N_ring5_neigh[third_ring]; index4++)
                                        {
                                            if (My_neigh_ring5[third_ring][index4] != first_ring &&
                                                My_neigh_ring5[third_ring][index4] != second_ring &&
                                                My_neigh_ring5[third_ring][index4] != third_ring &&
                                                compare(ring6[first_ring], ring5[My_neigh_ring5[third_ring][index4]], 6, 5) &&
                                                compare_adjacant(ring5[third_ring], ring5[My_neigh_ring5[third_ring][index4]], 5, 5, ring6[first_ring])
                                                )
                                            {
                                                fourth_ring = My_neigh_ring5[third_ring][index4];
                                                if (N_ring5_neigh[fourth_ring] > 3)
                                                {
                                                    for (int index5 = 0 ; index5 < N_ring5_neigh[fourth_ring]; index5++)
                                                    {
                                                        if (My_neigh_ring5[fourth_ring][index5] != first_ring &&
                                                            My_neigh_ring5[fourth_ring][index5] != second_ring &&
                                                            My_neigh_ring5[fourth_ring][index5] != third_ring  &&
                                                            My_neigh_ring5[fourth_ring][index5] != fourth_ring &&
                                                            compare(ring6[first_ring], ring5[My_neigh_ring5[fourth_ring][index5]], 6, 5) &&
                                                            compare_adjacant(ring5[fourth_ring], ring5[My_neigh_ring5[fourth_ring][index5]], 5, 5, ring6[first_ring] )
                                                            )
                                                        {
                                                            fifth_ring = My_neigh_ring5[fourth_ring][index5];
                                                            if (N_ring5_neigh[fifth_ring] != 0 && N_ring5_neigh[fifth_ring] > 3)
                                                            {
                                                                for (int index6 = 0 ; index6 < N_ring5_neigh[fifth_ring]; index6++)
                                                                {
                                                                    if (My_neigh_ring5[fifth_ring][index6] != first_ring  &&
                                                                        My_neigh_ring5[fifth_ring][index6] != second_ring &&
                                                                        My_neigh_ring5[fifth_ring][index6] != third_ring  &&
                                                                        My_neigh_ring5[fifth_ring][index6] != fourth_ring &&
                                                                        My_neigh_ring5[fifth_ring][index6] != fifth_ring  &&
                                                                        compare(ring6[first_ring], ring5[My_neigh_ring5[fifth_ring][index6]], 6, 5) &&
                                                                        compare_adjacant(ring5[fifth_ring], ring5[My_neigh_ring5[fifth_ring][index6]], 5, 5, ring6[first_ring])
                                                                        )
                                                                    {
                                                                        sixth_ring = My_neigh_ring5[fifth_ring][index6];
                                                                        
                                                                        if (N_ring5_neigh[sixth_ring] != 0 && N_ring5_neigh[sixth_ring] > 3)
                                                                        {
                                                                            for (int index7 = 0 ; index7 < N_ring5_neigh[sixth_ring]; index7++)
                                                                            {
                                                                                if (My_neigh_ring5[sixth_ring][index7] != first_ring  &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != second_ring &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != third_ring  &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != fourth_ring &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != fifth_ring  &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != sixth_ring  &&
                                                                                    compare(ring5[second_ring], ring5[My_neigh_ring5[sixth_ring][index7]], 5, 5) &&
                                                                                    compare_adjacant(ring5[sixth_ring], ring5[My_neigh_ring5[sixth_ring][index7]], 5, 5, ring6[first_ring])
                                                                                    )
                                                                                {
                                                                                    seventh_ring = My_neigh_ring5[sixth_ring][index7];
 
                                                                                    for (int index8 = 0 ; index8 < N_ring6_ring5_neigh[first_ring]; index8++)
                                                                                    {
                                                                                        if(My_neigh_ring6_ring5[first_ring][index8] == seventh_ring)
                                                                                        {
                                                                                            if (compare(ring6[first_ring], ring5[seventh_ring], 6, 5))
                                                                                            {
                                                                                                
                                                                                                temp_ring1.clear();
                                                                                                temp_ring2.clear();
                                                                                                temp_ring3.clear();
                                                                                                temp_ring4.clear();
                                                                                                temp_ring5.clear();
                                                                                                temp_ring6.clear();
                                                                                                temp_ring7.clear();
                                                                                                
                                                                                                //put all the rings in temp_ring vectors.
                                                                                                for (int i_1 = 0; i_1 < 6; i_1++)
                                                                                                {
                                                                                                    temp_ring1.push_back(ring6[first_ring ][i_1]);
                                                                                                }
                                                                                                for (int i_1 = 0; i_1 < 5; i_1++)
                                                                                                {
                                                                                                    temp_ring2.push_back(ring5[second_ring][i_1]);
                                                                                                    temp_ring3.push_back(ring5[third_ring ][i_1]);
                                                                                                    temp_ring4.push_back(ring5[fourth_ring][i_1]);
                                                                                                    temp_ring5.push_back(ring5[fifth_ring ][i_1]);
                                                                                                    temp_ring6.push_back(ring5[sixth_ring ][i_1]);
                                                                                                    temp_ring7.push_back(ring5[seventh_ring][i_1]);
                                                                                                }
                                                                                                
                                                                                                repeated_atoms.clear();
                                                                                                set_intersection(temp_ring1.begin(),temp_ring1.end(),temp_ring2.begin(),temp_ring2.end(),back_inserter(repeated_atoms));    //check ring1 and ring2 for overlaps
                                                                                                
                                                                                                if (repeated_atoms.size() <= 2)
                                                                                                {
                                                                                                    repeated_atoms.clear();
                                                                                                    set_intersection(temp_ring2.begin(),temp_ring2.end(),temp_ring3.begin(),temp_ring3.end(),back_inserter(repeated_atoms));
                                                                                                    if (repeated_atoms.size() <= 2)
                                                                                                    {
                                                                                                        repeated_atoms.clear();
                                                                                                        set_intersection(temp_ring3.begin(),temp_ring3.end(),temp_ring4.begin(),temp_ring4.end(),back_inserter(repeated_atoms));
                                                                                                        if (repeated_atoms.size() <= 2)
                                                                                                        {
                                                                                                            
                                                                                                            repeated_atoms.clear();
                                                                                                            set_intersection(temp_ring4.begin(),temp_ring4.end(),temp_ring5.begin(),temp_ring5.end(),back_inserter(repeated_atoms));
                                                                                                            if (repeated_atoms.size() <= 2)
                                                                                                            {
                                                                                                                repeated_atoms.clear();
                                                                                                                set_intersection(temp_ring5.begin(),temp_ring5.end(),temp_ring6.begin(),temp_ring6.end(),back_inserter(repeated_atoms));
                                                                                                                if (repeated_atoms.size() <= 2)
                                                                                                                {
                                                                                                                    repeated_atoms.clear();
                                                                                                                    set_intersection(temp_ring6.begin(),temp_ring6.end(),temp_ring7.begin(),temp_ring7.end(),back_inserter(repeated_atoms));
                                                                                                                    
                                                                                                                    if (repeated_atoms.size() <= 2)
                                                                                                                    {
                                                                                                                        repeated_atoms.clear();
                                                                                                                        set_intersection(temp_ring7.begin(),temp_ring7.end(),temp_ring1.begin(),temp_ring1.end(),back_inserter(repeated_atoms));
                                                                                                                        if (repeated_atoms.size() <= 2)
                                                                                                                        {
                                                                                                                            
                                                                                                                            
                                                                                                                            temp_vec6[0] = first_ring ;
                                                                                                                            temp_vec6[1] = second_ring;
                                                                                                                            temp_vec6[2] = third_ring ;
                                                                                                                            temp_vec6[3] = fourth_ring;
                                                                                                                            temp_vec6[4] = fifth_ring ;
                                                                                                                            temp_vec6[5] = sixth_ring ;
                                                                                                                            temp_vec6[6] = seventh_ring ;
                                                                                                                            
                                                                                                                            
                                                                                            cup62512.push_back(temp_vec6);
                                                                                                                            
                                                                                                                            temp_ring1.clear();
                                                                                                                            temp_ring2.clear();
                                                                                                                            temp_ring3.clear();
                                                                                                                            temp_ring4.clear();
                                                                                                                            temp_ring5.clear();
                                                                                                                            temp_ring6.clear();
                                                                                                                            temp_ring7.clear();
                                                                                                                        }
                                                                                                                    }
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }           //end of finding half of 6(2)5(12) cages---------------------------------------------------------------------------
    
}

// New Function--------------------------------------------------------------------------------------

// New Function--------------------------------------------------------------------------------------

//Find Cages of 5(12) from cups of 5(6).

int cage_Finder(vector<vector<int>> cups, unsigned long count_cups, vector<vector<int>>& My_neigh_ring5, vector<vector<int>>& cage, vector<vector<int>>& cage_rings, string time )
{
    
    int count = 0;      //Counter for number of cages.
    vector<int> current_cup;
    int current_ring = 0 ;
    int N = 0;      //counter for number of common rings between side rings of current cup and side rings of rest of the cups.
    int M = 0;      //counter for number of side rings of current cup which have 2 common rings with side rings of rest of the cups.
    unsigned long int cup_size=0;
    if( cups.size() != 0 ){cup_size = cups[0].size();}
    vector<int> temp_vec, temp_vec_rings;
    
    
    for (int i = 0  ; i < count_cups ; i++)     //Loop over all the cups.
    {
        for (int j = 0 ; j < cups[i].size() ; j++ ) {current_cup.push_back(cups[i][j]);}        //Write current cup into a new variable.
        for (int k = i+1 ; k < count_cups ; k++)        //Loop over the rest of the cups to compare with current cup.
        {
            for(int l = 1 ; l < current_cup.size() ; l++)                           //Loop over side rings of current cup.
                                                                                    //"l" is the counter for side rings of current cup. It starts from 1
                                                                                    //because first ring in each cup is the base ring.
            {
                current_ring = current_cup[l];
                for(int m = 1 ; m < current_cup.size() ; m++)                          //Loop over side rings of the rest of the cups.
                {
                    for(int p = 0 ; p < My_neigh_ring5[current_ring].size() ; p++)      //Loop over neighbour rings of current ring.
                                                                                        //"p" is counter for all the neighbour rings of
                                                                                        //side rings of current cup.
                    {
                        if( My_neigh_ring5[current_ring][p] == cups[k][m] )
                        {N++;}
                    }
                }
                if(N==2){ M++; N = 0;}      //If current side ring has 2 neighbors in other cups, move on
                                            //to next side ring of the current cup.
                else {N = 0 ; break;}       //If current side ring has less than 2 neighbors, move on to the next cup.
            }
            if(M == cup_size -1)
            {
                
                for(int h = 0 ; h < cup_size  ; h++){temp_vec_rings.push_back(cups[i][h]);}
                for(int h = 0 ; h < cup_size  ; h++){temp_vec_rings.push_back(cups[k][h]);}
                cage_rings.push_back(temp_vec_rings);
                temp_vec_rings.clear();
                
                temp_vec.push_back(i);
                temp_vec.push_back(k);
                cage.push_back(temp_vec);
                temp_vec.clear();
                count++;
            }
            M = 0;
        }
        current_cup.clear();
    }
    return count;
    
}

// New Function--------------------------------------------------------------------------------------
// New Function--------------------------------------------------------------------------------------

void print_vmd_cage_frings_test(vector<vector<int>> cups, vector<vector<int>> cages, int cage_count, vector<vector<int>> cage_rings, vector<vector<int>> ring5, vector<vector<int>> ring6, vector<vector<double>> atom_Pos, string time, string rawFilename , string box_size_xyz, vector<vector<double>> solutes, size_t & meth_counter, string solute1, int topSolute, string solute2, int count_solute2, int frameCounter)

{
    
    unsigned long int N_rings = 0;      //Number of rings in a cage. cage 512 has 12 and cage 62512 has 14.
    if(cage_rings.size() != 0){N_rings = cage_rings[0].size();}
    
    unsigned long int N_cage = cage_count;        //This is the number of cages after removing the duplicates. This should be the correct number to use, however, using it will reduce the number of cages printed in gro file output. Something is not right here!
    int atoms=0;
    int p;       //"i" is cup number. "p" is ring number.
    int ring_size;
    int atom_counter=1;        //"atom_counter" is a counter for atom number in gro file. It should reset to one for each frame.
    int sol_counter=1;         //counter for number of SOL in each frame.
    meth_counter=0;      //counter for number of Methanes/Methanol inside cages.
    unsigned long int tot_atoms=0;            // total number of atoms in gro file.
    string out6, out5, out;
    set<int> atoms_set;
    //Different name for different cage types are written out based on the cage type recognized from N_cup_col.
    
    long double cage_sum_x=0, cage_sum_y=0, cage_sum_z=0;
    long double cage_com_x=0, cage_com_y=0, cage_com_z=0;
    int oxygen_counter=0;
    vector<vector<long double>> all_cages_com;
    vector<long double> current_cage_com;
    vector<long double> temp_cage_com;
    double dx, dy, dz, boxX, boxY, boxZ, dist;
    
    out5 = rawFilename + "_cage512-temp.gro";
    out6 = rawFilename + "_cage62512-temp.gro";
    
    istringstream streamB(box_size_xyz);
    streamB >> boxX >> boxY >> boxZ ;
    
    if(N_rings == 14){
        tot_atoms = N_cage * 24 * 4 ;    //each 62512 cage has 24 water molecule. TIP4P has 4 atoms per molecule.
                                         //This number is not neccesarily correct, e.x., if cages share a ring with each other(That's why
                                         //the re-writting of output file corrects it at the end).
        out = out6;
    }
    else if (N_rings == 12){
        tot_atoms = N_cage * 20 * 4 ;    //each 512 cage has 20 water molecule. TIP4P has 4 atoms per molecule.
                                         //This number is not neccesarily correct, e.x.,if cages share a ring with each other(That's why
                                         //the re-writting of output file corrects it at the end).
        out = out5;
    }
    
    std::ofstream outFile (out, std::ofstream::app);     //Open file and append to it.
    outFile.setf(ios::fixed, ios::floatfield);           //Make output numbers look neat
    outFile << "Generated by GROMACS : t= " << time << "\n";      //First line of gro file, comment and frame time.
    outFile << tot_atoms << "\n";
    
    for(int l = 0 ; l < N_cage ; l++)       //Loop over all the cages.
    {
        
        current_cage_com.clear();
        for(int m = 0 ; m < N_rings; m++)         //Loop over 6 or 7 rings of each cage.
        {
            p = cage_rings[l][m];
            {
                if( (N_rings == 14 && m == 0) || (N_rings == 14 && m == 7) )  //This part is for base ring of 6(2)5(12) cages which is a 6-ring.
                {
                    ring_size = 6;
                    for (int q = 0 ; q < ring_size ; q++)       //Loop over all atoms of a ring.
                    {
                        atoms = ring6[p][q];
                        if(!atoms_set.count(atoms))
                        {                               //If atoms is not present in atoms_set, count will be zero, condition is one and
                                                        //the rest of this block will be executed.
                                                        //If atoms already exists in the set, it will not print out.
                            atoms_set.insert(atoms);
                            
                            outFile << setprecision(3);
                            outFile << setw(5) << sol_counter << "SOL" << setw(7) << "OW" << setw(5) << atom_counter;
                            outFile << setw(8) << atom_Pos[atoms][0];
                            outFile << setw(8) << atom_Pos[atoms][1];
                            outFile << setw(8) << atom_Pos[atoms][2];
                            outFile << "\n";
                            atom_counter++;
                            
                            outFile << setw(5) << sol_counter << "SOL" << setw(7) << "HW1" << setw(5) << atom_counter;
                            outFile << setw(8) << atom_Pos[atoms+1][0];
                            outFile << setw(8) << atom_Pos[atoms+1][1];
                            outFile << setw(8) << atom_Pos[atoms+1][2];
                            outFile << "\n";
                            atom_counter++;
                            
                            outFile << setw(5) << sol_counter << "SOL" << setw(7) << "HW2" << setw(5) << atom_counter;
                            outFile << setw(8) << atom_Pos[atoms+2][0];
                            outFile << setw(8) << atom_Pos[atoms+2][1];
                            outFile << setw(8) << atom_Pos[atoms+2][2];
                            outFile << "\n";
                            atom_counter++;
                            
                            outFile << setw(5) << sol_counter << "SOL" << setw(7) << "MW" << setw(5) << atom_counter;
                            outFile << setw(8) << atom_Pos[atoms+3][0];
                            outFile << setw(8) << atom_Pos[atoms+3][1];
                            outFile << setw(8) << atom_Pos[atoms+3][2];
                            outFile << "\n";
                            atom_counter++;
                            
                            sol_counter++;
                            
                        }
                    }
                }
                else            //This part is for all the rings of 5(12) cages which are all 5-rings.
                {
                    ring_size = 5;
                    for (int q = 0 ; q < ring_size ; q++)       //Loop over all atoms of a ring.
                    {
                        atoms = ring5[p][q];
                        if(!atoms_set.count(atoms))
                        {                               //If atoms is not present in atoms_set, count will be zero, condition is one and
                                                        //the rest of this block will be executed.
                                                        //If atoms already exists in the set, it will not print out.
                            
                            atoms_set.insert(atoms);
                            
                            outFile << setprecision(3);
                            outFile << setw(5) << sol_counter << "SOL" << setw(7) << "OW" << setw(5) << atom_counter;
                            outFile << setw(8) << atom_Pos[atoms][0];
                            outFile << setw(8) << atom_Pos[atoms][1];
                            outFile << setw(8) << atom_Pos[atoms][2];
                            outFile << "\n";
                            atom_counter++;
                            
                            outFile << setw(5) << sol_counter << "SOL" << setw(7) << "HW1" << setw(5) << atom_counter;
                            outFile << setw(8) << atom_Pos[atoms+1][0];
                            outFile << setw(8) << atom_Pos[atoms+1][1];
                            outFile << setw(8) << atom_Pos[atoms+1][2];
                            outFile << "\n";
                            atom_counter++;
                            
                            outFile << setw(5) << sol_counter << "SOL" << setw(7) << "HW2" << setw(5) << atom_counter;
                            outFile << setw(8) << atom_Pos[atoms+2][0];
                            outFile << setw(8) << atom_Pos[atoms+2][1];
                            outFile << setw(8) << atom_Pos[atoms+2][2];
                            outFile << "\n";
                            atom_counter++;
                            
                            outFile << setw(5) << sol_counter << "SOL" << setw(7) << "MW" << setw(5) << atom_counter;
                            outFile << setw(8) << atom_Pos[atoms+3][0];
                            outFile << setw(8) << atom_Pos[atoms+3][1];
                            outFile << setw(8) << atom_Pos[atoms+3][2];
                            outFile << "\n";
                            atom_counter++;
                            
                            sol_counter++;
                            
                        }
                    }
                }
            }
        }
    }
    
    
    
    //Find methanes closeest to Center of Mass of each cage and write them to -temp.gro file output.
    
    //---------------------------------------------------------------------------------------------------------
    
    oxygen_counter=1;
    cage_sum_x=0;
    cage_sum_y=0;
    cage_sum_z=0;
    temp_cage_com.clear();
    current_cage_com.clear();
    
    for(int l = 0 ; l < N_cage ; l++)       //Loop over all the cages.
    {
        current_cage_com.clear();
        for(int m = 0 ; m < N_rings; m++)         //Loop over all rings of each cage.(12 or 14)
        {
            
            p = cage_rings[l][m];
            
            {
                if( (N_rings == 14 && m == 0) || (N_rings == 14 && m == 7) ) //This part is for base ring of 6(2)5(12) cages which is a 6-ring.
                {
                    ring_size = 6;
                    for (int q = 0 ; q < ring_size ; q++)       //Loop over all atoms of a ring.
                    {
                        atoms = ring6[p][q];
                        
                        oxygen_counter++;
                        cage_sum_x += atom_Pos[atoms][0];       //Find sum of x coordinate of all "O" atoms, for COM calcultaion.
                        cage_sum_y += atom_Pos[atoms][1];       //find sum of y coordinate of all "O" atoms, for COM calcultaion.
                        cage_sum_z += atom_Pos[atoms][2];       //Find sum of z coordinate of all "O" atoms, for COM calcultaion.
                        
                    }
                }
                else            //This part is for all the rings of 5(12) cages which are all 5-rings.
                {
                    ring_size = 5;
                    for (int q = 0 ; q < ring_size ; q++)       //Loop over all atoms of a ring.
                    {
                        atoms = ring5[p][q];
                        
                        
                        oxygen_counter++;
                        cage_sum_x += atom_Pos[atoms][0];       //Find sum of x coordinate of all "O" atoms, for COM calcultaion.
                        cage_sum_y += atom_Pos[atoms][1];       //find sum of y coordinate of all "O" atoms, for COM calcultaion.
                        cage_sum_z += atom_Pos[atoms][2];       //Find sum of z coordinate of all "O" atoms, for COM calcultaion.
                    }
                    
                }
                
            }
            
        }
        
        cage_com_x = cage_sum_x / oxygen_counter;
        cage_com_y = cage_sum_y / oxygen_counter;
        cage_com_z = cage_sum_z / oxygen_counter;
        
        current_cage_com.push_back(cage_com_x);
        current_cage_com.push_back(cage_com_y);
        current_cage_com.push_back(cage_com_z);
        
        
        if(temp_cage_com != current_cage_com)
        {
            temp_cage_com.clear();
            all_cages_com.push_back(current_cage_com);
            temp_cage_com= current_cage_com;
            
            current_cage_com.clear();
        }
        cage_sum_x=0;
        cage_sum_y=0;
        cage_sum_z=0;
        oxygen_counter=0;
        
    }
    
    //---------------------------------------------------------------------------------------------------------
    //***********************************************************************************************************
    //This part removes the duplicates from all_cages_com variable.
    
    unsigned long int old_size = all_cages_com.size();              //Number of cups before removing the duplicates.
    int count=0;                                                   //Number of com.
    vector<double> current;                                        //Cup which is being read and compared at each loop.
    map<vector<double>,int> mymap;                                 //A map which counts the frequencty of each com( [com,frequency] )
    int prev = 0 ;
    
    for (int i = 0 ; i < old_size ; i++)      //Loop over all the com.
    {
        // number "3" in next loop over "j" is for x,y,z coordinates of com.
        for(int j = 0 ; j < 3 ; j++){current.push_back(all_cages_com[i][j]);}      //Loop to put the read COM into "current".
                                                                                   //sort(current.begin(), current.end());
        if(!mymap[current])                  //For each read com, check if it exists in the map. If it does NOT, write it in the map.
            
        {
            mymap[current] = 1;
            count++;
            for (int k = 0 ; k < 3 ; k++)
            {
                all_cages_com[prev][k] = all_cages_com[i][k];   //If the com does not exist in map, write it to first line of original com variable.
            }
            prev++;
        }
        else {mymap[current]++;}        //If the com already exists in the map, jump the writing step and add one to frequencty of that com.
        current.clear();
    }
    //*_*_*_*_*_*_*_*_*_*_*_*_*_
    /*Write all the solutes ( solute1 and solute2 ) at the end of the gro file. This part is added in v1.16. */
    size_t solutes_size = solutes.size();
    int s1 = 7, s2 = 7;
    int meth_counter2=0;
    for(int i = 1 ; i <= solute1.size(); ++i){s1 = 10 - i;}
    for(int i = 1 ; i <= solute2.size(); ++i){s2 = 10 - i;}
    
    
    for (int i = 0 ; i < solutes_size ; i++)
    {
        if(i < topSolute )
        {
            
            outFile << setprecision(3);
            outFile << setw(5) << i+1 << solute1  << setw(s1) << "CB" << setw(5) << ++meth_counter2;
            outFile << setw(8) << solutes[i][0];
            outFile << setw(8) << solutes[i][1];
            outFile << setw(8) << solutes[i][2];
            outFile << "\n";
            atom_counter++;
        }
        if(i >= topSolute)
        {
            outFile << setprecision(3);
            outFile << setw(5) << i+1 << solute2  << setw(s2) << "CB" << setw(5) << ++meth_counter2;
            outFile << setw(8) << solutes[i][0];
            outFile << setw(8) << solutes[i][1];
            outFile << setw(8) << solutes[i][2];
            outFile << "\n";
            atom_counter++;
        }
    }
    //*_*_*_*_*_*_*_*_*_*_*_*_*_
    
    //Write Methane molecules inside cages to the -temp.gro file.
    //Find methanes closeest to Center of Mass of each cage and write them to -temp.gro file output.
    vector<int> temp_meth_number_holder;
    temp_meth_number_holder.clear();
    
    for ( int i = 0; i < count ; i++)        //Calculate distance between com of cages and each solute molecule.
    {
        for ( int j = 0 ; j < solutes.size() ; j++)
        {
            dx = solutes[j][0] - all_cages_com[i][0];
            dy = solutes[j][1] - all_cages_com[i][1];
            dz = solutes[j][2] - all_cages_com[i][2];
            
            if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
            if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
            if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
            
            dist = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (dist <= 0.2)        //If distance between COM of cage and any Methane molecule is less than cut-off, append it to gro file.
            {
                
                if(find(temp_meth_number_holder.begin(), temp_meth_number_holder.end(), j) == temp_meth_number_holder.end()) //If the value of "j" (methane index) does not exist in the temp_meth_number_holder, then add it to the vector and write its coordinates to gro file. Otherwise, do not write it in gro output file. The find(v.begin(),v.end(),"value") gives the iterator number of found value in the vector. If the value is not in the vector, iterator will be equal to v.end().
                {
                    temp_meth_number_holder.push_back(j);
                    
                    outFile << setprecision(3);
                    outFile << setw(5) << j+1 << "CH3" << setw(7) << "CB" << setw(5) << ++meth_counter;
                    outFile << setw(8) << solutes[j][0];
                    outFile << setw(8) << solutes[j][1];
                    outFile << setw(8) << solutes[j][2];
                    outFile << "\n";
                    atom_counter++;
                }
            }
            
        }
    }
    
    outFile.close();
    
    
    
    //***********************************************************************************************************
    //This part re-writes the output gro file to correct the number of atoms in the gro file. It will do that, only if N_cage > 0 .
    if(N_cage>0)
    {
        string line;        // Reads each line.
        string out_name, raw_filename;
        int line_number = 0;
        ofstream fout;
        ifstream fin;
        size_t found_ext = out.find_last_of("-");
        string frame_number = std::to_string(frameCounter);
        raw_filename = out.substr(0,found_ext);
        out_name = raw_filename + "-" + frame_number + ".gro";
        
        
        fout.open(out_name,ofstream::app);
        fin.open(out);
        //Error Check
        if(fin.fail()){
            cerr << "Error Reading temp file for re-writing with correct atom numbers in gro file" << "\n";
            exit (1);}
        while ( !fin.eof() )
        {
            getline(fin,line);
            line_number++;
            
            if(line_number==2)fout << atom_counter-1 << "\n"  ;
            
            else if(!line.empty() )fout << line << "\n"  ;
        }
        
        fout << box_size_xyz << "\n" ;
        
        fin.close();
        fout.close();
        
        atom_counter = 0;        //reset atom number for each frame.(to be used in gro file output)
        sol_counter = 0;         //reset SOL number for each frame.(to be used in gro file output)
        
        remove(out.c_str());
        
    }
    }
// New Function--------------------------------------------------------------------------------------
// New Function--------------------------------------------------------------------------------------


double calc_F4(int count_solvent, int count_solute, vector<vector<int>>& My_neigh, vector<vector<double>>& atom_Pos, double& boxX, double& boxY, double& boxZ, vector<int>& Nneigh, int Natoms, int topSolute, string time, double HBOND_DIST ){
    
    vector<int> currentPair;        //Current pair of axygen atoms to be considered.
    vector<double> temp_vec;        //This is used to sort the 4 H pair distances.
    vector<double> A, B, C, D, AB, BC, CD, ABC, BCD;
    map<vector<int>,int> pairMap;   //A map to exclude double calculation of torsion angle between atom 'i-j' or 'j-i'.
    double dx=0, dy=0, dz=0 ;       //Distances for each pair.
    double dist_OA_H1B, dist_OA_H2B, dist_OB_H1A, dist_OB_H2A;      //Total distance for each pair.
    double phi=0, phi_avg=0, sum=0, magnitude_ABC=0, magnitude_BCD=0;          //Parameters used in cross product calculation.
    int dih_A_index = 0, dih_B_index = 0, dih_C_index = 0, dih_D_index = 0;         //Four atoms of the dihedral.

    for ( int i = topSolute + 1 ; i < topSolute + count_solvent  ; i+=4 )   //Loop over oxygen atoms, starting after solutes.
    {
        if (My_neigh[i].size() > 0 )                        //If current oxygen has any neighbors, do the following
        {
            for (int j = 1 ; j < My_neigh[i].size() ; j++)      //Loop over all the neighbors of oxygen 'i'.
            {

                currentPair.push_back(i);
                currentPair.push_back(My_neigh[i][j]);
                
                sort(currentPair.begin(), currentPair.end());   //Sort items in map.
              
                if (!pairMap[currentPair])              //If the currentPair does not exist in the map, do the following(all calculations):
                {
                    pairMap[currentPair] = 1;
                    
                    //Find the outermost hydrogen atoms in the two current water molecules.
                    
                    //Find distance between each H and other O. Ther order is OA_H*B and OB_H*A.
                    
                    //OA_H1B,  OA=i, H1B=My_neigh[i][j]+1
                    dx = atom_Pos[i][0] - atom_Pos[ My_neigh[i][j]+1 ][0];
                    dy = atom_Pos[i][1] - atom_Pos[ My_neigh[i][j]+1 ][1];
                    dz = atom_Pos[i][2] - atom_Pos[ My_neigh[i][j]+1 ][2];
                    if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                    if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                    if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                    dist_OA_H1B = sqrt(dx*dx + dy*dy + dz*dz);
                    
                    //OA_H2B, OA=i, H2B=My_neigh[i][j]+2
                    dx = atom_Pos[i][0] - atom_Pos[ My_neigh[i][j]+2 ][0];
                    dy = atom_Pos[i][1] - atom_Pos[ My_neigh[i][j]+2 ][1];
                    dz = atom_Pos[i][2] - atom_Pos[ My_neigh[i][j]+2 ][2];
                    if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                    if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                    if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                    dist_OA_H2B = sqrt(dx*dx + dy*dy + dz*dz);
                    
                    //OB_H1A, OB=My_neigh[i][j], H1A=i+1
                    dx = atom_Pos[i+1][0] - atom_Pos[ My_neigh[i][j] ][0];
                    dy = atom_Pos[i+1][1] - atom_Pos[ My_neigh[i][j] ][1];
                    dz = atom_Pos[i+1][2] - atom_Pos[ My_neigh[i][j] ][2];
                    if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                    if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                    if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                    dist_OB_H1A = sqrt(dx*dx + dy*dy + dz*dz);

                    
                    //OB_H2A, OB=My_neigh[i][j], H2A=i+2
                    dx = atom_Pos[i+2][0] - atom_Pos[ My_neigh[i][j] ][0];
                    dy = atom_Pos[i+2][1] - atom_Pos[ My_neigh[i][j] ][1];
                    dz = atom_Pos[i+2][2] - atom_Pos[ My_neigh[i][j] ][2];
                    if (abs(dx) >= boxX * 0.5){ dx = boxX - abs(dx) ;}
                    if (abs(dy) >= boxY * 0.5){ dy = boxY - abs(dy) ;}
                    if (abs(dz) >= boxZ * 0.5){ dz = boxZ - abs(dz) ;}
                    dist_OB_H2A = sqrt(dx*dx + dy*dy + dz*dz);
                    
                    //case1: dihedral atoms are : H1A-OA-OB-H1B or H1A-OA-OB-H2B
                    if ( min(dist_OA_H1B, min(dist_OA_H2B, min(dist_OB_H1A, dist_OB_H2A))) == dist_OB_H2A ) //If H2A is the H closest to OB, then first 3 atoms of dihedral are : H1A, OA and OB. The last H will be the furthest H from OA.
                    {
                        dih_A_index = i+1;                       //First atom in dihedral is H1A=i+1.
                        dih_B_index = i ;                        //Second atom in dihedral is OA=i.
                        dih_C_index = My_neigh[i][j];            //Third atom in dihedral is OB=My_neigh[i][j].
                        if( max(dist_OA_H1B, max(dist_OA_H2B, max(dist_OB_H1A, dist_OB_H2A))) == dist_OA_H1B)
                            dih_D_index = My_neigh[i][j]+1;         //If dist_OA_H1B is greatest, fourth atom in dihedral is H1B=My_neigh[i][j]+1.
                        else dih_D_index = My_neigh[i][j]+2;        //If dist_OA_H1B is not the greatest, fourth atom is H2B=My_neigh[i][j]+2.
                    }
                    
                    //case2: dihedral atoms are : H2A-OA-OB-H1B or H2A-OA-OB-H2B
                    if ( min(dist_OA_H1B, min(dist_OA_H2B, min(dist_OB_H1A, dist_OB_H2A))) == dist_OB_H1A ) //If H1A is the H closest to OB, then first 3 atoms of dihedral are : H2A, OA and OB. The last H will be the furthest H from OA.
                    {
                        dih_A_index = i+2;                       //First atom in dihedral is H2A=i+2.
                        dih_B_index = i ;                        //Second atom in dihedral is OA=i.
                        dih_C_index = My_neigh[i][j];            //Third atom in dihedral is OB=My_neigh[i][j].
                        if( max(dist_OA_H1B, max(dist_OA_H2B, max(dist_OB_H1A, dist_OB_H2A))) == dist_OA_H1B)
                            dih_D_index = My_neigh[i][j]+1;         //If dist_OA_H1B is greatest, fourth atom in dihedral is H1B=My_neigh[i][j]+1.
                        else dih_D_index = My_neigh[i][j]+2;        //If dist_OA_H1B is not the greatest, fourth atom is H2B=My_neigh[i][j]+2.
                    }
                    
                    //case3: dihedral atoms are : H2A-OA-OB-H2B OR H1A-OA-OB-H2B. This gives the same dihedral angles as case2-b and case1-b, however, it needs to be in another loop. Otherwise, the results will be wrong.
                    if (min(dist_OA_H1B, min(dist_OA_H2B, min(dist_OB_H1A, dist_OB_H2A))) == dist_OA_H1B )
                    {
                        dih_A_index = My_neigh[i][j]+2 ;                        //First atom in dihedral is H2B=My_neigh[i][j]+2.
                        dih_B_index = My_neigh[i][j] ;                          //Second atom in dihedral is OB=My_neigh[i][j].
                        dih_C_index = i ;                                       //Third atom in dihedral is OA=i.
                        if( max(dist_OA_H1B, max(dist_OA_H2B, max(dist_OB_H1A, dist_OB_H2A))) == dist_OB_H1A)
                            dih_D_index = i+1;                                  //If dist_OB_H1A is greatest, fourth atom in dihedral is H1A=i+1.
                        else dih_D_index = i+2;                                 //If dist_OB_H1A is not the greatest, fourth atom is H2A=i+2.
                    }
                    
                    //case4: dihedral atoms are : H2A-OA-OB-H1B OR H1A-OA-OB-H1B. This gives the same dihedral angles as case2-a and case1-a, however, it needs to be in another loop. Otherwise, the results will be wrong.
                    if (min(dist_OA_H1B, min(dist_OA_H2B, min(dist_OB_H1A, dist_OB_H2A))) == dist_OA_H2B )
                    {
                        dih_A_index = My_neigh[i][j]+1 ;                        //First atom in dihedral is H1B=My_neigh[i][j]+1.
                        dih_B_index = My_neigh[i][j] ;                          //Second atom in dihedral is OB=My_neigh[i][j].
                        dih_C_index = i ;                                       //Third atom in dihedral is OA=i.
                        if( max(dist_OA_H1B, max(dist_OA_H2B, max(dist_OB_H1A, dist_OB_H2A))) == dist_OB_H1A)
                            dih_D_index = i+1;                                  //If dist_OB_H1A is greatest, fourth atom in dihedral is H1A=i+1.
                        else dih_D_index = i+2;                                 //If dist_OB_H1A is not the greatest, fourth atom is H2A=i+2.
                    }
                    
                    
                    //Find the angle:(begin)----------------------
                    
                    //Put atom positions in vectors A to C.
                    A.push_back(atom_Pos[ dih_A_index ][0]);
                    A.push_back(atom_Pos[ dih_A_index ][1]);
                    A.push_back(atom_Pos[ dih_A_index ][2]);
                    
                    B.push_back(atom_Pos[ dih_B_index ][0]);
                    B.push_back(atom_Pos[ dih_B_index ][1]);
                    B.push_back(atom_Pos[ dih_B_index ][2]);
                    
                    C.push_back(atom_Pos[ dih_C_index ][0]);
                    C.push_back(atom_Pos[ dih_C_index ][1]);
                    C.push_back(atom_Pos[ dih_C_index ][2]);
                    
                    D.push_back(atom_Pos[ dih_D_index ][0]);
                    D.push_back(atom_Pos[ dih_D_index ][1]);
                    D.push_back(atom_Pos[ dih_D_index ][2]);
         
                    //Create vectors AB, BC, CD.
                    AB.push_back(B[0]-A[0]);
                    AB.push_back(B[1]-A[1]);
                    AB.push_back(B[2]-A[2]);
                    
                    BC.push_back(C[0]-B[0]);
                    BC.push_back(C[1]-B[1]);
                    BC.push_back(C[2]-B[2]);
                    
                    CD.push_back(D[0]-C[0]);
                    CD.push_back(D[1]-C[1]);
                    CD.push_back(D[2]-C[2]);
                    
                    
                    //Calculate cross products: ABxBC=ABC, BCxCD=BCD
                    
                    ABC.push_back(AB[1]*BC[2] - AB[2]*BC[1]);
                    ABC.push_back(AB[2]*BC[0] - AB[0]*BC[2]);
                    ABC.push_back(AB[0]*BC[1] - AB[1]*BC[0]);
                    
                    BCD.push_back(BC[1]*CD[2] - BC[2]*CD[1]);
                    BCD.push_back(BC[2]*CD[0] - BC[0]*CD[2]);
                    BCD.push_back(BC[0]*CD[1] - BC[1]*CD[0]);
                    
                    
                    //Calculate the angle between ABC and BCD from their dot product. call it phi.
                    
                    for (int i = 0 ; i < 3 ; i++)
                    {
                        sum += ABC.at(i) * BCD.at(i);
                    }
                    
                    magnitude_ABC = sqrt(pow(ABC.at(0), 2) + pow(ABC.at(1), 2) + pow(ABC.at(2), 2));
                    magnitude_BCD = sqrt(pow(BCD.at(0), 2) + pow(BCD.at(1), 2) + pow(BCD.at(2), 2));
                    
                    phi = acos( sum / (magnitude_ABC * magnitude_BCD) );                      //Keep the result in radians for cos calculation.

                    phi_avg += cos(3*phi);
                    
                    //Find the angle:(end)----------------------
                    
                    //Clear the vectors for next loop.
                    A.clear();
                    B.clear();
                    C.clear();
                    D.clear();
                    AB.clear();
                    BC.clear();
                    CD.clear();
                    ABC.clear();
                    BCD.clear();
                    sum = 0;
                 }
                else    pairMap[currentPair]++ ;        //If currentPair already exists in the map, increase the frequency.
                currentPair.clear();                    
            }
        }
    }
    
    cout << "phi_avg= " << phi_avg/count_solvent << "\n\n";
    return phi_avg/count_solvent;
}

