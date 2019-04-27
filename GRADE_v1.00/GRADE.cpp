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
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <unistd.h>


#include "MyFunctions.hpp"

using namespace std;


// Main function ----------------------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    double const Version=1.00 ;
    
    //output program information.
    cout << std::fixed ;
    cout << std::setprecision(2);
    cout << "\n\n\t\t** GRADE - VERSION " << Version << " **\n\n" ;
    
    string inputFilename, s1, s2 ;
    string outputFilename, rawFilename;
    int DT=1, FR=1, THETA=45;
    int in_theta=0, in_fr=0, in_dt=0, in_r=0, in_F4=0, in_d1=0, in_d2=0, in_s1=0, in_s2=0;         //This parameter is for whether theta is given as input parameter (1) or taken as default value (0).
    double HBOND_DIST = 0.35;   //Command line input parameter for Hbond_distance cutoff, taken from '-r' flag.
    double delta_p = 0.18, delta_h = 0.26;   //Command line input parameter for delta constraints for pentagon and hexagons.
    string F4 = "no";
    
    //Parsing command line parameters (input taken from "-i" and output taken from "-o")
    std::string arg;
    std::string inName, outName;
    
    
    
    //-------------------------------------------------------------------------------------
    
    int argi=1;
    
    while ( argi < argc )
    {
        const char* args = argv[ argi++ ] ;
        
        switch ( *args )
        {
            case '-':
                if( strcmp(args, "-i") == 0 )
                {
                    inputFilename = argv[argi];
                }
                else if (strcmp(args, "-o") == 0 )
                {
                    outputFilename = argv[argi];
                }
                else if (strcmp(args, "-dt") == 0)
                {
                    DT = atoi(argv[argi]);
                    in_dt = 1;
                }
                else if (strcmp(args, "-fr") == 0)
                {
                    FR = atoi(argv[argi]);
                    in_fr = 1 ;
                }
                else if (strcmp(args, "-theta") == 0)
                {
                    THETA = atoi(argv[argi]);
                    in_theta=1;
                }
                else if (strcmp(args, "-r") == 0)
                {
                    HBOND_DIST = atof(argv[argi]);
                    in_r = 1 ;
                }
                else if (strcmp(args, "-f4") == 0 || strcmp(args, "-F4") == 0)
                {
                    in_F4 = 1;
                    F4 = "YES" ;
                    if(argi < argc )F4 =argv[argi];
                    if(F4 == "no" || F4 == "No" || F4 == "NO") {in_F4=0; F4 = "NO";}
                }
                else if (strcmp(args, "-d1") == 0)
                {
                    delta_p = atof(argv[argi]);
                    in_d1 = 1 ;
                }
                else if (strcmp(args, "-d2") == 0)
                {
                    delta_h = atof(argv[argi]);
                    in_d2 = 1 ;
                }
                
                else if (strcmp(args, "-s1") == 0)
                {
                    s1 = argv[argi];
                    in_s1 = 1 ;
                }
                else if (strcmp(args, "-s2") == 0)
                {
                    s2 = argv[argi];
                    in_s2 = 1 ;
                }
                
                else cout << "Skipped unknown option(s): '" << args << " " <<argv[argi] << "'" << "\n\n";
                break;
                
        }
    }
    if(inputFilename.empty())
    {
        cout << "**No Input Provided!**" << "\n\n";
        print_usage();
        cout << "**No Input Provided!**" << "\n\n";
        return 0 ;
    }
    else{print_usage();}
    
    if (outputFilename.empty())
    {
        outputFilename = inputFilename;
    }
    
    cout << "Command line:\n" ;
    cout << argv[0] << " -i " << inputFilename << " -o " << outputFilename ;
    if (in_dt == 1)  {cout << " -dt " << DT ;}
    if (in_fr == 1)  {cout << " -fr " << FR ;}
    if (in_theta==1) {cout << " -theta " << THETA;}
    if (in_r == 1)   {cout << " -r " << HBOND_DIST ;}
    if (in_F4 == 1)  {cout << " -f4 " << "YES" ;}
    if (in_d1 == 1)  {cout << " -d1 " << delta_p ;}
    if (in_d2 == 1)  {cout << " -d2 " << delta_h ;}
    if (in_s1 == 1)  {cout << " -s1 " << s1 ;}
    if (in_s2 == 1)  {cout << " -s1 " << s2 ;}
    
    cout << "\n\n" ;
    //-------------------------------------------------------------------------------------
    
    vector<vector<int>> My_neigh ;
    int Natoms, count_solute=0 , count_solute2=0;
    vector<int> Nneigh;
    vector<vector<double>> atom_Pos;
    double F4_value=0;                                  //Value of F4 order parameter for each frame.
    vector<double> current_F4_line;
    vector<vector<double>> time_vs_F4;                 //This variable holds Time/frame_counter in first column and value of F4 order parameter in second column.
    
    
    
    string line="NONE", str1, str2, str3;
    int int1;
    double x , y , z ;
    double boxX = 0.0, boxY = 0.0, boxZ = 0.0;
    vector<double> temp_vect;
    vector<int> temp_vect2;
    int lineNumber=0;
    int firstSOL=0;
    int frameCounter=0;
    int topSolute =0;
    size_t methane_512 = 0, methane_62512 = 0, methane_64512 = 0;
    string time;
    Natoms=0;
    string solute1="AAA";
    string solute2="BBB";
    
    
    ofstream outFile;
    size_t found_ext = outputFilename.find_last_of(".");
    rawFilename = outputFilename.substr(0,found_ext);
    outputFilename = rawFilename + ".xvg";
    string temp1, temp2, solute1_atom_finder, solute2_atom_finder;
    temp1 = rawFilename + "_cage62512.gro" ;
    temp2 = rawFilename + "_cage512.gro" ;
    remove(temp1.c_str());
    remove(temp2.c_str());
    
    //Create the header for outputfile.
    outFile.open(outputFilename, ofstream::app);
    outFile << "# ------------------------------------------------------------------------------------------------------ " << endl;
    outFile << "#|(Frame) Time(ps)\t\t|cage\t|filled_cage\t|cage\t|filled_cage\t|cage\t|filled_cage\t|" << endl ;
    outFile << "#|\t\t\t\t|5¹²\t|5¹²\t\t|6²5¹²\t|6²5¹²\t\t|6⁴5¹²\t|6⁴5¹²\t\t|" << endl ;
    outFile << "# ------------------------------------------------------------------------------------------------------ " << endl;
    
    ifstream fileIN;
    fileIN.open(inputFilename);
    
    //Error Check
    if(fileIN.fail()){
        cerr << "Error Reading File" << "\n";
        exit (1);}
    
    ofstream outFile_F4;
    if(in_F4 == 1)                          //If F4 flag option is on, open a file for F4 as a function of time.
    {
        remove("F4.xvg");       //Remove any existing F4.xvg file and create a new one. 
        outFile_F4.open("F4.xvg", ofstream::app);
        outFile_F4 << "# -------------------------------------- \n" ;
        outFile_F4 << "#|Frame\t|F4\t\t|Time(ps)\t|" << endl;
        outFile_F4 << "# -------------------------------------- \n" ;

    }
    
    //Start reading the input file.
    while (!fileIN.eof())
    {
        getline(fileIN, line);
        lineNumber++;
        //Following if statement makes sure the "Natoms" always has the correct number, even if the code skips the first frame due to FR being greater than 1. (Added in v1.18)
        if(lineNumber == 2 || (lineNumber == (2 + 3*frameCounter + frameCounter*Natoms) && !fileIN.eof()))
        {
            istringstream streamA(line);
            
            streamA >> Natoms ;
        }
        size_t found =0;
        size_t found_time=0;
        if( line.find("t=") ) {found_time = line.find("t=");}
        if (lineNumber == 1 || (lineNumber == (1 + 3*frameCounter + frameCounter*Natoms) && !fileIN.eof()) )        //This is to find the first line of gro file.
        {
            frameCounter++;
            if(( (frameCounter-1) % FR) != 0 ) continue;        //If frameCounter-1 is not a multiple of FR, continue to next frame.
                                                                //"-1" is to ignore the frame with t=0.000. This ensures that FR=20 reads frames with t=0,20,40,... .
            if(found_time != string::npos)              //If found_time is not null, do the following.
            {
                time = line.substr(found_time+3);
                outFile << frameCounter << "\t " << time << "\t\t "  ;
            }
            else outFile << frameCounter << "\t\t\t " ;
            
            cout << " frame#: " << frameCounter << ", "  ;
            if(found_time != string::npos) cout << line.substr(found_time) << " ps\n" ;
            else cout << "\n" ;
            
            if(found == string::npos)cout << line << endl;
            
            getline(fileIN, line);      //Read 2nd line (number of atoms)
            lineNumber++;
            
            
            istringstream streamA(line);
            
            streamA >> Natoms ;
            
            vector<int> temp_vec = {0,0,0};
            
            int  count_solvent=0;
            vector<vector<double>> solutes;
            count_solute=0;
            temp_vect={0,0,0};
            atom_Pos.clear();
            atom_Pos.push_back(temp_vect);
            
            
            while (atom_Pos.size() <= Natoms && atom_Pos.size() <= 9999)
            {
                
                temp_vect.clear();
                
                getline(fileIN, line);
                lineNumber++;
                
                
                istringstream streamA(line);
                
                streamA >> str1 >> str2 >> int1 >> x >> y >> z ;
                
                temp_vect.push_back(x);
                temp_vect.push_back(y);
                temp_vect.push_back(z);
                
                
                atom_Pos.push_back(temp_vect);
                
                if(line.find("SOL") != string::npos)                        //If line includes "SOL", then add to number of count_solvent.          Else, add to number of count_solute.
                {
                    count_solvent++;
                    if(count_solvent == 1 ){firstSOL = lineNumber;}
                }
                else
                {
                    if(count_solute == 0)
                    {
                        solute1 = line.substr(5, 7);                        //Get the name of first solute in system.
                        size_t found_space = solute1.find_first_of(" ");
                        solute1 = solute1.substr(0, found_space);
                    }
                    
                    if(line.substr(5,7) != solute1)
                    {
                        
                        solute2 = line.substr(5,7);                         //Get the name of second solute in system.
                        size_t found_space = solute2.find_first_of(" ");
                        solute2 = solute2.substr(0, found_space);
                        
                        count_solute2++;
                        
                    }
                    count_solute++;
                    solutes.push_back(temp_vect);
                }
                
            }
            
            while (atom_Pos.size() <= Natoms && atom_Pos.size() > 9999)
            {
                
                
                temp_vect.clear();
                
                getline(fileIN, line);
                lineNumber++;
                
                istringstream streamA(line);
                
                streamA >> str1 >> str2 >> x >> y >> z ;
                
                temp_vect.push_back(x);
                temp_vect.push_back(y);
                temp_vect.push_back(z);
                
                
                atom_Pos.push_back(temp_vect);
                
                if(line.find("SOL") != string::npos)        //If line includes "SOL"
                {
                    count_solvent++;
                }
                else
                {
                    count_solute++;
                    solutes.push_back(temp_vect);
                    
                    if(line.substr(5,7) != solute1)
                    {
                        solute2 = line.substr(5,7);
                        size_t found_space = solute2.find_first_of(" ");
                        solute2 = solute2.substr(0, found_space);
                        count_solute2++;
                    }
                    
                }
                
            }
            
            getline(fileIN, line);      //get the box size from last line of frame
            //End of reading each frame.
            
            string box_size_xyz;
            box_size_xyz = line;
            lineNumber++;
            
            istringstream streamB(line);
            streamB >> boxX >> boxY >> boxZ ;
            
            if(frameCounter == 1 )
            {
                topSolute = firstSOL - 3;
                if(in_s1 == 0){cout << "solute1: " << solute1 << " " << topSolute << " atoms " ;}
                else {cout << "solute1: " << s1 << " " << topSolute << " atoms " ;}
                //if( (in_s1==0 && in_s2==1) || (in_s1==1 && in_s2==1) || (in_s1==0 && in_s2==0))
                {
                    
                    if(strcmp(solute1.c_str(), solute2.c_str()) != 0)
                    {
                       if(in_s2==0) cout << ", solute2: " << solute2 << " " << count_solute2-topSolute <<" atoms\n";
                        else cout << ", solute2: " << s2 << " " << count_solute2-topSolute <<" atoms\n";
                    }
                    else cout << "\n" ;
                    
                }
                
            }

            //Start of calculations for each frame
            
            calc_Distance(count_solvent, count_solute, My_neigh, atom_Pos, boxX, boxY, boxZ, Nneigh, Natoms, topSolute, time, HBOND_DIST);
            
            vector<vector<int>> ring5, ring6, ring5_temp, ring6_temp;
            vector<vector<int>> My_neigh_ring6, My_neigh_ring5, My_neigh_ring6_ring5;
            
            
            ring_Finder(count_solute, Natoms, Nneigh, My_neigh, ring5_temp, ring6_temp, topSolute ,count_solvent, atom_Pos, boxX, boxY, boxZ, HBOND_DIST,delta_p, delta_h);
            
            
            if( ring5_temp.size() > 0 ) coplanar_Points(ring5_temp, atom_Pos, time, ring5,THETA);        //Find the 5-rings which form a plane and get rid of the rest.
            
            int count_ring5 = remove_duplicates_map_rings(ring5);     //Remove duplicate lines from 5-rings.
            
            cout << "ring[5]: " << count_ring5 << "\n";
            
            if( ring6_temp.size() > 0 )coplanar_Points(ring6_temp, atom_Pos, time, ring6, THETA);
            
            
            int count_ring6 = remove_duplicates_map_rings(ring6);
            
            cout << "ring[6]: " << count_ring6 << "\n";
            
            vector<unsigned long int> N_ring5_neigh;      //Vector of ring-neighbours of all rings.
            
            
            find_shared_edges_ring5(count_ring5, ring5, My_neigh_ring5, N_ring5_neigh);
            
            vector<unsigned long int> N_ring6_neigh;
            
            vector<unsigned long int> N_ring6_ring5_neigh;
            find_shared_edges_ring6_ring5(count_ring6, count_ring5, ring5, ring6, My_neigh_ring6_ring5, N_ring6_ring5_neigh);
            
            
            if (count_ring5 == 0 && count_ring6 == 0 ){cout << "\n**NO RINGS FOUND!**\n\n" ;
                
            }
            
            vector<vector<int>> cup512;
            cup_512_Finder(ring5, count_ring5, N_ring5_neigh, My_neigh_ring5, cup512);
            int count_512_cups = 0;
            if(cup512.size() != 0){count_512_cups = remove_duplicates_map(cup512);}        //Remove duplicate lines from 512 cups.
            
            if (count_512_cups == 0){cout << "\n**NO 5⁶ CUPS FOUND**\n\n" ;}
            else {cout << "# 5⁶ \t cup: " << count_512_cups << "\n";}
            
            vector<vector<int>> cup62512;
            cup_62512_Finder(ring6, count_ring6, My_neigh_ring6_ring5, N_ring6_ring5_neigh, ring5, N_ring5_neigh, My_neigh_ring5, cup62512);
            int count_62512_cups=0;
            if(cup62512.size() != 0){count_62512_cups = remove_duplicates_map(cup62512);}
            
            if (count_62512_cups == 0 ){cout << "\n**NO 6¹5⁶ CUPS FOUND**\n\n" ;}
            else {cout << "# 6¹5⁶ \t cup: " << count_62512_cups << "\n";}
            
            
            vector<vector<int>> cage_512, cage_62512, cage_64512;
            vector<vector<int>> cage_512_rings, cage_62512_rings, cage_64512_rings;
            
            int cage_512_count = 0;
            
            if(count_512_cups > 0)
            {
                
                cage_512_count = cage_Finder(cup512, count_512_cups, My_neigh_ring5, cage_512, cage_512_rings, time);
                
                cage_512_count = remove_duplicates_map(cage_512_rings);
                
                //Write output gro file only if frameCounter is a multiple of DT. if -dt is not provided, every frame will be written.
                if(cage_512_count>0  && (frameCounter % DT == 0))print_vmd_cage64512_frings(cup512, cage_512_count, cage_512_rings, ring5, ring6, atom_Pos, time, rawFilename, box_size_xyz, solutes, methane_512, solute1, topSolute, solute2, count_solute2, frameCounter);
                
            }
            
            cout << "# 5¹²\tcage: " << cage_512_count << "\n";
            
            int cage_62512_count = 0;
            if(count_62512_cups > 0)
            {
                cage_62512_count = cage_Finder(cup62512, count_62512_cups, My_neigh_ring5, cage_62512, cage_62512_rings, time);
                
                cage_62512_count = remove_duplicates_map(cage_62512_rings);
                
                //Write output gro file only if frameCounter is a multiple of DT. if -dt is not provided, every frame will be written.
                if(cage_62512_count>0 && (frameCounter % DT) == 0)print_vmd_cage64512_frings(cup62512, cage_62512_count, cage_62512_rings, ring5, ring6, atom_Pos, time, rawFilename, box_size_xyz, solutes, methane_62512, solute1, topSolute, solute2, count_solute2, frameCounter);
            }
            cout << "# 6²5¹²\tcage: " << cage_62512_count << "\n";
            /**********************************************/
            /*Finding Cage 64512*/
            int cage_64512_count = 0;
            
            cage_64512_count = cage_Finder_64512(cup62512, count_62512_cups, cage_64512_rings);
            
            cout << "# 6⁴5¹²\tcage: " << cage_64512_count << "\n\n";

            print_vmd_cage64512_frings(cup62512, cage_64512_count, cage_64512_rings, ring5, ring6, atom_Pos, time, rawFilename, box_size_xyz, solutes, methane_64512, solute1, topSolute, solute2, count_solute2, frameCounter);
            
            /**********************************************/

            outFile << cage_512_count << "\t " << methane_512 << "\t\t " << cage_62512_count << "\t " << methane_62512 << "\t\t " << cage_64512_count << "\t " << methane_64512 << endl ;

            if(in_F4 == 1)      //If F4 flag is provided, calculate F4.
            {
                F4_value = calc_F4(count_solvent, count_solute, My_neigh, atom_Pos, boxX, boxY, boxZ, Nneigh, Natoms, topSolute, time, HBOND_DIST) ;
                if(F4_value > 0) outFile_F4 << frameCounter << "\t " << F4_value << "\t " << time << endl;
                else outFile_F4 << frameCounter << "\t" << F4_value << "\t " << time << endl;        //The if-else condition takes care of the extra space needed for "-" sign and alligns the output(lazy way!).
                
            }
            
            
            
        }
    }   //End of reading the input file.
    
    
    fileIN.close();
    outFile.close();
    
    if( in_F4 == 1 )outFile_F4.close();
    
    return 0;
    
} // Close main function







