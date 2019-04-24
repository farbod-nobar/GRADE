---------------------------------------------------------------------------------------------------------------------
GRADE 
**version 1.00**
GRADE analyzes atomic positions of oxygen atoms of water to compute the number of 512, 62512 and 64512cages and account for their three-dimensional structures. The latter can be used for visualization using software such as VMD (Visual Molecular Dynamics). GRADE stands of “cages” in Portuguese. F4 order parameter can also be calculated for trajectories.
---------------------------------------------------------------------------------------------------------------------
Prerequisites: 
-GNU Compiler Collection  version 6.1.0 or newer.
---------------------------------------------------------------------------------------------------------------------
Compilation:
GRADE is written in C++ and is made up of a main program file (GRADE.cpp) and two supporting resource files (MyFunctions.hpp and MyFunctions.cpp). Use the Makefile to compile GRADE by typing: 
$ make
$ make clean
---------------------------------------------------------------------------------------------------------------------
Usage: 
If the gro file containing atomic positions of water molecules is named “test.gro”, to run type: 
$ ./GRADE -i test.gro 
This will generate following files by default (if at least one cage is found in test.gro): test.xvg, test_cage512-frame.gro, test_cage62512-frame.gro (if 62512 cages exist). A separate gro file is generated for each time-frame of the test.gro. 
---------------------------------------------------------------------------------------------------------------------
Options:
Full list of options can be printed on terminal by using flag ‘-h’. These options as of version 1.00 are:
-i 	[<.gro>] 	(input)
	Trajectory in gro format
-theta 	<int> 	(45) 	(degree)
	Angle cut-off for planarity constraint
-r 	<real> 	(0.35) 	(nm)
	Hydrogen bond cutoff radius 	(nm, Oxygen-Oxygen distance)
-d1 	<real> 	(0.18) 	(nm)
	Minimum length of Pentagon diameter
-d2 	<real> 	(0.26) 	(nm)
	Minimum length of Hexagon diameter
-o 	[<.gro/.xvg>] 	(output) 	 (Opt.)
	(output name to be used in ~.xvg, ~_cage512.gro, ~_cage62512.gro)
-dt 	<int> 	(1) 	(Opt.)
	Read all input file, write output gro files every dt frame(s)
-fr 	<int>	(1)	(Opt.)
	Read input file every fr frame(s)
-[no]f4 	(yes)
	Compute four-body order parameter F4=<cos3ф>


---------------------------------------------------------------------------------------------------------------------

Authors: 
Farbod Mahmoudinobar, fm59@njit.edu
Cristiano L. Dias, cld@njit.edu
---------------------------------------------------------------------------------------------------------------------
License & Copyright:
© Cristiano L. Dias, New Jersey Institute of Technology, Physics
Licensed under the GNU GPL-3.0-or-later



