# GRADE 

![C++](https://img.shields.io/badge/Language-C%2B%2B-blue)

**version 1.00**

## Description
GRADE analyzes atomic positions of oxygen atoms of water to compute the number of $5^{12}$, $6^{2}5^{12}$, and $6^{4}5^{12}$ cages and account for their three-dimensional structures. The latter can be used for visualization using software such as VMD (Visual Molecular Dynamics). GRADE stands for "cages" in Portuguese. *`F4`* order parameter can also be calculated for trajectories.

To understand how the code works and for examples of the output files, please refer to the published paper.

If you have used this code, please cite the following publication:
Mahmoudinobar, Farbod, and Cristiano L. Dias. [*GRADE: A code to determine clathrate hydrate structures*](https://doi.org/10.1016/j.cpc.2019.06.004), Computer Physics Communications 244 (2019): 385-391



---

# Using Docker to run the code (recommended)

You can use the docker image of this code to skip compiling the c++ source code by first pulling the docker image from `dockerhub` and then run the code using the docker image. You can follow these steps:

1. Download and install [Docker](https://www.docker.com)

2. Pull Grade's docker image from [Dockerhub](https://hub.docker.com/repository/docker/farbodnobar/grade/general):
```.sh
docker pull farbodnobar/grade:v1.0
```
3. Go to the directory where you store your input .gro files:
```
cd ~/YOUR_WORKING_DIR/
```
4. Find the `IMAGE ID` corresponding to Grade's docker image:
```
image_id=$(docker images -q farbodnobar/grade:v1.0)
```
5. Execute Grade through the docker image:
```
docker run --rm -v $PWD:/data -w /data $image_id grade -i test.gro 
```
The last command runs a Docker container using Grade's docker image, attaches it to your current working directory, and sets the container's working directory to `/data`. It runs the Grade code and after completion of the command, removes the container on exit. The output files of Grade will be written in your current working directory. 

# Compiling from source code
GRADE is written in C++ and is made up of a main program file [GRADE.cpp](./GRADE.cpp) and two supporting resource files [MyFunctions.hpp](./MyFunctions.hpp) and [MyFunctions.cpp](./MyFunctions.cpp). Use the provided [Makefile](./Makefile) to compile GRADE by typing:



## Prerequisites:
GNU Compiler Collection  version 6.1.0 or newer.



## Compilation:

GRADE is written in C++ and is made up of a main program file (GRADE.cpp) and two supporting resource files (MyFunctions.hpp and MyFunctions.cpp). Clone the repo and change directory to `./PATH_TO_GRADE/GRADE/`. Use the Makefile to compile GRADE by typing: 

```.sh
$ make
$ make clean
```

## Usage: 

If the gro file containing atomic positions of water molecules is named “test.gro”, to run the code with all default values: 

```.sh
$ ./GRADE -i test.gro 
```

This will generate following files by default (if at least one cage is found in test.gro): test.xvg, test_cage512-frame.gro, test_cage62512-frame.gro (if 62512 cages exist).

A separate gro file is generated for each time-frame of the test.gro.

---

# Options

A full list of options can be printed on the terminal by using the flag `-h`. The available options as of version 1.00 are:

| Option             | Default Value | Description                                                      |
|--------------------|---------------|------------------------------------------------------------------|
| `-i <.gro>`        |  	         | Trajectory in gro format                             |
| `-theta <int>`     | 45 (degree)	 | Angle cut-off for planarity constraint                 |
| `-r <real>`        | 0.35 (nm)     | Hydrogen bond cutoff radius (nm, Oxygen-Oxygen distance)    |
| `-d1 <real>`       | 0.18 (nm)     | Minimum length of Pentagon diameter                         |
| `-d2 <real>`       | 0.26 (nm)     | Minimum length of Hexagon diameter                          |
| `-o <.gro/.xvg>`   | 			     | Output name to be used as prefix in ~.xvg, ~_cage512.gro, ~_cage62512.gro files|
| `-dt <int>`        | 1 (Opt.)      | Read all input file, write output gro files every dt frame(s) |
| `-fr <int>`        | 1  (Opt.)     | Read input file every fr frame(s)                         |
| `-[no]f4`          | yes           | Compute four-body order parameter F4=<cos3ф>. Use `-f4 no` to disable.                     |



## Authors
- Farbod Mahmoudinobar, [fm59@njit.edu](mailto:fm59@njit.edu), [farbodmahmoudi@gmail.com](mailto:farbodmahmoudi@gmail.com)

- Cristiano L. Dias, [cld@njit.edu](mailto:cld@njit.edu)

## License & Copyright
- © Cristiano L. Dias, New Jersey Institute of Technology, Physics
Licensed under the GNU GPL-3.0-or-later


License & Copyright:

- © Cristiano L. Dias, New Jersey Institute of Technology, Physics

Licensed under the GNU GPL-3.0-or-later

---

## Notes:
- Grade is optimized to work with .gro files written by GROMACS. Trajectories converted from other MD packages may encounter issues with reading the data. 
- Currently, Grade recognizes 4-atom water models.

---
## Citation
```
@article{Nobar2019grade,
	title={GRADE: A code to determine clathrate hydrate structures},
	author={Mahmoudinobar, Farbod and Dias, Cristiano L.},
	journal={Computer Physics Communications},
	volume={244},
	pages={385--391},
	year={2019},
	doi={10.1016/j.cpc.2019.06.004},
	publisher={Elsevier}
```
