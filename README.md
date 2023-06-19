# GRADE 

![C++](https://img.shields.io/badge/Language-C%2B%2B-blue)

**version 1.00**

## Description
GRADE analyzes atomic positions of oxygen atoms of water to compute the number of $5^{12}$, $6^{2}5^{12}$, and $6^{4}5^{12}$ cages and account for their three-dimensional structures. The latter can be used for visualization using software such as VMD (Visual Molecular Dynamics). GRADE stands for "cages" in Portuguese. *`F4`* order parameter can also be calculated for trajectories.

To understand how the code works and for examples of the output files, please refer to the published paper.

If you have used this code, please cite the following publication:
Mahmoudinobar, Farbod, and Cristiano L. Dias. [*GRADE: A code to determine clathrate hydrate structures*](https://doi.org/10.1016/j.cpc.2019.06.004), Computer Physics Communications 244 (2019): 385-391



---


# Prerequisites
GNU Compiler Collection version 6.1.0 or newer.

# Compilation
GRADE is written in C++ and is made up of a main program file [GRADE.cpp](./GRADE.cpp) and two supporting resource files [MyFunctions.hpp](./MyFunctions.hpp) and [MyFunctions.cpp](./MyFunctions.cpp). Use the provided [Makefile](./Makefile) to compile GRADE by typing:




## Prerequisites:
GNU Compiler Collection  version 6.1.0 or newer.



## Compilation:

GRADE is written in C++ and is made up of a main program file (GRADE.cpp) and two supporting resource files (MyFunctions.hpp and MyFunctions.cpp). Clone the repo and change directory to `./PATH_TO_GRADE/GRADE/`. Use the Makefile to compile GRADE by typing: 

```.sh
$ make
$ make clean
```

# Usage: 

If the gro file containing atomic positions of water molecules is named “test.gro”, to run the code with all default values: 

```.sh
$ ./GRADE -i test.gro 
```

This will generate following files by default (if at least one cage is found in test.gro): test.xvg, test_cage512-frame.gro, test_cage62512-frame.gro (if 62512 cages exist).

A separate gro file is generated for each time-frame of the test.gro.

# Options

A full list of options can be printed on the terminal by using the flag `-h`. The available options as of version 1.00 are:

- `-i <.gro>` (input): Trajectory in gro format
- `-theta <int>` (45) (degree): Angle cut-off for planarity constraint
- `-r <real>` (0.35) (nm): Hydrogen bond cutoff radius (nm, Oxygen-Oxygen distance)
- `-d1 <real>` (0.18) (nm): Minimum length of Pentagon diameter
- `-d2 <real>` (0.26) (nm): Minimum length of Hexagon diameter
- `-o <.gro/.xvg>` (output) (Opt.): Output name to be used in ~.xvg, ~_cage512.gro, ~_cage62512.gro
- `-dt <int>` (1) (Opt.): Read all input file, write output gro files every dt frame(s)
- `-fr <int>` (1) (Opt.): Read input file every fr frame(s)
- `-[no]f4` (yes): Compute four-body order parameter F4=<cos3ф>


## Authors
- Farbod Mahmoudinobar, [fm59@njit.edu](mailto:fm59@njit.edu), [farbodmahmoudi@gmail.com](mailto:farbodmahmoudi@gmail.com)

## PI
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
<div style="position: relative;">
  <button class="copy-button" onclick="copyCode(this)">Copy</button>
  <pre><code>
	@article{Nobar2019grade,
		title={GRADE: A code to determine clathrate hydrate structures},
		author={Mahmoudinobar, Farbod and Dias, Cristiano L.},
		journal={Computer Physics Communications},
		volume={244},
		pages={385--391},
		year={2019},
		doi={10.1016/j.cpc.2019.06.004},
		publisher={Elsevier}
	}
	</code></pre>
</div>
