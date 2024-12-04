# README
Here, I have several programs that solve the 2D diffusion equation on a fixed Cartesial grid with a constant diffusion coefficient. The goal of this repo is to compare different algorithms and different hardware/parallelisation strategies. A brief overview of the different programs:

## Expl_Diff_Serial:

Solves the diffusion equation via a simple explicit FTCS algorithm. Compile the code serially via:
 
	`gfortran Expl_Diff_Serial.f90`

Or, compile the code with better vectorisation:

	`gfortran -O3 Expl_Diff_Serial.f90`

## Expl_Diff_OMP

Uses the same algorithm, but is complemented with openMP (OMP) pragma's for shared memory parallelisation. Compile via:

	`gfortran -fopenmp Expl_Diff_OMP.f90`

## Expl_Diff_Acc

Uses again the same algorithm, but makes use of openACC for GPU acceleration. To compile this program, first connect to a gpu-node, load the nvidia hpc module:

	`module load nvhpc`

And compile the code with the nvidia compiler:

	'nvfortran -acc Expl_Diff_Acc.f90


## Relax_Diff_Serial

Makes use of an implicit jacobi iteration algorithm. Compile serially via:

	`gfortran Relax_Diff_Serial.f90`

## 2g_Diff_Serial

Makes use of the implicit jacobi iteration algorithm, but on two spatial grids, accelerating the error diffusion of larger wavelengths slightly. Compile via

	`gfortran 2q_Diff_Serial`

## Mg_Diff_Serial

Implements the multi-grid algorithm for solving the diffusion equation. Compile via

	`gfortran Mg_Diff_Serial'
