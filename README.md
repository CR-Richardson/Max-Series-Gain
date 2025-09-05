# Max-Series-Gain
This respository was first developed to acompany the paper [Strengthened Circle and Popov Criteria for the stability analysis of feedback systems with ReLU neural networks](https://eprints.soton.ac.uk/478495/) where the experimental setup is detailed in the *numerical examples* section. The code computes the maximum series gain for which global asymptotic stability is verified using various criteria. Furthermore, the number of decision variables used by the criteria is also returned to compare the complexity. The criteria are tested on a number of example Lurie systems assumed to have repeated ReLU nonlinearities. The code associated with this paper was developed using MATLAB's *Robust Control Toolbox* and the LMIs were solved using *LMI Lab*. This code is contained in the `LMI_Toolbox` directory.

Since then, a second, more efficient implementation has been developed using *YALMIP*, to pose the problems, and *MOSEK*, to solve the LMIs. This implementation is contained in the `YALMIP_MOSEK` directory and additionally demonstrates the criteria on some higher dimensional Hopfield network examples, which were intractable with the previous implementation.

### Authors:
* Carl R Richardson (cr2g16@soton.ac.uk)
* Matthew C Turner (m.c.turner@soton.ac.uk)

## Prerequisites
All the code is written in MATLAB. The `LMI_Toolbox` implementation requires the *Robust Control Toolbox* to be installed as a MATLAB add-on whereas the `YALMIP_MOSEK` implementation requires [YALMIP](https://yalmip.github.io/download/) and [MOSEK](https://www.mosek.com/) to be installed.

## Overview
The `Docs` directory contains associated documents explaining the code at varying levels of detail. The following file structure is used within the `LMI_Toolbox` and `YALMIP_MOSEK` directories:
- `Max_Series_Gain.m` The master script which loops through each example, computing the maximum series gain (and # of decision variables) according to each criterion.
- `Examples.m` Defines the (A,B,C,D) matrices of the example systems.
- `Circle.m` Implementation of the Circle Criterion - See Theorem 1 and Remark 3.
- `Circle_Like.m` Implementation of the Circle-like Criterion - See Theorem 1.
- `Popov.m` Implementation of the Popov Criterion - See Theorem 2 and Remark 4.
- `Popov_Like1.m` Implementation of the Relaxed Popov-like Criterion - See Corollary 1.
- `Popov_Like2.m` Implementation of the Relaxed Popov-like Criterion - See Corollary 2.
- `Park.m` Implementation of the Park Criterion - See Theorem 2 from Reference 11.
- `ZF.m` Implementation of the Zames-Falb Criterion - See Reference 31.
- `Aizerman.m` Calculates the Nyquist gain for each example.

## Getting Started
Run `Max_Series_Gain.m` to repeat the experiments in the paper or select a subset of the examples by defining them in the *Ex_array* variable.  
