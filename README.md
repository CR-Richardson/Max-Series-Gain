# Max-Series-Gain
This code acompanies the paper [Strengthened Circle and Popov cirteria for the stability analysis of feedback systems with ReLU neural networks](10.1109/LCSYS.2023.3287494) where the experimental setup is detailed in the *numerical examples* section. The code computes the maximum series gain for which global asymptotic stability is maintained using various criteria. Furthermore, the number of decision variables used by the criteria is also returned to compare the complexity. This is tested on a number of example Lurie systems assumed to have repeated ReLU nonlinearities.  

### Authors:
* Carl R Richardson (cr2g16@soton.ac.uk)
* Matthew C Turner (m.c.turner@soton.ac.uk)

## Prerequisites
All the code is written in MATLAB. The LMI's are solved using the *Robust Control Toolbox* which must be installed as an add-on.

## Overview
The repository is organised as follows:
- `Max_Series_Gain.m` The master script. It loops through each example, computing the maximum series gain (and # of decision variables) according to each criterion,  and displays the results.
- `Examples.m` Defines the (A,B,C,D) matrices of the example systems.
- `ZF_Parameters.m` Defines the Zames-Falb parameters used for each example.
- `Circle.m` Implementation of the Circle criterion - See Theorem 1 and Remark 3.
- `Circle_Like.m` Implementation of the Circle-Like criterion - See Theorem 1.
- `Popov.m` Implementation of the Popov criterion - See Theorem 2 and Remark 4.
- `Popov_Like1.m` Implementation of the Relaxed Popov-Like criterion - See Corollary 1.
- `Popov_Like2.m` Implementation of the Relaxed Popov-Like criterion - See Corollary 2.
- `Park.m` Implementation of the Park criterion - See Theorem 2 from Reference 11.
- `ZF.m` Implementation of the Zames-Falb criterion - See Reference 31.

## Getting Started
Run `Max_Series_Gain.m` to repeat the experiments in the paper or select a subset of the examples by defining them in the *Ex_array* variable.  
