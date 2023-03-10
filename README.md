# Max-Series-Gain
This code acompanies the paper *Strengthened Circle and Popov cirteria for the stability analysis of feedback systems with ReLU neural networks* where the experimental setup is detailed in the *numerical examples* section. The code computes the maximum series gain for which global asymptotic stability is maintained using various criteria. This is tested on a number of example Lurie systems assumed to have repeated ReLU nonlinearities.  

### Authors:
* Carl R Richardson (cr2g16@soton.ac.uk)
* Matthew C Turner (m.c.turner@soton.ac.uk)

## Prerequisites
All the code is written in MATLAB. The LMI's are solved using the *Robust Control Toolbox* which must be installed as an add-on.

## Overview
The repository is organised as follows:
- `Examples.m` The master script. It defines the example systems and parameters of the various criteria. It then loops through each example, computing the maximum series gain according to each criterion,  and displays the results.
- `Circle.m` Implementation of the Circle criterion - See Theorem 1 and Remark 4.
- `Circle_Like.m` Implementation of the Circle-Like criterion - See Theorem 1.
- `Popov.m` Implementation of the Popov criterion - See Theorem 2 and Remark 6.
- `Popov_Like1.m` Implementation of the Relaxed Popov-Like criterion - See Corollary 1.
- `Popov_Like2.m` Implementation of the Relaxed Popov-Like criterion - See Corollary 2.
- `Park.m` Implementation of the Park criterion - See Theorem 2 from Reference 13.
- `ZF.m` Implementation of the Zames-Falb criterion - See References 24 and 25.

## Getting Started
Run `Examples.m` to repeat the experiments in the paper or select a subset of the examples by defining them in the *Ex_array* variable.  
