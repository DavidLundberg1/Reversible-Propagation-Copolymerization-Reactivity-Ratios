# Accurate Determination of Reactivity Ratios for Copolymerization Reactions with Reversible Propagation Mechanisms
This repository supports the following manuscript:

David J. Lundberg, Landon J. Kilgallon, Julian C. Cooper, Francesca Starvaggi, Yan Xia, Jeremiah A. Johnson, "Accurate Determination of Reactivity Ratios for Copolymerization Reactions with Reversible Propagation Mechanisms", _Macromolecules_, 2024. [Link](https://doi.org/10.1021/acs.macromol.4c00835)

In this work we compare three distinct methods for the characterization of copolymerization reactions with reversible propagation mechanisms, and assess their general accuracy, usability, and consistency. 

This repository is intended to allow polymer chemists and engineers to employ the methods developed in this work to their own data.  This repository contains general MatLab scripts for running stochastic simulations, population balance models, and numerical copolymer equation integrations given inputs of rate constants, monomer concentrations, final polymerization conversion, and chain length. Additionally, an example script used to fit the endo-NB/iPrSi8 dataset from the manuscript using copolymer equation integration is included. 

### List of Scripts
_IzuCPE.m_ <br /> 
Numerical integration of the Izu copolymer equation.<br /> 

_BinaryReversibleStochasticSimulation.m_<br /> 
Stochastic simulation of binary copolymerization with reversible propagation.<br /> 

_PopulationBalanceODEs.m_<br /> 
Numerical integration of population balance equations.<br /> 

_Example_iPrSi8_EndoNB_Fitting.m_<br /> 
Example of non-linear least squares fitting of experimental data using numerical integration of the Izu copolymer equation. 

## Contact
  David Lundberg, PhD 
  
  Email: DavidJL@mit.edu
  
  GitHubID: DavidLundberg1
