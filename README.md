# composites_calculators by PC

Calculators related to composite materials

The repository contains Python scripts that use Classical Lamination Theory for Composite Materials to find different parameters of laminated composites.

composites_equivalent_sitffness_calculator.py

File composites_equivalent_sitffness_calculator.py is a Python script written using OOP to determine equivalent parameters of laminated composite.
The input for this file are Young's moduli in directions of fibres and matrix, E11 and E22, respectively, Poisson's ratio v12 and shear modulus of lamina G12,
list of angles "thetas" (fibers orientations) and list of ply thicknesses "thickness". 

The script computes ABD matrices of the laminate and equivalent Young's moduli, Poisson's ratio and shear modulus G12. 

Current limitation of the script is that only one material can be used (which in 80% of composites is sufficient), therefore no hybrid materials are considered yet. 
