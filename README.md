# Computational Data for the Lefschetz Property of Graphic Matroids

[![DOI](https://zenodo.org/badge/1166160242.svg)](https://doi.org/10.5281/zenodo.18764513)

## Overview
This dataset contains computational results regarding the Strong Lefschetz Property (SLP) for Artinian Gorenstein algebras associated with graphic matroids. These data support the findings presented in the paper: "Failure of the Lefschetz property for the Graphic Matroid".

## Environment and Software
- System: SageMath Version 9.5

## Directory Structure
- /basis: CSV files containing basis elements for the algebra A_3, represented as 3-tuples of variable indices.
- /F: Text files containing the components of the polynomial vector F, which generates the kernel of the Hessian matrix.
- /graphs: Text files providing edge lists (vertex set {0..7}).
- /H: CSV files of the higher Hessian matrices H_{B_3}(f_G).
- /polynomial: Text files of the basis generating functions f_G.
- table.csv: A comprehensive index of all 152 graphs, including the number of edges, Hilbert function values (dim A_2, dim A_3), degrees of F, edge sets in compressed binary format, and graph6 representations.
- verify.sage: A SageMath script to re-verify the algebraic relations (H * F = 0, etc.) and data consistency.

## Data Field Descriptions (table.csv)
- index: Serial number of the graph (1 to 152).
- num_edges: Number of edges in the graph (corresponds to dim A_1).
- dim_a2, dim_a3: Dimensions of the algebra in degrees 2 and 3, respectively.
- deg_f_1, deg_f_2: Degrees of the polynomial components in the kernel vector F. 
- edge_set_binary: 28-bit binary string representing the edges (subset of all pairs in {0..7}).
- graph6: Standard graph6 representation.

## How to Verify
To verify the data using SageMath, run:
$ sage verify.sage

This script will check the consistency between the various graph representations, re-calculate f_G and H from the matroid structure and the basis B_3, and confirm that H * F = 0.
