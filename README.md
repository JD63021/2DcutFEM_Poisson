# 2D CutFEM Poisson Solver: High-Order P3 Implementation with Circular Inclusions

## Description

This repository contains an advanced MATLAB implementation of the Cut Finite Element Method (CutFEM) for solving the 2D Poisson/Laplace equation. The solver is specifically designed to handle multiple circular Dirichlet boundaries that are unfitted, meaning they do not align with the underlying mesh edges. By using higher-order P3 (cubic) elements and a true-arc Nitsche integration scheme, the code achieves high accuracy even on relatively coarse meshes.

## Mathematical Problem

The solver addresses the steady-state diffusion equation:

-kappa * (div(grad(u))) = f

In this specific setup (Laplace):

1. Inner Boundary 1 (Hot): A circle at (0,0) with Radius = 1.0 and u = 50.
2. Inner Boundary 2 (Cold): A circle at (2.5,0) with Radius = 0.5 and u = 0.
3. Outer Boundary: Default "do-nothing" Neumann boundary condition.

The problem effectively simulates the temperature distribution or potential flow between two cylinders.

## Dependent Files

To run this solver, the following files must be present in your MATLAB path:

1. mesh5_gmsh.m: A custom utility function used to parse and load high-order Gmsh files into the MATLAB environment.
2. circlep3r.m: The mesh definition file (typically generated via Gmsh) containing the P3 triangle (10-node) and edge (4-node) connectivity.

## Technical Features

1. P3 Isoparametric Elements: The solver utilizes 10-node cubic triangles (P3) for high-order spatial accuracy.
2. Volume Masking: Integration points (Gauss points) located inside the circular inclusions are identified and skipped during the assembly of the global stiffness matrix.
3. True-Arc Nitsche Integration: Instead of approximating the circular boundary as a straight chord, the code integrates along the actual circular arc within each cut element. This maintains the geometric fidelity of the circles.
4. Newton-Based Point Inversion: Since P3 mapping is non-linear, the code uses a Newton-Raphson scheme to map physical boundary points back to the isoparametric (xi, eta) space of the element.
5. Strong Interior Clamping: Nodes located strictly inside the circles are "clamped" to the Dirichlet value. This prevents the global system matrix from becoming singular due to disconnected degrees of freedom inside the "cut" regions.
6. Analytical Validation: The code includes a closed-form analytical solution based on bipolar coordinates (Greenberg formulation) for direct error comparison.

## Neighborhood Visualization

The script includes a sophisticated visualization tool that extracts the "interface neighborhood." It identifies cut elements and their immediate neighbors (k-ring neighbors) and renders them with:

* Shrinkage effects to visualize the gaps between elements.
* Sub-element partitioning showing how cut triangles are split into smaller polygons for visualization.
* Distinct color palettes for inside vs. outside sub-regions.

## Configuration Parameters

Key variables in the Controls section:

* gmshFile: The name of the P3 mesh file.
* order: Polynomial order (set to 3 for P3).
* kappa: Thermal conductivity or diffusion coefficient.
* gammaN: Nitsche penalty parameter (set to 30 for stability).
* C(1), C(2): Structures defining the center, radius, and boundary value for the circles.

## How to Run

1. Ensure mesh5_gmsh.m and circlep3r.m are in your directory.
2. Run cutfem2d_P3_two_circles_ARC_mask_clamp.m.
3. The script will generate a tiled layout comparing the CutFEM solution to the analytical solution, followed by contour plots and a detailed interface neighborhood view.

---
