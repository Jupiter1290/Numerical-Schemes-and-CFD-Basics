# Introduction to CFD using MATLAB and OpenFOAM

This repository contains MATLAB programs developed as part of the course **“Introduction to CFD using MATLAB and OpenFOAM.”**  
Each folder represents a self-contained exercise focused on understanding the mathematical foundations of Computational Fluid Dynamics (CFD) through numerical methods and algorithmic implementation.  
The `README` files within each folder include the corresponding problem statements.

---

## 📂 Repository Structure and Topics

| Folder | Description |
|--------|--------------|
| `2DHeatConduction_TransientAndSteady` | Solves the 2D steady and transient heat conduction equation using point iterative techniques — Jacobi, Gauss–Seidel, and Successive Over-Relaxation (SOR). Boundary conditions are fixed-temperature (Dirichlet). |
| `4thOrder_Derivative_Schemes` | Derives and implements fourth-order central and skewed finite difference schemes for second-order derivatives, verifying the expected order of accuracy using analytical comparison. |
| `Discretisation_Basics` | Compares first-, second-, and fourth-order finite difference approximations for the first derivative of *f(x) = sin(x)/x³*, analyzing numerical vs analytical errors through bar plots. |
| `LinearConvection_ExplicitMethod` | Simulates the 1D linear convection equation using first-order explicit schemes. Includes grid refinement and time step studies to demonstrate stability and CFL dependence. |
| `ODE45_Routine_Pendulum` | Models nonlinear pendulum motion using MATLAB’s ODE45 solver. Demonstrates time integration of coupled ODEs with damping effects. |
| `Quasi1D_Super_SubSonicNozzle_ConsAndNonCons` | Simulates quasi-one-dimensional compressible flow through a converging–diverging nozzle using the MacCormack scheme in both conservative and non-conservative forms. Compares steady-state and transient solutions. |
| `StabilityInLinearSystems` | Solves a system of linear equations using Jacobi, Gauss–Seidel, and SOR methods. Computes iteration matrices, eigenvalues, and spectral radii to analyze convergence. |
| `Backward facing step-Ganymede.docx` | Technical report on a backward-facing step flow simulation performed in OpenFOAM. |

---

## ⚙️ Numerical Methods and Tools
- **Programming:** MATLAB  
- **Techniques:** Finite Difference Method (FDM), Explicit Time Stepping, Iterative Solvers (Jacobi, Gauss–Seidel, SOR), Stability and Convergence Analysis  
- **Applications:** Heat Conduction, Linear Convection, Compressible Nozzle Flow, and Matrix Systems  

---

## 🧠 Learning Outcomes
- Developed an understanding of numerical discretization and stability.  
- Implemented key CFD solvers from scratch in MATLAB.  
- Transitioned from foundational numerical schemes to OpenFOAM, leading into the course *“Introduction to OpenFOAM Development.”*

---

## 🧑‍💻 Author
**Thomas Abraham**  
B.Tech in Mechanical Engineering, IIT (ISM) Dhanbad  
*(Batch of 2026)*  
Specialization: Computational Fluid Dynamics (CFD)
