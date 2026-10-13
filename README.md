# ROV State-Space Modeling & Controllability Analysis

## Project Overview
This repository contains a MATLAB-based mathematical framework for analyzing the controllability and dynamics of a Remotely Operated Vehicle (ROV). The project leverages linear control theory to formulate the state-space representation of the ROV, evaluate its controllability under normal and degraded conditions, and simulate its state evolution using optimal control inputs.

> **Note:** [Add a screenshot of the MATLAB plots showing the X, Y, Yaw, and velocities tracking the reference state]

## Mathematical Architecture

### State-Space Representation
The ROV's dynamics are modeled as a linear time-invariant (LTI) system with a 6-dimensional state vector $X = [x, y, \theta, v_x, v_y, \omega]^T$:
$$\dot{X} = A X + B U$$

Where:
- $A \in \mathbb{R}^{6 \times 6}$ is the system matrix incorporating mass ($m$), inertia ($J$), and hydrodynamic drag coefficients ($\alpha, \beta, \gamma$).
- $B \in \mathbb{R}^{6 \times 4}$ is the input matrix derived from the Thruster Allocation Matrix (TAM) and DC-gain.

### Controllability Gramian
The system's ability to reach a desired final state $X_f$ within time $T$ is evaluated by computing the Controllability Gramian $W_c$:
$$W_c = \int_0^T e^{A(T-\tau)} B B^T e^{A^T(T-\tau)} d\tau$$

A full-rank Gramian ($rank(W_c) = 6$) guarantees complete state controllability.

## Core Features & Analysis

1. **Optimal Control Input Computation:**
   Calculates the exact voltage commands $U(\tau)$ required to drive the ROV from the origin to a reference state $X_f = [-1, 2, 3, 0, 0, 0]^T$ in $T=1$ second using the inverted Gramian.
2. **State Evolution Simulation:**
   Implements Euler integration ($\Delta t = 0.001s$) to simulate and plot the real-time evolution of the ROV's pose and velocities.
3. **Hardware Redundancy & Failure Mode Analysis:**
   - Evaluates system controllability when specific thrusters fail by nullifying respective columns in the Thruster Allocation Matrix (TAM).
   - Proves hardware redundancy: The ROV remains fully controllable if Thruster 2 fails, but loses controllability if both Thrusters 2 and 3 fail simultaneously.
4. **Subspace Analysis (SVD):**
   Applies Singular Value Decomposition (SVD) on the Gramian during uncontrollable failure modes to mathematically isolate the reachable and unreachable subspaces.

## How to Run

**Prerequisites:**
- MATLAB (Symbolic Math Toolbox required for exact Gramian integration)

**Execution:**
1. Clone the repository.
2. Open `Lab2_ex1.m` in MATLAB.
3. Run the script to execute the Gramian integrals, generate the trajectory plots, and output the thruster failure analysis to the command window.

---
*Developed by Hussein Elnaggar.*
