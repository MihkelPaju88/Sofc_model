# **SOFC Digital Twin: High-Fidelity Multiphysics Coupling**

A spatially resolved, fully coupled digital twin for Intermediate-Temperature Solid Oxide Fuel Cells (IT-SOFC). This model escapes the limitations of 0D lumped-parameter approximations by explicitly coupling 2D fluid dynamics, 2D solid-state heat and charge transport, and localized Triple Phase Boundary (TPB) surface chemistry.

## **Multiphysics Architecture**

This digital twin utilizes a staggered, explicit coupling architecture mapped across three specialized open-source solvers via **preCICE**.

> * **Fluid Domain (OpenFOAM):** Simulates multi-component gas advection and diffusion across distinct Cartesian Anode ($H\_2$, $H\_2O$) and Cathode ($O\_2$, $N\_2$) channels.  
> * **Solid Domain (FEniCSx):** Solves the non-linear Poisson equations for Ohmic charge transport and Fourier heat conduction through the dense GDC electrolyte.  
> * **Thermodynamics & Kinetics (Cantera):** Acts as the point-wise electrochemical engine. It calculates strictly localized Nernst driving potentials, Butler-Volmer activation overpotentials, and exact half-cell reversible enthalpies at the gas-solid interfaces.

### **The Coupling Loop**

> 1. **PreCICE** extracts localized fuel and oxidant species fractions from the OpenFOAM boundary patches.  
> 2. **FEniCS** maps parallel 2D unstructured fluid nodes to its solid electrolyte boundaries via nearest-neighbor distance matrices.  
> 3. **Cantera** evaluates the local Nernst potential based on the dual gas streams and computes the exact interfacial heat flux (W/m²), strictly splitting the endothermic $H\_2$ oxidation from the highly exothermic $O\_2$ reduction.  
> 4. **FEniCS** applies these galvanic driving potentials to solve the internal Ohmic current ($J$), integrates the interfacial heat via Neumann boundary measures (ds), and solves the thermal field.  
> 5. **PreCICE** pushes the resulting heat fluxes back to OpenFOAM to dictate the convective fluid temperature rise for the next time step.

## **Repository Structure**

&nbsp;

&nbsp;

&nbsp;

Plaintext

├── fluid-openfoam/       \# OpenFOAM 2D Cartesian gas channels  
│   ├── 0/                \# Initial boundary and internal fields (p, U, T, species)  
│   ├── constant/         \# Thermophysical and species transport properties  
│   └── system/           \# blockMeshDict, controlDict, and preCICE fluid adapter  
├── solid-fenics/         \# FEniCSx solid electrolyte solver  
│   ├── sofc\_solid.py     \# Main finite element thermal/electrical solver  
│   └── precice-config.xml\# Coupling definition and data exchange mapping  
└── kinetics-cantera/     \# Cantera 0D Surface Thermodynamics  
&nbsp;&nbsp;&nbsp;&nbsp;├── sofc\_0d.py        \# Python API for localized nodal kinetics  
&nbsp;&nbsp;&nbsp;&nbsp;├── ni\_gdc.yaml       \# Custom TPB Anode/Cathode heterogeneous mechanism  
&nbsp;&nbsp;&nbsp;&nbsp;└── gri30.yaml        \# Standard gas-phase combustion mechanism

## **Dependencies**

To run the coupled simulation, your environment must have the following installed:

> * **FEniCSx** (dolfinx) with PETSc linear solver backends.  
> * **OpenFOAM** (v2312 or newer) with the openfoam-adapter for preCICE compiled.  
> * **Cantera** (v3.0.0+) for Python.  
> * **preCICE** (v3.4.1+) and pyprecice bindings.  
> * **MPI** (mpi4py) for parallel execution.

## **Installation & Setup**

> 1. **Clone the repository:**  
>    Bash  
>    git clone https://github.com/your-username/sofc-digital-twin.git  
>    cd sofc-digital-twin

> 2. **Initialize the Fluid Mesh:**  
>    Navigate to the OpenFOAM directory and generate the strictly segregated Anode/Cathode channels.  
>    Bash  
>    cd fluid-openfoam  
>    blockMesh  
>    setFields

>    *Note: setFields is strictly required to safely initialize the fuel and air channels and prevent Courant-violating thermodynamic advection shocks at $t=0$.*

## **Running the Simulation**

Because this is a dynamically coupled simulation, the solid and fluid solvers must be executed simultaneously in separate terminal windows. preCICE will automatically synchronize the time steps via socket communication.

**Terminal 1 (Fluid Solver):**

&nbsp;

&nbsp;

&nbsp;

Bash

cd fluid-openfoam  
reactingFoam

**Terminal 2 (Solid Solver):**

&nbsp;

&nbsp;

&nbsp;

Bash

cd solid-fenics  
python3 sofc\_solid.py

## **Post-Processing**

> * **Solid-State Fields:** FEniCS exports results to solid-fenics/solid\_results.xdmf. Use ParaView to visualize the dynamic spatial gradients of Temperature, Electrical Potential ($\\Phi$), and Current Density Magnitude ($\\vert{}J\\vert{}$).  
> * **Global Metrics:** A lightweight CSV logger (sofc\_metrics.csv) tracks maximum cell temperature, integrated cell current (A/m), and operational power (W/m) in real-time.  
> * **Fluid Fields:** Standard OpenFOAM time directories are generated in fluid-openfoam/. Use paraFoam to visualize the depletion of $H\_2$ and the convective transport of heat out of the cell.

## **Recent Architectural Upgrades**

> * **Neumann Surface Integration:** Eliminated mesh-dependent volumetric heat scaling by applying strictly interfacial Cantera heat fluxes via ufl.Measure("ds").  
> * **Under-Relaxation:** Implemented explicit $\\omega \= 0.5$ dampening on the galvanic boundary conditions to stabilize highly non-linear Nernst/Butler-Volmer oscillations.  
> * **Two-Sided Thermodynamics:** Rewrote ni\_gdc.yaml to enforce true spatial half-reactions, isolating Anode endothermic cooling from Cathode exothermic heating.