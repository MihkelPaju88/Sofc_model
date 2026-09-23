import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile
import ufl
import precice
from dolfinx.fem import form, assemble_scalar
import csv
import sys
import os

# --- NEW: Link Cantera API ---
# Append the sibling directory to Python's system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'kinetics-cantera')))
import sofc_0d as kinetics

# --- 1. MESH GENERATION ---
domain = mesh.create_rectangle(
    MPI.COMM_WORLD,
    [np.array([0.0, 0.001]), np.array([0.05, 0.00101])],
    [50, 10],
    cell_type=mesh.CellType.quadrilateral
)


# --- 3. THERMAL & ELECTRICAL SOLVER SETUP ---
V = fem.functionspace(domain, ("Lagrange", 1))

# Thermal Variables
T = ufl.TrialFunction(V)
v_T = ufl.TestFunction(V)
T_n = fem.Function(V)
T_n.x.array[:] = 873.15
T_n.name = "Temperature"
T_next = fem.Function(V)

# --- NEW: Explicit Boundary Flux Arrays ---
q_anode_func = fem.Function(V)
q_anode_func.name = "Anode_Heat_Flux"
q_cathode_func = fem.Function(V)
q_cathode_func.name = "Cathode_Heat_Flux"

# Thermal Parameters
k = fem.Constant(domain, default_scalar_type(2.0))
rho = fem.Constant(domain, default_scalar_type(6000.0))
Cp = fem.Constant(domain, default_scalar_type(400.0))
dt = fem.Constant(domain, default_scalar_type(0.001))



# Electrical Variables
Phi = ufl.TrialFunction(V)
v_Phi = ufl.TestFunction(V)
Phi_n = fem.Function(V)
Phi_n.name = "Electrical_Potential"


# Electrical Parameters (GDC Arrhenius)
sigma_0 = fem.Constant(domain, default_scalar_type(3.34e4))
E_a = fem.Constant(domain, default_scalar_type(62000.0))
R_gas = fem.Constant(domain, default_scalar_type(8.314))

# Temperature-dependent Conductivity Field
sigma_GDC = (sigma_0 / T_n) * ufl.exp(-E_a / (R_gas * T_n))

# Electrical Weak Form (Poisson Equation: div(sigma * grad(Phi)) = 0)
a_Phi = sigma_GDC * ufl.inner(ufl.grad(Phi), ufl.grad(v_Phi)) * ufl.dx
L_Phi = fem.Constant(domain, default_scalar_type(0.0)) * v_Phi * ufl.dx

# Boundary Conditions
def inlet_boundary(x): return np.isclose(x[0], 0.0)
def bottom_boundary(x): return np.isclose(x[1], 0.001)
def top_boundary(x): return np.isclose(x[1], 0.00101)

fdim = domain.topology.dim - 1

# Inlet facets (kept exclusively for the Thermal solver's temperature boundary)
inlet_facets = mesh.locate_entities_boundary(domain, fdim, inlet_boundary)
inlet_dofs = fem.locate_dofs_topological(V, fdim, inlet_facets)

# Transverse facets (for the Electrical solver)
bottom_facets = mesh.locate_entities_boundary(domain, fdim, bottom_boundary)
bottom_dofs = fem.locate_dofs_topological(V, fdim, bottom_facets)

top_facets = mesh.locate_entities_boundary(domain, fdim, top_boundary)
top_dofs = fem.locate_dofs_topological(V, fdim, top_facets)

# --- NEW: Create Boundary Measure for Integration ---
marked_facets = np.hstack([bottom_facets, top_facets])
marked_values = np.hstack([np.full(len(bottom_facets), 1, dtype=np.int32), 
                           np.full(len(top_facets), 2, dtype=np.int32)])
sorted_facets = np.argsort(marked_facets)
facet_tags = mesh.meshtags(domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets])

ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
n = ufl.FacetNormal(domain)

# Integral of J dot n over the top boundary (tag 2)
I_form = form(ufl.inner(-sigma_GDC * ufl.grad(Phi_n), n) * ds(2))
# --------------------------------------------------

# --- NEW: Neumann Boundary Integration ---
# FEniCS requires ds(1) for bottom facets and ds(2) for top facets (mapped later in the code)
# The flux functions will be integrated correctly as W/m^2 over the surface boundaries.
F_T = rho * Cp * ((T - T_n) / dt) * v_T * ufl.dx \
  + k * ufl.inner(ufl.grad(T), ufl.grad(v_T)) * ufl.dx \
  - q_anode_func * v_T * ds(1) \
  - q_cathode_func * v_T * ds(2)
a_T, L_T = ufl.system(F_T)

# --- NEW: Dynamic Electrode Boundaries ---
V_operating = 0.6  # The external voltage constraint

bc_T_inlet = fem.dirichletbc(default_scalar_type(873.15), inlet_dofs, V)
bc_Phi_anode = fem.dirichletbc(default_scalar_type(0.0), bottom_dofs, V)

# Create a spatially varying function to hold the local driving forces
V_top_func = fem.Function(V)
V_top_func.name = "Cathode_Driving_Potential"
V_top_func.x.array[:] = 0.5  # Initial guess to start the solver

bc_Phi_cathode = fem.dirichletbc(V_top_func, top_dofs)

problem_T = LinearProblem(a_T, L_T, bcs=[bc_T_inlet], u=T_next, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}, petsc_options_prefix="T_solve")
problem_Phi = LinearProblem(a_Phi, L_Phi, bcs=[bc_Phi_anode, bc_Phi_cathode], u=Phi_n, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}, petsc_options_prefix="Phi_solve")

# Current Density Magnitude Space (Scalar for safe memory alignment)
J_ufl = -sigma_GDC * ufl.grad(Phi_n)
J_mag_ufl = ufl.sqrt(ufl.inner(J_ufl, J_ufl))

J_mag_expr = fem.Expression(J_mag_ufl, V.element.interpolation_points)
J_mag = fem.Function(V)
J_mag.name = "Current_Density_Magnitude"

# --- 4. PRECICE COUPLING ---
print("Initializing preCICE...")
participant = precice.Participant("Solid", "../precice-config.xml", MPI.COMM_WORLD.Get_rank(), MPI.COMM_WORLD.Get_size())

vertices = V.tabulate_dof_coordinates()
vertex_ids = participant.set_mesh_vertices("Solid-Mesh", vertices)
participant.initialize()
print("preCICE handshake complete!")

num_vertices = len(vertices)
write_flux = np.zeros(num_vertices)
T_checkpoint = np.copy(T_n.x.array)

xdmf = XDMFFile(domain.comm, "solid_results.xdmf", "w")
xdmf.write_mesh(domain)
t_solid = 0.0

# Initialize CSV logger
csv_file = None
csv_writer = None
if domain.comm.rank == 0:
    csv_file = open("sofc_metrics.csv", "w", newline="")
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(["Time_s", "Max_Temp_K", "Current_A", "Power_W"])

current_A = 0.1  # Initialize current to prevent errors on step 1

try:
    while participant.is_coupling_ongoing():
        
        if participant.requires_writing_checkpoint():
            T_checkpoint[:] = T_next.x.array[:]
        
        dt_precice = participant.get_max_time_step_size()

        # Read Data from OpenFOAM (PreCICE requires these reads to clear the buffer)
        read_temp = participant.read_data("Solid-Mesh", "Temperature", vertex_ids, dt_precice)
        read_h2   = participant.read_data("Solid-Mesh", "H2", vertex_ids, dt_precice)
        read_o2   = participant.read_data("Solid-Mesh", "O2", vertex_ids, dt_precice)
        read_h2o  = participant.read_data("Solid-Mesh", "H2O", vertex_ids, dt_precice)

        P_atm = 101325.0  # 1 atm in Pascals

        # --- NEW: Reordered Nodal Point-Wise Kinetics Loop ---
        # 1. Solve the Electrical Field FIRST so we have accurate physical currents
        problem_Phi.solve()
        J_mag.interpolate(J_mag_expr)

        # --- NEW: Robust Spatial Mapping & Under-Relaxation Loop ---
        coords = V.tabulate_dof_coordinates()
        x_anode = coords[bottom_dofs, 0]
        x_cathode = coords[top_dofs, 0]
    
        # 1. Robust Unstructured Mapping: Find the closest Anode node for every Cathode node
        closest_anode_indices = np.argmin(np.abs(x_cathode[:, None] - x_anode), axis=1)
        mapped_anode_dofs = bottom_dofs[closest_anode_indices]
    
        q_anode_func.x.array[:] = 0.0
        q_cathode_func.x.array[:] = 0.0
        write_flux[:] = 0.0
    
        omega = 0.5  # Under-relaxation factor for numerical stability
   
        # 2. Iterate through the dynamically mapped parallel nodes
        for idx, d_cathode in enumerate(top_dofs):
          d_anode = mapped_anode_dofs[idx]
        
          T_anode = max(T_n.x.array[d_anode], 300.0)
          T_cathode = max(T_n.x.array[d_cathode], 300.0)
        
          local_h2 = read_h2[d_anode]
          local_h2o = read_h2o[d_anode]
          local_o2 = read_o2[d_cathode]
        
          j_local = J_mag.x.array[d_cathode]
        
          E_nernst, eta_act, q_anode_vol, q_cathode_vol, _ = kinetics.calculate_nodal_kinetics(
            T_anode, T_cathode, P_atm, local_h2, local_h2o, local_o2, j_local
          )
        
          # 3. Under-Relaxation of the Galvanic Boundary
          V_calc = E_nernst - eta_act - V_operating
          V_current = V_top_func.x.array[d_cathode]
          V_top_func.x.array[d_cathode] = (1.0 - omega) * V_current + (omega * V_calc)
        
          # 4. Convert Cantera's volumetric heat (W/m^3) back to true interfacial flux (W/m^2)
          # This enables the mathematically correct Neumann ds integration.
          q_anode_flux = q_anode_vol / kinetics.A_v
          q_cathode_flux = q_cathode_vol / kinetics.A_v
        
          q_anode_func.x.array[d_anode] = q_anode_flux
          write_flux[d_anode] = q_anode_flux
        
          q_cathode_func.x.array[d_cathode] = q_cathode_flux
          write_flux[d_cathode] = q_cathode_flux

        # 5. Solve the Electrical Field with the relaxed boundary conditions
        problem_Phi.solve()
    
        # 6. Interpolate the new physical current
        J_mag.interpolate(J_mag_expr)
    
        #  7. Solve the Thermal Field using Neumann ds measures
        problem_T.solve()
        # -------------------------------------------
        # -------------------------------------------
        # -------------------------------------------
        # -------------------------------------------
      
        participant.write_data("Solid-Mesh", "Heat-Flux", vertex_ids, write_flux)
        participant.advance(dt_precice)
        
        if participant.requires_reading_checkpoint():
            T_next.x.array[:] = T_checkpoint[:]
        else:
            T_n.x.array[:] = T_next.x.array[:]
            t_solid += dt_precice
            
            # Export to XDMF (J_mag was already interpolated at the start of the time step)
            xdmf.write_function(T_n, t_solid)
            xdmf.write_function(q_anode_func, t_solid)
            xdmf.write_function(q_cathode_func, t_solid)
            xdmf.write_function(Phi_n, t_solid)
            xdmf.write_function(J_mag, t_solid)
            
            # --- NEW: Calculate Current and Power ---
            current_A = abs(domain.comm.allreduce(assemble_scalar(I_form), op=MPI.SUM))
            power_W = current_A * V_operating
            
            print(f"Time: {t_solid:.3f} s | Max Temp: {np.max(T_next.x.array):.2f} K | Current: {current_A:.2f} A/m | Power: {power_W:.2f} W/m")

# --- NEW: Write to CSV ---
            if domain.comm.rank == 0:
                max_T = np.max(T_next.x.array)
                csv_writer.writerow([t_solid, max_T, current_A, power_W])
                csv_file.flush()

finally:
    xdmf.close()
    participant.finalize()
    if domain.comm.rank == 0 and csv_file is not None:
        csv_file.close()
    print("Data files safely closed.")
