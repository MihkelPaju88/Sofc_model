import cantera as ct
import math

# ==========================================
# 1. CANTERA INITIALIZATION (RESTORED TPB MECHANISM)
# ==========================================
gas = ct.Solution('gri30.yaml', 'gri30')
ni = ct.Solution('ni_gdc.yaml', 'Ni_Bulk')
lscf = ct.Solution('ni_gdc.yaml', 'LSCF_Bulk')
gdc = ct.Solution('ni_gdc.yaml', 'GDC_Bulk')

anode_tpb = ct.Interface('ni_gdc.yaml', 'Anode_TPB', [gas, ni, gdc])
cathode_tpb = ct.Interface('ni_gdc.yaml', 'Cathode_TPB', [gas, lscf, gdc])

A_v = 1.0e5
F = ct.faraday
R = ct.gas_constant
n_elec = 2.0

# ==========================================
# 2. FULLY SPATIAL TPB KINETICS
# ==========================================
def calculate_nodal_kinetics(T_anode, T_cathode, P, X_H2, X_H2O, X_O2, j_local, sigma_0=3.34e4):
    # 1. Nernst Potential across the cell (using average temp for global Delta G)
    T_avg = (T_anode + T_cathode) / 2.0
    gas.TPX = T_avg, P, f"H2:{max(X_H2, 1e-10)}, O2:{max(X_O2, 1e-10)}, H2O:{max(X_H2O, 1e-10)}"
    
    mu = gas.chemical_potentials
    delta_G = mu[gas.species_index('H2O')] - mu[gas.species_index('H2')] - 0.5 * mu[gas.species_index('O2')]
    E_nernst = -delta_G / (n_elec * F)
    
    # 2. Activation Overpotentials (Inverse Butler-Volmer)
    j_safe = max(abs(j_local), 1e-6)
    j0_anode = 8000.0 * math.exp(-120000 / (R * T_anode))
    j0_cathode = 5000.0 * math.exp(-110000 / (R * T_cathode))
    
    eta_anode = (2.0 * R * T_anode) / (2.0 * F) * math.asinh(j_safe / (2.0 * j0_anode))
    eta_cathode = (2.0 * R * T_cathode) / (2.0 * F) * math.asinh(j_safe / (2.0 * j0_cathode))
    eta_act_total = eta_anode + eta_cathode
    
    # 3. Two-Sided TPB Thermodynamics (Reversible Enthalpy + Irreversible Activation Heat)
    rate = j_safe / (n_elec * F)
    
    # Anode Heat 
    gas.TPX = T_anode, P, f"H2:{max(X_H2, 1e-10)}, H2O:{max(X_H2O, 1e-10)}"
    anode_tpb.TP = T_anode, P
    dh_anode = anode_tpb.delta_enthalpy[0]
    q_anode = A_v * (rate * -dh_anode) + A_v * (j_safe * eta_anode)
    
    # Cathode Heat 
    gas.TPX = T_cathode, P, f"O2:{max(X_O2, 1e-10)}"
    cathode_tpb.TP = T_cathode, P
    dh_cathode = cathode_tpb.delta_enthalpy[0]
    q_cathode = A_v * (rate * -dh_cathode) + A_v * (j_safe * eta_cathode)
    
    # 4. GDC Ionic Conductivity 
    sigma_GDC = (sigma_0 / T_avg) * math.exp(-62000.0 / (8.314 * T_avg))
    
    return E_nernst, eta_act_total, q_anode, q_cathode, sigma_GDC
