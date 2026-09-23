import os

def write_field(filename, f_class, dims, internal, b_inlet, b_outlet, b_wall):
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {f_class};
    object      {filename};
}}
dimensions      {dims};
internalField   uniform {internal};
boundaryField
{{
    fuel_inlet        {{ {b_inlet} }}
    fuel_outlet       {{ {b_outlet} }}
    anode_interface   {{ {b_wall} }}
    bottom_wall       {{ {b_wall} }}
    
    air_inlet         {{ {b_inlet} }}
    air_outlet        {{ {b_outlet} }}
    cathode_interface {{ {b_wall} }}
    top_wall          {{ {b_wall} }}
    
    front             {{ type empty; }}
    back              {{ type empty; }}
}}
"""
    with open(f"0/{filename}", "w") as f:
        f.write(content)

# 1. Velocity (U): Drop to 0.5 m/s to ensure Courant < 1.0. Set internal field to match.
write_field(
    "U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0.5 0 0)",
    "type fixedValue; value uniform (0.5 0 0);",  # Inlets
    "type zeroGradient;",                         # Outlets
    "type noSlip;"                                # Walls/Interfaces
)

# 2. Temperature (T): Match 873.15 K everywhere to prevent thermal shock
write_field(
    "T", "volScalarField", "[0 0 0 1 0 0 0]", "873.15",
    "type fixedValue; value uniform 873.15;",     # Inlets
    "type zeroGradient;",                         # Outlets
    "type zeroGradient;"                          # Walls (preCICE overrides this)
)

# 3. Pressure (p): 1 atm (101325 Pa) everywhere.
write_field(
    "p", "volScalarField", "[1 -1 -2 0 0 0 0]", "101325",
    "type zeroGradient;",                         # Inlets (velocity is fixed, pressure must float)
    "type fixedValue; value uniform 101325;",     # Outlets
    "type zeroGradient;"                          # Walls
)
