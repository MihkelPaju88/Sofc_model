import os

def write_specie(filename, internal, f_in, a_in):
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      {filename};
}}
dimensions      [0 0 0 0 0 0 0];
internalField   uniform {internal};
boundaryField
{{
    fuel_inlet        {{ type fixedValue; value uniform {f_in}; }}
    fuel_outlet       {{ type zeroGradient; }}
    anode_interface   {{ type zeroGradient; }}
    bottom_wall       {{ type zeroGradient; }}
    
    air_inlet         {{ type fixedValue; value uniform {a_in}; }}
    air_outlet        {{ type zeroGradient; }}
    cathode_interface {{ type zeroGradient; }}
    top_wall          {{ type zeroGradient; }}
    
    front             {{ type empty; }}
    back              {{ type empty; }}
}}
"""
    with open(f"0/{filename}", "w") as f:
        f.write(content)

# Initialize the internal field with N2 to ensure sum=1.0 at time 0. 
# The H2 and O2 will flush into the channels immediately as the solver runs.
write_specie("H2",  0.0, 0.97, 0.0)
write_specie("H2O", 0.0, 0.03, 0.0)
write_specie("O2",  0.0, 0.0,  0.233)
write_specie("N2",  1.0, 0.0,  0.767)
