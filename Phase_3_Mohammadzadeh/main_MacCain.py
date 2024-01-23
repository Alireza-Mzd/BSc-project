# import math
# import numpy as np
from EOS import z_cal
from sat_fun import f_saturation
from equilibrium import equil
# from McCain_results import result_calculator

# Please be asware that this code, calculates the z factor, phi_i, fugacity and fs with 5 different EOSs (vdW, RK, SRK, PR, SW)

# #Data
TEMP=160+460 # Rankine
PRESS=1000 #psi
RGAS=10.73

CompName=['C1',	'NC4',	'C10']
COMP_liq=[0.2408, 0.1517, 0.6075] #liquid compostion
COMP_vap=[0.9613, 0.0366, 0.0021] #vapor composition
COMPOSITION=[COMP_liq, COMP_vap]

TCRIT=[343, 765.3, 1111.7] #from McCain, Rankine
PCRIT=[666.4, 550.6, 305.2] #from McCain, psi
AF=[0.0104,0.1995,0.4898] # ACENTRIC FACTOR

NCOMP=len(COMP_liq)

#BIC from mcCain
BIC=[[0, 0.02, 0.04],
     [0.02, 0, 0],
     [0.04, 0, 0]
    ]

z_liq=[]
phi_liq=[]
z_vap=[]
phi_vap=[]
fil=[] #fugacity of liquid of component i
fiv=[] #fugacity of vapor of component i


# NOTE: THE ABOVE LISTS (i.e. phi_liq, phi_vap, fil, fiv) CONTAIN 5 ROWS EACH HAVNIG THREE(=NCOMP) ELEMENTS.
#    EACH ROW IS FOR AN EOS AND EACH ELEMENT WITHTIN A ROW IS FOR THE COMPONENT i.
# NOTE: FOR Z FACTORS, EACH ROW IS FOR ONE EOS (VD, RK, SRK, PR, SW)

# NOTE: the results of the following loop are mainly for demonstration
for phase_type_counter in range(0,2): #liquid and vapor
    for EOS_num_counter in range(0,5): # only calculates based on PR and SW EoS 
        results=z_cal(EOS_num_counter, phase_type_counter, COMPOSITION[phase_type_counter], TEMP, PRESS, TCRIT, PCRIT, NCOMP, BIC, RGAS, AF)
        if phase_type_counter==0:
            z_liq.append(results[0])
            phi_liq.append(results[1])
            fil.append(results[2]) #psi
        else:
            z_vap.append(results[0])
            phi_vap.append(results[1])
            fiv.append(results[2]) #psi


fs=[]
ki_equil=[]
x_equil=[]
y_equil=[]
fil_equil=[]
fiv_equil=[]
variable_P_T=[]
for EQUTYPE in range(0,2):
    for VARTYPE in range(0,2):
        #just to demonstrate that 'f_saturation' works correctly
        fs.append(f_saturation(COMPOSITION, EQUTYPE, 3, TEMP, PRESS, TCRIT, PCRIT, RGAS, AF, NCOMP,BIC))
        
        ki_equil, x_equil, y_equil, fil_equil, fiv_equil, variable= equil(3, EQUTYPE, VARTYPE, COMPOSITION, TEMP, PRESS, TCRIT, PCRIT, NCOMP, BIC, RGAS, AF)
        variable_P_T.append(variable)
print('output results= ',variable_P_T)
        