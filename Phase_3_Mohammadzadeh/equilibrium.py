import math
from EOS import z_cal
from sat_fun import f_saturation
def equil(EOS_num_equil, EQUTYPE, VARTYPE, COMPOSITION, TEMP, PRESS, TCRIT, PCRIT, NCOMP, BIC, RGAS, AF):
    fil=[]
    fiv=[]
    xi=[]
    yi=[]
    ki=[]
    if VARTYPE==0: #pressure is the variable and temperature is known
        T=TEMP
        #for now
        P=950
        # P=float(input('Please enter your initial guess for PRESSURE(psia): '))
        var=P
    else: #temperature is the variable and pressure is known
        P=PRESS
        #for now
        T=600
        # T=float(input('Please enter your initial guess for TEMPERATURE(Rankine): '))
        var=T
        
    #intitial k value using WILSON eq.
    ki=[a/PRESS*math.exp(5.37*(1+b)*(1-c/TEMP)) for a,b,c in zip(PCRIT, AF, TCRIT)]
        
    ERRORsat=0.1
    while ERRORsat>10**-5:
        ERRORequil=0.1
        while ERRORequil>10**-13:
            if EQUTYPE==0: # bubble calculation
                nv=0
                xi=COMPOSITION[0] #sum(xi)=1
                yi=[a*b for a,b in zip(ki,xi)] #in such case sum(yi/ki)=1
            else: # dew calculation
                nv=1
                yi=COMPOSITION[1]
                xi=[a/b for a,b in zip(yi,ki)]
            Z_l, phi_liq, fil=z_cal(EOS_num_equil, 0, xi, T, P, TCRIT, PCRIT, NCOMP, BIC, RGAS, AF)
            Z_v, phi_vap, fiv=z_cal(EOS_num_equil, 1, yi, T, P, TCRIT, PCRIT, NCOMP, BIC, RGAS, AF)
            # ERRORequil = abs(sum([1-a/b for a,b in zip(phi_liq,phi_vap)]))
            ERRORequil = abs(sum([1-a/b for a,b in zip(fil,fiv)]))
            ki=[a/b for a,b in zip(phi_liq,phi_vap)] # according to McCain formula for Ki
            
        sat_COMP=[xi,yi]
        step=0.02
        fs =f_saturation(sat_COMP, EQUTYPE, EOS_num_equil, T, P, TCRIT, PCRIT, RGAS, AF, NCOMP,BIC)
        
        if VARTYPE==0:
            fs_plus =f_saturation(sat_COMP, EQUTYPE, EOS_num_equil, T, P+step, TCRIT, PCRIT, RGAS, AF, NCOMP,BIC)
            fs_minus=f_saturation(sat_COMP, EQUTYPE, EOS_num_equil, T, P-step, TCRIT, PCRIT, RGAS, AF, NCOMP,BIC)
        else:
            fs_plus =f_saturation(sat_COMP, EQUTYPE, EOS_num_equil, T+step, P, TCRIT, PCRIT, RGAS, AF, NCOMP,BIC)
            fs_minus=f_saturation(sat_COMP, EQUTYPE, EOS_num_equil, T-step, P, TCRIT, PCRIT, RGAS, AF, NCOMP,BIC)
        
        derivative=(fs_plus-fs_minus)/(2*step)
        var=var-fs/derivative
        
        #updating T or P
        if VARTYPE==0: #pressure is the variable and temperature is known
            P=var
        else: #temperature is the variable and pressure is known
            T=var
        
        ERRORsat = f_saturation(sat_COMP, EQUTYPE, EOS_num_equil, T, P, TCRIT, PCRIT, RGAS, AF, NCOMP,BIC)
        # ERRORsat=10**-13
    if VARTYPE==0:
        return ki, xi, yi, fil, fiv, P
    else:
        return ki, xi, yi, fil, fiv, T
           
