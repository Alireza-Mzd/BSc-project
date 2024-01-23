from EOS import z_cal

def f_saturation(COMPOSITION, EQUTYPE, EOS_num, TEMP, PRESS, TCRIT, PCRIT, RGAS, AF, NCOMP, BIC):
    
    phi_liq=[]
    phi_vap=[]

    for phase_type_counter in range(0,2): #liquid and vapor
        results = z_cal(EOS_num, phase_type_counter, COMPOSITION[phase_type_counter], TEMP, PRESS, TCRIT, PCRIT, NCOMP, BIC, RGAS, AF)
        if phase_type_counter==0:
            phi_liq=results[1]

        else:
            phi_vap=results[1]
            
    if EQUTYPE==0: #fs=fB
        COMP=COMPOSITION[0] #liquid
        fs1=sum([a*b/c for a,b,c in zip(COMP,phi_liq, phi_vap)])-1
    elif EQUTYPE==1:
        COMP=COMPOSITION[1] #vapor
        fs1=sum([a*b/c for a,b,c in zip(COMP,phi_vap, phi_liq)])-1
    else:
        print('check EQUTYPE value')

    return fs1
