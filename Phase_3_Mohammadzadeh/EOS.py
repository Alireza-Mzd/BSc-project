import numpy as np
import math

# this function calculates z factor, phi_i, and fugacity of a liquid or vapor (depending on phase type) and \
    #based on 1 EOS specified by EOS_num
def z_cal(EOS_num, phase_type, COMP, TEMP, PRESS, TCRIT, PCRIT, NCOMP, BIC, RGAS, AF):
    
 
    #     vdW        RK      SRK       PR     SW is calculated differently
    omegaA=[0.421875, 0.42747, 0.42747, 0.457235, 0] # Danesh
    omegaB=[0.125, 0.08664, 0.08664, 0.07780, 0] # Danesh
    z1=[]
    SI=[0]*NCOMP
    bmix=0
    amix=0
    umix=0
    wmix=0
    aci=[]
    ai= []
    bi=[]
    alpha=[]
    ui=[]
    wi=[]
    qi=[]
    eta=[]
    m=0
    m0=0
    m1=0
    m2=0
    omegaaci=[]
    omegaBi=[]
    # COMP=COMPOSITION[phase_type]
    
    for ICOMP in range (0,NCOMP): # EOS initialization (aci, bi, ai, alpha, ui, wi etc.)
        Tr=TEMP/TCRIT[ICOMP]
        if EOS_num!=4:
            aci.append(omegaA[EOS_num]*RGAS**2*TCRIT[ICOMP]**2/PCRIT[ICOMP])
            bi.append(omegaB[EOS_num]*RGAS*TCRIT[ICOMP]/PCRIT[ICOMP])
        
        if EOS_num==0: #vdW
             alpha.append(1)
             ui.append(0)
             wi.append(0)
             myType='vdW'
        elif EOS_num==1: #RK
             alpha.append(Tr**(-0.5))
             ui.append(bi[ICOMP])
             wi.append(0)
             myType='RK'
        elif EOS_num==2: #SRK
             m=0.48+1.574*AF[ICOMP]-0.176*AF[ICOMP]**2
             alpha.append((1+m*(1-Tr**0.5))**2)
             ui.append(bi[ICOMP])
             wi.append(0)
             myType='SRK'
        elif EOS_num==3:#PR
             m=0.37464+1.54226*AF[ICOMP]-0.26992*AF[ICOMP]**2 
             alpha.append((1+m*(1-Tr**0.5))**2)
             ui.append(2*bi[ICOMP])
             wi.append(bi[ICOMP])
             myType='PR'
        else: #SW
            myType='SW'
            q_coefficients=[6*AF[ICOMP]+1, 3, 3, -1]
            q_roots=np.roots(q_coefficients)
            qi.append(np.real(q_roots[np.isreal(q_roots)][0]) if np.any(np.isreal(q_roots)) else None)
            eta.append(1/(3*(1+qi[ICOMP]*AF[ICOMP])))
            omegaBi.append(eta[ICOMP]*qi[ICOMP])
            omegaaci.append((1-eta[ICOMP]*(1-qi[ICOMP]))**3)
            bi.append(omegaBi[ICOMP]*RGAS*TCRIT[ICOMP]/PCRIT[ICOMP])
            aci.append(omegaaci[ICOMP]*RGAS**2*TCRIT[ICOMP]**2/PCRIT[ICOMP])
            ui.append(1+3*AF[ICOMP])
            wi.append(math.sqrt(3*AF[ICOMP]*bi[ICOMP]))
            # for i in  range (0,NCOMP):
            if Tr<1: #sub-critical fluid
            
                if  AF[ICOMP]<= 0.3671:
                    m0=0.465+1.347*AF[ICOMP]-0.528*AF[ICOMP]**2
                else:
                    m0=0.5361+0.9591*AF[ICOMP]
                if AF[ICOMP]<= 0.4:
                    m1=m0+0.01429*(5*Tr-3*m0-1)**2
                    m=m1
                elif AF[ICOMP]>=0.55:
                    m2=m0+0.71*(Tr-0.779)**2
                    m=m2
                else:
                    m1=m0+0.01429*(5*Tr-3*m0-1)**2
                    m2=m0+0.71*(Tr-0.779)**2
                    m=((0.55-AF[ICOMP])/0.15)*m1+((AF[ICOMP]-0.4)/0.5)*m2
                alpha.append((1+m*(1-Tr**0.5))**2)
            else: # if Tr>=1 supercritical fluid
                alpha.append(1-(0.4774+1.328*AF[ICOMP])*math.log(Tr))   
            
        ai.append(aci[ICOMP]*alpha[ICOMP])
    
    for ICOMP in range(0,NCOMP): # amix calculation
        for JCOMP in range(0,NCOMP):
            amix += (ai[ICOMP]*ai[JCOMP])**0.5*COMP[ICOMP]*COMP[JCOMP]*(1-BIC[ICOMP][JCOMP])
            
    for ICOMP in range(0,NCOMP): # SI calculation
    
        for JCOMP in range(0,NCOMP):
           SI[JCOMP] +=(2*ai[JCOMP]**0.5*((COMP[ICOMP] * ai[ICOMP]**0.5)*(1-BIC[ICOMP][JCOMP]))) / (amix)
            
    umix=sum([a*b for a,b in zip(COMP,ui)])
    wmix=sum([a*b for a,b in zip(COMP,wi)])
    bmix=sum([a*b for a,b in zip(COMP,bi)])
    A=amix*PRESS/(RGAS*TEMP)**2
    B=bmix*PRESS/(RGAS*TEMP)
    U=umix*PRESS/(RGAS*TEMP)
    W=wmix*PRESS/(RGAS*TEMP)
    
    a1 = 1
    a2 = -(1+B-U)
    a3 = A-B*U-U-W**2
    a4 = -(A*B-B*W**2-W**2)
    coefficients=[a1, a2, a3, a4]
    ROOTS=np.roots(coefficients)
    z1=[np.real(element) for element in ROOTS if np.isreal(element)]
    if len(z1)!=1:
        if phase_type==0:
            z1=min(z1)
        else:
            z1=max(z1)
        
    z1=z1[0]
    
    LnPhi=[0]*NCOMP
    PHII=[0]*NCOMP
    fug=[0]*NCOMP
    #u2=U**2+4*W**2
    
    # print(myType)
    # print('z1= ',z1)
    # print('B= ',B)
    for i in range(0,NCOMP):
        if wmix==0:
            wmix=1
            
        if umix==0:
            umix=1
            
        W=wmix*PRESS/(RGAS*TEMP)
        U=umix*PRESS/(RGAS*TEMP)
        t5=(U**2 + 4*W**2)
        t1= - math.log(z1-B) + (B*(bi[i]/bmix))/(z1-B)
        t2=A/math.sqrt(t5)
        t3=SI[i]
        t4=(ui[i]/umix)*U**2 + 4*(wi[i]/wmix)*W**2
        
        t6=2*z1+U - math.sqrt(t5)
        t7=2*z1+U + math.sqrt(t5)
        t8=2*(2*z1+U) * (wi[i]/wmix)*W**2
        t9=(U*z1 - 2*W**2)*(ui[i]/umix)*U
        t10=(z1**2+U*z1-W**2)*(t5)
        
        LnPhi[i]= t1 + (t2 * (t3 - t4/t5) * math.log(t6/t7)) - A*((t8 + t9)/t10)
      
        PHII[i]=math.exp(LnPhi[i])
        fug[i]=PRESS*COMP[i]*PHII[i]

    return z1, PHII, fug