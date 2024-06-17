$title "CAMD of Novel Solvent for C02 Capture"
$onEolCom
* COURSE TITLE:  **ADVANCED PROCESS OPTIMIZATION (CENG 70003)**

********Unit of measurement for properties*****
* Melting Point Temperature (Tmp)= K          *
* Boiling Point Temperature (Tbp)= K          *
* Critical Temperature (Tct)= K               *
* Critical Pressure (Pct)= bar                *
* Molar volume (Vt)= m3/kmol                  *
* Density(Rhou)= g/cm3                        *
* RED                                         *
* Liquid Heat Capacity (cpla)= cal/mol.K      *
* Vapour Pressure (Vap)= bar                  *
* Surface Tension (Sigma_h)= mN/m             *
* Viscosity (eta_h)= mPa.s                    *
* Molecular weight (Mwt) = kg/kmol            *
*********************************************** 

Sets
i functional groups to be considered in our design space
/CH3, CH2,CH, C, CH-NH2, CH2-NH2,CH3-NH, CH2-NH , CH2-N, CH3-N, OH, CH-NH, C-NH2/
j property contibution of each functional groups 
/tb, tc, REDrhod, REDrhop, REDrhoh, tm, cpa, cpb, cpc, cpd, mm, vm, sigma, eta, OW, ai, bi, ci, di, pc/
p properties for use in equations
/ Tmp, Tbp, Tct, Pct, Vt, MWt, Rhou, Ra, RED, Omega, cpla, Vap, sigma_h, eta_h /

x define the number of integer cuts / 1*10/
dyn(x) dynamic set of x for integer cuts
;
dyn(x) = NO;

parameters
Tb0 Default Boiling Point
                         /244.5165/
Tc0 Default Critical Temperature
                               /181.6716/
Pc1 Default critical Pressure constant 1
                                       /0.0519/
Pc2 Default critical pressure constant 2
                                       /0.1347/
Tm0 Default Melt Temp
                    /143.5706/
*Lquid molar volume constant by Hukkerika et al  (m3/kmol)                  
Vm0 Default Molar Volume
                       /0.0160/

*RED property constant using Hukkerikar approach
D0  Dispersion cohesive energy density of solvent
                                               /15.7/
P0  Dispersion polar cohesive energy density
                                           /5.2/
H0  Dispersion hydrogen bonding cohesive energy density
                                                     /5.8/
R0  Default RED of component to be captured
                                          /3.3/
*T   Temperature for project
Tavg  column temperature
                       /353/

Tabs  average temperature of the absorption column
                                                   /313/
Tdes  average temperature of the desorption column
                                                  /393/
*conversion factor for properties                                                  
CFr conversion factor for molecular weight from kgm-3 to gcm-3
                                                          /1e-3/
*Property 
pL(p) Defines lower bound for the property of interest
      /RED 10e-5, cpla 50, Rhou 0.6, Vap 10e-4/

pU(p) Defines upper bound for the property of interest
    /RED 6.5, cpla 110, Rhou 1.5, Vap 0.1/
*Scaling range
pmax upper scale 
   /1/
pmin lower scale
   /0/
                                         
v(i) Defines the valency & number of bonds on group i of the building block molecular groups
    /CH3 1, CH2 2,  CH 3, C 4, CH-NH2 2, CH2-NH2 1, CH3-NH 1, CH2-NH 2, CH2-N 3, CH3-N 2, OH 1, CH-NH 3, C-NH2 3/
nv(i,x)     Molecule ID values
OFv(x)      Objective function values 
kv(p,x)     Property values
;

Table Propertyset(i,j) Table  of constants for the building-block molecular groups
               tb         tc         REDrhod       REDrhop     REDrhoh      tm          cpa          cpb        cpc          cpd            mm         vm           sigma      eta          OW         ai           bi           ci          di         pc                         
    CH3       0.8853     1.0893       7.5697        1.9996       2.2105    0.6699      1.95e+1   -8.08E-03    1.53E-4     -9.67E-08       15.034     0.0241        8.0328    -1.0278      0.0077    0.0570      -0.2383e-2   0.7556e4    -0.1765     0.0052        
    CH2       0.5815     3.4607      -0.0018       -0.1492      -0.215     0.2992     -9.09e-1    9.50E-02   -5.44E-05     1.19E-08       14.026     0.0165        0.6213     0.2125      0.0047   -0.1497       0.0060e-2   1.4157e4     0.0751     0.0087
    CH       -0.0039     4.667       -7.7208       -2.7099      -2.6826   -0.2943     -2.30e+1    2.04E-1    -2.65E-04     1.20E-07       13.019     0.0086       -7.7843     1.3180      0.0005   -2.2942       0.4028e-2   4.5094e4     0.6679     0.0123
    C        -0.4985     6.6186      -15.4498      -4.7191      -6.4821   -0.043      -6.62e+1    4.27E-1    -6.41E-4      3.01E-07       12.011     0.0007       -16.3927    2.8147     -0.0069    1.0031      -0.3677e-2  -6.0316e4     1.1972     0.0150
    CH-NH2    1.5629     11.2234     -0.3287        0.6603       2.8953    37.8806     3.9        0.1628     -1.01e-4      2.24e-8        29.041     0.0207        4.4131     1.0108      0.0192   -3.4287       0.1902e-2   11.5638e4    0.8015     0.0057
    CH2-NH2   2.3212     12.9223      8.1717        5.2964       6.7984    3.4368      25.991     0.0538      1.096e-4    -8.57e-8        30.049     0.0281        16.395     0.2902      0.0286   -1.2842      -0.2066e-2   8.4701e4     0.2087    -0.0043
    CH3-NH    1.9861     10.6027      8.1301        3.4132       7.2551    2.7205      18.29      0.06812     2.6e-4      -8.62e-8        30.049     0.0282        18.3072   -0.0637      0.0206   -6.8919      -0.4106e-2   6.5360e4     1.4702     0.0030                   
    CH2-NH    1.3838     9.7537       0.2374        0.1072       1.4183    2.0673     -2.119      0.1712     -1.03e-4      2.24e-8        29.041     0.0260        4.5002     1.0512      0.0250   -7.0986      -0.1663e-2   7.1961e4     1.7218     0.0154
    CH2-N     0.4199     7.5227      -6.6418        0.0847      -7.3014   -0.1982     -32.009     0.322      -3.744e-4     1.579e-7       28.033     0.0187       -4.8682     1.4376      0.0077   -2.2900       0.4902e-2   7.6050e4     0.5469     0.0383
    CH3-N     1.0505     8.1028       0.4260       -0.5166       2.4585    0.9396     -11.6       0.21892    -1.67e-4      4.93e-8        29.041     0.0259       -2.6557     0.8715      0.0113   -2.0833       0.2459e-2   6.9449e4     0.2953     0.0152                       
    OH        2.1385     10.1673      8.0236        4.9598       11.8005   3.2702      2.57e+1   -6.91e-2     1.77e-4     -9.88e-8        17.007     0.0044        16.0184    1.3057      0.0502    1.4351      -1.0010e-2   13.8366e4    0.3418    -0.0146
    CH-NH     0.7116     8.7789      -7.7581       -3.5886      -2.2824    1.6571     -2.421e+1   0.2802     -3.136e-4     1.305e-7       28.033     0.0209       -5.207      1.8378      0.0136   -9.2431       0.2305e-2   10.2898e4    2.3146     0.0199
    C-NH2     1.0473     11.0799     -6.9473        0.7950       2.5957    11.4243    -39.3       0.3858     -4.77e-4      2.034e-7       28.033     0.0179       -7.2285     2.8880      0.0085   -0.1314      -0.5803e-2   1.0228e4     1.3308     0.0198

       
                                                                              

        
Variables
OF objective value
Rhou density
Tmp Melting point of solvent
Vt Molar Volume
MWt Molecular weight
Tbp boiling point
Ra red value which determines the ability of solvent to dissolve in solute
RED final red value
tau
cpoa ideal gas cp
cpla liquid heat capacity
sigma_h surface tension
eta_h  viscosity
Vap vapour pressure
Tavgr reduced average temperature
Tct  critical temperature
omega accentric factor
f1 Ambrose & Walton coefficient
f2 Ambrose & Walton coefficient
f3 Ambrose & Walton coefficient
Pct Crtitical pressure
k(p) Bringing the property values into a variable
Alpha
Tbr
Betta
;

Integer Variable
n(i) number per molecular groups;

Positive Variable
MWt
Vt
Tavgr
Tct
RED
Vap
Pct
Rhou;


Equations
eq1  estimation of melting point temperature of solvent 'Tmp' using Hukkerika et al apporach
eq2  estimation of boiling point temperature of solvent 'Tbp' using Hukkerika et al approach
eq3  estimation of critical temperature of solvent Tct using Hukkerika et al approach
eq3b estimation of critical pressure of solvent Pct using Hukkerika et al approach
eq4a estimation of total molar volume Vt using Hukkerika et al approach
eq4b estimation of total molecular weight MWt
eq4  estimation of density of solvent 'Rhou' using Hukkerika et al approach
eq5  estimation of RED of solvent using Hukkerika et al apporach
eq5a estimation of ability of solvent to dissolve solute Ra
eq6  estimation of heat capacity of solvent 'Cpla'  
*eq6a estimation of accentric factor
eq6b estimation of ideal gas heat capacity using Joback and Reid
eq7  estimation of vapour pressure of solvent 'Vap' using Ambrose and Walton coefficients
eq7a estimation of vapour pressure of solvent 'Vap' using Ambrose and Walton coefficients       
eq7b estimation of vapour pressure of solvent 'Vap' using Ambrose and Walton coefficients       
eq7c estimation of vapour pressure of solvent 'Vap' using Ambrose and Walton coefficients      
eq7d estimation of vapour pressure of solvent 'Vap' using Ambrose and Walton coefficients      
eq8  estimation of surface tension of solvent 'sigma' using Conte et al approach
*eq9a etimation of viscosity of solvent 'eta_h' using conte et al approach*
eq9b estimation of viscosity of solvent 'eta_h' using Hsu et al approach
eq10 constraint to ensure Tm is less than the average temperature of the absorption column 𝑇Abs (313 K) to avoid solvent solidification
eq11 constraint to ensure the solvent’s normal boiling temperature 𝑇b is higher than the average temperature of the desorption column 𝑇des (393 K)to avoid excess vaporization of the pure solvent
eq13 constraint on RED to ensure higher dissolution of solute
eq15 to estimate reduced average temperature Tavgr
eq16 constraint to ensure octect rule is satisfied
eq17 constraint to ensure bounding rule is satisfied
eq18 objective function 'OF'
eq19 Constraint on reduced average temperature Tavgr
eq20 Accentic factor calculation
eqAux1 constraint to avoid trivial solution
eq20(x) Integer cut
eq21 Estimation of Accentric factor using Sahinidis Approach
eq21c Estimation of Accentric factor using Sahinidis Approach
eq21b Estimation of Accentric factor using Sahinidis Approach
eq22  Estimation of Accentric factor using Sahinidis Approach
;
*Objective function 
*eq18..  OF =e=  0.25*(power(((k('RED')-REDT)/REDT),2)) + 0.25*(power(((k('cpla')-cplaT)/cplaT ),2)) + 0.25*(power(((k('Vap')-VapT)/VapT ),2)) + 0.25*((power(((k('Rhou')-RhouT)/RhouT ),2)))  ;
eq18..  OF =e=  0.25*(((k('RED')-pL('RED'))/(pU('RED')- pL('RED')))*(pmax-pmin)) + 0.25*((((k('cpla')-pL('cpla'))/(pU('cpla')- pL('cpla')))*(pmax-pmin))) + 0.25*((((k('Vap')-pL('Vap'))/(pU('Vap')- pL('Vap')))*(pmax-pmin))) - 0.25*((((k('Rhou')-pL('Rhou'))/(pU('Rhou')- pL('Rhou')))*(pmax-pmin)));
*estimation of melting point Tm
eq1.. k('Tmp') =e= Tm0 *log(sum(i,n(i)*Propertyset(i,'tm')));

*estimation of boiling point Tb
eq2.. k('Tbp') =e= Tb0 *log(sum(i,Propertyset(i,'tb')*n(i))); 

*estimation of critical temperature of solvent
eq3.. k('Tct') =e= Tc0 *  log(sum(i,Propertyset(i,'tc')*n(i)));

*estimation of critical pressure
eq3b.. rPower((k('Pct')-Pc1),-0.5) =e=(sum(i,Propertyset(i,'pc')*n(i))+Pc2);

*estimation of density of solvent
eq4a.. k('Vt')=e=(sum(i,Propertyset(i,'vm')*n(i))+ Vm0);
eq4b.. k('MWt')=e=sum(i,Propertyset(i,'mm')*n(i));
eq4..  k('Rhou') =e= (k('MWt')/k('Vt'))*CFr;

*estimation of RED
eq5a.. k('Ra') =e= (sqrt((4 * (power((15.7 - (sum(i,Propertyset(i,'REDrhod')* n(i)))), 2))) + (power((5.2 - (sum(i,Propertyset(i,'REDrhoP')* n(i)))), 2)) + (power((5.8 - (sum(i,Propertyset(i,'REDrhoH')* n(i)))),2))));
eq5..  k('RED') =e= k('Ra')/R0;

*etsimation of accentric factor omega using hukkerika et al
*eq6a.. k('omega') =e= Wa*((log(sum(i, Propertyset(i,'OW')*n(i))+ Wc))**1/Wb); 
*estimation of ideal gas heat capacity
eq6b..cpoa =e= (sum(i,Propertyset(i,'cpa')*n(i)))- 37.93 +(( ( sum(i,Propertyset(i,'cpb')*n(i))) + 0.21)* Tavg ) + (((sum(i,Propertyset(i,'cpc')*n(i))) - 3.91e-4) * (Tavg**2)) + (( ( sum(i,Propertyset(i,'cpd')*n(i))) + 2.06e-7)*(Tavg**3) );
*Estimation of liquid heat capacity of solvent
eq6..k('cpla') =e= (1/4.1868)*(cpoa + 8.314*(1.45 + (0.45 /(1-Tavgr)) + 0.25 * k('omega')*(17.11 + 25.2 *(rpower((1-Tavgr),(1/3))/Tavgr) + 1.742/(1-Tavgr))));

*estimation of vapour pressure, Vap
eq7a.. f1 =e= ((-5.97616*tau)+(1.29874*(rpower(tau,1.5)))-(0.60394*(rpower(tau,2.5)))-(1.06841*(rpower(tau,5))))/Tavgr;
eq7b.. f2 =e= ((-5.03365*tau)+(1.11505*(rpower(tau,1.5)))-(5.41217*(rpower(tau,2.5)))-(7.46628*(rpower(tau,5))))/Tavgr;
eq7c.. f3 =e= ((-0.64771*tau)+(2.41539*(rpower(tau,1.5)))-(4.26979*(rpower(tau,2.5)))+(3.25259*(rpower(tau,5))))/Tavgr;
eq7d.. Tau + Tavgr -1 =e= 0;
eq7.. k('Vap') =e= (exp(f1 + k('omega')*f2 + k('omega')*k('omega')*f3))*k('Pct');

*estimation of surface tension
eq8.. k('sigma_h') =e= sum(i,Propertyset(i,'sigma')*n(i));

*estimation of viscosity by conte et al
*eq9a.. k('eta_h')  =e= exp(sum(i, Propertyset(i,'eta')*n(i)));

*estimation of viscosity by Hsu et al
eq9b.. k('eta_h') =e= exp(sum(i,Propertyset(i,'ai')*n(i))+ sum(i,Propertyset(i,'bi')*n(i))*Tavg +  sum(i,Propertyset(i,'ci')*n(i))/(Tavg**2) + sum(i,Propertyset(i,'di')*n(i))*log(k('Pct')));


*Constraint on melting point temperature
eq10.. k('Tmp')  =l= Tabs ;

*Constraint on boiling point of solvent
eq11.. k('Tbp')  =g=  Tdes;

*RED constraint
eq13.. k('RED') =l= 1;

*Tavgr reduced average temperature
eq15..  Tavgr =e= Tavg/k('Tct');

*Valence rule(octet rule)
eq16.. sum(i,n(i)*(2-v(i))) =e= 2;

*bonding rule  n(i)*(Propertyset(i,'valency')-1) + 2 - sum(i, n(i)) =l=0;
ALIAS (i, ii);
eq17(i)..  n(i)*((v(i)) - 1)+ 2 - sum(ii, n(ii))=l=0 ;

eq19.. Tavgr =l=1;

*To avoid trivial solutions
eqAux1..         sum(i,n(i))=g=1;
*Estimation of accentric factor by Sahinidis approach
eq21..  Tbr=e= k('Tbp')/k('Tct');
eq21b..Alpha =e= -5.97214-log(k('pct')/1.013) + (6.09648/Tbr)+ 1.28862 *log(Tbr)-0.169347*((Tbr)**6);
eq21c..Betta =e= 15.2518-(15.6875/Tbr)- 13.4721*log(Tbr) + 0.43577*((Tbr)**6);
eq22.. k('omega') =e= Alpha/(Betta);


*Evaluation of integer cut using integer variable
eq20(x)$ (dyn(x)).. sum(i,abs(nv(i,x)-n(i)))=g=1;
*eq20(x)$ (dyn(x)).. sum(i,(rpower((nv(i,x)-n(i)),2)))=g=1;

*Boundary conditions
*bound on surface tension (Lee et al)
k.lo('sigma_h') = 25;
k.up('sigma_h') = 60;
*boundary on viscosity (Lee et al)
k.lo('eta_h') = 10e-5;
k.up('eta_h') = 60;
*boundary on Tb (Lee et al)
k.lo('Tbp') = 393;
k.up('Tbp') = 550;
*boundary on Tm (Lee et al)
k.lo('Tmp') = 273;
k.up('Tmp') = 313;
*boundary on group occurrence
n.lo(i)=0;
n.up(i)=5;
*boundary on RED
k.lo('RED') = pL('RED');
k.up('RED') = pU('RED');
*boundary on molecular weight
k.lo('MWt')=30;
k.up('MWt')=1000;
*boundary on vapour pressure
k.lo('Vap')=pL('Vap');
k.up('Vap')=pU('Vap');
*boundary on liquid heat capacity
k.lo('cpla')=pL('cpla');
k.up('cpla')=pU('cpla');
*boundary on critical pressure
k.lo('Pct')=10;
k.up('Pct')=100;
*boundary on density
k.lo('Rhou')= pL('Rhou');
k.up ('Rhou')= pU('Rhou');
Tbr.lo=0.001;
Tbr.up=1;
Betta.lo=-16000;
Betta.up=-1e-10;

*Intial guesses
k.l('Vt')= 50;
Tavgr.l = 0.01;
k.l('Tct') = 273;
k.l('Pct') = 10;
k.l('Vap') = 0.01;
k.l('Tmp') = 280;
n.l(i)=5;
k.l('Rhou') = 1.0;
k.l('eta_h') = 2.0;
k.l('MWt')= 350;
k.l('Tbp') = 394;
*Ra.l= 0.1;
k.l('RED') = 0.4;
Tavgr.l= 0.8;
Tau.l= 0.5;
nv(i,x) = 0;

OPTION SYSOUT = ON
option MIP=cplex;
option MINLP=BARON;
option optcr=0;
option optca=1E-5;

Model modelPM /all/;
modelPM.OPTFILE=1;

*Integer cut
nv(i,x)=0;
ALIAS(x,xx);
LOOP(xx,
         SOLVE modelPM MINIMIZE OF USING MINLP;
         nv(i,xx)=n.l(i);
         OFv(xx)=OF.l;
         kv(p,xx) = k.l(p);
         dyn(xx)=YES;
);
Display nv,OFv, kv;