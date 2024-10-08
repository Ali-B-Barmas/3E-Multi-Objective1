%% Economic analysis
i = 0.12;                                    % Constant
n=20;                                        % year
CRF= (i* ( 1+i )^n)/ ( (1+i)^n-1 );
phi= 1.06;                                   % Constant
N=8000;                                      % hours in 1 year

%PTSC
    Aap=(6-0.15)*350;
    Z_PTSC= 354*(Aap);   
    Z_dot_PTSC= (Z_PTSC1*CRF*phi)/(N*3600)*60*60;                     %$/h
    
%Pump I
    W_Pump1=
    Z_Pump1= 3540*(W_Pump1);   
    Z_dot_Pump1= (Z_Pump1*CRF*phi)/(N*3600)*60*60;                    %$/h
    
%Pump II
    W_pump2=
    Z_Pump2= 3540*(W_pump2);   
    Z_dot_Pump2= (Z_Pump2*CRF*phi)/(N*3600)*60*60;                    %$/h
    
%BC Heat Exchanger
    A_BCHEX=
    EZ_BCHEX= 4.6656-0.1557*log(A_BCHEX)+0.1547*(log(A_BCHEX))^2;
    Z_BCHEX= 10^EZ_BCHEX;
    Z_dot_BCHEX= (Z_BCHEX*CRF*phi)/(N*3600)*60*60;                    %$/h
    
%BC Compressor
    W_Comp1=
    Z_Comp1= 10167*((W_Comp1)^0.46);   
    Z_dot_Comp1= (Z_Comp1*CRF*phi)/(N*3600)*60*60;                    %$/h
    
%Gas Cooler
    A_GC=
    EZ_GC= 4.6656-0.1557*log(A_GC)+0.1547*(log(A_GC))^2;
    Z_GC= 10^EZ_GC;
    Z_dot_GC= (Z_GC*CRF*phi)/(N*3600)*60*60;                          %$/h    
   
%BC Turbine
    W_T1=
    EZ_T1= 2.6259-1.4398*log(W_T1)-0.1776*(log(W_T1))^2;
    Z_T1= 10^EZ_T1;
    Z_dot_T1= (Z_T1*CRF*phi)/(N*3600)*60*60;                          %$/h
    
%RC Heat Exchanger
    A_RCHEX=
    EZ_RCHEX= 4.6656-0.1557*log(A_RCHEX)+0.1547*(log(A_RCHEX))^2;
    Z_RCHEX= 10^EZ_RCHEX;
    Z_dot_RCHEX= (Z_RCHEX*CRF*phi)/(N*3600)*60*60;                    %$/h
    
%RC Turbine
    W_T2=
    EZ_T2= 2.6259-1.4398*log(W_T2)-0.1776*(log(W_T2))^2;
    Z_T2= 10^EZ_T2;
    Z_dot_T2= (Z_T2*CRF*phi)/(N*3600)*60*60;                          %$/h
    
%Condenser
    mr=
    mb=
    Z_Cond= 1773*(mr+mb);   
    Z_dot_Cond= (Z_Cond*CRF*phi)/(N*3600)*60*60;                      %$/h
    
%VC Compressor
    W_Comp2=
    Z_Comp2= 10167*((W_Comp2)^0.46);   
    Z_dot_Comp2= (Z_Comp2*CRF*phi)/(N*3600)*60*60;                    %$/h
    
%Evaporator
    A_Ev=
    EZ_Ev= 4.6656-0.1557*log(A_Ev)+0.1547*(log(A_Ev))^2;
    Z_Ev= 10^EZ_Ev;
    Z_dot_Ev= (Z_Ev*CRF*phi)/(N*3600)*60*60;                          %$/h
    
%TEG I & II
    W_TEG1=
    W_TEG2=
    Z_TEG1=2000*W_TEG1;
    Z_TEG2=2000*W_TEGII;
    Z_dot_TEG1= (Z_TEG1*CRF*phi)/(N*3600)*60*60;                      %$/h
    Z_dot_TEG2= (Z_TEG2*CRF*phi)/(N*3600)*60*60;                      %$/h

Z_dot_total=Z_dot_PTSC+Z_dot_Pump1+Z_dot_BCHEX+Z_dot_Comp1+Z_dot_T1+...
    +Z_dot_GC+Z_dot_T2+Z_dot_RCHEX+Z_dot_Cond+Z_dot_Pump2+Z_dot_Comp2+...
    +Z_dot_Ev+Z_dot_TEG1+Z_dot_TEG2;