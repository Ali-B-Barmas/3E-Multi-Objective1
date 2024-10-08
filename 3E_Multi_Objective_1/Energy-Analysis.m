clear all
close all
clc

%% inputs 

Ta=298.15;                              % [K]    (Ambient Temperature)
T0=298.15;                              % [K]    (Dead State Temperature)
P0=101.3;                               % [kPa]  (Dead State Pressure)
Tsun=5739;                              % [k]    (Temperature Of Sun)

% PTSC
eta_ps=0.90;                            % Pump Isentropic Efficiency
fluid_S='T66';
Di=0.08;                                % Reciever Pipe Inlet Diameter [m]
Do=0.09;                                % Reciever Pipe Outlet Diameter [m]
Dc=0.15;                                % Glass Cover Main Diameter [m]
L=350;                                  % Collector Length [m]
epsc=0.87;                              % Emissivity Of Glass Cover
epsr=0.92;                              % Emissivity Of Reciever
S=900*0.64;                             % Solar Radiation [W/m^2]
Tin=422;                                % [K]
Tout=600;                               % [K]
zigma=5.67*10^(-8);                     % [J/(s.m^2.K^4)]
K_air=0.025;                            % [N/m.K]
K_f=0.1;                                % [N/m.K]
w=6;                                    % Width Of Collector [m]

% Brayton Cycle
fluid='CO2';
eta_tb=0.93;                            % Turbine Isentropic Efficiency
eta_cb=0.90;                            % Compressor Isentropic Efficiency
T04=32+273.15;                          % Compressor Inlet Temperature [k]
P04=8000;                               % Compressor Inlet Pressure [kPa]
P05=20000;                              % Compressor Outlet Pressure [kPa]
P06=P05;                         
P07=8000;                               % Turbine Outlet Pressure [kPa]
P08=8000;                          
mb=2.1;                                 % Mass Flow Rate Of CO2 [kg/s]
 
% Rankine Cycle
eta_tr=0.93;                            % Turbine Isentropic Efficiency
eta_pr=0.90;                            % Pump Isentropic Efficiency
P09=6130;                               % Pump Inlet Pressure [kPa]
P10=9000;                               % Pump Outlet Pressure [kPa]
P11=9000;                         
P12=6130;                               % Pressure Of Condenser [kPa]
mr=1.7;                                 % Mass Flow Rate Of CO2 [kg/s]

% Vapor Compression Refrigeration 
eta_cv=0.95;                            % Compressor Isentropic Efficiency
P13=3500;                               % Compressor Inlet Pressure [kPa]
P14=6130;                               % Compressor Outlet Pressure [kPa]
P15=6130;                               % Pressure Of Condenser [kPa]
P16=3500;                               % Pressure Of Evaporator [kPa]
mv=1;                                   % Mass Flow Rate Of CO2 [kg/s]

%% PTSC

Tr_av=(Tin+Tout)/2;         
Aap=(w-Dc)*L;                           % Unshaded Area
Ar=pi*(Di+Do)*L/2;                 
Ac=pi*(Dc)*L/2;
    
% Nusselt Number
Pr=0.74;                                % Prandtl Number
V=4;                                    % air velocity [m/s]
v=1.5*10^(-5);                          % Kinematic viscosity [m2/s]
Re= V*Dc/v;                             % Reynolds Number
Nu_f= 0.023*(Re^0.8)*(Pr^0.4);          % Nusselt Number Of Fluid
Nu_a=0.4;                               % Air Reynolds -> 0
hc_ga= Nu_a*K_air/Dc;                   % Convection (Glass & Air)
hc_f= Nu_f*K_f/Di;                      % Convection (Reciever)
    
Tg=500;                                 % [K]
 for n=1:5   

hr_rg=(zigma*(Tg+Tr_av)*(Tg^2+Tr_av^2))/((1/epsr)+(Ar/Ac)*((1/epsc)-1));
hr_ga=epsc*zigma*(Tg+Ta)*(Tg^2+Ta^2);
Tg=((hr_rg*Tr_av*Ar)+(Ac*(hc_ga+hr_ga)*Ta))/((hr_rg*Ar)+(Ac*(hc_ga+hr_ga)));

 end

UL= (Ar/((hc_ga+hr_ga)*Ac)+(1/hr_rg))^(-1);
U0= ((1/UL)+(Do/(hc_f*Di))+((Do/(2*K_f))*log(Do/Di)))^(-1);
F1= U0/UL;
    
% Calculating Q_solar & m_dot_solar (by guessing m_dot_solar )

P02=140;                               % [kPa]
P03=P02*0.90;                          % Pressure Drop Of 5% in PTSC [kPa]
P01=P03*0.95;                        
P_bar=(P03+P02)/2;
T_bar=(Tout+Tin)/2;
% specific heat of Therminal 66 at 504.3 K
C_T66= Props('C','P',P_bar,'T',T_bar,fluid_S)*1000;

ms=6;                                  % [kg/s]
for i=1:5
    FR= ((ms*C_T66)/(Ar*UL))*(1-exp((-Ar*UL*F1)/(ms*C_T66)));
    Q_solar= FR*((S*Aap)-Ar*UL*(Tin-Ta));                
    ms= Q_solar/(C_T66*(Tout-Tin));
end  
    ms;                                % [kg/s]
    Q_solar=Q_solar/1000;              % [kj/s]
%% Solar Cycle

% Pump I (state 1 & 2)
T01=Tin; 
T03=Tout;
h01= Props('H','P',P03,'T',T01,fluid_S);       % [kj/kg]
s01= Props('S','P',P03,'T',T01,fluid_S);       % [kj/(kg.K)]
D01= Props('D','P',P03,'T',T01,fluid_S);       % [kg/m^3]
v01= 1/D01;                                    % [m^3/Kg]
w_p1_s= v01 * (P02-P01);                       % [kj/kg]
w_p1_a= w_p1_s / eta_ps;                       % [kj/kg]
h02= h01 + w_p1_a;                             % [kj/kg]
s02= Props('S','P',P02,'H',h02,fluid_S);       % [kj/(kg.K)]
T02= Props('T','P',P02,'H',h02,fluid_S);       % [K]

% Solar Heat (state 3)

h03= h02 + (Q_solar/ms);                         % [kj/kg]
s03= Props('S','P',P03,'H',h03,fluid_S);       % [kj/(kg.K)]

%% Bryton Cycle

% Compressor I (state 4 & 5)k
h04= Props('H','T',T04,'P',P04,fluid);          % [kj/kg]
s04= Props('S','T',T04,'P',P04,fluid);          % [kj/kg.K]
s5s=s04;              
h5s= Props('H','P',P05,'S',s5s,fluid);          % [kj/kg]
w_c1_s = h5s - h04;                             % [kj/kg]
w_c1_a = w_c1_s / eta_cb;                       % [kj/kg]
h05= h04 + w_c1_a;                              % [kj/kg]
s05= Props('S','P',P05,'H',h05,fluid);          % [kj/kg.K]
T05= Props('T','P',P05,'H',h05,fluid);          % [K]

% Heat Exchanger (state 6)
h06= (ms*(h03-h01) + (mb*h05))/(mb);            % [kj/kg]
s06= Props('S','P',P06,'H',h06,fluid);          % [kj/kg.K]
T06= Props('T','P',P06,'H',h06,fluid);          % [K]
 
% Turbine I (state 7)
s7s=s06;              
h7s= Props('H','P',P07,'S',s06,fluid);          % [kj/kg]
w_t1_s= h06 - h7s;                              % [kj/kg]
w_t1_a= w_t1_s * eta_tb;                        % [kj/kg]
h07= h06 - w_t1_a;                              % [kj/kg]
s07= Props('S','P',P07,'H',h07,fluid);          % [kj/kg.K]

%% Vapor Compression Cycle

% Compressor (state 13 & state 14)
h13= Props('H','P',P13,'Q',1,fluid);            % [kj/kg]
s13= Props('S','P',P13,'Q',1,fluid);            % [kj/kg.K]
T13= Props('T','P',P13,'Q',1,fluid);            % [K]
s14s=s13;              
h14s= Props('H','P',P14,'S',s14s,fluid);        % [kj/kg]
w_c2_s = h14s - h13;                            % [kj/kg]
w_c2_a = w_c2_s / eta_cv;                       % [kj/kg]
h14= h13 + w_c2_a;                              % [kj/kg]
s14= Props('S','P',P14,'H',h14,fluid);          % [kj/kg.K]
T14= Props('T','P',P14,'H',h14,fluid);          % [K]


% Condenser (state 15)
h15= Props('H','P',P15,'Q',0,fluid);            % [kj/kg]
s15= Props('S','P',P15,'Q',0,fluid);            % [kj/kg.K]

% Expansion Valve (state 16)
h16 = h15;
s16 = Props('S','P',P16,'H',h16,fluid);         % [kj/kg.K]

%% Rankine Cycle

% Turbine (state 11 & 12)
w_t2_a= mv * w_c2_a / mr;
w_t2_s= w_t2_a / eta_tr;
  
h11=440;                                        % [kj/kg]
for i=1:5      
    h12s= h11 - w_t2_s;
    s12s= Props('S','P',P12,'H',h12s,fluid);
    s11=s12s;
    h11= Props('H','P',P11,'S',s11,fluid);
end

h12= h11 - w_t2_a;                              % [kj/kg]
T12= Props('T','P',P12,'H',h12,fluid);          % [K]
s12= Props('S','P',P12,'H',h12,fluid);     



% Condenser (state 9)
h09= Props('H','P',P09,'Q',0,fluid);           % [kj/kg]
s09= Props('S','P',P09,'Q',0,fluid);           % [kj/kg.K]
D09= Props('D','P',P09,'Q',0,fluid);           % [kg/m^3]
v09= 1/D09;                                    % [m^3/Kg]
T09= Props('T','P',P09,'Q',0,fluid);           % [K]


% Pump II (state 10)
w_p2_s = v09 * (P10-P09);                     % [kj/kg]
w_p2_a = w_p2_s / eta_pr;                     % [kj/kg]
h10= h09 + w_p2_a;                            % [kj/kg]
s10= Props('S','P',P10,'H',h10,fluid);        % [kj/(kg.K)]

% Heat Exchanger (state 8)
h08= (mb*h07 + mr*h10 - mr*h11)/mb;
s08= Props('S','P',P08,'H',h08,fluid);        % [kj/(kg.K)]
T08= Props('T','P',P08,'H',h08,fluid);        % [K]
%% TEG I 
T20=T08-10;
T27=T0;
P20=P0;
P27=P0; 
TL=0.5*(T20+T27);
TH=0.5*(T08+T04);
TM=0.5*(TL+TH);
ZTM=0.9;

eta_TEG1=(1-TL/TH)*((1+ZTM)^0.5-1)/((1+ZTM)^0.5+TL/TH);

h20=Props('H','T',T20,'P',P20,fluid);
s20=Props('S','P',P20,'T',T20,fluid);
h27=Props('H','T',T27,'P',P27,fluid);
s27=Props('S','P',P27,'T',T27,fluid);

mcool1=(mb*(h08-h04))/((1+eta_TEG1)*(h20-h27));
Qin1=mb*(h08-h04);
Q_ELEGANT1=mcool1*(h27-h20);

W_TEG1=eta_TEG1*Q_ELEGANT1;

%% TEG II 
T22=T09-10;
T23=T0;
P22=P0;
P23=P0; 
TL=0.5*(T22+T23);
TH=0.5*(T12+T09);
TM=0.5*(TL+TH);
ZTM=0.9;

eta_TEG2=(1-TL/TH)*((1+ZTM)^0.5-1)/((1+ZTM)^0.5+TL/TH);

h22=Props('H','T',T22,'P',P22,fluid);
s22=Props('S','P',P22,'T',T22,fluid);
h23=Props('H','T',T23,'P',P23,fluid);
s23=Props('S','P',P23,'T',T23,fluid);

mcool2=(mr*(h12-h09)+mv*(h14-h15))/((1+eta_TEG1)*(h23-h22));
Qin2=mr*(h12-h09)+mv*(h14-h15);
Q_ELEGANT2=mcool2*(h23-h22);

W_TEG2=eta_TEG2*Q_ELEGANT2;
%% Performance

W_turbine= mb*w_t1_a;                         % [kj/s]
W_compressor= mb*w_c1_a;                      % [kj/s]
W_pump1= ms * w_p1_a;                         % [kj/s]
W_pump2= mr*w_p2_a;                           % [kj/s]

W_net= W_turbine + W_TEG1 + W_TEG2 - W_compressor - W_pump1 - W_pump2;

Q_cond= mr*(h12-h09)+mv*(h14-h15);
Q_gc= mb*(h08-h04);
Q_eva= mv*(h13-h16);

eta_sys= (W_net/Q_solar);

%% Exergy Analysis

h0= Props('H','P',P0,'T',T0,fluid);
s0= Props('S','P',P0,'T',T0,fluid);

h0_S=Props('H','P',P0,'T',T0,fluid_S);
s0_S=Props('S','P',P0,'T',T0,fluid_S);

ex01= h01 - h0_S -T0*(s01 - s0_S); 
Ex01= ms*ex01;                                 % kw

ex02= h02 - h0_S -T0*(s02 - s0_S); 
Ex02= ms*ex02;                                 % kw

ex03= h03 - h0_S -T0*(s03 - s0_S); 
Ex03= ms*ex03;                                 % kw

ex04= h04 - h0 -T0*(s04 - s0); 
Ex04= mb*ex04;                                 % kw

ex05= h05 - h0 -T0*(s05 - s0); 
Ex05= mb*ex05;                                 % kw

ex06= h06 - h0 -T0*(s06 - s0); 
Ex06= mb*ex06;                                 % kw

ex07= h07 - h0 -T0*(s07 - s0); 
Ex07= mb*ex07;                                 % kw

ex08= h08 - h0 -T0*(s08 - s0); 
Ex08= mb*ex08;                                 % kw

ex09= h09 - h0 -T0*(s09 - s0); 
Ex09= mr*ex09;                                 % kw

ex10= h10 - h0 -T0*(s10 - s0); 
Ex10= mr*ex10;                                 % kw

ex11= h11 - h0 -T0*(s11 - s0); 
Ex11= mr*ex11;                                 % kw

ex12= h12 - h0 -T0*(s12 - s0); 
Ex12= mr*ex12;                                 % kw

ex13= h13 - h0 -T0*(s13 - s0); 
Ex13= mv*ex13;                                 % kw

ex14= h14 - h0 -T0*(s14 - s0); 
Ex14= mv*ex14;                                 % kw

ex15= h15 - h0 -T0*(s15 - s0); 
Ex15= mv*ex15;                                 % kw

ex16= h16 - h0 -T0*(s16 - s0); 
Ex16= mv*ex16;                                 % kw

ex22= h22 - h0 - T0*(s22-s0);
Ex22=mcool2*ex22;

ex20= h20 - h0 - T0*(s20-s0);
Ex20=mcool1*ex20;

%% Exergy Detruction

% Pump I
Ex_des_p1= ms*T0*(s02-s01); 
eta_ex_p1= 1-(Ex_des_p1/ms*w_p1_a);

%PTSC
Exsun= ((0.9*900*Aap*(1+(T0/Tsun)^4-((4*T0)/(3*Tsun))))/1000)-300;
Ex_des_PTSC= 0.64*Exsun -(Ex03 - Ex01)+150;
eta_ex_PTSC= 1-(Ex_des_PTSC/Exsun);

% BC Heat Exchanger
Ex_des_BCHEX=T0*(ms*(s01-s03)+mb*(s06-s05));
eta_ex_BCHEX=1-(Ex_des_BCHEX/(Ex03-Ex01));

% Compressor I
Ex_des_c1= mb*T0*(s05-s04); 
eta_ex_c1= 1-(Ex_des_c1/W_compressor);

% Turbine I
Ex_des_t1= mb*T0*(s07-s06);
eta_ex_t1= 1-(Ex_des_t1/(Ex06-Ex07));

% Gas Cooler
Ex_des_gc= -T0*(mb*(s04-s08)+ mcool1*(s20-s0));
eta_ex_gc= 1-(Ex_des_gc/(Ex08-Ex04));

% Rc Heat Exchanger
Ex_des_RCHEX= T0*(mb*(s08-s07)+mr*(s11-s10));
eta_ex_RCHEX= 1-(Ex_des_BCHEX/(Ex07-Ex07));

% Turbine II
Ex_des_t2= mr*T0*(s12-s11);
eta_ex_t2= 1-(Ex_des_t2/(Ex11-Ex12));

% Pump II
Ex_des_p2= T0*mr*(s10-s09); 
eta_ex_p2= 1-(Ex_des_p2/mr*w_p2_a);

% Compressor II
Ex_des_c2= T0*mv*(s14-s13); 
eta_ex_c2= 1-(Ex_des_c2/(mv*w_c2_a));

% Evaporator
Ex_des_eva=mv*T0*(s13-s16-(Q_eva/288));
eta_ex_eva=1-(Ex_des_eva/(Ex16-Ex13));

% Condenser
Ex_des_cond= T0*((mr*(s09-s12))+(mv*(s15-s14))+Q_cond/295);
eta_ex_cond= 1-(Ex_des_cond/(Ex12+Ex14-Ex09-Ex15));

% Expansion Valve
Ex_des_ev= T0*mv*(s16-s15);
eta_ex_ev=1-Ex_des_ev/Ex15;

% TEG 1
Ex_des_TEG1= T0*(mcool1*(s20-s0))/100;
Ex_des_TEG2= T0*(mcool2*(s0-s22))/100;

% Exergy Performance 
Ex_des_total= Ex_des_p1 + Ex_des_p2 + Ex_des_BCHEX + Ex_des_t1 + Ex_des_gc+...
    +Ex_des_c1+Ex_des_t2+Ex_des_RCHEX+Ex_des_cond+...
    Ex_des_eva+Ex_des_c2+Ex_des_ev+Ex_des_PTSC+Ex_des_TEG1+Ex_des_TEG2;
eta_ex_total = 1 - ( (Ex_des_total)/ (Exsun) );




