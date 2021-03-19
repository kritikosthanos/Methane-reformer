# Methane-reformer

function h=enthalpy(y,T) 

clc 

% A user defined function to calculate 

% the specific enthalpy of the gases of 

% the SR/SOFC/UB systems 

% 

% Input arguments : y = vector of mole fractions of the gas mixture [-] 

%                   T = temperature [K] 

% Output argument : h = specific enthalpy of gas mixture [kJ/kmol] 

 

% Components 

%     1        2      3      4    5      6      7 

%     CH4      H2O    CO     H2   CO2    O2     N2 

 

h=0; 

href=[-74870   -241830 -110530 0 -393520 0 0]; % enthalpy of formation of the ideal das at 298.15 K in [kJ/kmol] 

Tref=298.15;                                   % temperature at which enthapy equals enthalpy of formation of ig 

% cp = C + B * T + A * T^2 

% constants A,B,C are fitted values to raw data from NIST Chemistry webbook 

% 

A=[0           0         2.5359e-6     1.5807e-6   -2.2051e-05  -7.1534e-06    0      ]; 

B=[0.05346     0.01133   0.00300006    -0.00050303  0.05228      0.016494    0.0053629]; 

C=[19.676      29.704    27.774        29.058       23.853       24.991         27.237]; 

for j=1:7 

    h= h + y(j) * ((A(j)/3) * (T^3-Tref^3) + (B(j)/2) * (T^2-Tref^2) + C(j) * (T-Tref) + href(j)); 

end 

end 

 

Κώδικας αντιδραστήρα: 

clear all 

clc 

M=652.3/3600;%kg/s 

p=0.3915;%kg/m3 

kg=7.715e-2;%W/m*K 

vi=2.373e-5;%kg/m*s 

 

dt=[0.10226 0.127 0.1524];%m 

dp=dt/17;%m  

Nt=[6 7 8 9 10]; 

%    i=2; 

%    j=4; 

%%%%%%%% CH4=5.98880110637327e-004  %%%%%%%%%%%%% 

% Components 

%     1        2      3      4    5      6      7 

%     CH4      H2O    CO     H2   CO2    O2     N2 

%revma 1+2 

ni=[0.002414841330850 0.007244523992000 0 0 0 0 0.000426148470150]; 

Ni=sum(ni);%kmol/s 

yi=ni/Ni; 

%revma 7 

no=[0 0.226404287975066 0 0 4.52808686679159e-002 7.97242856909162e-002 0.648590557666102]; 

Nout=sum(no);%kmol/s 

yo=no/Nout; 

 

%Temperatures-pressure 

Ti=773.15;%C->K 

Tout=1839.15;% K  

 

P=1.4;%bar 

n=0; 

 

%Length 

zi=0;zf=5; 

z=linspace(zi,zf,10000); 

zz=zeros(length(z),length(dt)*length(Nt)*length(Tout)); 

s=zeros(length(z),11,length(dt)*length(Nt)*length(Tout)); 

%Mass matrix 

MASS=eye(11); MASS(8,8)=0; MASS(9,9)=0; 

 

 

 for i=1:length(dt) 

       for j=1:length(Nt) 

          for k=1:length(Tout) 

n=n+1 

 

At(n)=Nt(j)*((pi*(dt(i))^2)/4);%m^2 

ut(n)=M/(p*At(n));%m/s 

Re(n)=(dp(i)*ut(n)*p)/vi; 

e(n)=0.38+0.073*(1+(1-(2*dp(i))/dt(i))^2); 

f(n)=((1-e(n))/e(n)^3)*(1.75 + 150*((1-e(n))/Re(n))); 

 

dpdz(n)=(f(n)*p*(ut(n))^2)/(dp(i)*100*1000);%bar m- 

hi(n)=(kg*0.813*(Re(n)^0.9)*exp(-(6*dp(i))/dt(i)))/dt(i);%W/m^2 K 

U(n)=hi(n)/1000;%kJ/m^2 K s 

%epipleon 

pcat=2380;%kg/m3 

pbed(n)=pcat*(1-e(n));%kg/m3 

A(n)=At(n)*pbed(n);%kg/m 

B(n)=Nt(j)*pi*dt(i)*U(n);%kJ/m K s 

 

%systhma diaforikwn-algevrikwn-arxikes 

h=Ni*enthalpy(yi,Ti);%kJ/s 

ho(n)=Nout*enthalpy(yo,Tout(k));%kJ/s 

 

%arxikes 

x0=[]; 

x0(1:7)=ni;%kmol/s 

x0(8)=Ti;%K 

x0(9)=Tout(k);%K 

x0(10)=h;%kJ/s 

x0(11)=ho(n);%kJ/s 

 

%epilush systhmatos  

options =odeset('reltol',1e-4,'abstol',1e-15,'mass',MASS); 

         

[z,x]=ode15s(@reformer_sunartisi,z,x0,options,A(n),B(n),no,P); 

s(1:length(x(:,1)),:,n)=x; 

zz(1:length(z),n)=z; 

           end 

       end 

  end 

%Vcat=At*z*(1-e); 

 

 

Κώδικας επίλυσης ode: 

function f = reformer_sunartisi(z,x,A,B,no,P,dpdz) 

 

 

ni=x(1:7);%kmol/s 

Tr=x(8); 

Tb=x(9); 

H=x(10); 

Ho=x(11); 

% Pi=[]; 

% Pi=P-dpdz 

y=ni/sum(ni); 

srr=srreaction(P,y,Tr); 

 

f(1)=-A*srr(1)/3600;%CH4 kmol/s m 

f(2)=+A*(-srr(1)-srr(2))/3600;%H2O kmol/s m 

f(3)=+A*(srr(1)-srr(2))/3600;%CO kmol/s m 

f(4)=+A*(3*srr(1)+srr(2))/3600;%H2 kmol/s m 

f(5)=+A*srr(2)/3600;%CO2 kmol/s m 

f(6)=0;%O2 kmol/s m 

f(7)=0;%N2 kmol/s m 

f(10)=B*(Tb-Tr);%kJ/s m 

f(11)=-B*(Tb-Tr);% kJ/s m 

%algebraic 

f(8)=H-enthalpy(ni,Tr); 

f(9)=Ho-enthalpy(no,Tb); 

f=f(:); 

 

end 

 

Κώδικας srreactor: 

function srr=srreaction(P,y,T) 

clc 

% A user defined function to calculate the 

% reaction rates of the methane steam reformer 

% according to the paper of Xu&Froment  

% 

% CH4 + H2O <-> CO  + 3 H2   reaction 1 

% CO  + H2O <-> CO2 +   H2   reaction 2 

% 

% RETURNS 

% 

% r : reaction rates  [kmol/(h kg_cat)] 

% 

% INPUT ARGUMENTS 

% 

% P : pressure        [bar] 

% y : molar fractions [-] 

% T : temperature     [K] 

 

%components 

%      1    2    3   4   5   6   7 

%      CH4  H2O  CO  H2  CO2 O2  N2 

 

% gas constant 

R=8.314; %            [kJ/(kmol K)] 

 

% partial pressures   [bar] 

p=P*y; 

p(4)=max(p(4),eps); 

% denominator 

K_CH4 = 6.65E-4*exp(+38280/(R*T));% [bar^-1] 

K_CO  = 8.23E-5*exp(+70650/(R*T));% [bar^-1] 

K_H2  = 6.12E-9*exp(+82900/(R*T));% [bar^-1] 

K_H2O = 1.77E+5*exp(-88680/(R*T));% [-] 

DEN   = 1+K_CH4*p(1)+K_CO*p(3)+K_H2*p(4)+K_H2O*p(2)/max(p(4),eps); % [-] 

 

% equilibrium 

K1=7.85e12*exp(-220200/(R*T)); % [bar^2]  

% K2=exp(4400/T -4.036); 

K2=0.017668*exp(36581.6/(R*T)); % [-] 

 

% reaction rate constants 

  

k1=4.225e+15*exp(-240100/(R*T)); % [kmol bar^0.5/(h kgcat)] 

k2=1.955e+06*exp( -67130/(R*T)); % [kmol/(h kgcat bar)] 

 

% reaction rates 

 

srr(1) = ( k1 / p(4)^2.5 ) * ( p(1)*p(2) - p(4)^3*p(3)/K1 ) / DEN^2; 

srr(2) = ( k2 / p(4)     ) * ( p(3)*p(2) - p(4)  *p(5)/K2 ) / DEN^2; 

 
