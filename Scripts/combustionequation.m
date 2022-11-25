% THIS SCRIPT IS FORMED TO CALCULATE THE MOLE NUMBERS (RELATIVE & REAL) OF
% BOTH REACTANTS AND PRODUCTS 
a1f=1.13711;a2f=1.4553*10^-2;a3f=-2.95876*10^-6;a4f=0;
a1o2=3.25304;a2o2=6.5235*10^-4;a3o2=-1.49524*10^-7;a4o2=1.5389*10^-11;
a1n2=3.3443;a2n2=2.9426*10^-4;a3n2=1.953*10^-9;a4n2=-6.574*10^-12;
a1co2=3.09590;a2co2=2.73114*10^-3;a3co2=-7.8854*10^-7;a4co2=8.66002*10^-11;
a1h2o=3.74292;a2h2o=5.65590*10^-4;a3h2o=4.95240*10^-8;a4h2o=-1.81802*10^-11;
a1co=3.317;a2co=3.76970*10^-4;a3co=-3.22080*10^-8;a4co=-2.19450*10^-12;

%___________given values______________________
n_cum=3;            %CnHm of fuel
m_cum=8;
lamda_cum=1;        %amout of excees air of more air
%___________constant calculated vlues_________
X_air=(n_cum+(m_cum/4));    %amount of theortical air
b_h2o=m_cum/2;          %amount of H2O
%___________changing values___________________
a_co2=0;            %CO2 constant
e_o2=0;            %excess O2 constant 
d_co=0;            %CO constant in case of incomblite comdustion
N_r_r=0;        %relative mole number of the reactants
N_p_r=0;        %relative mole number of the products
N_r_re=0;       %real mole number of the reactants
N_p_re=0;       %real mole number of the products
xf=0;           %fuel mole fraction
xo2_r=0;        %oxcgen mole fraction of the reactants
xn2_r=0;        %nitrogen mole fraction of the reactants
xa=0;           %CO2 mole fraction of the products
xb=0;           %H2O mole fraction of the products
xd=0;           %CO mole fraction of the products
xe=0;           %oxcgen mole fraction of the products
xn2_p=0;        %nitrogen mole fraction of the products
%___________calculating other varaibles_______

if(lamda_cum==1)
    %this is a theortical mixture (e=d=0)
    e_o2=0;d_co=0;
    a_co2=n_cum;
elseif(lamda_cum>1)
    %this is excess air combustion(d=0)
    d_co=0;
    a_co2=n_cum;
    e_o2=(X_air*lamda_cum)-a_co2-(b_h2o/2);
else
    %this is an incomblite combustion(e=0)
    e_o2=0;
    d_co=2*(n_cum+(b_h2o/2)-(X_air*lamda_cum));
    a_co2=n_cum-d_co;
end
N_r_r=1+(4.76*(X_air*lamda_cum));
N_p_r=a_co2+b_h2o+d_co+e_o2+(3.76*(X_air*lamda_cum));
xf=1/N_r_r;
xo2_r=(X_air*lamda_cum)/N_r_r;
xn2_r=(X_air*lamda_cum*3.76)/N_r_r;
xa=a_co2/N_p_r;
xb=b_h2o/N_p_r;
xd=d_co/N_p_r;
xe=e_o2/N_p_r;
xn2_p=(3.76*X_air*lamda_cum)/N_p_r;
disp('X=');
disp(X_air);
disp('Lamda=');
disp(lamda_cum);
disp('a=');
disp(a_co2);
disp('b=');
disp(b_h2o);
disp('e=');
disp(e_o2);
disp('d=');
disp(d_co);
disp('relative mole number of reactants=');
disp(N_r_r);
disp('relative mole number of products=');
disp(N_p_r);

a_equ=((a1f*xf)+(a1o2*xo2_r)+(xn2_r*a1n2))
b_equ=((a2f*xf)+(a2o2*xo2_r)+(a2n2*xn2_r))
c_equ=((a3f*xf)+(a3o2*xo2_r)+(a3n2*xn2_r))
d_equ=((a4f*xf)+(a4o2*xo2_r)+(a4n2*xn2_r))

a_t=((a1f*xf)+(a1o2*xo2_r)+(xn2_r*a1n2));
a_tr=((a_t*N_r_r)-1);
b_tr=((a2f*xf)+(a2o2*xo2_r)+(a2n2*xn2_r))*N_r_r;
c_tr=((a3f*xf)+(a3o2*xo2_r)+(a3n2*xn2_r))*N_r_r;
d_tr=((a4f*xf)+(a4o2*xo2_r)+(a4n2*xn2_r))*N_r_r;

ap_tr=((a1co2*xa)+(a1h2o*xb)+(a1n2*xn2_p));
bp_tr=((a2co2*xa)+(a2h2o*xb)+(a2n2*xn2_p));
cp_tr=((a3co2*xa)+(a3h2o*xb)+(a3n2*xn2_p));
dp_tr=((a4co2*xa)+(a4h2o*xb)+(a4n2*xn2_p));

e0_r=N_r_r*xf*-9-0.90510*10^8;
e0_p=(N_p_r*xa*-3.9364*10^8)+(N_p_r*xb*-2.39225*10^8);

disp('Xf=');
disp(xf);
disp('xo2_r=');
disp(xo2_r);
disp('xn2_r=');
disp(xn2_r);
disp('xco2=');
disp(xa);
disp('xh2o=');
disp(xb);
disp('xn2 pro=');
disp(xn2_p);
e_r=(8314*((a_tr*300)+(b_tr*300^2)+(c_tr*300^3)+(d_tr*300^4)-300));
e_r_t=e_r-e0_r;
z=(e_r_t+e0_p)/8314;


