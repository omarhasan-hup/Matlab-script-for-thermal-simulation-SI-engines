clc;clear;
% syms x;      %declare the symbolic functions with variable x
%_____________________componant constants______________
% a1f=1.13711;a2f=1.4553*10^-2;a3f=-2.95876*10^-6;a4f=0;
% a1o2=3.25304;a2o2=6.5235*10^-4;a3o2=-1.49524*10^-7;a4o2=1.5389-10^-11;
% a1n2=3.3443;a2n2=2.9426*10^-4;a3n2=1.953*10^-9;a4n2=-6.574*10^-12;
% a1co2=3.09590;a2co2=2.73114*10^-3;a3co2=-7.8854*10^-7;a4co2=8.66002*10^-11;
% a1h2o=3.74292;a2h2o=5.65590*10^-4;a3h2o=4.95240*10^-8;a4h2o=-1.81802*10^-11;
% a1co=3.317;a2co=3.76970*10^-4;a3co=-3.22080*10^-8;a4co=-2.19450*10^-12;
%_____________________Constants__________________________
cd=0.6;                                    %cd of inlet area of intake valve
ro_air_s=1.18;                             %standard air density
P_atm=101325;                              %standard atmospheric pressure
d_theta=0.5;
%_____________________scripts calls____________________
piston_subroutine;                       %call this to get the instantanous volume
new_cam;                                       %call this to get the instantanous valve area
combustionequation;       %call this to get all the coumbustion prcentages
matrix_end=length(a_crank)
%__________________________Variables and arays to store data________________
m_new=zeros(1,matrix_end);                            %this variable to store the instantanous intake air mass
m_dot=zeros(1,matrix_end);                        %this variable to store the instantanous intake mass flow rate 
mole_test=zeros(1,matrix_end);
P_cy=zeros(1,matrix_end);
T_cy=zeros(1,matrix_end);
H_in=zeros(1,matrix_end);
df_H_in=zeros(1,matrix_end);
work=zeros(1,matrix_end);
df_work=zeros(1,matrix_end);
E_r_t=zeros(1,matrix_end);                      %this vector to store the total reactant energy
df_E_r_t=zeros(1,matrix_end);  
e1=zeros(1,matrix_end);                         % for the first term in E_r_t
df_e1=zeros(1,matrix_end);
e2=zeros(1,matrix_end);                         % for the second term in E_r_t
df_e2=zeros(1,matrix_end);
e3=zeros(1,matrix_end);                         % for the third term in E_r_t
df_e3=zeros(1,matrix_end);
e4=zeros(1,matrix_end);                         % for the fourth term in E_r_t
df_e4=zeros(1,matrix_end);
itter_count=zeros(1,matrix_end);
et_nolumetric=zeros(1,matrix_end);
dv_test_cse=zeros(1,matrix_end);
%__________intial values___________________________
P0=100000;       %at the start we will assume intial pressure of 10 pascal
t0=500;           %at the beginging the inlet tempreture will be 320K
% T_cy(1)=390;
N=2000;           % the engine runs at 2000rpm
omega=(2*pi*N)/60;
%m0=0;            % assume there was no mass in the cylinder at the start
v0=7.136530000000000e-05;
eps=10;
itter=0;
Ru=8314         %universal gas constant in J/mole.K
dv=7.136530000000000e-05;           % delta V at the start is equal to the clearance volume
N_r=2.899161664*10^-6 %number of real moles =(p*v)/(Ru*t)
N_r_re=N_r;
P_cy(1)=P0;
dE=0;            % delta internal enegry
sss=(N_r_re*Ru);
x0=300;                        %this X is a symbol fot the cylinder temperature
xn=0;
count_test=0;
ii=0;
xf
E0= N_r_re*xf *(-0.905*10^8)
E_prev=0;
dE_prev=0;
E_prev_n=0;
dE_prev_n=0;
df_E=0;
MW_rec=(xf*44000)+(xo2_r*32000)+(xn2_r*28000);
v_max=max(g_crank);
th_mass=1.16*v_max;
% a_equ=0;
% b_equ=0;
% c_equ=0;
% d_equ=0;
Cp=1.2;
%___________________________________SOLVING______________________
for i=1:719
%     ss=((N_r_re*Ru)/(g_crank(i)));
sss=(N_r_re*Ru);
sss_check=sss;
%  H_in(i)=(276.989*(av_cam(i))*(P_atm-P_cy(i))^0.5) /1000;
 H_in(i)=(276.989*(av_cam(i))*((P_atm-P_cy(i))^0.5 ))*1000;
 work(i)=(dv*6*N*(P_cy(i)));
 test_case_h=H_in(1);
 test_case_w=work(1);

 e1=@(x)(x * ((a_equ*N_r_re)-1));
 e2=@(x)(x^2 *N_r_re*b_equ);
 e3=@(x)(x^3 *N_r_re*c_equ);
 e4=@(x)(x^4 *N_r_re*d_equ);
 ee=@(x)(e1(x)+e2(x)+e3(x)+e4(x));
 E_r_t=@(x)(6*N)*((Ru * ee(x))-E0);
 dE=E_r_t(x0)-E_prev;
 F= dE + work(i)-H_in(i);
 energy= E_r_t(x0);
 energy_difference = dE;
 total_equ= F;
 
%  E_prev=E_r_t(i);
%  df_H_in(i)=276.989*(av_cam(i))*(P_atm-P_cy(i))^0.5 *1000;    %this term
%  will be zero after the differentiaition 
%  df_work(i)=dv*6*N*(ss);
 df_e1=@(x)(((a_equ*N_r_re)-1));
 df_e2=@(x)(2*x*b_equ*N_r_re);
 df_e3=@(x)(3*x^2*c_equ*N_r_re);
 df_e4=@(x)(4*x^3*d_equ*N_r_re);
 df_ee=@(x)(df_e1(x)+df_e2(x)+df_e3(x)+df_e4(x));
 df_E_r_t=@(x) (6*N)*((Ru * df_ee(x)));
 df_E=df_E_r_t(x0)-dE_prev;
 energy_diff=df_E_r_t(x0);
 dF=df_E;
 dF; 
 
%  dF=inline(dF);
 xn=x0-(F/dF);
%  xn=abs(xn);
 tol=abs(xn-x0);
tol

 ii=i;
count_test=count_test+1;

% the beggining of newton raphson method
e_iterate=zeros(1,matrix_end);
de_iterate=zeros(1,matrix_end);
e_iterate(1)=xn;
de_iterate(1)=xn;
j=1;
x0=xn;
e_new_real=zeros(1,50);
e_prev_real=zeros(1,50);
de_new_real=zeros(1,50);
de_prev_real=zeros(1,50);
count_iterate=zeros(1,50);
ii=0;
%  for omar=1:1000
while tol > (0.08)
   itter =itter+1
   i;
   ii=ii+1;
   count_iterate(ii)=ii;
%    x0=abs(x0);
   %________________________________
%  H_in(ii)=(276.989*(av_cam(ii))*((P_atm-P_cy(ii))^0.5 ))/1000;
 H_in(ii)=(276.989*(av_cam(ii))*((P_atm-P_cy(ii))^0.5 ))*1000;
 work(ii)=dv*6*N*(P_cy(ii));
 e1=@(x)(x * ((a_equ*N_r_re)-1));
 e2=@(x)(x^2 *N_r_re*b_equ);
 e3=@(x)(x^3 *N_r_re*c_equ);
 e4=@(x)(x^4 *N_r_re*d_equ);
 ee=@(x)(e1(x)+e2(x)+e3(x)+e4(x));
 E_r_t=@(x)(6*N)*((Ru * ee(x))-E0);
 dE=E_r_t(x0)-E_prev_n;
%  e_new_real(omar)=E_r_t(x0);
%  e_prev_real(omar)=E_prev_n;
 E_prev_n=E_r_t(x0);
 e_previous=E_prev_n;     %test
 e_now=E_r_t(x0);
 d_energy=dE;
 F= dE+ work(ii)-H_in(ii);
 total_equ=F;
%  F=inline(F);
 df_e1=@(x)(((a_equ*N_r_re)-1));
 df_e2=@(x)(2*x*b_equ*N_r_re);
 df_e3=@(x)(3*x^2*c_equ*N_r_re);
 df_e4=@(x)(4*x^3*d_equ*N_r_re);
 df_ee=@(x)(df_e1(x)+df_e2(x)+df_e3(x)+df_e4(x));
 df_E_r_t=@(x) ((6*N))*((Ru * df_ee(x)));
 df_E=df_E_r_t(x0)-dE_prev_n;
%  de_new_real(omar)=df_E_r_t(x0);
%  de_prev_real(omar)=dE_prev_n;
de_previous=dE_prev_n;
de_now=df_E_r_t(x0);
dd_energy=df_E;

 dF=df_E;
 dE_prev_n=df_E_r_t(x0);
%  E_prev_n=E_r_t(x0);
 %_______________________________________
   xn=x0-(F/dF);
%    xn=abs(xn);
   tol=abs(xn-x0);
   fprintf('In this iteration the temp is  %0.8f\n',xn)
   fprintf('de_prev_n is  %0.8f\n',dE_prev_n)
   fprintf('E_prev_n is  %0.8f\n',E_prev_n)
   x0=xn;
 end
 %End of newton raphson method
%  figure
% plot(count_iterate,e_new_real,'g','LineWidth',3);
% title('iteration VS now energy')
% xlabel('iteration');
% ylabel('now energy');
%  figure
% plot(count_iterate,e_prev_real,'g','LineWidth',3);
% title('iteration VS prev energy')
% xlabel('iteration');
% ylabel('prev energy');
itter_count(i)=itter;
itter
xn


T_cy(i)=xn;
% P_cy(i)=((N_r_re * Ru *T_cy(i))/(g_crank(i)));
% ppp=P_cy(i);
mole_test(i)=N_r_re;
dv=(g_crank(i+1)-g_crank(i));
% m_dot(i)=0.6*av_cam(i)*1.53622915*(P_atm-P_cy(i))^0.5;
% m_dot(i)=(0.6*av_cam(i)*1.56524)*((abs(P_atm-P_cy(i))/((abs(P_atm-P_cy(i)))^0.5)));
m_dot(i)=(0.6*av_cam(i)*1.56524)*(((P_atm-P_cy(i))/((abs(P_atm-P_cy(i)))^0.5)));
m_new(i+1)=m_new(i)+((1/6*N)*m_dot(i));
N_r_re=m_new(i+1)/MW_rec;
P_cy(i)=((N_r_re * Ru *T_cy(i))/(g_crank(i)));
mole_test(i)=N_r_re;
% test_cse(i)=sss;
et_nolumetric(i)=m_new(i+1)/th_mass;
E_prev = E_r_t(xn);
df_E = df_E_r_t(xn);
x0=xn;    
E0= N_r_re*xf *(-0.905*10^8);    
    
      
end

count_test;
figure
plot(a_crank,P_cy,'g','LineWidth',3);
title('Crank angle VS Cylinder pressure')
xlabel('crank angle');
ylabel('cylinder pressure in PASCAL');
figure
plot(a_crank,T_cy,'r','LineWidth',3);
title('Crank angle VS Cylinder Temp')
xlabel('crank angle');
ylabel('cylinder Temp in K');
figure
plot(mole_test,'b','LineWidth',3);
title('Crank angle VS real mole numbrs')
xlabel('crank angle');
ylabel('real mole number');
figure
plot(a_crank,m_new,'c','LineWidth',3);
title('Crank angle VS instantanous mass in cylinder')
xlabel('crank angle');
ylabel('instantanous mass in cylinder');
figure
plot(a_crank,et_nolumetric,'c','LineWidth',3);
title('Crank angle VS volumetric efficiency')
xlabel('crank angle');
ylabel('volumetric efficiency');
% max(et_nolumetric)
% figure
% plot(a_crank,E_r_t,'c','LineWidth',3);
% title('Crank angle VS internal energy')
% xlabel('crank angle');
% ylabel('internal energy');

% a_equ=((a1f*xf)+(a1o2*xo2_r)+(xn2_r*a1n2)-(xf+xo2_r+xn2_r))
% b_equ=((a2f*xf)+(a2o2*xo2_r)+(a2n2*xn2_r))
% c_equ=((a3f*xf)+(a3o2*xo2_r)+(a3n2*xn2_r))
% d_equ=((a4f*xf)+(a4o2*xo2_r)+(a4n2*xn2_r))
















