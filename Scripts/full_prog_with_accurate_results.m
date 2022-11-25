clc;clear;
%_____________________Constants__________________________
cd=0.6;                                    %cd of inlet area of intake valve
ro_air_s=1.18;                             %standard air density
P_atm=101325;                              %standard atmospheric pressure
t_in=320;                                  %at the beginging the inlet tempreture will be 320K
d_theta=0.1;
%_____________________scripts calls____________________
piston_subroutine;         %call this to get the instantanous volume
new_cam;                  %call this to get the instantanous valve area
combustionequation;       %call this to get all the coumbustion prcentages
matrix_end=length(a_crank);
volume_length=length(g_crank);
%__________________________Variables and arays to store data________________
m_new=zeros(1,matrix_end);           %this variable to store the instantanous intake air mass
m_dot=zeros(1,matrix_end);           %this variable to store the instantanous intake mass flow rate 
mole_test=zeros(1,matrix_end);
P_cy=zeros(1,matrix_end);
T_cy=zeros(1,matrix_end);
et_nolumetric=zeros(1,matrix_end);
testing_mass=zeros(1,matrix_end);
%__________intial values___________________________
P0=100000;       %at the start we will assume intial pressure of 10 pascal
N=2000;           % the engine runs at 2000rpm
omega=(2*pi*N)/60;
v0=7.136530000000000e-05;
x0=320;                        %this X is a symbol fot the cylinder temperature
eps=10;
itter=0;
Ru=8314;                          %universal gas constant in J/Kmole.K
N_r=((P0*v0)/(Ru*x0));
N_r_re=N_r;
P_cy(1)=P0;
xn=0;
count_test=0;
T_CYLINDER=0;
P_CYLINDER=P0;
xf;
E0= N_r_re*xf *(-0.905*10^8);
E_prev=19.04;
MW_rec=(xf*44)+(xo2_r*32)+(xn2_r*28);
v_max=max(g_crank)+v0;
th_mass=ro_air_s*v_max;
iter=0;
H_in_newton=0;
work_newton=0;
Cp=1;
dt=(d_theta/(6*N));
m_new(1)=((P0*v0)/(287*x0));
volume_length=length(g_crank);
 fprintf('____________Beginning of the new code_______________\n')

%___________________________________SOLVING______________________
for i=1:(719*(1/d_theta))
    fprintf('__________loop number %0.8f______________________\n',i)
    m_dot(i)=(cd*(av_cam(i))*((2*ro_air_s)^0.5))*(((P_atm-P_cy(i))/((abs(P_atm-P_cy(i)))^0.5)));
    m_new(i+1)=m_new(i)+((dt)*m_dot(i));
   testing=m_new(i+1);
   N_r_re=m_new(i+1)/MW_rec;
    E0
    T_CYLINDER=x0
    P_CYLINDER
    E_prev
    dv=g_crank(i)-v0
    work=((P_CYLINDER)*dv)
    H_dot=(m_dot(i))*(Cp)*(t_in)*(1000)
    H=H_dot*dt;
    e1=@(x)((a_equ-1)*x);
    e2=@(x)((b_equ)*x^2);
    e3=@(x)((c_equ)*x^3);                            
    e4=@(x)((d_equ)*x^4);
    etot_prev=@(x)(Ru*N_r_re*((e1(x))+(e2(x))+(e3(x))+(e4(x))))+(E0);
    E_equ=((-E0+E_prev+H-work)/(N_r_re*Ru));
    etot=@(x)((e1(x))+((e2(x))+((e3(x))+((e4(x))+(-E_equ)))));
    de1=@(x)((a_equ-1));
    de2=@(x)(2*(b_equ)*x);
    de3=@(x)(3*(c_equ)*x^2);
    de4=@(x)(4*(d_equ)*x^3);
    detot=@(x)((de1(x))+((de2(x))+((de3(x))+((de4(x))))));
    F=etot(x0);
    dF=detot(x0);
    xn=x0-(F/dF);
    tol=abs(xn-x0);
    x0=xn;
while tol > (0.001)
    F=etot(x0);
    dF=detot(x0);
    xn=x0-(F/dF);
    tol=abs(xn-x0);
    x0=xn;
    x0;
    iter=iter+1;
end
iter
T_cy(i)=xn;
fprintf('The temperature of this loop is %0.8f\n',T_cy(i));
mole_test(i)=N_r_re;
v0=g_crank(i);
P_cy(i)=((N_r_re * Ru *xn)/(g_crank(i)));
P_CYLINDER=P_cy(i);
mole_test(i)=N_r_re;
et_nolumetric(i)=m_new(i+1)/th_mass;
x0=xn; 
E_prev=etot_prev(x0);
E0= N_r_re*xf *(-0.905*10^8);    
testing_mass(i)=(m_dot(i)*(1/dt))/(th_mass);   
end
et_nolumetric(matrix_end)=et_nolumetric(matrix_end-1);
mole_test(matrix_end)=mole_test(matrix_end-1);
T_cy(matrix_end)=T_cy(matrix_end-1);
P_cy(matrix_end)=P_cy(matrix_end-1);
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
title('Crank angle VS mdot with mtotal')
xlabel('crank angle');
ylabel('volumetric efficiency');