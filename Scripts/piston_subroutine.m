% THIS SCRIPT IS USED TO CONTRACT A RELATION BETWEEN THE CRANK ANGLE AND
% THE PISTON POSITION
% d_theta=1;                             % the incremental step of theta 
bo_crank=86;                             % assuming the engine is square engin so stroke= bore
Vc_crank=7.13653e-05;                    % clearance volume in m^3
pa_crank=(pi/4)*((bo_crank*10^-3)^2)    %this is the piston area in m^2
a_crank=[1:d_theta:720];                           %theta in degree
b_crank=[1:length(a_crank)];                           % theta in radians
c_crank=[1:length(a_crank)];                           %piston position
lc_crank=215;                            % conecting rod length
rc_crank=43;                             %crank radious       
sum=lc_crank+rc_crank;                         %sumation of conecting rod lengt and the crand radious
d_crank=[1:length(a_crank)];                           % cos(crank angle)
e_crank=[1:length(a_crank)];                           % (sin(crank angle))^2
f_crank=[1:length(a_crank)];                           
g_crank=[1:length(a_crank)];                           %this vector is to store the clinder volume
for i=1:length(a_crank)
    b_crank(i)=deg2rad(a_crank(i));      % to change the crank angle from degrees to rad 
    d_crank(i)=(cos(b_crank(i)));
    e_crank(i)=((sin(b_crank(i)))^2);
    f_crank(i)=sqrt((lc_crank^2)-((rc_crank^2)*e_crank(i)));
    c_crank(i)=sum-(rc_crank*(d_crank(i)))-f_crank(i);
    g_crank(i)=((c_crank(i)*10^-3*pa_crank))+Vc_crank;
end
test=min(g_crank)
figure
plot(a_crank,c_crank,'linewidth',3);
title('Piston position VS Crank angle');
xlabel('Crank angle (degree)');
ylabel('piston position (mm)');
figure
plot(a_crank,g_crank,'linewidth',3)
title('Stroke volume VS Crank angle');
xlabel('Crank angle (degree)');
ylabel('Stroke volume (m^3)');
v_max=max(g_crank)
th_mass=1.16*v_max
% g_crank(1)
% dv1=g_crank(1)-Vc_crank
% dv_test=g_crank(2)-g_crank(1)
% pdv=(10^5)*dv1