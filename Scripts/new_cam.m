% clear all;clc;
max_area=0;
max_vlift=0;
p_area=0.0058; %in m^2
% d_theta=0.1;
%____________________assumpations_____________________________________%
rb_cam=10;%radious of base circle
rn_cam=8; %radious of nosie circle
r_cam=8;  %radious of cam follower
dp_cam=30; %radious of intake port
%____________________martsises________________________________________%
a_cam=1:d_theta:720;%this matrix is uesed to store the theta value in degree
b_cam=1:d_theta:720;%this matrix is uesed to store the theta value in radians
% l=1:d_theta:720;
ll_cam=zeros(1,(length(a_cam)));% this matrix to store the valve lift
aiv_cam=zeros(1,(length(a_cam)));%this matrix to store the instantions valve area
av_cam=zeros(1,(length(a_cam)));%this matrix to store the intake port area with respect to valve lift
test_b_cam_length=length(b_cam)
test_ll_cam_cam_length=length(ll_cam)
test_value_1=(89*(1/d_theta))
test_value_2=(90*(1/d_theta))
test_value_3=(132*(1/d_theta))
test_value_4=round(45*(1/d_theta))
test_case=zeros(1,(length(a_cam)));
%__________________helpful variables__________________________________%
gama_cam=deg2rad(110); %maximum cam angle 
D_cam=((rb_cam-rn_cam)/(cos(gama_cam)));
d2_cam=D_cam^2;
xr_cam=(r_cam+rn_cam)^2;
Ap_cam=((pi/4)*(dp_cam^2))*10^-6 % intake port area(constant)
for i=1:test_value_2
b_cam(i)=deg2rad(a_cam(i));
% l(i)=1/(cos(b(i)));
% ll(i)=((r+rb)*(l(i)-1));
ll_cam(i)=(sin(b_cam(i)))*10;
aiv_cam(i)=((pi*dp_cam*(ll_cam(i)))-(i*2*(d_theta)))*10^-6;
if(aiv_cam(i)< aiv_cam((test_value_4)))
   aiv_cam(i)=aiv_cam(i-1); 
end
% av_cam(i)=min(Ap_cam,aiv_cam(i));
av_cam(i)=(min(Ap_cam,aiv_cam(i)))/(5);
test_case(i)=av_cam(i)/p_area;
end
% ll(91)=15;
for i=test_value_2:test_value_3
        th=deg2rad(110-a_cam(i));
        sinth=sin(th);
        costh=cos(th);
        ll_cam(i)=((xr_cam-(d2_cam*(sinth^2)))^0.5)+(D_cam*costh)-(r_cam+rb_cam)-8;
        ll_cam(i)=(1/ll_cam(i))+7.365+2.764;
        aiv_cam(i)=(pi*dp_cam*ll_cam(i))*10^-6;
%         av_cam(i)=min(Ap_cam,aiv_cam(i));
        av_cam(i)=(min(Ap_cam,aiv_cam(i)))/(5);
        test_case(i)=av_cam(i)/p_area;
end
% ll(133)=15;
for i=1:test_value_1
      ll_cam(test_value_3+i)=ll_cam(test_value_2-i);
      aiv_cam(test_value_3+i)=(pi*dp_cam*ll_cam(test_value_2-i));
%       av_cam(test_value_3+i)=min(Ap_cam,aiv_cam(test_value_2-i));
      av_cam(test_value_3+i)=(min(Ap_cam,aiv_cam(test_value_2-i)))/(5);
      test_case(i)=av_cam(i)/p_area;
end
max_area=max(av_cam)
max_vlift=max(ll_cam)
figure
subplot(2,1,1);
 plot(a_cam,ll_cam,'g','LineWidth',3);
 title('valve lift VS cam angle');
 xlabel('cam angle');
 ylabel('valve lift in (mm)');
 subplot(2,1,2);
 plot(a_cam,av_cam,'r','LineWidth',3);
 title('intake area VS cam angle');
 xlabel('cam angle');
 ylabel('intake area in m^2');
 figure
 plot(a_cam,test_case,'g','LineWidth',3);
 title('valve area / piston area VS cam angle');
 xlabel('cam angle');
 ylabel('valve area / piston area');
 av_cam(1)