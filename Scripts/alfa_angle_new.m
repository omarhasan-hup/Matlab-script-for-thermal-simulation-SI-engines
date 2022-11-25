d_theta=0.1;
%________________assumptions__________________
theta_ign=300;                 %crank angle begining of igntion
delta_theta_ign=20;            %delata crank angle of igntion delay
boc=theta_ign+delta_theta_ign; %crank angle beginging of combustion
delta_theta_comb=60;           %delata crank angle of combustion
eoc=boc+delta_theta_comb;      %crank angle end of combustion
%______________usful variables__________________




%____________variables to store data_________________
theta_crank=[1:d_theta:720]; % this matrix is used to store the values of crank angle
a=zeros(1,length(theta_crank)); % this matrix is used to store the alfa angle value 
% theta_crank_length=length(theta_crank)
% a_length=length(a)

j=(boc*(1/d_theta));
%__________________main program_______________________
for i=(boc*(1/d_theta)):(eoc*(1/d_theta))
   b=deg2rad(((abs(theta_crank(j)-boc))/delta_theta_comb)*90);
   a(i)=sin(b);
   j=j+1;
end

for i=((eoc)*(1/d_theta)):((720*(1/d_theta))-(1/d_theta))
   a(i)=1; 
end
a(a_length)=1;
% theta_crank_length=length(theta_crank)
% a_length=length(a)

plot(theta_crank,a,'r','LineWidth',4);
title('alfa_combustion VS crank angle');
xlabel('Crank_angle');
ylabel('Alfa');
