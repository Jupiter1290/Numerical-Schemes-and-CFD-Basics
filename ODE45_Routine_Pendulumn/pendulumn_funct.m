function [angdot]= pendulumn_funct(t,ang,b,m,l)
g=9.81;
angdot=zeros(2,1);
%The rate of change of angular velocity=angular acceleration
angdot(1)=ang(2);
angdot(2)=-((b/m)*ang(2)+(g/l)*sin(ang(1)));
end