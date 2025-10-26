clear
close all
clc
x=pi/3;
dx=linspace(pi/3000,pi/3,20);
%given function: f(x)=exp(x)*cos(x);
%analytical_derivative=-exp(x)*(cos(x)+sin(x))
analytical_derivative=-2*exp(x)*sin(x);
for i=1:20
%central_difference=(-0.0833f(x-2*dx)+1.3333*f(x-dx)-2.5*f(x)+1.3333*f(x+dx)-0.0833f(x+2*dx))/dx^2
%right_skewed=(3.75*f(x)-12.8333*f(x+dx)+17.8333f(x+2*dx)-13*f(x+3*dx)+5.0833*f(x+4*dx)-0.8333f(x+5*dx))/dx^2
%left_skewed=(3.75*f(x)-12.8333*f(x-dx)+17.8333f(x-2*dx)-13*f(x-3*dx)+5.0833*f(x-4*dx)-0.8333f(x-5*dx))/dx^2
%Coefficients for central_difference:
A=[1 1 1 1 1;-2 -1 0 1 2;2 1/2 0 1/2 2;-4/3 -1/6 0 1/6 4/3;2/3 1/24 0 1/24 2/3];
B=[0 0 1 0 0 ]';
%By method of matrix multiplication:
X=inv(A)*B;
%Coefficients for right_skewed
A1=[1 1 1 1 1 1;0 1 2 3 4 5;0 1/2 2 9/2 8 25/2;0 1/6 4/3 9/2 32/3 125/6;0 1/24 2/3 27/8 32/3 625/24;0 1/120 4/15 81/40 128/15 625/24];
B1=[0 0 1 0 0 0]';
%By method of matrix multiplication:
X1=inv(A1)*B1;
%Coefficients for left_skewed
A2=[1 1 1 1 1 1;0 -1 -2 -3 -4 -5;0 1/2 2 9/2 8 25/2;0 -1/6 -4/3 -9/2 -32/3 -125/6;0 1/24 2/3 27/8 32/3 625/24;0 -1/120 -4/15 -81/40 -128/15 -625/24];
B2=[0 0 1 0 0 0]';
%By method of matrix multiplication:
X2=inv(A2)*B2;
central_difference(i)=(X(1)*exp(x-2*dx(i))*cos(x-2*dx(i))+X(2)*exp(x-dx(i))*cos(x-dx(i))+X(3)*exp(x)*cos(x)+X(4)*exp(x+dx(i))*cos(x+dx(i))+X(5)*exp(x+2*dx(i))*cos(x+2*dx(i)))/dx(i)^2;
right_skewed(i)=(X1(1)*exp(x)*cos(x)+X1(2)*exp(x+dx(i))*cos(x+dx(i))+X1(3)*exp(x+2*dx(i))*cos(x+2*dx(i))+X1(4)*exp(x+3*dx(i))*cos(x+3*dx(i))+X1(5)*exp(x+4*dx(i))*cos(x+4*dx(i))+X1(6)*exp(x+5*dx(i))*cos(x+5*dx(i)))/dx(i)^2;
left_skewed(i)=(X2(1)*exp(x)*cos(x)+X2(2)*exp(x-dx(i))*cos(x-dx(i))+X2(3)*exp(x-2*dx(i))*cos(x-2*dx(i))+X2(4)*exp(x-3*dx(i))*cos(x-3*dx(i))+X2(5)*exp(x-4*dx(i))*cos(x-4*dx(i))+X2(6)*exp(x-5*dx(i))*cos(x-5*dx(i)))/dx(i)^2;
end
error_central_difference=abs(central_difference-analytical_derivative);
error_right_skewed=abs(right_skewed-analytical_derivative);
error_left_skewed=abs(left_skewed-analytical_derivative);
figure(1)
title('loglog plot for error-vs-dx')
loglog(dx,error_central_difference,'r')
hold on
loglog(dx,error_right_skewed,'g')
hold on
loglog(dx,error_left_skewed,'b')
legend('Red-> Central Difference','Green->Right Skewed','Blue->Left Skewed')
xlabel('dx');
ylabel('absolute error')
