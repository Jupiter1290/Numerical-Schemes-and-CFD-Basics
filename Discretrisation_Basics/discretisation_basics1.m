clear
close all
clc
x=pi/3; dx=linspace(pi/4000,pi/40,20);
%given function: f(x)=sin(x)/x^3;
%f'(x)=(cos(x)*x^3)-3*x^2*sin(x))/x^6 (analytically computed);
analytical_derivative=(cos(x)*x^3-3*x^2*sin(x))/x^6;

for k=1:20
%forward_differencing(first order approximation):forward_differencing=(f(x+dx)-f(x))/dx
numerical_derivative(1,k)=((sin(x+dx(k))/(x+dx(k))^3)-(sin(x)/x^3))/(dx(k));
%central differnecing(second order approximation) central_difference=(f(x+dx)-f(x-dx))/2*dx;
numerical_derivative(2,k)=((sin(x+dx(k))/(x+dx(k))^3)-(sin(x-dx(k))/(x-dx(k))^3))/(2*dx(k));
%fourth_order_approximation:fourth_order_approx=(f(x-2*dx)-8*f(x-dx)+8*f(x+dx)-f(x+2*dx))/12*dx;
numerical_derivative(3,k)=((sin(x-2*dx(k))/(x-2*dx(k))^3)-8*(sin(x-dx(k))/(x-dx(k))^3)+8*(sin(x+dx(k))/(x+dx(k))^3)-(sin(x+2*dx(k))/(x+2*dx(k))^3))/(12*dx(k));
err=abs(numerical_derivative-analytical_derivative);
end
%plots for "error-vs-dx"
figure(1)
bar(err)
set(gca,'XTickLabel', {'1st Order','2nd Order','4th Order'})
ylabel('error')
title('Comparison of Numerical Errors in Derivative Approximations')
grid on
%calculating slopes;
slope_firstO=diff(err(1,:))./diff(dx);
slope_seconderO=diff(err(2,:))./diff(dx);
slope_thirdO=diff(err(3,:))./diff(dx);
%loglog;
figure(2)
loglog(dx,err(1,:),'r')
title('Log–Log Plot: Error vs dx')
xlabel('dx')
ylabel('error')
hold on
loglog(dx,err(2,:),'g')
xlabel('dx')
ylabel('error')
hold on
loglog(dx,err(3,:),'b')
xlabel('dx')
ylabel('error')
legend('Red->First Order','Green->Second Order','Blue->Fourth Order')


