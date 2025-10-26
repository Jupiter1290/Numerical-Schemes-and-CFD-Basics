function [rho,v,T]=non_cons(rho,v,T,a,x,nt,dt,n)
dx=x(2)-x(1);
rho_old=rho;
v_old=v;
T_old=T;
gamma=1.4;

for k=1:nt
rho_i=rho;
v_i=v;
T_i=T;
for i=2:n-1
%Predictor Step:
dvdx(i)=(v(i+1)-v(i))/(dx);
dlnadx(i)=(log(a(i+1))-log(a(i)))/(dx);
drhodx(i)=(rho(i+1)-rho(i))/(dx);
dTdx(i)=(T(i+1)-T(i))/(dx);

drhodt_p(i)=-rho(i)*dvdx(i)-rho(i)*v(i)*dlnadx(i)-v(i)*drhodx(i);
dvdt_p(i)=-v(i)*dvdx(i)-(1/gamma)*(dTdx(i)+(T(i)/rho(i))*drhodx(i));
dTdt_p(i)=-v(i)*dTdx(i)-(gamma-1)*T(i)*(dvdx(i)+v(i)*dlnadx(i));

%Updating the values:
rho_u(i)=rho_old(i)+drhodt_p(i)*dt;
v_u(i)=v_old(i)+dvdt_p(i)*dt;
T_u(i)=T_old(i)+dTdt_p(i)*dt;

%Corrector step (using the updated values):
dvdx_u(i)=(v(i)-v(i-1))/(dx);
dlnadx_u(i)=(log(a(i))-log(a(i-1)))/(dx);
drhodx_u(i)=(rho(i)-rho(i-1))/(dx);
dTdx_u(i)=(T(i)-T(i-1))/(dx);

drhodt_c(i)=-rho_u(i)*dvdx_u(i)-rho_u(i)*v_u(i)*dlnadx_u(i)-v_u(i)*drhodx_u(i);
dvdt_c(i)=-v_u(i)*dvdx_u(i)-(1/gamma)*(dTdx_u(i)+(T(i)/rho_u(i))*drhodx_u(i));
dTdx_c(i)=-v_u(i)*dTdx_u(i)-(gamma-1)*T_u(i)*(dvdx_u(i)+v(i)*dlnadx_u(i));

%Average Time Derivative that will be used for the final values of the
%time-step;
drhodt_av(i)=(drhodt_p(i)+drhodt_c(i))/2;
dvdt_av(i)=(dvdt_p(i)+dvdt_c(i))/2;
dTdt_av(i)=(dTdt_p(i)+dTdx_c(i));

%Final Values for the time-step:
rho(i)=rho_i(i)+drhodt_av(i)*dt;
v(i)=v_i(i)+dvdt_av(i)*dt;
T(i)=T_i(i)+dTdt_av(i)*dt;
end

%Boundary Conditions:
v(1) = 2*v(2)-v(3);

v(n) = 2*v(n-1)-v(n-2);
rho(n) = 2*rho(n-1)-rho(n-2);
T(n) = 2*T(n-1)-T(n-2);
end
end

