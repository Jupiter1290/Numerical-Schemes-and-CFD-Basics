function [rho,v,T,mfr,rho_th,v_th,T_th]=non_cons2(rho,v,T,a,x,nt,dt,n)
dx=x(2)-x(1);
gamma=1.4;

for k=1:nt
rho_old=rho;
v_old=v;
T_old=T;
%rho_i=rho;
%v_i=v;
%T_i=T;
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
rho(i)=rho_old(i)+drhodt_p(i)*dt;
v(i)=v_old(i)+dvdt_p(i)*dt;
T(i)=T_old(i)+dTdt_p(i)*dt;
end

for i=2:n-1
%Corrector step (using the updated values):
dvdx(i)=(v(i)-v(i-1))/(dx);
dlnadx(i)=(log(a(i))-log(a(i-1)))/(dx);
drhodx(i)=(rho(i)-rho(i-1))/(dx);
dTdx(i)=(T(i)-T(i-1))/(dx);

drhodt_c(i)=-rho(i)*dvdx(i)-rho(i)*v(i)*dlnadx(i)-v(i)*drhodx(i);
dvdt_c(i)=-v(i)*dvdx(i)-(1/gamma)*(dTdx(i)+(T(i)/rho(i))*drhodx(i));
dTdx_c(i)=-v(i)*dTdx(i)-(gamma-1)*T(i)*(dvdx(i)+v(i)*dlnadx(i));
end
 
%Average Time Derivative that will be used for the final values of the
%time-step;
drhodt_av=(drhodt_p+drhodt_c)/2;
dvdt_av=(dvdt_p+dvdt_c)/2;
dTdt_av=(dTdt_p+dTdx_c);

%Final Values:
for i=2:n-1
rho(i)=rho_old(i)+drhodt_av(i)*dt;
v(i)=v_old(i)+dvdt_av(i)*dt;
T(i)=T_old(i)+dTdt_av(i)*dt;
end

%Boundary Conditions:
v(1) = 2*v(2)-v(3);

v(n) = 2*v(n-1)-v(n-2);
rho(n) = 2*rho(n-1)-rho(n-2);
T(n) = 2*T(n-1)-T(n-2);
%Mass flow rate:
mfr(k)=rho(16)*a(16)*v(16);
v_th(k)=v(16);
rho_th(k)=rho(16);
T_th(k)=T(16);
end

end