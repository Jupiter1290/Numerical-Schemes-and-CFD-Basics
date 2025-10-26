function [rho,v,T,mfr,rho_th,v_th,T_th]=cons(rho,v,T,a,x,nt,dt,n)
dx=x(2)-x(1);
rho_old=rho;
v_old=v;
T_old=T;
gamma=1.4;
%Solution vectors:
u1=rho.*a; u2=rho.*a.*v; u3=rho.*(T/(gamma-1)+(gamma/2)*v.^2).*a;
for k=1:nt
%Flux vectors:
f1=rho.*a.*v; f2=rho.*a.*v.^2+(1/gamma)*rho.*T.*a;
f3=rho.*(T/(gamma-1)+(gamma/2)*v.^2).*a.*v+rho.*T.*a.*v;
%Source Terms:
for l=2:n-1
    dadx(l)=(a(l+1)-a(l))/dx;
    j2(l)=(1/gamma)*rho(l)*T(l)*dadx(l);
end

u1_i=u1; u2_i=u2; u3_i=u3;
for i=2:n-1
%Predictor Step:
du1dt_p(i)=-(f1(i+1)-f1(i))/dx;
du2dt_p(i)=-(f2(i+1)-f2(i))/(dx)+j2(i);
du3dt_p(i)=-(f3(i+1)-f3(i))/dx;
%Updating the values:
u1(i)=u1(i)+du1dt_p(i)*dt;
u2(i)=u2(i)+du2dt_p(i)*dt;
u3(i)=u3(i)+du3dt_p(i)*dt;
end

%Updating rho,v,T:
rho=u1./a; v=u2./u1; T=(gamma-1)*(u3./u1-(gamma/2)*v.^2);
%Flux vectors:
f1=rho.*a.*v; f2=rho.*a.*v.^2+(1/gamma)*rho.*T.*a;
f3=rho.*(T/(gamma-1)+(gamma/2)*v.^2).*a.*v+rho.*T.*a.*v;
%Source Terms:
for l=2:n-1
    dadx(l)=(a(l)-a(l-1))/dx;
    j2(l)=(1/gamma)*rho(l)*T(l)*dadx(l);
end

%Predictor Step, averaged derivative and updating the final values:
for i=2:n-1
du1dt_c(i)=-(f1(i)-f1(i-1))/dx;
du2dt_c(i)=-(f2(i)-f2(i-1))/(dx)+j2(i);
du3dt_c(i)=-(f3(i)-f3(i-1))/dx;

du1dt_av(i)=(du1dt_p(i)+du1dt_c(i))/2;
du2dt_av(i)=(du2dt_p(i)+du2dt_c(i))/2;
du3dt_av(i)=(du3dt_p(i)+du3dt_c(i))/2;

u1(i)=u1_i(i)+du1dt_av(i)*dt;
u2(i)=u2_i(i)+du2dt_av(i)*dt;
u3(i)=u3_i(i)+du3dt_av(i)*dt;
end
%Boundary Conditions:
u1(1)=rho(1)*a(1);
u2(1)=2*u2(2)-u2(3);
u3(1)=u1(1)*(T(1)/(gamma-1)+gamma/(2)*v(1)^2);

u1(n)=2*(u1(n-1))-(u1(n-2));
u2(n)=2*(u2(n-1))-(u2(n-2));
u3(n)=2*(u3(n-1))-(u3(n-2));
%plot(x,v)
%pause(0.1)
%Evaluating the primitive variables:
rho=u1./a; v=u2./u1; T=(gamma-1)*(u3./u1-(gamma/2)*v.^2);
mfr(k)=rho(16)*a(16)*v(16);
v_th(k)=v(16);
rho_th(k)=rho(16);
T_th(k)=T(16);
end
end










    