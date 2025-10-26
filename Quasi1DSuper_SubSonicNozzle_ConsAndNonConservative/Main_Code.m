close all
clear 
clc
%defining the 1-D mesh:
L=3; n=31; x=linspace(0,3,n); dx=x(2)-x(1); nt=1400;gamma=1.4;R=8.3145;
%Initial Profiles:
rho=1-0.3146*x; rho_i=rho;
T=1-0.2314*x; T_i=T;
v=(0.1+1.09*x).*T.*0.5; v_i=v;
a=1+2.2*(x-1.5).^2;
a_s=(gamma*R*T).^(1/2);
rho_old=rho;
v_old=v;
T_old=T;
for i=1:n
    if x(i)<=0.5
        rho_p(i)=1;
        T_p(i)=1;
    elseif x(i)>=0.5 && x(i)<=1.5
        rho_p(i)=1-0.366*(x(i)-0.5);
        T_p(i)=1-0.167*(x(i)-0.5);
    elseif x(i)>=1.5 && x(i)<=3.5
        rho_p(i)=0.634-0.3879*(x(i)-1.5);
        T_p(i)=0.833-0.3507*(x(i)-1.5);
    end
end
v_p=0.59*(rho_p.*a).^(-1);
%Incorporating the CFL criteria:
C=0.85;
del_t=((a_s+v).^-1)*dx*0.8;
dt=min(del_t);
%Calling the functions:
[rho_nc2, v_nc2, T_nc2,mfr_nc2,rho_th_nc2,v_th_nc2,T_th_nc2]=non_cons2(rho,v,T,a,x,nt,dt,n);
[rho_c, v_c, T_c,mfr_c,rho_th_c,v_th_c,T_th_c]=cons(rho_p,v_p,T_p,a,x,nt,dt,n);
%Plotting the steady state distribution of primitive variables:
subplot(3,1,1)
plot(x,rho_nc2,'--',x,rho_c,':')
title('Density along the nozzle')
xlabel('x') 
ylabel('Density')
legend('Non-conservative','Conservative')
subplot(3,1,2)
plot(x,v_nc2,'--',x,v_c,':')
title('Velocity along the nozzle')
xlabel('x') 
ylabel('Velocity')
legend('Non-conservative','Conservative')
subplot(3,1,3)
plot(x,T_nc2,'--',x,T_c,':')
title('Temperature along the nozzle')
xlabel('x') 
ylabel('Temperature')
legend('Non-conservative','Conservative')
%Plotting the time-stepwise variation of primitive variables:
tstep=1:1400;
figure(2)
subplot(3,1,1)
plot(tstep,rho_th_nc2,'--',tstep,rho_th_c,':')
xlabel('Time-Step')
ylabel('Throat Density')
title('Throat Density vs Time-Step')
legend('Non-conservative','Conservative')
subplot(3,1,2)
plot(tstep,v_th_nc2,'--',tstep,v_th_c,':')
xlabel('Time-Step')
ylabel('Throat Velocity')
title('Throat Velocity vs Time-Step')
legend('Non-conservative','Conservative')
subplot(3,1,3)
plot(tstep,T_th_nc2,'--',tstep,T_th_c,':')
xlabel('Time-Step')
ylabel('Throat Temperature')
title('Throat Temperature vs Time-Step')
legend('Non-conservative','Conservative')
%Plotting the mass flow rate at the throats vs time-step:
figure(3)
plot(tstep,mfr_nc2,'--',tstep,mfr_c,':')
title('Mass Flow Rate at Different Time-Steps')
xlabel('Time-Step')
ylabel('MFR at the Throat')
legend('Non-conservative','Conservative')
figure(4)
mfr_c_ss=rho_c.*a.*v_c;
mfr_nc2_ss=rho_nc2.*a.*v_nc2;
plot(x,mfr_nc2_ss,'--',x,mfr_c_ss,':')
title('Mass Flow Rate along the Nozzle')
xlabel('x')
ylabel('Normalised Mass Flow Rate')
