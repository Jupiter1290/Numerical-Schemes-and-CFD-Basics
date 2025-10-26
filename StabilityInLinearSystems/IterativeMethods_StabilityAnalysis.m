close all
clear
clc
A=[5 1 2; -3 9 4; 1 2 -7];
b=[10;-14;33]; mf=0.5:0.1:1.5;
%For loop to run the code for each value of the magnification factor
for i=1:length(mf)
    %The iterative schemes are called as functions that are user-defined:
    [x_jacobi,jacobi_iter(i),sr_jacobi(i)]=jacobi(A,b,mf(i));
    [x_gs,gs_iter(i),sr_gs(i)]=gs(A,b,mf(i));
    [x_sor,sor_iter(i),sr_sor(i)]=sor(A,b,mf(i));
end

%Plotting
figure(1)
subplot(3,1,1)
plot(mf,jacobi_iter)
xlabel('Magnifying Factor')
ylabel('Iterations')
title('Jacobi')

subplot(3,1,2)
plot(mf,gs_iter)
xlabel('Magnifying Factor')
ylabel('Iterations')
title('Gauss-Seidl')

subplot(3,1,3)
plot(mf,sor_iter)
xlabel('Magnifying Factor')
ylabel('Iterations')
title('SOR')

figure(2)
subplot(3,1,1)
plot(mf,sr_jacobi);
xlabel('Magnifying Factor');
ylabel('Spectral Radius');
title('Jacobi')

subplot(3,1,2)
plot(mf,sr_gs);
xlabel('Magnifying Factor');
ylabel('Spectral Radius');
title('Gauss-Seidl')

subplot(3,1,3)
plot(mf,sr_sor);
xlabel('Magnifying Factor');
ylabel('Spectral Radius');
title('SOR')
   
figure(3)
subplot(3,1,1)
plot(sr_jacobi,jacobi_iter);
xlabel('Spectral Radius');
ylabel('Iterations');
title('Jacobi')

subplot(3,1,2)
plot(sr_gs,gs_iter);
xlabel('Spectral Radius');
ylabel('Iterations');
title('Gauss-Seidl')

subplot(3,1,3)
plot(sr_sor,sor_iter);
xlabel('Spectral Radius');
ylabel('Iterations');
title('SOR')





