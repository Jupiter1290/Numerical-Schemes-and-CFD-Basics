clear
close all
clc
%all units in SI
%Length of domain(L),number of nodes(n),constant term in the equation(c) and timestep(dt) defined
L=1; n=[20 40 80 160]; c=1; dt=0.01;
for l=1:4
%Domain
x=linspace(0,L,n(l));
dx=x(2)-x(1);
mesh_length(l)=x(2)-x(1);
%Specifying the intial velocity profile:
% for calculating n_start and n_end:
arr1=abs(x-0.1); arr2=abs(x-0.3); n_start=0; n_end=0;
v=1;v1=1;
for i=1:length(arr1)
    for k=1:length(arr1)
        if arr1(i)<arr1(k)
            v=v+1;
        end
       
    end
     if v==length(arr1)
            n_start=i;
        else 
            v=1;
    end
end

for i=1:length(arr2)
    for k=1:length(arr2)
        if arr2(i)<arr2(k)
            v1=v1+1;
        end
       
    end
     if v1==length(arr2)
            n_end=i;
        else 
            v1=1;
    end
end
u=0.1*ones(1,n(l));
u(n_start:n_end)=0.3;
%Performing Time Marching Integration:
uold=u;
u_initial=u;
%Time loop (each loop corresponds to one time step(ts))
%Duration:0.4s, Size of each time-step=0.01s,thus, number of time-steps=0.4/0.01=40: 
ts=40;
for e=1:ts
    for p=2:n(l)
    u(p)=uold(p)-(c*dt/dx)*(uold(p)-uold(p-1));
    end
    %Updating the velocities for the next time step
    uold=u;
    %plotting the velocities
    figure(l);
    plot(x,u_initial,'r')
    hold on
    plot(x,u,'b')
    axis([0 1 0 0.35])
    xlabel('X')
    ylabel('Velocity')
    hold off
    pause(0.1)
end
end