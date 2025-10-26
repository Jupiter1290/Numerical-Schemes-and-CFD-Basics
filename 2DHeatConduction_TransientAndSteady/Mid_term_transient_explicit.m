clear
close all
clc

%Defining the domain:(number of nodes along x and y directions are 10)
nx=10; ny=10;
x=linspace(0,1,nx);y=linspace(1,0,ny);
dx=x(2)-x(1); dy=y(2)-y(1);
K=dx^2*dy^2/(2*(dx^2+dy^2));
T=300*ones(nx,ny);
Ttop=600; Tbottom=900; Tl=400; Tr=800;
T(1,:)=Ttop; T(ny,:)=Tbottom; T(:,1)=Tl; T(:,nx)=Tr;
T(1,1)=(Ttop+Tl)/2; T(nx,1)=(Ttop+Tr)/2; T(nx,ny)=(Tbottom+Tr)/2;
T(ny,1)=(Tbottom+Tl)/2;
Told=T;
T_i=Told;
time_simulation=54; dt=0.001;
%time-step taken such that 'k1+k2</=0.5' to meet stability criterion;
alpha=0.5;
k1=alpha*dt/dx^2; k2=alpha*dt/dy^2;
%dt is the time step size; calculating the number of time-steps=t;
nt=time_simulation/dt;
for e=1:nt
    for i=2:nx-1
        for j=2:ny-1
    T(i,j)=Told(i,j)*(1-2*k1-2*k2)+k1*(Told(i+1,j)+Told(i-1,j))+k2*(Told(i,j+1)+Told(i,j-1));
        end
    end
    Told=T;
end

Text=sprintf('Explicit Scheme');

contourf(x,y,T)
colorbar
xlabel('x-axis');
ylabel('y-axis');
title(Text);


