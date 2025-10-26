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
%Defining the schemes:
err=1;
tol=0.0001;
j_iter=1; gs_iter=1; sor_iter=1; sorf=1.5;
%Input for scheme:- 1->jacobi, 2->gauss-seidl, 3->SOR;
disp('Choose the numerical scheme to use:')
disp('  1 → Jacobi')
disp('  2 → Gauss-Seidel')
disp('  3 → Successive Over-Relaxation (SOR)')
scheme = input('Enter your choice (1/2/3): ');

if scheme==1
while err>tol
    for i=2:nx-1
        for j=2:ny-1
            T(i,j)=K*((Told(i+1,j)+Told(i-1,j))/dx^2+(Told(i,j+1)+Told(i,j-1))/dy^2);
        end
    end
    err=max(max(abs(T-Told)));
    Told=T;
    j_iter=j_iter+1;
end
 Text=sprintf('Jacobi Scheme\nNumber of Iterations=%d',j_iter);
end

if scheme==2
    while err>tol
    for i=2:nx-1
        for j=2:ny-1
            T(i,j)=K*((T(i+1,j)+T(i-1,j))/dx^2+(T(i,j+1)+T(i,j-1))/dy^2);
        end
    end
    err=max(max(abs(T-Told)));
    Told=T;
    gs_iter=gs_iter+1;
    end
    Text=sprintf('Gauss-Seidl Scheme\nNumber of Iterations=%d',gs_iter);
end

if scheme==3
    while err>tol
    for i=2:nx-1
        for j=2:ny-1
            T(i,j)=T(i,j)+sorf*(K*((T(i+1,j)+T(i-1,j))/dx^2+(T(i,j+1)+T(i,j-1))/dy^2)-T(i,j));
        end
    end
    err=max(max(abs(T-Told)));
    Told=T;
    sor_iter=sor_iter+1;
    end
    Text=sprintf('SOR Scheme\nNumber of Iterations=%d',sor_iter);
end
%plotting results
contourf(x,y,T)
colorbar
xlabel('x-axis');
ylabel('y-axis');
title(Text);

