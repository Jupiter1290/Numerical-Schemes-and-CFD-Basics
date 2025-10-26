function [x_jac,jacobi_iter,sr]=jacobi(A,b,mf)
%L,U,D decomposition of the inputted matrix:
L=tril(A); U=triu(A); d=mf*diag(diag(A));d1=diag(diag(A));
l=L-d1;u=U-d1;
T_jac=-inv(d)*(l+u);
C_jac=inv(d)*b;
tol=0.0001;
err_jac=1;
x_old=ones(3,1);
x_jac=ones(3,1);
jacobi_iter=1;
while err_jac>tol
    x_jac=T_jac*x_old+C_jac;
    err_jac=max(max(abs(x_jac-x_old)));
    x_old=x_jac;
    %Introducing an "if" statement to break the code when the solution
    %diverges:
    if jacobi_iter>10000
        break
    end
    jacobi_iter=jacobi_iter+1;
end
lambda=eig(T_jac);
sr=max(abs(lambda));
end    