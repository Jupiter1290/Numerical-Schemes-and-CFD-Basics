function [x_jac,jacobi_iter,sr]=jacobiX(A,b,mf)
L=tril(A); U=triu(A); d=mf*diag(diag(A));
l=L-d;u=U-d;
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
    jacobi_iter=jacobi_iter+1;
end
syms lambda;
p=det(T_jac-lambda*eye(3));
e=solve(p==0,lambda);
sr=max(abs(double((e))));
end    