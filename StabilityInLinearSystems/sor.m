function [x_sor,sor_iter,sr]=sor(A,b,mf)
L=tril(A); U=triu(A); d=mf*diag(diag(A));d1=diag(diag(A));
%L,U,D decomposition of the inputted matrix:
l=L-d1; u=U-d1;
sorf=0.9;
T_sor=(eye(3)-sorf*(inv((l+d))*u+eye(3)));
C_sor=sorf*(inv(l+d)*b);
tol=0.0001;
err_sor=1;
x_old=ones(3,1);
x_sor=ones(3,1);
sor_iter=1;

while err_sor>tol
    x_sor=T_sor*x_old+C_sor;
    err_sor=max(max(abs((x_sor-x_old))));
    x_old=x_sor;
     %Introducing an "if" statement to break the code when the solution
    %diverges:
     if sor_iter>10000
        break
     end
    sor_iter=sor_iter+1;
end
lambda=eig(T_sor);
sr=max(abs(lambda));
end



