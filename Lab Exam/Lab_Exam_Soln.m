%% Problem
x_old = [-1;-2;1];
D_inv = jac_i(x_old);
Fx = F(x_old);
x = x_old - D_inv*Fx; %Update as per Newton's Method for first time
k=1; %Set iteration counter
Fx = F(x); %Evaluate function at new variable
fprintf("k\t x1\t\t x2\t\t x3\t\t ||dx||_inf\t ||F(x)||_inf\n");
fprintf("%d\t %.6f\t %.6f\t %.6f\t %.2e\t %.2e\n",...
        k,x(1),x(2),x(3),norm(x-x_old,"inf"),norm(Fx,"inf"));
while(norm(x-x_old,"inf")>=1e-8) %Keep iterating while error is more than 1e-8
    x_old = x;
    D_inv = jac_i(x_old);
    Fx = F(x_old);
    x = x_old - D_inv*Fx;
    Fx = F(x);
    k = k+1;
    fprintf("%d\t %.6f\t %.6f\t %.6f\t %.2e\t %.2e\n",...
        k,x(1),x(2),x(3),norm(x-x_old,"inf"),norm(Fx,"inf"));
end
fprintf("x = [%.6f %.6f %.6f]^T\n",x(1),x(2),x(3));
%fprintf("%.16f %.16f %.16f",Fx(1),Fx(2),Fx(3));

%%
function D_inv = jac_i(x) %Creates Jacobian Matrix evaluated at a point
    D(:,1)= [3*x(1)^2+2*x(1)*x(2)-x(3); exp(x(1));-2*x(3)];
    D(:,2)= [x(1)^2;exp(x(2));2*x(2)];
    D(:,3)= [-x(1);-1;-2*x(1)];
    D_inv = inv(D);
end
%%
function y = F(x) %Evaluates given function at some arbitrary vector say x
    y = [x(1)^3+(x(1)^2)*x(2)-x(1)*x(3)+6;exp(x(1))+exp(x(2))-x(3);x(2)^2-2*x(1)*x(3)-4];
end