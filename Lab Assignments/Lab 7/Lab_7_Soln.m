clc
clear all
close all
%% Problem 1
%plots the graph of x^2*abs(sin x)-4 for x in [0,4]; uses the Secant method
%to find the smallest positive zero accurate to within 10^-6 using p0 = 3.6
%and p1 = 3.7 (also p0 = 2.8 and p1 = 2.9)

%function for which we want to approximate root
f = @(x) x.^2.*abs(sin(x))-4;
%plots function on the interval [0,4]
figure(1)
plot(0:0.01:4, f(0:0.01:4),'b-','LineWidth',2)
xlabel('x'); ylabel('y')

p0 = 3.6; %initial guesses (two required for secant method)
p1 = 3.7;
tolexp = -6; %tolerance exponent
tolerance = 10^tolexp; %stopping criteria tolerance
N = 50; %max number of iterations

%uses secant method to produce approximation for root
Secant(p0, p1, f, N, tolexp)

p0 = 2.8; 
p1 = 2.9;

Secant(p0, p1, f, N, tolexp) %converges to another zero

%% Problem 2
% Use Newton’s method to approximate the zero of the function 
% $f(x) = x^2 - 2e^{-x}x + e^{-2x}$ starting from $p_0 = 1$ using accurate 
% to within $10^{-8}$. What do you notice about the convergence of Newton’s
% method? 
% ==> Note that here the root has multiplicity > 1, that is f(r) = f'(r)=0
% and hence Newton's method may not convergence quadratically. 

% Function
f = @(x) x.^2 - 2*x.*exp(-x) + exp(-2*x);

% Derivative
fprime = @(x) 2*(x - exp(-x)).*(1 + exp(-x));

% Plot f and f'
x = 0:0.01:2;
figure(2)
plot(x, f(x), 'b-', 'LineWidth', 2); hold on
plot(x, fprime(x), 'r-', 'LineWidth', 2)
legend('f(x)', 'f''(x)')
xlabel('x'); ylabel('y')
title('Function and its derivative')
grid on

% Parameters
tolexp = -10;          % tighter tolerance for better "true root"
maxiterations = 50;
p0 = 1;

[p_vals, n] = Newton_Scalar(f, fprime, tolexp, maxiterations, p0);

r = p_vals(n);

% True errors
errors_true = abs(p_vals(1:n) - r);
% Instead of measuring |pn-r| we could also measure |pn-p_(n+1)|, since
% near the root |pn-p_(n+1)| = |pn-r|

kmax = length(errors_true) - 10;   % ignore last 3 points

e_n   = errors_true(1:kmax-1);
e_np1 = errors_true(2:kmax);

figure
loglog(e_n, e_np1, 'o-','LineWidth',2); hold on
xlabel('e_n')
ylabel('e_{n+1}')
title('Order of convergence (restricted range)')
grid on

% Rate of convergence

% pick a small middle segment
i1 = floor(length(e_n)/2);
i2 = i1 + 5;

x_ref = e_n(i1:i2);

% slope 1: y = C*x (match locally)
C = e_np1(i1)/e_n(i1);
y_ref = 1e+1*C * x_ref;
figure(3)
loglog(x_ref, y_ref, 'k--','LineWidth',2)
legend('Data','Slope = 1','Location','best')

%% Problem 3
% Use Newton’s method to find a solution to the following nonlinear system.
% Iterate until ||x^(k)-x^(k-1)|| < 10^(-6).
% 6*x1 - 2*cos(x2*x3) - 1 = 0
% 9*x2 + \sqrt(x1^2 + sin(x3) + 1.06) + 0.9 = 0,
% 60*x3 + 3e^{-x1*x2} + 10*pi - 3 = 0.
 

x0 = [0;0;0]; tol = 1e-6; maxit = 50;
x = Newton_system(@f_fun, @Jacobian, x0, tol, maxit);
fprintf('\nFinal solution:\n');
disp(x)


f_fun(x)
%------------------------------Functions----------------------------------%
function Secant(p0, p1, f, maxiterations, tolexp)
%function uses the Secant method to approximate zero of a function

%INPUT
%  p0: first initial guess
%  p1: second initial guess
%  f: function for which we want to approximate a zero
%  maxiterations: max number of iterations
%  tolexp: exponent of stopping criteria tolerance
  
    tolerance = 10^tolexp;
    
    n = 2; %iteration counter  
    fp0 = f(p0); %computes function value of initial guesses
    fp1 = f(p1);

    fprintf(['n: 1\tp%i: %.' num2str(abs(tolexp)+1) 'f\t|f(p%i)|: %.8f\t|(p%i-p%i)/p%i|: %.' num2str(abs(tolexp)+1) 'f\n'], n-1, p1, n-1, abs(fp1),n-1, n-2, n-2, abs((p1-p0)/p0)) 

    while n <= maxiterations
        %produces new approximation of root
        p = p1-fp1*(p1-p0)/(fp1-fp0);
    
        fprintf(['n: %i\tp%i: %.' num2str(abs(tolexp)+1) 'f\t|f(p%i)|: %.8f\t|(p%i-p%i)/p%i|: %.' num2str(abs(tolexp)+1) 'f\n'],...
            n, n, p, n, abs(f(p)), n, n-1, n-1, abs((p-p1)/p1)) 
    
        if abs((p-p1)/p1) < tolerance
             break
        end
        
        %updates the iteration counter
        n = n+1;
    
        %updates first guess with the second guess
        p0 = p1;
        fp0 = fp1;
        %updates second guess with the new approximation
        p1 = p;
        fp1 = f(p);
    end

    if n < maxiterations+1
        fprintf('Tolerance: 10^%i, Approximation: %.10f, Iterations: %i\n', tolexp, p, n);
    else
        fprintf('The max number of iterations %i was reached before convergence\n', maxiterations);
    end

    %insert carriage return to separate output
    fprintf('\n')
end


function [p_vals, n] = Newton_Scalar(f, fprime, tolexp, maxiterations, p0)

    p_vals = zeros(maxiterations+1,1);
    p_vals(1) = p0;

    fprintf(['n: 0\tp0: %.' num2str(abs(tolexp)+1) 'f\n'], p0)

    for n = 1:maxiterations

        p = p0 - f(p0)/fprime(p0);
        p_vals(n+1) = p;

        relerr = abs((p - p0)/p);

        fprintf(['n: %d\tp%d: %.' num2str(abs(tolexp)+1) ...
            'f\t|(p%d-p%d)/p%d|: %.' num2str(abs(tolexp)+1) 'f\n'], ...
            n, n, p, n, n-1, n, relerr)

        if relerr < 10^tolexp
            break
        end

        p0 = p;
    end

    fprintf('\nTolerance: 10^%d, Approximation: %.*f, Iterations: %d\n\n', ...
        tolexp, abs(tolexp)+1, p, n);

end
%
function x = Newton_system(Ffun, Jfun, x0, tol, maxit)
    x = x0;
    err = 1.0;
    it = 0;
    fprintf('k\t   x1\t\t   x2\t\t   x3\t\t   ||dx||_inf\n');
    while err > tol && it < maxit
        Fx = Ffun(x);
        Jx = Jfun(x);
        dx = Jx \ Fx;
        x = x - dx;
        err = max(abs(dx));
        it = it + 1;

        fprintf('%d\t% .6f\t% .6f\t% .6f\t% .2e\n', ...
            it, x(1), x(2), x(3), err);
    end
end
%
function F = f_fun(x)
    F = zeros(3,1);
    F(1) = 6*x(1) - 2*cos(x(2)*x(3)) - 1;
    F(2) = 9*x(2) + sqrt(x(1)^2 + sin(x(3)) + 1.06) + 0.9;
    F(3) = 60*x(3) + 3*exp(-x(1)*x(2)) + 10*pi - 3;
end
%
function J = Jacobian(x)
    J = zeros(3,3);
    % Row 1
    J(1,1) = 6;
    J(1,2) = 2*x(3)*sin(x(2)*x(3));
    J(1,3) = 2*x(2)*sin(x(2)*x(3));
    % Row 2
    denom = sqrt(x(1)^2 + sin(x(3)) + 1.06);
    J(2,1) = x(1)/denom;
    J(2,2) = 9;
    J(2,3) = cos(x(3))/(2*denom);
    % Row 3
    J(3,1) = -3*x(2)*exp(-x(1)*x(2));
    J(3,2) = -3*x(1)*exp(-x(1)*x(2));
    J(3,3) = 60;
end
%-------------------------------------------------------------------------%

