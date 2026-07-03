clear
close all
clc

%% Problem 1 using a while loop
%uses the Bisection method to find a solution accurate to within 10^-8
%for x-2^{-x}=0 on the interval [0,1]
fprintf('\n -------Solution of Problem 1-------\n\n');
f = @(x) x-2.^(-x);
a = 0;
b = 1;
tolexp = -8;
tolerance = 10^tolexp;
maxiterations = 30;

k = 1;

while k <= maxiterations
    
    x = (a+b)/2;
    
    if f(x) == 0 || (b-a)/2 < tolerance
        break
    end
        
    if sign(f(a))*sign(f(x)) == -1
        b = x;
    else 
        a = x;
    end
    
    k = k+1;
end

if k < maxiterations+1
    fprintf('Tolerance: 10^%i, Approximation: %.8f, Iterations: %i\n', tolexp, x, k);
else
    fprintf('The max number of iterations (%i) was reached before convergence\n', maxiterations);
end


%% Problem 2
%uses the Bisection method to find a solution accurate to within 10^-12
%for x-2^{-x}=0 on the interval [0,1]
fprintf('\n -------Solution of Problem 2-------\n\n');
f = @(x) x-2.^(-x);
a = 0;
b = 1;
tolexp = -12;
tolerance = 10^tolexp;
maxiterations = 30;

k = 1;

while k <= maxiterations
    
    x = (a+b)/2;
    
    if f(x) == 0 || (b-a)/2 < tolerance
        break
    end
        
    if sign(f(a))*sign(f(x)) == -1
        b = x;
    else 
        a = x;
    end
    
    k = k+1;
end

if k < maxiterations+1
    fprintf('Tolerance: 10^%i, Approximation: %.12f, Iterations: %i\n', tolexp, x, k);
else
    fprintf('The max number of iterations (%i) was reached before convergence\n', maxiterations);
end


%% Problem 3
%plots y = x and y = 2 sin(x) and finds the first positive solution of
%x = 2 sin(x)
fprintf('\n -------Solution of Problem 3-------\n\n');
xp = 0:0.01:2.5;

figure
plot(xp, xp,'b','LineWidth',2)
hold on
plot(xp,2*sin(xp),'k','LineWidth',2)
axis([0 2.5 -0.5 2.5])
legend({'x','2sin(x)'},'Location','NorthWest')

f = @(x) x-2*sin(x);

a = 1;
b = 3;

tolexp = -8;
tolerance = 10^tolexp;
maxiterations = 100;

xold = 999;
k = 1;

while k <= maxiterations
    
    xnew = (a+b)/2;
    
    if f(xnew) == 0 || abs((xnew-xold)/xnew) < tolerance
        break
    end
        
    if sign(f(a))*sign(f(xnew)) == -1
        b = xnew;
    else
        a = xnew;
    end
    
    k = k+1;
    xold = xnew;
end

if k < maxiterations+1
    fprintf('Tolerance: 10^%i, Approximation: %.8f, Iterations: %i\n', tolexp, xnew, k);
else
    fprintf('The max number of iterations (%i) was reached before convergence\n', maxiterations);
end

%% Problem 4
%uses four different fixed-point iteration functions to approximate 
%21^(1/3) accurate to within 10^{-10}; compares the speed of 
%convergence based on the number of iterations for the four methods
fprintf('\n -------Solution of Problem 4-------\n\n');
x0 = 1; 
tolexp = -10;
maxiterations = 150;

%first fixed-point iteration function
g = @(x) (20*x+21/x^2)/21;

fpiteration(g, x0, maxiterations, tolexp);

%second fixed-point iteration function
g = @(x) x-(x^3-21)/(3*x^2);

fpiteration(g, x0, maxiterations, tolexp);

%third fixed-point iteration function
g = @(x) x-(x^4-21*x)/(x^2-21);

fpiteration(g, x0, maxiterations, tolexp);

%fourth fixed-point iteration function
g = @(x) (21/x)^(1/2);

fpiteration(g, x0, maxiterations, tolexp);

fprintf('From fastest to slowest: b, d, a (c does not converge). \n')

%---------------------------Function--------------------------------------%
function fpiteration(g, x0, maxiterations, tolexp)
%uses fixed-point iteration to generate approximations

%INPUT
%   g: function used for fixed-point iteration
%   maxiterations: max number of iterations
%   x0: initial guess
%   tolexp: exponent ofstopping criteria tolerance

    %initializes counter for first iteration
    k = 1;

    %outputs initial guess
    fprintf('k: 0\tx0: %.10f\n', x0) 

    while  k <= maxiterations
        %performs fixed-point iteration
        x = g(x0); 
    
        %outputs iteration number, approximation, and difference from last
        %approximation
        fprintf('k: %i \tx_%i: %.10f\t|x_%i-x_%i|: %.10f\n', k, k, x,...
            k, k-1, abs(x-x0)) 
    
        %checks stopping criteria
        if abs(x-x0) < 10^tolexp
            break
        end
        
        %increments iteration count
        k = k+1;
    
        %updates 'initial guess' for next iteration
        x0 = x;
    end
    
    if k < maxiterations+1
        fprintf(['Tolerance: 10^%i, Approximation: %.' num2str(abs(tolexp)+1) 'f, Iterations: %i\n\n'],...
        tolexp, x, k);
    else
        fprintf('The max number of iterations %i was reached before convergence\n\n', maxiterations);
    end
end
%-------------------------------------------------------------------------%