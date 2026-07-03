%  Lab 10: Comparison of Quadrature Rules
%  Clenshaw–Curtis vs Gauss–Legendre

Nmax = 50;

% Exact integrals
Itrue = [ ...
    1/2;                                        % I1 = ∫ |x|^3 dx
    2*(exp(-1) + sqrt(pi)*(erf(1)-1));          % I2 = ∫ exp(-x^{-2}) dx
    pi/2;                                       % I3 = ∫ 1/(1+x^2) dx
    2/11                                        % I4 = ∫ x^{10} dx
];

labels = {'|x|^3','exp(-x^{-2})','1/(1+x^2)','x^{10}'};

% Storage for errors
E_cc = zeros(4, Nmax);
E_ga = zeros(4, Nmax);

%  Main Computation Loop

for N = 1:Nmax
    
    % Get nodes and weights
    [x_cc, w_cc] = clencurt(N);
    [x_ga, w_ga] = gauss(N);
    
    % Evaluate functions safely
    F_cc = test_functions(x_cc);
    F_ga = test_functions(x_ga);
    
    % Compute errors
    for k = 1:4
        E_cc(k,N) = abs(w_cc * F_cc(k,:)' - Itrue(k));
        E_ga(k,N) = abs(w_ga * F_ga(k,:)' - Itrue(k));
    end
end

%  Plot: Clenshaw–Curtis Errors
plot_errors(E_cc, labels, Nmax, 'Clenshaw–Curtis Quadrature Errors');

%  Plot: Gauss–Legendre Errors
plot_errors(E_ga, labels, Nmax, 'Gauss–Legendre Quadrature Errors');

%  Degree of Exactness Test (Gauss–Legendre)
fprintf('\nDegree of Exactness Test for Gauss–Legendre:\n')

for N = 1:6
    [x,w] = gauss(N);
    approx = w * (x.^10);
    error = abs(approx - 2/11);
    
    fprintf('N = %d, Error for x^{10} = %.2e\n', N, error);
end


%  Function: Test Functions
function F = test_functions(x)

    F = zeros(4,length(x));
    
    % (a) |x|^3
    F(1,:) = abs(x).^3;
    
    % (b) exp(-x^{-2}) (handle x=0 safely)
    xtmp = x;
    xtmp(xtmp==0) = eps;
    F(2,:) = exp(-xtmp.^(-2));
    
    % (c) 1/(1+x^2)
    F(3,:) = 1./(1 + x.^2);
    
    % (d) x^{10}
    F(4,:) = x.^10;
end

%  Function: Plotting
function plot_errors(E, labels, Nmax, title_str)

    figure('Position',[100 100 900 900])
    sgtitle(title_str,'FontSize',18,'FontWeight','bold')

    for k = 1:4
        subplot(2,2,k)
        semilogy(1:Nmax, E(k,:) + 1e-100,'o-','LineWidth',1.2)
        grid on
        
        xlim([1 Nmax])
        ylim([1e-18 1e+1])
        
        xlabel('N','FontSize',14)
        ylabel('Error','FontSize',14)
        title(labels{k},'FontSize',15,'FontWeight','bold')
    end
end

%  Clenshaw–Curtis Quadrature
function [x,w] = clencurt(N)

    theta = pi*(0:N)'/N;
    x = cos(theta);

    w = zeros(1,N+1);
    ii = 2:N;
    v = ones(N-1,1);

    if mod(N,2)==0
        w(1) = 1/(N^2 - 1);
        w(N+1) = w(1);
        for k = 1:(N/2 - 1)
            v = v - 2*cos(2*k*theta(ii)) / (4*k^2 - 1);
        end
        v = v - cos(N*theta(ii))/(N^2 - 1);
    else
        w(1) = 1/N^2;
        w(N+1) = w(1);
        for k = 1:((N-1)/2)
            v = v - 2*cos(2*k*theta(ii)) / (4*k^2 - 1);
        end
    end

    w(ii) = 2*v/N;
end

%  Gauss–Legendre Quadrature
function [x,w] = gauss(N)

    beta = .5 ./ sqrt(1 - (2*(1:N-1)).^(-2));
    T = diag(beta,1) + diag(beta,-1);

    [V,D] = eig(T);
    x = diag(D);
    [x,ind] = sort(x);

    w = 2 * V(1,ind).^2;
end