% Lab 2 Solutions

%% Testing code - you can copy and paste this section into your MATLAB file, and run it to check if your functions give the desired outputs. You can match them with the results shared (i.e. the output of this script)

% Problem 1
test1 = roots(6, 5, -4);
if (test1(1) == 0)
    type1 = 'real';
else
    type1 = 'complex';
end
fprintf("Solutions to 6x^2 + 5x - 4 = 0 are %s and are equal to %.5f and %.5f.\n \n", type1, test1(2), test1(3));

test2 = roots(6*10^154, 5*10^154, -4*10^154);
if (test2(1) == 0)
    type2 = 'real';
else
    type2 = 'complex';
end
fprintf("Solutions to (6*10^154)x^2 + (5*10^154)x - (4*10^154) = 0 are %s and are equal to %.5f and %.5f.\n \n", type2, test2(2), test2(3));

test3 = roots(0,1,1);
if (test3(1) == 0)
    type3 = 'real';
else
    type3 = 'complex';
end
fprintf("Solutions to x + 1 = 0 are %s and are equal to %.5f and %.5f.\n \n", type3, test3(2), test3(3));


test4 = roots(1,-10^5,1);
if (test4(1) == 0)
    type4 = 'real';
else
    type4 = 'complex';
end
fprintf("Solutions to x^2 - (10^5)x + 1 = 0 are %s and are equal to %.5f and %.5f.\n \n", type4, test4(2), test4(3));

test5 = roots(1, -4, 3.999999);
if (test5(1) == 0)
    type5 = 'real';
else
    type5 = 'complex';
end
fprintf("Solutions to x^2 - 4x + 3.999999 = 0 are %s and are equal to %.5f and %.5f.\n \n", type5, test5(2), test5(3));

test6 = roots(10^(-155), 10^155, 10^155);
if (test6(1) == 0)
    type6 = 'real';
else
    type6 = 'complex';
end
fprintf("Solutions to (10^-155)x^2 + (10^155)x + (10^155) = 0 are %s and are equal to %.5f and %.5f.\n \n", type6, test6(2), test6(3));

test7 = roots(1, 1, 1);
if (test7(1) == 0)
    type7 = 'real';
else
    type7 = 'complex';
end
fprintf("Solutions to x^2 + x + 1 = 0 are %s and have real and imaginary parts- %.5f and %.5f.\n \n", type7, test7(2), test7(3));

% Problem 2
x = 1;                      % point of evaluation
K = 16;                     % determines step sizes
plot_error(x,K);

%% Problem 1

% Function to find roots of a quadratic equation

function y = roots(a, b, c)
    % Scaling to avoid overflow and underflow
    s = max([abs(a), abs(b), abs(c)]);
    
    % Handle the trivial case where all coefficients are zero
    if s == 0
        error('All coefficients are zero; equation is indeterminate.');
    end

    as = a / s;
    bs = b / s;
    cs = c / s;

    % Handle the Linear Case (a is zero or relatively negligible)
    % Using 'as' here is safer for numerical precision, this also includes
    % underflow.
    if as == 0
        % Equation is bs*x + cs = 0 => x = -cs/bs
        root = -cs / bs;
        y = [0, root, root];
        return;
    end

    % Calculate Discriminant of scaled coefficients
    d_sq = bs^2 - 4*as*cs;

    if d_sq < 0
        % CASE: COMPLEX ROOTS
        type = 1;
        % Real part: -bs / (2*as)
        % Imaginary part: sqrt(-d_sq) / (2*as)
        real_part = -bs / (2 * as);
        imag_part = sqrt(-d_sq) / (2 * as);
        y = [type, real_part, imag_part];
        
    else
        % CASE: REAL ROOTS
        type = 0;
        % Use the sign of bs to avoid catastrophic cancellation error
        if bs > 0
            Q = -bs-sqrt(d_sq);
            x1 = 2*cs/Q; x2 = Q/(2*as);
        else
            Q = -bs+sqrt(d_sq);
            x1 = Q/(2*as); x2 = 2*cs/Q;
        end
        y = [type, x1, x2];
    end
end

%% Problem 2

function y = fin_diff_tan(x,h)
    y = (tan(x+h) - tan(x))/h;
end

function e = error_fd_tan(x,h)
    e = abs(sec(x)^2 - fin_diff_tan(x,h));
end

function y = cent_diff_tan(x,h)
    y = (tan(x+h) - tan(x-h))/(2*h);
end

function e = error_cd_tan(x,h)
    e = abs(sec(x)^2 - cent_diff_tan(x,h));
end

function [] = plot_error(x,K)
    H = 10.^-(0:K);            % step sizes
    
    error_fd = zeros(size(H));
    error_cd = zeros(size(H));
    
    for i = 1:length(H)
        error_fd(i) = error_fd_tan(x, H(i));
        error_cd(i) = error_cd_tan(x, H(i));
    end
    
    [~, fd_idxmin] = min(error_fd);
    [~, cd_idxmin] = min(error_cd);
    
    % Plot errors
    loglog(H, error_fd, 'r-o', ...
           'LineWidth', 2, 'MarkerSize', 7);
    hold on
    loglog(H(2:7), 1e2*H(2:7).^1,'k --', 'LineWidth',2)
    loglog(H, error_cd, 'b--*', ...
           'LineWidth', 2, 'MarkerSize', 7);
    loglog(H(3:6), 1e2*H(3:6).^2 - 10^(-5),'k -', 'LineWidth',2)
    loglog(H(fd_idxmin), error_fd(fd_idxmin), 'ko', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'y','LineWidth', 2);
    loglog(H(cd_idxmin), error_cd(cd_idxmin), 'ko', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'm','LineWidth', 2);

    xlabel('h');
    ylabel('Error');
    title('Finite Difference Error');
    
    legend( ...
        'Forward difference error', ...
        '$\mathcal{O}(h)$', ...
        'Central difference error', ...
        '$\mathcal{O}(h^2)$', ...
        'Minimum FD error', ...
        'Minimum CD error', ...
        'Location', 'north', ...
        'Interpreter', 'latex' ...
    );
    set(gca,'Fontsize',10)
    xlim([H(end) H(1)])
    grid on
end