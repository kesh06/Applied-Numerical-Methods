%% Testing Code

% Problem 1
%implements Gaussian elimination with partial pivoting

clc
A = [ 2  0  1  -1;
      6  3  2  -1;
      4  3 -2   3;
     -2 -6  2 -14];
 
b = [6; 15; 3; 12];

x = Gaussian_elimination_partial_pivoting(A,b);

%display the solution vector
fprintf('Solution vector\n');
fprintf('%8.4f\n', x');

% Problem 2
%implements Gaussian elimination with scaled partial pivoting

clc 
A = [     pi  -exp(1)   sqrt(2) -sqrt(3);
          pi   exp(1) -exp(1)^2      3/7;
     sqrt(5) -sqrt(6)         1 -sqrt(2);
        pi^3 exp(1)^2  -sqrt(7)     1/9];

b = [sqrt(11); 0; pi; sqrt(2)];

x = Gaussian_elimination_scaled_partial_pivoting(A,b);

%display the solution vector
fprintf('Solution vector\n');
fprintf('%8.4f\n', x');


%% Problem 1
function x = Gaussian_elimination_partial_pivoting(A,b)
    A = [A b];
     
    n = size(A,1);
    
    for i = 1:n-1
        %find indices for max value in column of submatrix of interest
        maxindices  = find(abs(A(i:n, i)) == max(abs(A(i:n, i))));
    
        %performs row swap if necessary
        if maxindices(1) ~= 1  % maxindicies could be a vector 
            temp = A(i, :);    % as there may be two maximums
            A(i, :) = A(maxindices(1)+i-1, :);
            A(maxindices(1)+i-1, :) = temp;
            fprintf('Step %i: Swapped row %i and row %i\n', i, i, maxindices(1)+i-1)
        else
            fprintf('Step %i: No row swap\n', i)
        end
    
        %performs row ops
        for j = i+1:n
            A(j, i+1:end) = A(j, i+1:end)-A(j, i)/A(i, i)*A(i, i+1:end);
        end
    
        %sets values of elements below pivot element to zero
        A(i+1:end, i) = zeros(n-i, 1);   
    end
    
    fprintf('\n')
    
    %creates a vector to store the solution
    x = zeros(n, 1);
    
    %solves for x_n
    x(n) = A(n, n+1)/A(n, n);
    
    %uses back substitution to solve for x_i
    for i = n-1:-1:1
        x(i) = (A(i, n+1)-dot(A(i, i+1:n), x(i+1:n)))/A(i, i);
    end
end

%% Problem 2
function x = Gaussian_elimination_scaled_partial_pivoting(A,b)
    %creates column vector for scaling values
    scale = max(abs(A'))';
    
    A = [A b];
     
    n = size(A,1);
    
    for i = 1:n-1
        %find indices for max value of scaled column elements in submatrix of interest
        maxindices = find(abs(A(i:n, i)./scale(i:n)) == max(abs(A(i:n, i))./scale(i:n)));
    
        %performs row swap (and scale element swap!) if necessary
        if maxindices(1) ~= 1
            temprow = A(i, :);
            A(i, :) = A(maxindices(1)+i-1, :);
            A(maxindices(1)+i-1, :) = temprow;
            tempscale = scale(i);
            scale(i) = scale(maxindices(1)+i-1);
            scale(maxindices(1)+i-1) = tempscale;
            fprintf('Step %i: Swapped row %i and row %i\n', i, i, maxindices(1)+i-1)
        else
            fprintf('Step %i: No row swap\n', i)
        end
    
        %performs row ops
        for j = i+1:n
            A(j, i+1:end) = A(j, i+1:end)-A(j, i)/A(i, i)*A(i, i+1:end);
        end
    
        %sets values of elements below pivot element to zero
        A(i+1:end, i) = zeros(n-i, 1);   
    end
    
    fprintf('\n')
    
    %creates a vector to store the solution
    x = zeros(n, 1);
    
    %solves for x_n
    x(n) = A(n, n+1)/A(n, n);
    
    %uses back substitution to solve for x_i
    for i = n-1:-1:1
        x(i) = (A(i, n+1)-dot(A(i, i+1:n), x(i+1:n)))/A(i, i);
    end
end
