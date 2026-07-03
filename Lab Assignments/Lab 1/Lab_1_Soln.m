% Lab 1 Solutions

%% Testing code - you can copy and paste this section into your MATLAB file, and run it to check if your functions give the desired outputs. You can match them with the results shared (i.e. the output of this script)
elements = system_elements(2,3,-1,2);
fprintf("The elements of the system F(2,3,-1,2) is:\n");
disp(elements);
fprintf("The number of elements in the system is: %d\n", no_of_elements(2,3,-1,2));
if length(elements) == no_of_elements(2,3,-1,2)
    fprintf("This matches the number of elements in the list displayed above, which we got from the system_elements function!\n");
end
fprintf("The Underflow Level of the system is: %f\n", UFL(2,-1));
fprintf("The Overflow Level of the system is: %f\n", OFL(2,3,2));

% Other test cases
% Case 1
elements = system_elements(10,1,1,1); 
% We can easily see that this system must consist of integers from -9 to 9.
% Let us check if our functions reflect the same.
fprintf("\nThe elements of the system F(10,1,1,1) is:\n");
disp(elements);
fprintf("The number of elements in the system is: %d\n", no_of_elements(10,1,1,1));
if length(elements) == no_of_elements(10,1,1,1)
    fprintf("This matches the number of elements in the list displayed above, which we got from the system_elements function!\n");
end
fprintf("The Underflow Level of the system is: %f\n", UFL(10,1));
fprintf("The Overflow Level of the system is: %f\n", OFL(10,1,1));

fprintf("\nVisualizing the elements of both systems on the number line, we get:\n");
plot_system(2,3,-1,2);
plot_system(10,1,1,1);

%% Problem 1

% Function to generate all representable floating-point numbers in a system

function y = system_elements(beta, t, L, U)
% let n be the total number of elements in the system
n = 2*(beta-1)*(beta^(t-1))*(U-L+1)+1;
y = zeros(n,1); % pre-allocating memory for y
m = zeros(beta^(t-1), t-1); % m will store the different possible mantissas (excluding d1)
% filling in the matrix, m
for k = 1:t-1
    for times = 0:beta^(t-k-1)-1
        for i = 0:beta-1
            for j = 1:beta^(k-1)
                m(times*(beta^k)+i*(beta^(k-1))+j,k) = i;
            end
        end
    end
end
ind = 2;
for e = L:U
    for i = 1:beta^(t-1)
        for j = 1:(beta-1)
            s = 0;
            for k = 1:t-1
                s = s + ((1/beta)^(k+1))*m(i,k);
            end
            s = s+ (1/beta)*j;
            s = s*(beta^e);
            y(ind,1) = s;
            y(ind+1,1) = -1*s;
            ind = ind + 2;
        end
    end
end
y = sort(y);
end

% Functions to return the UFL, OFL, and number of elements in a system
function y = UFL(beta, L)
y = beta^(L-1);
end

function y = OFL(beta, t, U)
y = (1 - beta^(-t))*beta^U;
end

function y = no_of_elements(beta, t, L, U)
y = 2*(beta-1)*(beta^(t-1))*(U-L+1)+1;
end

%% Problem 2

% Function to plot all the points in a given floating point system
function [] = plot_system(beta, t, L, U)
y = system_elements(beta, t, L, U);
plot(y,0, 'b.');
end