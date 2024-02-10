clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked')

% Setting x and y boundaries of the solution
nx = 100;
ny = 100;
ni = 7000; % Number of steps
V = zeros(nx, ny);

%Boundary Conditions
V(1, :) = 1; % Left side BC is 1V
V(nx, :) = 1; % Right side BC is 0V

for m = 1 : ni
    for i = 2 : nx - 1
        for j = 2 : ny - 1
            V(i, j) = (V(i + 1, j) + V(i - 1, j) + V(i, j + 1) + V(i, j - 1)) / 4; % Update Equation
        end 
    end
    
    if((V(1, :) == 1) & (V(nx, :) == 0))
        V(:, 1) = V(:, 2); % Bottom BC
        V(:, ny) = V(:, ny - 1); % Top BC
    else 
        V(:, 1) = 0; % Bottom BC
        V(:, ny) = 0; % Top  BC
    end


    if(mod(m, 50)==0)
        figure(1)
        surf(V)
        title("Voltage by Solving Laplace's Equation using Iteration")
        xlabel("x")
        ylabel("y")
        pause(0.01);
    end
end

[Ex, Ey] = gradient(V);

figure(2)
quiver(-Ex, -Ey, 8)
title("Electric Field Vector Plot")
xlabel("x")
ylabel("y")

boxFilter = imboxfilt(V, 3);
figure(3)
surf(boxFilter)