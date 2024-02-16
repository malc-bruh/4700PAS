clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked')

nx = 50;
ny = 100;
nxp = 0;
nxm = 0;
nyp = 0;
nym = 0;
nmodes = 9;
G = sparse(nx * ny, nx * ny);
V = zeros(nx, ny);

for i = 1 : nx
    for j = 1 : ny
        n = j + (i - 1) * ny;
        nxm = j + (i - 2) * ny;
        nxp = j + (i) * ny;
        nym = (j - 1) + (i - 1) * ny;
        nyp = (j + 1) + (i -  1) * ny;
        
        if i == 1
            G(n, n) = 1;
        elseif i == nx
            G(n, n) = 1;
        elseif j == 1
            G(n, n) = 1;
        elseif j == ny
            G(n, n) = 1;
        elseif (i > 10 && i < 20 && j > 10 && j < 20)
            G(n, n) = -2;
            G(n, nxm) = 1;
            G(n ,nxp) = 1;
            G(n, nym) = 1;
            G(n, nyp) = 1;
        else
            G(n, n) = -4;
            G(n, nxm) = 1;
            G(n ,nxp) = 1;
            G(n, nym) = 1;
            G(n, nyp) = 1;
        end
    end
end
figure(1)
spy(G)
title("Sparsity Pattern of Matrix G");

[E, D] = eigs(G, nmodes, 'SM'); % Eigenvalues and Eigenvectors
plot(diag(D))

np = ceil(sqrt(nmodes));
for k = 1 : nmodes
    M = E(:, k);
    for i = 1 : nx
        for j = 1 : ny
            n = i + (j - 1) * nx;
            V(i, j) = M(n);
        end
        subplot(np, np, k)
        surf(V)
        title(['EV = ', num2str(D(k, k))])
    end
end