% This example shows how to calculate and plot both the
% fundamental TE and TM eigenmodes of an example 3-layer ridge
% waveguide using the full-vector eigenmode solver.  

% Refractive indices:
n1 = 3.34;          % Lower cladding
n2 = 3.44;          % Core
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding

% Horizontal dimensions:
rh = 1.1;           % Ridge height
rw = 1.0;           % Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 1;         % number of modes to compute

r = (3.44 - 3.305) / 10;
fig = 1;
neffArray = [];
for i = 3.305 : r : 3.44
    [x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,i,n3],[h1,h2,h3], ...
                                            rh,rw,side,dx,dy); 
    [Hx,Hy,neff] = wgmodes(lambda,i,nmodes,dx,dy,eps,'000A');
    neffArray = [neffArray, neff];
    
    %fprintf(1,'neff = %.6f\n', neff);
    
    figure(fig);
    subplot(211);
    contourmode(x,y,Hx);
    title('Hx (TE mode)'); xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end
    
    subplot(212);
    contourmode(x,y,Hy);
    title('Hy (TE mode)'); xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end
    fig = fig + 1;
end

figure(12)
subplot(1,1,1)
plot(rr, neffArray)
title('n_{eff} as the Ridge Index Increases'); xlabel('Ridge Index'); ylabel('n_{eff}'); 