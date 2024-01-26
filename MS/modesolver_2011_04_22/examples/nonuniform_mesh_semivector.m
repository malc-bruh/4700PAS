% This example shows how to use the stretchmesh function to expand the
% size of the computation window at the edges without increasing the
% number of gridpoints.  It uses the same example waveguide
% structure considered in 'basic-semivector.m'

% Refractive indices:
n1 = 3.34;          % Lower cladding
n2 = 3.44;          % Core
n3 = 1.00;          % Upper cladding (air)

% Vertical dimensions:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding
rh = 1.1;           % Ridge height

% Horizontal dimensions:
rw = 1.0;           % Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 1;         % number of modes to compute

fprintf (1,'generating index mesh...\n');
[x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,n2,n3],[h1,h2,h3], ...
                                            rh,rw,side,dx,dy); 

% Stretch the 40 gridpoints on the south and east edges of the
% computational window by factors of 4 and 3, respectively.

[x,y,xc,yc,dx,dy] = stretchmesh(x,y,[0,40,40,0],[1,4,3,1]);

[Ex,neff] = svmodes(lambda,n2,nmodes,dx,dy,eps,'000S','EX');

fprintf(1,'neff = %.6f\n',neff);

figure(1);
contourmode(x,y,Ex);
title('Ex (TE Mode)'); xlabel('x'); ylabel('y'); 
for v = edges, line(v{:}); end
