function [Ez_t, Ez_r, Ez, Hx, Hy] = ...
    Yee2DEM(nx,ny,epi,mu,sigma,sigmaH,xMax,tsim,nSteps,Ez0,Hx0,Hy0,...
    bc,pml,Plot,...
    Reg,movie)
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global mainTitle;
global MetaSurf;
global Pzt_total;
global Pzr_total;
global Ez_output_t;
global Ez_output_r;
global E0_input;

global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight;

% subplot orientations
if Plot.ori == '13'
    np = 1;
    mp = 3;
else
    np = 3;
    mp = 1;
end

% standard yee 2d Cell code
ds = xMax{1}/nx{1};

%Courant criteria
for k = 1:Reg.n
    Maxdt(k) = min(min(sqrt(epi{k}.*mu{k}/c_eps_0/c_mu_0)))*0.5*ds/c_c;
end

Maxdt = max(Maxdt);

dt = tsim/nSteps;
if dt > Maxdt
    dt = Maxdt;
    nSteps = round(tsim/dt);
    display(['Courant criteria not satisfied reseting time step: ' num2str(nSteps)]);
end


% for every region set up the paramaters

for k = 1:Reg.n
    
    yMax{k} = ny{k}*ds;
    xMax{k} = nx{k}*ds;
    
    % boundary and center positions of yee cells (no offset)
    xB{k} = [0:nx{k}]*ds; % Hy
    yB{k} = [0:ny{k}]*ds; % Hx
    
    %     if surf.Mode == 't' || surf.Mode == 'i'
    %         xB{k}(surf.XOffset+1:end) = xB{k}(surf.XOffset+1:end)-ds;
    %     end
    
    xC{k} = [0.5:nx{k}-0.5]*ds; % Ez/Hx
    yC{k} = [0.5:ny{k}-0.5]*ds; % Ez/Hy
    
    %     if surf.Mode == 't' || surf.Mode == 'i'
    %         xC{k}(surf.XOffset:end) = xC{k}(surf.XOffset:end)-ds;
    %         xB{k}(surf.XOffset+1:end) = xB{k}(surf.XOffset+1:end)-ds;
    %     end
    
    
    % Initial fields
    Ez{k} = Ez0{k};
    Hx{k} = Hx0{k};
    Hy{k} = Hy0{k};
    
    %% pml
    
    pmlMax{k}=(0.8*(pml.m+1))./(sqrt(c_mu_0/c_eps_0)*ds*sqrt(mu{k}.*epi{k}./c_eps_0./c_mu_0));
    
    pmlTh{k} = pml.width+1;
    pmlS{k} = ([1:pmlTh{k}]./pml.width).^pml.m;
    pmlSm{k} = fliplr(pmlS{k});
    pmlX{k} = repmat(pmlS{k},ny{k},1)';
    pmlXm{k} = repmat(pmlSm{k},ny{k},1)';
    
    pmlY{k} = repmat(pmlS{k}',1,nx{k})';
    pmlYm{k} = repmat(pmlSm{k}',1,nx{k})';
    
    pmlSigmax{k} = zeros(nx{k},ny{k});
    
    if bc{k}.xm.type == 'a'
        pmlSigmax{k}(1:pmlTh{k},:) = pmlMax{k}(1:pmlTh{k},:).*pmlXm{k};
    end
    
    if bc{k}.xp.type == 'a'
        pmlSigmax{k}(nx{k}-pmlTh{k}+1:nx{k},:) = pmlMax{k}(nx{k}-pmlTh{k}+1:nx{k},:).*pmlX{k};
    end
    
    pmlSigmay{k} = zeros(nx{k},ny{k});
    
    if bc{k}.ym.type == 'a'
        pmlSigmay{k}(:,1:pmlTh{k}) = pmlMax{k}(:,1:pmlTh{k}).*pmlYm{k};
    end
    
    if bc{k}.yp.type == 'a'
        pmlSigmay{k}(:,ny{k}-pmlTh{k}+1:ny{k}) = pmlMax{k}(:,ny{k}-pmlTh{k}+1:ny{k}).*pmlY{k};
    end
    
    byHx{k} = exp(-(pmlSigmay{k}.*dt/c_eps_0));
    bxHy{k} = exp(-(pmlSigmax{k}.*dt/c_eps_0));
    bxEz{k} = exp(-(pmlSigmax{k}.*dt/c_eps_0));
    byEz{k} = exp(-(pmlSigmay{k}.*dt/c_eps_0));
    cyHx{k} = (1-byHx{k})/ds;
    cxHy{k} = (1-bxHy{k})/ds;
    cxEz{k} = (1-bxEz{k})/ds;
    cyEz{k} = (1-byEz{k})/ds;
    
    dHydx{k} = zeros(nx{k},ny{k});
    dHxdy{k} = zeros(nx{k},ny{k});
    
    QHxy{k} = zeros(nx{k},ny{k});
    QHyx{k} = zeros(nx{k},ny{k});
    QEyz{k} = zeros(nx{k},ny{k}-1);
    QExz{k} = zeros(nx{k}-1,ny{k});
    
    %% coeff
    Ca{k} = (1 - sigma{k}*dt./(2*epi{k}))./(1 + sigma{k}*dt./(2*epi{k}));
    Cb{k} = (dt./(epi{k}*ds))./(1 + sigma{k}*dt./(2*epi{k}));
    
    Da{k} = (1 - sigmaH{k}*dt./(2*mu{k}))./(1 + sigmaH{k}*dt./(2*mu{k}));
    Db{k} = (dt./(mu{k}*ds))./(1 + sigmaH{k}*dt./(2*mu{k}));
    
    figure
    subplot(1,2,1),surface((xC{k}+Reg.xoff{k})/1e-6,(yC{k}+Reg.yoff{k})/1e-6,epi{k}.','linestyle','none'); title('\epsilon');
    subplot(1,2,2),surface((xC{k}+Reg.xoff{k})/1e-6,(yC{k}+Reg.yoff{k})/1e-6,mu{k}.','linestyle','none');title('\mu');
    
%% Simulation

    if ~Plot.off || (movie || simulationStopTimes(1) ~= 0)
        figure
        subplot(np,mp,1),s1(k) = surface(xC{k}+Reg.xoff{k},yC{k}+Reg.yoff{k},real(Ez{k}.'),'linestyle','none');
        axis([Plot.reglim -Plot.MaxEz Plot.MaxEz])
        %         axis equal
        view(Plot.pv)
        subplot(np,mp,2),s2(k) = surface(xC{k}+Reg.xoff{k},yB{k}+Reg.yoff{k},real(Hx{k}.'),'linestyle','none');
        axis([Plot.reglim -Plot.MaxH Plot.MaxH])
        %         axis equal
        view(Plot.pv)
        subplot(np,mp,3),s3(k) = surface(xB{k}+Reg.xoff{k},yC{k}+Reg.yoff{k},real(Hy{k}.'),'linestyle','none');
        axis([Plot.reglim -Plot.MaxH Plot.MaxH])
        %         axis equal
        view(Plot.pv)
    end
end



Ez_t = zeros(nSteps,ny{1});
Ez_r = zeros(nSteps,ny{1});

for i = 1:nSteps
    
    for k = 1:Reg.n
        
        % basic propagation with PML
        
        % Feb 5/17 tjs reversed update
        %          so H at n+1/2 (t+dt/2) and E at n+1 (t+dt)
        
        % Hx/y update at t + dt/2
        t = i*dt - dt/2; % n + 1/2 update
        th(i) = t;
        
        dEzdx{k} = Ez{k}(2:end,:) - Ez{k}(1:end-1,:);
        dEzdy{k} = Ez{k}(:,2:end) - Ez{k}(:,1:end-1);
        
        Hx{k}(:,2:end-1) = Da{k}(:,1:end-1).*Hx{k}(:,2:end-1) - Db{k}(:,1:end-1).*(dEzdy{k});
        Hy{k}(2:end-1,:) = Da{k}(1:end-1,:).*Hy{k}(2:end-1,:) + Db{k}(1:end-1,:).*(dEzdx{k});
        
        QExz{k} = bxHy{k}(1:end-1,:).*QExz{k} - cxHy{k}(1:end-1,:).*dEzdx{k};
        QEyz{k} = byHx{k}(:,1:end-1).*QEyz{k} - cyHx{k}(:,1:end-1).*dEzdy{k};
        
        Hx{k}(:,2:end-1) = Hx{k}(:,2:end-1) - dt./mu{k}(:,1:end-1).*QEyz{k};
        Hy{k}(2:end-1,:) = Hy{k}(2:end-1,:) + dt./mu{k}(1:end-1,:).*QExz{k};
        
        
        % boundary conditions on H
        
        if bc{k}.ym.type == 's'
            [Hx{k}(:,1), bc{k}.ym.paras] = bc{k}.ym.fct(xC{k},i,t,bc{k}.ym.paras);
        elseif bc{k}.ym.type == 'e'
            Hx{k}(:,1) = Hx{k}(:,2);
        elseif bc{k}.ym.type == 'm'
            Hx{k}(:,1) = 0;
        elseif bc{k}.ym.type == 'p'
            Hx{k}(:,1) = Da{k}(:,1).*Hx{k}(:,1) - Db{k}(:,1).*(Ez{k}(:,1) - Ez{k}(:,end));
        end
        
        if bc{k}.yp.type == 's'
            [Hx{k}(:,end),bc{k}.yp.paras] = bc{k}.yp.fct(xC{k},i,t,bc{k}.yp.paras);
        elseif bc{k}.yp.type == 'e'
            Hx{k}(:,end) = Hx{k}(:,end-1);
        elseif bc{k}.yp.type == 'm'
            Hx{k}(:,end) = 0;
        elseif bc{k}.yp.type == 'p'
            Hx{k}(:,end) = Hx{k}(:,1);
        end
        
        if bc{k}.xm.type == 's'
            [H, bc{k}.xm.paras] = bc{k}.xm.fct(yC{k},i,t,bc{k}.xm.paras);
            Hy{k}(1,:) = H.';
        elseif bc{k}.xm.type == 'e'
            Hy{k}(1,:) = Hy{k}(2,:);
        elseif bc{k}.xm.type == 'm'
            Hy{k}(1,:) = 0;
        elseif bc{k}.xm.type == 'p'
            Hy{k}(1,:) = Da{k}(1,:).*Hy{k}(:,1) + Db{k}(1,:).*(Ez{k}(1,:) - Ez{k}(end,:));
        end
        
        if bc{k}.xp.type == 's'
            [H,bc{k}.xp.paras] = bc{k}.xp.fct(yC{k},i,t,bc{k}.xp.paras);
            Hy{k}(end,:) = H.';
        elseif bc{k}.xp.type == 'e'
            Hy{k}(end,:) = Hy{k}(end-1,:);
        elseif bc{k}.xp.type == 'm'
            Hy{k}(end,:) = 0;
        elseif bc{k}.xp.type == 'p'
            Hy{k}(end,:) = Hy{k}(1,:);
        end
        
        % Ez update`
        t = i*dt; % n + 1 update
        te(i) = t;
        
        dHydx{k} = Hy{k}(2:end,:) - Hy{k}(1:end-1,:);
        dHxdy{k} = Hx{k}(:,2:end) - Hx{k}(:,1:end-1);
        
        %         Ez_p = Ez{k};
        Ez{k} = Ca{k}.*Ez{k} + Cb{k}.*(-dHxdy{k} + dHydx{k});
        
        QHxy{k} = bxEz{k}.*QHxy{k} - cxEz{k}.*dHydx{k}; %Update of the PML-Matrices
        QHyx{k} = byEz{k}.*QHyx{k} - cyEz{k}.*dHxdy{k};
        
        Ez{k} = Ez{k} - dt./epi{k}.*QHyx{k} + dt./epi{k}.*QHxy{k};
        
        if bc{k}.NumS > 0
            % Source 1
            Einc = bc{k}.s(1).fct(yC{k},i,t,bc{k}.s(1).paras);
            Eincp = bc{k}.s(1).fct(yC{k},i,t-ds/c_c/2,bc{k}.s(1).paras);
            
            Hy{k}(bc{k}.s(1).xpos,:) = Hy{k}(bc{k}.s(1).xpos,:) + dt/ds/c_mu_0*Einc;
            Ez{k}(bc{k}.s(1).xpos-1,:) = Ez{k}(bc{k}.s(1).xpos-1,:) - dt/ds/c_eps_0/c_eta_0*Eincp;

            % Source 2
            Einc = bc{k}.s(2).fct(yC{k},i,t,bc{k}.s(2).paras);
            Eincp = bc{k}.s(2).fct(yC{k},i,t-ds/c_c/2,bc{k}.s(2).paras);
            
            Hy{k}(bc{k}.s(2).xpos,:) = Hy{k}(bc{k}.s(2).xpos,:) + dt/ds/c_mu_0*Einc;
            Ez{k}(bc{k}.s(2).xpos-1,:) = Ez{k}(bc{k}.s(2).xpos-1,:) - dt/ds/c_eps_0/c_eta_0*Eincp;
            
        end
        
        
    end
    
    if mod(i,Plot.N) == 0 && movie
        
        %         i/nsteps*100
        
        delete(s1);
        delete(s2);
        delete(s3);
        
        for k = 1:Reg.n
            
            if Plot.pl
                
                colormap(hot);
                subplot(np,mp,1),plot(xC{k}+Reg.xoff{k},Ez{k}(:,Plot.y0),'-');hold off
                ylim([-Plot.MaxEz Plot.MaxEz]);
                subplot(np,mp,2),plot(xB{k}+Reg.xoff{k},Hy{k}(:,Plot.y0),'-');hold off
                ylim([-Plot.MaxH Plot.MaxH]);
            else
                subplot(np,mp,1),s1(k) = surface(xC{k}+Reg.xoff{k},yC{k}+Reg.yoff{k},real(Ez{k}.'),'linestyle','none');
                axis([Plot.reglim -Plot.MaxEz Plot.MaxEz])
                %            axis equal
                caxis([-Plot.MaxEz Plot.MaxEz]);
                view(Plot.pv)
                xlabel('x');
                ylabel('y');
                title('Ez');
                
                subplot(np,mp,2),s2(k) = surface(xC{k}+Reg.xoff{k},yB{k}+Reg.yoff{k},real(Hx{k}.'),'linestyle','none');
                axis([Plot.reglim -Plot.MaxH Plot.MaxH])
                %             axis equal
                caxis([-Plot.MaxH Plot.MaxH]);
                view(Plot.pv)
                xlabel('x');
                ylabel('y');
                title('Hx');
                
                subplot(np,mp,3),s3(k) = surface(xB{k}+Reg.xoff{k},yC{k}+Reg.yoff{k},real(Hy{k}.'),'linestyle','none');
                axis([Plot.reglim -Plot.MaxH Plot.MaxH])
                %             axis equal
                caxis([-Plot.MaxH Plot.MaxH]);
                view(Plot.pv)
                xlabel('x');
                ylabel('y');
                title('Hy');
            end
        end
        pause(0.01);
    end
    
    
    
    for p = 1 : length(simulationStopTimes)
        timeStop = simulationStopTimes(p);
        if t >= timeStop && (t-dt) < timeStop
            %             delete(s1);
            %             delete(s2);
            %             delete(s3);
            
            timeTitle = strcat(num2str(round(timeStop / 1e-12 * 100)/100), 'ps');
            savingTitle = strcat(mainTitle, '(');
            savingTitle = strcat(savingTitle, 'ny=', num2str(ny{1}));
            if surf.Mode == 1 && MetaSurf.dc > 0
                savingTitle = strcat(savingTitle, ',fm=', strrep(num2str(round(MetaSurf.wm/MetaSurf.w*100)/100),'.',','), 'f');
            end
            savingTitle = strcat(savingTitle, ',', timeTitle, ').mat');
            save(savingTitle, 'Ez', 'Hx', 'Hy');
            
            
        end
    end
    
end

% figure
% subplot(2,1,1),plot(Ezs0); hold on
% subplot(2,1,2),plot(Hys0); hold on
%
% subplot(2,1,1),plot(Ezss); hold on
% subplot(2,1,2),plot(Hyss); hold on


if Plot.pl && ~Plot.off
    colormap(hot);
    subplot(np,mp,1),plot(xC{k}+Reg.xoff{k},Ez{k}(:,Plot.y0),'-x');hold off
    ylim([-Plot.MaxEz Plot.MaxEz]);
    subplot(np,mp,2),plot(xB{k}+Reg.xoff{k},Hy{k}(:,Plot.y0),'-x');hold off
    ylim([-Plot.MaxH Plot.MaxH]);
    
elseif ~Plot.off
    subplot(np,mp,1),s1(1) = surface((xC{1}+Reg.xoff{1})/1e-6,(yC{1}+Reg.yoff{1})/1e-6,real(Ez{1}.'),'linestyle','none');
    axis([Plot.reglim/1e-6 -Plot.MaxEz Plot.MaxEz])
    %             axis equal
    caxis([-Plot.MaxEz Plot.MaxEz]);
    view(Plot.pv)
    set(gca, 'XTick', [0 25 50 75 100 125 150]);
    set(gca, 'YTick', [0 25 50]);
    % xlabel('x');
    % ylabel('y');
    title('Ez');
    
    subplot(np,mp,2),s2(1) = surface((xC{1}+Reg.xoff{1})/1e-6,(yB{1}+Reg.yoff{1})/1e-6,real(Hx{1}.'),'linestyle','none');
    axis([Plot.reglim/1e-6 -Plot.MaxH Plot.MaxH])
    %             axis equal
    caxis([-Plot.MaxH Plot.MaxH]);
    view(Plot.pv)
    set(gca, 'XTick', [0 25 50 75 100 125 150]);
    set(gca, 'YTick', [0 25 50]);
    % xlabel('x');
    ylabel('distance y (um)');
    title('Hx');
    
    subplot(np,mp,3),s3(1) = surface((xB{1}+Reg.xoff{1})/1e-6,(yC{1}+Reg.yoff{1})/1e-6,real(Hy{1}.'),'linestyle','none');
    axis([Plot.reglim/1e-6 -Plot.MaxH Plot.MaxH])
    %             axis equal
    caxis([-Plot.MaxH Plot.MaxH]);
    view(Plot.pv)
    set(gca, 'XTick', [0 25 50 75 100 125 150]);
    set(gca, 'YTick', [0 25 50]);
    xlabel('distance x (um)');
    % ylabel('y');
    title('Hy');
end

end

