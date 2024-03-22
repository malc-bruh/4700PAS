
function [ output_args ] = PlotVars(c, Limits)
global Vx Vy L W x y Fx Fy C
global Phi nAtoms time Mass0 Mass1 Mass2 Pty0in Pty1in Pty2in
global LJEpsilon Phi0 PhiTot KETot MinX MaxX MinY MaxY
global T T0 T1 T2 ScaleV MarkerSize doPlotImage PlotCount
global PlotFig map im PlotSize ScaleF

if isempty(Limits)
    Limits = Limits;
end

if doPlotImage == 1
    if PlotCount == 0
        close all;
        set(0,'DefaultFigureWindowStyle','normal')
        PlotFig = figure('Position', PlotSize);
        PlotCount = 1;
    else
        PlotCount = PlotCount+1;
    end
end

V2 = Vx.*Vx + Vy.*Vy;
MaxV = max(sqrt(V2));
dx = MaxX - MinX;
dy =  MaxY - MinY;

if c==0 && ~ScaleV
    ScaleV = sqrt(dx*dx + dy*dy) / MaxV * 0.15;
end

subplot(3, 2, 1), plot(x(Pty0in), y(Pty0in), 'bo', 'markers',...
    MarkerSize,'MarkerFaceColor', 'b');
hold on
subplot(3, 2, 1), plot(x(Pty1in), y(Pty1in), 'go', 'markers',...
    MarkerSize,'MarkerFaceColor', 'g');
subplot(3, 2, 1), plot(x(Pty2in), y(Pty2in), 'ro', 'markers',...
    MarkerSize,'MarkerFaceColor', 'r');
quiver(x, y, Fx * ScaleF, Fy * ScaleF, 0, 'm', 'linewidth', 2);
hold off
axis(Limits);
title('Atoms and Forces')
xlabel('X')
ylabel('Y')

subplot(3, 2, 2), plot(time, KETot, 'linewidth', 2);
hold on
subplot(3, 2, 2), plot(time, PhiTot + LJEpsilon * nAtoms, 'g', 'linewidth',2);
subplot(3, 2, 2), plot(time, PhiTot + KETot + LJEpsilon * nAtoms, 'k', 'linewidth',2);
%     axis([0 TStop 0e-22 3e-23]);
xlabel('time')
ylabel('Energy')
title('Total, Kinetic and Potential Energy')

subplot(3,2,3), plot(x(Pty0in), y(Pty0in), 'bo', 'markers',...
    MarkerSize,'MarkerFaceColor','b');
hold on
subplot(3,2,3), plot(x(Pty1in), y(Pty1in), 'go', 'markers',...
    MarkerSize,'MarkerFaceColor','g');
subplot(3,2,3), plot(x(Pty2in), y(Pty2in), 'ro', 'markers',...
    MarkerSize,'MarkerFaceColor','r');
subplot(3,2,3),quiver(x, y, Vx * ScaleV, Vy * ScaleV, 0, 'r', 'linewidth', 2);
hold off
axis(Limits);
title('Atoms and Velocities')
xlabel('X')
ylabel('Y')

AE(Pty0in) = 1 / 2 * Mass0 * V2(Pty0in) + Phi(Pty0in) - Phi0;
AE(Pty1in) = 1 / 2 * Mass1 * V2(Pty1in) + Phi(Pty1in) - Phi0;
AE(Pty2in) = 1 / 2 * Mass1 * V2(Pty2in) + Phi(Pty2in) - Phi0;
% if AddParticle, AE(end) = 0;end
% AE = Phi;
subplot(3, 2, 4),scatter3(x,y,AE,ones(1,nAtoms)*300,AE,'fill');
axis(Limits);
view(2);
title('Atom Energy (KE + Pot)')

subplot(3, 2, 5), plot(x(Pty0in), y(Pty0in), 'bo', 'markers',...
    MarkerSize, 'MarkerFaceColor', 'b');
hold on
subplot(3, 2, 5), plot(x(Pty1in), y(Pty1in), 'go', 'markers',...
    MarkerSize, 'MarkerFaceColor', 'g');
subplot(3, 2, 5), plot(x(Pty2in), y(Pty2in), 'ro', 'markers',...
    MarkerSize, 'MarkerFaceColor', 'r');
hold off
title('Atoms (All)')
xlabel('X')
ylabel('Y')

subplot(3, 2, 6), plot(time, T, 'k', 'linewidth', 2);
hold on
subplot(3, 2, 6), plot(time, T0, 'b', 'linewidth', 2);
subplot(3, 2, 6), plot(time, T1, 'g', 'linewidth', 2);
%
xlabel('time')
ylabel('Temperature')
title('Temperature: Each type and all')

% if doPlotImage
%     if PlotCount == 1
%         f = getframe(PlotFig);
%         [im, map] = rgb2ind(f.cdata, 256, 'nodither');
%         im(1, 1, 1, 2) = 0;
% 
%     else
%         f = getframe(PlotFig);
%         im(:, :, 1, PlotCount) = rgb2ind(f.cdata, map, 'nodither');
%     end
% end

end
