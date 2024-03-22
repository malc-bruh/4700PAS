
function [ output_args ] = PlotVars(c, Limits)
global Vx Vy L W x y Fx Fy C
global Phi nAtoms time Mass0 Mass1 Mass2 Pty0in Pty1in Pty2in
global LJEpsilon Phi0 PhiTot KETot MinX MaxX MinY MaxY
global T T0 T1 T2 ScaleV MarkerSize doPlotImage PlotCount
global PlotFig map im PlotSize ScaleF

if isempty(Limits)
    Limits = Limits;
end

plot(x(Pty0in), y(Pty0in), 'bo', 'markers',...
    MarkerSize,'MarkerFaceColor', 'b');
hold on
plot(x(Pty1in), y(Pty1in), 'go', 'markers',...
    MarkerSize,'MarkerFaceColor', 'g');
plot(x(Pty2in), y(Pty2in), 'go', 'markers',...
    MarkerSize,'MarkerFaceColor', 'r');
quiver(x, y, Fx * ScaleF, Fy * ScaleF, 0, 'm', 'linewidth', 2);
hold off
axis(Limits);
title('Atoms and Forces')
xlabel('X')
ylabel('Y')


end
