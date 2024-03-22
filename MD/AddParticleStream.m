function [ output_args ] = AddParticleStream(num, x0, y0, PartAng, Type, Ep, Seper)
global AtomSpacing x y AtomType Vx Vy Mass0 Mass1 Mass2 nAtoms

if Type == 0
    Mass = Mass0;
elseif Type == 1
    Mass = Mass1;
else
    Mass = Mass2;
end

for p = 0:num - 1
    nAtoms = nAtoms + 1;
    x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng);
    y(nAtoms) = y0 * AtomSpacing - Seper * p * AtomSpacing * sin(PartAng);
    AtomType(nAtoms) = Type;
end

V = sqrt(2 * Ep / Mass);

for p = 1:num
    Vx(nAtoms - num + p) = V * cos(PartAng);
    Vy(nAtoms - num + p) = V * sin(PartAng);
end

end
