clear;
close all;

G1 = 1;
G2 = 1/2;
G3 = 1/10;
G4 = 1/0.1;
G5 = 1/1000;

C1 = 0.25;
L = 0.2;
a = 100;
omega = linspace(0, 100, 100);
Vout = zeros(1, 100);
GV5 = zeros(1, 100);

Vin = 1;

G = [1,0,0,0,0,0,0;
     G1,-G1,1,0,0,0,0;
     -G1,G1+G2,0,1,0,0,0;
     0,0,0,-1,G3,0,0;
     0,0,0,0,-a*G3,1,0;
     0,1,0,0,-1,0,0;
     0,0,0,0,0,-G4,G4+G5];
C =[0,0,0,0,0,0,0;
    C1,-C1,0,0,0,0,0;
    -C1,C1,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,-L,0,0,0;
    0,0,0,0,0,0,0];
F = [Vin;
    0;
    0;
    0;
    0;
    0;
    0];

VIN = linspace(-10,10,100);

for i = 1:100
    F(1) = VIN(i);
    Vop = G\F;
    V5(i) = Vop(7);
    V3(i) = Vop(5);

end

figure(1);
plot(VIN,V5)
title('Plot of V_O with a DC Sweep of the Input Voltage');
xlabel('V_{in}');
ylabel('V_O');

figure(2);
plot(VIN,V3)
title('Plot of V_3 with a DC Sweep of the Input Voltage');
xlabel('V_{in}');
ylabel('V_3');

for i = 1:100
    F(1) = VIN(i);
    Vop = (G + 1i * omega(i) * C)\F;
    V5(i) = Vop(7);
    V3(i) = Vop(5);
    Vout(i) = Vop(7);
    GV5(i) = Vout(i) / F(1);
end

% Vout as a function of omega
figure(3);
plot(omega, real(Vout));
title('Plot of V_O over \omega');
xlabel('Frequency (\omega)');
ylabel('V_O (v)');

% Magnitude in dB of Vout over Vin
figure(4);
semilogx(omega, 20*log10(real(GV5)));
title('Plot of the Magnitude of the V_O / V_{in} in dB');
xlabel('Frequency (\omega)');
ylabel('Magnitude (dB)');

omegaPartE = pi;

for i = 1 : 1000
    C(2, 1) = normrnd(C1, 0.05);
    C(2, 2) = normrnd(C1, 0.05);
    C(3, 1) = normrnd(C1, 0.05);
    C(3, 2) = normrnd(C1, 0.05);
    C(6, 4) = normrnd(L, 0.05);

    Vop = (G + 1i * omegaPartE * C)\F;
    VoutR(i) = Vop(7);
    GV5R(i) = VoutR(i) / F(1);
end

figure(5);
histogram(real(GV5R));
title('Histograms of the Gain With Normally Distributed Perturbations in C Matrix')
xlabel('V_O / V_{in}')
ylabel('Number')