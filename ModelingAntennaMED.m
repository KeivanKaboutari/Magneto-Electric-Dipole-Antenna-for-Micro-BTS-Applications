clear all;
% close all;
clc;

%%

% % resonance frequency [Hz]
% fr = 905e6;

% % Angular frequency [rad. Hz]
% Omega = 2 * pi * fr;

% frequency range (GHz)
StartingFreq = 0.5;
EndingFreq = 1.5;
Steps = 0.005;
FrRange = StartingFreq : Steps : EndingFreq;

% Input impedance of a loaded transmission line
% BL = Beta * l
Zin = @(Z_L, Zw, BL) Zw * (Z_L + 1i * Zw .* tan(BL)) / (Zw + 1i * Z_L .* tan(BL));

% Sine and cosine integrals (for Schelkunoff's biconical antenna model):
Si = @(Z) sinint(Z);
Ci = @(Z) cosint(Z);
Cin = @(Z) 0.57722 + log(Z) - cosint(Z);

% Characteristic impedance of the biconical anteanna in Schelkunoff's model (Alpha is the cone half-angle, in rad)
Rho = @(Alpha) 120 * log(cot(Alpha / 2));

% ???????????????????????????????????????????????
Kl = 2e9 * pi * 0.1 * FrRange / 3e8;

% Effective cone half-angle
Alpha = 20.0 * pi / 180.0;

% Schelkunoff's model for the input impedance of a biconical antenna (Kl = (Omega / c) * dipole_half_length)
Ra = 60 * Cin(2 * Kl) + 30 * (2 * Cin(2 * Kl) - Cin(4 * Kl)) .* cos(2 * Kl) + 30 * (Si(4 * Kl) - 2 * Si(2 * Kl)) .* sin(2 * Kl);
Xa = 60 * Si(2 * Kl) - 30 * (Cin(4 * Kl) - log(4.0)) .* sin(2 * Kl) - 30 * Si(4 * Kl) .* cos(2 * Kl);

% Impedance of the dipole antenna
% Calculate real and imaginary part of the Zdip
% Correction factor takes into account the effect of the ground
Zdip = @(Kl, Alpha, Ra, Xa, Correction) Zin(Correction * Ra + 1i * Xa, Rho(Alpha), Kl - 0.5 * pi);
RealZd = zeros(1, length(FrRange));
ImaginaryZd = zeros(1, length(FrRange));
for Counter = 1 : 1 : length(Kl)
    % Calculate real part of impedance
    RealZd(Counter) = real(Zdip(Kl(Counter), Alpha, Ra(Counter), Xa(Counter), 1));
    % Calculate imaginary part of impedance
    ImaginaryZd(Counter) = imag(Zdip(Kl(Counter), Alpha, Ra(Counter), Xa(Counter), 1));
end

% for Counter = 1 : 1 : length(Kl)
%     Zd = Zin(Correction * Ra(Counter) + 1i * Xa(Counter), Rho(Alpha), Kl(Counter) - 0.5 * pi);
%     % Calculate real part of impedance
%     RealZd(Counter) = real(Zd);
%     % Calculate imaginary part of impedance
%     ImaginaryZd(Counter) = imag(Zd);
% end

% Biconical antenna Impedance plot
% Set axes color for figures
left_color = [0 0 0];
right_color = [0 0 0];
Fig = figure('Name', 'Impedance of the dipole (Ohm.) using Schelkunoff''s model', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 24, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
plot(FrRange, RealZd, 'LineWidth', 3, 'Color', 'black');
hold on;
plot(FrRange, ImaginaryZd, 'LineWidth', 3, 'Color', 'blue', 'LineStyle', '-.');
xlim([-inf inf]);
ylim([-inf inf]);
grid on;
% xlabel('Frequency (GHz)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 20, 'Interpreter','latex');
% ylabel('Impedance ($\Omega$)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 20, 'Interpreter','latex');
xlabel('Frequency (GHz)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 24);
ylabel('Impedance (\Omega)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 24);
axis tight;
hold off;
% LegProp = legend('Real part', 'Imaginary part', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'white', 'Location', 'best', 'FontSize' , 20, 'Interpreter','latex');
LegProp = legend('Real part', 'Imaginary part', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'white', 'Location', 'best', 'FontSize' , 24);

%% ????????????????

% Width of the microstripline in [m]
w = 8e-3;

% Dielectric hight in [m]
h = 3.5e-3;

% WperH = W / H
WperH = w / h;

% Relative dielectric of the substrate
e_r =  1;

% Approximate formula for the microstrip line impedance (assuming epsr = 1, WperH > 1)
if (WperH) > 1
    % Effective dielectric of the substrate
    e_eff = (e_r + 1) / 2 + (e_r - 1) / 2 * sqrt(1 + 12 * h / w);
    % Characteristic impedance of the microstripline [Ohm.]
    Zw = @(WperH) 120 * pi / (sqrt(e_eff) * (WperH + 1.393 + 2 / 3 * log(WperH + 1.444)));
else
    % Effective dielectric of the substrate
    e_eff = (e_r + 1) / 2 + (e_r - 1) / 2 * (1 / sqrt(1 + 12 * h / w) + 0.04 * power(1 - WperH, 2));
    % Characteristic impedance of the microstripline [Ohm.]
    Zw = @(WperH) 60 / sqrt(e_eff) * log(8 * h / w + 0.25 * WperH);
end

% % Speed of light [m / s]
% c = physconst('LightSpeed');

% % Wavenumber
% Beta = 2 * pi * fr * sqrt(e_r) / c;

% ???????????????????
k0 = 2e9 * pi * FrRange / 3e8;

% Characteristic impedance of the antenna feed [Ohm.]
Zw0 = 50.0;

% Characteristic impedance of the 1st microstrip line [Ohm.]
Zw1 = Zw(WperH);
% Length of the 1st microstrip [mm]
L1 = 93e-3;

% Characteristic impedance of the 2nd microstrip line [Ohm.]
Zw2 = Zw(WperH);
% Length of the 2nd microstrip [mm]
L2 = 81e-3;

% Width of the microstripline for loop in [m]
wl = 21e-3;
% Dielectric hight in for loop [m]
hl = 21e-3 / 2;
% WperH = W / H
WperHl = wl / hl;
% Characteristic impedance of the loop
Zw_loop = 2 * Zw(WperHl);
% Length of the loop
l_loop = L1;
% Width of the loop
h_loop = 21e-3;

% Effective halflength of the biconical antenna [mm]
% Leff = 170e-3;
Leff = 195e-3;
% Effective cone half-angle
% Alpha = 25.0 * pi / 180;
Alpha = 25.0 * pi / 180;

Rrad_loop = zeros(1, length(FrRange));
% Correction factor takes into account the effect of the ground
Correction = zeros(1, length(FrRange));
% Input impedance of the complete antenna obtained from the equivalent network
Z1 = zeros(1, length(FrRange));
Za = zeros(1, length(FrRange));
for Counter = 1 : 1 : length(Kl)
    Rrad_loop(Counter) = 2 * (80 * power(k0(Counter), 2) * power(k0(Counter) * l_loop * h_loop * 2 / pi, 2));
    
    Correction(Counter) = 2 * power(sin(k0(Counter) * L1), 2);
    
    Ra(Counter) = 60 * Cin(2 * k0(Counter) * Leff) + 30 * (2 * Cin(2 * k0(Counter) * Leff) - Cin(4 * k0(Counter) * Leff)) .* cos(2 * k0(Counter) * Leff) + 30 * (Si(4 * k0(Counter) * Leff) - 2 * Si(2 * k0(Counter) * Leff)) .* sin(2 * k0(Counter) * Leff);
    Xa(Counter) = 60 * Si(2 * k0(Counter) * Leff) - 30 * (Cin(4 * k0(Counter) * Leff) - log(4.0)) .* sin(2 * k0(Counter) * Leff) - 30 * Si(4 * k0(Counter) * Leff) .* cos(2 * k0(Counter) * Leff);
    
    Z1(Counter) = 1 / (1 / Zdip(k0(Counter) * Leff, Alpha, Ra(Counter), Xa(Counter), Correction(Counter)) + 1 / (Rrad_loop(Counter) + 1i * Zw_loop * tan(k0(Counter) * l_loop))) - 1i * Zw2 * cot(k0(Counter) * L2);
    
    Za(Counter) = Zin(Z1(Counter), Zw1, k0(Counter) * L1);
end

% Reflection Coefficient (Gamma) is equal to (Zl - Z0) / (Zl + Z0), then, Standing Wave Ratio (SWR) is
Gamma = @(Zl, Z0) (Zl - Z0) / (Zl + Z0);
SWR = @(Zl, Z0) (1 + abs(Gamma(Zl, Z0))) / (1 - abs(Gamma(Zl, Z0)));

VSWR = zeros(1, length(FrRange));
for Counter = 1 : 1 : length(FrRange)
    VSWR(Counter) = SWR(Za(Counter), Zw0);
end

Fig = figure('Name', 'Impedance of the dipole (Ohm.) using Schelkunoff''s model', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 24, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
plot(FrRange, real(Za), 'LineWidth', 3, 'Color', 'black');
hold on;
plot(FrRange, imag(Za), 'LineWidth', 3, 'Color', 'blue', 'LineStyle', '-.');
plot(FrRange, repmat(50, length(FrRange)), 'LineWidth', 3, 'Color', [0.6350 0.0780 0.1840], 'LineStyle', ':');
xlim([-inf inf]);
ylim([-inf inf]);
grid on;
% xlabel('Frequency (GHz)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 20, 'Interpreter','latex');
% ylabel('Impedance ($\Omega$)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 20, 'Interpreter','latex');
xlabel('Frequency (GHz)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 24);
ylabel('Impedance (\Omega)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 24);
axis tight;
hold off;
% LegProp = legend('Real part', 'Imaginary part', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'white', 'Location', 'best', 'FontSize' , 20, 'Interpreter','latex');
LegProp = legend('Real part', 'Imaginary part', '50 \Omega', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'white', 'Location', 'best', 'FontSize' , 24);

% Plot VSWR
Fig = figure('Name', 'VSWR', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 24, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
plot(FrRange, VSWR, 'LineWidth', 3, 'Color', 'black');
xlim([0.4 1.4]);
ylim([-0 14]);
grid on;
xlabel('Frequency (GHz)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 24);
ylabel('VSWR', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 24);
% axis tight;

RL = -20 * log10((VSWR - 1) ./ (VSWR + 1));
S11Lin = (VSWR - 1) ./ (VSWR + 1);
S11dB = 20 * log10(S11Lin);

% Plot S11
Fig = figure('Name', 'S11', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 24, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
plot(FrRange, S11dB, 'LineWidth', 3, 'Color', 'black');
xlim([0.4 1.4]);
ylim([-40 0]);
grid on;
xlabel('Frequency (GHz)', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 24);
ylabel('S_{\rm 11} [dB]', 'FontWeight' , 'bold' , 'FontName' , 'Times New Roman' , 'Color' , 'black', 'FontSize' , 24);
% axis tight;
