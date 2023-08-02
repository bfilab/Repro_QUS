%% Plot Attenuation Functions
% Values are given in dB/cm/Mhz
% 5/12/2020 (THL): Created

figure(1); hold on;

% 18 um Phantom
f = linspace(10,40,100);                                                    % MHz
plot(f, 4*(0.0279 * f .^ 1.862 + 0),'r','LineWidth',2);                     % T15, T17, T35
f = linspace(10,20,100);                                                    % MHz
plot(f, 4*( 0.52 * f .^1.0 - 3.31),'--','LineWidth',2);                     % T15
f = linspace(12,22,100);                                                    % MHz
plot(f, 4*( 0.77 * f .^1.0 - 8.16),'--','LineWidth',2);                     % T17
f = linspace(30,40,100);                                                    % MHz
plot(f, 4*( 0.95 * f .^1.0 - 12.59),'--','LineWidth',2);                    % T35

f = linspace(20,30,100);                                                    % MHz
plot(f, 4*( 0.8 * f .^1.0 - 8.5),'LineWidth', 2);                             % 25 MHz, single element
f = linspace(15,25,100);                                                    % MHz
plot(f, 4*( 0.67 * f .^1.0 - 6.25),'LineWidth',2);                             % P6 annular array center element (~20 MHz) F/6?
f = linspace(5,15,100);                                                     % MHz
plot(f, 4*( 0.37 * f .^1.0 - 1.5),'LineWidth',2);                             % 10 MHz F/5 single element

legend('T15, T17, T35',...
    'T15',...
    'T17',...
    'T35',...
    '25 MHz, single element',...
    'P6 annular array center element (~20 MHz) F/6?',...
    '10 MHz F/5 single element');

xlabel('MHz');
ylabel('dB/cm');
title('18 um Phantom');

hold off;

% 15 um Phantom
figure(2);
hold on;

f = linspace(10,40,100);                                                    % MHz
plot(f, 4*(0.0254 * f .^ 1.842 + 0),'b','LineWidth',2);                     % T15, T17, T35
f = linspace(10,20,100);                                                    % MHz
plot(f, 4*( 0.43 * f .^1.0 - 2.62),'--','LineWidth',2);                     % T15
f = linspace(12,22,100);                                                    % MHz
plot(f, 4*( 0.65 * f .^1.0 - 6.79),'--','LineWidth',2);                     % T17
f = linspace(30,40,100);                                                    % MHz
plot(f, 4*( 0.81 * f .^1.0 - 10.74),'--','LineWidth',2);                    % T35

f = linspace(15,25,100);                                                    % MHz
plot(f, 4*( 0.55 * f .^1.0 - 4.5),'LineWidth',2);                             % 20 MHz, single element

legend('T15, T17, T35',...
    'T15',...
    'T17',...
    'T35',...
    '20 MHz, single element');

xlabel('MHz');
ylabel('dB/cm');
title('15 um Phantom');

hold off;

% 60 um phantom

figure(3);
hold on;

f = linspace(10,40,100);                                                    % MHz
plot(f, 4*(0.55 * f .^ 1.0 - 1.0),'g','LineWidth',2);

hold off;
