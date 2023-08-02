%% Plot Fluid-Filled Sphere Correlation Function
% 11/03/2020 (THL): Created

% Reference: Insana, M. F., R. Wagner, D. Brown and T. Hall (1990).
% "Describing small-scale structure in random media using pulse-echo
% ultrasound." The Journal of the Acoustical Society of America 87:
% 179-192.


r = -100:100;
a = 50;
lambda = 100;
k = (2*pi)/lambda;
% max_a = max(abs(delta_r))/2;

b = 1 - ( (3*abs(r))./(4*a) ) .* ( (abs(r).^3)./(16*a^3) );
figure(1);
plot(r,b); axis tight;
xlabel('\Delta r');

x = b.*sin(2*k*r).*r;
figure(2);
plot(r,x); axis tight;
xlabel('\Delta r');