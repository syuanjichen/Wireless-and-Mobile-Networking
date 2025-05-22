%% Task 1. Radio Propagation: Only Path Loss is Considered
temp = 27 + 273.15; % Temperature (Kelvin)
bandwidth = 10e6; % Channel Bandwidth (Hz)
power_base = 33; % Power of Base station (dBm, the dimension of mW)
transmitter_gain = 14; % Transmitter Gain (dB)
receiver_gain = 14; % Receiver Gain (dB)
base_height = 50 + 1.5; % Base station height from the ground (m)
mobile_height = 1.5; % Mobile Device height from the ground (m)
max_distance = 1000; % The maximum distance considered (m)
n = 1000; % The number of samples
k = physconst("Boltzmann"); % Boltzmann's Constant (J/K)
thermal_noise_power = (k * temp * bandwidth) / 10^(-3); % Calculate Thermal Noise Power (mW)
rng('default')

% Formula: P_R = P_T * g(d) * 10^(x/10) * G_T * G_R
% Take logarithm on both sides yields 
% log(P_R) = log(P_T) + log(g(d)) + x/10 + log(G_T) + log(G_R)

distance = linspace(1, max_distance, n); % Distance samples
path_loss = zeros(1, n); % For calculating path loss g(d)
power_received = zeros(1, n); % Store received power with only path loss
power_received_shadowing = zeros(1, n); % Store received power with path loss and shadowing
SINR = zeros(1, n); % Store SINR with only path loss
SINR_shadowing = zeros(1, n); % Store SINR with path loss and shadowing

for j = 1 : n
    path_loss(j) = ((base_height * mobile_height)^2) / ((distance(j))^4);
    power_received(j) = power_base + 10*log10(path_loss(j)) + transmitter_gain + receiver_gain;
    SINR(j) = power_received(j) - 10*log10(thermal_noise_power);
end

plot(distance, power_received) % Plot Received Power with respect to distance
title("Received Power to Distance using the Two-Ray-Ground Model")
xlabel("Distance (m)")
ylabel("Received Power (dBm)")

figure;
plot(distance, SINR) % Plot SINR with respect to distance
title("SINR to Distance using the Two-Ray-Ground Model")
xlabel("Distance (m)")
ylabel("SINR (dB)")

%% Task 2. Radio Propagation: Both Path Loss and Shadowing are Considered
x_mean = 0; % Mean of x, in dB
x_stdev = 6; % Standard Deviation of x, in dB

for j = 1:n
    x = normrnd(x_mean, x_stdev);
    power_received_shadowing(j) = power_received(j) + x;
    SINR_shadowing(j) = power_received_shadowing(j) - 10*log10(thermal_noise_power);
end

figure;
plot(distance, power_received_shadowing) % Plot Received Power with respect to distance
title("Received Power to Distance, with Two-Ray-Ground Path Loss and Log-Normal Shadowing")
xlabel("Distance (m)")
ylabel("Received Power (dBm)")

figure;
plot(distance, SINR_shadowing) % Plot SINR with respect to distance
title("SINR to Distance, with Two-Ray-Ground Path Loss and Log-Normal Shadowing")
xlabel("Distance (m)")
ylabel("SINR (dB)")
