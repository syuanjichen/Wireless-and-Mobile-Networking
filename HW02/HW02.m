clear all;
clc;

%% Task 1. Downlink
n_base_station = 19; % Number of base stations
temperature = 27 + 273.15; % Temperature in Kelvin
isd = 500; % Inter-site Distance
bandwidth = 10e6; % Channel Bandwidth
power_base_station = 33; % Power of base station in dBm
power_mobile = 23; % Power of mobile device in dBm
gain_transmitter = 14; % Transmitter gain in dB
gain_receiver = 14; % Receiver gain in dB
height_base_station = 1.5 + 50; % Height of each base station
height_mobile = 1.5; % Height of Mobile Device
rng default

n_mobile_device = 50; % Number of mobile devices in the CBS
CBS = [0, 0]; % Location of Central Base Station
location_mobile_device_x = 500/sqrt(3) * 2 * (rand(n_mobile_device, 1) - 0.5);
location_mobile_device_y = 250 * 2 * (rand(n_mobile_device, 1) - 0.5);

% Scatter: BS and mobile devices
hold on
scatter(0, 0, 'red', 'filled')
scatter(location_mobile_device_x, location_mobile_device_y, 'blue', 'filled')
hold off
xlabel("x (m)")
ylabel("y (m)")
title("Base Station (0, 0) and the Distribution of Mobile Devices")
legend('Base Station', 'Mobile Device')

% For computing the received power
dist_mobile_device = zeros(n_mobile_device, 1);
path_loss = zeros(n_mobile_device, 1);
power_received_mobile = zeros(n_mobile_device, 1);

for i = 1:n_mobile_device
    dist_mobile_device(i) = sqrt((location_mobile_device_x(i))^2 + (location_mobile_device_y(i))^2);
    path_loss(i) = (height_base_station * height_mobile)^2 / (dist_mobile_device(i))^4;
    power_received_mobile(i) = power_base_station + gain_transmitter + 10*log(path_loss(i)) + gain_receiver;
end

figure;
scatter(dist_mobile_device, power_received_mobile, 'filled')
xlabel("Distance (m)")
ylabel("Received Power of Mobile Device (dBm)")
title("Received Power of Mobile Device and Distance")

% For Calculating SINR
power_thermal = physconst("Boltzmann") * temperature * bandwidth;
path_loss_near1 = zeros(n_mobile_device, 1);
path_loss_D1 = zeros(n_mobile_device, 1);
path_loss_far1 = zeros(n_mobile_device, 1);
path_loss_nearsq3 = zeros(n_mobile_device, 1);
path_loss_Dsq3 = zeros(n_mobile_device, 1);
path_loss_farsq3 = zeros(n_mobile_device, 1);
path_loss_near2 = zeros(n_mobile_device, 1);
path_loss_D2 = zeros(n_mobile_device, 1);
path_loss_far2 = zeros(n_mobile_device, 1);
interference_mobile = zeros(n_mobile_device, 1);
SINR_mobile = zeros(n_mobile_device, 1);

for i = 1:n_mobile_device
    path_loss_near1(i) = (height_base_station * height_mobile)^2 / (isd - dist_mobile_device(i))^4;
    path_loss_D1(i) = (height_base_station * height_mobile)^2 / (isd)^4;
    path_loss_far1(i) = (height_base_station * height_mobile)^2 / (isd + dist_mobile_device(i))^4;
    path_loss_nearsq3(i) = (height_base_station * height_mobile)^2 / (isd * sqrt(3) - dist_mobile_device(i))^4;
    path_loss_Dsq3(i) = (height_base_station * height_mobile)^2 / (isd * sqrt(3))^4;
    path_loss_farsq3(i) = (height_base_station * height_mobile)^2 / (isd * sqrt(3) + dist_mobile_device(i))^4;
    path_loss_near2(i) = (height_base_station * height_mobile)^2 / (isd * 2- dist_mobile_device(i))^4;
    path_loss_D2(i) = (height_base_station * height_mobile)^2 / (isd * 2)^4;
    path_loss_far2(i) = (height_base_station * height_mobile)^2 / (isd * 2+ dist_mobile_device(i))^4;
    interference_mobile(i) = 2 * (10^((power_base_station + gain_transmitter + 10*log(path_loss_near1(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_D1(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_far1(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_nearsq3(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_Dsq3(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_farsq3(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_near2(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_D2(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_far2(i)) + gain_receiver)/10));
    SINR_mobile(i) = power_received_mobile(i) - 10*log((power_thermal/1e-3) + interference_mobile(i));
end

figure;
scatter(dist_mobile_device, SINR_mobile, 'filled')
xlabel("Distance (m)")
ylabel("SINR of Mobile Device(dB)")
title("SINR of Mobile Device and Distance")

%% Task 2. Uplink

% Scatter: BS and mobile devices
figure;
hold on
scatter(0, 0, 'red', 'filled')
scatter(location_mobile_device_x, location_mobile_device_y, 'blue', 'filled')
hold off
xlabel("x (m)")
ylabel("y (m)")
title("Base Station (0, 0) and the Distribution of Mobile Devices")
legend('Base Station', 'Mobile Device')

% For calculating the received power
power_received_bs = zeros(n_mobile_device, 1);

for i = 1:n_mobile_device
    power_received_bs(i) = power_mobile + gain_transmitter + 10*log(path_loss(i)) + gain_receiver;
end

figure;
scatter(dist_mobile_device, power_received_bs, 'filled')
xlabel("Distance (m)")
ylabel("Received Power of Base Station (dBm)")
title("Received Power of Base Station and Distance")

% For calculating SINR
SINR_bs = zeros(n_mobile_device, 1);
power_to_bs = 0;

for i = 1:n_mobile_device
    power_to_bs = power_to_bs + 10^(power_received_bs(i)/10);
end

for i = 1:n_mobile_device
    SINR_bs(i) = power_received_bs(i) - 10*log((power_thermal/1e-3) + power_to_bs - 10^(power_received_bs(i)/10));
end

figure;
scatter(dist_mobile_device, SINR_bs, 'filled')
xlabel("Distance (m)")
ylabel("SINR of Base Station(dB)")
title("SINR of Base Station and Distance")
