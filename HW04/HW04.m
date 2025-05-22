clear all;
clc;

%% Task 1. Downlink with Constant Bit Rate
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
BS_traffic_buffer = 6e6; % Size of BS traffic buffer
time_simulation = 1000; % Simulation Time

CBR_low = 5e6; % Low traffic load
CBR_medium = 10e6; % Medium traffic load
CBR_high = 20e6; % High traffic load
CBR = [CBR_low, CBR_medium, CBR_high];

rng shuffle

n_mobile_device = 50; % Number of mobile devices in the CBS
BS = zeros(2, n_base_station);
BS(:, 1) = [0, 0]; % Location of Central Base Station
for i = 2:7
    BS(:, i) = 500 * [cos((pi * i / 3) - pi/2), sin((pi*i/3)-pi/2)];
end
CBS = [BS(1, 1), BS(2, 1)]; % Location of Central Base Station
mobile_device_dist = 500/sqrt(3) * (rand(n_mobile_device, 1));
mobile_device_angle = 2 * pi * (rand(n_mobile_device, 1));
mobile_device_x = mobile_device_dist .* cos(mobile_device_angle);
mobile_device_y = mobile_device_dist .* sin(mobile_device_angle);

% 1.1 Plot Central BS and 50 Mobile Devices
hold on
axis equal
pgon = nsidedpoly(6, 'Center', CBS, 'SideLength', 500/sqrt(3));
plot(pgon)

scatter(0, 0, 'red', 'filled')
scatter(mobile_device_x, mobile_device_y, 'blue', 'filled')
hold off
xlabel("x (m)")
ylabel("y (m)")
title("Base Station (0, 0) and the Distribution of Mobile Devices")
legend('Cell Range', 'Base Station', 'Mobile Device')

% 1.2 Plot Shannon Capacity

% For computing the received power;
path_loss = zeros(n_mobile_device, 1);
mobile_received_power = zeros(n_mobile_device, 1);
mobile_capacity = zeros(n_mobile_device, 1);

for i = 1:n_mobile_device
    path_loss(i) = (height_base_station * height_mobile)^2 / (mobile_device_dist(i))^4;
    mobile_received_power(i) = power_base_station + gain_transmitter + 10*log(path_loss(i)) + gain_receiver; % in dBm
end

%{
figure;
scatter(mobile_device_dist, mobile_received_power, 'filled')
xlabel("Distance (m)")
ylabel("Received Power of Mobile Device (dBm)")
title("Received Power of Mobile Device and Distance")
%}

% For Calculating SINR
power_thermal = physconst("Boltzmann") * temperature * bandwidth;
path_loss_neighbor = zeros(n_mobile_device, 6);
dist_to_neighbor = zeros(n_mobile_device, 6);
% path_loss_near1 = zeros(n_mobile_device, 1);
% path_loss_D1 = zeros(n_mobile_device, 1);
% path_loss_far1 = zeros(n_mobile_device, 1);
% path_loss_nearsq3 = zeros(n_mobile_device, 1);
% path_loss_Dsq3 = zeros(n_mobile_device, 1);
% path_loss_farsq3 = zeros(n_mobile_device, 1);
% path_loss_near2 = zeros(n_mobile_device, 1);
% path_loss_D2 = zeros(n_mobile_device, 1);
% path_loss_far2 = zeros(n_mobile_device, 1);
interference_mobile = zeros(n_mobile_device, 1);
mobile_SINR = zeros(n_mobile_device, 1); % in dB

%{
for i = 1:n_mobile_device
    path_loss_near1(i) = (height_base_station * height_mobile)^2 / (isd - mobile_device_dist(i))^4;
    path_loss_D1(i) = (height_base_station * height_mobile)^2 / (isd)^4;
    path_loss_far1(i) = (height_base_station * height_mobile)^2 / (isd + mobile_device_dist(i))^4;
    path_loss_nearsq3(i) = (height_base_station * height_mobile)^2 / (isd * sqrt(3) - mobile_device_dist(i))^4;
    path_loss_Dsq3(i) = (height_base_station * height_mobile)^2 / (isd * sqrt(3))^4;
    path_loss_farsq3(i) = (height_base_station * height_mobile)^2 / (isd * sqrt(3) + mobile_device_dist(i))^4;
    path_loss_near2(i) = (height_base_station * height_mobile)^2 / (isd * 2 - mobile_device_dist(i))^4;
    path_loss_D2(i) = (height_base_station * height_mobile)^2 / (isd * 2)^4;
    path_loss_far2(i) = (height_base_station * height_mobile)^2 / (isd * 2 + mobile_device_dist(i))^4;
    interference_mobile(i) = 2 * (10^((power_base_station + gain_transmitter + 10*log(path_loss_near1(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_D1(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_far1(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_nearsq3(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_Dsq3(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_farsq3(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_near2(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_D2(i)) + gain_receiver)/10) + ...
        10^((power_base_station + gain_transmitter + 10*log(path_loss_far2(i)) + gain_receiver)/10));
    mobile_SINR(i) = mobile_received_power(i) - 10*log((power_thermal/1e-3) + interference_mobile(i));
    mobile_capacity(i) = bandwidth * log2(1 + 10^(mobile_SINR(i)/10));
end
%}

for i = 1:n_mobile_device
    for j = 2:7
        dist_to_neighbor(i, j - 1) = sqrt((mobile_device_x(i) - BS(1, j))^2 + (mobile_device_y(i) - BS(2, j))^2);
        path_loss_neighbor(i, j - 1) = ((height_base_station * height_mobile)^2 / (dist_to_neighbor(i, j - 1))^4);
        interference_mobile(i) = interference_mobile(i) + 10^((power_base_station + gain_transmitter + 10*log(path_loss_neighbor(i, j - 1)) + gain_receiver)/10);
    end
    mobile_SINR(i) = mobile_received_power(i) - 10*log((power_thermal/1e-3) + interference_mobile(i));
    mobile_capacity(i) = bandwidth * log2(1 + 10^(mobile_SINR(i)/10));
end

figure;
scatter(mobile_device_dist, mobile_capacity, 'filled')
xlabel("Distance (m)")
ylabel("Shannon Capacity of Mobile Device (bits / s)")
title("SINR of Mobile Device and Distance")

% 1.3 CBR transmission
CBR_bitloss = zeros(3, 1);

% Low traffic
for k = 1:3
    mobile_bits = zeros(n_mobile_device, 1);
    mobile_buffer = zeros(n_mobile_device, 1);
    mobile_loss = zeros(n_mobile_device, 1);
    mobile_bitlossprob = zeros(n_mobile_device, 1);
    for t = 1:time_simulation
        for i = 1:n_mobile_device
            mobile_bits(i) = 0;
            if mobile_bits(i) > mobile_buffer(i)
                mobile_bits(i) = mobile_bits(i) + mobile_buffer(i);
                mobile_buffer(i) = 0;
            else
                mobile_buffer(i) = mobile_buffer(i) - mobile_capacity(i);
                mobile_bits(i) = mobile_capacity(i);
            end
    
            if mobile_bits(i) + CBR(k) < mobile_capacity(i)
                mobile_bits(i) = mobile_bits(i) + CBR_low;
            else
                remain = CBR(k) - (mobile_capacity(i) - mobile_bits(i));
                mobile_bits(i) = mobile_capacity(i);
                if remain < BS_traffic_buffer
                    mobile_buffer(i) = remain;
                else
                    mobile_loss(i) = mobile_loss(i) + (remain - BS_traffic_buffer);
                    mobile_buffer(i) = BS_traffic_buffer;
                end
            end
        end
    end
    
    for i = 1:n_mobile_device
        mobile_bitlossprob(i) = mobile_loss(i) / (CBR(k) * time_simulation);
    end

    CBR_bitloss(k) = mobile_bitlossprob(i);
end

figure;
CBR_x = ["Low", "Medium", "High"];
bar(CBR_x, CBR_bitloss, 0.5)
title("Comparison of CBR low, medium, high traffic")
xlabel("Type")
ylabel("Bit Loss Probability")
