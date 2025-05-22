clear all;
clc;

%% Task 1. 1 Moving Mobile Device, Downlink
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
rng shuffle

minSpeed = 1;
maxSpeed = 15;
minT = 1;
maxT = 6;

MS_direction = 2*pi*rand;
MS_velocity = (maxSpeed - minSpeed)*rand + minSpeed;
MS_duration = (maxT - minT)*rand + minT;
MS_pos_init = [250, 0];

BS = zeros(2, n_base_station);
BS(:, 1) = [0, 0]; % Location of Central Base Station
hold on
axis equal
pgon = nsidedpoly(6, 'Center', [0, 0], 'SideLength', 500/sqrt(3));
plot(pgon)
txt = "1";
text(50, 0, txt)

for i = 2:7
    BS(:, i) = 500 * [cos((pi * i / 3) - pi/2), sin((pi*i/3)-pi/2)];
    pgon = nsidedpoly(6, 'Center', 500 * [cos((pi * i / 3) - pi/2), sin((pi*i/3)-pi/2)], 'SideLength', 500/sqrt(3));
    plot(pgon)
    txt = [num2str(i)];
    text(500 * cos((pi * i / 3) - pi/2) + 50, 500 * sin((pi*i/3)-pi/2), txt)
end

for i = 1:6
    BS(:, 6 + 2*i) = 500 * sqrt(3) * [cos(pi*i/3 - pi/3), sin(pi*i/3 - pi/3)];
    pgon = nsidedpoly(6, 'Center', 500 * sqrt(3) * [cos(pi*i/3 - pi/3), sin(pi*i/3 - pi/3)], 'SideLength', 500/sqrt(3));
    plot(pgon)
    txt = [num2str(6 + 2*i)];
    text(500 * sqrt(3) * cos(pi*i/3 - pi/3) + 50, 500 * sqrt(3) * sin(pi*i/3 - pi/3), txt)
end

for i = 1:6
    BS(:, 7 + 2*i) = 1000 * [cos((pi * i / 3) - pi/6), sin((pi*i/3)-pi/6)];
    pgon = nsidedpoly(6, 'Center', 1000 * [cos((pi * i / 3) - pi/6), sin((pi*i/3)-pi/6)], 'SideLength', 500/sqrt(3));
    plot(pgon)
    txt = [num2str(7 + 2*i)];
    text(1000 * cos((pi * i / 3) - pi/6) + 50, 1000* sin((pi*i/3)-pi/6), txt)
end

scatter(BS(1, :), BS(2, :), 'blue', 'filled')


SIMULATION_TIME = 900;


t = 0; % Time
dt = 0.1; % Time step
tstamp = 0;
dist = zeros(1, n_base_station);

for i = 1 : n_base_station
    dist(i) = sqrt((BS(1, i) - MS_pos_init(1))^2 + (BS(2, i) - MS_pos_init(2))^2);
end

min = dist(1);
min_temp = min;
id = 1;
id_temp = id;
n_handoff = 0;

for i = 1 : n_base_station
    if dist(i) < min
        min = dist(i);
        id = i;
    end
end

while t < SIMULATION_TIME   
    if t - tstamp <= MS_duration
        t = t + dt;
        MS_pos_init(1) = MS_pos_init(1) + MS_velocity * cos(MS_direction) * dt;
        MS_pos_init(2) = MS_pos_init(2) + MS_velocity * sin(MS_direction) * dt;

        for i = 1:n_base_station
            dist(i) = sqrt((BS(1, i) - MS_pos_init(1))^2 + (BS(2, i) - MS_pos_init(2))^2);
        end

        min_temp = dist(1);
        for i = 1 : n_base_station
            if dist(i) < min_temp
                min_temp = dist(i);
                id_temp = i;
            end
        end

        if min_temp > 500/sqrt(3)
            if(id_temp + 6 > 19)
                id_temp = id_temp - 6;
            else
                id_temp = id_temp + 6;
            end
            MS_pos_init(1) = -MS_pos_init(1);
            MS_pos_init(2) = -MS_pos_init(2);
        end

        if id_temp ~= id
            fprintf("Time: %.2f s, Source cell ID: %d, Destination cell ID: %d\n", t, id, id_temp);
            id = id_temp;
            n_handoff = n_handoff + 1;
        end
    else
        tstamp = t;
        scatter(MS_pos_init(1), MS_pos_init(2), 'red', 'filled')
        txt = num2str(t);
        text(MS_pos_init(1) + 5, MS_pos_init(2), txt)
        MS_direction = 2*pi*rand;
        MS_velocity = (maxSpeed - minSpeed)*rand + minSpeed;
        MS_duration = (maxT - minT)*rand + minT;
    end
end

fprintf("Number of handoffs: %d\n", n_handoff);

%% Task 2. 100 Moving Mobile Devices, Uplink
n_mobile_device = 100; % Number of mobile devices in the CBS
CBS = [0, 0]; % Location of Central Base Station
location_mobile_device_x = 500/sqrt(3) * 2 * (rand(n_mobile_device, 1) - 0.5);
location_mobile_device_y = 250 * 2 * (rand(n_mobile_device, 1) - 0.5);

% Scatter: BS and mobile devices
%{
figure;
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

%}
