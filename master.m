close all;
clear;
clc;
format short;

source_speed = 0.05; % strength of source in hood, ~0.05 m3/s
vent_speed = 0.200; % strengh of vent in hood, 0.1572-0.498 m3/s
show_graphs = false; %show intermediate graphs or not
threshold = 0.0002; % amount of gas (m3) at entrance before gas has escaped

heights = 0.1:0.01:0.65; % height of sash, 0.1-0.65m
speeds_hi = 0.20:0.01:0.50; % range of vent speeds to check
speeds_lo = 0:0.01:0.20; % range of vent speeds to check


volumes = zeros(1, length(heights));
times = zeros(1, length(heights));
for i = 1:length(heights)
    height = heights(i);
    volumes(i) = stationary(height, source_speed, vent_speed, show_graphs);
    if volumes(i) < threshold
        times(i) = -1;
    else
        times(i) = projecttube(height, source_speed, vent_speed, show_graphs);
    end
    disp([num2str(height), ',', num2str(volumes(i)), ',', num2str(times(i))]);
end
figure
plot(heights, volumes);
figure
plot(heights, times);

volumes = zeros(1, length(speeds_lo));
height = 0.1;

for i = 1:length(speeds_lo)
    vent_speed = speeds_lo(i);
    volumes(i) = stationary(height, source_speed, vent_speed, show_graphs);
    disp([num2str(vent_speed), ',', num2str(volumes(i))]);
end
figure
plot(speeds_lo, volumes);

volumes = zeros(1, length(speeds_hi));
height = 0.65;

for i = 1:length(speeds_hi)
    vent_speed = speeds_hi(i);
    volumes(i) = stationary(height, source_speed, vent_speed, show_graphs);
    disp([num2str(vent_speed), ',', num2str(volumes(i))]);
end
figure
plot(speeds_hi, volumes);


