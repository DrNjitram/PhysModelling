close all;
clear;
clc;
format short;

source_speed = 0.05; % strength of source in hood, ~0.05 m3/s
vent_speed = 0.200; % strengh of vent in hood, 0.1572-0.498 m3/s
show_graphs = false; %show intermediate graphs or not
threshold = 0.0002; % amount of gas (m3) at entrance before gas has escaped

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