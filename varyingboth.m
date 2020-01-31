close all;
clear;
clc;
format short;

source_speed = 0.05; % strength of source in hood, ~0.05 m3/s
show_graphs = false; %show intermediate graphs or not
threshold = 0.0002; % amount of gas (m3) at entrance before gas has escaped

speeds = 0:0.01:0.50; % range of vent speeds to check
heights = 0.1:0.01:0.65; % height of sash, 0.1-0.65m

% this will take >1h to run
volumes = zeros(length(speeds), length(heights));
for i = 1:length(speeds)
    for j = 1:length(heights)
        vent_speed = speeds(i);
        height =  heights(j);
        
        volumes(i, j) = stationary(height, source_speed, vent_speed, show_graphs);
        
        disp([num2str(vent_speed), ',', num2str(height), ',', num2str(volumes(i, j))]);
    end
end

figure
surf(heights, speeds, volumes)
xlabel('vent speed (m3/s)')
ylabel('vent speed (m3/s)')
xlabel('height of sash (m)')
zlabel('amount of gas escaped (m3)')