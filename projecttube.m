close all;
clear;
clc;
format short;

%% physical parameters
length_hood = 0.9; % height of entire hood
size_opening = 0.3; % size of vent
height = 0.1; % 0.1 - 0.65 height of opening
width = 0.695; % depth of hood
hood_thick = 0.025; % thickness of hood

blockage_thick = 0.01; % thickness of the blockage
blockage_offset = [0.638 0.250]; % x y offset of the blockage


width_source = 0.01;
diffusion = 10 * 0.000176; % m2/s (oxygen in air)
density = 1.2041; % air at room temperature
production_s = 0.1; % production of oxygen (L/s)

% combined value must not exceed 0.215
vx = 0.0; % air speed in x in m/s
vy = 0.0; % air speed in y in m/s

%% Model Parameters
max_mesh = 0.01; % maximum mesh size ( scaled by x/length)
sim_length = 50; % duration (s)
sim_steps = 100; % time steps

%% create a model container
model = createpde();

%% set geometry
fume = [2; 10;... % Fume hood
    0; width; width;       width - 0.057; width - 0.057;      width - size_opening - 0.057;width - size_opening - 0.057; hood_thick;  hood_thick;  0;...
    0; 0;     length_hood; length_hood;   length_hood + 0.05; length_hood+0.05;length_hood;                  length_hood;      height;       height];
        
source = [3; 4;... % source
    width * 0.8; width * 0.8 + width_source; width * 0.8 + width_source; width * 0.8;...
    0;           0;                  0.08;               0.08;];

blocking = [2; 6;... % blockage
    blockage_offset(1); blockage_offset(1) + blockage_thick; blockage_offset(1) + blockage_thick - 0.143; blockage_offset(1) + blockage_thick * 2.5 - 0.343; blockage_offset(1) - 0.321; blockage_offset(1) - 0.143;...
    blockage_offset(2); blockage_offset(2)                 ; blockage_offset(2) + 0.558                 ; blockage_offset(2) + 0.605                       ; blockage_offset(2) + 0.595; blockage_offset(2) + 0.55;];


source = [source; zeros(length(fume)-length(source),1)];
blocking = [blocking; zeros(length(fume)-length(blocking),1)];


gd = [fume source blocking];
ns = char('P1', 'SQ1', 'P2')';
sf = 'P1 - SQ1 - P2';


g = decsg(gd, sf, ns);
geometryFromEdges(model, g);

%% Generate Mesh
generateMesh(model, 'Hmax', max_mesh);
%pdegplot(model, 'EdgeLabels', 'on');

%% Boundary Conditions
source_strength = production_s / width_source; 
applyBoundaryCondition(model,'neumann','Edge', 1:20 ,'q',0,'g',0); % static walls
applyBoundaryCondition(model,'neumann','Edge', [1 5] ,'q',vx,'g',0); % 'open' walls
applyBoundaryCondition(model,'dirichlet','Edge', 7, 'h', 1, 'r', source_strength); % source edge


%% set coefficients in PDE
% use location.x/y/x and state.u/ux/uy/uz
f =  @(location, state) -1*( vx*state.ux + vy*state.uy);

specifyCoefficients(model, 'm', 0, 'd', density, 'c', diffusion, 'a', 0, 'f', f); % outer pipe

setInitialConditions(model, 0);

%% solve
tlist = linspace(0, sim_length, sim_steps);
results = solvepde(model, tlist);

%% show results
figure;
% subplot(3, 1, 1);
% t = 1;
% ux = results.XGradients(:, t);
% uy = results.YGradients(:, t);
% pdeplot(model,'XYData',results.NodalSolution(:, t), 'Contour', 'on','ColorMap','hot');
% title(['Spread of gas at t=', num2str(sim_length * (t - 1)/sim_steps), 's']);
% xlabel("x (m)");
% ylabel('y (m)');
% c = colorbar;
% c.Label.String = 'gas concentration (???)';
% 
% subplot(3, 1, 2);
% t = round(sim_steps/2);
% ux = results.XGradients(:, t);
% uy = results.YGradients(:, t);
% pdeplot(model,'XYData',results.NodalSolution(:, t), 'Contour', 'on','ColorMap','hot');
% title(['Spread of gas at t=', num2str(sim_length * t/sim_steps), 's']);
% xlabel("x (m)");
% ylabel('y (m)');
% c = colorbar;
% c.Label.String = 'gas concentration (???)';
% 
% subplot(3, 1, 3);
t = sim_steps;
ux = results.XGradients(:, t);
uy = results.YGradients(:, t);
pdeplot(model,'XYData',results.NodalSolution(:, t), 'Contour', 'on','ColorMap','hot');
title(['Spread of gas at t=', num2str(sim_length * t/sim_steps), 's']);
xlabel("x (m)");
ylabel('y (m)');
c = colorbar;
c.Label.String = 'gas concentration (???)';