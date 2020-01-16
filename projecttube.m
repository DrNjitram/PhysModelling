close all;
clear;
clc;
format short;

%% physical parameters
length = 0.5; % length in m
ratio = 0.33; % ratio of height to length
source_radius = 0.01; % radius of source in m
position_source = [0.1 0.5]; % position of source (ratio of length and height)

diffusion = 2 * 0.000176; % m2/s (oxygen in air)
density = 1.2041; % air at room temperature
production_s = 1; % production of oxygen (L/s)

% combined value must not exceed 0.215
vx = 0.515; % air speed in x in m/s
vy = 0.0; % air speed in y in m/s

%% Model Parameters
max_mesh = 0.005; % maximum mesh size ( scaled by x/length)
sim_length = 1; % duration (s)
sim_steps = 100; % time steps

%% create a model container
model = createpde();

%% set geometry
height = ratio*length;
poly = [3; 4;... 
    0; length; length;  0;...
    0; 0;     height; height;];

source = [1; length*position_source(1); height*position_source(2); source_radius; 0; 0; 0; 0; 0; 0];

gd = [poly source];
ns = char('P1', 'C1')';
sf = 'P1 - C1';

g = decsg(gd, sf, ns);
geometryFromEdges(model, g);

%% Generate Mesh
generateMesh(model, 'Hmax', max_mesh/length);

%% Boundary Conditions
source_strength = production_s / (2*pi*source_radius); 
applyBoundaryCondition(model,'neumann','Edge', 1:4 ,'q',vx,'g',0); % static walls
applyBoundaryCondition(model,'dirichlet','Edge', 5:8, 'h', 1, 'r', source_strength/4); % source edge


%% set coefficients in PDE
% use location.x/y/x and state.u/ux/uy/uz
f =  @(location, state) -1*( vx*state.ux + vy*state.uy);

specifyCoefficients(model, 'Face', 1, 'm', 0, 'd', density, 'c', diffusion, 'a', 0, 'f', f); % outer pipe

setInitialConditions(model, 0, 'Face', 1);

%% solve
tlist = linspace(0, sim_length, sim_steps);
results = solvepde(model, tlist);

%% show results
figure;
subplot(3, 1, 1);
t = 1;
ux = results.XGradients(:, t);
uy = results.YGradients(:, t);
pdeplot(model,'XYData',results.NodalSolution(:, t), 'Contour', 'on','ColorMap','hot');

subplot(3, 1, 2);
t = round(sim_steps/2);
ux = results.XGradients(:, t);
uy = results.YGradients(:, t);
pdeplot(model,'XYData',results.NodalSolution(:, t), 'Contour', 'on','ColorMap','hot');

subplot(3, 1, 3);
t = sim_steps;
ux = results.XGradients(:, t);
uy = results.YGradients(:, t);
pdeplot(model,'XYData',results.NodalSolution(:, t), 'Contour', 'on','ColorMap','hot');