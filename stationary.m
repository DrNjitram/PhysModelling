close all;
clear;
clc;
format short;


%% physical parameters
length_hood = 0.9; % height of entire hood
width = 0.695; % depth of hood

size_vent = 0.3; % size of vent
height_vent = 0.025; %height of vent

height = 0.65; % 0.1 - 0.65 height of sash 0.302 is crossover
hood_thick = 0.025; % thickness of sash

blockage_thick = 0.01; % thickness of the blockage
blockage_offset = [0.638 0.250]; % x y offset of the blockage


width_source = 0.00635; % width of source
speed_s = 5; % speed of oxygen (L/s) min 3.5
source_pos = [0.556 0.08]; % x y pos of source

% natural constants
diffusion = 10 * 0.000176; % m2/s (oxygen in air)


%% Model Parameters
max_mesh = 0.05; % maximum mesh size ( scaled by x/length) 0.05 works, below crashes?
trip_gas = 0.001; % amount of gas (m3) at entrance before gas has escaped


%% create a model container
model = createpde(1);

%% set geometry
fume = [2; 10;... % Fume hood
    0; width; width;       width - 0.057; width - 0.057;      width - size_vent - 0.057;width - size_vent - 0.057; hood_thick;  hood_thick;  0;...
    0; 0;     length_hood; length_hood;   length_hood + height_vent; length_hood+height_vent;length_hood;                  length_hood;      height;       height];
        
source = [3; 4;... % source
    source_pos(1) - width_source/2; source_pos(1) + width_source/2; source_pos(1) + width_source/2; source_pos(1) - width_source/2;...
    0;                              0;                              source_pos(2);                  source_pos(2);];

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


%% Retrieve flow model
flowresults = airflow(height, speed_s, blockage_offset, false); % retrieve flow model with hood at a certain height and show the results in a graph


%% Boundary Conditions
source_strength = (speed_s/1000)/(pi * (width_source/2)^2);
edge_cond = @(location, state) sum(evaluateGradient(flowresults, [location.x; location.y]));

applyBoundaryCondition(model,'neumann','Edge', 1:20 ,'q', 0,'g', 0); % static walls
applyBoundaryCondition(model,'neumann','Edge', [1 5] ,'q', edge_cond,'g', 0); % 'open' walls
applyBoundaryCondition(model,'neumann','Edge', 7,'g', source_strength,'q',1); % source


%% set coefficients in PDE
f = @(location, state) sum(evaluateGradient(flowresults, [location.x; location.y])' .* [state.ux; state.uy] * -1); % implement spread of gas based on flow model

specifyCoefficients(model, 'm', 0, 'd', 0, 'c', diffusion, 'a', 0, 'f', f); % Coefficients for model


%% solve
results = solvepde(model);

%% determine when front is reached
y = linspace(0, height, 101);
x = zeros(1, length(y));

entrance_values = interpolateSolution(results, x, y);
volume = trapz(y, entrance_values); % integrate amount over entrance

%% show results at the exit time
figure;
pdeplot(model,'XYData',results.NodalSolution, 'Contour', 'on','ColorMap','hot', 'Levels', 20);
if volume > trip_gas
    title(['Gas escaped:', num2str(volume)]);
else
    title('Gas never escaped');
end