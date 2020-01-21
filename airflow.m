close all;
clear;
clc;
format short;

% physical parameters
% source http://www.conditionaire.com.au/school-type.html
dynamic_visc = 1.825E-5; % dynamic viscocity of air at 20C
length_hood = 0.9; % height of entire hood
size_opening = 0.3; % size of vent
height = 0.50; % 0.1 - 0.65 height of opening
width = 0.695; % depth of hood
hood_thick = 0.025; % thickness of hood
blockage_thick = 0.01; % thickness of the blockage
blockage_offset = [0.638 0.250]; % x y offset of the blockage

vent_strength = 200; % 200-317 L/s
source_strength = 0.025;

% Model Parameters
max_mesh = 0.01;

% create a model container
model = createpde();

% set geometry (see documentation)
fume = [2; 10;... % Fume hood
    0; width; width;       width - 0.057; width - 0.057;      width - size_opening - 0.057;width - size_opening - 0.057; hood_thick;  hood_thick;  0;...
    0; 0;     length_hood; length_hood;   length_hood + 0.05; length_hood+0.05;length_hood;                  length_hood;      height;       height];
        
source = [3; 4;... % source
    width * 0.8; width * 0.8 + 0.01; width * 0.8 + 0.01; width * 0.8;...
    0;           0;                  0.08;               0.08;];

blocking = [2; 6;... % blockage
    blockage_offset(1); blockage_offset(1) + blockage_thick; blockage_offset(1) + blockage_thick - 0.143; blockage_offset(1) + blockage_thick * 2.5 - 0.343; blockage_offset(1) - 0.321; blockage_offset(1) - 0.143;...
    blockage_offset(2); blockage_offset(2)                 ; blockage_offset(2) + 0.558                 ; blockage_offset(2) + 0.605                       ; blockage_offset(2) + 0.595; blockage_offset(2) + 0.55;];

outside = [3; 4;...
    -1; 0; 0; -1;...
    -0.5; -0.5; 1.2; 1.2;];

source = [source; zeros(length(fume)-length(source),1)];
blocking = [blocking; zeros(length(fume)-length(blocking),1)];
outside = [outside; zeros(length(fume)-length(outside),1)];

gd = [fume source blocking outside];
ns = char('P1', 'SQ1', 'P2', 'SQ2')';
sf = 'P1 - SQ1 - P2 + SQ2';


g = decsg(gd, sf, ns);
geometryFromEdges(model, g);


% Generate Mesh
generateMesh(model, 'Hmax', max_mesh);

%pdeplot(model)
pdegplot(model, 'EdgeLabels', 'on');

% Boundary Conditions
applyBoundaryCondition(model,'neumann','Edge', 1:25,'q',0,'g',0); % static walls
applyBoundaryCondition(model,'dirichlet','Edge', 1, 'h', 1, 'r', vent_strength); % vent
applyBoundaryCondition(model,'dirichlet','Edge', 6, 'h', 1, 'r', source_strength); % source
applyBoundaryCondition(model,'dirichlet','Edge', [15 14], 'h', 1, 'r', 1); % opening

% Specify model parameters
specifyCoefficients(model,'m',0,'d', 0,'c', -dynamic_visc ,'a',0,'f',0);
setInitialConditions(model, 1);

% solve
results = solvepde(model);

% show results
figure;
ux = results.XGradients;
uy = results.YGradients;
pdeplot(model,'XYData',results.NodalSolution, 'FlowData', [ux uy], 'Contour', 'on');
xlabel("x (m)");
ylabel('y (m)');
c = colorbar;
c.Label.String = 'air concentration (???)';



% Use [vx, vy] =  evaluateGradient(results,X,Y) to get gradient values
% is decently fast, show using f = @() evaluateGradient(results,0.5, 0.5);
% timeit(f)