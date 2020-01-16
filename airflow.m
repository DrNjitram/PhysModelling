close all;
clear;
clc;
format short;

% physical parameters
% source http://www.conditionaire.com.au/school-type.html
diffusion = 1;
length = 0.9; % height of entire hood
size_opening = 0.3; % size of vent
height = 0.1; % 0.1 - 0.65 height of opening
width = 0.695; % depth of hood
thick = 0.055; % thickness of hood

vent_strength = 200; % 200-317 L/s
source_strength = 0.025;

% Model Parameters
max_mesh = 0.01;

% create a model container
model = createpde();

% set geometry (see documentation)
fume = [2; 7;... % Fume hood
    0; width; width;  width - size_opening; thick;  thick;  0;...
    0; 0;     length; length;               length; height; height];

source = [3; 4;... % source
    width * 0.8; width * 0.8 + 0.01; width * 0.8 + 0.01; width * 0.8;...
    0;           0;                  0.08;               0.08;...
    0; 0; 0; 0; 0; 0]; % padding

blocking = [2; 6;...
    ];

gd = [fume source];
ns = char('P1', 'SQ1')';
sf = 'P1 - SQ1';

g = decsg(gd, sf, ns);
geometryFromEdges(model, g);


% Generate Mesh
generateMesh(model, 'Hmax', max_mesh);

%pdegplot(model, 'Edgelabels', 'on')
%pdeplot(model)

% Boundary Conditions
applyBoundaryCondition(model,'neumann','Edge', 1:10,'q',0,'g',0); % static walls
applyBoundaryCondition(model,'dirichlet','Edge', 9, 'h', 1, 'r', vent_strength); % vent
applyBoundaryCondition(model,'dirichlet','Edge', 7, 'h', 1, 'r', source_strength); % source
applyBoundaryCondition(model,'dirichlet','Edge', 4, 'h', 1, 'r', 1); % opening

% Specify model parameters
specifyCoefficients(model,'m',0,'d', 0,'c', 0.1 ,'a',0,'f',0);
setInitialConditions(model, 1);

% solve
results = solvepde(model);

% show results
figure;
ux = results.XGradients;
uy = results.YGradients;
pdeplot(model,'XYData',results.NodalSolution, 'FlowData', [ux uy], 'Contour', 'on');