function results = airflow(height, source_str, blockage, showgraph)
    % physical parameters
    dynamic_visc = 1.825E-5; % dynamic viscocity of air at 20C
     
    length_hood = 0.9; % height of entire hood
    width = 0.695; % depth of hood
    size_opening = 0.3; % size of vent

    hood_thick = 0.025; % thickness of sash

    blockage_thick = 0.01; % thickness of the blockage
    blockage_offset = blockage; % x y offset of the blockage


    width_source = 0.01; % width of source
    source_strength =(source_str/1000)/(pi * (width_source/2)^2); % strength, roughly 0.025 for testing
    source_pos = [0.556 0.08]; % x y pos of source
    
    
    % Model Parameters
    max_mesh = 0.01;
    vent_strength = 200; % 200-317 L/s
    vent_strength = (vent_strength/1000)/(pi * (size_opening/2)^2); 

    % create a model container
    model = createpde(1);

    % set geometry source http://www.conditionaire.com.au/school-type.html
    fume = [2; 10;... % Fume hood
        0; width; width;       width - 0.057; width - 0.057;      width - size_opening - 0.057;width - size_opening - 0.057; hood_thick;  hood_thick;  0;...
        0; 0;     length_hood; length_hood;   length_hood + 0.05; length_hood+0.05;length_hood;                  length_hood;      height;       height];

    source = [3; 4;... % source
        source_pos(1) - width_source/2; source_pos(1) + width_source/2; source_pos(1) + width_source/2; source_pos(1) - width_source/2;...
        0;                              0;                              source_pos(2);                  source_pos(2);];

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
    %pdegplot(model, 'EdgeLabels', 'on');

    % Boundary Conditions
    applyBoundaryCondition(model,'neumann','Edge', 1:25,'q',0,'g',0); % static walls
    applyBoundaryCondition(model,'neumann','Edge', 1,'q',vent_strength,'g',1); % vent
    applyBoundaryCondition(model,'neumann','Edge', 6,'q',source_strength,'g',1); % source
    applyBoundaryCondition(model,'neumann','Edge', [15 14],'q',vent_strength - source_strength,'g',1); % opening

    % Specify model parameters
    specifyCoefficients(model,'m',0,'d', 0,'c', -dynamic_visc ,'a',0,'f',0);
    setInitialConditions(model, 1);

    % solve
    results = solvepde(model);
    
    if showgraph
        % show results
        figure;
        ux = results.XGradients;
        uy = results.YGradients;
        pdeplot(model,'XYData',results.NodalSolution, 'FlowData', [ux uy], 'Contour', 'on');
        xlabel("x (m)");
        ylabel('y (m)');
        c = colorbar;
        c.Label.String = 'air pressure (P)';
    end
end


% Use [vx, vy] =  evaluateGradient(results,X,Y) to get gradient values
% is decently fast, show using f = @() evaluateGradient(results,0.5, 0.5);
% timeit(f)