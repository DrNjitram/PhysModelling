function exit_time = projecttube(h, src, vent, showgr)
    if nargin < 1
        h = 0.18;
        showgr = true;
        src = 0.05;
        vent = 0.200;
    end

    %% physical parameters
    length_hood = 0.9; % height of entire hood
    width = 0.695; % depth of hood

    size_vent = 0.3; % size of vent
    height_vent = 0.025; %height of vent
    vent_str = vent; % 0.200-0.317 m3/s

    height = h; % 0.1 - 0.65 height of sash 0.302 is crossover
    hood_thick = 0.025; % thickness of sash

    blockage_thick = 0.01; % thickness of the blockage
    blockage_offset = [0.638 0.250]; % x y offset of the blockage

    width_source = 0.00635; % width of source
    speed_s = src; % speed of oxygen (m3/s) min 0.0035
    source_pos = [0.556 0.08]; % x y pos of source

    % natural constants
    diffusion = 100 * 0.000176; % m2/s (oxygen in air)
    density = 1.2041; % air at room temperature


    %% Model Parameters
    max_mesh = 0.05; % maximum mesh size ( scaled by x/length)0.01 for good results, but slow, 0.05 for fast testing
    sim_length = 10; % duration (s)
    sim_steps = sim_length*100; % time steps
    trip_gas = 0.0002; % amount of gas (m3) at entrance before gas has escaped


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
    flowresults = airflow(height, speed_s, vent_str, blockage_offset, showgr); % retrieve flow model with hood at a certain height and show the results in a graph


    %% Boundary Conditions
    edge_cond = @(location, state) sum(evaluateGradient(flowresults, [location.x; location.y]));

    applyBoundaryCondition(model,'neumann','Edge', 1:20 ,'q', 0,'g', 0); % static walls
    applyBoundaryCondition(model,'neumann','Edge', [1 5] ,'q', edge_cond,'g', 0); % 'open' walls
    applyBoundaryCondition(model,'dirichlet','Edge', 7,'r', speed_s,'h',1); % source


    %% set coefficients in PDE
    f = @(location, state) sum(evaluateGradient(flowresults, [location.x; location.y])' .* [state.ux; state.uy] * -1); % implement spread of gas based on flow model

    specifyCoefficients(model, 'm', 0, 'd', density, 'c', diffusion, 'a', 0, 'f', f); % Coefficients for model

    setInitialConditions(model, 0); % Initial condition is no gas anywhere


    %% solve
    tlist = linspace(0, sim_length, sim_steps);
    results = solvepde(model, tlist);


    %% determine when front is reached
    y = linspace(0, height, 101);
    x = zeros(1, length(y));

    exit_time = 0;
    for i = 1:sim_steps
        entrance_values = interpolateSolution(results, x, y, i); % measure values at entrance
        volume = trapz(y, entrance_values); % integrate amount over entrance
        if volume > trip_gas % compare check because float maths and i dont trust matlab
            exit_time = i;
            break
        end
    end

    if exit_time ~= 0 && showgr % the gas escapes
        figure;
        plot(y, entrance_values);
        title(['Gas exited at t=', num2str(sim_length * exit_time/sim_steps), 's']);
    elseif showgr % the gas never escapes
        disp(['The gas never escapes within the given time range of ', num2str(sim_length), 's']);
        exit_time = sim_steps;
    end

    if showgr
        %% show results
        figure;
        subplot(1, 3, 1);
        t = 1;
        pdeplot(model,'XYData',results.NodalSolution(:, t), 'Contour', 'on','ColorMap','hot', 'Levels', 20);
        title(['Spread of gas at t=', num2str(sim_length * (t - 1)/sim_steps), 's']);
        xlabel("x (m)");
        ylabel('y (m)');
        c = colorbar;
        c.Label.String = 'gas concentration (???)';

        subplot(1, 3, 2);
        t = round(sim_steps/2);
        pdeplot(model,'XYData',results.NodalSolution(:, t), 'Contour', 'on','ColorMap','hot', 'Levels', 20);
        title(['Spread of gas at t=', num2str(sim_length * t/sim_steps), 's']);
        xlabel("x (m)");
        ylabel('y (m)');
        c = colorbar;
        c.Label.String = 'gas concentration (???)';

        subplot(1, 3, 3);
        t = sim_steps;
        pdeplot(model,'XYData',results.NodalSolution(:, t), 'Contour', 'on','ColorMap','hot', 'Levels', 20);
        title(['Spread of gas at t=', num2str(sim_length * t/sim_steps), 's']);
        xlabel("x (m)");
        ylabel('y (m)');
        c = colorbar;
        c.Label.String = 'gas concentration (???)';

        %% show results at the exit time
        figure;
        pdeplot(model,'XYData',results.NodalSolution(:, exit_time), 'Contour', 'on','ColorMap','hot', 'Levels', 20);
        title(['Spread of gas at t=', num2str(sim_length * exit_time/sim_steps), 's']);
        xlabel("x (m)");
        ylabel('y (m)');
        c = colorbar;
        c.Label.String = 'gas concentration (???)';
    end
end