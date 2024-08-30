% Download OptimalTransportNetworkToolbox from https://github.com/SebKrantz/OptimalTransportNetworkToolbox
% and set the path to the lib folder here. Possibly also need to add a path to Ipopt if not detected. See docs/User Guide.pdf
addpath("code/OptimalTransportNetworkToolbox/lib")
% Read helper functions
addpath("code/matlab_helpers")

% Read Undirected Graph
graph_orig = readtable('data/transport_network/csv/graph_orig.csv');
graph_orig.add = false(height(graph_orig), 1);
graph_orig.Properties.VariableNames{'ug_cost_km'} = 'cost_per_km';

% Read Additional Links 
graph_add = readtable('data/transport_network/csv/graph_add.csv');
graph_add.add = true(height(graph_add), 1);
graph_add.Properties.VariableNames{'duration_100kmh'} = 'duration';
graph_add.Properties.VariableNames{'cost_km'} = 'cost_per_km';

% Combining
graph = [graph_orig(:, {'add', 'from', 'to', 'sp_distance', 'distance', 'duration', 'cost_per_km', 'border_dist'}); ...
         graph_add(:, {'add', 'from', 'to', 'sp_distance', 'distance', 'duration', 'cost_per_km', 'border_dist'})];

% Adjusting 
graph.distance = graph.distance / 1000;          % Convert to km
graph.sp_distance = graph.sp_distance / 1000;    % Convert to km
graph.border_dist = graph.border_dist / 1000;    % Convert to km
graph.duration = graph.duration / 60;            % Convert to hours
graph.cost_per_km = graph.cost_per_km / 1000000; % Convert to millions
graph.total_cost = graph.cost_per_km .* graph.distance

n = max([max(graph.from); max(graph.to)]);

% Create Adjacency Matrix
adj_matrix = false(n, n);
for i = 1:height(graph)
    adj_matrix(graph.from(i), graph.to(i)) = true;
    adj_matrix(graph.to(i), graph.from(i)) = true;
end

% Read Nodes Data
nodes = readtable('data/transport_network/csv/graph_nodes.csv');
nodes.population = nodes.population / 1000; % Convert to thousands
nodes.outflows = nodes.outflows / 1000;     % Adjust in line with population
% Describe function is not directly available in MATLAB. Instead, you can use summary.
summary(nodes);

% Create Infrastructure Matrix: Following Graff (2024)
infra_matrix = zeros(n, n);
for i = 1:height(graph)
    if graph.add(i)
        infra_matrix(graph.from(i), graph.to(i)) = 0.0;
        infra_matrix(graph.to(i), graph.from(i)) = 0.0;
    else
        speed = graph.distance(i) / graph.duration(i);
        infra_matrix(graph.from(i), graph.to(i)) = speed;
        infra_matrix(graph.to(i), graph.from(i)) = speed;
    end
end

% Create Iceberg Trade Cost Matrix: Following Graff (2024)
% graph.distance = graph.distance + graph.border_dist; % uncomment to enable frictions
iceberg_matrix = zeros(n, n);
for i = 1:height(graph)
    iceberg_matrix(graph.from(i), graph.to(i)) = 0.1158826 * log(graph.distance(i) / 1.609);
    iceberg_matrix(graph.to(i), graph.from(i)) = 0.1158826 * log(graph.distance(i) / 1.609);
end
iceberg_matrix(iceberg_matrix < 0) = 0;

% Create Infrastructure Building Cost Matrix: Following Graff (2024)
infra_building_matrix = zeros(n, n);
for i = 1:height(graph)
    infra_building_matrix(graph.from(i), graph.to(i)) = graph.total_cost(i);
    infra_building_matrix(graph.to(i), graph.from(i)) = graph.total_cost(i);
end


% Basic characteristics of the economy
population = nodes.population;
population = population + (population == 0) * 1e-6; % this is needed as ipopt crashes otherwise (need to have population > 0 always), so I add an infinitisimal person

% Productivity Matrix
with_ports = true;
productivity = zeros(n, 4); 
if with_ports
    for i = 1:n
        if nodes.outflows(i) > 0
            ind = 4;
            productivity(i, ind) = (37 * nodes.outflows(i)) / population(i); % Productivity of ports
        else
            if population(i) >= 1000
                ind = 3;
            elseif population(i) >= 200
                ind = 2;
            else
                ind = 1;
            end
        end
        productivity(i, ind) = productivity(i, ind) + nodes.IWI(i);
    end
else
    for i = 1:n
        if nodes.outflows(i) > 0
            ind = 4;
        elseif population(i) >= 1000
            ind = 3;
        elseif population(i) >= 200
            ind = 2;
        else
            ind = 1;
        end
        productivity(i, ind) = nodes.IWI(i);
    end
end

size(productivity)
[min(productivity(:)), max(productivity(:))]
for i = 1:size(productivity,2)
    [min(productivity(:,i)), max(productivity(:,i))]
end

[J, N] = size(productivity);

%% Infrastructure Bounds
min_mask = infra_matrix;                        
max_mask = max(infra_matrix, adj_matrix * 100);   

% Cost of all Road Activities (162 billion)
K_inv = sum(sum(infra_building_matrix)) / 2; 
% Total implied costs
infra_building_matrix = infra_building_matrix ./ (max_mask - infra_matrix);
infra_building_matrix(isinf(infra_building_matrix)) = 0;
infra_building_matrix(isnan(infra_building_matrix)) = 0;
K_base = sum(sum(infra_building_matrix .* infra_matrix)) / 2;
% Now setting budget
K = (K_base + 50e3) * 2 % 50 billion USD'15 investment volume (*2 because symmetric)

% Check Budget
K / sum(sum(infra_building_matrix .* max_mask))
K / sum(sum(infra_building_matrix .* min_mask))
if K > sum(sum(infra_building_matrix .* max_mask))
    error("Infrastructure budget is too large for the network")
end

% Parameters
alpha = 0.7; % Spending share on traded goods in utility (curvature parameter)
gamma = 0.946; % F&S: 0.1 ; TG: 0.946; Parameter governing intensity of congestion in transport
beta = 1.2446 * gamma; % F&S: 0.13 ; TG: 1.2446 * gamma; Parameter governing returns to scale in infrastructure investment
% gamma = beta^2/gamma; % IRS case
sigma = 1.5; %5; % Elasticity of substitution parameter
a = 1; % F&S: 1; TG: 0.7; Returns to scale to labor in production function Zn * Ln^a
rho = 0; % inequality aversion: not possible to solve model if enabled (= 2) -> set alpha = 0.1 instead and resolve allocation afterwards

% Initialise geography
param = init_parameters('Annealing', 'on', 'ADiGator', 'off', 'LaborMobility', 'off', 'CrossGoodCongestion', 'on', ...
                        'a', a, 'sigma', sigma, 'N', N, 'alpha', alpha, 'beta', beta, 'gamma', gamma, 'rho', rho, ...
                        'verbose', 'on', 'K', K, 'TolKappa', 1e-5);
[param, g] = create_graph(param,[],[],'X',nodes.lon,'Y',nodes.lat,'Type','custom','Adjacency',adj_matrix);

param.Lj = population; 
param.Zjn = productivity;
% param.omegaj = weights; % In case want to give more weight to specific locations
param.Hj = population .* (1-alpha); % TG: I normalise this because the general utility function has a (h_j/(1-alpha))^(1-alpha) thing with it

g.delta_i = infra_building_matrix;
g.delta_tau = iceberg_matrix;

% Naming conventions: if IRS -> 'cgc_irs'; if alpha = 0.1 -> add '_alpha01'; if with_ports = false -> add '_noport'; if frictions -> add '_bc' (border cost)
filename = '4g_50b_fixed_cgc_sigma15'; % adjust if sigma != 1.5
fprintf('File extension: %s\n', filename)

% Solve allocation from existing infrastructure
strcat("Started P_stat on ", datestr(datetime('now')))
res_stat = solve_allocation(param, g, infra_matrix, true);

% Solve Optimal Network
strcat("Started P_opt on ", datestr(datetime('now')))
res_opt = optimal_network(param, g, infra_matrix, min_mask, max_mask, false);
% res_ann = annealing(param,g,res.Ijk,'Il',min_mask,'Iu',max_mask); % Annealing is enabled by default if gamma > beta

% Plot Network
plot_graph(param, g, res_opt.Ijk, 'Margin', 0.01);
axis equal;
plot_graph(param, g, res_opt.Ijk - infra_matrix, 'Margin', 0.01);
axis equal;

% Saving: Nodes
res_nodes = nodes;
res_nodes.uj_orig = res_stat.uj;
res_nodes.Lj_orig = res_stat.Lj;
res_nodes.Cj_orig = res_stat.Cj;
res_nodes.Dj_orig = res_stat.Dj;
res_nodes.PCj_orig = res_stat.PCj;
% res_nodes.welfare = res_opt.welfare; % Same as: sum(res_opt.Lj .* res_opt.uj)
res_nodes.uj = res_opt.uj;
res_nodes.Lj = res_opt.Lj;
res_nodes.Cj = res_opt.Cj;
res_nodes.Dj = res_opt.Dj;
res_nodes.PCj = res_opt.PCj;
for n=1:N
   res_nodes = setfield(res_nodes, ['Lj_', num2str(n)], res_opt.Ljn(:,n));
   res_nodes = setfield(res_nodes, ['Dj_', num2str(n)], res_opt.Djn(:,n));
   res_nodes = setfield(res_nodes, ['Yj_', num2str(n)], res_opt.Yjn(:,n));
   res_nodes = setfield(res_nodes, ['Pj_', num2str(n)], res_opt.Pjn(:,n));
end
writetable(res_nodes, sprintf('results/transport_network/regional/nodes_results_%s.csv', filename))

% Saving: Graph / Edges
res_graph = graph;
res_graph.Ijk_orig = res_to_vec(infra_matrix, graph);
res_graph.Ijk = res_to_vec(res_opt.Ijk, graph);
for n=1:N
   res_graph = setfield(res_graph, ['Qjk_', num2str(n)], res_to_vec(res_opt.Qjkn(:,:,n), graph));
end
writetable(res_graph, sprintf('results/transport_network/regional/edges_results_%s.csv', filename));

