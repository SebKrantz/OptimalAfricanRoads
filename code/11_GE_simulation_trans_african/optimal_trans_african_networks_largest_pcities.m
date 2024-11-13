% Download OptimalTransportNetworkToolbox from https://github.com/SebKrantz/OptimalTransportNetworkToolbox
% and set the path to the lib folder here. Possibly also need to add a path to Ipopt if not detected. See docs/User Guide.pdf
addpath('code/OptimalTransportNetworkToolbox/lib');
% Read helper functions
addpath('code/helpers');

% NOTE: can use fastest_routes or all_routes (incl. shortest routes -> larger graph)
% With all_routes, need to set beta = 1 and CrossGoodCongestion 'off' to use efficient dual solution

% Read Undirected Graph
edges = readtable("data/transport_network/trans_african/fastest_routes_graph_edges.csv");
% histogram(edges.ug_cost, 'BinWidth', 100);
% Adjusting 
edges.distance = edges.distance / 1000;       % Convert to km
edges.border_dist = edges.border_dist / 1000; % Convert to km
edges.duration = edges.duration / 60;         % Convert to hours
edges.total_cost = edges.ug_cost / 1e6;       % Convert to millions

n = max(max(edges.from), max(edges.to));

% Create Adjacency Matrix
adj_matrix = false(n, n);
for i = 1:height(edges)
    adj_matrix(edges.from(i), edges.to(i)) = true;
    adj_matrix(edges.to(i), edges.from(i)) = true;
end

% Read Nodes Data
nodes = readtable("data/transport_network/trans_african/fastest_routes_graph_nodes.csv");
nodes.population = nodes.population / 1000; % Convert to thousands
nodes.outflows = nodes.outflows / 1000;     % Adjust in line with population
summary(nodes);

% Create Infrastructure Matrix: Following Graff (2024) = average speed in km/h: length of route is accounted for in cost function
infra_matrix = zeros(n, n);
for i = 1:height(edges)
    speed = edges.distance(i) / edges.duration(i);
    infra_matrix(edges.from(i), edges.to(i)) = speed;
    infra_matrix(edges.to(i), edges.from(i)) = speed;
end
summary(edges.distance ./ edges.duration);

% Create Iceberg Trade Cost Matrix: Following Graff (2024)
% edges.distance = edges.distance + edges.border_dist; % uncomment to enable frictions
iceberg_matrix = zeros(n, n);
for i = 1:height(edges)
    iceberg_matrix(edges.from(i), edges.to(i)) = 0.1158826 * log(edges.distance(i) / 1.609);
    iceberg_matrix(edges.to(i), edges.from(i)) = 0.1158826 * log(edges.distance(i) / 1.609);
end
iceberg_matrix(iceberg_matrix < 0) = 0;

% Create Infrastructure Building Cost Matrix: Following Graff (2024)
infra_building_matrix = zeros(n, n);
for i = 1:height(edges)
    infra_building_matrix(edges.from(i), edges.to(i)) = edges.total_cost(i);
    infra_building_matrix(edges.to(i), edges.from(i)) = edges.total_cost(i);
end

% Basic characteristics of the economy
population = nodes.population;
population = population + (population == 0) * 1e-6;

% Productivity Matrix 
% Check largest pcities
assert(nnz(population > 2000 | nodes.outflows > 1000) == 47);
productivity = zeros(n, max(nodes.product));
with_ports = true;
for i = 1:n
    productivity(i, nodes.product(i)) = nodes.IWI(i);
    if with_ports && nodes.outflows(i) > 0
        productivity(i, nodes.product(i)) = productivity(i, nodes.product(i)) + (37 * nodes.outflows(i)) / population(i); % Productivity of ports
    end
end
disp(minmax(productivity));
assert(all(sum(productivity > 0, 2) == 1)); % Check for Armington parameterization

[J, N] = size(productivity);

% Infrastructure Bounds
min_mask = infra_matrix;                          
max_mask = max(infra_matrix, adj_matrix * 100);

% Cost of all Road Activities (33 billion)
K_inv = sum(sum(infra_building_matrix)) / 2;
% Total implied costs
infra_building_matrix = infra_building_matrix ./ (max_mask - infra_matrix);
infra_building_matrix(isinf(infra_building_matrix)) = 0;
infra_building_matrix(isnan(infra_building_matrix)) = 0;
K_base = sum(sum(infra_building_matrix .* infra_matrix)) / 2;
% Now setting budget
K = (K_base + 10e3) * 2; % 10 or 20 billion USD'15 investment volume (*2 because symmetric)

disp(K / sum(sum(infra_building_matrix .* max_mask)));
disp(K / sum(sum(infra_building_matrix .* min_mask)));
if K > sum(sum(infra_building_matrix .* max_mask))
    error('Infrastructure budget is too large for the network');
end

% Parameters
alpha = 0.7; % Spending share on traded goods in utility (curvature parameter)
gamma = 0.946; % F&S: 0.1 ; TG: 0.946; Parameter governing intensity of congestion in transport
beta = 1.2446 * gamma; % F&S: 0.13 ; TG: 1.2446 * gamma; Parameter governing returns to scale in infrastructure investment
% gamma = beta^2/gamma; % IRS case
sigma = 3.8; %5; % Elasticity of substitution parameter (3.8 = Armington)
a = 1; % F&S: 1; TG: 0.7; Returns to scale to labor in production function Zn * Ln^a
rho = 0; % 2 % inequality aversion

% Initialise geography
param = init_parameters('Annealing', 'on', 'ADiGator', 'off', 'LaborMobility', 'off', 'CrossGoodCongestion', 'on', ...
                        'a', a, 'sigma', sigma, 'N', N, 'alpha', alpha, 'beta', beta, 'gamma', gamma, 'rho', rho, ...
                        'verbose', 'on', 'K', K, 'TolKappa', 1e-5);

[param, graph] = create_graph(param,[],[],'X',nodes.lon,'Y',nodes.lat,'Type','custom','Adjacency',adj_matrix); 

param.Lj = population; 
param.Zjn = productivity;
% param.omegaj = weights; % In case want to give more weight to specific locations
param.Hj = population .* (1-alpha); % TG: I normalise this because the general utility function has a (h_j/(1-alpha))^(1-alpha) thing with it
                          
graph.delta_i = infra_building_matrix;
graph.delta_tau = iceberg_matrix;

% Naming conventions: if IRS -> 'cgc_irs' or 'irs_na' without annealing; if alpha = 0.1 -> add '_alpha01'; if with_ports = false -> add '_noport'; if frictions -> add '_bc' (border cost)
filename = '22g_10b_fixed_cgc_sigma38_rho0'; % adjust if sigma != 1.5
fprintf('File extension: %s\n', filename);

% Solve allocation from existing infrastructure
strcat("Started P_stat on ", datestr(datetime('now')))
res_stat = solve_allocation(param, graph, infra_matrix, true);

% Solve Optimal Network (this can take long - up to 48h)
strcat("Started P_opt on ", datestr(datetime('now')))
res_opt = optimal_network(param, graph, infra_matrix, min_mask, max_mask, false);

% Check: should be 1
K_ratio = K / sum(sum(res_opt.Ijk .* graph.delta_i));
disp(K_ratio);

% Plot Network
plot_graph(param, graph, res_opt.Ijk, 'Margin', 0.01);
axis equal;
plot_graph(param, graph, res_opt.Ijk - infra_matrix, 'Margin', 0.01);
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
for n = 1:N
   res_nodes = setfield(res_nodes, ['Lj_', num2str(n)], res_opt.Ljn(:,n));
   res_nodes = setfield(res_nodes, ['Dj_', num2str(n)], res_opt.Djn(:,n));
   res_nodes = setfield(res_nodes, ['Yj_', num2str(n)], res_opt.Yjn(:,n));
   res_nodes = setfield(res_nodes, ['Pj_', num2str(n)], res_opt.Pjn(:,n));
end
writetable(res_nodes, sprintf('results/transport_network/GE/trans_african/nodes_results_%s.csv', filename));

% Saving: Graph / Edges
res_edges = edges;
res_edges.Ijk_orig = res_to_vec(infra_matrix, edges);
res_edges.Ijk = res_to_vec(res_opt.Ijk, edges);
for n = 1:N
    res_edges = setfield(res_edges, ['Qjk_', num2str(n)], res_to_vec(res_opt.Qjkn(:,:,n), edges));
end
writetable(res_edges, sprintf('results/transport_network/GE/trans_african/edges_results_%s.csv', filename));
