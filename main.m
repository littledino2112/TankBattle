[G, node_to_coordinate, gas_list, heart_list, nrow, ncol] = buildgraph();
root_node_coordinate = [1 1];
root_node_node_id = coordinate_to_node_idx(root_node_coordinate,nrow);
root_node_name = num2str(root_node_node_id);
global gas_tank_max;
gas_tank_max = 25;
gas_tank = gas_tank_max;


% Build 2 cell arrays containing node names for heart and gas station
gas_list_name = cell(length(gas_list),1);
heart_list_name = cell(length(heart_list),1);
for i=1:length(gas_list)
    gas_list_name(i) = {['G_' num2str(gas_list(i))]};
end
for i=1:length(heart_list)
    heart_list_name(i) = {['H_' num2str(heart_list(i))]};
end
% Build a new graph T with all gas stations
% Create edges to connect stations that are at most gas_tank apart - This
% is already done above by setting all element of the adjecancy matrix
% larger than gas_tank to 0
% Compute distances btw gas stations
distances_gas = distances(G,gas_list,gas_list);
[row, col] = size(distances_gas);
for i=1:row
    for j=i+1:col
        if distances_gas(i,j) > gas_tank
            distances_gas(i,j) = 0;
            distances_gas(j,i) = 0;
        end
    end
end
% T = graph(distances_gas,node_names(gas_list));
T = graph(distances_gas,gas_list_name);
% Variables to keep track if connection from start node to gas stations and
% hearts have been done
root_node_to_gas_stations = 0; 
root_node_to_hearts = 0;

% Create edges from heart position to all gas stations within (gas_tank-dx)
% Add root node (starting position) and create connection
for i=1:length(heart_list)
    dis_to_gas_stations = distances(G,heart_list(i),gas_list);
    [dx, closest_gas_station] = min(dis_to_gas_stations);
    radius = gas_tank - dx;
    for j=1:length(gas_list)
        if dis_to_gas_stations(j) < radius
            T = addedge(T,heart_list_name(i),gas_list_name(j),dis_to_gas_stations(j));
        end
        if ~root_node_to_gas_stations
            distance_root_to_gas = distances(G,gas_list(j),root_node_node_id);
            if  distance_root_to_gas <= gas_tank
                T = addedge(T,root_node_name,gas_list_name(j),distance_root_to_gas);
            end
        end
    end
    root_node_to_gas_stations = 1;
end

% Create edges from heart to heart if distance btw them (in original graph
% G) is at most (gas_tank - dx - dy) with dx,dy is the distance to the
% closest gas station of x and y, respectively
for i=1:length(heart_list)-1
    if ~findnode(T,heart_list_name(i))
           continue; 
    end
    dx = [];
    dy = [];
    for j=i+1:length(heart_list)
        if ~findnode(T,heart_list_name(j))
           continue; 
        end
        dx = distances(T,heart_list_name(i),gas_list_name);
        dx = min(dx);
        dy = distances(T,heart_list_name(j),gas_list_name);
        dy = min(dy);
        dxy = distances(G,heart_list(i),heart_list(j));
        if dxy <= (gas_tank - dx(1) - dy(1))
            T = addedge(T,heart_list_name(i),heart_list_name(j),dxy);
        end
    end
    if ~root_node_to_hearts
        dis_root_to_heart = distances(G,heart_list(i),root_node_node_id);
        if dis_root_to_heart < (gas_tank - dx(1))
            T = addedge(T,heart_list_name(i),root_node_name,dis_root_to_heart);
        end
    end
end
root_node_to_hearts = 1;

% Find Heart position names in graph T since not all hearts in G end up in
% T
% Find the indexes of Hearts in node array
node_names_in_T = T.Nodes.Name;
heart_list_in_T = [];
for idx=1:length(node_names_in_T)
    if node_names_in_T{idx}(1) == 'H'
        heart_list_in_T = [heart_list_in_T idx];  %#ok<AGROW>
    end
end

% Create a new graph to contain all the hearts and connection btw them
nodes_in_final_graph = [node_names_in_T(heart_list_in_T); root_node_name];
heart_distances = distances(T,nodes_in_final_graph,nodes_in_final_graph);
H = graph(heart_distances,nodes_in_final_graph);
global travel_map;
travel_map = {};
mintree = minspantree(H);
directed_mintree = shortestpathtree(mintree,root_node_name);
traverse_mst(directed_mintree,root_node_name);

%% Go back from map H to map T
travel_map_T = {};
for i=1:length(travel_map)-1
    p = shortestpath(T,travel_map{i},travel_map{i+1});
    travel_map_T = [travel_map_T p(1:end-1)]; %#ok<AGROW>
end
travel_map_T(1) = {root_node_node_id};
% for i=2:length(travel_map_T)
%     travel_map_T(i) = {str2double(travel_map_T{i}(3:end))}; %#ok<SAGROW>
% end

%% Traverse the final map and print out coordinates and gas tank
step_and_fuel = [];
for i=1:length(travel_map_T)-1
    % If current location is a heart, visit its nearset gas station to
    % refill before moving on
    destination = 0;
    source = 0;
    if travel_map_T{i}(1) == 'H'
        d = distances(T,travel_map_T(i),gas_list_name);
        [min_val, nearest_station_idx] = min(d);
        nearest_station = gas_list_name{nearest_station_idx};
        [path,fuel_left] = find_steps_and_fuel(G,travel_map_T{i},nearest_station,gas_tank);
        step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]]; %#ok<*AGROW>
        gas_tank = fuel_left(end);
        [path,fuel_left] = find_steps_and_fuel(G,nearest_station,travel_map_T{i},gas_tank);
        step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]];
        gas_tank = fuel_left(end);
    end
    [path,fuel_left] = find_steps_and_fuel(G,travel_map_T{i},travel_map_T{i+1},gas_tank);
    gas_tank = fuel_left(end);
    step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]];
end


%%
function traverse_mst(graph, node)
    global travel_map;
    if isempty(find(strcmp(travel_map,node), 1))
        travel_map = [travel_map node]; 
    end
    child_list = successors(graph,node);
    for i=1:length(child_list)
        traverse_mst(graph,child_list{i});
    end
end

%% Function to find what step to take and fuel left from node A to B
function [steps, fuel] = find_steps_and_fuel(graph,node_1,node_2,current_gas_tank)
    global gas_tank_max;
    node_id_1 = node_1;
    node_id_2 = node_2;
    if ischar(node_1)
        node_id_1 = str2double(node_1(3:end));
    end
    if ischar(node_2)
        node_id_2 = str2double(node_2(3:end));
    end
    [steps,d] = shortestpath(graph,node_id_1,node_id_2);
    steps = steps.';
    fuel = current_gas_tank:-1:current_gas_tank-d;
    fuel = fuel.';
    if node_2(1) == 'G'
        current_gas_tank = gas_tank_max;
        fuel(end) = current_gas_tank;
    end
end

function node_idx = coordinate_to_node_idx(coordinate,nrow)
    node_idx = nrow*(coordinate(1)-1) + coordinate(2); 
end
