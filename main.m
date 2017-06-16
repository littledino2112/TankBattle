global gas_tank_max;
input_file = 'map50.data';
root_node_coordinate = dlmread(input_file,' ',[0 0 0 1]);
gas_tank_max = dlmread(input_file,' ',[1 0 1 0]);
gas_tank = gas_tank_max;
raw_map = dlmread(input_file,' ',3,0);
[G, node_to_coordinate, gas_list, heart_list, nrow, ncol] = buildgraph(raw_map);
root_node_node_id = coordinate_to_node_idx(root_node_coordinate,nrow);
root_node_name = num2str(root_node_node_id);
large_graph_limit = 400*400;



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
max_edges = nrow*ncol;
s = cell(max_edges,1);
t = cell(max_edges,1);
w = zeros(max_edges,1);
edge_idx = 1;
distances_hearts_to_gas_stations = distances(G,heart_list,gas_list);
distances_root_to_gas_stations = distances(G,root_node_node_id,gas_list);
correspondent_nearest_gas_stations = zeros(length(heart_list),2); % station_idx | dx
for i=1:length(heart_list)
    distance_to_gas_stations = distances_hearts_to_gas_stations(i,:);
    [dx, closest_gas_station] = min(distance_to_gas_stations);
    correspondent_nearest_gas_stations(i,:) = [closest_gas_station, dx];
    radius = gas_tank - dx;
    for j=1:length(gas_list)
        if distance_to_gas_stations(j) < radius
            s(edge_idx) = heart_list_name(i);
            t(edge_idx) = gas_list_name(j);
            w(edge_idx) = distance_to_gas_stations(j);
            edge_idx = edge_idx + 1;
        end
        if ~root_node_to_gas_stations
            distance_root_to_gas = distances_root_to_gas_stations(root_node_node_id,j);
            if  distance_root_to_gas <= gas_tank
                s(edge_idx) = {root_node_name};
                t(edge_idx) = gas_list_name(j);
                w(edge_idx) = distance_root_to_gas;
                edge_idx = edge_idx + 1;
            end
        end
    end
    root_node_to_gas_stations = 1;
end
s = s(1:edge_idx-1);
t = t(1:edge_idx-1);
w = w(1:edge_idx-1);
T = addedge(T,s,t,w);

% Create edges from heart to heart if distance btw them (in original graph
% G) is at most (gas_tank - dx - dy) with dx,dy is the distance to the
% closest gas station of x and y, respectively
% If there're too many hearts, bypass this step since there're too many
% connections to create
if (nrow*ncol < large_graph_limit)
    s = {};
    t = {};
    w = [];
    for i=1:length(heart_list)-1
        if ~findnode(T,heart_list_name(i))
               continue; 
        end
        dx = correspondent_nearest_gas_stations(i,2);
        for j=i+1:length(heart_list)
            if ~findnode(T,heart_list_name(j))
               continue; 
            end
            dy = correspondent_nearest_gas_stations(j,2);
            dxy = distances(G,heart_list(i),heart_list(j));
            if dxy <= (gas_tank - dx - dy)
                s = [s, heart_list_name(i)];
                t = [t, heart_list_name(j)];
                w = [w, dxy];
            end
        end

        if ~root_node_to_hearts
            distance_root_to_heart = distances(G,heart_list(i),root_node_node_id);
            if distance_root_to_heart < (gas_tank - dx)
                s = [s, heart_list_name(i)];
                t = [t, {root_node_name}];
                w = [w, distance_root_to_heart];            
            end
        end
    end
    root_node_to_hearts = 1;
    T = addedge(T,s,t,w);    
end


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
    if (i == length(travel_map) - 1)
        travel_map_T = [travel_map_T p(1:end)]; %#ok<AGROW>
    else
        travel_map_T = [travel_map_T p(1:end-1)]; %#ok<AGROW>
    end
    
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
    current_location = travel_map_T{i};
    if travel_map_T{i}(1) == 'H'
        %In case the next location need to go is Heart, so we need to check the
        %left fuel is enough or not
        if travel_map_T{i+1} (1) == 'H'
            [steps,length_2_nearest_H] = shortestpath(G,str2double(travel_map_T{i}(3:end)),str2double(travel_map_T{i+1}(3:end)));
            idx = find(strcmp(heart_list_name,travel_map_T(i+1)));
            nearest_H_2_nearest_gas = correspondent_nearest_gas_stations(idx, 2);
            if ((length_2_nearest_H + nearest_H_2_nearest_gas +1) < gas_tank)
                [path,fuel_left] = find_steps_and_fuel(G,current_location,travel_map_T{i+1},gas_tank);
                gas_tank = fuel_left(end);
                step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]];
            else
                %goto nearest gas station
                temp = travel_map_T(i);
                idx = find(strcmp(heart_list_name,temp));
                nearest_station_idx = correspondent_nearest_gas_stations(idx,1);
                nearest_station = gas_list_name{nearest_station_idx};
                [path,fuel_left] = find_steps_and_fuel(G,travel_map_T{i},nearest_station,gas_tank);
                step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]]; %#ok<*AGROW>
                gas_tank = fuel_left(end);
                current_location = nearest_station;
                %Goto heart
                [path,fuel_left] = find_steps_and_fuel(G,current_location,travel_map_T{i+1},gas_tank);
                gas_tank = fuel_left(end);
                step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]];
            end
        else
            %If next stop is gas station:
            %   -Check the fuel left before go to gas station
            %   -If fuel left is not enough so go to the nearest gas staion
            %   then go to the index + 1
            [steps,length_2_next_stop] = shortestpath(G,str2double(travel_map_T{i}(3:end)),str2double(travel_map_T{i+1}(3:end)));
            if (length_2_next_stop < gas_tank)
                [path,fuel_left] = find_steps_and_fuel(G,current_location,travel_map_T{i+1},gas_tank);
                gas_tank = fuel_left(end);
                if (i == length(travel_map_T)-1)
                    step_and_fuel = [step_and_fuel; [path(1:end),fuel_left(1:end)]];
                else
                    step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]];       
                end
            else
                 %goto nearest gas station
                temp = travel_map_T(i);
                idx = find(strcmp(heart_list_name,temp));
                nearest_station_idx = correspondent_nearest_gas_stations(idx,1);
                nearest_station = gas_list_name{nearest_station_idx};
                [path,fuel_left] = find_steps_and_fuel(G,travel_map_T{i},nearest_station,gas_tank);
                step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]]; %#ok<*AGROW>
                gas_tank = fuel_left(end);
                current_location = nearest_station;
                %Goto next index
                [path,fuel_left] = find_steps_and_fuel(G,current_location,travel_map_T{i+1},gas_tank);
                gas_tank = fuel_left(end);
                step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]];
            end
        end
    else
        [path,fuel_left] = find_steps_and_fuel(G,current_location,travel_map_T{i+1},gas_tank);
        gas_tank = fuel_left(end);
        if (i == length(travel_map_T)-1)
            step_and_fuel = [step_and_fuel; [path(1:end),fuel_left(1:end)]];
        else
            step_and_fuel = [step_and_fuel; [path(1:end-1),fuel_left(1:end-1)]];       
        end
    end
 
end

result = node_to_coordinate(step_and_fuel(:,1));
result = cell2table(result);
writetable(result,'result.txt','Delimiter',' ','WriteVariableNames',0);

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
