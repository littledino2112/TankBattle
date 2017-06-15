function [G, node_to_coordinate, gas_station_list, heart_list, nrow, ncol] = buildgraph(input_data)
% This builds a graph G from input matrix
[nrow,ncol] = size(input_data);
max_edges = nrow*ncol*2; % Each node can connect to 2 others nodes (down and right)
max_nodes = nrow*ncol;
s = zeros(1,max_edges);
d = zeros(1,max_edges);
node_to_coordinate = cell(max_nodes,1); % Format: node_idx -> {[coordicate]}
edge_idx = 1;
heart_list = [];
gas_station_list = [];
for i=1:nrow
   for j=1:ncol
      if (input_data(i,j) ~= 0)
        node_idx = nrow*(i-1)+j;
        node_to_coordinate(node_idx,:) = {[i,j]};
        if (input_data(i,j) == 2)
            gas_station_list = [gas_station_list node_idx]; %#ok<AGROW>
        elseif (input_data(i,j) == 3)
            heart_list = [heart_list node_idx]; %#ok<AGROW>
        end
        if (j < ncol)
            if (input_data(i,j+1) ~= 0)
            %                 node_name_1 = build_node_name(i,j+1,input_data(i,j+1));
            %                 G = addedge(G,node_name,node_name_1);
                % Create connection btw node (i,j) to (i,j+1)
                s(edge_idx) = nrow*(i-1)+j;
                d(edge_idx) = nrow*(i-1)+(j+1);
                edge_idx = edge_idx + 1;
            end
        end
        if (i < nrow)
            if (input_data(i+1,j) ~= 0)
            %                 node_name_2 = build_node_name(i+1,j,input_data(i+1,j));
            %                 G = addedge(G,node_name,node_name_2);
                % Create connection btw node (i,j) to (i+1,j)            
                s(edge_idx) = nrow*(i-1)+j;
                d(edge_idx) = nrow*(i)+j;
                edge_idx = edge_idx + 1;
            end
        end
      end
   end
end

s = s(1:edge_idx-1);
d = d(1:edge_idx-1);
G = graph(s,d);
end

function node_name = build_node_name(row, col, value)
    node_name = sprintf('%d_%d',row,col);
    if (value == 2)
        node_name = ['G_' node_name];   % Gas station
    elseif (value == 3)
        node_name = ['H_' node_name];   % Heart
    end
end

function distance = compute_distance_btw_hearts(graph, heart_list)
    d_heart = distances(graph,heart_list,heart_list);
    distance = d_heart;
end