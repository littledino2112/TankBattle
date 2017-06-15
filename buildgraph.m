function [G, gas_station_list, heart_list] = buildgraph()
% This builds a graph G from input matrix
input_data = load('input_data.csv');
input_size = size(input_data);
row_len = input_size(1,1);
col_len = input_size(1,2);
G = graph; % Create empty graph
for i=1:row_len
   for j=1:col_len
      if (input_data(i,j) ~= 0)
          node_name = build_node_name(i,j,input_data(i,j));
          if (j < col_len)
            if (input_data(i,j+1) ~= 0)
                node_name_1 = build_node_name(i,j+1,input_data(i,j+1));
                G = addedge(G,node_name,node_name_1);
            end
          end
          if (i < row_len)
            if (input_data(i+1,j) ~= 0)
                node_name_2 = build_node_name(i+1,j,input_data(i+1,j));
                G = addedge(G,node_name,node_name_2);
            end
          end
      end
   end
end

% Find the indexes of Hearts in node array
node_names = G.Nodes.Name;
heart_list = [];
gas_station_list = [];
for idx=1:length(node_names)
    if node_names{idx}(1) == 'H'
        heart_list = [heart_list idx]; %#ok<AGROW>
    elseif node_names{idx}(1) == 'G'
        gas_station_list = [gas_station_list idx]; %#ok<AGROW>
    end
end
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