input_data = load('/Users/NHua/Documents/01_Work/07_Misfit/10_Misc/12_ShortestPath/input_data.csv');
input_size = size(input_data);
row_len = input_size(1,1);
col_len = input_size(1,2);
% node_list = {};
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

function node_name = build_node_name(row, col, value)
    node_name = sprintf('%d_%d',row,col);
    if (value == 2)
        node_name = ['G_' node_name];   % Gas station
    elseif (value == 3)
        node_name = ['H_' node_name];   % Heart
    end
end