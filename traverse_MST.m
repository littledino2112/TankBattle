start_node = 'H_2_49';
global travel_map;
travel_map = {};
mintree = minspantree(H);
directed_mintree = shortestpathtree(mintree,start_node);
traverse_mst(directed_mintree,start_node);

% Go back from map to original graph G
travel_map_T = {};
for i=1:length(travel_map)-1
    p = shortestpath(T,travel_map{i},travel_map{i+1});
    travel_map_T = [travel_map_T p(1:end-1)];
end

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