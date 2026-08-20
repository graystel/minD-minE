function [branchData, adj, edgeMatrix] = loopremover(branchData, adj, edgeMatrix)
    [numNodes, ~, ~] = size(adj);
    todelete = [];
    for i=1:numNodes
        for j=1:i-1 %between each pair of nodes
            edges = squeeze(adj(i, j, :)); %track the branches in the third dim for this particular pair
            edges = edges(edges > 0); %we want to ensure that the edges are nonzero

            if length(edges) > 1
                lengths = zeros(length(edges), 1); %make a list of the lengths of all the branches connecting these two nodes
                for k=1:length(edges)
                    lengths(k) = size(branchData{edges(k)}, 1);
                end

                [~, sortIdx] = sort(lengths, "descend"); %since we wanna keep either the shortest or the longest (will experiment to see
                %which one leads to better results)
                sortededges = edges(sortIdx);   
                todelete = [todelete; sortededges(2:end)]; %keep the first, append the rest (loop/redundant) to a list of branches 
                % to delete all at the same time 
            end
        end
    end

    todelete = unique(todelete); %no repeats

    %first clear these from adj
    for k = 1:length(todelete)
        adj(adj == todelete(k)) = 0;
    end

    %delete from branchData in descending order for preservation
    alltodelete = sort(todelete, 'ascend');
    for k = 1:length(alltodelete)
        delIdx = alltodelete(k);
        branchData(delIdx) = [];
        edgeMatrix(delIdx, :) = [];
        %branches later than the current branch get decreased by 1 in adj to keep branches in branchdata and adj synchronized
        adj(adj > delIdx) = adj(adj > delIdx) - 1;
    end
end
    
