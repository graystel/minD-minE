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

                [~, sortIdx] = sort(lengths, "ascend"); %since we wanna keep either the shortest or the longest (will experiment to see
                %which one leads to better results)
                sortededges = edges(sortIdx);   
                todelete = [todelete; sortededges(2:end)]; %keep the first, append the rest (loop/redundant) to a list of branches 
                % to delete all at the same time 
            end
        end
    end

    %for singlenode loops
    for i=1:numNodes
        selfloops = squeeze(adj(i, i, :));
        selfloops = selfloops(selfloops > 0);
        if ~isempty(selfloops)
            todelete = [todelete; selfloops(:)];
        end
    end

    todelete = unique(todelete); %no repeats

    %first clear these from adj
    for k = 1:length(todelete)
        adj(adj == todelete(k)) = 0;
    end

    %delete from branchData in descending order for preservation
    alltodelete = sort(todelete, 'descend');
    for k = 1:length(alltodelete)
        delIdx = alltodelete(k);
        branchData(delIdx) = [];
        edgeMatrix(delIdx, :) = [];
        %branches later than the current branch get decreased by 1 in adj to keep branches in branchdata and adj synchronized
        adj(adj > delIdx) = adj(adj > delIdx) - 1;
    end

    edgeMatrix(:, 1) = (1:length(branchData));

    %using minspantree to remove loops of multiple nodes
    % branches = edgeMatrix(:, 1);
    % starts   = edgeMatrix(:, 2);
    % ends     = edgeMatrix(:, 3);
    % 
    % branchlengths = zeros(length(branches), 1);
    % for i = 1:length(branches)
    %     branchlengths(i) = size(branchData{branches(i)}, 1);
    % end
    % 
    hasNode72 = sum(edgeMatrix(:, 2) == 72 | edgeMatrix(:, 3) == 72);
    fprintf('Node 72 appears %d times in edgeMatrix before MST.\n', hasNode72);
    % originaltable = table([starts, ends], branchlengths, branches, 'VariableNames', {'EndNodes', 'Weight', 'branches'});
    % untrimmedgraph = graph(originaltable);
    % 
    % trimmedgraph = minspantree(untrimmedgraph);
    % startvertices = trimmedgraph.Edges.EndNodes(:,1);
    % endvertices = trimmedgraph.Edges.EndNodes(:, 2);
    % keptbranches = trimmedgraph.Edges.branches;
    % 
    % newbranchData = cell(length(keptbranches), 1);
    % newadj = zeros(length(adj), length(adj));
    % newedgeMatrix = zeros(length(keptbranches), 3);
    % 
    % for i=1:length(keptbranches)
    %     a = startvertices(i);
    %     b = endvertices(i);
    %     branchnumber = keptbranches(i);
    % 
    %     newbranchData{i} = branchData{branchnumber};
    %     newadj(a, b) = i;
    %     newadj(b, a) = i;
    %     newedgeMatrix(i,:) = [i a b];
    % 
    % end

end
    
