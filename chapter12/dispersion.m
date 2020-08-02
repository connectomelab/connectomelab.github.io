function [d, dv] = dispersion(matrix, groups)
% [d, dv] = dispersion(matrix,groups)
% Calculates how many groups (e.g. brain regions)
% a single node (e.g. ROI) can on average reach.
% This measure is different from the degree of that node.
% Input:
%    matrix: adjacency matrix (has to be binary and undirected!) 
%    groups: list to which group each node belongs to
% Output:
%    d: average dispersal (mean of dv) 
%    dv: dispersal for all individual nodes
%        (connected groups divided by the number of all groups)
% Author: Marcus Kaiser  Date: 9 August 2011

N = length(matrix);
G = max(groups);

for i=1:N
    dv(i) = length(unique(nonzeros( matrix(i,:) .* groups ))) / G;
end;

d = mean(dv);

return
