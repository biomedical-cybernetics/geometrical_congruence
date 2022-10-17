function [ptsp, gsp, pgrp] = compute_paths_patent(x, geo)

% code to compute the geometrical paths in networks with patent geometry
%
%%% INPUT %%%
% x - adjacency matrix of the network (unweighted, symmetric, zero-diagonal)
% geo - matrix of pairwise geodesics between the nodes (nonnegative, symmetric, zero-diagonal)
%
%%% OUTPUT %%%
% ptsp - mean projection of topological shortest paths
% gsp - geometrical shortest paths
% pgrp - projection of greedy routing paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
validateattributes(x, {'numeric','logical'}, {'square','binary'});
if ~issymmetric(x); error('The input matrix ''x'' must be symmetric.'); end
if any(x(speye(size(x))==1)); error('The input matrix ''x'' must be zero-diagonal.'); end
x = double(sparse(x));
validateattributes(geo, {'numeric'}, {'square','ncols',length(x),'finite','nonnegative'});
if ~issymmetric(geo); error('The input matrix ''geo'' must be symmetric.'); end
if any(geo(speye(size(geo))==1)); error('The input matrix ''geo'' must be zero-diagonal.'); end

% compute mean projection of topological shortest paths
tsp = graphallshortestpaths(x, 'Directed', 0);
[~, order] = sort(sum(tsp), 'descend');
ptsp = compute_pTSP_mex(x, geo, tsp, order);
clear tsp order;

% compute geometrical shortest paths
gsp = graphallshortestpaths(x.*geo, 'Directed', 0);

% compute projection of greedy routing paths
pgrp = compute_pGRP_mex(x, geo);