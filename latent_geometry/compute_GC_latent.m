function GC_gsp = compute_GC_latent(x, weights)

% code to compute the geometrical congruence (GC) in networks with latent geometry
%
% Authors:
% Alessandro Muscoloni, 2022-02-07
%
% Reference:
% "Geometrical congruence and efficient greedy navigability of complex networks"
% C. V. Cannistraci, A. Muscoloni, arXiv:2005.13255, 2020
% https://arxiv.org/abs/2005.13255
%
% Released under MIT License
% Copyright (c) 2022, C. V. Cannistraci, A. Muscoloni

%%% INPUT %%%
% x - adjacency matrix of the network (unweighted or weighted, symmetric, zero-diagonal)
% weights - string indicating the type of weights to use:
%   'original' - the adjacency matrix in input is assumed weighted and the associated weights are used,
%       they should represent distances between the node pairs (and not similarities)
%   'RA' - Repulsion-Attraction rule is used to weight the adjacency matrix (weights given in input are discarded)
%   'EBC' - Edge-Betweenness-Centrality rule is used to weight the adjacency matrix (weights given in input are discarded)
%
%%% OUTPUT %%%
% GC_gsp - geometrical congruence with respect to geometrical shortest paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
validateattributes(x, {'numeric'}, {'square','finite','nonnegative'});
if ~issymmetric(x); error('The input matrix ''x'' must be symmetric.'); end
if any(x(speye(size(x))==1)); error('The input matrix ''x'' must be zero-diagonal.'); end
x = double(sparse(x));
validateattributes(weights, {'char'}, {'scalartext'});
if ~ismember(weights,{'original','RA','EBC'})
    error('Possible weights options: ''original'', ''RA'', ''EBC''.');
end

% set weights
if strcmp(weights,'original')
    xw = x;
    x = double(x>0);
elseif strcmp(weights,'RA')
    x = double(x>0);
    xw = RA_weighting(x);
elseif strcmp(weights,'EBC')
    x = double(x>0);
    xw = EBC_weighting(x);
end

% compute mean projection of topological shortest paths
tsp = graphallshortestpaths(x, 'Directed', 0);
[~, order] = sort(sum(tsp), 'descend');
ptsp = compute_pTSP_mex(x, full(xw), tsp, order);
clear tsp order;

% compute geometrical shortest paths
gsp = graphallshortestpaths(xw, 'Directed', 0);

% compute geometrical congruence
mask = triu(full(x==0),1);
GC_gsp = mean(gsp(mask)./ptsp(mask));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xw = RA_weighting(x)

n = size(x,1);
cn = full(x*x);
ext = repmat(full(sum(x,2)),1,n) - cn - x;
xw = x .* (1 + ext + ext') ./ (1 + cn);

function xw = EBC_weighting(x)

[~, xw] = betweenness_centrality(x);
