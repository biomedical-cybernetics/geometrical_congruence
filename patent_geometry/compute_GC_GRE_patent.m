function [GC_geo, GC_gsp, GRE_geo, GRE_gsp] = compute_GC_GRE_patent(x, geo)

% code to compute the geometrical congruence (GC) and greedy routing efficiency (GRE)
% in networks with patent geometry
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
% x - adjacency matrix of the network (unweighted, symmetric, zero-diagonal)
% geo - matrix of pairwise geodesics between the nodes (nonnegative, symmetric, zero-diagonal)
%
%%% OUTPUT %%%
% GC_geo - geometrical congruence with respect to geodesics
% GC_gsp - geometrical congruence with respect to geometrical shortest paths
% GRE_geo - greedy routing efficiency with respect to geodesics
% GRE_gsp - greedy routing efficiency with respect to geometrical shortest paths

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

% compute geometrical congruence
mask = triu(full(x==0),1);
[GC_geo, GC_gsp] = meanratio2_mask_mex(geo, gsp, ptsp, mask);
clear ptsp mask;

% compute projection of greedy routing paths
pgrp = compute_pGRP_mex(x, geo);

% compute greedy routing efficiency
mask = full(x==0);
mask(speye(size(mask))==1) = 0;
[GRE_geo, GRE_gsp] = meanratio2_mask_mex(geo, gsp, pgrp, mask);
