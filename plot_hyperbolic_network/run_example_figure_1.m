% The file "nPSO_example_data.mat" contains data related to a nPSO network,
% generated with parameters N=100, m=4, T=0.1, gamma=2.5, C=5.
% Variables:
% x - adjacency matrix (NxN) of the network
% coords - polar [theta,radius] hyperbolic coordinates (Nx2) of the nodes in the native space
% geo - hyperbolic geodesics NxN
% comm - community memberships Nx1

load('nPSO_example_data.mat', 'x', 'coords', 'geo', 'comm')

% set colors for node communities
colors = NaN(length(comm),3);
colors(comm==1,:) = repmat([255 0 204]./255,sum(comm==1),1);
colors(comm==2,:) = repmat([204 0 255]./255,sum(comm==2),1);
colors(comm==3,:) = repmat([102 0 255]./255,sum(comm==3),1);
colors(comm==4,:) = repmat([0 102 255]./255,sum(comm==4),1);
colors(comm==5,:) = repmat([0 204 255]./255,sum(comm==5),1);

figure('color','white')

% plot network
subplot(1,2,1)
plot_hyperbolic_network(x, coords, 'native', colors)

% plot geometrical congruence histogram
subplot(1,2,2)
[ptsp, gsp, pgrp] = compute_paths_patent(x, geo);
mask = triu(full(x==0),1);
y = geo(mask)./ptsp(mask);
xi = [0,0.05:0.1:0.95,1.01];
yi = NaN(1,length(xi)-1);
for i = 1:length(xi)-1
   yi(i) = mean(y>=xi(i) & y<xi(i+1)); 
end
xi = [0,0.05:0.1:0.95,1];
histogram('BinEdges',xi,'BinCounts',yi,'FaceColor',[179 179 179]./255)
hold on
xlim([-0.05 1.05])
ylim([0 0.35])
xlabel('GEO/pTSP')
ylabel('probability')
set(gca,'xtick',0:0.1:1,'ytick',0:0.05:0.35,'FontSize',8)
box on