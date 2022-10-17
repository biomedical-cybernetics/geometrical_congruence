% The file "nPSO_example_data.mat" contains data related to two nPSO networks,
% generated with parameters N=100, m=6, T=0.1, C=4,
% gamma=2 for the network #1 and gamma=3 for the network #2.
% Variables:
% x1, x2 - adjacency matrix NxN
% coords1, coords2 - hyperbolic polar coordinates Nx2 (angle, radius)
% geo1, geo2 - hyperbolic geodesics NxN
% comm1, comm2 - community memberships Nx1

% Test #1: gamma=2
load('nPSO_example_data.mat', 'x1', 'geo1')
[GC_geo1, GC_gsp1, GRE_geo1, GRE_gsp1] = compute_GC_GRE_patent(x1, geo1);

% Test #2: gamma=3
load('nPSO_example_data.mat', 'x2', 'geo2')
[GC_geo2, GC_gsp2, GRE_geo2, GRE_gsp2] = compute_GC_GRE_patent(x2, geo2);

% Plot results: the test highlights the large decrease of GC_geo and GRE_geo
% when going from gamma=2 to gamma=3.
figure('color','white')
measures = {'GC_geo', 'GRE_geo', 'GC_gsp', 'GRE_gsp'};
subtitles = {'GC($\overline{pTSP},GEO$)', 'GRE($pGRP,GEO$)', 'GC($\overline{pTSP},GSP$)', 'GRE($pGRP,GSP$)'};
for i = 1:length(measures)
    subplot(2,2,i)
    eval(['y1 = ' measures{i} '1; y2 = ' measures{i} '2;'])
    hold on
    bar(1,y1,'k');
    bar(2,y2,'r');
    set(gca,'ylim',[0.5,1],'ytick',0.5:0.1:1)
    ylabel(subtitles{i},'Interpreter','Latex')
    set(gca,'xlim',[0,3],'xtick',[1,2],'xticklabels',{'\gamma = 2','\gamma = 3'})
    text(1,y1+0.05,sprintf('%.2f',y1),'horizontalalignment','center')
    text(2,y2+0.05,sprintf('%.2f',y2),'horizontalalignment','center')
end
suptitle('nPSO networks (N=100, m=6, T=0.1, C=4)')
