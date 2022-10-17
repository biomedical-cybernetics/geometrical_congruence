function run_example

% The file "brain_connectomes_data.mat" contains data related to brain connectomes,
% described in the Methods section of the paper.
% Variables:
% matrices_NOS - cell array of 614 adjacency matrices 114x114 weighted by NOS
% coords_3D - cell array of 614 Euclidean 3D cartesian coordinates of the nodes 114x3
% info - table containing ID, age and gender of the 614 subjects

% load brain connectomes data
load('brain_connectomes_data.mat', 'matrices_NOS', 'coords_3D', 'info')

% compute geometrical congruence (GC)
GC_gsp_NOS = NaN(length(matrices_NOS),1);
GC_gsp_3D = NaN(length(matrices_NOS),1);
for i = 1:length(matrices_NOS)
    % GC-NOS
    xw = matrices_NOS{i};
    xw(xw>0) = 1 ./ (1 + xw(xw>0));
    GC_gsp_NOS(i) = compute_GC_latent(xw, 'original');
    
    % GC-3D
    x = double(matrices_NOS{i}>0);
    geo_3D = squareform(pdist(coords_3D{i}));
    xw = x .* geo_3D;
    GC_gsp_3D(i) = compute_GC_latent(xw, 'original');
end

% gender analysis
labels_gender = strcmp(info.gender,'MALE');
[GC_gsp_NOS_gender_AUPR, GC_gsp_NOS_gender_trustworthiness] = compute_AUPR_trustworthiness(GC_gsp_NOS, labels_gender);
[GC_gsp_3D_gender_AUPR, GC_gsp_3D_gender_trustworthiness] = compute_AUPR_trustworthiness(GC_gsp_3D, labels_gender);

% age analysis
mask_age = info.age<=30 | info.age>=55;
labels_age = info.age(mask_age)<=30;
[GC_gsp_NOS_age_AUPR, GC_gsp_NOS_age_trustworthiness] = compute_AUPR_trustworthiness(GC_gsp_NOS(mask_age), labels_age);
[GC_gsp_3D_age_AUPR, GC_gsp_3D_age_trustworthiness] = compute_AUPR_trustworthiness(GC_gsp_3D(mask_age), labels_age);

% Plot results: the test highlights the better performance using NOS weights with respect to 3D weights
figure('color','white')
tests = {'gender','age'};
subtitles = {'gender (male,female)', 'age-range ([7-30],[55-85])'};
for i1 = 1:length(tests)
    subplot(1,2,i1)
    eval(sprintf('y1 = GC_gsp_NOS_%s_AUPR; t1 = GC_gsp_NOS_%s_trustworthiness;', tests{i1}, tests{i1}))
    eval(sprintf('y2 = GC_gsp_3D_%s_AUPR; t2 = GC_gsp_3D_%s_trustworthiness;', tests{i1}, tests{i1}))
    hold on
    bar(1,y1,'k');
    bar(2,y2,'r');
    set(gca,'ylim',[0.5,0.8],'ytick',0.5:0.1:0.8)
    ylabel('AUPR')
    set(gca,'xlim',[0,3],'xtick',[1,2],'xticklabels',{'GC-NOS', 'GC-3D'})
    if t1<=0.01
        text(1,y1+0.01,'*','FontSize',12)
    end
    if t2<=0.01
        text(2,y2+0.01,'*','FontSize',12)
    end
    axis square
    text(0.5,1.1,subtitles{i1},'units','normalized','horizontalalignment','center')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AUPR, trustworthiness] = compute_AUPR_trustworthiness(scores, labels)

%%% INPUT %%%
% scores - numerical scores for the samples
% labels - binary labels indicating the positive and negative samples

%%% OUTPUT %%%
% AUPR - area under precision-recall
% trustworthiness - p-value of trustworthiness for the AUPR

validateattributes(scores, {'numeric'}, {'vector','finite'})
validateattributes(labels, {'numeric','logical'}, {'vector','binary','numel',length(scores)})
if ~any(labels==1) || ~any(labels==0)
    error('labels cannot be all ones or all zeros')
end
if isrow(scores); scores = scores'; end
if isrow(labels); labels = labels'; end

AUPR = max(compute_AUPR(scores, labels), compute_AUPR(scores, ~labels));
trustworthiness = compute_trustworthiness(scores, labels, AUPR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AUPR = compute_AUPR(scores, labels)

n = length(scores);
[scores,idx] = sort(-scores, 'ascend');
labels = labels(idx);
[~,ut,~] = unique(scores);
ut = [ut(2:end)-1; n];
tp = full(cumsum(labels));
recall = tp ./ sum(labels);
precision = tp ./ (1:n)';
recall = recall(ut);
precision = precision(ut);
if all(recall==1)
    AUPR = precision(1);
else
    AUPR = trapz(recall,precision) / (1-recall(1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trustworthiness = compute_trustworthiness(scores, labels, AUPR)

rstr = RandStream('mt19937ar','Seed',1);
iters = 1000;
AUPR_rand = zeros(iters,1);
for i = 1:iters
    scores_perm = scores(randperm(rstr,length(scores)));
    AUPR_rand(i) = max(compute_AUPR(scores_perm, labels), compute_AUPR(scores_perm, ~labels));
end
trustworthiness = (sum(AUPR_rand >= AUPR) + 1) / (iters + 1);
