%clear workspace
clear
clc
%cluster_from_matfile(var)
var = 'nit';
year = '2015';
close all

%load data
ts = strcat('./datamats/',var,'_',year,'hind.mat')
load(ts)

signalmat = datamat;

if strcmp(var,'bio')
    signalmat = datamat_linear;
end

dend_tit = 'Dendrogram';
Zn1 = linkage(signalmat,'ward','euclidean');
figure
dendrogram(Zn1,580)
title(dend_tit, 'FontSize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%part 2
% 
cl_n = 10; 
clustermap = cluster(Zn1, 'maxclust', cl_n);

%clearvars -except clustermap cl_n datamat datamat_linear

