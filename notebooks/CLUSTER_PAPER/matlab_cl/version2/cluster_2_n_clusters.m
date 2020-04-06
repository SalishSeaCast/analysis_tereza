%clear workspace
clear
clc
var = 'ved';
year = '2016';
close all

%load data
ts = strcat('./datamats/',var,'_',year,'hind.mat')
%for VED 
ts = strcat('./datamats/',var,'_',year,'hind_999MASK.mat')
load(ts)

signalmat = datamat;

%for VED
%signalmat = dm_MASK

if strcmp(var,'bio')
    signalmat = datamat_linear;
end

dend_tit = 'Dendrogram';
Zn1 = linkage(signalmat,'ward','euclidean');
n_stn = 580;

figure
dendrogram(Zn1,n_stn)
title(dend_tit, 'FontSize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%part 2
%n = max amount of clusters
n = 100;

%for VED
%n_stn = 570
clusters = zeros(n,n_stn);
w = 'walrus'
for i = 1:100
    cl_n = i
	
    clustermap = cluster(Zn1, 'maxclust', cl_n);

    clusters(i,:) = clustermap';
    
    
end

clearvars -except signalmat datamat clusters

