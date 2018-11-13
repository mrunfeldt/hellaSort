
% % % USE a Gaussian Mixture model to cluster data
% input: "theData": observation x variable (e.g. [spikes in PC2 ; spikes in PC2]
% "nClusters" 
% MJRunfeldt 03_02_2015

% theData = [u(nodes,1),u(nodes,2)] ;
% nClusters=3;

function [clusterID, clusterz] = gausEMcluster(pcData,nClusters)
%nDim = 4;  % % % OPTION: Restrict Number of dimensions to consider
%dataRedu= pcData(:,1:nDim) ; % % %truncate dimensions

dataRedu = pcData;

gm=fitgmdist(dataRedu,nClusters) ; % fit
clusterID = cluster(gm,dataRedu) ; % cluster

for k = 1:nClusters % assign neurons to clusters
   clusterz{k} = clusterID == k; 
end

end