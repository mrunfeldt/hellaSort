

function [newNeurons] = hellaSort_mergeOption(nrns,spiketimes)
% Cluster 1 - neuron A, B, C, or X
newNeurons = cell(1,4)  ;
isiBins= logspace(-2,4.5,80);

interBin = 0.5; % ms, bin rate
rateEdges = 0:interBin:max(spiketimes) ; %linspace(0,max(spiketimes)*1e-3,300); % bins for rate over time

ff=figure; xwidth = 16e2; ywidth=4e2; ccol = 'mbrg' ;
set(ff, 'Position', [0 0 xwidth ywidth]);
movegui(ff,'north');tHand = title([]); drawnow

nSubs = length(nrns)+1 ; % # of subplots (1summary, 1 autoCorr, rest xCorr)
pearCount = 1:length(nrns);




for w = 1:length(nrns) % for each cluster
    
        thiscluster = nrns{w} ; mySpikes = spiketimes(thiscluster) ; 
        rate = histc(mySpikes,rateEdges) ./ interBin ;% firing rate++
% % % % % % % % Plot ISI histogram % % % % % % % % % % % % % %        
subplot(1,nSubs,2);hold on;
isi = diff(mySpikes);isi1 = 100 * sum(isi<=1)./length(isi);
bins=histc(isi,isiBins) ; leg{w} = [ 'Cluster ', num2str(w)];
plot(log10(isiBins),bins./sum(bins),'k','linewidth',2,'color',ccol(w));
xlim([log10(isiBins(1)) log10(isiBins(end-1))])
binPlot = isiBins(2:10:end) ; xlabel('ISI (ms) - logscale');ylabel('Count')
set(gca,'xtick',log10(binPlot),'xticklabel',num2cell(round(binPlot.*10)./10))%'XScale','log')
msMark = isiBins(isiBins >1); msMark = log10(msMark(1));
line([msMark msMark],[0 max(bins./sum(bins))],'linewidth',3)
ylabel('Probability'); xlabel('ISI (ms)','fontsize',16,'FontWeight','bold')
drawnow        
        
        
        allOthers = pearCount; allOthers(w)=[];
        % % % CrossCorr % % % 
        for q = 1:length(allOthers) % xcorr for all pairings
            clu1 = spiketimes(nrns{w}) ;  clu2 = spiketimes(nrns{allOthers(q)}) ;  
            rate1 = histc(clu1,rateEdges) ./ interBin ;% firing rate
            rate2 = histc(clu2,rateEdges) ./ interBin ;% firing rate
            [crxCorr,lags] = xcorr(rate1,rate2,20);
            subplot(1,nSubs,2+q,'replace');hold on
            plot(lags.*interBin,crxCorr,'linewidth',3,'color',ccol(q))
            tx=text(-10,median(crxCorr)*0.7,['Cluster ',num2str(w),' spikes sooner']);
            set(tx,'fontweight','bold','fontsize',11)
            
            tx=text(3,median(crxCorr)*1.15,['Cluster ',num2str(allOthers(q)),' spikes sooner']);
            set(tx,'fontweight','bold','fontsize',11)
            
            xlabel(['CrossCorr with Cluster ',num2str(allOthers(q)),' (ms)'],'fontsize',16,'FontWeight','bold')
            ylim([-3 max(crxCorr)+2]); drawnow
        end
        
subplot(1,nSubs,3);
title(['Cluster ',num2str(w),' belongs to which neuron?  Hit "a", "b", "c", or "x" to toss' ],'fontsize',25,'color',[0 0.6 0.2]);
if size(nrns{w},1) ~=1 ; nrns{w} = nrns{w}';end

waitforbuttonpress; sVal=double(get(ff,'CurrentCharacter'))  ;
switch sVal
    case 97 % "a"
        num = 'Neuron A';
        newNeurons{1}  = [newNeurons{1} nrns{w}] ;
    case 98 % "b"
        num = 'Neuron B';
        newNeurons{2}  = [newNeurons{2} nrns{w}] ;
    case 99 % c
        num = 'Neuron C';
        newNeurons{3}  = [newNeurons{3} nrns{w}] ;
    case 120 % x
        num = 'NOT a neuron';
end


subplot(1,nSubs,1);hold on
text(0.2,w, ['Cluster ',num2str(w),' is ',num ],'fontsize',14,'fontweight','bold')
set(gca,'xtick',[],'ytick',[]);xlim([0 1]);ylim([0.5 nSubs+0.5])
hold off
end

empties = find(cell2mat(cellfun(@(x) isempty(x),newNeurons,'un',0))==1) ;
newNeurons(empties)=[];
end % END function 


