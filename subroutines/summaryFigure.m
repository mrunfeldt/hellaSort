

% events - original waveforms
% spiketimes - original spiketimes (in ms)
% time - "T" for waveform Plots
% nrns - cell, where each cell is a list of indicies that assign spikes to
% neuron identities
% MJRunfeldt 03_11_2015

function [exitVal,nrns,confusion] = summaryFigure(events,spiketimes,time,nrns)
%time = T; % for devel

pickN = 20 ; % plot a subset of neurons in cluster
durBins = linspace(min(spiketimes),max(spiketimes),40).*1e-3; % bins for waveforms over time
edges = linspace(0,max(spiketimes)*1e-3,40); % bins for rate over time
interBin = diff(edges); interBin=interBin(1);

tryAgain = 1;
while tryAgain == 1 
nNodes = size(events,1) ; myColors = 'crbg' ;
sumFig = figure;xwidth = 15e2; ywidth=9e2;
set(sumFig, 'Position', [0 0 xwidth ywidth]);movegui(sumFig,'center');

try% may or may not need to transpose
weAreNeurons = unique(cell2mat(nrns'))';
catch;  weAreNeurons = unique(cell2mat(nrns))';
end

unwanted = setdiff(1:nNodes,weAreNeurons) ; % events not sorted


% % % Plot MEan +/- STD waveforms % % 
summary=[]; legHand=cell(1,length(nrns)) ; cHand=legHand;
midHand=subplot(3,3,[2,5]);hold on
for sf = 1:length(nrns) % for each cluster
    summary = [summary length(nrns{sf})] ; %cHand{sf}) = 
    legHand{sf} = ['Cluster ', num2str(sf)];
    waveMean = mean(events(nrns{sf},:)) ; waveSTD = std(events(nrns{sf},:)); % mean and STD of waveform
    clusterForms{sf} = events(nrns{sf},:) ; % all waveforms for cluster
    shadedLines(time,waveMean, waveSTD,myColors(sf))
end
set(gca,'xlim',[time(1) time(end)]);xlabel('Time (sec)');


% % % PLOT PIE CHART % % % 
summary = [summary length(unwanted)]; legHand{sf+1} = 'Not included' ;
subplot(3,3,1); pHand = pie(summary,legHand); 
for k = 1:length(nrns)
set(pHand(k*2-1),'facecolor',myColors(k));
set(pHand(k*2),'String',['Cluster ', num2str(k)],'fontsize',18);
end

% % % Plot randomSampling of waveforms % % ONLY ALLOW UP TO 3
if length(nrns)> 3 ; plotSub = 3;else plotSub = length(nrns);end
for n = 1:plotSub % for each cluster
    thiscluster = nrns{n} ; % indicies
    if length(thiscluster) < pickN;pickNN = length(thiscluster);else pickNN = pickN; end
    whyNotUs = randperm(length(thiscluster),pickNN) ; % random selection of spikes
    subplot(3,3,n*3);hold on;plot(time,events(thiscluster(whyNotUs),:))
    plot(time,mean(events(thiscluster,:)),'color',myColors(n),'linewidth',3)
    set(gca,'xlim',[time(1) time(end)]);xlabel('Time (sec)');
    title(['Cluster ',num2str(n)],'color',myColors(n),'fontsize',18)
    %ylim([min(min(events)) max(max(events))-(1*mean(std(events)))])
end


% % % PLOT RATES % % %      
newLeg = [];
subplot(3,3,4);hold on;
for ss = 1:sf
    newLeg{ss} = ['Cluster ', num2str(ss)];
plot(edges,histc(spiketimes(nrns{ss}).*1e-3,edges)./interBin,myColors(ss),'linewidth',3);
end
xlabel('Time (sec)');title('Binned Spike Rate');ylabel('Rate (Hz)');
set(gca,'xlim',[0 edges(end-1)]);legend(newLeg)

set(get(midHand,'title'),'String','...Be patient while quality metrics run...','fontsize',17,'color',[0 0.8 0.2])
drawnow; 


% % % Plot WaveForm Descriptor % % % 
subplot(3,3,7,'replace');hold on; 
ylabel('Max Amplitude');set(gca,'xlim',[0 durBins(end-1)]); xlabel('Time (sec)');
for w = 1:length(nrns) % for each cluster
        thiscluster = nrns{w} ;  ourWaves = events(thiscluster,:);
        ourSpikes = spiketimes(thiscluster).*1e-3 ;
duration=[] ;maxAmp=[];
for g = 1:length(ourSpikes)
    [duration(g),maxAmp(g)] = waveForm_descript(ourWaves(g,:),time); 
end
%[duration(g),maxAmp(g),timeToPeak(g),preMin(g),postMin(g),,firstMin(g),secondMin(g)] = waveForm_descript(ourWaves(g,:),time);


%plot(ourSpikes,yvar,'.','color',myColors(w),'linewidth',8) % all
yvar = maxAmp ; [~,binDur] = histc(ourSpikes,durBins);wine = [];
for f = 1:length(durBins);wine(f) = mean(yvar(binDur ==f)) ;end
plot(durBins,wine,'-','color',myColors(w),'linewidth',3)  % binned mean

%plot(maxAmp,duration,'.','color',myColors(w))
end


drawnow

% % % % Quality Metrics % % % 
% row1 = false positives; row2 = false negatives. Columns == clusters
texLocs = 0;pears = combnk(1:length(nrns),2) ; % all pairwise combinations of clusters
subplot(3,3,8,'replace'); set(gca,'color',[0.2 0.2 0.2])
if length(nrns)>1 & length(nrns)<4
    title('Quality Metrics');hold on;
    
    for a = 1:size(pears,1)
        confusion{a} = gaussian_overlap( clusterForms{pears(a,1)} , clusterForms{pears(a,2)} ); 
        % w1 = clusterForms{pears(a,1)} ; w2 = clusterForms{pears(a,2)}  ;
        texLocs=texLocs+7;
    text(0.3,texLocs+4.5, ['False Positive, Cluster ',num2str(pears(a,2)),' assigned to  ',num2str(pears(a,1)),' :   ', ...
        num2str(confusion{a}(1,1))],'color',myColors(a),'fontsize',12,'FontWeight','bold') 
    text(0.3,texLocs+3, ['False Negative, Cluster ',num2str(pears(a,1)),' assigned to ',num2str(pears(a,2)),' :  ', ...
        num2str(confusion{a}(1,2))],'color',myColors(a),'fontsize',12,'FontWeight','bold')  
    text(0.3,texLocs+1.5, ['False Positive, Cluster ',num2str(pears(a,1)),' assigned to  ',num2str(pears(a,2)),' :   ', ...
        num2str(confusion{a}(2,1))],'color',myColors(a),'fontsize',12,'FontWeight','bold')  
    text(0.3,texLocs, ['False Negative, Cluster ',num2str(pears(a,2)),' assigned to ',num2str(pears(a,1)),' :  ', ...
        num2str(confusion{a}(2,2))],'color',myColors(a),'fontsize',12,'FontWeight','bold') 
    end
    ylim([4 texLocs+7]);xlim([0.25 0.9])
else confusion=[];
end
%confusion=[];
set(get(midHand,'title'),'String','Hit "y" to keep, "x" to toss or merge, "t" to start again, "n" to cluster further','fontsize',20,'color',[0.7 0.2 0.4] )

waitforbuttonpress; exitVal=double(get(sumFig,'CurrentCharacter'))  ;

if exitVal == 120 %to toss or merge
            movegui(sumFig,'south'); drawnow
            [newNeurons] = hellaSort_mergeOption(nrns,spiketimes) ;
            nrns = newNeurons ;
            tryAgain = 1; close
    
elseif exitVal == 121 | exitVal ==116  | exitVal==110 % close exit
    tryAgain = 0; close
else set(get(midHand,'title'),'String','What?!?!? Hit RETURN','fontsize',20,'color',[0.7 0.2 0.4] )
    pause; tryAgain = 1;
end

end % END while tryAgain == 1

end

