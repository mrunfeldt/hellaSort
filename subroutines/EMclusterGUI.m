
% % % "EMclusterGUI" - GUI for using GAussian Mixture model to cluster data
% % % % % calls "gausEMcluster" [[can be replaced with different cluster method % % 
% updated: MJRunfeldt 2015_03_10

% "clustVal": waitforbuttonpress; clustVal=double(get(f2H,'CurrentCharacter')) ;
% "theData": spikes/nodes in PC space (nNodes x nPCs)
function [inClust,decision] = EMclusterGUI(clustVal,waveForms,spikeTimes,T) 

%waveForms = events(keepers,:) ;  spikeTimes = spiketimes(keepers);% for devel

% % %  PARAMETERS % % % 
isiBins= logspace(-2.5,4.5,60); 
perVariance = 60 ;
durBins = linspace(min(spikeTimes),max(spikeTimes),50).*1e-3;
rateBins = linspace(0,max(spikeTimes)*1e-3,100); 

% % % OBTAIN WaveForm Descriptors % % % 
for g = 1:length(spikeTimes)
[duration(g),maxAmp(g),timeToPeak(g),preMin(g),postMin(g)] = waveForm_descript(waveForms(g,:),T);
end

 % % GENERATE NEW PCA SPACE based on waveforms % % % %
[~,pcData,eigVal] = princomp(waveForms); eigVal = 100.* eigVal ./ sum(eigVal) ;
%figure;bar(eigVal);xlabel('Eigenvector');ylabel('Percent variance')

x=find(cumsum(eigVal) > perVariance); nPCs = x(1) ; % # of PC's to capture XX% variance
input = pcData(:,1:nPCs); % reduced dimensionality
% % % AUGMENT PC Space with waveform descriptors % % %     
input = [duration' maxAmp' preMin' postMin' input]; % construct Gaus Mix Model from this

keepTrying = 1; 
while keepTrying == 1  
        
inClust =[];clusterIDs =[]; nClusters = double(clustVal-48) ; % convert key into number
% % initialize figure % %        
emH = figure;xwidth = 18e2; ywidth=8e2;set(emH, 'Position', [0 0 xwidth ywidth])
movegui(emH,'center')

try % TRY Gaussian mixture model - if it doesn't fail, plot outcome

    clusterIDs =  gausEMcluster(input,nClusters) ; 
    tossMerge = 0;myColors = 'kcmgbry';
        for j = 1:nClusters
        inClust{j} = find(clusterIDs==j);
        end
    
catch 
        clusterIDs{1} = 1; nClusters = 1; % give it something to not crash
        disp('gauss mix model fail')
        tH=subplot(2,4,1,'replace');
        set(get(tH,'title'),'String','Gaus Model Fail,"n" to try another # of clusters ,"t" to move on','fontsize',15,'color','r');
        tossMerge=1;drawnow
        clustVal = 49;
end % END try/catch for Gaus mix model    


while tossMerge == 0    
% % % PLOT IN PCA SPACE % % % 
tH=subplot(2,4,1,'replace');hold on;set(gca,'Color',[0.8 0.8 0.8])
legger = [];
    for i = 1:nClusters
        %inClust{i} = find(clusterIDs==i);
        plot(pcData(inClust{i},1),pcData(inClust{i},2),'.','color',myColors(i),'markersize',6);
        legger{i} = ['Cluster ', num2str(i)];
    end;
    
% % % PLOT WAVEFORMS+/- STD% % % 
    hTwo=subplot(2,4,2,'replace');hold on
    for ii = 1:nClusters
    waveMean = mean(waveForms(inClust{ii},:)) ; waveSTD = std(waveForms(inClust{ii},:)); 
    shadedLines(T,waveMean,waveSTD,myColors(ii));
    end ;set(gca,'yticklabel',[]);drawnow
        
% % % PLOT MEAN WAVEFORMS % % % 
    midH=subplot(2,4,3,'replace');hold on
    for x = 1:nClusters
    waveMean = mean(waveForms(inClust{x},:)) ; plot(T,waveMean,myColors(x),'linewidth',4)
    end ; 
    xlim([T(1) T(end)]);leg=legend(legger);set(leg,'fontsize',15,'color',[1 1 1]);drawnow
          
% % % PLOT ISI DISTRIBUTIONS % % %     
    rtH=subplot(2,4,4,'replace'); hold on; title('Inter Spike Intervals');
    for v = 1:nClusters
        isi = diff(spikeTimes(inClust{v}));bins=histc(isi,isiBins) ; 
        plot(log10(isiBins),bins./sum(bins),'color',myColors(v),'linewidth',2);
    end
    xlim([log10(isiBins(1)) log10(isiBins(end))])
    binPlot = isiBins(2:10:end) ; xlabel('ISI (ms) - logscale');ylabel('Count')
    set(gca,'xtick',log10(binPlot),'xticklabel',num2cell(round(binPlot.*10)./10))%'XScale','log')
    msMark = isiBins(isiBins >1); msMark = log10(msMark(1));
    line([msMark msMark],[0 max(bins./sum(bins))],'linewidth',3)
        
% % % PLOT RATES % % %      
    subplot(2,4,5,'replace');hold on;interBin = diff(rateBins); interBin=interBin(1);rates=[];
    for g = 1:nClusters
        rates{g} = histc(spikeTimes(inClust{g}).*1e-3,rateBins)./interBin ;
    plot(rateBins,rates{g},'color',myColors(g),'linewidth',3);
    end;set(gca,'xlim',[0 rateBins(end-2)]);
    xlabel('Time (sec)');title('Binned Spike Rate');ylabel('Rate (Hz)');     
    
%     set(get(tH,'title'),'String','Hit "z" for corrs and amplitudes');
%     waitforbuttonpress; more=double(get(emH,'CurrentCharacter')) ;    
%    if more== 122

% % % SPIKE RATE CORR % % %     
subplot(2,4,6,'replace');hold on;    
pears = combnk(1:nClusters,2); cnt = 0;
    for j = 1:size(pears,1)
        n1 = rates{pears(j,1)} ; n2 = rates{pears(j,2)} ;
        [r,p]=corrcoef(n1,n2); r=r(2); p=p(2);cnt = cnt +1;
        text(0.2, cnt , ['N ',num2str(pears(j,1)),' vs N ',num2str(pears(j,2)),...
            ', r = ', num2str(r),', p = ', num2str(p)],'fontsize',14,'fontweight','bold')
    end
ylim([0.1 cnt+1]);xlim([0 1.3])


% % % PLOT Waveform Descriptors over time % % %      
subplot(2,4,7,'replace');hold on; pickMe = maxAmp; ylabel('Max Amplitude')
    for d = 1:nClusters
        spikeRound = spikeTimes(inClust{d}).*1e-3;yvar = pickMe(inClust{d});
         [~,binDur] = histc(spikeRound,durBins);wine = [];
        for f = 1:length(durBins);wine(f) = mean(yvar(binDur ==f)) ;end
        plot(durBins,wine,'-','color',myColors(d),'linewidth',3)    
    end ; set(gca,'xlim',[0 durBins(end-1)]); xlabel('Time (sec)');
    
subplot(2,4,8,'replace');hold on; pickMe = duration; ylabel('Spike Duration (ms)')
    for d = 1:nClusters
        spikeRound = spikeTimes(inClust{d}).*1e-3;yvar = pickMe(inClust{d});
         [~,binDur] = histc(spikeRound,durBins);wine = [];
        for f = 1:length(durBins);wine(f) = mean(yvar(binDur ==f)) ;end
        plot(durBins,wine,'-','color',myColors(d),'linewidth',3)    
    end ; set(gca,'xlim',[0 durBins(end-1)]); xlabel('Time (sec)');
 
%    end % END option for corrs and waveform amps

set(get(hTwo,'title'),'String','Hit "y" to accept, "n" try another # of clusters, ','fontsize',17,'color','r');        
set(get(rtH,'title'),'String','"m" to manually cluster, "t" start from the top','fontsize',17,'color','b');    
set(get(midH,'title'),'String','Hit "x" to toss or merge','fontsize',18,'color','k');

waitforbuttonpress; decision=double(get(emH,'CurrentCharacter')) ;


switch decision 
    case 110 % "n" try again 
        set(get(hTwo,'title'),'String','How Many Clusters in model ?','fontsize',20,'color',[0.7 0 0.3])
        waitforbuttonpress; clustVal=double(get(emH,'CurrentCharacter')) ;
        set(get(midH,'title'),'String','Be patient while EM runs','fontsize',18,'color',[0.2 0.8 0.2]);drawnow
        keepTrying = 1; tossMerge = 1; close; inClust=[];
    case 116 %"t" exit (try again from the top)
        inClust=[]; keepTrying=0; tossMerge = 1;% break from function
    case 109 % "m" exit (move on to manual clustering
        inClust=[]; keepTrying=0; tossMerge = 1;% break from function
    case 121 %"yes" keep "nrns" assignment and exit
        keepTrying =0; tossMerge = 1; movegui(emH,'northwest')
    
    case 120 % "x" Remove or merge a cluster
      
set(get(hTwo,'title'),'String','Hit "x" to remove a cluster, "m" to merge, ','fontsize',20,'color','r');        
set(get(midH,'title'),'String','"y" to exit','fontsize',17,'color','k');


set(get(rtH,'title'),'String', []); set(get(tH,'title'),'String',[]);
waitforbuttonpress; kval=double(get(emH,'CurrentCharacter')) ;
        
        switch kval
            case 120 % "x" Remove or merge a cluster and exit
                set(get(midH,'title'),'String',[]);set(get(rtH,'title'),'String', []); set(get(tH,'title'),'String', []);
                set(get(midH,'title'),'String','Which cluster to toss?','fontsize',18,'color',[0.7 0 0.3])
                waitforbuttonpress; kval=double(get(emH,'CurrentCharacter')) ;
                kToss = double(kval-48) ; % convert key into number
                
                if kToss <= nClusters
                inClust(kToss)=[]; % remove cluster from output
                nClusters = nClusters - 1; 
                myColors(kToss) = [];
                end
                
                
            case 109 % "m" merge 
                set(get(midH,'title'),'String',[]);
                set(get(midH,'title'),'String','Type # of First cluster of Merge','fontsize',18,'color',[0.7 0 0.3])
                waitforbuttonpress; kval=double(get(emH,'CurrentCharacter')) ;
                mOne = double(kval-48) ;
                
                set(get(midH,'title'),'String','Type # of Second cluster of Merge','fontsize',18,'color',[0.3 0 0.7])
                waitforbuttonpress; kval=double(get(emH,'CurrentCharacter')) ;
                mTwo = double(kval-48);
                
                inClust{mOne} = ([inClust{mOne}' inClust{mTwo}'])' ;inClust(mTwo) = [] ;
                nClusters = nClusters - 1;
                myColors(mTwo) = []; 
        end % ENd kvalswitch [toss or merge]
        
        
end % END "decision" - switch for [keep,try again,manual,toss/merge]

end % END "tossMerge - switch 

end % END EM cluster loop

     movegui(emH,'northwest')   
end % END EMclusterGui function
