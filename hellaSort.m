
% % % % % % HELLA-SORT : SPIKESORTING GUI of AWESOMENESS % % % % %
% % % Created by MJRunfeldt, March 2015 % % % 

% INPUT: "nev": struct  of neural events and spiketimes
% () nev.events - [N x T], event waveforms (in sampled time)
% () nev.spiketimes - [ 1 x N ], spiketimes in samples
% () nev.spiketimesMS - [ 1 x N ], spiketimes in ms
% () nev.stepBack - [ 1 ], used to align spikes at t=0 in plots 
% OUTPUT: one cell per neuron/cluster
% (1) "IDs" - indicies of spikes in a cluster, referring back to original input
% % % indicies. e.g. IDs{1}(1:5) - original spikes 1-5 correspond to neuron 1
% (2) "spike_samp" - spiketimes at original sampling rate
% (3) "spike_ms" - spiketimes in ms
% (4) "forms" - waveforms


% % % % % Functions Called % % % % %   
% #EMclusterGUI #gausEMcluster #gaussian_overlap #hellaSort_mergeOption
% #pcasvd #parseparam #plotTwoWaves #shadedLines #summaryFigure
% #waveForm_descript


function [IDs, spike_samp, spike_ms, forms ] = hellaSort(nev)
close all

% % % % INPUT FORMAT % % % % % % % % % % % % % % % % % % % %
samprate = 24414.0625; % assumes the data is from TDT/brainware
events = nev.events; % spike x timePoint
spiketimes = nev.spikeTimesMS ;
spiketimesBins = nev.spikeTimes; 
stepBack = nev.stepBack ;% samples (not ms): Parameters from event extraction
% % % % % % % % % % % % % % % % % % % % 
nEvents = size(events,1); % total # of spikes
eventDur = size(events,2); % duration of events (in bins) 


% % %  VISUALIZATION PARAMETERS % % % % 
plotN = 1:20 ; % # of waveforms plotted during initial exploration of pc space
pickN = 25 ; % how many waveforms to plot after initial selection
fractionClust = 20 ;% what fraction of nodes to visualize (manual cluster explore)

edges = linspace(0,max(spiketimes)*1e-3,40); % for plotting binned rates
interBin = diff(edges); interBin=interBin(1);% for plotting binned rates
isiBins= logspace(-2,4,80);
T = (-stepBack:eventDur-stepBack-1)/samprate*1000; %  provide time axis for waveform plotting


% PCA
[~,uu] = princomp(events) ;

% % % % % % % % % % % % WHEN "terminate ==1", SORT IS DONE % % % % % %
terminate = 0 ;keepExploring = 1;
while terminate == 0

% % % % % % % % PLOT SPIKES AS POINTS/NODES IN PCA SPACE % % % % % % % % % 
% % CLICK AROUND TO CHECK OUT WAVEFORMS AND MARK SOEM DOTS % % % % %
% keepExploring = 1;
while keepExploring == 1
    
hFig = figure ;xwidth = 12e2; ywidth=10e2;
set(hFig, 'Position', [0 0 xwidth ywidth]);movegui(hFig,'center');
pcH=subplot(2,4,[1 2]);set(gca,'Color',[0.7 0.7 0.7])
title('Click around','fontsize',16,'color','r');

% set axis limits
if min(uu(:,1))<0; xxAxis=([1.2*min(uu(:,1)) 1.2*max(uu(:,1))]) ;
elseif min(uu(:,1))>0;xxAxis = ([min(uu(:,1))-.2*min(uu(:,1))   1.2*max(uu(:,1))]);end
if min(uu(:,2))<0;yyAxis=([1.2*min(uu(:,2)) 1.2*max(uu(:,2))]) ;
elseif min(uu(:,1))<0;yyAxis=([min(uu(:,2))-0.2.*min(uu(:,2)) 1.2*max(uu(:,2))]);end
xborder = max(abs(xxAxis));yborder = max(abs(yyAxis));

% % % % EXPLORE PC SPACE AND SPIKE WAVEFORMS, PLOTS SOME DOTS % % % 
mousepoint(1) = xborder/2;mousepoint(2) = yborder/2;
ptCount=0; val = 0; nodes = []; nodes2=[]; starz=[];starz2=[];
while val ~=113 % Until you hit "q"
    set(gca,'xlim',xxAxis,'ylim', yyAxis); % maintain axis
    
subplot(2,4,[1 2]);hold on;plot(uu(:,1),uu(:,2),'y.','markersize',1); % plot all pts
plot(uu(nodes,1),uu(nodes,2),'r.','markersize',5); % red dots
plot(uu(nodes2,1),uu(nodes2,2),'c.','markersize',5); % blue dots
set(gca,'Color',[0.7 0.7 0.7])

[mousepoint]=ginput(1) ;ptCount =ptCount+1 ; % CROSSHAIR/arrow 
% % % % Sort indicies (spikes) according to nearest to crosshair % % % 
xDiff = abs(uu(:,1)-mousepoint(1)) ; yDiff = abs(uu(:,2)-mousepoint(2)) ; 
tDiff=sqrt(xDiff.^2+yDiff.^2)  ; [~,nearest] = sort(tDiff,'ascend') ; % sort inidicies within "nodes"
plot(uu(nearest(plotN),1),uu(nearest(plotN),2),'m.','markersize',2);drawnow

% Plot waveforms for "plotN" number of nodes nearestnea to mousepoint
wvH=subplot(2,4,[3 4],'replace');hold on;line([0 0],[min(min(events)) max(max(events))],'color','k');
set(gca,'xlim',[T(1) T(end)]); plot(T,events(nearest(plotN)',:),'linewidth',2);

set(get(pcH,'title'),'String','Press "1" to mark red, "2" blue, "4" to skip, "x" if there are no real spikes','fontsize',17,'color',[0.1 0.8 0.2]);
xlabel('time (sec)');title('Hit "q" to move on.','fontsize',18,'color','k');hold off

% % PRESS "1" TO MAKE RED DOTS, "2". PRESS ANY OTHER KEY TO SKIP. IF NO KEY IS
waitforbuttonpress;val=double(get(hFig,'CurrentCharacter')) ;
switch val
    case 49  % if "1" is pressed
    starz = [starz; mousepoint(1), mousepoint(2)] ; % location of all those starz (to use if polygon needs to be cleared
    subplot(2,4,[1 2]);hold on;plot(mousepoint(1),mousepoint(2),'x','color','r','markersize',3,'linewidth',2);hold off
    nodes=[nodes nearest(plotN)'];
    case 50       % if 2 is pressed
    starz2 = [starz2; mousepoint(1), mousepoint(2)] ; % location of all those starz (to use if polygon needs to be cleared
    subplot(2,4,[1 2]);hold on;plot(mousepoint(1),mousepoint(2),'x','color','c','markersize',3,'linewidth',2);hold off
    nodes2=[nodes2 nearest(plotN)'];
    case 120 % EXIT entire sorter
        IDs = []; timesBins=[]; forms = []; timesMS =[]; close; return
        
end
end % END ["while val ~=q"] click-around and make dots loop



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % SELECT POLYGONAL AREA OF SPIKES THAT YOU WANT TO KEEP % % % 
val=0;
while val ~=13 % Until you hit RETURN

    % % % SELECT ALL OR CLICK POLYGON % % % 
set( get(wvH,'title'), 'String','Type "a" to select all OR CLICK To make polygon')
waitforbuttonpress;pcv=double(get(hFig,'CurrentCharacter'))  ;

if pcv == 97 % "a" - select all
    keepers = 1:nEvents ; val = 13 ; 
else
    
% % % % % DRAW A POLYGON AROUND ALL SPIKES YOU WANT TO KEEP % % % %     
set( get(pcH,'title'), 'String', 'Click polygon outline and then hit RETURN','fontsize',20,'color',[0.1 0.8 0.2]);drawnow
       
% % DRAW POLYGON % % 
subplot(2,4,[1 2]);hold on;
startNode = ginput(1); % beginnning of poly (for visualization)
plot(startNode(1) ,startNode(2),'x','linewidth',2,'markersize',5,'Color',[0.5 0 0.5]);drawnow
polyDum = ginput ;polygon = [startNode; polyDum; startNode]; 
plot(polygon(:,1) ,polygon(:,2),'Color',[0.5 0 0.5]);drawnow % % press "return

keepers = find(inpolygon(uu(:,1),uu(:,2),polygon(:,1),polygon(:,2))); % indicies within selected area
plot(uu(keepers,1),uu(keepers,2),'.','markersize',1,'Color',[1 0 0.5]);
        
% % %"keepers" ==> indices of spikes that are kept (i.e. within polygon) % % % 

% % % % % % % % % Plot MEAN+/- STD of waveform % % % % % % % % % % % 
group1 = keepers; group2=setdiff(1:nEvents,keepers);
sh=subplot(2,4,[7 8],'replace'); waveDum =  events(keepers,:);
waveMean = mean(waveDum) ;  waveSTD = std(waveDum) ;%./ length(keepers) ; 
shadedLines(T,waveMean, waveSTD,'r')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % % Plot Binned Rate % % % % % % % % % % % 
subplot(2,4,[5 6],'replace');hold on; 
plot(edges,histc(spiketimes(group1).*1e-3,edges)./interBin,'m','linewidth',2);
plot(edges,histc(spiketimes(group2).*1e-3,edges)./interBin,'b','linewidth',2); 
leg=legend('Selected','All Others');set(leg,'fontsize',14) ;
xlabel('Time (sec)');title('Binned Spike Rate');ylabel('Rate (Hz)');
set(gca,'xlim',[0 edges(end-1)]);
 
% % % % PLOT SUBSET OF SELECTED WAVEFORMS % % % % 
if length(keepers) < pickN;pickN = length(keepers);end
whyNotUs = randperm(length(keepers),pickN) ; % random selection of spikes
subplot(2,4,[1 2]),title([]);subplot(2,4,[3 4],'replace');wvH=subplot(2,4,[3 4]);plot(T,events(whyNotUs,:)');
set(gca,'xlim',[T(1) T(end)],'ytick',[]);xlabel('time (sec)')
title('Hit RETURN if polygon is a keeper, "q" to try again','fontsize',14)
waitforbuttonpress; val=double(get(hFig,'CurrentCharacter')) ;

if val== 113 % if "q" is pressed, clear polygon. Re-plot nodes and starz
% % % clear plot: redraw nodes and starz keep red and blue dot % % %
    subplot(2,4,[1 2],'replace');pcH=subplot(2,4,[1 2]);hold on
    set(gca,'Color',[0.7 0.7 0.7]);set(gca,'xlim',xxAxis,'ylim', yyAxis);
    hold on;plot(uu(:,1),uu(:,2),'y.','markersize',1); 
    if ~isempty(starz);plot(starz(:,1),starz(:,2),'x','color','r','markersize',3,'linewidth',3);end
    if ~isempty(starz2);plot(starz2(:,1),starz2(:,2),'x','color','c','markersize',3,'linewidth',3);end
    hold off
end
end % END if select all OR click polygon
end % END ["while val ~= RETURN"] draw-polygon


% % % % % % % % Plot ISI histogram % % % % % % % % % % % % % %
isi = diff(spiketimes(keepers));isi1 = 100 * sum(isi<=1)./length(isi);
bins=histc(isi,isiBins) ; 
subplot(2,4,[5 6],'replace');plot(log10(isiBins),bins,'k','linewidth',2);xlim([log10(isiBins(1)) log10(isiBins(end-1))])
binPlot = isiBins(2:10:end) ; xlabel('ISI (ms) - logscale');ylabel('Count')
set(gca,'xtick',log10(binPlot),'xticklabel',num2cell(round(binPlot.*10)./10))%'XScale','log')
msMark = isiBins(isiBins >1); msMark = log10(msMark(1));
line([msMark msMark],[0 max(bins)],'linewidth',3)
title(['Percent ISI less than 1ms : ',num2str(round(isi1.*10)./10)],'fontsize',14)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


set( get(pcH,'title'), 'String', 'Hit "1" if single neuron. Hit "k" to use EM on selected,','fontsize',16,'color','b');
set( get(wvH,'title'), 'String', 'Hit "n" to manually cluster','fontsize',16,'color','b');drawnow

waitforbuttonpress; Kval=double(get(hFig,'CurrentCharacter')) ;
switch Kval
    case 49 % "1" only one neuron, plot SUMARY FIGURE and allow to accept, try again, or move on % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % EXIT SUMMARY % % % % % % % % % % % % % % % %       
close; nrns=[]; nrns{1} = keepers'; exitVal = 121;
        switch exitVal
            case 110 %"n" move on to manual cluster
                keepExploring = 0 ; nrns=[]; nextMod=1;
            case 121 % "y" save and exit
                 keepExploring = 0 ; terminate = 1; tryAgain = 0 ;nextMod=0;
            case 116 % "t" restart loop (i.e. make dots and draw polygon)
                close; keepExploring = 1; nrns=[];
        end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    case 110 % "n" Move on to manually cluster
        keepExploring = 0 ; nrns=[]; nextMod = 1; 
        
% % % % % % GAUSSIAN MIXTURE MODEL CLUSTERNG % % % %         
    case 107 %"k" USE EM (RESTRICT TO ONLY SELECTED SPIKES) % % %        
        set( get(pcH,'title'), 'String', 'How Many Clusters?','fontsize',16,'color',[0.7 0.2 0.3]);
        set( get(wvH,'title'), 'String', []);drawnow
        waitforbuttonpress; clustVal=double(get(hFig,'CurrentCharacter')) ;       

        [gmClusters,decision] = EMclusterGUI(clustVal,events(keepers,:),spiketimes(keepers),T) ; % Gaus Mix Model

        switch decision
            case 121 % "y" save and exit
                nrns=[];
                for i = 1:length(gmClusters); 
                    nrns{i} = keepers(gmClusters{i});end % CONVERT BACK TO ORIGINAL INDICIES
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % EXIT SUMMARY % % % % % % % % % % % % % % % % 
             [exitVal,nrns,confusion] = summaryFigure(events,spiketimes,T,nrns) ; 
                switch exitVal; 
                    case 110; keepExploring = 0 ; nrns=[]; nextMod=1;%"n" move on to manual cluster
                    case 121; keepExploring = 0 ; terminate = 1; nextMod=0;% "y" save and exit
                    case 116; close; keepExploring = 1; nrns=[];% "t" restart loop (i.e. make dots and draw polygon)
                end 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            case 109  %"m" move on to manual cluster
                keepExploring = 0 ; nrns=[]; nextMod=1;
            case 116  %"t" try again from the top
                close; keepExploring = 1 ; nrns=[];
                
        end % END switch - EM/Gauss mix loop
end % END switch - decision post polygon cluster
        
end % END while-loop ["while keepExploring == 1", select subsets manually]
    



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % MANUAL CLUSTERING ON SELECTED SPIKES IN NEW PCA SPACE% % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if nextMod == 1;
newInd = keepers; newWaves = events(keepers,:);
newSpiketimes = spiketimes(keepers); %% Restrict to only the kept spikes. VIP keep track of your indicies!!!!!
nClust = round((1/fractionClust)*length(keepers)); % determines what fraction of nodes to select at a time

% % % NEW FIGURE %  % % PLOT "kept" spikes in new PC space % % % % % % 
[~,u] = princomp(newWaves) ;
myColors = [0.7 0 0.7; 0 0.8 0.5; 0 0.3 0.8; 1 0 0];
end


% % % % % % % % % %  % % % % NEW FIGURE % % % % % % % % % % % % % % 
% % % %  % % % % MANUAL/POLYGON CLUSTERING % % % %  % % % % 
if nextMod == 1; keepItUp = 1; else keepItUp=0;end
while keepItUp == 1 
    
mkFig = figure; xwidth = 14e2; ywidth=6e2;set(mkFig, 'Position', [0 0 xwidth ywidth]) % main figure
movegui(mkFig,'center');
    mPCh=subplot(2,3,[2,5]);hold on; title('How Many Clusters (max == 4) ?','fontsize',16); 
    plot(u(:,1),u(:,2),'k.','markersize',1);set(gca,'Color',[0.8 0.8 0.8]) % PC space
    if min(u(:,1))<0; xxAxis=([1.2*min(u(:,1)) 1.2*max(u(:,1))]) ;elseif min(u(:,1))>0;xxAxis = ([min(u(:,1))-.2*min(u(:,1))   1.2*max(u(:,1))]);end
    if min(u(:,2))<0;yyAxis=([1.2*min(u(:,2)) 1.2*max(u(:,2))]) ;elseif min(u(:,1))<0;yyAxis=([min(u(:,2))-0.2.*min(u(:,2)) 1.2*max(u(:,2))]);end
    xborder = max(abs(xxAxis));yborder = max(abs(yyAxis));
    
    waitforbuttonpress;mClusters=double(get(mkFig,'CurrentCharacter')) ; % HOW MANY CLUSTERS?
    nClusters = double(mClusters-48) ; % convert key into number

% % % DRAW EACH CLUSTER AND WAVEFORMS/RATES/ISI % % %     
myClusters=[]; cnt=0;  
    for ck = 1:nClusters
        set(mPCh,'xlim',xxAxis,'ylim', yyAxis); % maintain axis
        set(get(mPCh,'title'),'String',['Select cluster # ', num2str(ck), '(click poly, then RETURN)'],'fontsize',18,'color',myColors(ck,:))
        
% % %  DRAW POLYGON % % %         
        pcHand=subplot(2,3,[2 5]);hold on; 
        startNode = ginput(1); % beginnning of poly (for visualization)
        plot(startNode(1) ,startNode(2),'x','linewidth',2,'markersize',5,'Color',myColors(ck,:));drawnow
        polyDum = ginput ;polygon = [startNode; polyDum; startNode]; 
        plot(polygon(:,1) ,polygon(:,2),'Color',myColors(ck,:));drawnow % % press "return" when full polygon selected
        thisCluster = find(inpolygon(u(:,1),u(:,2),polygon(:,1),polygon(:,2))); % indicies within selected area
        
        if ~isempty(thisCluster) % if the cluster is not empty
            cnt=cnt+1;myClusters{cnt} = thisCluster;
% % if nodes already belong to another cluster, assign to the most-recently-drawn % % 
        if ck>1; for v=1:cnt-1; [crossOver,traitors] = intersect(myClusters{v},thisCluster); 
                if ~isempty(crossOver); myClusters{v}(traitors)=[];end
        end;end
             
        plot(u(thisCluster,1),u(thisCluster,2),'.','markersize',2,'Color',myColors(ck,:)); % plot cluster pts
% % % WaveForm of Selected cluster  % % % %        
        atH= subplot(2,3,1,'replace');
         waveMean = mean(newWaves(thisCluster,:)) ; waveSTD = std(newWaves(thisCluster,:));
         shadedLines(T,waveMean,waveSTD,myColors(ck,:));legend('Selected')

% % % Selected Cluster vs Other clusters % % %        
        subplot(2,3,3);hold on;title('All Clusters','fontsize',15)% PLOT WAVEFORMS
        waveMean = mean(newWaves(thisCluster,:)) ; waveSTD = std(newWaves(thisCluster,:));
        shadedLines(T,waveMean,waveSTD,myColors(ck,:));legend('Selected')
        %set(gca,'xlim',[time(1) time(end)],'color',[0.9 0.9 0.9]);xlabel('Time (sec)');
        
        subplot(2,3,6);hold on; % PLOT BINNED RATES
        plot(edges,histc(newSpiketimes(thisCluster).*1e-3,edges)./interBin,'color',myColors(ck,:),'linewidth',3);
        xlabel('Time (sec)');title('Binned Spike Rate');ylabel('Rate (Hz)');
        set(gca,'xlim',[0 edges(end-1)]);
               
% % % % % % % % Plot ISI histogram % % % % % % % % % % % % % %
        isi = diff(newSpiketimes(thisCluster));isiBins= logspace(-2,3.5,70); bins=histc(isi,isiBins) ; 
        subplot(2,3,4);hold on;plot(log(isiBins),bins,'color',myColors(ck,:),'linewidth',2);xlim([log(isiBins(1)) log(isiBins(end))])
        binPlot = isiBins(2:10:end) ; xlabel('ISI (ms) - logscale');ylabel('Count')
        set(gca,'xtick',log(binPlot),'xticklabel',num2cell(round(binPlot.*10)./10))%'XScale','log')
        line([log(1) log(1)],[0 max(bins)],'linewidth',3);title('InterSpikeIntervals')
        end % END if statement for cluster comtaining at least one spike
    end % END (ck) per manual cluster
    
set(get(mPCh,'title'),'String',[]); set(get(atH,'title'),'String',[]);
subplot(2,3,1);title('Type "y" to keep, "n" to try again, "t" to start from the top','fontsize',16);
waitforbuttonpress; anotherVal=double(get(mkFig,'CurrentCharacter'));subplot(2,3,1);title([])

switch anotherVal
    case 121 % "y" save and exit
        inItsRightPlace = unique(cell2mat(myClusters'))';
        if length(inItsRightPlace) ~= length(newSpiketimes) ;  % if all nodes in NOT in unique group
         subplot(2,3,3);title([]);subplot(2,3,[2 5]);hold on;
            title('Assign stragglers to nearest node? "y" = yes, "n" = toss','color','k','fontsize',18);
            stragglers = setdiff(1:length(newSpiketimes),inItsRightPlace) ; % IDs of stragglers
            waitforbuttonpress; dah=double(get(mkFig,'CurrentCharacter')) ; cAss = double(dah-48) ; % convert key into number
            switch dah
                case  121 % "y" assign stragglers to nearest node center
                    tDiff=[];for kk=1:ck; ccenter = [mean(u(myClusters{kk},1)), mean(u(myClusters{kk},2))];
                        subplot(2,3,[2,5]);hold on; plot(ccenter(1),ccenter(2),'x','color',myColors(kk,:),'markersize',20,'linewidth',3);
                        stragCent = u(stragglers,:) ; % x,y
                        stragDistX =  abs(stragCent(:,1)-ccenter(1)) ; stragDistY =  abs(stragCent(:,2)-ccenter(2)) ;
                        tDiff(kk,:) = sqrt(stragDistX.^2+stragDistY.^2) ;end
                    [~,nearest]=min(tDiff) ; clear tDiff% inidicies correspond to clusters;
                    for kk=1:ck; myClusters{kk} = [myClusters{kk}' stragglers(nearest==kk)];
                        plot(u(myClusters{kk},1),u(myClusters{kk},2),'.','markersize',2,'Color',myColors(kk,:));
                        nrns{kk} = keepers(myClusters{kk}) ;
                    end % assign to cluster

                case 110 % "n" % ignore unassigned nodes
                    for kk=1:nClusters; nrns{kk} = keepers(myClusters{kk}) ;end
            end % END switch - what to do with stragglers (unassigned spikes)
        end % END if stragglers exist
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % EXIT SUMMARY % % % % % % % % % % % % % % % % 
    [exitVal,nrns,confusion] = summaryFigure(events,spiketimes,T,nrns) ; 
        switch exitVal; 
            case 121; close; keepItUp = 0 ; terminate = 1; % "y" save and exit
            case 110; keepItUp = 0 ; nrns=[]; keepExploring=0;
                movegui(mkFig,'northwest');terminate = 0;%"n" move on to manual cluster
            case 116; keepItUp = 1; nrns=[];% "t" restart loop (i.e. make dots and draw polygon)
                set(get(mPCh,'title'),'String','Restart from very beginning "y", or just the clustering "n" ?','fontsize',16);
                waitforbuttonpress; doom=double(get(mkFig,'CurrentCharacter'));
                switch doom
                    case 110 % "n" try again - re-draw clusters (new fig)
                        oldfig = mkFig; movegui(oldfig,'northwest'); % move old GUI to northwest             
                        keepItUp = 1; nrns = [];
                    case 121 % "y" restart from the very beginning
                        close all; keepItUp = 0 ; terminate = 0; keepExploring = 1; keepers = []; nrns=[];
                end
        end 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %              
    
    case 110 % "n" try again - re-draw clusters (new fig)
        keepItUp = 1; nrns = []; oldfig = mkFig; movegui(oldfig,'northwest'); % move old GUI to northwest    
        
    case 116 % "t" clear it all and start from the top
        keepItUp = 0; nrns = []; keepers = []; terminate = 0 ;  
        
        set(get(atH,'title'),'String','Restart from very beginning "y", or just the clustering "n" ?','fontsize',16);
        waitforbuttonpress; doom=double(get(mkFig,'CurrentCharacter'));
        
        switch doom
            case 110 % "n" try again - re-draw clusters (new fig)
                oldfig = mkFig; movegui(oldfig,'northwest'); % move old GUI to northwest             
                keepItUp = 1; nrns = [];
            case 121 % "y" restart from the very beginning
                close all; keepItUp = 0 ; terminate = 0; keepExploring = 1; keepers = []; nrns=[];
        end
        
end % END switch - what to do after assigning clusters with polygon

end % END ["while keepItUp==1] - manual/polygon clustering

end % END master loop ["while terminate == 0]

% % % % Rename variables for output. Extract SpikeTimes from sorted neurons % % % 
IDs = nrns; spike_samp =[]; spike_ms=[]; forms =[] ;
for z = 1:length(nrns)
   spike_samp{z} = spiketimesBins(nrns{z}) ; forms{z} = events(nrns{z},:) ; 
   spike_ms{z} = spiketimes(nrns{z}) ;
end



close all
figure;spy;set(gca,'ytick',[],'xtick',[],'xlabel',[]);title('Good job, human. Hit any key to continue','fontsize',20)
waitforbuttonpress; close 

end % END hellSort function



