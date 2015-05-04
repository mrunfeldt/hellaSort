

% % % 

clear all

% % % % % % PARAMETERS (you set these) % % % % % %
dPath = '/Users/melissarunfeldt/GitHub/hellaSort_version_1' ; % path were raw data was saved
samprate = 24414.0625 ; % sampling rate in Hz (samples/second)
nevThresh = 3.5 ; % Threshold of extracting neural events (std above median)

visRaw = 1; % Do you want to plot raw data? (0=no, 1=yes)
chunk = [1 1.5] ; % Plot a subset of raw data to save time (set parameter in SECONDS)
chunkInSamp = round(chunk.*samprate) ; % convert to samples/


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % Step 1: extract neural events from raw voltage trace % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

load([dPath,'/exampleData']) % load data

% % % % % % % If desired, plot a subset of raw data % % % % % % % % %
if visRaw==1;timeTick = linspace(0,chunkInSamp(end),10); 
ff=figure; hold on;title('...extracting neural events from threshold crossings')
plot(rawV (chunkInSamp(1):chunkInSamp(end)) ); xlabel('Time (seconds)'); ylabel('Voltage')
set(gca,'xtick',timeTick,'xticklabel',num2cell( round(((timeTick+chunkInSamp(1))/samprate)*1e2).*1e-2 ) )
xlim([0 chunkInSamp(end)-chunkInSamp(1)]);drawnow
else disp('...extracting neural events from threshold crossings')
end


nev = raw2nev(rawV,nevThresh) ; % EXTRACT neural events


% % % % % % % If desired, indicate times of neural events % % % % % % % % %
if visRaw==1; tCross = median(rawV) + (nevThresh * std(rawV)); % plot threshold as y-intercept
    inWindow = nev.spikeTimes(nev.spikeTimes > chunkInSamp(1) & nev.spikeTimes < chunkInSamp(end)); % ID events inchucnk window
    plot(inWindow-chunkInSamp(1),tCross,'r*','markersize',5);
    title('Neural events extracted from threshold crossings');
    movegui(ff,'northwest');xlim([0 chunkInSamp(end)-chunkInSamp(1)]);drawnow
else
disp('Neural events extracted from threshold crossings') 
end 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % Step 2: sort those spikes - hellaSort % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
[IDs, timesBins, forms, timesMS ] = hellaSort(nev) ;

disp('Solve brain') 
