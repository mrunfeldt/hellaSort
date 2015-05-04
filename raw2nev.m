
% Input: (1) single vector of raw voltage trace (2) threshold
% event is when : rawVoltage >  median + (threshold *std(rawVoltage))

% Output: (1) nev.waveform (nEvents x time (samples), 
% (2) nev.spikeTimes (time in samples)
% (3) nev.spikeTimesMS (times in ms) 
% created MJRunfeldt 2015_02_29


function [nev] = raw2nev(input,threshold,samprate)
% input = raw; threshold = 3.5 ;

% % % HARDCODED PARAMETERS % % % % % % % % % % % % % %  % % % %
minPeakdurMS = 0.08 ; % raw Voltage must be above threshold for this many ms

% % % Params for peak search % % %
stepBackSearch = round(0.5 *1e-3 *samprate) ; % (ms) step Back from peak [for waveform window]
stepForSearch = round(0.5 *1e-3 *samprate) ; % (ms) step Forward from peak

% % % Params for window extraction (larger than for peak search)
stepBackMS = 0.8 ; % (ms) step Back from peak [for waveform window]
stepForMS = 1.2 ; % (ms) step Forward from peak

% % % % % % % % % % % % % % % % % % % % % % % % % %  % % % %
minPeakdur = ceil(minPeakdurMS *1e-3 *samprate) ;
stepBack = round(stepBackMS *1e-3 *samprate) ; % (samples) step Back from peak
stepFor = round(stepForMS *1e-3 *samprate) ; % (samples) step Back from peak
nev.stepBack = stepBack; nev.stepFor=stepFor; % output parameters (in samples)
nev.stepBackMS = stepBackMS; nev.stepForMS=stepForMS; % output parameters (in samples)

% % % pad onset (usually has transient) % % %
input(1:5) = median(input);

% % % % IDENTIFY TIME POINTS WHEN VOLT RISES ABOVE THRESHOLD AND IS
% MAINTAINED FOR "minPeakdur" UNTIL DROPPING BELOW THRESHOLD % % %
threshPos = median(input) + threshold * std(input); % voltage cross threshold
peakDT = diff((input>threshPos)); peakDT = [peakDT 0] ; % add extra value to make same length as "peakBin"
ups = find(peakDT>0) ;%indicies of positive slope
downs = find(peakDT<0)+1; % indicies of neg slope

%figure;hold on; plot(input); plot(ups,threshPos,'rx');plot(downs,threshPos,'bx')

% % correct if DT vectors are not the same length
if length(ups) > length(downs); downs = [downs ups(end)+minPeakdur]; 
elseif length(ups) < length(downs) ; downs(1)=[];
end
nevs = ups( (downs - ups) >= minPeakdur )+1 ; % TIMES OF EVENTS

if median(downs-ups) < 1 ; disp('check line33 of raw2nev. Might have shift');end

% % % opt: investigate short ISI situations % % %
%ISIs = diff(nevs) / (samprate*1e-3) ; % TIME BETWEEN PEAKS (in ms)
%isiMin=1;closeCall = find(ISIs<isiMin);

% % % Remove events that are too close to beginning or end of recording % %
ons = nevs-stepBackSearch ; offs = nevs+stepForSearch ;
tooClose = find(ons-stepBack<2 |offs+stepFor>length(input));
if ~isempty(tooClose);nevs(tooClose)=[];ons(tooClose)=[];offs(tooClose)=[];end

spikeTimes=zeros(1,length(nevs));
waveForms = zeros(length(nevs),stepBack+stepFor+1) ; % event x time
for n = 1:length(nevs)
    window = input(ons(n):offs(n)); wPeak = find(window==max(window)) ; wPeak = wPeak(1) ;
    spikeTimes(n)= ons(n)+wPeak-1; % identify peak within window 
    newWindow = input(ons(n) + wPeak-stepBack -1: ons(n) + stepFor + wPeak -1 ) ; % re-align by peak
    waveForms(n,:) = newWindow ; 
    % figure;hold on; plot(window,'linewidth',4);plot(waveForms(n,:),'r','linewidth',2);legend('Original','Re-aligned')
    % line([stepBack stepBack],[min(newWindow) max(newWindow)+1],'color','k');pause;close
end


nev.events = waveForms ;
nev.spikeTimes = spikeTimes ; 
nev.spikeTimesMS = spikeTimes ./ (samprate*1e-3) ; 
nev.shortisidx = 1:length(spikeTimes);

%%%
% subset = 50; 
% for k = 1:20
% whyNotUs = randperm(length(nevs),subset) ; % random selection of spikes
% figure;hold on; plot(waveForms(whyNotUs,:)'); 
% line([stepBack stepBack],[min(min(waveForms(whyNotUs,:))) max(max(waveForms(whyNotUs,:)))+1],'color','k')
% xlabel('Time (samples)')
% pause;close
% end
% 
% for i = 1:subset
%     figure;hold on; plot(waveForms(whyNotUs(i),:)'); 
% line([stepBack stepBack],[min(min(waveForms(whyNotUs,:))) max(max(waveForms(whyNotUs,:)))+1],'color','k')
% xlabel('Time (samples)')
% pause;close
% end
% 
% figure; hold on; xlim([1 1e3]); plot(raw(1:1e3));plot(spikeTimes(1:1e3),max(raw(1:1e3)),'r*')