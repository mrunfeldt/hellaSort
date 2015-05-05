
% [duration,maxAmp,timeToPeak, preMin, postMin ] = waveForm_descript(wave,time)

function [duration,maxAmp,timeToPeak, preMin, postMin,firstMin,secondMin ] = waveForm_descript(wave,time)
%wave = waveForms(20,:);time = T;

pre_halfMS = find(time >= -0.5) ;
post_halfMS = find(time >= 0.65) ;
tZero = find(time >= 0) ; tZero = tZero(1);

if ~isempty(pre_halfMS) & ~isempty(post_halfMS) 
    pre_halfMS = pre_halfMS(1) ; post_halfMS=post_halfMS(1);    
    
else % in case time is shorter than "pre_halfMS" parameter
    pre_halfMS = 1 ; post_halfMS=length(time);
end
    
    first = wave(pre_halfMS:tZero); % Chunk of wave prior to peak
    firstMin= find(first==min(first))+pre_halfMS-1 ; % time at min(samples)
    firstMin = firstMin(end);
    
    second = wave(tZero:post_halfMS); % % Chunk of wave post peak
    secondMin= find(second==min(second))+tZero-1 ; % time at min(samples)
    secondMin = secondMin(1);
    
duration = time(secondMin)-time(firstMin) ; % duration of spike (ms)
preMin = min(first) ;postMin = min(second) ;
peakAmp = wave(tZero) ;
maxAmp = peakAmp - preMin ;
timeToPeak = abs(time(firstMin)) ;


% figure;hold on; plot(time,wave,'k','linewidth',3); 
% plot(time(firstMin),preMin,'b*','markersize',15);
% plot(time(secondMin),postMin,'b*','markersize',15)
% line([time(firstMin) time(secondMin)],[min(first) min(first)],'color','r','linewidth',4)
% pause;close
end
