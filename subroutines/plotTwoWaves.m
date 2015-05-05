
function [] = plotTwoWaves(waveFormz,group1,group2,time)

waveKeeps = waveFormz(group1,:); otherz = waveFormz(group2,:);
waveMean = mean(waveKeeps) ; waveSTD = std(waveKeeps); % mean and STD of waveform
if length(otherz)<1;otherMean= zeros(size(waveMean));otherSTD=zeros(size(waveSTD));
else otherMean = mean(otherz) ; otherSTD = std(otherz); end% mean and STD of waveform
boundedline(time,waveMean,waveMean+waveSTD,'-xm',time,otherMean,otherMean+waveSTD,'-ob','alpha');
set(gca,'xlim',[time(1) time(end)],'Color',[0.9 0.9 0.9]);xlabel('Time (sec)');%title('Mean+STD waveform')

end