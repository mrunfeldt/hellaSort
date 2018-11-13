
% % % Draws middle line and +/- value, transparent patch
% MJRunfeldt March 2015

function [] = shadedLines(T,mid,wing,mainColor)

bottom = mid-wing; top = mid+wing;

x= [1:length(top), fliplr(1:length(top))];  % samples
x = T(x);

y = [top, fliplr(bottom)] ;

hold on; set(gca,'color',[0.9 0.9 0.9])
patch(x,y,mainColor,'facealpha',0.3,'edgecolor','none')
plot(T,mid,'color',mainColor,'linewidth',4); xlabel('Time (ms)');
xlim([T(1) T(end)])

end


%mid = waveMean; wing = waveSTD ;