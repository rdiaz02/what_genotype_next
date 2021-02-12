clc
close all
clearvars

data = table2array(readtable('data_for_m.txt'));
col = table2array(readtable('colors_for_m.txt'));

fig = figure;

for i = 1:max(data(:,1))
    plot3(data(data(:,1)==i,2),...
          data(data(:,1)==i,3),...
          data(data(:,1)==i,4),...
          'o',...
          'MarkerFaceColor',col(i,:)./255,...
          'MarkerEdgeColor','none',...
          'MarkerSize',4)
      hold all
end

% aesthetics
axis square

d = 0.05;
set(gca,'XLim',[-d 1+d],'YLim',[-d 1+d],'ZLim',[-d 1+d])

set(gca,'XTick',[-d 0 0.25 0.5 0.75 1 1+d],...
    'XTickLabel',{'' '0' '' '0.5' '' '1' ''})
set(gca,'YTick',[-d 0 0.25 0.5 0.75 1 1+d],...
    'YTickLabel',{'' '0' '' '0.5' '' '1' ''})
set(gca,'ZTick',[-d 0 0.25 0.5 0.75 1 1+d],...
    'ZTickLabel',{'' '0' '' '0.5' '' '1' ''})
set(gca,'TickLength',[0 0.025])

set(gca,'XDir','reverse','YDir','reverse')

grid on

view(-15,20)

%{
xlabel({'Distance between';'CE predictions,';'{\it JS} (CBN, MHN)'},...
    'interpreter','tex')
ylabel({'Distance between';'TD predictions,';'{\it JS} (CBN\_td, MHN\_td)'},...
    'interpreter','tex')
zlabel({'Distance between';'CE and TD predictions,';'{\it JS} (MHN, MHN\_td)'},...
    'interpreter','tex')
%}