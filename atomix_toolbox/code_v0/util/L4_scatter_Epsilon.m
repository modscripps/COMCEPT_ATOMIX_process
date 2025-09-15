function [] = L4_scatter_Epsilon(filePI,fileTEST,fileL3TEST,pi,tester)
%L4_EPSILON_SCATTER plots of epsilon from the PI and the tester 
% input:
%   dataset = path and prefix data set
%   pi = suffix PI data set (without .nc ending)
%   tester = suffix tester data set (without .nc ending)
%   tester_name = ID of tester for figure legend and file name

% filePI = [dataset pi '.nc'];
% fileTEST = [dataset tester '.nc'];

% read in epsi
epsiPI = ncread(filePI,'/L4_dissipation/EPSI');
epsiTEST = ncread(fileTEST,'/L4_dissipation/EPSI');
epsiL3TEST = ncread(fileL3TEST,'/L4_dissipation/EPSI');

n = size(epsiPI,2); % number shear probes
if size(epsiTEST,1)~=size(epsiPI,1)
    epsiTEST=epsiTEST.';
    epsiL3TEST=epsiL3TEST.';
end

if size(epsiTEST,1)~=size(epsiPI,1)
   disp('Error: Number of epsi records does not match!')
   disp('Stop plotting.')
   return
end

mineps=min(min([epsiPI epsiTEST]));
maxeps=max(max([epsiPI epsiTEST]));

axlim = [ ((mineps)) ((maxeps)) ];

patchx = [axlim fliplr(axlim) axlim(1)];
patchy = [axlim*(1/sqrt(2)) fliplr(axlim)*sqrt(2) axlim(1)*(1/sqrt(2))];

figure('rend','painters','pos',[10 10 1000 500])

for ii=1:n
    subplot(1,n,ii)
    hold on
    patch(patchx,patchy,.7*[1 1 1 ],'LineStyle','none')
    line(axlim,axlim,'Color','k')
    %line(axlim,axlim*sqrt(2),'Color','k','LineStyle','--')
    %line(axlim,axlim*(1/sqrt(2)),'Color','k','LineStyle','--')
    p1=scatter(epsiPI(:,ii),epsiTEST(:,ii),[],'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    p2=scatter(epsiPI(:,ii),epsiL3TEST(:,ii),[],'filled',...
        'MarkerEdgeColor','m','MarkerFaceColor','m');
    ylim(axlim)
    xlim(axlim)
    axis equal
    hold off
    set(gca,'xscale','log','yscale','log')
    legend([p1,p2],"L2","L3")
    box on
    grid on
    xlabel([char(949) ' ' pi ' (W kg^-^1)'])
    ylabel([char(949) ' ' tester ' (W kg^-^1)'])
    text(0,1,['shear ' int2str(ii) ' (W kg^-^1)'],'units','normalized','VerticalAlignment','bottom')
    set(gca,'layer','top')
end

end

