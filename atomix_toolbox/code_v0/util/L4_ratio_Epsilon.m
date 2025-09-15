function [] = L4_ratio_Epsilon(filePI,fileTEST,fileL3TEST,pi,tester)
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
patchy = [1/sqrt(2)*[1 1] sqrt(2)*[1 1] (1/sqrt(2))];

figure('rend','painters','pos',[10 10 1000 500])

for ii=1:n
    subplot(1,n,ii)
    hold on
    patch(patchx,patchy,.7*[1 1 1 ],'LineStyle','none')
    line(axlim,[1 1],'Color','k','LineWidth',1.5)
    p1=plot(epsiPI(:,ii),epsiTEST(:,ii)./epsiPI(:,ii),'k.','MarkerSize',16);
    p2=plot(epsiPI(:,ii),epsiL3TEST(:,ii)./epsiPI(:,ii),'m.','MarkerSize',16);

    xlim(axlim)
    hold off
    set(gca,'xscale','log')
    box on
    grid on
    ylim([0 3])
    legend([p1,p2],'L2','L3')
    xlabel([char(949) '' pi ' (W kg^-^1)'])
    ylabel([char(949) ' ' tester ':' char(949) ' PI'])
    text(0,1,['shear ' int2str(ii)],'units','normalized','VerticalAlignment','bottom')
    set(gca,'layer','top')
end

end

