function [] = L4_timeseries_Epsilon(filePI,fileTEST,fileL3TEST,pi,tester)
%L4_TIMESERIES_EPSILON  plots time series of epsilon from the PI and the 
% tester 
% input:
%   dataset = path and prefix data set
%   pi = suffix PI data set (without .nc ending)
%   tester = suffix tester data set (without .nc ending)
%   tester_name = ID of tester for figure legend and file name

% filePI = [dataset pi '.nc'];
% fileTEST = [dataset tester '.nc'];

% read in time
timePI = ncread(filePI,'/L4_dissipation/TIME');
timeTEST = ncread(fileTEST,'/L4_dissipation/TIME');
timeL3TEST = ncread(fileL3TEST,'/L4_dissipation/TIME');

if (floor(1e4.*(timeTEST(1)))~=floor(1e4.*timePI(1)))
    disp('Warning: Mismatch in starting time of records.')
end


% convert time to minutes 
timeL3TEST = (timeL3TEST - timeL3TEST(1))*24*60*60;
timeTEST = (timeTEST - timeTEST(1))*24*60*60;
timePI = (timePI - timePI(1))*24*60*60;

% read in epsi
epsiPI = ncread(filePI,'/L4_dissipation/EPSI');
epsiTEST = ncread(fileTEST,'/L4_dissipation/EPSI');
epsiL3TEST = ncread(fileL3TEST,'/L4_dissipation/EPSI');

n = size(epsiPI,2); % number shear probes

for ii=1:n
    subplot(n,1,ii)
    hold on
    semilogy(timePI,epsiPI(:,ii),'k.-','Color',[0.8 0.38 0.08],'LineWidth',3,...
        'MarkerSize',20)
    try
    semilogy(timeTEST,epsiTEST(:,ii),'k.-','Color','k','LineWidth',1.5,...
        'MarkerSize',10)
    semilogy(timeL3TEST,epsiL3TEST(:,ii),'k.-','Color','m','LineWidth',1.5,...
        'MarkerSize',10)
    catch
    semilogy(timeTEST,epsiTEST(ii,:),'k.-','Color','k','LineWidth',1.5,...
        'MarkerSize',10)
    semilogy(timeL3TEST,epsiL3TEST(ii,:),'k.-','Color','m','LineWidth',1.5,...
        'MarkerSize',10)
    end
    hold off
    xlim(timePI([1 end]))
    ylim([1e-12 1e-1])
    box on
    grid on
    ylabel([char(949) ' shear ' int2str(ii) ' (W kg^-^1)'])
    if ii==n
       xlabel('time elapsed (seconds)') 
       legend(pi,tester,[tester '_L3'],'Location','best')
    end
    set(gca,"YScale",'log')
end

end

