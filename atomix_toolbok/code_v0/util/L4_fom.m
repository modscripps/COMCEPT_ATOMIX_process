function [] = L4_fom(filePI,fileTEST,fileL3TEST,pi,tester)
%L4_EPSILON_SCATTER plots of epsilon from the PI and the tester 
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

if (timeTEST(1)~=timePI(1))
    disp('Warning: Mismatch in starting time of records.')
end


% convert time to minutes 
timeL3TEST = (timeL3TEST - timeL3TEST(1))*24*60*60;
timeTEST = (timeTEST - timeTEST(1))*24*60*60;
timePI = (timePI - timePI(1))*24*60*60;

% read in epsi
epsiPI = ncread(filePI,'/L4_dissipation/FOM');
epsiTEST = ncread(fileTEST,'/L4_dissipation/FOM');
epsiL3TEST = ncread(fileL3TEST,'/L4_dissipation/FOM');

n = size(epsiPI,2); % number shear probes

figure('rend','painters','pos',[10 10 1000 500])

for ii=1:n
    subplot(n,1,ii)
    hold on
    plot(timePI,epsiPI(:,ii),'k.-','Color',[0.8 0.38 0.08],'LineWidth',3,...
        'MarkerSize',20)
    try
    plot(timeTEST,epsiTEST(:,ii),'k.-','Color','k','LineWidth',1.5,...
        'MarkerSize',10)
    plot(timeL3TEST,epsiL3TEST(:,ii),'m.-','Color','m','LineWidth',1.5,...
        'MarkerSize',10)
    catch
    plot(timeTEST,epsiTEST(ii,:),'k.-','Color','k','LineWidth',1.5,...
        'MarkerSize',10)
    plot(timeL3TEST,epsiL3TEST(ii,:),'m.-','Color','m','LineWidth',1.5,...
        'MarkerSize',10)
    end
    hold off
    xlim(timePI([1 end]))
    if(ii==1)
        title('FOM','FontSize',15,'FontName','times new roman')
    end

    box on
    grid on
    ylabel([char(949) ' shear ' int2str(ii) ' (fom)'])
    if ii==n
       xlabel('time elapsed (seconds)') 
       legend(pi,tester,[tester '_L3'],'Location','best')
    end
end

end
