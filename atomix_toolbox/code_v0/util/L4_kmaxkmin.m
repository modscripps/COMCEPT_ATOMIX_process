function [] = L4_kmaxkmin(filePI,fileTEST,pi,tester)
%L4_EPSILON_SCATTER plots of epsilon from the PI and the tester 
% input:
%   dataset = path and prefix data set
%   pi = suffix PI data set (without .nc ending)
%   tester = suffix tester data set (without .nc ending)
%   tester_name = ID of tester for figure legend and file name

% filePI = [dataset pi '.nc'];
% fileTEST = [dataset tester '.nc'];

% read in epsi
% read in time
timePI = ncread(filePI,'/L4_dissipation/TIME');
timeTEST = ncread(fileTEST,'/L4_dissipation/TIME');

if (timeTEST(1)~=timePI(1))
    disp('Warning: Mismatch in starting time of records.')
end


% convert time to minutes 
timeTEST = (timeTEST - timeTEST(1))*24*60*60;
timePI = (timePI - timePI(1))*24*60*60;

% read in epsi
try
    kmaxPI = ncread(filePI,'/L4_dissipation/KMAX');
    novar=0;
catch
    novar=1;
    kmaxPI=timePI.*0;
end
try
    kminPI = ncread(filePI,'/L4_dissipation/KMIN');
catch
    kminPI=timePI.*0;
end


if novar==0

    kmaxTEST = ncread(fileTEST,'/L4_dissipation/KMAX');
    try
    kminTEST = ncread(fileTEST,'/L4_dissipation/KMIN');
    catch
    kminTEST = timeTEST.*0;;
    end

    n = size(kmaxPI,2); % number shear probes

    figure('rend','painters','pos',[10 10 1000 500])
    for k=1:2 % kmax kmin
        if k==1
            wh_kPI=kmaxPI;
            wh_kTEST=kmaxTEST;
            markerstyle='s-'
        else
            wh_kPI=kminPI;
            wh_kTEST=kminTEST;
            markerstyle='p-'
        end
    for ii=1:n

        subplot(n,1,ii)
        hold on
        plot(timePI,wh_kPI(:,ii),markerstyle,'Color',[0.8 0.38 0.08],'LineWidth',1,...
            'MarkerSize',5)
        try
            plot(timeTEST,wh_kTEST(:,ii),markerstyle,'Color','k','LineWidth',.5,...
                'MarkerSize',5)
        catch
            plot(timeTEST,wh_kTEST(ii,:),markerstyle,'Color','k','LineWidth',.5,...
                'MarkerSize',5)
        end
        if(ii==1 && k==1)
            title('Kmax','FontSize',15,'FontName','times new roman')
        end
        if(ii==1 && k==2)
            title('Kmax - Kmin','FontSize',15,'FontName','times new roman')
        end
        hold off
        xlim(timePI([1 end]))
        box on
        grid on
        ylabel([char(949) ' K ' int2str(ii) ' (VAR_RESOLVED\)'])
        if ii==n
            xlabel('time elapsed (seconds)')
            legend(pi,tester,'Location','best')
        end
    end
    end
end

end

