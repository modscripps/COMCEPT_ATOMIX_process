function [] = L3_compare_spectra_plot(filePI,fileTEST,fileL3TEST,pi,tester,fig_dir,section)
%L3_COMPARE_SPECTRA_PLOT compares PI and tester spectra 
% input:
%   dataset = path and prefix data set
%   pi = suffix PI data set (without .nc ending)
%   tester = suffix tester data set (without .nc ending)
%   tester_name = ID of tester for figure legend and file name
%   section = section in record to plot

% filePI = [dataset pi '.nc'];
% fileTEST = [dataset tester '.nc'];

ncdisp(fileTEST)

% read in eps value and spectra
epsiPI = ncread(filePI,'/L4_dissipation/EPSI',[section 1],[1 Inf]);
specRawPI = squeeze(ncread(filePI,'/L3_spectra/SH_SPEC',[section 1 1],[1 Inf Inf]));
try
    specCleanPI = squeeze(ncread(filePI,'/L3_spectra/SH_SPEC_CLEAN',[section 1 1],[1 Inf Inf]));
catch
    specCleanPI = specRawPI.*nan;
end
kcycPI = ncread(filePI,'/L3_spectra/KCYC',[section 1],[1 Inf]);

epsiTEST      = ncread(fileTEST,'/L4_dissipation/EPSI',[section 1],[1 Inf]);
specRawTEST   = squeeze(ncread(fileTEST,'/L3_spectra/SH_SPEC',[section 1 1],[1 Inf Inf]));
try
    specCleanTEST = squeeze(ncread(fileTEST,'/L3_spectra/SH_SPEC_CLEAN',[section 1 1],[1 Inf Inf]));
catch
    specCleanTEST=specRawTEST.*nan;
end
kcycTEST      = ncread(fileTEST,'/L3_spectra/KCYC',[section 1],[1 Inf]);
kmaxTEST      = ncread(fileTEST,'L4_dissipation/KMAX',[section 1],[1 Inf]);

epsiL3TEST      = ncread(fileL3TEST,'/L4_dissipation/EPSI',[section 1],[1 Inf]);
specRawL3TEST   = squeeze(ncread(fileL3TEST,'/L3_spectra/SH_SPEC',[section 1 1],[1 Inf Inf]));
try
    specCleanL3TEST = squeeze(ncread(fileL3TEST,'/L3_spectra/SH_SPEC_CLEAN',[section 1 1],[1 Inf Inf]));
catch
    specCleanL3TEST = specRawL3TEST.*nan;
end
kcycL3TEST      = ncread(fileL3TEST,'/L3_spectra/KCYC',[section 1],[1 Inf]);
kmaxL3TEST      = ncread(fileL3TEST,'L4_dissipation/KMAX',[section 1],[1 Inf]);

try
    kvisc = ncread(filePI,'L4_dissipation/KVISC',section,1);
catch
    kvisc =1e-6;
end
n = size(epsiPI,2); % number shear probes

% wvn, min and max wvn
try
    kminPI = ncread(filePI,'L4_dissipation/KMIN',[section 1],[1 Inf]);
catch
    kminPI=0*ones(1,n);
end
kmaxPI = ncread(filePI,'L4_dissipation/KMAX',[section 1],[1 Inf]);

try
    kminTEST = ncread(fileTEST,'L4_dissipation/KMIN',[section 1],[1 Inf]);
catch 
   kminTEST = 0*ones(1,n);
   disp('No KMIN specified.')
end
try
    kminL3TEST = ncread(fileL3TEST,'L4_dissipation/KMIN',[section 1],[1 Inf]);
catch
    kminL3TEST = 0*ones(1,n);
    disp('No KMIN specified.')
end


% plotting
figure('rend','painters','pos',[10 10 1000 1000])
FS = 12;

%raw spectra
for ii=1:n  
    % fitlines
     ktPI = linspace(kminPI(ii),kmaxPI(ii),20);
     PtPI = nasmyth(epsiPI(ii),kvisc,ktPI);
     ktTEST = linspace(kminTEST(ii),kmaxTEST(ii),20);
     PtTEST = nasmyth(epsiTEST(ii),kvisc,ktTEST);
     PtL3TEST = nasmyth(epsiL3TEST(ii),kvisc,ktTEST);
    
    subplot(2,n,ii)
%     fNas=create_nasmyth_baseplot; % create Nasmyth background from CB
    create_nasmyth_baseplot; % create Nasmyth background from CB
    hold on
    loglog(kcycPI,specRawPI(:,ii),'Color',[0.8 0.38 0.08])
    loglog(kcycL3TEST,specRawL3TEST(:,ii),'Color','m')    
    loglog(kcycTEST,specRawTEST(:,ii),'Color','k')    
    loglog(ktPI,PtPI,'Color','g')    
    loglog(ktTEST,PtTEST,'Color','c')    
    loglog(ktTEST,PtL3TEST,'Color','y')    
    hold off
    ylabel('Raw shear spectra (s^-^2 cpm^-^1)')
    box on
    grid on
    xlim([1e-1 1e3])
    ylim([1e-10 1e0])   
    set(gca,'FontSize',FS)
end

% clean spectra
for ii=1:n  
    % fitlines
%     ktPI = linspace(kminPI(ii),kmaxPI(ii),20);
%     PtPI = nasmyth(epsiPI(ii),kvisc,ktPI);
%     ktTEST = linspace(kminTEST(ii),kmaxTEST(ii),20);
%     PtTEST = nasmyth(epsiTEST(ii),kvisc,ktTEST);
    
    subplot(2,n,ii+n)
    create_nasmyth_baseplot; % create Nasmyth background from CB
    hold on
    p1 = loglog(kcycPI,specCleanPI(:,ii),'Color',[0.8 0.38 0.08]);
    p2 = loglog(kcycL3TEST,specCleanL3TEST(:,ii),'Color','m');
    p3 = loglog(kcycTEST,specCleanTEST(:,ii),'Color','k');
    hold off
    ylabel('Clean shear spectra (s^-^2 cpm^-^1)')
    legend([p1 p2, p3],{[pi ,char(949),'=',num2str(epsiPI(ii),'%3.1d')],...
        [tester,' L3,',char(949),'=',num2str(epsiL3TEST(ii),'%3.1d')],...
        [tester,',',char(949),'=',num2str(epsiTEST(ii),'%3.1d')]},...
        'Location','SouthWest','FontSize',FS) 
    box on
    grid on
    xlim([1e-1 1e3])
    ylim([1e-10 1e0])   
    set(gca,'FontSize',FS)
end



end

