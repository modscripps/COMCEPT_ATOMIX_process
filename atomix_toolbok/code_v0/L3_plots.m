function L3_plots(filePI,fileTEST,fileL3TEST,pi,tester,fig_dir,flags)

% clear all
% close all
% clc
% 
% %% set path, data set and PI / tester
% 
% path = '/scratch/kirstin/ATOMIX/data/'; % path to the nc data files
% fig_dir = './figures/'; % output directory to store figures
% dataset = 'VMP250_TidalChannel_024'; % prefix of the nc files
% pi_suffix = ''; % suffix of PI nc file
% %tester_suffix = '_fromL3_IF'; tester = 'IF'; % suffix of test nc file (yours)
% tester_suffix = '_fromL3_ALB'; tester = 'ALB'; % suffix of test nc file (yours)


%% select which figures to plot

plot_L3_spectrum_low10percentile = flags(1);
plot_L3_spectrum_high10percentile = flags(2);
plot_L3_spectrum_midrange = flags(3);
plot_L3_spectrum_any = flags(4);

%% pick segment and run plotting routines

if plot_L3_spectrum_low10percentile==1
    section = lower10(filePI);
    L3_compare_spectra_plot(filePI,fileTEST,fileL3TEST,pi,tester,fig_dir,section)
    print(gcf,'-dpng',[fig_dir 'L3_Spectrum_low10percentile_PI_' tester '.png'],'-r300')
end

if plot_L3_spectrum_high10percentile==1
    section = upper10(filePI);
    L3_compare_spectra_plot(filePI,fileTEST,fileL3TEST,pi,tester,fig_dir,section)    
    print(gcf,'-dpng',[fig_dir 'L3_Spectrum_high10percentile_PI_' tester '.png'],'-r300')
end

if plot_L3_spectrum_midrange==1
    section = mid_range(filePI);
    L3_compare_spectra_plot(filePI,fileTEST,fileL3TEST,pi,tester,fig_dir,section)        
    print(gcf,'-dpng',[fig_dir 'L3_Spectrum_midrange_PI_' tester '.png'],'-r300')
end

if plot_L3_spectrum_any==1
    section = 12;
    % test if within range
    tmpvar = ncread(filePI,'/L4_dissipation/EPSI_FINAL');
    tmpvar = length(tmpvar);
    if section>tmpvar
        disp(['Error: segment must be between 1 and ' num2str(tmpvar)])
        return
    end
    clearvars tmpvar
    L3_compare_spectra_plot(filePI,fileTEST,fileL3TEST,pi,tester,fig_dir,section)        
    print(gcf,'-dpng',[fig_dir 'L3_Spectrum_section' int2str(section) '_PI_' tester '.png'],'-r300')
end