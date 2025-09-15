function L4_plots(filePI,fileTEST,fileL3TEST,pi,tester,fig_dir,flags)

%% set path, data set and PI / tester

% filePI = [dataset_vertical_structure '.nc'];
% fileTEST = [dataset_tester '.nc
% 
% plot stuff
% [1 1 1 1 1 1];

%% select which figures to plot
close all
plot_L4_epsi_timeseries = flags(1);
plot_L4_epsi_scatter    = flags(2);
plot_L4_ratio_epsi      = flags(3);
plot_L4_fom_epsi        = flags(4);
plot_L4_mad_epsi        = flags(5);
plot_L4_varres_epsi     = flags(6);


%% call plotting routines

if plot_L4_epsi_timeseries==1
   L4_timeseries_Epsilon(filePI,fileTEST,fileL3TEST,pi,tester);
   print(gcf,'-dpng',[fig_dir filePI '_L4_timeseries_Epsilon_' pi '_' tester '.png'],'-r300')
end

if plot_L4_epsi_scatter==1
   L4_scatter_Epsilon(filePI,fileTEST,fileL3TEST,pi,tester);
   print(gcf,'-dpng',[fig_dir filePI 'L4_Epsilon_scatter_' pi '_' tester '.png'],'-r300')
end

if plot_L4_ratio_epsi==1
   L4_ratio_Epsilon(filePI,fileTEST,fileL3TEST,pi,tester);
   print(gcf,'-dpng',[fig_dir filePI 'L4_Epsilon_ratio_' pi '_' tester '.png'],'-r300')
end

if plot_L4_fom_epsi==1
   L4_fom(filePI,fileTEST,fileL3TEST,pi,tester);
   print(gcf,'-dpng',[fig_dir filePI 'L4_fom_' pi '_' tester '.png'],'-r300')
end

if plot_L4_mad_epsi==1
   L4_mad(filePI,fileTEST,pi,tester);
   print(gcf,'-dpng',[fig_dir filePI 'L4_mad_' pi '_' tester '.png'],'-r300')
end

if plot_L4_varres_epsi==1
   L4_varres(filePI,fileTEST,fileL3TEST,pi,tester);
   print(gcf,'-dpng',[fig_dir filePI 'L4_varres_' pi '_' tester '.png'],'-r300')
end

if plot_L4_varres_epsi==1
   L4_kmaxkmin(filePI,fileTEST,pi,tester);
   print(gcf,'-dpng',[fig_dir filePI 'L4_kmaxkmin_' pi '_' tester '.png'],'-r300')
end

if plot_L4_varres_epsi==1
   L4_epsiflag(filePI,fileTEST,pi,tester);
   print(gcf,'-dpng',[fig_dir filePI 'L4_epsiflag_' pi '_' tester '.png'],'-r300')
end




