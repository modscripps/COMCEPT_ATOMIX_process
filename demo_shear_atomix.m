% Example demo to read, process and write ATOMIX shear data.
% Load the entire path onto your path
tic
addpath(genpath('../atomix_toolbok'))
% addpath /Volumes/DataDrive4T/TOOLBOXES/seawater/
addpath '~/aleboyer@ucsd.edu - Google Drive/My Drive/TOOLBOXES/seawater/'
%clear all; close all

% FieldData: 
% Structure filled up with our data where each field will be written to its own NetCDF
% group. Each group can have its own dimensions defined via the yml files below.
%
% The routine follows the following steps
% 1/ Global Attibute are defined with the instrument/ dataset yml (e.g., 'epsilometer_metadata.yml')
% % Meta:structure with global metadata. All fields will be written into the
%   NetCDF 

% 2/ The netcdf group are defined in the (hopefully) common specific atomix .yml. 
% This file will be used to defined Group Name, Group Dimensions, Group Variable Name and Dimension
% GroupMeta:structure with Group metadata. All fields will be written into the
%   NetCDF
% Note that: NFFT, DOF ,disslength... are defined in the group attribute.
% 
% 3/ FieldData=process_ATOMIX_data(Meta,GroupMeta) : This is were you
% process your data. THe current file process_ATOMIX_data.m use the Meta
% and GroupMeta to fill up the different Levels using the names defined by the yml. 
% 
% 4/ TODO: Deal with the flags a-la-ATOMIX. As of now: This part of the code is copy/paste of Cynthia's routines.   
% 
% 5/ create_netcdf_fielddata(FieldData,Meta,'shear_atomix_test1'):
% FieldData is a structure with pro!cessed data with names defineds by the
% yml. Write data and attributes in the netcdf file.
% Note the create_netcdf_fielddata read the CSV file 'atomix_netcdf.csv'
% directly extracted from the .xlsx file defiend by the ATOMIX group on
% TEAM. 
% (WARNING: I had to addd a couple of raws to match my yml names.
%           This should be cleared out by the ATOMIX group. )
% 
%% add YML library
addpath(genpath('/Volumes/DataDrive4T/TOOLBOXES/YAMLMatlab_0.4.3'))

%% path to the data
hdir='./TestData'; % Relative path to where the data is located.

%% Read global and group level metadata via yaml input files

Meta = ReadYaml(fullfile(hdir,'epsilometer_metadata.yml'),[],1);
GroupMeta = ReadYaml(fullfile(hdir,'shear_atomix_metada.yml'),[],1);
GroupMeta.L3_spectra.spectral_dof=Meta.spectral_dof;
Meta.disslength_sample=Meta.fftlength_sample*((Meta.spectral_dof+1)/2);
Meta.spectral_disslength_overlap=Meta.disslength_sample.*Meta.spectral_disslength_percent_overlap/100;
Meta.dissoverlap_sec=Meta.spectral_disslength_overlap/Meta.instrument_sample_rate;
Meta.date_update=datestr(now,'dd-mmm-yyyy');
GroupMeta.L3_spectra.spectral_disslength=Meta.disslength_sample;
GroupMeta.L3_spectra.spectral_fft_length=Meta.fftlength_sample;
Meta.num_fft_segments=Meta.spectral_dof;
%% Load field data and global metadata  
% This function should be stand alone so you could actually run A.Le Boyer
% processing routines 

% Prep first group here L1_converted inside FieldDasta. This level is
% instrument dependent and should be customized for every platform/PI
FieldData=prep_L1_epsilometer(Meta,GroupMeta);
%%
FieldData=process_L1_L2_L3_L4_ATOMIX_ALB(Meta,GroupMeta,FieldData);

%% Load optional flags that explain meangings of flags in variables_QC
% You can edit the yml file to suit your definitions. This YAML can also
% contain extra variable attributes that aren't in the CSV databases.
% Anything in this custom yaml will override other databases...
ShearFlags = convert_flags_yaml('variables_flags_databases/qcflags_databases/shear_qc_flags.yml');

% Add flags to FieldData
vF=fieldnames(ShearFlags);
gF=fieldnames(FieldData);
mF=fieldnames(GroupMeta);
nGrp=length(gF); % nbre Groups
for ii=1:nGrp
    if any(strcmp(gF{ii},vF))
        FieldData.(gF{ii}).Flags=ShearFlags.(gF{ii});
    end
    % add Meta structure to the group
    if any(strcmp(gF{ii},mF))
        FieldData.(gF{ii}).Meta=GroupMeta.(gF{ii});
    end
end

%% Write the file
timeRef=datenum(Meta.origin_of_time);

Meta.date_updated=date;
create_netcdf_fielddata(FieldData,Meta,fullfile(hdir,'epsifish_epsilometer_northatlantic_final.nc'));
toc

%%
% cd TestData/ 
% pi='DISS3072';
% tester='DISS4096';
% fig_dir='./figures/';
% filename='epsifish_epsilometer_blt_north_atl.nc';
% ALB_nc_filename='epsifish_epsilometer_northatlantic_1.nc';
% 
% plot_flags=[1 1 1 1 1 1];
% L4_plots(filename,ALB_nc_filename,ALB_nc_filename,pi,tester,fig_dir,plot_flags)

%%
% profile=2
% semilogx(FieldData.L4_dissipation.EPSI(FieldData.L4_dissipation.SECTION_NUMBER==profile),FieldData.L4_dissipation.PRES(FieldData.L4_dissipation.SECTION_NUMBER==profile),'b')
% hold on
% semilogx(ncdata.L4_dissipation.EPSI(ncdata.L4_dissipation.SECTION_NUMBER==profile),ncdata.L4_dissipation.PRES(ncdata.L4_dissipation.SECTION_NUMBER==profile),'r')




