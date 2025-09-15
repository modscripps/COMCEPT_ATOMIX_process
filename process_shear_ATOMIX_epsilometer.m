%process shear ATOMIX : Epsilometer

% function process_shear_ATOMIX_epsilometer(filename)
%
% This routine intend to load a NETCDF file provided by ATOMIX
% through https://atomix.app.uib.no/Benchmark_datasets_for_shear_probes

% get the data, process them from L1 to L4, and compare the data at
% the different levels. 
% 

%%
%add the needed libraries
addpath(genpath('./atomix_toolbok/'))
%%

dataroot='DATA/';
filename='epsifish_epsilometer_blt_north_atl.nc';
%%

attributes = load_netcdf_attributes(fullfile( ...
    dataroot,filename));
f_attribute= fieldnames(attributes);

Meta       = attributes.Global;
if length(f_attribute)>1
    for f=2:length(f_attribute)
        wh_attribute=f_attribute{f};
        GroupMeta.(wh_attribute)  = attributes.(wh_attribute);
    end
else
    %build Group Meta from the NCfile
    % f_attribute= fieldnames(ncdata);
    A=ncinfo(fullfile(dataroot,filename));
    for g=1:length(A.Groups)
        wh_group=A.Groups(g).Name;
        group_dim=[A.Groups(g).Variables(:).Dimensions];
        group_dim=unique({group_dim(:).Name});
        group_dim=cellfun(@(x) x(2:end),group_dim,'un',0);

        GroupMeta.(wh_group).dimensions = group_dim;%dimension name
        for v=1:length(A.Groups(g).Variables)
            wh_var=(A.Groups(g).Variables(v).Name);
            wh_dim={A.Groups(g).Variables(v).Dimensions.Name};
            GroupMeta.(wh_group).variables.(wh_var)  = cellfun(@(x) x(2:end),wh_dim,'un',0); %dimension name
        end
    end
end

%% ready to load the data
ncdata=load_netcdf_data(fullfile(dataroot,filename));
%% Copy data to process them on your own
FieldData.L1_converted=ncdata.L1_converted;

%% Process the data with MOD routines. 
%  These routine are for the epsilometer. 
%  Here I use ACC (accelerometer) to correct the data. 
%  TODO: Figure out a way to identify which variable should be used for the
%  correction
%  TODO: create a flag so we could choose different processing if someone
% wants to test their own "recipie" 

%  
tic
[FieldData]=process_L1_L2_L3_L4_ATOMIX_ALB(Meta,GroupMeta,FieldData);
toc

%%
