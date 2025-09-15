function create_netcdf_fielddata(FieldData,Meta,oname)
% create_netcdf_fielddata(FieldData,Meta,oname,tgrp,timeRef)
%   Accomodates Structures in FieldData where each field is a separate
%   group! e.g., FieldData.summer and FieldData.winter for summer and
%   winter cruises, or these could be different instruments? The name of
%   the structure fields dictate the name of the groups!
% Required inputs:
%	oname: desired  filename for output NetCDF
%	FieldData: struct that I typically save in matlab .mat that I want to convert to netcdf.
%		e.g. FieldData may contain fields with data: PRES,  TIME, TEMP if flat NetCDF
%		or fields that are structures i.e., groups (Level1_raw,
%		Level2_quality_controlled)
%   tgrp=0 or 1 number of levels if the FieldData should be sorted in groups according to
%       each field. tgrp=0 means the NetCDF is "flat", while tgrp=1 means
%       we have one hierarchy level. For ATOMIX data, use tgrp=1.
%
%	Meta: structure with global attributes/metadata that will be writtn to file.
%           This could be nominal depth of water, nominal height of instrument abbove seabed, time deployed and recovered.
%
%   Special fields that should be present at global or group level:
%       Meta.featureType OR Meta.dimensions: dictates the global dimensions via a cell array with names of
%           dimension variables e.g.,Meta.featureType={'TIME','HEIGHT_AB'}
%           featureType can also be a string  (e.g., timeSeries) that will be queried by get_coordinate_type.m to determine the expected dimensions.
%           See get_coordinate_type for different definitions.
%      The global featureType (or dimensions) can be overiden at any group level (e.g.,
%               FieldData.group1) if it contains a field called "featureType".
%       Meta.Flags (purely optional): follows similar principle as Meta.featureType.
%
%   Any group in FieldData with field "Flags" (FieldData.group1.Flags) will overide global Meta.Flags
%           The Flags is structure with fields named "variable_FLAGS" which
%           will append flag attribute information (flag_meanings/flag_values/flag_masks) to the variable
%           FieldData.group1.variable_FLAGS
%           The Flags' fields include only attributes to write into the
%           NetCDF if the variable is found in any of the groups
%
%           For example: "FieldData.group1.UCUR_FLAGS" will have flags that are explained in
%           Meta.Flags.UCUR_FLAGS, which can be overriden by flags explained in
%           FieldData.group1.Flags.UCUR_FLAGS.
%
% Optional input:
%   timeRef: datenum(1950,1,1,0,0,0) is the default time ref applied to ALL
%   variables that contain the word TIME (TIME_HPR, TIME_BNDS).

%% Data checks and initialisations
if isfield(Meta,'Flags')
    DataFlags=Meta.Flags;
    Meta=rmfield(Meta,'Flags');
else
    DataFlags=struct;
end

%% Write global attributes
% [~,oname]=fileparts(oname);
fid = netcdf.create(oname,'netcdf4');
write_netcdf_global(fid,Meta); % Write the global metadata (field experiment, author, etc) to the netcdf file

%% Write groups
vF=fieldnames(FieldData);
nGrp=length(vF);
childId=zeros(nGrp,1);

for ii=1:nGrp
    disp(['Group: ', vF{ii}])
    childId(ii) = netcdf.defGrp(fid,vF{ii}); % define group

    FieldData.(vF{ii}) = set_dimensions(FieldData.(vF{ii}));

    write_group(fid,childId(ii),FieldData.(vF{ii}),Meta);
    
end
% close netcdf
netcdf.close(fid)
close all
end % end of Main

function write_group(fid,childId,GroupData,Meta)


attribute_table=readtable('variables_flags_databases/variables_databases/atomix_netcdf.csv');
Nattribute_csv=size(attribute_table,2);
attribute_name=fieldnames(attribute_table);

% write global group attribute
mF=fieldnames(GroupData.Meta);

groupVarID = netcdf.getConstant('GLOBAL'); % as in the current group

for kk=1:length(mF)
    wh_attribute=mF{kk};
%     netcdf.putAtt(childId,groupVarID, mF{kk},GroupData.Meta.(mF{kk})); % or moored
    switch wh_attribute
        case{'dimensions','variables','dims'}
        otherwise
            netcdf.putAtt(childId,groupVarID, mF{kk},GroupData.Meta.(mF{kk})); % or moored
    end
end



idx_type=find(cellfun(@(x) strcmp(x,"NetcDFType"),fieldnames(attribute_table)));

% Defining dimensions
Ndim=length(GroupData.Meta.dimensions);
% dimension IDs in nc-file
dimid=zeros(Ndim,1);
for i=1:Ndim
    wh_dim=GroupData.Meta.dimensions{i}(1:end);
    wh_dim_size=GroupData.Meta.dims{i};
    % if dim already exisxt in the netcdf defDIM will choke
    % using inqDimID to get the dimID
    try
        dimid(i) = netcdf.inqDimID(childId,wh_dim);
    catch
        dimid(i) = netcdf.defDim(childId,wh_dim,wh_dim_size);
    end

end

% Defining variables
var_names=fieldnames(GroupData.Meta.variables);
Nvar=length(var_names);
lvl_var_id=zeros(Nvar,1);
for i=1:Nvar
    % get the variable name
    wh_var  = var_names{i};
    % get the tyoe of data for the variable (e.g., double,int...)
    idx_attribute_var=find(cellfun(@(x) strcmp(x,wh_var),attribute_table{:,2}));
    wh_type = attribute_table{idx_attribute_var,idx_type};
    if ~isempty(wh_type)
        wh_type=wh_type{1};
        % get the dimensions of the variable (can be a matrix)
        wh_var_dim=GroupData.Meta.variables.(wh_var);
        % initialize the index-of-dimension array.
        wh_dims=zeros(length(wh_var_dim),1);
        for ii=1:length(wh_var_dim)
            idx_dim=find(cellfun(@(x) strcmp(wh_var_dim{ii},x), ...
                GroupData.Meta.dimensions,'un',1));
            wh_dims(ii)= dimid(idx_dim);
        end
        % Define the variable in the NetCDF.
        lvl_var_id(i) = netcdf.defVar(childId,wh_var,wh_type,wh_dims);
    else
        fprintf("No variable %s, pass.\r\n",wh_var)
    end
end

% writing variables to the netcdf
for i=1:Nvar
    % get the variable name
    wh_var  = var_names{i}; 
    idx_attribute_var=find(cellfun(@(x) strcmp(x,wh_var),attribute_table{:,2}));
    if ~isempty(GroupData.(wh_var)) && ~isempty(idx_attribute_var)
        % write the variable inside the NetCDF
        size_data=size(GroupData.(wh_var));
%         if ~isvector(GroupData.(wh_var))
%             GroupData.(wh_var)=permute(GroupData.(wh_var),length(size_data):-1:1);
%         end

        netcdf.putVar(childId,lvl_var_id(i),GroupData.(wh_var));
        % add the attributes for the variable
        for ii=2:Nattribute_csv
            wh_attribute=attribute_table{idx_attribute_var,ii};
            if iscell(wh_attribute)
                wh_attribute=wh_attribute{1};
            end
            if strcmp(wh_attribute,'days, or Days since Meta.origin_of_time')
                wh_attribute=sprintf('days since %s',Meta.origin_of_time);
            end

            sprintf("Write %s-%s\r\n",wh_var,attribute_name{ii})
            netcdf.putAtt(childId,lvl_var_id(i),attribute_name{ii},wh_attribute)
        end
    else
        fprintf("Warning - %s is empty\r\n", wh_var)
    end
end
    
end

function GroupData=set_dimensions(GroupData)
% returns the numerical values of the dimension
disp("Get the dimension size with the data and set them in Meta.dims field\r\n")
Ndims=length(GroupData.Meta.dimensions);
variable_names=fieldnames(GroupData.Meta.variables);
Nvar=length(variable_names);
clear dims
dims=cell(Ndims,1);
for i=1:Nvar
    wh_variable= variable_names{i};
    if isfield(GroupData,wh_variable)
        var_dims=GroupData.Meta.variables.(variable_names{i});
        for ii=1:length(var_dims)
            idx_dims=cell2mat(cellfun(@(y) (strcmp(var_dims{ii},y)), ...
                GroupData.Meta.dimensions,'un',0));
            dims{idx_dims}=size(GroupData.(variable_names{i}),ii);
        end
    end
end
GroupData.Meta.dims=dims;
end





%% Do I need this



%% Sub fcts
% function [dims,newCell]=writegroups(StnDat,childId,featureType,DataFlags,timeRef,Glob)
% if nargin==6
%     StnDat=append_struc(StnDat,Glob);% needed if dimensions are stored/written "global" group...
%     % there's a glitch though since they will get written anyways :( This
%     % code only works for storing global VARIABLES, but not global dim
% end
% 
% 
% [StnDat,varTim]=adjust_time_for_reference(StnDat,timeRef);
% cellData=convert_struc_cell_netcdf(StnDat,DataFlags);
% cellData=correct_time_units_attributes(cellData,varTim,StnDat,timeRef);
% [dims,newCell]=identify_dimension(cellData,featureType);
% 
% if nargin==6
%     %dims=remove_global_cells(dims,Glob);
%     newCell=remove_global_cells(newCell,Glob);
% end
% 
% write_netcdf_var(childId,newCell,dims); 
% end

% %% fct to remove redundant global attributes
% function dims=remove_global_cells(dims,Glob)
% % Remove cells/dims that are already written via glob variable
% vF=fieldnames(Glob);
% cc=0;
% ind=[];
% for ii=1:length(dims)
%     if ~strcmp(dims{ii}.short_name,vF)
%         cc=cc+1;
%         ind(cc)=ii;
%     end
% end
% 
% if ~isempty(ind)
%     dims={dims{ind}};
% end
% end

% %% Fct to correct time attributes units as per timeRef
% %950-01-01T00:00:00Z
% function cellData=correct_time_units_attributes(cellData,varTim,StnDat,timeRef)
% vF=fieldnames(StnDat);
% 
% for ii=1:length(varTim)
%     ind=find(strcmp(varTim{ii},vF));
%     if ~isempty(ind)
%         cellData{ind}.Atts.units=['Days since ',datestr(timeRef,'yyyy-mm-ddTHH:MM:SSZ')];
%     end
% end
% 
% 
% end

%%
% function NewField=add_group_metadata(FieldData,childId)
% % assumes group ID has been called/defined before calling this subfct
% % e.g., childId = netcdf.defGrp(fid,vF{ii});
% %   Looks for structure "Meta", and assigns all fields to netcdf as
% %   metadata at the group level.
% %Output: NewField which is the same as FieldData but with Meta structure removed.
% aF={'Meta','META'}; % so
% 
% for ii=1:length(aF)
%     if isfield(FieldData,aF{ii})
%         
%         [~,TmpMeta]=check_dimension_feature(FieldData.(aF{ii}));
%         
%         
%         mF=fieldnames(TmpMeta);
%         groupVarID = netcdf.getConstant('GLOBAL'); % as in the current group
%         
%         for kk=1:length(mF)
%             netcdf.putAtt(childId,groupVarID, mF{kk},TmpMeta.(mF{kk})); % or moored
%         end
%         FieldData=rmfield(FieldData,aF{ii});
%     end
% end
% 
% NewField=FieldData;
% end

%%

% function dataType=check_featuretype(Data,defaultType)
% % Checks for Meta field in each group, so that different feature types may
% % Input: Data Structure for the group
% %       Meta: global metadata structure to set default
% 
% aF={'Meta','META'}; % sometimes tools will uppercase all variable names
% dataType=defaultType; % default type is set globally
% 
% for ii=1:length(aF)
%     if isfield(Data,aF{ii})
%         [dataType,~,tdef]=check_dimension_feature(Data.(aF{ii}));
%         if tdef==1
%             dataType=defaultType;
%         end
%         % if strcmp(Data.(aF{ii}),'featureType')
%         %    dataType=Data.(aF{ii}).featureType;
%         %end
%     end
% end
% 
% end
%%
%
% function [FieldData,childId]=addnew_group(FieldData,fid,tgrp,lev)
% % lev is the depth we're processing
% vF=fieldnames(FieldData);
% nGrp=length(vF);
%
% for ii=1:nGrp
%     disp(['Group: ', vF{ii}])
%     childId(ii) = netcdf.defGrp(fid,vF{ii}); % define group
%     if tgrp>0 % global one have already been placecd
%         FieldData.(vF{ii})=add_group_metadata(FieldData.(vF{ii}),childId(ii));
%     end
%     %if lev==tgrp
%     writegroups(FieldData.(vF{ii}),childId(ii),'timeSeries',DataFlags)
%     %end
% end
% end

% function [featureType,Meta,tdef]=check_dimension_feature(Meta)
% % FUnction checks if a featureType or dimension was specified
% % Input:
% %   Meta: structure of attributes... the returned structure will habe
% %   dimensions variable removed...
% %   ALB: Adding the "variable" struct for backward compatibility 
% vF={'dimensions','dimensionType','featureType','variables'}; % for backward compatibility
% tDim=isfield(Meta,vF);
% tdef=0;
% if any(tDim)
%     ii=find(tDim);
%     for iii=ii
%         wh_vF=vF{iii};
%         switch wh_vF
%             case 'dimensions'
%             featureType=Meta.(wh_vF);
%         end
%         Meta=rmfield(Meta,wh_vF);
%     end
%     disp("check_dimension_feature(ALB: I removed the existing error msg: Please specify the dimensions field. TODO: check if/when it is needed")
% %     if length(ii)==1
% %         featureType=Meta.(vF{ii});
% %         Meta=rmfield(Meta,vF{ii});
% %     else
% %         error('Please specify the dimensions field within Meta structure. e.g. Meta.dimensions={"Time","Depth"}')
% %     end
% else
%     tdef=1; % using defaults
%     warning('Specify your  dimensions within Meta structure. Will attempt using global defaults or  TIME')
%     %warning('Supply a feature type, see get_coordinate_type.m (e.g., timeSeries)')
%     featureType={'TIME'};%'timeseries';
% end
% 
% 
% end

% function StruR=append_struc(StruR,StrutoAppend,vF)
% % function grabs all the fields in StrutoAppend, and rewrites/appends them
% % to structure StruR.
% % Optional
% %   vF: cell array of fields to append. If noneprovided or empty, all the
% %   fields in StruToAppend will be added
% if nargin<3
%     vF=fieldnames(StrutoAppend);
% end
% 
% if isempty(vF)
%     vF=fieldnames(StrutoAppend);
% end
% 
% %%
% for ii=1:length(vF)
%     if isfield(StruR,vF{ii})
%         disp(['Field ',vF{ii},' exists, so skipping'])
%     else
%         if isfield(StrutoAppend,vF{ii})
%             StruR.(vF{ii})=StrutoAppend.(vF{ii});
%         else
%             disp(['Field ',vF{ii},' doesnt exists'])
%         end
%     end
% end
% end
