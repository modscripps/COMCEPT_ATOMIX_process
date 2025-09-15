function write_netcdf_global(fid,Meta)
% write_netcdf_global describing dataset and author of it...
% Essentially wanted a template of info that should always be included,
% I need to improve the code to check if ALL required variables/metadata have been supplied


Meta.disclaimer= 'Data provided as is without warranty';
Meta.license= 'http://creativecommons.org/licenses/by/4.0/';
Meta.date_created=datestr(now,'yyyy/mm/dd HH:MM:SS');
Meta.creation_time=datestr(now,'yyyy/mm/dd HH:MM:SS');
%Meta.Conventions='CF-1.6';

%Required metadata to pass on for featureType: 'timeSeries' i.e., bottom moored

%    instrument: 'Nortek ADV'
%        instrument_serial_number: '03724'
%      instrument_sample_interval: 12
%              geospatial_lat: []
%              geospatial_lon: []
%                instrument_nominal_height: []; OR
%              site_nominal_depth: []
%           time_deployment_start: []
%             time_deployment_end: []

% If featureType='profiler'

%% Optional metadata
%					 abstract: ''
%                        citation: [1x128 char]
%                acknowledgement: [1x434 char]
%      instrument_average_interval: []
%       instrument_burst_interval: []
%    instrument_burst_duration: []

%fid = netcdf.open(strcat(fname,'.nc'),'WRITE');  % Rq: here we use 'WRITE' because the file already exists
NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');

vf=fieldnames(Meta);


for ii=1:length(vf)
    disp(vf{ii})
    if strfind(vf{ii},'time')
        
        tmp=Meta.(vf{ii});
        if isnumeric(tmp)==0 & ischar(tmp)==0
            %disp('stop')
            netcdf.putVar(fid,NC_GLOBAL,vf{ii},datestr(Meta.(vf{ii}),'yyyy-mm-ddTHH:MM:SSZ'));
        else
            netcdf.putAtt(fid,NC_GLOBAL,vf{ii},Meta.(vf{ii}));
        end
        
    else
        netcdf.putAtt(fid,NC_GLOBAL,vf{ii},Meta.(vf{ii})); % or moored
    end
end
disp(' Add input checks of passed metadata');

%netcdf.close(fid)

end


