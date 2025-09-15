function [AllData,varName] = load_netcdf_data(filename,nlist)
%function [Data,varName] = load_netcdf(filename,nlist)
%  Load the variable data into a matlab struct, and variable names from a netcdf
%  file into a cellstr "varName"
% Optional Input:
%   nlist: cell array of desired variables if known... Best to use this
%       option for large datasets!
% Improvements
%   At the moment it handles the situations where this is no variables, and
%    only 1 group (netCDF4). It's a matter of adding a for loop for
%    multiple groups but adding the ability to grab variables stored
%    "outside" a group would be more painful
%
% Created by CBluteau in Oct 2017 when I couldn't open a netcdf file in
% ncview/ncbrowse
%%%
if nargin<2
    nlist=[];
end

nInfo=ncinfo(filename);
ncid = netcdf.open(filename,'NOWRITE');
gInfo=nInfo.Groups;


if isempty(gInfo)
    %vInfo=nInfo;
    gid=ncid;
    nG=1;
    gName{1}=[];
else
    %vInfo=nInfo.Groups;
    nG=length(gInfo);
    for kk=1:nG
        gName{kk}=gInfo(kk).Name;
        gid(kk) = netcdf.inqNcid(ncid,gName{kk});
        
    end
end



for kk=1:nG


    varids = netcdf.inqVarIDs(gid(kk));


    for ii=1:length(varids)
        varName{ii}=netcdf.inqVar(gid(kk),varids(ii));
    end

    if isempty(nlist)
        nVar=varName;
    else
       [nVar,varids]=getvarName(varName,varids,nlist);
    end


    for ii=1:length(varids)
        Data.(nVar{ii}) = netcdf.getVar(gid(kk),varids(ii));
        size_data=size(Data.(nVar{ii}));

         if ~isvector(Data.(nVar{ii}))
             Data.(nVar{ii})=permute(Data.(nVar{ii}),length(size_data):-1:1);
         end
    end

    if all(isempty([gName{:}])) & nG==1 %ILKER: inserted [] around gName
       AllData=Data; 
    else
        AllData.(gName{kk})=Data;
    end
    clear Data;
end
end % end of main


function [nVar,vIds]=getvarName(varName,varids,nlist)
    cc=0;
    for ii=1:length(nlist)
        ind=find(strcmp(nlist{ii},varName));
        
        if isempty(ind)
            warning(['Variable ', nlist{ii},'doesnt exist']);
        else 
            cc=cc+1;
            nVar{cc}=varName{ind};
            vIds(cc)=varids(ind);
            
        end
    end

end



