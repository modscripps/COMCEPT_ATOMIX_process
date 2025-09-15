function ADVFlags=convert_flags_yaml(fname)
% Reads the yaml files and processes the flag_values and flag_masks such
% that they have the correct format for writing to a NetCDF file.
% ReadYaml doesn't do a great job at recognising number vectors. It stores them as cell arrays.
% CBluteau manually fixed toolbox ReadYaml so that it would check the cell
% arrays (if they were numbers).  Code below may be redundant now
[~,fname]=fileparts(fname);
ADVFlags = ReadYaml([fname,'.yml'],[],1);
gF=fieldnames(ADVFlags);

flagF={'flag_masks','flag_values'};

for ii=1:length(gF)
    Flags=ADVFlags.(gF{ii}); % group's flags
    vF=fieldnames(Flags);
    nVar=length(vF);
    for jj=1:nVar
        
        for kk=1:2
            
            if isfield(Flags.(vF{jj}),flagF{kk})
                if iscell(Flags.(vF{jj}).(flagF{kk}))
                    dataTmp=cell2mat(Flags.(vF{jj}).(flagF{kk}));
                else
                    dataTmp=Flags.(vF{jj}).(flagF{kk});
                end
                
                switch flagF{kk}
                    case{'flag_masks'}
                        dataTmp=uint8(dataTmp);
                    case{'flag_values'}
                        dataTmp=int32(dataTmp);
                end
                Flags.(vF{jj}).(flagF{kk})=dataTmp;
            end
                 
        end
        ADVFlags.(gF{ii})=Flags;
    end
    
end