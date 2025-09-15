function Att = load_netcdf_attributes(filename)
%% Ilker Fer, 20220512
% Att = load_netcdf_attributes(filename)
% fetch all global attributes from an NC file to a structure Att. 
% Field % names are the attribute names.
% 
%% ALB, 20220522
% Add Groups attribute
%%

A=ncinfo(filename);
for I=1:length(A.Attributes)
    name = A.Attributes(I).Name;
    val = A.Attributes(I).Value;
    Att.Global.(name)=val;
    clear name val
end
%ALB add loop to include Group attribute.
if ~isempty(A.Groups)
    for II=1:length(A.Groups)
        wh_group=A.Groups(II).Name;
        %get attributesat
        if ~isempty(A.Groups(II).Attributes)
            for III=1:length(A.Groups(II).Attributes)
                wh_group_att=A.Groups(II).Attributes(III).Name;
                Att.(wh_group).(wh_group_att)= ...
                    A.Groups(II).Attributes(III).Value;
            end
            %get variable attributes
            for III=1:length(A.Groups(II).Variables)
                wh_var=A.Groups(II).Variables(III).Name;

                Att.(wh_group).variables.(wh_var)=...
                    {A.Groups(II).Variables(III).Dimensions.Name};
                %         Att.(wh_group).variables.(wh_var) = ...
                %             cellfun(@(x) x(2:end),Att.(wh_group).variables.(wh_var),...
                %                     'UniformOutput',false);
            end
            Att.(wh_group).dimensions={A.Groups(II).Dimensions.Name};
        end
    end%end of empty group attribute
end%end of empty group


