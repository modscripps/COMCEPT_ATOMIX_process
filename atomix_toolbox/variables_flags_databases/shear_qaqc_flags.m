function [Data,Flags]=adv_qaqc_flags(Mask)
% function adv_qaqc_flags(Mask)
% Passes attributes corresponding to various masking flags
% Required inputs:
%		Mask: structure with each field representing 
%Output:
%		Data: structure of UCUR_QC, etc with variables being vector of numeric flags for bad data...
%		Flags: strcuture with fields flag_values and flag_meanings: attributes to be supplied to netcdf file
flag_values=[0 1 2 3 4 ];
flag_meanings='no_qc_performed value_passes_all_tests instrument_status_error ancilliary_outside_valid_range identified_spike' ;


vF=fieldnames(Mask);

maskVar=ones(size(Mask.(vF{1})));

for kk=1:length(vF)
    ind=find(Mask.(vF{kk})>0);
    qaqc=lower(vF{kk});
    switch qaqc
        case{'statmask'}
            maskVar(ind)=2;
        case{'qaqcmask'}
            maskVar(ind)=3;
        case{'despike'}
            maskVar(ind)=4;
    end
end




%% Now append to a structure
vF={'UCUR_QC','VCUR_QC','WCUR_QC'};
for ii=1:length(vF)
   Data.(vF{ii})=maskVar(:,ii); 
   Flags.(vF{ii}).flag_values=flag_values;
    Flags.(vF{ii}).flag_meanings=flag_meanings;
end