function [FieldData]=process_L3_L4_ATOMIX_ALB(Meta,GroupMeta,FieldData)
% Process levels L4 using FieldData.L3
% using my (Le Boyer) epsilon processing
% The variable names inside L2, L3, L4 are defined in GroupMeta.
% GroupMeta is defined by shear_atomix_metada.yml

% TODO: Figure out a way to desing what variable is used for acceleration correction

Levelnames=fieldnames(GroupMeta);
Levelnames=Levelnames(3:4);

time_slow_flag=0;
if isfield(GroupMeta.(Levelnames{1}),'TIME_CTD')
    % name of the vibration acceleration field
    F1='ACC';F1size=3;
else
    time_slow_flag=1;
    % name of the vibration acceleration field
    F1='VIB';F1size=2;
end


%% diss_length and fft_length criteria
FS         = Meta.instrument_sample_rate; %sampling freq. should be FS  = 1./nanmean(diff(L1_converted.TIME).*86400)
FS_CTD     = Meta.instrument_ctd_sample_rate;
nfft       = Meta.fft_length_sample;
dof        = Meta.spectral_dof; % dof has to be odd
disslength = 0.5 .* nfft * (dof+1) ; %if dof = 5 you only 3 nfft length with 50% overlap   
disslength_overlap=Meta.spectral_disslength_overlap/100; % percentage of disslength overlap
FPUMP      = Meta.instrument_ctd_pump_freq; % CTD pump 50Hz;
FN         = FS/2;
equivalent_dof= 8/3 * disslength./ nfft ./ 2;

[~, f]=pwelch(1:nfft,nfft,[],nfft,FS);



%% fix issues with netcdf
switch Meta.contact
    case 'rinoue@jamstec.go.jp'
        GroupMeta.L3_spectra.dimensions={'N_SHEAR_SENSORS','N_ACC_SENSORS','TIME_SPECTRA','WAVENUMBER'};
        GroupMeta.L3_spectra.variables.SH_SPEC{3}=GroupMeta.L3_spectra.dimensions{1}
end
%% get Transfer functions
switch Meta.instrument
    case "Epsilometer"
        H.samplingfilter=(sinc(f/(2*f(end)))).^4;
        ca_filter = load('TestData/FILTER/charge_coeffilt.mat'); %from network analysis
        H.epsi_ca   = interp1(ca_filter.freq,ca_filter.coef_filt ,f);
        % H.epsi_ca   = f.*0 +1;
        H.gain_ca   = 1;
        H.electshear= H.epsi_ca*H.gain_ca;%
        H.gainshear=1;
        H.adcshear=H.gainshear.* H.samplingfilter;
        H.shear=(H.electshear .* H.adcshear).^2;
        H.diff = (2*pi.*f./(4 .* FN .* sin (pi/2 .* f./FN))).^2;
        H.shear=(H.electshear .* H.adcshear).^2+H.diff;
    case "MSS"
        H.samplingfilter=f.*0+1;
        H.shear=f.*0+1;
    case {"Rockland","Microrider"}
        H.samplingfilter=f.*0+1;
        H.diff = (2*pi.*f./(4 .* FN .* sin (pi/2 .* f./FN))).^2;

        H.shear=f.*0+H.diff;
end

%% Process profile

% get group dims to initialize arrays
ndim=length(GroupMeta.(Levelnames{1}).dimensions);
varnames=fieldnames(GroupMeta.(Levelnames{1}).variables);
nvar=length(varnames);
for d=1:ndim
    wh_dim=GroupMeta.(Levelnames{1}).dimensions{d}(2:end);
    for v=1:nvar
        wh_var=varnames{v};
        idx_dim=cellfun(@(x) strcmp(x(2:end),(wh_dim)),GroupMeta.(Levelnames{1}).variables.(wh_var));
        if ~isempty(find(idx_dim))
            dims.(wh_dim)=size(FieldData.(Levelnames{1}).(wh_var),find(idx_dim));
        end
    end
end

%initialize array if field does not exist (in the case the nc file was created by me)
for l=1:length(Levelnames)
        if ~isfield(FieldData,Levelnames{l})
            FieldData.(Levelnames{l})=[];
        end
    varnames=fieldnames(GroupMeta.(Levelnames{l}).variables);
    nvar=length(varnames);
    for v=1:nvar
        wh_var=varnames{v};
        if ~isfield(FieldData.(Levelnames{l}),wh_var)
            ndim=length(GroupMeta.(Levelnames{l}).variables.(wh_var));
            var_size=ones(1,max(ndim,2));
            for d=1:ndim
                wh_dim=GroupMeta.(Levelnames{l}).variables.(wh_var){d}(2:end);
                try
                var_size(d)=dims.(wh_dim);
                catch
                    fprintf("%s - no dim %s, set it to 1 \r\n",(Levelnames{l}),wh_dim)
                    var_size(d)=1;
                end
            end
            FieldData.(Levelnames{l}).(wh_var)=zeros(var_size);
        end
    end
end

% process the spectra
n_disslength=dims.TIME_SPECTRA;
N_SHEAR_SENSOR=size(FieldData.L3_spectra.SH_SPEC,3);

try
    N_VIB_SENSOR=size(FieldData.L3_spectra.VIB_SPEC,3);
catch
    N_VIB_SENSOR=0;
end

try
    N_ACC_SENSOR=size(FieldData.L3_spectra.ACC_SPEC,3);
catch
    N_ACC_SENSOR=0;
end
try
    N_AP_SENSOR=size(FieldData.L3_spectra.AP_SPEC,3);
catch
    N_AP_SENSOR=0;
end


equivalent_dof=dof-1;
s    = 31; %arbitrary

for p=1:n_disslength

        w    = FieldData.L3_spectra.PSPD_REL(p);
        pres = FieldData.L3_spectra.PRES(p);
        try
            t = FieldData.L3_spectra.TEMP(p);
        catch
            t=20;
        end

        % compute kvis 
        kvis = nu(s,t,db2MPa(pres));
        k=FieldData.L3_spectra.KCYC(p,:);
        kmax      = FPUMP./w;


        for n=1:N_SHEAR_SENSOR
            SH_SPEC(:,n)=squeeze(FieldData.L3_spectra.SH_SPEC(p,:,n));
            try
                Ps_clean_f(:,n)=squeeze(FieldData.L3_spectra.SH_SPEC_CLEAN(p,:,n));
            catch
                Ps_clean_f(:,n)=SH_SPEC(:,n);
            end
            SH_SPEC_CLEAN(:,n)=Ps_clean_f(:,n);
            Psheark(:,n) = squeeze(SH_SPEC_CLEAN(:,n));

            try
            [epsilon(n),kmin(n),kc(n),method(n)] = eps1_mmp(k,Psheark(:,n),kvis,kmax);
            catch
                disp('bobo')
            end

            if epsilon(n)==1e-11 
                disp("no epsilon. Why?")
            end

            Psheark_intrange=Psheark(k>=kmin(n) & k<=kc(n),n);
            k_intrange=k(k>=kmin(n) & k<=kc(n));
            [fom(n),var_resolve(n),mad_spec(n)] = ...
                get_panchev_fom(epsilon(n),...
                                Psheark_intrange, ...
                                k_intrange, ...
                                kvis, ...
                                equivalent_dof);


        end % end of N_SHEAR_SENSOR loop
        for a=1:N_VIB_SENSOR
            AP_SPEC(:,a)=squeeze(FieldData.L3_spectra.VIB_SPEC(p,:,a));
        end
        for a=1:N_AP_SENSOR
            AP_SPEC(:,a)=squeeze(FieldData.L3_spectra.AP_SPEC(p,:,a));
        end
        for a=1:N_ACC_SENSOR
            AP_SPEC(:,a)=squeeze(FieldData.L3_spectra.ACC_SPEC(p,:,a));
        end
        if N_VIB_SENSOR==0 && N_AP_SENSOR==0
            AP_SPEC=squeeze(FieldData.L3_spectra.SH_SPEC(p,:,1)).*nan;
        end

%         SH_SPEC1=squeeze(FieldData.L3_spectra.SH_SPEC(p,:,1));
%         SH_SPEC2=squeeze(FieldData.L3_spectra.SH_SPEC(p,:,2));
% 
%         try
%             Ps1_clean_f=squeeze(FieldData.L3_spectra.SH_SPEC_CLEAN(p,:,1)).';
%             Ps2_clean_f=squeeze(FieldData.L3_spectra.SH_SPEC_CLEAN(p,:,2)).';
%         catch
%             Ps1_clean_f=SH_SPEC1.';
%             Ps2_clean_f=SH_SPEC2.';
%         end

%         % shear clean convert frequency spectra to wavenumber spectra
%         SH_SPEC_CLEAN1=Ps1_clean_f;
%         SH_SPEC_CLEAN2=Ps2_clean_f;
        

%         try
%             AP_SPEC1=squeeze(FieldData.L3_spectra.VIB_SPEC(p,:,1));
%             AP_SPEC2=squeeze(FieldData.L3_spectra.VIB_SPEC(p,:,2));
%         catch
%             try
%                 AP_SPEC1=squeeze(FieldData.L3_spectra.ACC_SPEC(p,:,1));
%                 AP_SPEC2=squeeze(FieldData.L3_spectra.ACC_SPEC(p,:,2));
%             catch
%                 try
%                 AP_SPEC1=squeeze(FieldData.L3_spectra.AP_SPEC(p,:,1));
%                 AP_SPEC2=squeeze(FieldData.L3_spectra.AP_SPEC(p,:,2));
%                 catch
%                     AP_SPEC1=SH_SPEC1.*nan;
%                     AP_SPEC2=SH_SPEC1.*nan;
%                 end
% 
%             end
%         end

        % compute epsilon and cutoff wavenumber
%         
%         Ps1sheark = squeeze(SH_SPEC_CLEAN1);
%         Ps2sheark = squeeze(SH_SPEC_CLEAN2);
%         k         = f./w;

        % compute epsilon estimates

%         kmax=FieldData.PI_KMAX(p,1);
%         kmin=FieldData.PI_KMIN(p,1);
%         [epsilon1,kmin1,kc1,method1] = eps1_mmp_ALB(k,Ps1sheark,kvis,kmax,kmin);
%     disp(p)
% if p==2
%     disp("issue")
% end
%         [epsilon1,kmin1,kc1,method1] = eps1_mmp(k,Ps1sheark,kvis,kmax);

%         kmax=FieldData.PI_KMAX(p,2);
%         kmin=FieldData.PI_KMIN(p,2);
%         [epsilon2,kmin2,kc2,method2] = eps1_mmp(k,Ps2sheark,kvis,kmax,kmin);
% try
%         [epsilon2,kmin2,kc2,method2] = eps1_mmp(k,Ps2sheark,kvis,kmax);
% catch
%     disp('bobo')
% end


%         if epsilon1==1e-11 || epsilon2==1e-11
%             disp("no epsilon. Why?")
%         end



        % qc flag
        % get_qc_flag(epsilon,fom,fom_limit,despike,despike_limit,varepsilon,diss_ratio_limit,var_resolve,var_resolve_limit)
        fom_limit=3;
        despike_limit=.4; %40%
        var_resolve_limit=.5; %50%
        varepsilon=nan; %TODO deal with varepsilon and epsi_std. Right now I do not know how.
        diss_ratio_limit=nan;
        despike=epsilon.*nan;
        
        qc=get_qc_flag(epsilon, ...
                       fom,fom_limit,...
                       despike,despike_limit,...
                       varepsilon,...
                       diss_ratio_limit,...
                       var_resolve,var_resolve_limit);

        % final epsilon product
        [epsi_final, ~,epsi_std]=get_final_epsilon(epsilon,qc);


        variable_names=fieldnames(GroupMeta.(Levelnames{1}).variables);
        for n=1:length(variable_names)
            wh_var=variable_names{n};
            switch wh_var
                case "dimension"
                case "TIME"
                    FieldData.(Levelnames{1}).TIME(p)    = FieldData.(Levelnames{1}).TIME(p);
                case {"FLOWPS","PSPD_REL"}
                    FieldData.(Levelnames{1}).PSPD_REL(p)  = w;
                case "PRES"
                    FieldData.(Levelnames{1}).PRES(p)    = pres;
                case "TEMP"
                    FieldData.(Levelnames{1}).TEMP(p)    = t;
                case "PSAL"
                    FieldData.(Levelnames{1}).PSAL(p)     = s;
                case "FREQ"
                    FieldData.(Levelnames{1}).FREQ(:)    = f;
                case "KCYC"
                    try
                        FieldData.(Levelnames{1}).KCYC(p,:)  = f./w;
                    catch
                        FieldData.(Levelnames{1}).KCYC(p,:)  = k;
                    end
                case "ACC_SPEC"
                    FieldData.(Levelnames{1}).ACC_SPEC(p,:,:)=AP_SPEC;
%                     try
%                     FieldData.(Levelnames{1}).ACC_SPEC(p,:,3)=AP_SPEC3;
%                     catch
%                     end
                case "VIB_SPEC"
                    FieldData.(Levelnames{1}).VIB_SPEC(p,:,:)=AP_SPEC;
%                     FieldData.(Levelnames{1}).VIB_SPEC(p,:,1)=AP_SPEC1;
%                     FieldData.(Levelnames{1}).VIB_SPEC(p,:,2)=AP_SPEC2;
                case "SH_SPEC"
                    FieldData.(Levelnames{1}).SH_SPEC(p,:,:)=SH_SPEC;
%                     FieldData.(Levelnames{1}).SH_SPEC(p,:,1)=SH_SPEC1;
%                     FieldData.(Levelnames{1}).SH_SPEC(p,:,2)=SH_SPEC2;
                case "SH_SPEC_CLEAN"
                    FieldData.(Levelnames{1}).SH_SPEC_CLEAN(p,:,:)=SH_SPEC_CLEAN;
%                     FieldData.(Levelnames{1}).SH_SPEC_CLEAN(p,:,1)=SH_SPEC_CLEAN1;
%                     FieldData.(Levelnames{1}).SH_SPEC_CLEAN(p,:,2)=SH_SPEC_CLEAN2;
                case "CA_TF"
                    FieldData.(Levelnames{1}).CA_TF=H.shear;
                case "SECTION_NUMBER"
                    FieldData.(Levelnames{1}).SECTION_NUMBER(p)=FieldData.(Levelnames{1}).SECTION_NUMBER(p);
                case "DOF"
                    FieldData.(Levelnames{1}).DOF(p)=equivalent_dof;
            end
        end

        %L4 level
        variable_names=fieldnames(GroupMeta.(Levelnames{2}).variables);
        for n=1:length(variable_names)
            wh_var=variable_names{n};
            switch wh_var
                case "dimension"
                case "TIME"
                    FieldData.(Levelnames{2}).TIME(p)   = FieldData.(Levelnames{1}).TIME(p);
                case {"FLOWPS","PSPD_REL"}
                    FieldData.(Levelnames{2}).PSPD_REL(p) = FieldData.(Levelnames{1}).PSPD_REL(p);
                case "PRES"
                    FieldData.(Levelnames{2}).PRES(p)   = FieldData.(Levelnames{1}).PRES(p);
                case "TEMP"
                    FieldData.(Levelnames{2}).TEMP(p)   = FieldData.(Levelnames{1}).TEMP(p);
                case "PSAL"
                    FieldData.(Levelnames{2}).PSAL(p)    = FieldData.(Levelnames{1}).PSAL(p);
                case "KVISC"
                    FieldData.(Levelnames{2}).KVISC(p)  = kvis;
                case "KMAX"
                        FieldData.(Levelnames{2}).KMAX(p,:,:)   = kc;
%                     try
%                         FieldData.(Levelnames{2}).KMAX(p,:)   = [kc1 kc2];
%                     catch
%                         FieldData.(Levelnames{2}).KMAX(p)   = kc1;
%                     end
                case "KMIN"
                    FieldData.(Levelnames{2}).KMIN(p,:,:)   = kmin;
%                     FieldData.(Levelnames{2}).KMIN(p,:)   = [kmin1 kmin2];
                case "SECTION_NUMBER"
                    FieldData.(Levelnames{2}).SECTION_NUMBER(p) = FieldData.(Levelnames{1}).SECTION_NUMBER(p);
                case "FOM"
                    FieldData.(Levelnames{2}).FOM(p,:,:)     = fom;
%                     FieldData.(Levelnames{2}).FOM(p,:)     = [fom1 fom2];
                case "MAD"
                    FieldData.(Levelnames{2}).MAD(p,:,:)     = mad_spec;
                    FieldData.(Levelnames{2}).VAR_RESOLVED(p,:,:)  = var_resolve;
%                     FieldData.(Levelnames{2}).MAD(p,:)     = [mad_spec1 mad_spec2];
%                     FieldData.(Levelnames{2}).VAR_RESOLVED(p,:)  = [var_resolve1 var_resolve2];
                case "METHOD"
                    FieldData.(Levelnames{2}).METHOD(p,:,:)  = method;
%                     FieldData.(Levelnames{2}).METHOD(p,:)  = [method1 method2];
                case "VAR_RESOLVED"
                    FieldData.(Levelnames{2}).VAR_RESOLVED(p,:,:)  = var_resolve;
%                     FieldData.(Levelnames{2}).VAR_RESOLVED(p,:)  = [var_resolve1 var_resolve2];
                case "EPSI"
                        FieldData.(Levelnames{2}).EPSI(p,:,:)    = epsilon;
%                     try
%                         FieldData.(Levelnames{2}).EPSI(p,:)    = [epsilon1 epsilon2];
%                     catch
%                         FieldData.(Levelnames{2}).EPSI(p)    = epsilon1;
%                     end
                case {"EPSI_QC","EPSI_FLAGS"}
                    FieldData.(Levelnames{2}).EPSI_FLAGS(p,:,:) = qc;
                case "EPSI_FINAL"
                    FieldData.(Levelnames{2}).EPSI_FINAL(p)    = epsi_final;
                case "EPSI_STD"
                    FieldData.(Levelnames{2}).EPSI_STD(p,:)    = epsi_std;
            end
        end

end % end for p=1:P

% transpose if needed
for l=1:length(Levelnames)
    variable_names=fieldnames(GroupMeta.(Levelnames{l}).variables);
    for n=1:length(variable_names)
        wh_var=variable_names{n};
        if isfield(FieldData.(Levelnames{l}),wh_var)
            if isrow(FieldData.(Levelnames{l}).(wh_var))
                FieldData.(Levelnames{l}).(wh_var)=....
                    FieldData.(Levelnames{l}).(wh_var)(:);
            end
        else
            fprintf('Warning - no field %s\r\n',wh_var)
        end
    end
end


end

function y = haf_oakey(f,w)
% haf_oakey
%   Usage: y = haf_oakey(f,w)
%	   f is frequency in a column vector, 
%      w is the fall rate
%   Function: Power transfer function of airfoil probe according
%   to Oakey

lc=0.02;
y = 1 ./ (1 + (lc .* (f' / w) ).^2 );
y=y';
end

function [epsilon,kmin,kc,method]=eps1_mmp(k,Psheark,kvis,kmax)
% eps1_mmp
%   Usage: epsilon=eps1_mmp(k,Psheark,kvis,w);
%      k is a vector array with wavenumber in cpm
%      Psheark is a vector array with the shear spectrum
%      kvis is the kinematic viscocity, in m^2/s
%      w is the vehicle speed, in m/s
%      dk is the elementary waveSD! 1enumber determined by eps1_mmp
%      epsilon is the estimated dissipation rate, in W/kg
%      kc is the wavenumber at which integration is cutoff, in cpm
%   Function: To integrate airfoil shear versus k spectra to
%      obtain epsilon.  The algorithm is similar to that of
%      Wesson & Gregg, 1994, JGR, 99, -9877, but uses Panchev's
%      universal spectrum for reference and stops integration to 
%      avoid vibration peaks.  This limit, kmax, is determined by eps2_mmp.

%    Epsilon is determined iteratively, in three stages.
%    (1) integration between 2 and 10 cpm, eps1
%          The observed spectrum is interpolated onto a wavenumber grid

%       between 2 and 10 cpm, with 0.2 cpm steps.  The integral, shear 10,
%       is compared with an integral of Panchev's universal spectrum over
%       the same grid.  If log10(shear10)>-3, 2 to 10 cpm lies solely in
%       the inertial subrange, which does not depend on viscosity.  Epsilon
%       is obtained directly from a polynomial fit of epsilon to shear10 for
%       the Panchev spectrum.  If log10(shear10)<=-3, 2 to 10 cpm contains at
%       least some of the viscous rolloff and epsilon/7.5 nu is obtained from
%       a polynomial fit to minimize errors due to viscosity.
%
%    (2) integration to the wavenumber containing 90% variance of Panchev's
%       spectrum evaluated with eps1 and nu, eps2
%         The upper limit for integration is reduced if it exceeds kmax, the limit
%       determined for noise-free spectra by script eps2_mmp.m.  The lower bound
%       is at 1 cpmk.  If fewer than 4 spectral estimates are in the wavenumber band, 
%       no integration is done and epsilon is set to 1e-10, taken as the base level.  
%       The estimate is raised by script epsilon_correct.m if the signal has been 
%       reduced by probe attenuation.  
%
%    (3) repeat of (2) with wavenumber determined from eps2    

method=0; % we are not in the inertial subrange. 

% kmin=2;
if k(1)==0
    kmin=k(2);
else
    kmin=k(1);
end
% dKI=0.2;
% KI=(2:dKI:10); % wavenumber array for interpolation
% 
% dk=nanmean(diff(k));
% first estimation of epsilon of the sum of the shear variance is too high
% and falls into the inertial subrange (no roll off). I think it assumes that
% the shear variance directly follows a Panchev spectrum and predict
% epsilon directly without influence of the viscosity. I have no idea how
% these coefficients are computed.
% 
eps_fit_shear10=[8.6819e-04, -3.4473e-03, -1.3373e-03, 1.5248, -3.1607];

% If the shear variance include a part of the roll-off, we use these coef
% which (I think) give a panchev spectrum for a given shear variance value 
shtotal_fit_shear10=[6.9006e-04, -4.2461e-03, -7.0832e-04, 1.5275, 1.8564];


% first estimate, using Psheark between 2 & 10 cpm
% Interpolate Psheark onto 0.2 cpm grid & integrate
% Only this estimate is interpolated, as it is the only one input to
% a polynomial integrated with the same grid
krange=find(k>=kmin & k<10); 
% P_interpolated=interp1(k(krange),Psheark(krange),KI);
% ALB change to nansum since coherence correction can introduces nansclear 
% shear10=nansum(P_interpolated)*0.2;
% shear10=nansum(P_interpolated)*dKI;
shear10=trapz(k(krange),Psheark(krange));
%
% estimate epsilon using poly fits to log10(shear10)
logshear10=log10(shear10);
if logshear10>-3 % 2-10 cpm lies entirely in inertial subrange
	log10eps=polyval(eps_fit_shear10,logshear10);
	eps1=10^log10eps;
    method=1;% we ARE in the inertial subrange. 
else
	log10_sheartotal=polyval(shtotal_fit_shear10,logshear10);
	eps1=7.5*kvis*10^log10_sheartotal;
end

% second estimate: we use the first estimate of epsilon to find the
% kolmogrov scale and re-do the integration.
% kc2 = 0.0816*( eps1  / kvis^3 )^(1/4);  % k for 90% variance of Panchev spectrum
kc_coef=.13;
kc2 = kc_coef*( eps1  / kvis^3 )^(1/4);  % k for 90% variance of Panchev spectrum
if kc2>kmax
	kc2=kmax; % limit set by noise spectrum
end
krange=find(k>=kmin & k<=kc2);
% eps2=7.5*kvis*nansum(Psheark(krange))*dk/.9; % .9 we want to get 90% of the shear variance
eps2=7.5*kvis*trapz(k(krange),Psheark(krange));

% third estimate: same as before.
% kc=0.0816*( eps2 / kvis^3 )^(1/4);
kc = kc_coef*( eps2  / kvis^3 )^(1/4);  % k for 90% variance of Panchev spectrum
if kc > kmax
	kc=kmax;
end
krange=find(k>=kmin & k<=kc);
if length(krange)>2
    eps3=7.5*kvis*trapz(k(krange),Psheark(krange));
else
    eps3=1e-12;
end

if eps3<1e-11 
	epsilon=1e-11;
else
  mf=epsilon2_correct_ALB(eps3,kvis,kmin,kc);
	epsilon=mf*eps3;
end

kc=k(find(k>kc,1,'first')-1);

end



function [epsilon,kmin,kc,method]=eps1_mmp_ALB(k,Psheark,kvis,kmax,kmin)
% eps1_mmp
%   Usage: epsilon=eps1_mmp(k,Psheark,kvis,w);
%      k is a vector array with wavenumber in cpm
%      Psheark is a vector array with the shear spectrum
%      kvis is the kinematic viscocity, in m^2/s
%      w is the vehicle speed, in m/s
%      dk is the elementary waveSD! 1enumber determined by eps1_mmp
%      epsilon is the estimated dissipation rate, in W/kg
%      kc is the wavenumber at which integration is cutoff, in cpm
%   Function: To integrate airfoil shear versus k spectra to
%      obtain epsilon.  The algorithm is similar to that of
%      Wesson & Gregg, 1994, JGR, 99, -9877, but uses Panchev's
%      universal spectrum for reference and stops integration to 
%      avoid vibration peaks.  This limit, kmax, is determined by eps2_mmp.

%    Epsilon is determined iteratively, in three stages.
%    (1) integration between 2 and 10 cpm, eps1
%          The observed spectrum is interpolated onto a wavenumber grid

%       between 2 and 10 cpm, with 0.2 cpm steps.  The integral, shear 10,
%       is compared with an integral of Panchev's universal spectrum over
%       the same grid.  If log10(shear10)>-3, 2 to 10 cpm lies solely in
%       the inertial subrange, which does not depend on viscosity.  Epsilon
%       is obtained directly from a polynomial fit of epsilon to shear10 for
%       the Panchev spectrum.  If log10(shear10)<=-3, 2 to 10 cpm contains at
%       least some of the viscous rolloff and epsilon/7.5 nu is obtained from
%       a polynomial fit to minimize errors due to viscosity.
%
%    (2) integration to the wavenumber containing 90% variance of Panchev's
%       spectrum evaluated with eps1 and nu, eps2
%         The upper limit for integration is reduced if it exceeds kmax, the limit
%       determined for noise-free spectra by script eps2_mmp.m.  The lower bound
%       is at 1 cpmk.  If fewer than 4 spectral estimates are in the wavenumber band, 
%       no integration is done and epsilon is set to 1e-10, taken as the base level.  
%       The estimate is raised by script epsilon_correct.m if the signal has been 
%       reduced by probe attenuation.  
%
%    (3) repeat of (2) with wavenumber determined from eps2    

method=0; % we are not in the inertial subrange. 

% kmin=0;
dKI=0.2;
KI=(2:dKI:10); % wavenumber array for interpolation

dk=nanmean(diff(k));

% third estimate: same as before.
kc=kmax;
if kc > kmax
	kc=kmax;
end
krange=find(k>=kmin & k<=kc);
eps3 = 7.5*kvis*nansum(Psheark(krange))*dk; 

if eps3<1e-11 || length(krange) < 4
	epsilon=1e-11;
else
  mf=epsilon2_correct_ALB(eps3,kvis,kc);
	epsilon=mf*eps3;
end

%selecting the closest k for gc;
kc=k(find(k>kc,1,'first')-1);

end


function [MPa] = db2MPa(db)
% [MPa] = db2MPa(db)
% Converts pressure in db to pressure in MPa
MPa = db/100;
end

function f=epsilon2_correct(eps_res,kvis)

pf=[-3.12067617e-05, -1.31492870e-03, -2.29864661e-02, -2.18323426e-01, -1.23597906, ...
    -4.29137352,-8.91987933, -9.58856889, -2.41486526];

ler=log10(eps_res);
	
if ler <= -6 % no correction needed
	lf_eps=0;
elseif ler>-6 && ler<=-1 % range of fit
	lf_eps=polyval(pf,ler);	
	if ler < -2 % apply viscosity correction
 		lf_eps=lf_eps+0.05*(ler+6)*(1.3e-6-kvis)/(0.3e-6);
	end
elseif ler>-1
	lf_eps=polyval(pf,-1);
end

f=10.^lf_eps;
end

function f=epsilon2_correct_ALB(eps_res,kvis,kmin,kmax)

pf=[-3.12067617e-05, -1.31492870e-03, -2.29864661e-02, -2.18323426e-01, -1.23597906, ...
    -4.29137352,-8.91987933, -9.58856889, -2.41486526];

ler=log10(eps_res);
	
if ler <= -6 % no correction needed
	lf_eps=0;

    [kpan,Ppan] = panchev(eps_res,kvis);
%     [Ppan,kpan]=nasmyth(eps_res,kvis);

    idx_kmax2=find(kpan>=kmax,1,'first');
    idx_kmax1=find(kpan>=kmin,1,'first');
    idx_range=idx_kmax1:idx_kmax2;
    tot_var=trapz(kpan(idx_kmax1:end),Ppan(idx_kmax1:end));

    var_explain=trapz(kpan(idx_range),Ppan(idx_range));
    lf_eps=log10(2-var_explain/tot_var);

elseif ler>-6 && ler<=-1 % range of fit
	lf_eps=polyval(pf,ler);	
	if ler < -2 % apply viscosity correction
 		lf_eps=lf_eps+0.05*(ler+6)*(1.3e-6-kvis)/(0.3e-6);
        %ALB 12/09/2022 I have no idea where the coefs come from but if want to match Rolf
        %Lueck results I need to do 
        lf_eps=1.6.*lf_eps;

	end
elseif ler>-1
	lf_eps=polyval(pf,-1);
end

f=10.^lf_eps;
end

function [fom,percent_var_explain,mad_spec] = get_panchev_fom(epsilon,Psheark,k,kvis,dof)
%sig_lnS=5/4*dof^(-7/9);
%here dof is number of segment use (our classic dof) 
% minus the number of channel used to correct for acceleration
% https://wiki.app.uib.no/atomix/index.php/Figure_of_merit_(FM)
% k should be the range used for the spectral integration kmin<=k<=kmax
% Psheark should be the spectrum(kmin<=k<=kmax)

sig_lnS=sqrt(5/4*dof^(-7/9));
Ns=length(k);
Tm=0.8+sqrt(1.56/Ns);
[kpan,Ppan] = panchev(epsilon,kvis);
interp_Ppan=interp1(kpan(~isnan(Ppan)),Ppan(~isnan(Ppan)),k);

fom=log(Psheark(:)./interp_Ppan(:));
mad_spec=mad(fom);
fom=mad_spec./sig_lnS./Tm;

tot_var=trapz(kpan,Ppan);
idx=find(kpan>=k(1) & kpan<=k(end));
var_explain=trapz(kpan(idx),Ppan(idx));

percent_var_explain=var_explain./tot_var;

end

function [k,Pxx] = panchev(epsilon,kvis,kin) 
% panchev
%   Usage: [k,Pxx]=panchev(epsilon,kvis,kin);
%     inputs
%      epsilon is the dissipation rate in W/kg
%      kvis is the kinematic viscosity in m^2/s
%      kin, vector of wavenumbers in cpm, is an optional argument
%     outputs
%      k, a vector of wavenumbers,  is kin if it is specified.  Otherwise
%         it is uniformly spaced between 1/(1000*eta) and 1/(5*eta)
%      Pxx is Panchev's theoretical spectrum of transverse shear 
%   Function: evaluate and return Panchev's universal turbulent
%      spectrum for a specified wavenumber vector, if specified, or
%      for one uniformly spaced to span the lower inertial subrange and
%      the dissipation range.

%   Rewritten in function form by kw 8/5/94, modified by mg 1/16/95, 4/25/95,
%     9/15/95
%   Revised to partially vectorize calculation, decreasing cpu
%     time 18-fold for no kin. 12/16/96

% set values for constants and zeta integration
eta=(kvis^3/epsilon)^(1/4); % Kolmogoroff length scale, in meters
a = 1.6;
delta = 0.1;
conv = 1/(eta*2*pi);
c32 = 3/2;
sc32 = sqrt(c32);
c23 = 2/3;
ac32 = a^c32;
c43 = 4/3;


% if kin is not given,set wavenumber range, in cpm, based on eta
if nargin < 3
	k0=1/(10000*eta);
	kmax=1/(5*eta);
	dk=(kmax-k0)/200;
	kin=k0:dk:kmax;
end

%  do the computation:
scale=2*pi*(epsilon*kvis^5)^0.25;
sum=zeros(size(kin));
kn=kin/conv;
zeta=delta/2:delta:1;
nzeta=length(zeta);
for i=1:nzeta
  z=zeta(i);
	sum=sum ...
	+ delta .* (1+z.^2) ...
	.* (a.*z.^c23 + (sc32.*ac32) .* kn.^c23) ...
	.* exp(-c32.*a.*(kn/z).^c43 - (sc32.*ac32.*(kn/z).^2));	
end
phi = 0.5.*kn.^(-5/3).*sum;
Pxx = scale.*(kn/eta).^2.*phi;

k=kin(:); Pxx=Pxx(:);

end

function y = nu(s,t,p)
%nu: kvis=nu(s,t,p) - kinematic viscosity
%  Inputs: s: salinity, in concentration units or psu
%          t: temperature, in deg C
%          p: pressure, in MPa
%  Output: y: kinematic viscosity in m^2/s
%
%	Check value	:	mu (.02891, 30) = 8.499e-4
%	check values	:	nu (0.035, 15.0, 0.0) = 1.19343e-6
%				nu (0.040, 40.0, 100.0) = 6.13058e-7

% M.Gregg: modified 21dec96

% make all inputs into column vectors and check lengths
s=s(:); t=t(:); p=p(:);
if (length(s)~= length(t)) || (length(t)~=length(p))
  disp('nu: inputs must be same length')
end

% check salinity magnitude
ig=find(~isnan(s));
s_avg=mean(s);
if s_avg > 1
  s=s/1000;
%   disp('nu: input s vector > 1, divided by 1000 to compute nu')
end
 
mu = (1.779e-3 - t.*(5.9319e-5 - t.*(1.2917e-6 - t * 1.3402e-8))) ...
      + s.*(2.8782e-3 - t.*(3.0553e-6 + t * 1.1835e-6));
rho = density(s,t,p);
y = mu./rho; 

end

function rho = density(s,t,p)
% density
%   Usage: rho = density(s,t,p)
%      s is salinity in concentration units
%      t is temperature in deg C
%      p is sea pressure in MPa
%      rho is density in kg/m^3
%
%
%	Coded from: Landolt-Bornstein,
%	"Numerical Data and Fundamental Relationships in Science
%	and Technology, New Series V/3a, Oceanography", pp 238-239.
%
%	Check Values:
%	Inputs:
%		S = .040 CU
%		T = 40 deg C
%		P = 100 MPa
%
%	Outputs:
%
%		Intermediate Values:
%
%		Bw   = -7.536450e-6 dbar^-1
%		A    =  3.446154
%		B    = -2.311202e-6 dbar^-1
%		A    =  3.465376
%		Kw   =  2.260440e5 dbar
%		K0   =  2.434421e5 dbar
%		K    =  2.778648e5 dbar
%
%		Return Values:
%
%		rhow = rho(0, T, 0)
%		rho0 = rho(S, T, 0)
%
%		rhow =  992.2204 kg m^-3
%		rho0 = 1021.6788 kg m^-3
%		rho  =  1059.8204 kg m^-3
%
% 
%------------------------------------------------------------------------------

	s = s*1000;	
	p = p*100;	

	t2 = t.*t;
	t3 = t2.*t;
	t4 = t3.*t;
	t5 = t4.*t;

	s2 = s.*s;
	s32 = s.*sqrt(s);
	rhow = 999.842594 + 6.793952e-2*t - 9.095290e-3*t2 +  ... 
	       1.001685e-4*t3 -1.120083e-6*t4 + 6.536332e-9*t5;

	A1 = 8.24493e-1 - 4.0899e-3*t + 7.6438e-5*t2 -8.2467e-7*t3 +  ...
	    5.3875e-9*t4;
	rho0 = rhow+(A1.*s) ...
	       + (-5.72466e-3 + 1.0227e-4*t - 1.6546e-6*t2).*s32 ...
	       + 4.8314e-4*s2;

	if (p < 0) 
	   !echo out of sea water
	   return
	end
	Kw = 1.965221e5 + 1484.206*t - 23.27105*t2 + 1.360477e-1*t3 ...
	     -5.155288e-4*t4;
	Aw = 3.239908 + 1.43713e-3*t + 1.16092e-4*t2 - 5.77905e-7*t3;
	Bw = 8.50935e-6 - 6.12293e-7*t + 5.2787e-9*t2;

	K0 = Kw+(546.746 - 6.03459*t + 1.09987e-1*t2 - 6.1670e-4*t3).*s ...
	        + (7.944e-1 + 1.6483e-1*t - 5.3009e-3*t2).*s32;
	A = Aw + (2.2838e-3 - 1.0981e-5*t -1.6078e-6*t2).*s ...
	       + 1.91075e-4*s32;
	B = Bw + (-9.9348e-8 +  2.0816e-9*t + 9.1697e-11*t2).*s;

	K = K0 + A.*p + B.*(p.*p);

	rho = rho0./(1.0-p./K);
end

function qc=get_qc_flag(epsilon,fom,fom_limit,despike,despike_limit,varepsilon,diss_ratio_limit,var_resolve,var_resolve_limit)
% 1	    Bit 0	if FOM > FOM_limit	2	0	0
% 2	    Bit 1	if despike_fraction > despike_fraction_limit	40%	0	0
% 4	    Bit 2	log(e_max)-log(e_min)|> diss_ratio_limit X \sigma_{\ln\varepsilon}	N/A	1	4
% 8	    Bit 3	if despike_iterations > despike_iterations_limit	To be confirmed	0	0
% 16	Bit 4	if variance resolved less than a threshold	50%	1	16
% 32	Bit 5	manual flag to be defined by user	N/A	0	0
% 64	Bit 6	manual flag to be defined by user	N/A	0	0
% 128	Bit 7	manual flag to be defined by user	N/A	0	0

qc=epsilon.*0; % good data
Depsi=(log(max(epsilon))-log(min(epsilon)));


fom_bit               = 0;
despike_bit           = 1;
varepsilon_bit        = 2;
% despike_iteration_bit = 3;
var_resol_bit         = 4;

for n=1:length(epsilon)
    if (fom(n)>fom_limit)
        qc(n)=qc(n) + bitshift(1,fom_bit);
    end
    if (despike(n)>despike_limit)
        qc(n)=qc(n) + bitshift(1,despike_bit);
    end
    if (Depsi>diss_ratio_limit*varepsilon)
        %     qc(1)=qc(1) + bitshift(1,varepsilon_bit);
        %     qc(1)=qc(1) + bitshift(1,varepsilon_bit);
    end
    if (var_resolve(n)>var_resolve_limit)
        qc(n)=qc(n)+  bitshift(1,var_resol_bit);
    end

end
end


function [epsilon_final, qc_final, epsi_std]=get_final_epsilon(epsilon,qc)
epsi_std=epsilon.*nan;
if length(unique(qc))==1
    %good epsilons we use the mean of both epsilons
    epsilon_final = mean(epsilon,'omitnan');
    qc_final=qc(1);    
else
    %otherwise the estimate with the small QC
    epsilon_final = epsilon(qc==min(qc));
    qc_final=min(qc);
end
end

