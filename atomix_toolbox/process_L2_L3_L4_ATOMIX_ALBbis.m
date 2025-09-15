function [FieldData]=process_L2_L3_L4_ATOMIX_ALBbis(Meta,GroupMeta,FieldData)
% Process levels L2 L3 and L4 using FieldData.L1
% using my (Le Boyer) epsilon processing
% The variable names inside L2, L3, L4 are defined in GroupMeta.
% GroupMeta is defined by shear_atomix_metada.yml
%
% 1/  Compute PSPD_REL (dP/dt) and smooth it with a 1s lowpass. 
% 2/  Select section with a speed cut out (0.2 m/s)
% 3/  Remove outliers within a sliding window of 10x 50ms
% 4/  Start a dor loop over a section and split the section into disslength
%     timeseries
% 5/  Compute the spectrum of each channel with pwelch(detrend(data),nfft,[],nfft,FS);
% 6/  Compute coherence spectrum between shear channels and acceleration
%     channels
% 7/  Correct shear spectra with a3 coherence spectrum. (a3 is aligned with the shear probes)
% 8/  Correct all spectra with the 320 Hz sin4 sampling filter
% 9/  Correct shear spectra with Oakey and 1st difference operator transfer function
% 10/ - if shear spectrum is not too high (function eps1_mmp):
%       3 iteration of epsilon estimates:
%           1st estimate eps1 is from int(shear)
%           between 2cpm and 10cpm
%           2nd estimated between 2cpm and K_kolmogorov (from eps1)
%           3rd estimated between 2cpm and K_kolmogorov (from eps2)
%      - if shear spectrum is high (function eps1_mmp): 
%           spectral fitting (polyfit) between 2 and 10 cpm 
%           
% 11/   KMAX is  K_kolmogorov (from eps3; eps1_mmp)
% 12/   add missing variance (epsilon2_correct)
%
%
% TODO: Figure out a way to desing what variable is used for acceleration correction

Levelnames=fieldnames(GroupMeta);
% 
% if isfield(FieldData.(Levelnames{1}),'TIME_CTD')
%     L_CTD=length(FieldData.(Levelnames{1}).TIME_CTD);
%     % name of the vibration acceleration field
%     F1='ACC';F1size=3;
% else
%     L_CTD=length(FieldData.(Levelnames{1}).TIME_SLOW);
%     % name of the vibration acceleration field
%     F1='VIB';F1size=2;
% end

if isfield(FieldData.(Levelnames{2}),'ACC')
%     L_CTD=length(FieldData.(Levelnames{1}).TIME_CTD);
    % name of the vibration acceleration field
    ACCsize=size(FieldData.(Levelnames{1}).ACC);
    F1='ACC';F1size=ACCsize(2);
elseif isfield(FieldData.(Levelnames{2}),'VIB')
%     L_CTD=length(FieldData.(Levelnames{1}).TIME_SLOW);
    % name of the vibration acceleration field
    F1='VIB';F1size=2;
else
    F1='';F1size=0;
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
% equivalent_dof= 8/3 * disslength./ nfft ./ 2
equivalent_dof= 2*dof;

[~, f]=pwelch(1:nfft,nfft,[],nfft,FS);


%% split in profiles
speed_limit=Meta.speed_cutout;
% find FLOWPS with the right sign and the right amplitude
time_slow_flag=0;
time_ctd_name='TIME_CTD';
if isfield(FieldData.(Levelnames{1}),'PRES_SLOW')
    time_slow_flag=1;
    time_ctd_name='TIME_SLOW';
end


%% get Transfer functions
if ~time_slow_flag
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
else
    H.samplingfilter=f.*0+1;
    H.shear=f.*0+1;
end

CTD_disslength=floor(disslength*FS_CTD/FS);
CTD_disslength=CTD_disslength-mod(CTD_disslength,2);
%% Process profile
n_disslength=0;

% for i=CTD_disslength*disslength_overlap+1:floor(CTD_disslength.*disslength_overlap):L_CTD-CTD_disslength*disslength_overlap
for tt=1:length(FieldData.PI_TIME)

    i    = find(FieldData.PI_TIME(tt)>FieldData.L1_converted.TIME_SLOW,1,'last');
    % get the ctd timestamp at the center of our disslength
    idx_ctd    = i+ (-CTD_disslength/2:CTD_disslength/2);
    % move the idx_ctd back if idx_ctd go over the length
    idx_ctd=idx_ctd-max(0,idx_ctd(end)-length(FieldData.L1_converted.TEMP));


    i    = find(FieldData.PI_TIME(tt)>FieldData.L2_cleaned.TIME,1,'last');
    idx_shear  = i-1+(-disslength/2:disslength/2-1);
    idx_shear  = idx_shear-max(0,idx_shear(end)-length(FieldData.L1_converted.SHEAR));

    %     fprintf("disslength=%f sec\r\n",disslength./Meta.fs_fast)

    if length(FieldData.(Levelnames{2}).PSPD_REL)==length(FieldData.(Levelnames{2}).SHEAR)
        idx_PSD_REL=idx_shear;
    else
        idx_PSD_REL=idx_ctd;
    end
    if(mean(FieldData.(Levelnames{2}).PSPD_REL(idx_PSD_REL))>speed_limit) && ...
            sum(isnan(FieldData.(Levelnames{2}).SHEAR(idx_shear,1)))==0   && ...
            sum(isnan(FieldData.(Levelnames{2}).SHEAR(idx_shear,2)))==0

        n_disslength=n_disslength+1;
        p=n_disslength;
        % evaluate despiking 
        %TODO add tilt computation and flag here
        despike1=FieldData.(Levelnames{2}).SHEAR(idx_shear,1)-FieldData.(Levelnames{1}).SHEAR(idx_shear,1);
        despike1=sum(abs(despike1)>0)./length(despike1); % percent of corrected samples 
        despike2=FieldData.(Levelnames{2}).SHEAR(idx_shear,2)-FieldData.(Levelnames{1}).SHEAR(idx_shear,2);
        despike2=sum(abs(despike2)>0)./length(despike2); % percent of corrected samples 


        % select the disslength shear timeseries
        s1 = FieldData.(Levelnames{2}).SHEAR(idx_shear,1);
        s2 = FieldData.(Levelnames{2}).SHEAR(idx_shear,2);
        switch F1
            case "ACC"
                a1 = FieldData.(Levelnames{2}).(F1)(idx_shear,1);
                a2 = FieldData.(Levelnames{2}).(F1)(idx_shear,2);
                a3 = FieldData.(Levelnames{2}).(F1)(idx_shear,3);
            case "VIB"
                a1 = FieldData.(Levelnames{2}).(F1)(idx_shear,1);
                a2 = FieldData.(Levelnames{2}).(F1)(idx_shear,2);
        end
        % compute the average variable of disslength
        if length(FieldData.(Levelnames{2}).PSPD_REL)==length(FieldData.(Levelnames{2}).SHEAR)
            w     = mean(FieldData.(Levelnames{2}).PSPD_REL(idx_shear));
        else
            w     = mean(FieldData.(Levelnames{2}).PSPD_REL(idx_ctd));
        end
        if length(FieldData.(Levelnames{1}).TEMP)==length(FieldData.(Levelnames{2}).SHEAR)
            t     = mean(FieldData.(Levelnames{1}).TEMP(idx_shear));
        else
            t     = mean(FieldData.(Levelnames{1}).TEMP(idx_ctd));
        end
        if length(FieldData.(Levelnames{1}).PRES)==length(FieldData.(Levelnames{2}).SHEAR)
            pres  = mean(FieldData.(Levelnames{1}).PRES(idx_shear));
        else
            pres  = mean(FieldData.(Levelnames{1}).PRES(idx_ctd));
        end

        if ~time_slow_flag
        if length(FieldData.(Levelnames{1}).PSAL)==length(FieldData.(Levelnames{2}).SHEAR)
            s = mean(FieldData.(Levelnames{1}).PSAL(idx_shear));
        else
            s = mean(FieldData.(Levelnames{1}).PSAL(idx_ctd));
        end
        else
            s = 31;
        end

        % compute the frequency spectra
        [Ps1_f, f] = pwelch(detrend(s1),nfft,[],nfft,FS);
        [Ps2_f, ~] = pwelch(detrend(s2),nfft,[],nfft,FS);
        [Pa1_f, ~] = pwelch(detrend(a1),nfft,[],nfft,FS);
        [Pa2_f, ~] = pwelch(detrend(a2),nfft,[],nfft,FS);
        if ~time_slow_flag
            [Pa3_f, ~] = pwelch(detrend(a3),nfft,[],nfft,FS);
        end

        %compute coherence spectra (shear probes are aligned with a3)
        if ~time_slow_flag
            [Cu1a,~] = mscohere(detrend(s1),detrend(a3),nfft,[],nfft,FS);
            [Cu2a,~] = mscohere(detrend(s2),detrend(a3),nfft,[],nfft,FS);
        else
            [Cu1a,~] = mscohere(detrend(s1),detrend(a1),nfft,[],nfft,FS);
            [Cu2a,~] = mscohere(detrend(s2),detrend(a2),nfft,[],nfft,FS);
        end
        % correct shear spectra with coherence
        Ps1_clean_f = Ps1_f.*(1-Cu1a);
        Ps2_clean_f = Ps2_f.*(1-Cu2a);

        % compute the transfer function of the mechanical probe response
        TFoakey = haf_oakey(f,w);

        % compute wavenumber axis

        % Acceleration convert frequency spectra to wavenumber spectra
        AP_SPEC1 = Pa1_f .* w ./  H.samplingfilter;
        AP_SPEC2 = Pa2_f .* w ./  H.samplingfilter;
        if ~time_slow_flag
            AP_SPEC3 = Pa3_f .* w ./  H.samplingfilter;
        end

        % shear convert frequency spectra to wavenumber spectra
        goodman_correction=1-(1.02/9);
        SH_SPEC1=Ps1_f .* w ./ TFoakey./H.shear/goodman_correction;
        SH_SPEC2=Ps2_f .* w ./ TFoakey./H.shear/goodman_correction;

        % shear clean convert frequency spectra to wavenumber spectra
        SH_SPEC_CLEAN1=Ps1_clean_f .* w ./ TFoakey./H.shear/goodman_correction;
        SH_SPEC_CLEAN2=Ps2_clean_f .* w ./ TFoakey./H.shear/goodman_correction;


        % compute kvis 
        kvis = nu(s,t,db2MPa(pres));
        % compute epsilon and cutoff wavenumber
        
        Ps1sheark = squeeze(SH_SPEC_CLEAN1);
        Ps2sheark = squeeze(SH_SPEC_CLEAN2);
        k         = f./w;
        kmax      = FPUMP./w;

        [epsilon1,kmin1,kc1,method1] = eps1_mmp(k,Ps1sheark,kvis,kmax);
        [epsilon2,kmin2,kc2,method2] = eps1_mmp(k,Ps2sheark,kvis,kmax);


        %compute  theoretical spectra panchev
        [fom1,var_resolve1] = get_panchev_fom(epsilon1,Ps1sheark,k,kvis,kc1,equivalent_dof);
        [fom2,var_resolve2] = get_panchev_fom(epsilon2,Ps2sheark,k,kvis,kc2,equivalent_dof);

        % qc flag
        % get_qc_flag(epsilon,fom,fom_limit,despike,despike_limit,varepsilon,diss_ratio_limit,var_resolve,var_resolve_limit)
        fom_limit=3;
        despike_limit=.4; %40%
        var_resolve_limit=.5; %50%
        varepsilon=nan; %TODO deal with varepsilon and epsi_std. Right now I do not know how.
        diss_ratio_limit=nan;

        qc=get_qc_flag([epsilon1 epsilon2], ...
                       [fom1,fom2],fom_limit,...
                       [despike1,despike2],despike_limit,...
                       varepsilon,...
                       diss_ratio_limit,...
                       [var_resolve1,var_resolve2],var_resolve_limit);

        % final epsilon product
        [epsi_final, ~,epsi_std]=get_final_epsilon(epsilon1,epsilon2,qc);


        variable_names=fieldnames(GroupMeta.(Levelnames{3}).variables);
        for n=1:length(variable_names)
            wh_var=variable_names{n};
            switch wh_var
                case "dimension"
                case "TIME"
                    FieldData.(Levelnames{3}).TIME(p)    = mean(FieldData.(Levelnames{1}).TIME(idx_shear));
                case {"FLOWPS","PSPD_REL"}
                    FieldData.(Levelnames{3}).PSPD_REL(p)  = w;
                case "PRES"
                    FieldData.(Levelnames{3}).PRES(p)    = pres;
                case "TEMP"
                    FieldData.(Levelnames{3}).TEMP(p)    = t;
                case "PSAL"
                    FieldData.(Levelnames{3}).PSAL(p)     = s;
                case "FREQ"
                    FieldData.(Levelnames{3}).FREQ(:)    = f;
                case "KCYC"
                    FieldData.(Levelnames{3}).KCYC(p,:)  = f./w;
                case "ACC_SPEC"
                    FieldData.(Levelnames{3}).ACC_SPEC(p,:,1)=AP_SPEC1;
                    FieldData.(Levelnames{3}).ACC_SPEC(p,:,2)=AP_SPEC2;
                    FieldData.(Levelnames{3}).ACC_SPEC(p,:,3)=AP_SPEC3;
                case "VIB_SPEC"
                    FieldData.(Levelnames{3}).VIB_SPEC(p,:,1)=AP_SPEC1;
                    FieldData.(Levelnames{3}).VIB_SPEC(p,:,2)=AP_SPEC2;
                case "SH_SPEC"
                    FieldData.(Levelnames{3}).SH_SPEC(p,:,1)=SH_SPEC1;
                    FieldData.(Levelnames{3}).SH_SPEC(p,:,2)=SH_SPEC2;
                case "SH_SPEC_CLEAN"
                    FieldData.(Levelnames{3}).SH_SPEC_CLEAN(p,:,1)=SH_SPEC_CLEAN1;
                    FieldData.(Levelnames{3}).SH_SPEC_CLEAN(p,:,2)=SH_SPEC_CLEAN2;
                case "CA_TF"
                    FieldData.(Levelnames{3}).CA_TF=H.shear;
                case "SECTION_NUMBER"
                    FieldData.(Levelnames{3}).SECTION_NUMBER(p)=mean(FieldData.(Levelnames{2}).SECTION_NUMBER(idx_shear));
                case "DOF"
                    FieldData.(Levelnames{3}).DOF(p)=equivalent_dof;
            end
        end

        %L4 level
        variable_names=fieldnames(GroupMeta.(Levelnames{4}).variables);
        for n=1:length(variable_names)
            wh_var=variable_names{n};
            switch wh_var
                case "dimension"
                case "TIME"
                    FieldData.(Levelnames{4}).TIME(p)   = FieldData.(Levelnames{3}).TIME(p);
                case {"FLOWPS","PSPD_REL"}
                    FieldData.(Levelnames{4}).PSPD_REL(p) = FieldData.(Levelnames{3}).PSPD_REL(p);
                case "PRES"
                    FieldData.(Levelnames{4}).PRES(p)   = FieldData.(Levelnames{3}).PRES(p);
                case "TEMP"
                    FieldData.(Levelnames{4}).TEMP(p)   = FieldData.(Levelnames{3}).TEMP(p);
                case "PSAL"
                    FieldData.(Levelnames{4}).PSAL(p)    = FieldData.(Levelnames{3}).PSAL(p);
                case "KVISC"
                    FieldData.(Levelnames{4}).KVISC(p)  = kvis;
                case "KMAX"
                    FieldData.(Levelnames{4}).KMAX(p,:)   = [kc1 kc2];
                case "KMIN"
                    FieldData.(Levelnames{4}).KMIN(p,:)   = [kmin1 kmin2];
                case "SECTION_NUMBER"
                    FieldData.(Levelnames{4}).SECTION_NUMBER(p) = FieldData.(Levelnames{3}).SECTION_NUMBER(p);
                case "FOM"
                    FieldData.(Levelnames{4}).FOM(p,:)     = [fom1 fom2];
                case "MAD"
                    FieldData.(Levelnames{4}).MAD(p,:)     = [nan nan];
                case "METHOD"
                    FieldData.(Levelnames{4}).METHOD(p,:)  = [method1 method2];
                case "VAR_RESOLVED"
                    FieldData.(Levelnames{4}).VAR_RESOLVED(p,:)  = [var_resolve1 var_resolve2];
                case "EPSI"
                    FieldData.(Levelnames{4}).EPSI(p,:)    = [epsilon1 epsilon2];
                case {"EPSI_QC","EPSI_FLAGS"}
                    FieldData.(Levelnames{4}).EPSI_FLAGS(p,:) = qc;
                case "EPSI_FINAL"
                    FieldData.(Levelnames{4}).EPSI_FINAL(p)    = epsi_final;
                case "EPSI_STD"
                    FieldData.(Levelnames{4}).EPSI_STD(p)    = epsi_std;
            end
        end

    end % end of (idx_shear_center-disslength/2>0)
end % end for p=1:P
% transpose if needed
for l=1:length(Levelnames)
    variable_names=fieldnames(GroupMeta.(Levelnames{l}).variables);
    if isfield(FieldData,Levelnames{l})
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
    else
        fprintf('Warning - no field %s\r\n',Levelnames{l})

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

kmin=2;
dKI=0.2;
KI=(2:dKI:10); % wavenumber array for interpolation

dk=nanmean(diff(k));
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
P_interpolated=interp1(k(krange),Psheark(krange),KI);
% ALB change to nansum since coherence correction can introduces nansclear 
% shear10=nansum(P_interpolated)*0.2;
shear10=nansum(P_interpolated)*dKI;
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
kc2 = 0.0816*( eps1  / kvis^3 )^(1/4);  % k for 90% variance of Panchev spectrum
if kc2>kmax
	kc2=kmax; % limit set by noise spectrum
end
krange=find(k>=kmin & k<=kc2);
% eps2=7.5*kvis*nansum(Psheark(krange))*dk/.9; % .9 we want to get 90% of the shear variance
eps2=7.5*kvis*nansum(Psheark(krange))*dk; %

% third estimate: same as before.
kc=0.0816*( eps2 / kvis^3 )^(1/4);
if kc > kmax
	kc=kmax;
end
krange=find(k>=kmin & k<=kc);
eps3 = 7.5*kvis*nansum(Psheark(krange))*dk; 

if eps3<1e-11 || length(krange) < 4
	epsilon=1e-11;
else
  mf=epsilon2_correct(eps3,kvis);
	epsilon=mf*eps3;
end

%selecting the closest k for gc;
kc=k(find(k>kc,1,'first'));

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

% kmin=2;
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
  mf=epsilon2_correct(eps3,kvis);
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

function [fom, var_resolve] = get_panchev_fom(epsilon,Psheark,k,kvis,kc,dof)
sig_lnS=5/4*dof^(-7/9);
[kpan,Ppan] = panchev(epsilon,kvis);
Ppan=interp1(kpan(~isnan(Ppan)),Ppan(~isnan(Ppan)),k);

fom=log(Psheark./Ppan);
fom=nanvar(fom(k>2 & k<kc))./sig_lnS;

var_resolve=100.* sum(Psheark).*nanmean(diff(k))  ./ ...
                 (sum(Ppan).*nanmean(diff(kpan)));

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
	k0=1/(1000*eta);
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

qc=[0 0]; % good data
epsilon1=epsilon(1);
epsilon2=epsilon(2);
fom1=fom(1);
fom2=fom(2);
depsike1=despike(1);
depsike2=despike(2);
var_resolve1=var_resolve(1);
var_resolve2=var_resolve(2);
Depsi=(log(max(epsilon))-log(min(epsilon)));


fom_bit               = 0;
despike_bit           = 1;
varepsilon_bit        = 2;
% despike_iteration_bit = 3;
var_resol_bit         = 4;

if (fom1>fom_limit)
    qc(1)=qc(1) + bitshift(1,fom_bit);
end
if (fom2>fom_limit)
    qc(2)=qc(2)+  bitshift(1,fom_bit);
end

if (depsike1>despike_limit)
    qc(1)=qc(1) + bitshift(1,despike_bit);
end
if (depsike2>despike_limit)
    qc(2)=qc(2)+  bitshift(1,despike_bit);
end

if (Depsi>diss_ratio_limit*varepsilon)
%     qc(1)=qc(1) + bitshift(1,varepsilon_bit);
%     qc(1)=qc(1) + bitshift(1,varepsilon_bit);
end

if (var_resolve1>var_resolve_limit)
    qc(1)=qc(1)+  bitshift(1,var_resol_bit);
end
if (var_resolve2>var_resolve_limit)
    qc(2)=qc(2)+  bitshift(1,var_resol_bit);
end
    

end


function [epsilon_final, qc_final, epsi_std]=get_final_epsilon(epsilon1,epsilon2,qc)
epsi_std=nan;
if qc(1)==qc(2)
    %good epsilons we use the mean of both epsilons
    epsilons=[epsilon1 epsilon2];
    epsilon_final = nanmean(epsilons);
    qc_final=1;    
else
    %otherwise the estimate with the small QC
    epsilons=[epsilon1 epsilon2];
    epsilon_final = epsilons(qc==min(qc));
    qc_final=min(qc);
end
end


