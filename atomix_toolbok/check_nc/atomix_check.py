"""
Python software to check Atomix netCDF files for consistency

author(s): Peter Holtermann (peter.holtermann@io-warnemuende.de)
version: 0.1

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or any
later version.

"""

import logging
import sys
import os
import netCDF4

# Setup logging module
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
logger = logging.getLogger('atomix')

# Structure from Atomix webpage https://wiki.app.uib.no/atomix/index.php (version from the 13.09.2023)
atomix_structure = {}
atomix_structure['L1_converted']   = {'DIMENSIONS':{},'VARIABLES':{}}
atomix_structure['L2_cleaned']     = {'DIMENSIONS':{},'VARIABLES':{}}
atomix_structure['L3_spectra']     = {'DIMENSIONS':{},'VARIABLES':{}}
atomix_structure['L4_dissipation'] = {'DIMENSIONS':{},'VARIABLES':{}}

# Level1 data
levelname = 'L1_converted'
atomix_structure['L1_converted']['DIMENSIONS']['TIME'] = {'description':'length of the record from turbulence (fast) data channels'}
atomix_structure['L1_converted']['DIMENSIONS']['TIME_XYZ'] = {'description':'length of the record from sensor XYZ'}
atomix_structure['L1_converted']['DIMENSIONS']['N_SHEAR_SENSORS'] = {'description':'number of shear channels (shear sensors)' }
atomix_structure['L1_converted']['DIMENSIONS']['N_ACC_SENSORS']   = {'description':'number of acc channels (acc sensors)'}
atomix_structure['L1_converted']['DIMENSIONS']['N_VIB_SENSORS']   = {'description':'number of vib channels (vib sensors)'}
atomix_structure['L1_converted']['DIMENSIONS']['N_GRADT_SENSORS']   = {'description':'number of channels/sensors of type GRADT'}
atomix_structure['L1_converted']['DIMENSIONS']['N_GRADC_SENSORS']   = {'description':'number of channels/sensors of type GRADC'}
atomix_structure['L1_converted']['DIMENSIONS']['N_CNDC_SENSORS']   = {'description':'number of channels/sensors of type C'}
atomix_structure['L1_converted']['DIMENSIONS']['N_XYZ_SENSORS']   = {'description':'number of channels/sensors of type XYZ'}


atomix_structure['L1_converted']['VARIABLES']['TIME'] = {'required':'required'}
atomix_structure['L1_converted']['VARIABLES']['TIME']['standard_name'] = 'time'
atomix_structure['L1_converted']['VARIABLES']['TIME']['units'] = 'CF-Conventions compatible offset'
atomix_structure['L1_converted']['VARIABLES']['TIME']['DIMENSIONS'] = [atomix_structure['L1_converted']['DIMENSIONS']['TIME']]

atomix_structure['L1_converted']['VARIABLES']['SHEAR'] = {'required':'required'}
atomix_structure['L1_converted']['VARIABLES']['SHEAR']['standard_name'] = 'sea_water_velocity_shear'
atomix_structure['L1_converted']['VARIABLES']['SHEAR']['units'] = 's-1'
atomix_structure['L1_converted']['VARIABLES']['SHEAR']['DIMENSIONS'] = [atomix_structure['L1_converted']['DIMENSIONS']['TIME'],atomix_structure['L1_converted']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L1_converted']['VARIABLES']['PSPD_REL'] = {'required':'highly-recommended'}
atomix_structure['L1_converted']['VARIABLES']['PSPD_REL']['standard_name'] = 'platform_speed_wrt_sea_water'
atomix_structure['L1_converted']['VARIABLES']['PSPD_REL']['units'] ='m s-1'
atomix_structure['L1_converted']['VARIABLES']['PSPD_REL']['DIMENSIONS'] = [atomix_structure['L1_converted']['DIMENSIONS']['TIME']]

atomix_structure['L1_converted']['VARIABLES']['ACC'] = {'required':'highly-recommended'}
atomix_structure['L1_converted']['VARIABLES']['ACC']['standard_name'] = 'platform_acceleration'
atomix_structure['L1_converted']['VARIABLES']['ACC']['units'] ='m s-2'
atomix_structure['L1_converted']['VARIABLES']['ACC']['DIMENSIONS'] = [atomix_structure['L1_converted']['DIMENSIONS']['TIME'],atomix_structure['L1_converted']['DIMENSIONS']['N_ACC_SENSORS']]

atomix_structure['L1_converted']['VARIABLES']['VIB'] = {'required':'highly-recommended'}
atomix_structure['L1_converted']['VARIABLES']['VIB']['standard_name'] = 'platform_vibration'
atomix_structure['L1_converted']['VARIABLES']['VIB']['units'] =''
atomix_structure['L1_converted']['VARIABLES']['VIB']['DIMENSIONS'] = [atomix_structure['L1_converted']['DIMENSIONS']['TIME'],atomix_structure['L1_converted']['DIMENSIONS']['N_VIB_SENSORS']]

atomix_structure['L1_converted']['VARIABLES']['PRES'] = {'required':'highly-recommended'}
atomix_structure['L1_converted']['VARIABLES']['PRES']['standard_name'] = 'water_pressure'
atomix_structure['L1_converted']['VARIABLES']['PRES']['units'] = 'dbar'
atomix_structure['L1_converted']['VARIABLES']['PRES']['DIMENSIONS'] = [atomix_structure['L1_converted']['DIMENSIONS']['TIME']]

atomix_structure['L1_converted']['VARIABLES']['TEMP'] = {'required':'highly-recommended'}
atomix_structure['L1_converted']['VARIABLES']['TEMP']['standard_name'] = 'water_pressure'
atomix_structure['L1_converted']['VARIABLES']['TEMP']['units'] = 'degree_Celsius'
atomix_structure['L1_converted']['VARIABLES']['TEMP']['DIMENSIONS'] = [atomix_structure['L1_converted']['DIMENSIONS']['TIME']]
varname = 'ACC'
atomix_structure[levelname]['VARIABLES'][varname] = {'required':'optional'}
atomix_structure[levelname]['VARIABLES'][varname]['standard_name'] = 'platform_acceleration'
atomix_structure[levelname]['VARIABLES'][varname]['units'] ='m s-2'
atomix_structure[levelname]['VARIABLES'][varname]['DIMENSIONS'] = [atomix_structure[levelname]['DIMENSIONS']['TIME'],atomix_structure[levelname]['DIMENSIONS']['N_ACC_SENSORS']]
varname = 'VIB'
atomix_structure[levelname]['VARIABLES'][varname] = {'required':'optional'}
atomix_structure[levelname]['VARIABLES'][varname]['standard_name'] = 'platform_vibration'
atomix_structure[levelname]['VARIABLES'][varname]['units'] =''
atomix_structure[levelname]['VARIABLES'][varname]['DIMENSIONS'] = [atomix_structure[levelname]['DIMENSIONS']['TIME'],atomix_structure[levelname]['DIMENSIONS']['N_VIB_SENSORS']]
varname = 'ROLL'
atomix_structure[levelname]['VARIABLES'][varname] = {'required':'optional'}
atomix_structure[levelname]['VARIABLES'][varname]['standard_name'] = 'platform_roll_angle'
atomix_structure[levelname]['VARIABLES'][varname]['units'] ='degree'
atomix_structure[levelname]['VARIABLES'][varname]['DIMENSIONS'] = [atomix_structure[levelname]['DIMENSIONS']['TIME']]
varname = 'PITCH'
atomix_structure[levelname]['VARIABLES'][varname] = {'required':'optional'}
atomix_structure[levelname]['VARIABLES'][varname]['standard_name'] = 'platform_pitch_angle'
atomix_structure[levelname]['VARIABLES'][varname]['units'] ='degree'
atomix_structure[levelname]['VARIABLES'][varname]['DIMENSIONS'] = [atomix_structure[levelname]['DIMENSIONS']['TIME']]
varname = 'GRADT'
atomix_structure[levelname]['VARIABLES'][varname] = {'required':'optional'}
atomix_structure[levelname]['VARIABLES'][varname]['standard_name'] = 'derivative_of_seawater_temperature_wrt_X'
atomix_structure[levelname]['VARIABLES'][varname]['units'] ='degree_Celcius m-1'
atomix_structure[levelname]['VARIABLES'][varname]['DIMENSIONS'] = [atomix_structure[levelname]['DIMENSIONS']['TIME'],atomix_structure[levelname]['DIMENSIONS']['TIME'],atomix_structure[levelname]['DIMENSIONS']['N_GRADT_SENSORS']]
varname = 'GRADC'
atomix_structure[levelname]['VARIABLES'][varname] = {'required':'optional'}
atomix_structure[levelname]['VARIABLES'][varname]['standard_name'] = 'derivative_of_seawater_conductivity_wrt_X'
atomix_structure[levelname]['VARIABLES'][varname]['units'] ='Typically not calibrated'
atomix_structure[levelname]['VARIABLES'][varname]['DIMENSIONS'] = [atomix_structure[levelname]['DIMENSIONS']['TIME'],atomix_structure[levelname]['DIMENSIONS']['TIME'],atomix_structure[levelname]['DIMENSIONS']['N_GRADC_SENSORS']]
varname = 'CNDC'
atomix_structure[levelname]['VARIABLES'][varname] = {'required':'optional'}
atomix_structure[levelname]['VARIABLES'][varname]['standard_name'] = 'sea_water_electrical_conductivity (water_electrical_conductivity)'
atomix_structure[levelname]['VARIABLES'][varname]['units'] ='S m-1'
atomix_structure[levelname]['VARIABLES'][varname]['DIMENSIONS'] = [atomix_structure[levelname]['DIMENSIONS']['TIME'],atomix_structure[levelname]['DIMENSIONS']['TIME'],atomix_structure[levelname]['DIMENSIONS']['N_CNDC_SENSORS']]


# Level2 data
atomix_structure['L2_cleaned']['DIMENSIONS']['TIME'] = {'description':'length of the record from turbulence (fast) data channels'}
atomix_structure['L2_cleaned']['DIMENSIONS']['N_SHEAR_SENSORS'] = {'description':'number of shear channels (shear sensors)' }
atomix_structure['L2_cleaned']['DIMENSIONS']['N_ACC_SENSORS']   = {'description':'number of acc channels (acc sensors)'}
atomix_structure['L2_cleaned']['DIMENSIONS']['N_VIB_SENSORS']   = {'description':'number of vib channels (vib sensors)'}
atomix_structure['L2_cleaned']['DIMENSIONS']['N_XYZ_SENSORS']   = {'description':'number of other channels/sensors'}

atomix_structure['L2_cleaned']['VARIABLES']['TIME'] = {'required':'required'}
atomix_structure['L2_cleaned']['VARIABLES']['TIME']['standard_name'] = 'time'
atomix_structure['L2_cleaned']['VARIABLES']['TIME']['units'] = 'CF-Conventions compatible offset'
atomix_structure['L2_cleaned']['VARIABLES']['TIME']['DIMENSIONS'] = [atomix_structure['L2_cleaned']['DIMENSIONS']['TIME']]

atomix_structure['L2_cleaned']['VARIABLES']['SECTION_NUMBER'] = {'required':'required'}
atomix_structure['L2_cleaned']['VARIABLES']['SECTION_NUMBER']['standard_name'] = 'unique_identifier_for_each_section_of_data_from_timeseries'
atomix_structure['L2_cleaned']['VARIABLES']['SECTION_NUMBER']['units'] = ''
atomix_structure['L2_cleaned']['VARIABLES']['SECTION_NUMBER']['DIMENSIONS'] = [atomix_structure['L2_cleaned']['DIMENSIONS']['TIME']]

atomix_structure['L2_cleaned']['VARIABLES']['SHEAR'] = {'required':'required'}
atomix_structure['L2_cleaned']['VARIABLES']['SHEAR']['standard_name'] = 'sea_water_velocity_shear'
atomix_structure['L2_cleaned']['VARIABLES']['SHEAR']['units'] = 's-1'
atomix_structure['L2_cleaned']['VARIABLES']['SHEAR']['DIMENSIONS'] = [atomix_structure['L2_cleaned']['DIMENSIONS']['TIME'],atomix_structure['L2_cleaned']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L2_cleaned']['VARIABLES']['PSPD_REL'] = {'required':'required'}
atomix_structure['L2_cleaned']['VARIABLES']['PSPD_REL']['standard_name'] = 'platform_speed_wrt_sea_water'
atomix_structure['L2_cleaned']['VARIABLES']['PSPD_REL']['units'] ='m s-1'
atomix_structure['L2_cleaned']['VARIABLES']['PSPD_REL']['DIMENSIONS'] = [atomix_structure['L2_cleaned']['DIMENSIONS']['TIME']]

atomix_structure['L2_cleaned']['VARIABLES']['ACC'] = {'required':'optional'}
atomix_structure['L2_cleaned']['VARIABLES']['ACC']['standard_name'] = 'platform_acceleration'
atomix_structure['L2_cleaned']['VARIABLES']['ACC']['units'] ='m s-2'
atomix_structure['L2_cleaned']['VARIABLES']['ACC']['DIMENSIONS'] = [atomix_structure['L2_cleaned']['DIMENSIONS']['TIME'],atomix_structure['L2_cleaned']['DIMENSIONS']['N_ACC_SENSORS']]

atomix_structure['L2_cleaned']['VARIABLES']['VIB'] = {'required':'optional'}
atomix_structure['L2_cleaned']['VARIABLES']['VIB']['standard_name'] = 'platform_vibration'
atomix_structure['L2_cleaned']['VARIABLES']['VIB']['units'] =''
atomix_structure['L2_cleaned']['VARIABLES']['VIB']['DIMENSIONS'] = [atomix_structure['L2_cleaned']['DIMENSIONS']['TIME'],atomix_structure['L2_cleaned']['DIMENSIONS']['N_VIB_SENSORS']]



# Level3 data
atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA']    = {'description':'length of the record of average times of spectral segments. This also equals time of dissipation estimates. '}
atomix_structure['L3_spectra']['DIMENSIONS']['N_WAVENUMBER']    = {'description':'length of the wavenumber array' }
atomix_structure['L3_spectra']['DIMENSIONS']['N_SHEAR_SENSORS'] = {'description':'number of shear channels (shear sensors)' }
atomix_structure['L3_spectra']['DIMENSIONS']['N_ACC_SENSORS']   = {'description':'number of acc channels (acc sensors)'}
atomix_structure['L3_spectra']['DIMENSIONS']['N_VIB_SENSORS']   = {'description':'number of vib channels (vib sensors)'}
atomix_structure['L3_spectra']['DIMENSIONS']['N_XYZ_SENSORS']   = {'description':'number of other channels/sensors'}
atomix_structure['L3_spectra']['DIMENSIONS']['N_SH_ACC_SPEC']   = {'description':'number of shear and acc sensor combinations: N_SHEAR_SENSORS * N_ACC_SENSORS'}
atomix_structure['L3_spectra']['DIMENSIONS']['N_SH_VIB_SPEC']   = {'description':'number of shear and vib sensor combinations: N_SHEAR_SENSORS * N_VIB_SENSORS'}
atomix_structure['L3_spectra']['DIMENSIONS']['N_GLOBAL_VALUES'] = {'description':'Dimension for 1 data point (for the entire analysis)'}


atomix_structure['L3_spectra']['VARIABLES']['TIME'] = {'required':'required'}
atomix_structure['L3_spectra']['VARIABLES']['TIME']['standard_name'] = 'time'
atomix_structure['L3_spectra']['VARIABLES']['TIME']['units'] = 'CF-Conventions compatible offset'
atomix_structure['L3_spectra']['VARIABLES']['TIME']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L3_spectra']['VARIABLES']['SECTION_NUMBER'] = {'required':'required'}
atomix_structure['L3_spectra']['VARIABLES']['SECTION_NUMBER']['standard_name'] = 'unique_identifier_for_each_section_of_data_from_timeseries'
atomix_structure['L3_spectra']['VARIABLES']['SECTION_NUMBER']['units'] = ''
atomix_structure['L3_spectra']['VARIABLES']['SECTION_NUMBER']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L3_spectra']['VARIABLES']['PSPD_REL'] = {'required':'required'}
atomix_structure['L3_spectra']['VARIABLES']['PSPD_REL']['standard_name'] = 'platform_speed_wrt_sea_water'
atomix_structure['L3_spectra']['VARIABLES']['PSPD_REL']['units'] ='m s-1'
atomix_structure['L3_spectra']['VARIABLES']['PSPD_REL']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L3_spectra']['VARIABLES']['SH_SPEC'] = {'required':'required'}
atomix_structure['L3_spectra']['VARIABLES']['SH_SPEC']['standard_name'] = 'shear_probe_spectrum'
atomix_structure['L3_spectra']['VARIABLES']['SH_SPEC']['units'] = 's-2 cpm-1'
atomix_structure['L3_spectra']['VARIABLES']['SH_SPEC']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L3_spectra']['DIMENSIONS']['N_WAVENUMBER'],atomix_structure['L3_spectra']['DIMENSIONS']['N_SHEAR_SENSORS'] ]

atomix_structure['L3_spectra']['VARIABLES']['KCYC'] = {'required':'required'}
atomix_structure['L3_spectra']['VARIABLES']['KCYC']['standard_name'] = 'cyclic_wavenumber'
atomix_structure['L3_spectra']['VARIABLES']['KCYC']['units'] = 'cpm'
atomix_structure['L3_spectra']['VARIABLES']['KCYC']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L3_spectra']['DIMENSIONS']['N_WAVENUMBER'] ]

atomix_structure['L3_spectra']['VARIABLES']['DOF'] = {'required':'optional'}
atomix_structure['L3_spectra']['VARIABLES']['DOF']['standard_name'] = 'degrees_of_freedom_of_spectrum'
atomix_structure['L3_spectra']['VARIABLES']['DOF']['units'] = ''
atomix_structure['L3_spectra']['VARIABLES']['DOF']['DIMENSIONS'] = [1]

atomix_structure['L3_spectra']['VARIABLES']['SH_SPEC_CLEAN'] = {'required':'required'}
atomix_structure['L3_spectra']['VARIABLES']['SH_SPEC_CLEAN']['standard_name'] = 'shear_probe_spectrum_clean'
atomix_structure['L3_spectra']['VARIABLES']['SH_SPEC_CLEAN']['units'] = 's-2 cpm-1'
atomix_structure['L3_spectra']['VARIABLES']['SH_SPEC_CLEAN']['DIMENSIONS']  = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L3_spectra']['DIMENSIONS']['N_WAVENUMBER'],atomix_structure['L3_spectra']['DIMENSIONS']['N_SHEAR_SENSORS'] ]
atomix_structure['L3_spectra']['VARIABLES']['SH_SPEC_CLEAN']['comment'] = 'SH_SPEC_CLEAN must be included in the file if these spectra were used to compute EPS'

atomix_structure['L3_spectra']['VARIABLES']['PRES'] = {'required':'optional'}
atomix_structure['L3_spectra']['VARIABLES']['PRES']['standard_name'] = 'water_pressure'
atomix_structure['L3_spectra']['VARIABLES']['PRES']['units'] = 'dbar'
atomix_structure['L3_spectra']['VARIABLES']['PRES']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L3_spectra']['VARIABLES']['ACC_SPEC'] = {'required':'optional'}
atomix_structure['L3_spectra']['VARIABLES']['ACC_SPEC']['standard_name'] = 'acceleration_sensor_spectrum'
atomix_structure['L3_spectra']['VARIABLES']['ACC_SPEC']['units'] = 'm2 s-4 cpm-1'
atomix_structure['L3_spectra']['VARIABLES']['ACC_SPEC']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L3_spectra']['DIMENSIONS']['N_WAVENUMBER'],atomix_structure['L3_spectra']['DIMENSIONS']['N_ACC_SENSORS'] ]

atomix_structure['L3_spectra']['VARIABLES']['VIB_SPEC'] = {'required':'optional'}
atomix_structure['L3_spectra']['VARIABLES']['VIB_SPEC']['standard_name'] = 'vibration_sensor_spectrum'
atomix_structure['L3_spectra']['VARIABLES']['VIB_SPEC']['units'] = ''
atomix_structure['L3_spectra']['VARIABLES']['VIB_SPEC']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L3_spectra']['DIMENSIONS']['N_WAVENUMBER'],atomix_structure['L3_spectra']['DIMENSIONS']['N_VIB_SENSORS'] ]

atomix_structure['L3_spectra']['VARIABLES']['SH_VIB_SPEC'] = {'required':'optional'}
atomix_structure['L3_spectra']['VARIABLES']['SH_VIB_SPEC']['standard_name'] = 'shear_and_vibration_cross-spectral_matrix'
atomix_structure['L3_spectra']['VARIABLES']['SH_VIB_SPEC']['units'] = ''
atomix_structure['L3_spectra']['VARIABLES']['SH_VIB_SPEC']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L3_spectra']['DIMENSIONS']['N_WAVENUMBER'],atomix_structure['L3_spectra']['DIMENSIONS']['N_SH_ACC_SPEC'] ]

atomix_structure['L3_spectra']['VARIABLES']['SH_ACC_SPEC'] = {'required':'optional'}
atomix_structure['L3_spectra']['VARIABLES']['SH_ACC_SPEC']['standard_name'] = 'shear_and_acceleration_cross-spectral_matrix'
atomix_structure['L3_spectra']['VARIABLES']['SH_ACC_SPEC']['units'] = ''
atomix_structure['L3_spectra']['VARIABLES']['SH_ACC_SPEC']['DIMENSIONS'] = [atomix_structure['L3_spectra']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L3_spectra']['DIMENSIONS']['N_WAVENUMBER'],atomix_structure['L3_spectra']['DIMENSIONS']['N_SH_VIB_SPEC'] ]

# Level4 data
atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'] = {'description':'length of the record of average times of spectral segments. This also equals time of dissipation estimates. '}
atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS'] = {'description':'length of the wavenumber array' }

atomix_structure['L4_dissipation']['VARIABLES']['TIME'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['TIME']['standard_name'] = 'time'
atomix_structure['L4_dissipation']['VARIABLES']['TIME']['units'] = 'CF-Conventions compatible offset'
atomix_structure['L4_dissipation']['VARIABLES']['TIME']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L4_dissipation']['VARIABLES']['SECTION_NUMBER'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['SECTION_NUMBER']['standard_name'] = 'unique_identifier_for_each_section_of_data_from_timeseries'
atomix_structure['L4_dissipation']['VARIABLES']['SECTION_NUMBER']['units'] = ''
atomix_structure['L4_dissipation']['VARIABLES']['SECTION_NUMBER']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L4_dissipation']['VARIABLES']['PSPD_REL'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['PSPD_REL']['standard_name'] = 'platform_speed_wrt_sea_water'
atomix_structure['L4_dissipation']['VARIABLES']['PSPD_REL']['units'] ='m s-1'
atomix_structure['L4_dissipation']['VARIABLES']['PSPD_REL']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L4_dissipation']['VARIABLES']['EPSI'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['EPSI']['standard_name'] = 'specific_turbulent_kinetic_energy_dissipation _in_water'
atomix_structure['L4_dissipation']['VARIABLES']['EPSI']['units'] ='W kg-1'
atomix_structure['L4_dissipation']['VARIABLES']['EPSI']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['EPSI_FINAL'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['EPSI_FINAL']['standard_name'] = 'specific_turbulent_kinetic_energy_dissipation _in_water'
atomix_structure['L4_dissipation']['VARIABLES']['EPSI_FINAL']['units'] ='W kg-1'
atomix_structure['L4_dissipation']['VARIABLES']['EPSI_FINAL']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L4_dissipation']['VARIABLES']['KMIN'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['KMIN']['standard_name'] = 'minimum_wavenumber_used_for_estimating_ turbulent_kinetic_energy_dissipation'
atomix_structure['L4_dissipation']['VARIABLES']['KMIN']['units'] ='cpm'
atomix_structure['L4_dissipation']['VARIABLES']['KMIN']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['KMAX'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['KMAX']['standard_name'] = 'maximum_wavenumber_used_for_estimating_ turbulent_kinetic_energy_dissipation'
atomix_structure['L4_dissipation']['VARIABLES']['KMAX']['units'] ='cpm'
atomix_structure['L4_dissipation']['VARIABLES']['KMAX']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['N_S'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['N_S']['standard_name'] = 'number_of_spectral_points_used_for_estimating_turbulent_kinetic_energy_dissipation'
atomix_structure['L4_dissipation']['VARIABLES']['N_S']['units'] ='-'
atomix_structure['L4_dissipation']['VARIABLES']['N_S']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['EPSI_FLAGS'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['EPSI_FLAGS']['standard_name'] = 'dissipation_qc_flags'
atomix_structure['L4_dissipation']['VARIABLES']['EPSI_FLAGS']['units'] ='-'
atomix_structure['L4_dissipation']['VARIABLES']['EPSI_FLAGS']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['METHOD'] = {'required':'required'}
atomix_structure['L4_dissipation']['VARIABLES']['METHOD']['standard_name'] = 'method_used_for_estimating_ turbulent_kinetic_energy_dissipation'
atomix_structure['L4_dissipation']['VARIABLES']['METHOD']['units'] ='-'
atomix_structure['L4_dissipation']['VARIABLES']['METHOD']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['PRES'] = {'required':'optional'}
atomix_structure['L4_dissipation']['VARIABLES']['PRES']['standard_name'] = 'water_pressure'
atomix_structure['L4_dissipation']['VARIABLES']['PRES']['units'] = 'dbar'
atomix_structure['L4_dissipation']['VARIABLES']['PRES']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L4_dissipation']['VARIABLES']['KVISC'] = {'required':'optional'}
atomix_structure['L4_dissipation']['VARIABLES']['KVISC']['standard_name'] = 'kinematic_viscosity_of_water'
atomix_structure['L4_dissipation']['VARIABLES']['KVISC']['units'] = 'm2 s-1'
atomix_structure['L4_dissipation']['VARIABLES']['KVISC']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA']]

atomix_structure['L4_dissipation']['VARIABLES']['FOM'] = {'required':'optional'}
atomix_structure['L4_dissipation']['VARIABLES']['FOM']['standard_name'] = 'figure_of_merit'
atomix_structure['L4_dissipation']['VARIABLES']['FOM']['units'] = '-'
atomix_structure['L4_dissipation']['VARIABLES']['FOM']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['MAD'] = {'required':'optional'}
atomix_structure['L4_dissipation']['VARIABLES']['MAD']['standard_name'] = 'mean_absolute_deviation'
atomix_structure['L4_dissipation']['VARIABLES']['MAD']['units'] = '-'
atomix_structure['L4_dissipation']['VARIABLES']['MAD']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['VAR_RESOLVED'] = {'required':'optional'}
atomix_structure['L4_dissipation']['VARIABLES']['VAR_RESOLVED']['standard_name'] = 'variance_resolved'
atomix_structure['L4_dissipation']['VARIABLES']['VAR_RESOLVED']['units'] = '-'
atomix_structure['L4_dissipation']['VARIABLES']['VAR_RESOLVED']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['EPSI_STD'] = {'required':'optional'}
atomix_structure['L4_dissipation']['VARIABLES']['EPSI_STD']['standard_name'] = 'expected_standard_deviation_of_the_ logarithm_of_the_dissipation_estimate'
atomix_structure['L4_dissipation']['VARIABLES']['EPSI_STD']['units'] = '-'
atomix_structure['L4_dissipation']['VARIABLES']['EPSI_STD']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['DESPIKE_FRACTION_SH'] = {'required':'optional'}
atomix_structure['L4_dissipation']['VARIABLES']['DESPIKE_FRACTION_SH']['standard_name'] = 'fraction_of_shear_data_modified_by_despiking_algorithm'
atomix_structure['L4_dissipation']['VARIABLES']['DESPIKE_FRACTION_SH']['units'] = '-'
atomix_structure['L4_dissipation']['VARIABLES']['DESPIKE_FRACTION_SH']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

atomix_structure['L4_dissipation']['VARIABLES']['DESPIKE_PASS_COUNT_SH'] = {'required':'optional'}
atomix_structure['L4_dissipation']['VARIABLES']['DESPIKE_PASS_COUNT_SH']['standard_name'] = 'number_of_despike_passes_for_shear_probes'
atomix_structure['L4_dissipation']['VARIABLES']['DESPIKE_PASS_COUNT_SH']['units'] = '-'
atomix_structure['L4_dissipation']['VARIABLES']['DESPIKE_PASS_COUNT_SH']['DIMENSIONS'] = [atomix_structure['L4_dissipation']['DIMENSIONS']['TIME_SPECTRA'],atomix_structure['L4_dissipation']['DIMENSIONS']['N_SHEAR_SENSORS']]

class atomix():
    """An object to read atomix netCDF files

    """ 
    def __init__(self,filename,verbosity = logging.DEBUG):
        logger.setLevel(verbosity)
        logger.info(' Opening file: {:s}'.format(filename))
        self.filename = filename
        basename = os.path.basename(filename)
        nc = netCDF4.Dataset(filename)
        self.nc = nc


    def check_requirements(self,req_levels=['required','highly-recommended','optional']):
        """
        Checks if the netCDF file fulfills the requirements of an Atomix netCDF file.

        Returns:

        """
        #funcname = self.__class__.__name__ + '.check_requirements():'
        funcname = 'check:'
        groups = list(atomix_structure.keys())
        groups_found = []
        for g in self.nc.groups.keys():
            flag_atomix_standard = False
            if g in groups:
                flag_atomix_standard = True
                groups_found.append(g)

            logger.debug(funcname + ' Group {:s} which is ATOMIX standard: {:s}'.format(g,str(flag_atomix_standard)))
            
        logger.debug(funcname + ' Found {:d} ATOMIX level of data {:s}'.format(len(groups_found),str(groups_found)))
        for req_level in req_levels:
           logger.info(funcname + 'Checking variables for requirement: {:s}'.format(req_level))
           # Fill a dictionary with the requirements of variables
           variable_required = {}
           for g in groups:
               variable_required[g] = []
               for var in atomix_structure[g]['VARIABLES'].keys():
                   req = atomix_structure[g]['VARIABLES'][var]['required']
                   if req_level in req:
                       variable_required[g].append(var)

           # Check requirements
           variables_req_found = {}
           variables_req_notfound = {}
           for g in groups_found:
               if len(variable_required[g]) == 0:
                  logger.info(funcname + ' Group: {:s}: No {:s} variables specified by ATOMIX'.format(g, req_level))
               else:
                  variables_req_found[g] = []
                  variables_req_notfound[g] = []
                  for var in self.nc.groups[g].variables.keys():
                      flag_atomix_standard = False
                      if var in variable_required[g]:
                          flag_atomix_standard = True
                          variables_req_found[g].append(var)
                      logger.debug(
                          funcname + ' Group {:s}: Var {:s} which is ATOMIX standard: {:s}'.format(g, var, str(flag_atomix_standard)))

                  # Check the variable that have not been found
                  for var in variable_required[g]:
                      if var not in variables_req_found[g]:
                          variables_req_notfound[g].append(var)


                  #print('req',variables_req_found[g],variable_required[g])
                  #print('not found', variables_req_notfound[g])
                  if len(variables_req_notfound[g]) == 0:
                      logger.info(funcname + ' Group: {:s}: ALL {:s} variables found: {:s}'.format(g, req_level, str(variables_req_found[g])))
                  else:
                      logger.info(funcname + ' Group: {:s}: MISSING {:s} variables: {:s} '.format(g, req_level, str(variables_req_notfound[g])))                     
                      logger.info(funcname + ' Group: {:s}: {:s} variables found: {:s}'.format(g, req_level, str(variables_req_found[g])))               







