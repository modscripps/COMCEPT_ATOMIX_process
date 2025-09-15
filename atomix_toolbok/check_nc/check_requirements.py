import atomix_check
import logging
fnames = []
#fnames.append('benchmark_download_2023.10.19/Data/MSS_BalticSea/MSS_Baltic_0330.nc')
#fnames.append('./benchmark_download_2023.10.19/Data/Nemo_MR1000_Minas_Passage/Nemo_MR1000_Minas_Passage_InStream.nc')
#fnames.append('benchmark_download_2023.10.19/Data/VMP2000_FaroeBankChannel/VMP2000_FaroeBankChannel.nc')
#fnames.append('benchmark_download_2023.10.19/Data/VMP250_TidalChannel/VMP250_TidalChannel_024.nc')
fnames.append('../../atomix_create_compare_netcdf_ALB/TestData/epsifish_epsilometer_northatlantic_with_ns.nc')
#fnames.append('dengler_download_2023.10.25/Glider_MicroRider_Sal_eddy.nc')

for fname in fnames:
    print('---')
    print('Checking variables of {:s}'.format(fname.split('/')[-1]))
    ato_nc = atomix_check.atomix(fname,verbosity=logging.INFO)
    ato_nc.check_requirements(req_levels=['required'])
    print('---')


print('-------------------')
print('-------------------')
print('-------------------')
print('Detailed check of optional and highly recommended variables')
for fname in fnames:
    print('---')
    print('Checking variables of {:s}'.format(fname.split('/')[-1]))
    ato_nc = atomix_check.atomix(fname,verbosity=logging.INFO)
    ato_nc.check_requirements(req_levels=['highly-recommended','optional'])
    print('---')

