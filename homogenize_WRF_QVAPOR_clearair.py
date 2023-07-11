# The name homogenize_WRF_QVAPOR_clearair says it all I hope
# clearair is defined by ds.QCLOUD[0,ilev,:,:] <1e-7
# WRF seems to have actual QCLOUD zeros, but 1e-7 kg/kg is safe 
# Usage: python homogenize_WRF_QVAPOR_clearair.py filename

import sys
import xarray as xr

# Accessing the script's filename
script_name = sys.argv[0]

# Accessing the passed arguments
filename = sys.argv[1]

# Open the file
try: 
    ds = xr.open_dataset(filename)
except: 
    print("First argument should be filename of a WRF restart netCDF file: like this")
    print("python homogenize_WRF_QVAPOR_clearair.py filename")
    
    
# Loop over levels, compute the eddy component (deviation from spatial mean)
# for clear air only: where QCLOUD < 1e-7, in kg/kg units I think)

for ilev in range(ds.bottom_top.size):
    meanQV_clr = ds.QVAPOR[0,ilev,:,:].where(ds.QCLOUD[0,ilev,:,:] <1e-7).mean()
    eddy_clr = ds.QVAPOR[0,ilev,:,:] - meanQV_clr

# subtract the eddy component, again only where masked as clear   
    ds.QVAPOR[0,ilev,:,:] -= eddy_clr *(ds.QCLOUD[0,ilev,:,:]<1e-7)

# Write out the result
ds.to_netcdf(filename+'_homogQVAPOR.nc')