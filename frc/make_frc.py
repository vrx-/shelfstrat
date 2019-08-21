
import numpy as np
import xarray as xr
from datetime import datetime


def make_frc(output='../tests/shelfstrat_frc.nc',
             u=0.0, v=5.0, Tramp=1.0, Cd=1.5e-3, ndays=365.):

    sustr0 = Cd * np.sqrt(u**2 + v**2) * u
    svstr0 = Cd * np.sqrt(u**2 + v**2) * v
    t = xr.DataArray(np.linspace(0, ndays, ndays * 24 * 4), dims=['sms_time'])
    ramp = (1.0 - np.exp(-t / Tramp)) * np.sin(t * 2.0 * np.pi)
    sustr = xr.DataArray(sustr0 * ramp, dims=['sms_time'])
    svstr = xr.DataArray(svstr0 * ramp, dims=['sms_time'])

    # Create dataset

    ds = xr.Dataset({'sms_time': t, 'sustr': sustr, 'svstr': svstr})

    ds.attrs['Description'] = 'Initial conditions for ideal shelf'
    ds.attrs['Author'] = 'Vero and Lixin'
    ds.attrs['Created'] = datetime.now().isoformat()
    ds.attrs['type'] = 'ROMS FRC file'
    ds['sms_time'].attrs['units'] = 'days'
    ds['sustr'].attrs['units'] = 'Newton meter-2'
    ds['svstr'].attrs['units'] = 'Newton meter-2'

    print('Writing netcdf FRC file..')
    ds.to_netcdf(output)


if __name__ == '__main__':
    make_frc()
