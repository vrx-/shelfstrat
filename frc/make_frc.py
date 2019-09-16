
import numpy as np
import xarray as xr
from datetime import datetime


def make_frc(output='../tests/shelfstrat_frc.nc',
             uwind=0.0, vwind=5.0, Cd=1.5e-3, Rho0=1027.,
             ndays=30, dtw=1 / 24, Tramp=3.0, Tflat=3.):

    sustr0 = Cd * np.sqrt(uwind**2 + vwind**2) * uwind
    svstr0 = Cd * np.sqrt(uwind**2 + vwind**2) * vwind
    t = xr.DataArray(np.arange(0, ndays+dtw, dtw), dims=['sms_time'])
    ramp = (1.0 - np.exp(-t / Tramp)) * np.sin(t * 2.0 * np.pi)
    nt = int(Tflat / dtw)
    sustr = xr.DataArray(sustr0 * np.insert(ramp, 0, np.zeros(nt))[:-nt], dims=['sms_time'])
    svstr = xr.DataArray(svstr0 * np.insert(ramp, 0, np.zeros(nt))[:-nt], dims=['sms_time'])

    # Create dataset

    ds = xr.Dataset({'sms_time': t, 'sustr': sustr, 'svstr': svstr})

    ds.attrs['Description'] = 'Initial conditions for ideal shelf'
    ds.attrs['Author'] = 'Vero and Lixin'
    ds.attrs['Created'] = datetime.now().isoformat()
    ds.attrs['type'] = 'ROMS FRC file'
    ds['sms_time'].attrs['units'] = 'days'
    ds['sustr'].attrs['units'] = 'Newton meter-2'
    ds['svstr'].attrs['units'] = 'Newton meter-2'

    print('Writing netcdf FRC file: '+output)
    ds.to_netcdf(output)


if __name__ == '__main__':
    make_frc()
