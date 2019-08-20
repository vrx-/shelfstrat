import numpy as np
import xarray as xr
from datetime import datetime


def get_depths(h, hc, s, Cs, dim):
    """SEE Eq. (2) or (3) on https://www.myroms.org/wiki/Vertical_S-coordinate"""
    return (hc * s + h * Cs) / (hc + h) * h


def C(theta_s, theta_b, s):
    C = (1.0 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1.0)
    if theta_b > 0.0:
        return (np.exp(theta_b * C) - 1.0) / (1.0 - np.exp(-theta_b))
    else:
        return -C


def make_ini(grd_path='../tests/shelfstrat_grd.nc', output='../tests/shelfstrat_ini.nc',
             zlevs=30, theta_s=3.0, theta_b=0.4, hc=5.0,
             R0=1027.0, T0=25.0, S0=35.0, TCOEF=1.7e-4, SCOEF=7.6e-4,
             M20=1e-7, M2_yo=50e3, M2_r=5e3,
             N20=1e-4, N2_zo=50.0, N2_r=50.0,
             balanced_run=True):
    '''
    Create an initialization file.

    Horizontal stratification is controlled by salinity (only)
    Vertical stratification is controlled by temperature (only)

    Stratification properties are conservative through
    a linear equation of state:

       R0 == 1027.0d0                   ! kg/m3
       T0 == 25.0d0                      ! Celsius
       S0 == 35.0d0                     ! PSU
    TCOEF == 1.7d-4                     ! 1/Celsius
    SCOEF == 7.6d-4                     ! 1/PSU
    '''
    grd = xr.open_dataset(grd_path)
    g = 9.81
    dy = 1 / grd.pn
    s_w = xr.DataArray(np.linspace(-1., 0., zlevs + 1), dims=['s_w'])
    s_rho = np.linspace(-1., 0., zlevs + 1)
    s_rho = s_rho[:-1] + np.diff(s_rho) / 2
    s_rho = xr.DataArray(s_rho, dims=['s_rho'])
    Cs_r = C(theta_s, theta_b, s_rho)
    Cs_w = C(theta_s, theta_b, s_w)

    M2 = M20 * np.exp((M2_yo - grd.y_rho) / M2_r)
    M2 = M2.where(grd.y_rho > M2_yo, M20)
    salt = (M2 * dy / g / SCOEF).cumsum(axis=0)
    salt -= salt[-1] - S0
    salt = salt.expand_dims('s_rho') * np.ones((zlevs, 1, 1), 'd')
    salt.coords['s_rho'] = s_rho
    # (h, hc, s, Cs)
    z = get_depths(grd.h, hc, s_rho, Cs_r, 's_rho')
    Hz = get_depths(grd.h, hc, s_w, Cs_w, 's_w').diff('s_w').rename({'s_w': 's_rho'})
    Hz.coords['s_rho'] = s_rho
    N2 = N20 * np.exp(-(N2_zo - z) / N2_r)
    N2 = N2.where(z > N2_zo, N20)

    temp = xr.zeros_like(salt)
    for n in range(zlevs):
        temp[n] = T0 - np.trapz(N2[n:] / (g * TCOEF), x=z[n:], axis=0)

    #########################################
    # Create dataset

    ds = xr.Dataset({'temp': temp, 'salt': salt,
                     's_rho': s_rho, 'xi_rho': grd.xi_rho, 'eta_rho': grd.eta_rho})

    ds.attrs['Description'] = 'Initial conditions for ideal shelf'
    ds.attrs['Author'] = 'Vero and Lixin'
    ds.attrs['Created'] = datetime.now().isoformat()
    ds.attrs['type'] = 'ROMS FRC file'
    ds['ocean_time'] = xr.DataArray([0.0], dims=['ocean_time'])
    ds['ocean_time'].attrs['units'] = 'days'
    print('Writing netcdf INI file..')
    ds.to_netcdf(output)


if __name__ == '__main__':
    make_ini()
