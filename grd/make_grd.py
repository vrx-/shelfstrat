""" This script creates a grid file for an idealized near-field plume case to evaluate bottom traping by a bump
"""
import numpy as np
import xarray as xr


def make_grd(output='../tests/shelfstrat_grd.nc',
             Hmin=5.0, alpha=0.001,
             ho=20., dh=0., wdh=1e4,
             f=1e-4,
             dx=1e3, dy=1e3,
             shp=(131, 259),
             spherical=True,
             angle=0.):

    x = np.arange(shp[1], dtype='d') * dx
    y = np.arange(shp[0], dtype='d') * dy
    xvert, yvert = np.meshgrid(x, y)
    grd = make_CGrid(xvert, yvert)

    grd['xl'] = 1
    grd.xl.attrs['long_name'] = 'domain length in the XI-direction'
    grd.xl.attrs['units'] = 'meter'
    grd['el'] = 1
    grd.el.attrs['long_name'] = 'domain length in the ETA-direction'
    grd.el.attrs['units'] = 'meter'
    grd['spherical'] = spherical
    grd.spherical.attrs['long_name'] = 'Grid type logical switch'
    grd['f'] = f
    grd.f.attrs['long_name'] = 'Coriolis parameter at RHO-points'
    grd.f.attrs['units'] = 'second-1'
    grd.f.attrs['field'] = 'Coriolis, scalar'
    grd['angle'] = angle
    grd.f.attrs['long_name'] = 'angle between xi axis and east'
    grd.f.attrs['units'] = 'degree'

    # create a depth profile with slope alpha, a value of Hmin
    # and a bump localized approximately at a depth ho, as a function of step height (dh) and width (wdh).
    cff1 = alpha * grd.y_rho + dh * np.tanh((grd.y_rho - (ho - Hmin) / alpha) / wdh * 4) + 2 * Hmin
    cff1 += 0.01 * np.random.randn(*grd.y_rho.shape) * cff1
    grd['h'] = np.maximum(cff1, Hmin)
    grd.h.attrs['long_name'] = 'Final bathymetry at RHO-points'
    grd.h.attrs['units'] = 'meter'
    grd.h.attrs['field'] = 'bath, scalar'
    grd['mask_rho'] = grd.mask_rho
    grd.mask_rho.attrs['long_name'] = 'Final bathymetry at RHO-points'
    grd.mask_rho.attrs['units'] = 'meter'
    grd.mask_rho.attrs['field'] = 'bath, scalar'
    print('Writing netcdf GRD file..')
    grd.to_netcdf(output)


def make_CGrid(x, y):
    if np.any(np.isnan(x)) or np.any(np.isnan(y)):
        x = np.ma.masked_where((np.isnan(x)) | (np.isnan(y)), x)
        y = np.ma.masked_where((np.isnan(x)) | (np.isnan(y)), y)

    ds = xr.Dataset({'x_vert': (['eta_vert', 'xi_vert'], x),
                     'y_vert': (['eta_vert', 'xi_vert'], y)})
    ds['x_rho'] = (['eta_rho', 'xi_rho'], 0.25 * (x[1:, 1:] + x[1:, :-1] + x[:-1, 1:] + x[:-1, :-1]))
    ds['y_rho'] = (['eta_rho', 'xi_rho'], 0.25 * (y[1:, 1:] + y[1:, :-1] + y[:-1, 1:] + y[:-1, :-1]))
    ds['x_psi'] = (['eta_psi', 'xi_psi'], x[1:-1, 1:-1])
    ds['y_psi'] = (['eta_psi', 'xi_psi'], y[1:-1, 1:-1])
    ds['x_u'] = (['eta_u', 'xi_u'], 0.5 * (x[:-1, 1:-1] + x[1:, 1:-1]))
    ds['y_u'] = (['eta_u', 'xi_u'], 0.5 * (y[:-1, 1:-1] + y[1:, 1:-1]))
    ds['x_v'] = (['eta_v', 'xi_v'], 0.5 * (x[1:-1, :-1] + x[1:-1, 1:]))
    ds['y_v'] = (['eta_v', 'xi_v'], 0.5 * (y[1:-1, :-1] + y[1:-1, 1:]))
    ds['x_psi'] = (['eta_psi', 'xi_psi'], x[1:-1, 1:-1])
    ds['y_psi'] = (['eta_psi', 'xi_psi'], y[1:-1, 1:-1])
    x_temp = 0.5 * (ds.x_vert[1:, :] + ds.x_vert[:-1, :])
    y_temp = 0.5 * (ds.y_vert[1:, :] + ds.y_vert[:-1, :])
    dx = np.sqrt(np.diff(x_temp, axis=1)**2 + np.diff(y_temp, axis=1)**2)
    x_temp = 0.5 * (ds.x_vert[:, 1:] + ds.x_vert[:, :-1])
    y_temp = 0.5 * (ds.y_vert[:, 1:] + ds.y_vert[:, :-1])
    dy = np.sqrt(np.diff(x_temp, axis=0)**2 + np.diff(y_temp, axis=0)**2)

    ds['pm'] = (['eta_rho', 'xi_rho'], 1. / dx)
    ds['pn'] = (['eta_rho', 'xi_rho'], 1. / dy)

    try:
        ds['mask_rho'] = ds.pm.mask
    except:
        print('no mask')
    return ds


if __name__ == '__main__':
    make_grd(output='../tests/shelfstrat_grd.nc')
