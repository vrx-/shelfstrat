#!/bin/env python3.7
# encoding: utf-8
"""
run_case.py

Created by Rob Hetland on 2007-10-15.
Modified by Veronica Ruiz Xomchuk on 2019-08
Copyright (c) 2019 Texas A&M Univsersity. All rights reserved.
Release under MIT license.
"""
from grd import make_grd
from frc import make_frc
from ini import make_ini

import os
import warnings


class ROMS_in(object):
    """docstring for ROMS_IN"""

    def __init__(self, infile):
        self.infile = infile
        f = open(self.infile)

        self.variables = {}
        self._varlist = []       # to keep the order of variables

        for line in f.readlines():
            # strip comments
            s = line.split('!')[0].strip()
            if len(s) == 0:
                continue     # comment only
            vals = s.split('=')
            self.variables[vals[0].strip()] = vals[-1].strip()
            self._varlist.append(vals[0].strip())

    def write(self, filename):
        """docstring for write"""
        f = open(filename, 'w')
        for key in self._varlist:
            f.write('%s == %s\n' % (key, str(self.variables[key])))
        f.close()

    def __setitem__(self, key, val):
        """docstring for __setitem__"""
        if key not in self._varlist:
            self._varlist.append(key)
            warnings.warn('%s not previously in variable list.' % key)
        self.variables[key] = str(val)


def run_case(case, z0=0.003, dt=30.0, exec=False, rootdir='runs/'):

    print('BUILD case ID %s' % case['ID'])
    if not os.path.exists(rootdir):
        os.makedirs(rootdir)
    grd_name = rootdir + 'shelf_' + case['ID'] + '_grd.nc'
    frc_name = rootdir + 'shelf_' + case['ID'] + '_frc.nc'
    ini_name = rootdir + 'shelf_' + case['ID'] + '_ini.nc'

    make_grd(grd_name,
             Hmin=case['grd']['Hmin'],
             alpha=case['grd']['alpha'],
             ho=case['grd']['ho'],
             dh=case['grd']['dh'],
             wdh=case['grd']['wdh'],
             f=case['grd']['f'],
             dx=case['grd']['dx'],
             dy=case['grd']['dy'],
             shp=case['grd']['shp'])
    make_frc(frc_name,
             u=case['frc']['u'],
             v=case['frc']['v'],
             Tramp=case['frc']['Tramp'],
             Cd=case['frc']['Cd'])
    make_ini(ini_name, grd_name,
             case['ini']['zlevs'],
             case['ini']['theta_s'],
             case['ini']['theta_b'],
             case['ini']['hc'],
             case['ini']['R0'],
             case['ini']['T0'],
             case['ini']['S0'],
             case['ini']['TCOEF'],
             case['ini']['SCOEF'],
             case['ini']['M20'],
             case['ini']['M2_yo'],
             case['ini']['M2_r'],
             case['ini']['N20'],
             case['ini']['N2_zo'],
             case['ini']['N2_r'])

    infile = os.path.join(rootdir, 'ocean_shelf_' + case['ID'] + '.in')
    outfile = os.path.join(rootdir, 'ocean_shelf_' + case['ID'] + '.out')
    # run 3D case
    rin_3d = ROMS_in('./project/ocean_shelfstrait.in')
    rin_3d['GRDNAME'] = os.path.join(rootdir, grd_name)
    rin_3d['FRCNAME'] = os.path.join(rootdir, frc_name)
    rin_3d['HISNAME'] = os.path.join(rootdir, 'shelf_' + case['ID'] + '_his.nc')
    rin_3d['AVGNAME'] = os.path.join(rootdir, 'shelf_' + case['ID'] + '_avg.nc')
    rin_3d['DIANAME'] = os.path.join(rootdir, 'shelf_' + case['ID'] + '_dia.nc')
    rin_3d['ININAME'] = os.path.join(rootdir, ini_name)
    rin_3d['RSTNAME'] = os.path.join(rootdir, 'shelf_' + case['ID'] + '_rst.nc')
    rin_3d['VARNAME'] = './project/varinfo.dat'

    rin_3d['Lm'] = case['grd']['shp'][1]
    rin_3d['Mm'] = case['grd']['shp'][0]
    rin_3d['N'] = case['ini']['zlevs']

    rin_3d['NTIMES'] = int(86400 * case['frc']['ndays'] / dt)
    rin_3d['DT'] = dt
    rin_3d['NDTFAST'] = int(30)
    rin_3d['NHIS'] = int((86400.0 / dt) / 8.0)
    rin_3d['NAVG'] = int((86400.0 / dt) * 3.0)
    rin_3d['NDIA'] = int((86400.0 / dt) * 1.0)

    rin_3d['Zob'] = z0
    rin_3d.write(infile)
    print('RUN case ID %s' % case['ID'])
    print(infile)

    if exec:
        print(' ### Running 3D ROMS...')
        os.system('/usr/mpi/gcc/openmpi-1.4.3/bin/mpiexec -np 8 ./project/coawstM %s > %s &' % (infile, outfile))

    return 0


if __name__ == '__main__':

    case = {'ID': 'ho_5_wdh_1e4',
            'grd': {'Hmin': 5.0,
                    'alpha': 0.001,
                    'ho': 5.,
                    'dh': 0.,
                    'wdh': 1e4,
                    'f': 1e-4,
                    'dx': 1e3,
                    'dy': 1e3,
                    'shp': (131, 259),
                    },
            'frc': {'u': 0.,
                    'v': 5.,
                    'Tramp': 1.,
                    'Cd': 1.5e-3,
                    'ndays': 365,
                    },
            'ini': {'zlevs': 30,
                    'theta_s': 3.,
                    'theta_b': .4,
                    'hc': 5.,
                    'R0': 1027.,
                    'T0': 25.,
                    'S0': 35.,
                    'TCOEF': 1.7e-4,
                    'SCOEF': 7.6e-4,
                    'M20': 1e-7,
                    'M2_yo': 50e3,
                    'M2_r': 5e3,
                    'N20': 1e-4,
                    'N2_zo': 50.,
                    'N2_r': 50.,
                    },
            }

    run_case(case)
