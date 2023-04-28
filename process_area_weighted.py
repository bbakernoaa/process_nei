#!/usr/bin/env python
import xarray as xr
import xesmf as xe
import monetio as mio
import pyproj
from pyproj import Geod
from numpy import zeros, arange, ones
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
import pandas as pd
import numpy as np
import numba
import PseudoNetCDF as pnc
import warnings
warnings.filterwarnings('ignore')

def open_cmaq_file(fname):
    print('open with pnc')
    cmaqf = xr.open_dataset(fname,engine='pseudonetcdf',backend_kwargs=dict(format='ioapi', addcf=True),decode_cf=False,decode_times=False).set_coords(["latitude", "longitude", "time"]).rename({"COL": "x", "ROW": "y"}).squeeze().drop(["TFLAG",'layer','level','time_bounds'])
    cmaqf.coords["TSTEP"] = cmaqf.time
    cmaqf = cmaqf.drop(["time"]).rename({"TSTEP": "time"})
    print('getting corners')
    lon_b,lat_b = get_corners(cmaqf)
    cmaqf =  condition_cmaq_dset(cmaqf,lat_b,lon_b)
    print('returning dataset')
    return cmaqf.drop(['lambert_conformal_conic','latitude_bounds','longitude_bounds'])

def get_corners(dset):
    lon_b = zeros((dset.dims['y']+1,dset.dims['x']+1))
    lat_b = zeros((dset.dims['y']+1,dset.dims['x']+1))
    proj4_srs = "+proj=lcc +lat_1={} +lat_2={} +lat_0={} +lon_0{} +x_0={} +y_0={} +datum=WGS84 +units=m +a={} +b={}".format(
        dset.lambert_conformal_conic.attrs["latitude_of_projection_origin"],
        dset.lambert_conformal_conic.attrs["latitude_of_projection_origin"],
        dset.lambert_conformal_conic.attrs["latitude_of_projection_origin"],
        dset.lambert_conformal_conic.attrs["longitude_of_central_meridian"],
        dset.lambert_conformal_conic.attrs["false_easting"],
        dset.lambert_conformal_conic.attrs["false_northing"],
        dset.lambert_conformal_conic.attrs["semi_major_axis"],
        dset.lambert_conformal_conic.attrs["semi_minor_axis"]
        )
    p = pyproj.Proj(proj4_srs)
    lon_b,lat_b = zeros((dset.dims['y']+1,dset.dims['x']+1)), zeros((dset.dims['y']+1,dset.dims['x']+1))
    xx,yy = p(dset.longitude.values,dset.latitude.values)
    lon_b_upper,lat_b_upper = p(xx+dset.XCELL / 2., yy+dset.YCELL /2.,inverse=True)
    lon_b_lower,lat_b_lower = p(xx - dset.XCELL/2, yy-dset.YCELL /2.,inverse=True)
    lon_b[1:,1:], lat_b[1:,1:] = lon_b_upper,lat_b_upper
    lon_b[:-1,:-1],lat_b[:-1,:-1] = lon_b_lower,lat_b_lower
    # the upper left corner and lower right corner didn't get filled here
    # so fill them manually
    lon_b[0,-1],lat_b[0,-1] = p(xx[0,-1]+dset.XCELL/2.,yy[0,-1] - dset.YCELL/2.,inverse=True)
    lon_b[-1,0],lat_b[-1,0] = p(xx[-1,0]-dset.XCELL/2.,yy[-1,0] + dset.YCELL/2.,inverse=True)
    return lon_b,lat_b

def get_proj4_corners(dset, projstr):
     p = pyproj.Proj(dset.proj4_srs)
     lon_b,lat_b = zeros((dset.dims['y']+1,dset.dims['x']+1)), zeros((dset.dims['y']+1,dset.dims['x']+1))
     xx,yy = p(dset.longitude.values,dset.latitude.values)
     lon_b_upper,lat_b_upper = p(xx+dset.XCELL / 2., yy+dset.YCELL /2.,inverse=True)
     lon_b_lower,lat_b_lower = p(xx - dset.XCELL/2, yy-dset.YCELL /2.,inverse=True)
     lon_b[1:,1:], lat_b[1:,1:] = lon_b_upper,lat_b_upper
     lon_b[:-1,:-1],lat_b[:-1,:-1] = lon_b_lower,lat_b_lower
     # the upper left corner and lower right corner didn't get filled here
     # so fill them manually
     lon_b[0,-1],lat_b[0,-1] = p(xx[0,-1]+dset.XCELL/2.,yy[0,-1] - dset.YCELL/2.,inverse=True)
     lon_b[-1,0],lat_b[-1,0] = p(xx[-1,0]-dset.XCELL/2.,yy[-1,0] + dset.YCELL/2.,inverse=True)
     return lon_b,lat_b

def condition_cmaq_dset(dset,lat_b,lon_b):
    dset['lat_b'] = (('y_b','x_b'),lat_b)
    dset['lon_b'] = (('y_b','x_b'),lon_b)
    return dset.set_coords(['lat_b','lon_b']).rename({'latitude':'lat','longitude':'lon'})

def create_target_latlon_grid(dset,res):
    target = xe.util.grid_2d(dset.lon.min()-.05,dset.lon.max()+.05,res,dset.lat.min()-0.05,dset.lat.max()+0.05,res)
    return target

def get_target_area(dset):
    # Please see https://pyproj4.github.io/pyproj/stable/examples.html#geodesic-area
    area = zeros((dset.dims['y'],dset.dims['x']))
    for j in range(dset.dims['y']):
        for i in range(dset.dims['x']):
            geod = Geod(ellps='WGS84')
            lons = [dset.lon_b.values[j,i],dset.lon_b.values[j,i+1],dset.lon_b.values[j+1,i+1],dset.lon_b.values[j+1,i]]
            lats = [dset.lat_b.values[j,i],dset.lat_b.values[j,i+1],dset.lat_b.values[j+1,i+1],dset.lat_b.values[j+1,i]]
            area[j,i], _ = geod.polygon_area_perimeter(lons, lats)
    dset['area'] = (('y','x'),area)
    return dset

def create_regridder(source,target,method='conservative_normed',weights=None):
    if weights is not None:
        r = xe.Regridder(source,target,method=method,weights=weights)
    else:
        r = xe.Regridder(source,target,method=method)
    return r

def regrid_dataset(regrid_obj,source):
    return regrid_obj(source)

def convert_units(da,scale_factor=None):
    if 'g/s' in da.attrs['units']:
        return da.attrs['units'].replace('g/s','kg m-2 s-1')
    if 'moles/s' in da.attrs['units']:
        return da.attrs['units'].replace('moles/s','kg m-2 s-1')

def get_species_mw(filename='species_weights.txt'):
    import pandas as pd
    df = pd.read_csv(filename,names=['','species','long_name','molecular_weight','a','b'],sep='|')[['species','molecular_weight']]
    df['species'] = df.species.str.strip()
    d = df.set_index('species').to_dict()['molecular_weight']
    return d

def convert_to_COARDS_time(dset):
    import pandas as pd
    from numpy import timedelta64
    index_in_hours = (dset.time - dset.time[0]) / timedelta64(1,'h')
    hours_since = pd.Timestamp(dset.time.values[0]).strftime('hours since %Y-%m-%d 00:00:00')
    dset['time'] = index_in_hours
    dset.time.attrs['units'] = hours_since
    dset.time.attrs['long_name'] = 'time'
    dset.time.attrs['calendar'] = 'standard'
    dset.time.attrs['axis'] = 'T'
    return dset

def convert_to_COARDS_latlon(dset):
    dset.lat.attrs['long_name'] = 'latitude'
    dset.lat.attrs['units'] = 'degrees_north'
    dset.lon.attrs['long_name'] = 'longitude'
    dset.lon.attrs['units'] = 'degrees_east'
    lat = dset.lat[:,0]
    lon = dset.lon[0,:]
    dset['x'] = lon
    dset['y'] = lat
    return dset.drop(['lat','lon']).rename({'y':'lat','x':'lon'})

def convert_to_COARDS_lev(dset):
    import xarray as xr
    import numpy as np
    lev=xr.DataArray(np.arange(dset.z.ndim),dims="z")
    dset['lev'] = lev
    dset.lev.attrs['units'] = 'levels'
    dset=dset.rename(dict(z='lev'))
    return dset

def convert_to_COARDS_area(dset):
    dset.area.attrs['units'] = 'm**2'
    return dset

def write_ncf(dset,outfile):
    print('Output File:', outfile)
    encoding = {}
    for v in dset.data_vars:
        encoding[v] = dict(zlib=True, complevel=4)
    encoding["lat"] = dict(zlib=True, complevel=4)
    encoding["lon"] = dict(zlib=True, complevel=4)
    if 'lat_b' in dset:
        encoding["lat_b"] = dict(zlib=True, complevel=4)
        encoding["lon_b"] = dict(zlib=True, complevel=4)
    dset.load().to_netcdf(outfile, encoding=encoding)

def process(infile,outfile,target_file,weight_file,convert=False, area=None, verify=True, target_res=0.125):
    # first check if something is wrong with the time slice in the output file if exists
    if os.path.exists(outfile.replace('nc','nc4')):
        print('Intermediate file exists..... opening')
        out = xr.open_dataset(outfile.replace('nc','nc4'),decode_times=False,decode_cf=False)
        cmaq = mio.cmaq.open_dataset(infile)
    else:
        # open cmaq file
        print('Opening CMAQ file:', infile)
        c = open_cmaq_file(infile)
        orig = c.copy()
        # create area weighted cmaq file (units/m2/second)
        print('Getting Area')
        if area is None:
            print('Area not found.... generating')
            c = get_target_area(c)
        else:
            c['area'] = area['area']
        #print('creating area weighted variables')
        for v in c.data_vars:
           if v != 'area':
               print('          ',v)
               attrs = c[v].attrs
               c[v] = (c[v] / c.area.data).astype('float32')
               c[v].attrs = attrs
        # check if target file exists
        print('Getting target')
        if os.path.isfile(target_file):
            t = xr.open_dataset(target_file)

        else:
            # create boundary from minimum and maximum latitudes of the cmaq file
            t = create_target_latlon_grid(c, target_res)

            # calculate the area using the boundaries and pyproj (WSG84)
            print('Getting target area')
            t = get_target_area(t)

            #save target file for later use
            write_ncf(target_file)

        # regrid object
        print('Creating Regridder Object')
        if weight_file is not None and os.path.isfile(weight_file):
            print('.....creating using weight file', weight_file)
            r = create_regridder(c,t,method='conservative_normed',weights=weight_file)
        else:
            r = create_regridder(c,t,method='conservative_normed')
            r.to_netcdf()

        print('Regrid CMAQ file to target')
        out = regrid_dataset(r, c)
        out['area'] = t.area

        # add netcdf attributes back in
        for v in c.data_vars:
            out[v].attrs = c[v].attrs
            # convert to area weighted
#            if v != 'area':
#                out[v].data = (out[v].data/out['area'].data).astype('float32')

#        out.attrs['Convention'] = 'COARDS'
#        out.attrs['Format'] = 'NetCDF-4'
#        out = convert_to_COARDS_time(out)
#        out = convert_to_COARDS_latlon(out)
#        out = convert_to_COARDS_lev(out)
#        out = convert_to_COARDS_area(out)
        out = out.drop(['lat_b','lon_b'])

        
        write_ncf(out,outfile.replace('nc','nc4'))
        # if convert divide by target area and change units attr for each variable


    out = xr.open_dataset(outfile.replace('nc','nc4'),decode_times=False,decode_cf=False)
    print('Changing Units')
    mw = get_species_mw()
    # print(mw)
    for v in out.data_vars:
        #print(v,out[v].attrs)
        if 'units' in out[v].attrs:
            if (out[v].attrs['units'].strip() == 'g/s') | (out[v].attrs['units'].strip() == 'moles/s'):
                if out[v].attrs['units'].strip() == 'g/s':
                    out[v].attrs['units'] = 'kg m-2 s-1' # convert_units(out[v])
                    out[v].data[:] = out[v].data / 1000.
                    print('Convert {} from "g/s" to "kg m-2 s-1"'.format(v))
                else:
                    if v in mw:
                        print('{} in species_weights.txt'.format(v))
                        out[v].attrs['units'] = 'kg m-2 s-1'
                        out[v].data[:] = out[v].data * mw[v] / 1000.
                        print('Convert {} from "moles/s" to "kg m-2 s-1"'.format(v))
                    else:
                        print(' Check For {}'.format(v))
                        out[v].attrs['units'] = convert_units(out[v])
                # strip blank spaces from long_name and var_desc
                out[v].attrs['long_name'] = out[v].attrs['long_name'].strip()
                out[v].attrs['var_desc'] = out[v].attrs['var_desc'].strip()
        #print(v,out[v].attrs)
    # ensure units on latitude longitude and area
    #    write_ncf(out,'thisisatest.nc')
    out['x'] = out.lon[0,:].astype('float32')
    out['y'] = out.lat[:,0].astype('float32')
    out['x'].attrs['units'] = 'degrees_east'
    out['y'].attrs['units'] = 'degrees_north'
    out['x'].attrs['long_name'] = 'longitude'
    out['y'].attrs['long_name'] = 'latitude'
    out = out.drop(['lat','lon'])
    out = out.rename({'y':'lat','x':'lon'})
    out['area'].attrs['units'] = 'm**2'
    out['area'].attrs['long_name'] = 'area'
    current_date = (pd.to_timedelta(out.time.values,unit='seconds') + pd.Timestamp('1970-01-01'))[0].strftime('hours since %Y-%m-%d 00:00:00')
    out['time'] = np.arange(0,25)
    out['time'].attrs['units'] = current_date
    out['time'].attrs['long_name'] = 'time'
    out['time'].attrs['calendar'] = 'standard'
    out.attrs['conventions'] = 'COARDS'
    out.attrs['history'] = 'processed using process_nei'
    out.attrs['title'] = 'GMU NEMO Emissions'
#    out = out.drop(['lat','lon'])

    # output final file : outfile
    write_ncf(out.isel(time=slice(0,24)),outfile)

    if verify:
        print(' Variable    Remapped (mol/s | g/s)   Orig (mol/s | g/s)')
        if 'area' not in out.data_vars:
            print('AREA NOT IN ----->', outfile, 'Adding')
            t = xr.open_dataset(target_file)
            out['area'] = t.area
        for v in out.data_vars:
            if (v != 'area') & (v != 'lat') & (v != 'lon') & (v != 'time'):
                # print(v)
                if v not in mw:
                    n = float((out[v].isel(time=0) * out.area.data * 1000).sum())
                else:
                    try:
                        #print(v, orig[v].attrs)
                        if orig[v].attrs['units'].strip() == 'mol/s':
                            n = float((out[v].isel(time=0) * out.area.data / mw[v] * 1000).sum())
                        else:
                            n = float((out[v].isel(time=0) * out.area.data * 1000).sum())
                    except AttributeError:
                        print('AREA NOT IN -> ',outfile)
                o = float(orig[v].isel(time=0).sum().item(0))
                print(' {}          {:.5f}          {:.5f}'.format(v,n,o))
                        
                        # remove intermediate file
    os.remove(outfile.replace('nc','nc4'))

if __name__ == '__main__':
    parser = ArgumentParser(description='Regrid CMAQ IOAPI to WSG84 grid', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--infile', help='input ioapi file name', type=str, required=True)
    parser.add_argument('-o', '--outfile', help='output file', default='./', type=str,required=False)
    parser.add_argument('-t', '--target_file', help='predefined target domain file', type=str, required=False)
    parser.add_argument('-tr', '--target_resolution', help='predefined target domain file', type=float, default=None, required=False)
    parser.add_argument('-w', '--weight_file', help='filename of precomputed weight file', type=str, required=False)
    parser.add_argument('-c', '--change_units', help='change units to m-2 s-1 and divide by target area',required=False)
    parser.add_argument('-s', '--source_area', help='source area file',required=False,type=str,default=None)
    parser.add_argument('-v', '--verify', help='verify',required=False,type=bool, default=True)
    args = parser.parse_args()

    infile = args.infile
    # Check if input file exists
    outfile = args.outfile
    if args.target_file is not None:
        target = args.target_file
    else:
        target = 'target.nc'
    if args.weight_file	is not None:
        weights = args.weight_file
    else:
        weights= 'weights.nc'
    #    print(args.source_area)
    if args.source_area is not None:
        area = xr.open_dataset(args.source_area)
    else:
        area = None
    if args.target_resolution is not None:
        res = args.target_resolution
    else:
        res = 0.125
    #    print(area)
    if os.path.exists(infile):
        if os.path.exists(outfile):
            f = xr.open_dataset(outfile)
            if len(f.time) == 25:
                print('fixing time:',outfile)
                os.remove(outfile)
                f = f.isel(time=slice(0,24))
                write_ncf(f,outfile)

            elif 'area' in f.data_vars:
                print('file, {}, already exists... skipping'.format(outfile))
            else:
                os.remove(outfile)
                process(infile,outfile,target,weights,area=area,verify=args.verify, target_res=res)
            f.close()
        else:
            process(infile,outfile,target,weights,area=area,verify=args.verify, target_res=res)
        print('===============================================')
        print('--------------------SUCCESS--------------------')
        print('===============================================')
    else:
        print('input file, {} does not exits. Exiting'.format(infile))
