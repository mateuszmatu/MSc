#!/usr/bin/env python3.10
import argparse
import os
from ftle import FTLE
import matplotlib.pyplot as plt
import flow_map as fm

parser = argparse.ArgumentParser(description='LCS from ensemble')
group = parser.add_mutually_exclusive_group()

group.add_argument('-g','--generate', action='store_true', help='Generate .nc files.')
group.add_argument('-fl', '--flow_map', action='store_true', help='Show flow map.')
group.add_argument('-lcs', '--LCS', action='store_true', help='Calculate LCS field with FTLE')
group.add_argument('-dg', '--double_gyre', action='store_true', help='Run double gyre example')

#For --generate
parser.add_argument('-lon', '--longitude', nargs=2, metavar='', type=float, help='Longitude array')
parser.add_argument('-lat', '--latitude', nargs=2, metavar='', type=float, help='Latitude array')
parser.add_argument('-t', '--time_step', metavar='', type=int, help='Time step in seconds')
parser.add_argument('-del', '--delta', metavar='', type=float, help='Distance between seeded particles in x and y')
parser.add_argument('-dur', '--duration', metavar='', type=int, help='Duration of simulation')
parser.add_argument('-if', '--in_file', metavar='', type=str, help='Input file or threeds URL')
parser.add_argument('-of', '--out_file', metavar='', type=str, help='Outfile name')

#For --lcs
parser.add_argument('-path', '--path', metavar='', type=str, help='Path to files')
parser.add_argument('-p', '--plot', action='store_true', help='Plot LCS field.')

#For --double_gyre
parser.add_argument('-pf', '--particle_flow', action='store_true', help='Plot particle positions over LCS field')
args = parser.parse_args()

def run_LCS():
    if args.generate:
        if not args.longitude or len(args.longitude) != 2:
            parser.error('Longitude is not specified or shape is wrong. Shape shoud be 2.')
        if not args.latitude or len(args.latitude) != 2:
            parser.error('Latitude is not specified or shape is wrong. Shape shoud be 2.')
        if not args.time_step:
            parser.error('Time step is not specified')
        if not args.delta:
            parser.error('Delta is not specified')
        if not args.duration:
            parser.error('Duration is not specified')
        if not args.in_file:
            parser.error('Input file or URL is not specified')
        
        from generate_files import initiate_files
        initiate_files(args.longitude, args.latitude, args.time_step, args.delta, args.duration, args.in_file, args.out_file)
    
    if args.LCS:
        if not args.path:
            parser.error('Path to files is not specified')
        

    #Double gyre example
    if args.double_gyre:
        req_file1 = 'flow_maps/double_gyre/double_gyre_f.nc'
        req_file2 = 'flow_maps/double_gyre/double_gyre_b.nc'
        req_file3 = 'flow_maps/double_gyre/initial_values.npy'

        if not os.path.exists(req_file1) or not os.path.exists(req_file2) or not os.path.exists(req_file3):
            from example_double_gyre import example_double_gyre 
            example_double_gyre()
        vals = 'flow_maps/double_gyre/initial_values.npy'
        repelling = FTLE('flow_maps/double_gyre/double_gyre_f.nc', vals)
        attracting = FTLE('flow_maps/double_gyre/double_gyre_b.nc', vals)
        #attracting = attracting[::-1,::-1] #flipping for some reason

        fig, ax = plt.subplots(2)
        a = ax[0].imshow(repelling, interpolation='nearest', origin='lower', cmap='jet')
        b = ax[1].imshow(attracting, interpolation='nearest', origin='lower', cmap='jet')
        ax[0].set(title='Repelling LCS')
        ax[1].set(title='Attracting LCS')
        l = [a,b]
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        for i in range(2):
            divider = make_axes_locatable(ax[i])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(l[i], cax=cax)
        plt.show()

if __name__ == '__main__':
    run_LCS()