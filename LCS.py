#!/usr/bin/env python3.10
import argparse
import os
from generate_files import initiate_files as gen
from flow_map_all_members import flow_map as fl
from ftle import ftle as ftle

parser = argparse.ArgumentParser(description='LCS from ensemble')
group = parser.add_mutually_exclusive_group()

group.add_argument('-g','--generate', action='store_true', help='Generate .nc files.')
group.add_argument('-fl', '--flow_map', action='store_true', help='Show flow map.')
group.add_argument('-lcs', '--LCS', action='store_true', help='Calculate LCS field with FTLE')

#For --generate
parser.add_argument('-gr', '--grid_size', metavar='', type=int, help='Grid size.')
parser.add_argument('-lon', '--longitude', nargs=2, metavar='', type=float, help='longitude array.')
parser.add_argument('-lat', '--latitude', nargs=2, metavar='', type=float, help='latitude array.')
parser.add_argument('-t', '--time_step', metavar='', type=int, help='Time step.')
parser.add_argument('-of', '--o_file', metavar='', type=str, help='outfile name.')
parser.add_argument('-if', '--i_file', metavar='', type=str, help='input file.') #this is used as a general input file for other commands aswell

#For --lcs
parser.add_argument('-p', '--plot', action='store_true', help='Plot LCS field')

args = parser.parse_args()

def run_LCS():
    if args.i_file is None:
        parser.error('Input file is not specified')
    #Generate files
    if args.generate:
        if args.grid_size is None:
            parser.error('Grid Size is not specified')
        if args.longitude is None:
            parser.error('Longitude is not specified')
        if args.latitude is None:
            parser.error('Latitude is not specified')
        if args.time_step is None:
            parser.error('Time step is not specified')
        of = None
        if args.o_file is not None:
            of = args.o_file
        gen(args.grid_size, args.longitude, args.latitude, args.time_step, args.i_file, of)
    #Flow map
    if args.flow_map:
        fl(args.i_file)
    #LCS
    if args.LCS:
        ftle(args.i_file, args.plot)




if __name__ == '__main__':
    run_LCS()
