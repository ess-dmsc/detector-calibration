#!/usr/bin/env python3

from scipp import array, DataArray, ones_like
from argparse import ArgumentParser
import h5py, os
import matplotlib.pyplot as plt


# Convert Ring and FEN to numbers or if not set, to 'any'
def id2chr(id):
    if id == -1:
        return 'any'
    else:
        return f'{id}'


# Read EFU H5 file into scipp data structure. Adding positions and straw (as
# float values) calculated from the four amplitudes A,B,D and D to the structure
def readtoscipp(filename):

    f = h5py.File(filename, 'r')
    dat =  f['loki_readouts']

    tube = array(values=dat['TubeId'].astype('int'), dims=['event'])
    ring = array(values=dat['RingId'].astype('int'), dims=['event'])
    fen = array(values=dat['FENId'].astype('int'), dims=['event'])

    ampl_a = array(values=1.0 * dat['AmpA'].astype('int'), dims=['event'], unit='mV')
    ampl_b = array(values=1.0 * dat['AmpB'].astype('int'), dims=['event'], unit='mV')
    ampl_c = array(values=1.0 * dat['AmpC'].astype('int'), dims=['event'], unit='mV')
    ampl_d = array(values=1.0 * dat['AmpD'].astype('int'), dims=['event'], unit='mV')

    events = ones_like(1. * tube)
    events.unit = 'counts'

    pos = (ampl_a + ampl_b) / (ampl_a + ampl_b + ampl_c + ampl_d)
    straw = (ampl_b + ampl_d) / (ampl_a + ampl_b + ampl_c + ampl_d)

    return DataArray(data=events,
            coords={'pos': pos, 'straw': straw, # 'time': time,
                    'tube': tube, 'ring': ring, 'fen': fen,
                    'amplitude_a': ampl_a, 'amplitude_b': ampl_b,
                    'amplitude_c': ampl_c, 'amplitude_d': ampl_d})


# This is the 'main' program entry
def load_and_process_data(args):
    dat = readtoscipp(args.filename)

    rgrp = array(dims=['ring'], values=[args.ring])
    fgrp = array(dims=['fen'], values=[args.fen])

    fig, ax = plt.subplots(4,2, figsize=(16,16))

    for i in range(8): # 8 is the number of tubes on a FEN
        print(f'processing ring {id2chr(args.ring)}, fen {id2chr(args.fen)}, tube {i}')
        tgrp = array(dims=['tube'], values=[i])
        if args.ring == -1 and args.fen == -1:
            grp = dat.group(tgrp).bins.concat()
        elif args.ring == -1 and args.fen != -1:
            grp = dat.group(fgrp, tgrp).bins.concat()
        elif args.ring != -1 and args.fen == -1:
            grp = dat.group(rgrp, tgrp).bins.concat()
        else:
            grp = dat.group(rgrp, fgrp, tgrp).bins.concat()

        yi = i // 2
        xi = i % 2
        cax = ax[yi, xi]
        grp.hist(pos=200, straw=200).plot(aspect=1.,norm='log', ax=cax)
        cax.title.set_text(f'Tube {i}')
        cax.set_xlim(0, 1)
        cax.set_ylim(0, 1)
        cax.yaxis.tick_left()
        cax.yaxis.set_label_position('left')
        if i <= 5:
            cax.set(xlabel='', ylabel='pos')
        else:
            cax.set(xlabel='straw', ylabel='pos')

    plt.suptitle(f'Ring: {id2chr(args.ring)}, FEN: {id2chr(args.fen)}, Tubes 0 - 8', size='28')
    plt.savefig(os.path.join(args.outdir, f'strawpos_{id2chr(args.ring)}_{id2chr(args.fen)}.png'))



if __name__ == '__main__':
    parser = ArgumentParser(prog='dattoplot', description=__doc__)
    parser.add_argument('filename', type=str, nargs='?', default="",
                        help='.h5 file to load and plot')
    parser.add_argument('-o','--outdir', type=str, default="",
                        help='output directory')
    parser.add_argument('-r','--ring', type=int, default=-1, help='Ring Id (default all rings)')
    parser.add_argument('-f','--fen', type=int, default=-1, help='FEN Id (default all fens)')

    args = parser.parse_args()

    load_and_process_data(args)
