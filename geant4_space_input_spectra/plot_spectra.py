# geant4-space-input-spectra – Convert and create input spectra for the Geant4 GPS.
# Copyright (C) 2023  Christian Pommranz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


from astropy.io import ascii
import click
import itertools
import matplotlib.pyplot as plt
import os

from . import cli

# Use fixed markers for now.
markers = itertools.cycle(['.', '+', 'x', '1', '*', 'v', 's', '<', 'd', '>', '^'])


@cli.command()
@click.argument('outputfile', nargs=1)
@click.argument('spectra', type=click.Path(exists=True), nargs=-1)
@click.option('--particle', '-p', default='particles', help="Name of the particles for the y-axis.")
@click.option('--xlog/--no-xlog', default=True, help="Logarithmic xscale.")
@click.option('--ylog/--no-ylog', default=True, help="Logarithmic yscale.")
@click.option('--labels', '-l', default=None, help="Comma-separated list of labels.")
@click.option('--figsize', nargs=2, type=float, default=(6.4, 4.8), help="Size of Figure. Default: 6.4 4.8")
@click.option('--style', default=None, help="Custom matplotlib stylesheet.")
def plot_spectra(outputfile, spectra, particle, xlog, ylog, labels, figsize, style):
    if style is not None:
        plt.style.use(style)

    _, ax1 = plt.subplots(figsize=figsize, nrows=1, ncols=1)
    if xlog:
        ax1.set_xscale('log')
    if ylog:
        ax1.set_yscale('log')
    ax1.set_xlabel('Energy [MeV]')
    ax1.set_ylabel(f'Flux [{particle} / (cm² * s * sr * MeV)]')

    if labels is None:
        labels = [os.path.basename(x).removesuffix('.ecsv') for x in spectra]
    else:
        labels = labels.split(',')
    if len(labels) != len(spectra):
        click.echo("Number of labels must match number of spectra! Exiting.")
        exit(-1)
    for spectrum, label in zip(spectra, labels):
        data = ascii.read(spectrum)
        if 'yerrtot_lo' in data.colnames and 'yerrtot_hi' in data.colnames:
            ax1.errorbar(data['E'], data['flux'], yerr=[data['yerrtot_lo'], data['yerrtot_hi']],
                         label=label, marker=next(markers), linestyle='', capsize=2, elinewidth=1)
        else:
            ax1.plot(data['E'], data['flux'], label=label, marker=next(markers), linestyle='')

    ax1.legend()
    if outputfile == "-":
        plt.show()
    else:
        plt.savefig(outputfile, bbox_inches="tight")
