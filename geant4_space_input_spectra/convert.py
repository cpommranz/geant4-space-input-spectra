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

from . import cli


@cli.group()
def convert():
    pass


@convert.command()
@click.argument('datfile', type=click.Path(exists=True))
@click.argument('ecsvfile', type=click.Path(writable=True))
def crdb_dat_to_ecsv(datfile, ecsvfile):
    t = ascii.read(datfile, names=['Qty', 'E', 'E_lo', 'E_up', 'flux', 'ystat_lo', 'ystat_up',
                                   'ysyst_lo', 'ysyst_up', 'yerrtot_lo', 'yerrtot_hi'])
    t.sort('E')

    # GeV to MeV for energies
    for field in ['E', 'E_lo', 'E_up']:
        t[field] *= 1000

    # Conversion from counts / (GeV * m² * s * sr) to counts / (MeV * cm² * s * sr) -> 1 / 10000000
    for field in ['flux', 'ystat_lo', 'ystat_up', 'ysyst_lo', 'ysyst_up', 'yerrtot_lo', 'yerrtot_hi']:
        t[field] *= 1E-7

    t.write(ecsvfile, format='ascii.ecsv')


@convert.command()
@click.argument('txtfile', type=click.Path(exists=True))
@click.argument('nucleons', type=int)
@click.argument('ecsvfile', type=click.Path(writable=True))
def spenvis_txt_to_ecsv(txtfile, nucleons, ecsvfile):
    """Converts Spenvis output as a three-column txt file to ecsv.

    TXTFILE is the textfile containing a copy-and-paste version of the three columns Energy, IFlux, and DFlux from the
    Spenvis HTML output. Without headers.

    NUCLEONS is the number of nucleons. This is needed to convert from MeV/n to MeV.

    ECSVFILE is the name of the file to write to.
    """
    t = ascii.read(txtfile, names=['E', 'IFlux', 'flux'])
    t.sort('E')

    # MeV/n to MeV for energies
    t['E'] *= nucleons

    # Conversion from counts / (MeV/n * m² * s * sr) to counts / (MeV * cm² * s * sr)
    t['flux'] /= (nucleons * 10000)

    t.write(ecsvfile, format='ascii.ecsv')


@convert.command()
@click.argument('ecsvfile', type=click.Path(exists=True))
@click.argument('macfile', type=click.File('w'))
@click.option('--particle', '-p', default='geantino', help="Particle name.")
@click.option('--norm/--no-norm', default=True, help="Normalize spectrum.")
def ecsv_to_mac(ecsvfile, macfile, particle, norm):
    t = ascii.read(ecsvfile)

    if norm:
        total = sum(t['flux'])
        t['flux'] /= total

    macfile.write(f'/gps/particle {particle}\n\n')
    macfile.write("/gps/ene/type Arb\n/gps/hist/type arb\n\n")
    for row in t:
        macfile.write(f"/gps/hist/point    {row['E']:.4e}    {row['flux']:.4e}\n")
    macfile.write("/gps/hist/inter Log\n")
