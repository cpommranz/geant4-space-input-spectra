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


from astropy.table import Table
import click
from math import exp
import numpy as np

from . import cli

def gruber_spectrum(energy: float) -> float:
    """Spectrum from Gruber et al., 1999
    Note that energy is in keV!"""
    if energy <= 60:
        return 7.877 * pow(energy, -0.29) * exp(-energy / 41.13)
    return 0.0259 * pow(energy/60, -5.5) + \
        0.504 * pow(energy/60, -1.58) + \
        0.0288 * pow(energy/60, -1.05)


def gruber_flux(energy: float) -> float:
    """Calculate the flux in counts / (MeV * cm^2 * s * sr)
    Note that energy is expected MeV."""
    return 1000 * gruber_spectrum(energy * 1000) / (energy * 1000)


@cli.command()
@click.argument('ecsvfile', type=click.Path(writable=True))
@click.option('--number', '-n', type=int, default=1000, help="Number of data points. Default: 1000")
@click.option('--elo', '-e', type=float, default=3e-3, help="Minimum energy in MeV. Default: 3e-3")
@click.option('--ehi', '-E', type=float, default=1e5, help="Maximum energy in MeV. Default: 1e5")
def gruber1999(ecsvfile, number, elo, ehi):
    xs = np.logspace(np.log10(elo), np.log10(ehi), number)
    fluxes = np.array([gruber_flux(x) for x in xs])
    t = Table(data=[xs, fluxes], names=['E', 'flux'])
    t.write(ecsvfile, format='ascii.ecsv')
