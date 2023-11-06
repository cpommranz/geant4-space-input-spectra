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
import numpy as np

from . import cli

A_P = 1.7e5
A_HE = 1.0e4
GAMMA = 2.72
KAPPA_P = 4.64
KAPPA_HE_HCP = 3.26
EPSILON_0_P = 817  # MeV/nucl
EPSILON_0_HE_HCP = 576  # MeV/nucl
DELTA = 3.7


def kuznetsov_modulation_potential(particle: str, sunspot_number: int) -> float:
    """Calculate the modulation potential according to Formula 6 in (Kuznetsov et al., 2017).

    :param particle: 'helium' or 'proton'. Use 'helium' also for heavy charged particles (HCP).
    :param sunspot_number: Monthly smoothed sunspot number from https://www.sidc.be/silso/datafiles
    """
    if particle == 'proton':
        return EPSILON_0_P + (KAPPA_P * sunspot_number)
    elif particle == 'helium':
        return EPSILON_0_HE_HCP + (KAPPA_HE_HCP * sunspot_number)


def kuznetsov_modulation_function(energy: float, particle: str, sunspot_number: int) -> float:
    """Calculate the modulation potential according to Formula 4 in (Kuznetsov et al., 2017).

    :param energy: Energy in MeV/nucl.
    :param particle: 'helium' or 'proton'. Use 'helium' also for heavy charged particles (HCP).
    :param sunspot_number: Monthly smoothed sunspot number from https://www.sidc.be/silso/datafiles
    """
    psi = pow(energy / (energy + kuznetsov_modulation_potential(particle, sunspot_number)), DELTA)
    return psi


def kuznetsov_flux(energy: float, particle: str, sunspot_number: int, normalization_coefficient: float) -> float:
    """Calculate the Kuznetsov spectral flux at the given energy according to Formula 3 in (Kuznetsov et al., 2017).

    :param energy: Energy in MeV/nucl.
    :param particle: 'helium' or 'proton'. Use 'helium' also for heavy charged particles (HCP).
    :param sunspot_number: Monthly smoothed sunspot number from https://www.sidc.be/silso/datafiles
    :param normalization_coefficient: Normalization coefficient according to Table 2 in (Kuznetsov et al., 2017).
    """
    modulation = 1.
    if energy <= 20000:  # E<=20 GeV/nucl
        modulation = kuznetsov_modulation_function(energy, particle, sunspot_number)
    if particle == 'proton':
        a = A_P
    elif particle == 'helium':
        a = A_HE
    else:
        raise ValueError("Unknown particle.")
    flux = normalization_coefficient * a * energy**(-GAMMA) * modulation
    return flux


@cli.command()
@click.argument('ecsvfile', type=click.Path(writable=True))
@click.option('--number', '-n', type=int, default=1000, help="Number of data points. Default: 1000")
@click.option('--elo', '-e', type=float, default=80, help="Minimum energy in MeV/nucleon. Default: 80")
@click.option('--ehi', '-E', type=float, default=1e5, help="Maximum energy in MeV/nucleon. Default: 1e5")
@click.option('--particle', '-p', type=click.Choice(['proton', 'helium']), default='proton',
              help="Must be 'proton' or 'helium'. Use 'helium' for other heavy charged particles (HCP), too.")
@click.option('--sunspots', '-s', type=float, default=0,
              help="Monthly averaged sunspot number from https://www.sidc.be/silso/datafiles. Pay attention that for "
                   "odd solar cycles the sunspot number must be taken 15.5 months before the investigated time point, "
                   "for even solar cycles 5.5 months.")
@click.option('--norm', '-N', type=float, default=1.,
              help="Normalization coefficient according to Table 2 in (Kuznetsov et al., 2017).")
@click.option('--nucl', '-A', type=int, default=1, help="Number of nucleons.")
def kuznetsov2017(ecsvfile, number, elo, ehi, particle, sunspots, norm, nucl):
    xs = np.logspace(np.log10(elo), np.log10(ehi), number)
    fluxes = np.array([kuznetsov_flux(x, particle, sunspots, norm) / nucl for x in xs])
    t = Table(data=[xs * nucl, fluxes], names=['E', 'flux'])
    t.write(ecsvfile, format='ascii.ecsv')
