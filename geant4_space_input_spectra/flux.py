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
import numpy as np
import os
import scipy.integrate
import scipy.interpolate

from . import cli


@cli.command()
@click.argument('spectra', type=click.Path(exists=True), nargs=-1)
@click.option('--interpolate/--no-interpolate', default=True, help="Interpolate the spectrum.")
@click.option('--number', '-n', default=1000, help="Number of interpolated data points.")
@click.option('--iptype', '-t', type=click.Choice(choices=['log', 'lin'], case_sensitive=False), default='log',
              help="Type of interpolation. Default: log")
def flux(spectra, interpolate, number, iptype):
    for spectrum in spectra:
        data = ascii.read(spectrum)
        name = os.path.basename(spectrum).removesuffix('.ecsv')

        # From https://stackoverflow.com/questions/29346292/logarithmic-interpolation-in-python
        def log_interp1d(xx, yy, kind='linear'):
            logx = np.log10(xx)
            logy = np.log10(yy)
            lin_interp = scipy.interpolate.interp1d(logx, logy, kind=kind)
            log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
            return log_interp

        def lin_interp1d(xx, yy, kind='linear'):
            return scipy.interpolate.interp1d(xx, yy, kind=kind)

        energies = data['E']
        fluxes = data['flux']
        if interpolate:
            # Data points are always chosen in logspace
            energies = np.logspace(np.log10(data['E'][0]), np.log10(data['E'][-1]), number)
            if iptype == 'log':
                fluxes = log_interp1d(data['E'], data['flux'])(energies)
            elif iptype == 'lin':
                fluxes = lin_interp1d(data['E'], data['flux'])(energies)

        total_flux = scipy.integrate.simps(fluxes, energies)
        click.echo(f"Total flux of spectrum '{name}': {total_flux} particles / (cm² * s * sr) "
                   f"in [{data['E'][0]}, {data['E'][-1]}] MeV")
