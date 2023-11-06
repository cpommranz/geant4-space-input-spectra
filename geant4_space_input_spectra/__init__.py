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


import click


@click.group()
def cli():
    """Convert and create input spectra for radiation environments in space for the Geant4 GPS."""
    pass


from . import convert
from . import flux
from . import gruber_1999
from . import kuznetsov_2017
from . import plot_spectra
