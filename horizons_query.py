"""
This script queries NASA JPL Horizons for orbital elements of selected asteroids
and saves the results to a file called `Particles.el`.

It is intended for educational and research use, and the output format is suitable
for use in orbital simulations or post-processing tools.
"""

# Import the Horizons query interface from astroquery
from astroquery.jplhorizons import Horizons

# Import astropy utilities for working with tables of astronomical data
import astropy.table


# --------------------------------------------------------------------
# 1. Define the list of asteroid identifiers
# --------------------------------------------------------------------
# These identifiers are recognized by the JPL Horizons system.
# They may refer to provisional designations or named objects.
asteroid_ids = ['2025 PN7', '2025 SC', '2025 RB7']


# --------------------------------------------------------------------
# 2. Query JPL Horizons for orbital elements
# --------------------------------------------------------------------
# For each asteroid:
#   - `id` specifies the asteroid designation
#   - `location='500@10'` specifies Earth as the observer
#   - `epochs=2460800.5` is the Julian Date at which elements are requested
#
# The `.elements()` method returns a table of Keplerian orbital elements.
# All individual tables are stacked into one combined table.
all_results = astropy.table.vstack([
    Horizons(id=i, location='500@10', epochs=2460800.5).elements()
    for i in asteroid_ids
])


# --------------------------------------------------------------------
# 3. Write results to an output file
# --------------------------------------------------------------------
# The output file `Particles.el` will contain:
#   - The number of asteroids (first line)
#   - One line per asteroid with orbital elements in a fixed numeric format
#
# This replaces printing to the terminal, making the results reusable
# by other programs or simulations.
with open("Particles.el", "w") as f:

    # Write the total number of asteroids
    f.write(f"{len(asteroid_ids)}\n")

    # Loop over each asteroid's orbital elements
    for row in all_results:

        # Target name or designation (not written, but retained for clarity)
        designation = row['targetname']

        # Write the orbital elements in the following order:
        #   a      = semi-major axis (AU)
        #   e      = eccentricity
        #   incl   = inclination (degrees)
        #   Omega  = longitude of ascending node (degrees)
        #   w      = argument of perihelion (degrees)
        #   M      = mean anomaly (degrees)
        #
        # Numeric precision is chosen to preserve accuracy for simulations.
        f.write(
            f"{row['a']:.15f} "
            f"{row['e']:.16f} "
            f"{row['incl']:.14f} "
            f"{row['Omega']:.15f} "
            f"{row['w']:.15f} "
            f"{row['M']:.15f}\n"
        )
