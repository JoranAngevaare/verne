"""
verne.py - Code for calculating the Earth stopping effect, primarily for heavy Dark Matter.


Last updated: 13/12/2017
Contact: Bradley Kavanagh, bradkav@gmail.com

"""

import numba
import numpy as np
from scipy.integrate import quad, simps
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import odeint
import scipy.special
from verne import MaxwellBoltzmann as MB
import os

# --------------------
# Theta = 0 is directly from BELOW, angles in radians
# Gamma = 0 is mean DM flux directly from BELOW, angles in radians

# Densities are in number/cm^3
# Distances in m
# --------------------

print("*********************************************")
print("WARNING: SOME v^-4 FACTORS HAVE BEEN ADDED...")
print("*********************************************")

## Integrate options
MXSTEPS = 15
RTOL = 1.e-1

# MXSTEPS=1000, RTOL=1e-6
print(f'Be arware. The integrateion parameters are set to MXSTEPS, RTOL =  {MXSTEPS}, {RTOL}')

isotopes = None
dens_profiles = None
dens_interp = None
Avals = None
Niso = None
Niso_full = None

phi_interp = None
corr_interp = None
corr_Pb = None
corr_Cu = None

NEGLECT_FF = False

if (NEGLECT_FF):
    print(" ")
    print("*********************************************")
    print("WARNING: NEGLECTING FORM FACTORS...")
    print("*********************************************")

isoID = {"O": 0, "Si": 1, "Mg": 2, "Fe": 3, "Ca": 4, "Na": 5, "S": 6, "Al": 7, "O_A": 8, "N_A": 9}

h_A = 80e3  # Height of atmosphere in (m)
R_E = 6371.0e3  # Earth Radius in (m)


def loadIsotopes(path='.'):
    print("    Loading isotope data and density profiles...")

    global dens_profiles
    global isotopes
    global Avals
    global dens_interp
    global Niso
    global Niso_full

    rootdir = os.path.join(path, "..", "data")

    # Load in Earth isotopes
    Avals = np.loadtxt(os.path.join(rootdir, "isotopes.txt"), usecols=(1,))
    isotopes = np.loadtxt(os.path.join(rootdir, "isotopes.txt"), usecols=(0,))
    Niso = len(isotopes)  # Number of #arth isotopes
    Niso_full = Niso + 2  # Plus isotopes in atmosphere

    # Load density profiles for Earth isotopes
    r_list = np.loadtxt(os.path.join(rootdir, "dens_profiles", "n_1.dat"), usecols=(0,))
    r_list0 = 0.0 * r_list
    dens_profiles = [np.loadtxt(os.path.join(rootdir, "dens_profiles", "n_" + str(int(iso)) + ".dat"), usecols=(1,))
                     for iso in isotopes]

    # Grid of heights in the atmosphere
    h_list = np.linspace(0, h_A, 100) + 1e-5

    # Make a slight correction of the profile - truncate at 6371 km
    # This doesn't affect things on large scales, but we have to make
    # sure we've truncated the Earth at R_E exactly if we have very shallow
    # detectors
    r_list[-1] = R_E

    # Append atmospheric points to the list of radii and densities for the Earth
    r_list = np.append(r_list, R_E + h_list)
    for i, dens in enumerate(dens_profiles):
        dens[-1] = dens[-2]
        dens_profiles[i] = np.append(dens, 0.0 * h_list)

    # Add the atmospheric elements:
    Avals = np.append(Avals, [16, 14])

    # Load atmospheric parameters and
    # calculate the atmosphere density profiles...
    Hvals, Tvals, beta = np.loadtxt(os.path.join(rootdir, "ISA.txt"), unpack=True)
    Hvals *= 1e3
    beta *= 1e-3  # Get everything in m

    # Fraction of Oxygen and Nitrogen
    frac = [0.21, 0.78]
    dens = lambda x: atmos_density(x, Hvals, Tvals, beta)
    dens_vec = np.vectorize(dens)
    dens_atmos = [np.append(r_list0, 2 * f_n * dens_vec(h_list)) for f_n in frac]
    dens_profiles.extend(dens_atmos)

    # Generate interpolation functions for the density profiles
    dens_interp = [interp1d(r_list, dens, bounds_error=False, fill_value=0) for dens in
                   dens_profiles]


# Generate interpolation functions for the Form Factor corrections (C_i(v))
def loadFFcorrections(m_x, path='.'):
    global corr_interp
    global corr_Pb
    global corr_Cu

    # Check that the isotope list has been loaded
    if (Avals is None):
        loadIsotopes(path)

    print("    Calculating Form Factor corrections for m_x = ", m_x, " GeV...")
    corr_interp = [calcFFcorrection(m_x, Avals[ID]) for ID in range(Niso_full)]
    # Also need Lead + Copper, for the shielding
    corr_Pb = calcFFcorrection(m_x, 207)
    corr_Cu = calcFFcorrection(m_x, 63.5)


# International standard atmosphere, ISO 2533:1975
# https://www.iso.org/standard/7472.html
def atmos_density(h, Hvals, Tvals, beta):
    H = R_E * h / (R_E + h)
    if (H > 80000.0):
        return 0
    R = 287.05
    p0 = 1.01e5
    g = 9.807

    # determine the layer
    ib = np.digitize(H, Hvals, right=False) - 1

    if (ib == 0):
        return 0
    if (ib == 1):
        p = p0 * (1.0 + beta[1] * H / Tvals[1]) ** (-g / (beta[1] * R))
    else:
        p = p0 * (1.0 + beta[1] * (Hvals[2]) / Tvals[1]) ** (-g / (beta[1] * R))
        for i in range(2, ib):
            if (beta[i] < 1e-3):
                p *= np.exp(-g * (Hvals[i + 1] - Hvals[i]) / (Tvals[i] * R))
            else:
                p *= (1.0 + beta[i] * (Hvals[i + 1] - Hvals[i]) / Tvals[i]) ** (-g / (beta[i] * R))
        if (beta[ib] < 1e-3):
            p *= np.exp(-g * (H - Hvals[ib]) / (Tvals[ib] * R))
        else:
            p *= (1.0 + beta[ib] * (H - Hvals[ib]) / Tvals[ib]) ** (-g / (beta[ib] * R))
    n = 1e-6 * 6.022e23 * p / (8314.32e-3 * Tvals[ib])  # Air particles per cubic cm

    return n


# Path length [m], as measured from the top of the atmosphere to the detector
# (at 'depth' m underground)
@numba.jit()
def pathLength(depth, theta):
    r_det = R_E - depth
    return +np.cos(theta) * r_det + np.sqrt(
        (-np.cos(theta) * r_det) ** 2 - (r_det ** 2 - (R_E + h_A) ** 2))


# Path length [m], as measured from the Earth's surface to the detector
# (at 'depth' m underground)
def pathLength_Earth(depth, theta):
    r_det = R_E - depth
    return +np.cos(theta) * r_det + np.sqrt(
        (-np.cos(theta) * r_det) ** 2 - (r_det ** 2 - (R_E) ** 2))


# Calculate the Form Factor correction for a nucleus of mass-number A0
# See Eq. (10) of the paper

# Maximum recoil energy (in keV)
def ERmax(mX, mA, v):
    mu = mX * mA * 1.0 / (mX + mA)
    return (1e6 / (3e5 * 3e5)) * 2 * (mu * v) ** 2 / mA


# Calculate the spin-independent form factor
# for nucleon number A0 and recoil energy E
def calcSIFormFactor(E, A0):
    # Helm
    if (E < 1e-5):
        return 1.0

    # Define conversion factor from amu-->keV
    amu = 931.5 * 1e3

    # Convert recoil energy to momentum transfer q in keV
    q1 = np.sqrt(2 * A0 * amu * E)

    # Convert q into fm^-1
    q2 = q1 * (1e-12 / 1.97e-7)

    # Calculate nuclear parameters
    s = 0.9
    a = 0.52
    c = 1.23 * (A0 ** (1.0 / 3.0)) - 0.60
    R1 = np.sqrt(c * c + 7 * np.pi * np.pi * a * a / 3.0 - 5 * s * s)

    x = q2 * R1
    J1 = np.sin(x) / x ** 2 - np.cos(x) / x
    F = 3 * J1 / x
    return (F ** 2) * (np.exp(-(q2 * s) ** 2))


def calcFFcorrection(m_x, A0):
    v_vals = np.linspace(0, 1000, 200)
    corr_fact = v_vals * 0.0
    for i, v in enumerate(v_vals):
        corr_fact[i] = \
        quad(lambda x: 2.0 * x * calcSIFormFactor(x * ERmax(m_x, 0.9315 * A0, v), A0), 0, 1)[0]
    corr_fact[0] = 1.0

    if (NEGLECT_FF):
        return interp1d(v_vals, 1.0 + 0.0 * corr_fact, kind='linear', bounds_error=False,
                        fill_value=0.0)
    else:
        return interp1d(v_vals, corr_fact, kind='linear', bounds_error=False, fill_value=0.0)


# Calculate the DM-nucleus 'effective' cross section
# which takes into account the average energy loss
def effectiveXS(sigma_p, m_X, A, v=1.0):
    m_p = 0.9315  # Proton mass
    m_A = 0.9315 * A
    mu_A = m_A * m_X / (m_A + m_X)
    mu_p = m_p * m_X / (m_p + m_X)
    # (v**-4)
    return sigma_p * (1.0 / (m_X * m_A)) * (A ** 2) * (mu_A ** 4 / mu_p ** 2)


# Calculate the final speed distribution at the detector
def CalcF(vf, gamma, depth, sigma_p, m_x, target, vmax_interp):
    # Define a grid of values for theta which we sample over
    # theta = pi/2 is often problematic, so we sample more densely there
    tlist = np.linspace(0, np.pi, 51)
    tlist = np.append(tlist, (np.pi / 2) * (1 + np.logspace(-3, -0.01, 25)))
    tlist = np.append(tlist, (np.pi / 2) * (1 - np.logspace(-3, -0.01, 25)))
    tlist = np.sort(tlist)

    fint = tlist * 0.0
    for i in range(len(tlist)):
        # If maximum vf you can get this value of theta is greater
        # than the speed we're interested in, set to zero
        if (vmax_interp(tlist[i]) < vf):
            fint[i] = 0.0
        else:
            fint[i] = f_integrand_full(vf, tlist[i], gamma, depth, sigma_p, m_x, target)

    # Integrate with Simpson's rule
    return simps(fint, tlist)


# Integrand for calculating the final speed distribution at the detector
def f_integrand_full(vf, theta, gamma, depth, sigma_p, m_x, target):
    # Calculate the initial velocity corresponding to this final velocity vf
    dv = 1.5
    vi1 = calcVinitial_full(vf + dv / 2.0, theta, depth, sigma_p, m_x, target)
    vi2 = calcVinitial_full(vf - dv / 2.0, theta, depth, sigma_p, m_x, target)

    # Calculate the average and the numerical derivative
    vi = (vi1 + vi2) / 2.0
    dvi_by_dvf = np.abs(vi1 - vi2) * 1.0 / dv

    return (dvi_by_dvf) * np.sin(theta) * (vi ** 2) * MB.calcf_integ(vi, theta, gamma)


# Calculate the distance of a point from the centre of the Earth
# The point is defined by:
#   - theta, the angle of the trajectory
#   - depth,the detector depth
#   - D, the distance along the trajectory, starting at the top of the atmosphere
# @numba.jit()
def radius(D, theta, depth):
    # print(type(D), type(theta), type(depth))
    r_det = R_E - depth
    return np.sqrt(
        (R_E + h_A) ** 2 + D ** 2 + 2 * D * (r_det * np.cos(theta) - pathLength(depth, theta)))


# Derivative of DM speed along path length D
# To be used by the ODE integrator
# @numba.jit(nopython = True)
def dv_by_dD(v, D, params):
    theta, depth, sigma_p, m_x, target = params
    res = 0.0
    isovals = []
    if (target == "atmos"):
        isovals = [8, 9]
    elif (target == "earth"):
        isovals = range(Niso)
        # for i in range(Niso):
        #    isovals.append(i)
    else:
        isovals = range(Niso_full)
        # for i in range(Niso_full):
        #     isovals.append(i)

    r = radius(D, theta, depth)

    # Loop over the relevant isotopes
    # BJK!
    for i in isovals:
        res += dens_interp[i](r) * effectiveXS(sigma_p, m_x, Avals[i], v=v) * corr_interp[i](v)
    return -1e2 * v * res  # (km/s)/m


# Derivative for the case of Pb shielding
def dv_by_dD_Pb(v, D, params):
    # Pb density
    n_Pb = 3.3e22
    A_Pb = 207

    sigma_p, m_x = params
    res = n_Pb * effectiveXS(sigma_p, m_x, A_Pb, v=v) * corr_Pb(v)
    return -1e2 * v * res  # (km/s)/m


# Derivative for the case of Cu shielding
def dv_by_dD_Cu(v, D, params):
    # Cu density
    n_Cu = 8.5e22
    A_Cu = 63.5

    sigma_p, m_x = params
    res = n_Cu * effectiveXS(sigma_p, m_x, A_Cu, v=v) * corr_Cu(v)
    return -1e2 * v * res  # (km/s)/m


# Calculate the final velocity after propagating across 'target'
# Here, target = "atmos" or "earth"
def calcVfinal(vi, theta, depth, sigma_p, m_x, target="full"):
    params = [theta, depth, sigma_p, m_x, target]

    # Propagate across the atmosphere
    if (target == "atmos"):
        d1 = 0
        d2 = pathLength(depth, theta) - pathLength_Earth(depth, theta)

    # Propagate from the surface of the Earth to the detector
    if (target == "earth"):
        d1 = pathLength(depth, theta) - pathLength_Earth(depth, theta)
        d2 = pathLength(depth, theta)

    psoln = odeint(dv_by_dD, vi, [d1, d2], args=(params,), mxstep=MXSTEPS, rtol=RTOL)
    vf = psoln[1]
    return vf


# Calculate the final velocity after propagating from the top of the
# atmosphere to the detector, account for all the steps
# Recommend using target="MPI" or "SUF" depending on the detector
def calcVfinal_full(vi, theta, depth, sigma_p, m_x, target="full"):
    vf = 1.0 * vi
    if (target == "XENON" or target == "SNOLAB"):
        # No need to calculate any other contribution that the earth shielding effect at 1400m depth
        return calcVfinal(vf, theta, depth, sigma_p, m_x, target="earth")
    if (target in ["atmos", "full", "no_shield", "SUF", "MPI", "EDE"]):
        vf = calcVfinal(vf, theta, depth, sigma_p, m_x, target="atmos")
    if (target in ["earth", "full", "no_shield", "SUF", "MPI", "EDE"]):
        vf = calcVfinal(vf, theta, depth, sigma_p, m_x, target="earth")
    if (target == "MPI"):
        vf = calcVfinal_shield_MPI(vf, sigma_p, m_x)
    if (target == "SUF"):
        vf = calcVfinal_shield_SUF(vf, sigma_p, m_x)
    if (target == "EDE"):
        vf = calcVfinal_shield_EDE(vf, sigma_p, m_x)
    return vf


# Calculate the initial velocity (for a given final velocity) after propagating across 'target'
# Here, target = "atmos" or "earth"
def calcVinitial(vf, theta, depth, sigma_p, m_x, target="earth"):
    params = [theta, depth, sigma_p, m_x, target]

    # Propagate across the atmosphere
    if (target == "atmos"):
        d1 = 0
        d2 = pathLength(depth, theta) - pathLength_Earth(depth, theta)

    # Propagate from the surface of the Earth to the detector
    if (target == "earth"):
        d1 = pathLength(depth, theta) - pathLength_Earth(depth, theta)
        d2 = pathLength(depth, theta)

    psoln = odeint(dv_by_dD, vf, [d2, d1], args=(params,), mxstep=MXSTEPS, rtol=RTOL)
    return psoln[1]


# Calculate the initial speed at the top of the atmosphere for a
# given final speed at the detector
# Recommend using target="MPI" or "SUF" depending on the detector
def calcVinitial_full(vf, theta, depth, sigma_p, m_x, target="full"):
    vi = 1.0 * vf
    if (target == "XENON" or target == "SNOLAB"):
        # No need to calculate any other contribution that the earth shielding effect at 1400m depth
        return calcVinitial(vi, theta, depth, sigma_p, m_x, target="earth")
    if (target == "MPI"):
        vi = calcVinitial_shield_MPI(vi, sigma_p, m_x)
    if (target == "SUF"):
        vi = calcVinitial_shield_SUF(vi, sigma_p, m_x)
    if (target == "EDE"):
        vi = calcVinitial_shield_EDE(vi, sigma_p, m_x)
    if (target in ["earth", "full", "no_shield", "SUF", "MPI", "EDE", "XENON"]):
        vi = calcVinitial(vi, theta, depth, sigma_p, m_x, target="earth")
    if (target in ["atmos", "full", "no_shield", "SUF", "MPI", "EDE", "XENON"]):
        vi = calcVinitial(vi, theta, depth, sigma_p, m_x, target="atmos")

    return vi


# Calculate final (or initial) speed after crossing the Lead shielding at SUF
def calcVfinal_shield_SUF(v0, sigma_p, m_x):
    params = [sigma_p, m_x]
    # Propagate through 16cm of Lead
    psoln = odeint(dv_by_dD_Pb, v0, [0, 16.0e-2], args=(params,), mxstep=MXSTEPS, rtol=RTOL)
    return psoln[1]


def calcVinitial_shield_SUF(v0, sigma_p, m_x):
    params = [sigma_p, m_x]
    # Propagate through 16cm of Lead (backwards)
    psoln = odeint(dv_by_dD_Pb, v0, [16.0e-2, 0], args=(params,), mxstep=MXSTEPS, rtol=RTOL)
    return psoln[1]


# Calculate final (or initial) speed after crossing the Copper shielding at MPI
def calcVfinal_shield_MPI(v0, sigma_p, m_x):
    params = [sigma_p, m_x]
    # Propagate through 1mm Copper
    psoln = odeint(dv_by_dD_Cu, v0, [0, 1e-3], args=(params,), mxstep=MXSTEPS, rtol=RTOL)
    return psoln[1]


def calcVinitial_shield_MPI(v0, sigma_p, m_x):
    params = [sigma_p, m_x]
    # Propagate through 1mm Copper
    psoln = odeint(dv_by_dD_Cu, v0, [1e-3, 0], args=(params,), mxstep=MXSTEPS, rtol=RTOL)
    return psoln[1]


# Calculate final (or initial) speed after crossing the Lead shielding at EDE...
def calcVfinal_shield_EDE(v0, sigma_p, m_x):
    params = [sigma_p, m_x]
    # Propagate through 16cm of Lead
    psoln = odeint(dv_by_dD_Pb, v0, [0, 10.0e-2], args=(params,), mxstep=MXSTEPS, rtol=RTOL)
    return psoln[1]


def calcVinitial_shield_EDE(v0, sigma_p, m_x):
    params = [sigma_p, m_x]
    # Propagate through 16cm of Lead (backwards)
    psoln = odeint(dv_by_dD_Pb, v0, [10.0e-2, 0], args=(params,), mxstep=MXSTEPS, rtol=RTOL)
    return psoln[1]
