import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp2d
import scipy.special
import os

phi_interp = None


# ------- Velocity distribution stuff----------
# ----------------------------------------------


# Load an interpolation function for the integral
# of the Maxwell-Boltzmann velocity distribution
# over phi
# *This is called as soon as the module is loaded*
def loadPhiInterp(path='.'):
    global phi_interp
    fname = os.path.join(path, "..", "data", "PhiIntegrals.dat")
    os.makedirs(os.path.join(path, "..", "data",), exist_ok=True)

    xvals = np.arange(-7, 7.001, 0.05)
    phivals = np.arange(0, np.pi + 0.1, 0.05)
    xlen = len(xvals)
    ylen = len(phivals)
    # print xlen, ylen
    if (os.path.exists(fname)):
        data = np.loadtxt(fname, usecols=(2,))
        z = data.reshape((xlen, ylen))

    else:
        print("    File '../data/PhiIntegrals.dat' doesn't exist...")
        print("    Calculating from scratch...")
        z = np.zeros((xlen, ylen))
        for i, x in enumerate(xvals):
            for j, phi in enumerate(phivals):
                z[i, j] = quad(lambda y: np.exp(x * np.cos(y)), 0, phi)[0]

        xgrid, phigrid = np.meshgrid(xvals, phivals, indexing='ij')
        np.savetxt(fname, list(zip(xgrid.flatten(), phigrid.flatten(), z.flatten())))

    phi_interp = interp2d(xvals, phivals, z.T)


def IntegralOverPhi(x, phi_max):
    if (phi_max < 0):
        return 0
    if (phi_max >= np.pi):
        return 2.0 * np.pi * scipy.special.i0(x)
    else:
        return 2.0 * phi_interp(x, phi_max)


IntegralOverPhiVec = np.vectorize(IntegralOverPhi)


# Integrand for integrating over the velocity distribution
# The phi integral has already been performed, so all you have
# left is v and theta.
def calcf_integ(v, theta, gamma):
    if (np.sin(gamma) * np.sin(theta) <= 1e-10):
        return 2.0 * np.pi * VelDist(v, theta, 0, gamma)

    ve = np.sqrt(2.0) * sigmav
    delsq = v ** 2 + ve ** 2 - 2 * v * ve * np.cos(gamma) * np.cos(theta)

    cosmin = (v ** 2 + ve ** 2 - vesc ** 2) / (2 * v * ve * np.sin(gamma) * np.sin(theta)) \
             - (np.cos(gamma) * np.cos(theta)) / (np.sin(gamma) * np.sin(theta))

    x0 = np.sin(theta) * np.sin(gamma) * v * ve / (sigmav ** 2)
    phi_max = np.arccos(np.clip(cosmin, -1.0, 1.0))
    A = IntegralOverPhiVec(x0, phi_max) * np.exp(-delsq / (2.0 * sigmav ** 2))

    return A * 1.0 / NNORM


# Full 3-D velocity distribution (v, theta, phi)
def VelDist(v, theta, phi, gamma):
    cdel = np.sin(gamma) * np.sin(theta) * np.cos(phi) + np.cos(gamma) * np.cos(theta)
    ve = np.sqrt(2.0) * sigmav
    dsq = v ** 2 - 2 * v * ve * cdel + ve ** 2
    A = np.exp(-dsq / (2 * sigmav ** 2)) / NNORM
    if hasattr(A, "__len__"):
        A[np.where(dsq > vesc ** 2)] = A[np.where(dsq > vesc ** 2)] * 0.0
    else:
        if (dsq > vesc ** 2):
            A = 0
    return A


# Calculate the free MB speed distribution
# after integrating over all angles
def calcf_SHM(v):
    ve = np.sqrt(2) * sigmav
    beta = ve / (sigmav ** 2)
    N1 = 1.0 / (Nesc * sigmav ** 3 * np.sqrt(2 * np.pi))
    f = v * 0.0

    a = (v <= vesc - ve)
    f[a] = np.exp(-(v[a] ** 2 + ve ** 2) / (2.0 * sigmav ** 2)) * (
                np.exp(beta * v[a]) - np.exp(-beta * v[a]))

    b = (vesc - ve < v) & (v < vesc + ve)
    f[b] = np.exp(-(v[b] ** 2 + ve ** 2) / (2.0 * sigmav ** 2)) * (
                np.exp(beta * v[b]) - np.exp((v[b] ** 2 + ve ** 2 - vesc ** 2) / (2 * sigmav ** 2)))
    return f ** v


# Minimum velocity required for a recoil of energy E_R
def vmin(E, m_N, m_x):
    res = E * 0.0
    m_N2 = m_N * 0.9315
    mu = (m_N2 * m_x) / (m_N2 + m_x)
    res = 3e5 * np.sqrt((E / 1e6) * (m_N2) / (2 * mu * mu))
    return res


vesc = None
sigmav = None
Nesc = None
NNORM = None
if __name__ == '__main__':
    loadPhiInterp()
    vesc = 533.0
    sigmav = 156.0

    # Nesc - normalisation constant
    Nesc = (scipy.special.erf(vesc / (np.sqrt(2.0) * sigmav)) - np.sqrt(2.0 / np.pi) * (
                vesc / sigmav) * np.exp(
        -vesc ** 2 / (2.0 * sigmav ** 2)))
    NNORM = Nesc * sigmav ** 3 * np.sqrt(2.0 * np.pi) * 2.0 * np.pi
