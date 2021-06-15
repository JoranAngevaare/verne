import os
import numpy as np
from scipy.interpolate import interp1d
import argparse
import pandas as pd
import sys
import scipy
import time
from verne import MaxwellBoltzmann as MB
import verne


def write_calcveldist(m_x, sigma_p, loc, N_gamma, v_esc, v_0, save_as):
    df_avg = avg_calcveldist(m_x, sigma_p, loc, N_gamma, v_esc, v_0)
    if save_as is not None and type(save_as) == str:
        print(f'saving at {save_as}')
        fname_avg = os.path.abspath(save_as)
    else:
        raise ValueError
    save_aververage_df(df_avg, fname_avg)


def avg_calcveldist(m_x, sigma_p, loc, N_gamma, v_esc, v_0):
    path = verne.__path__[0]

    MB.loadPhiInterp(path)

    if (loc == "SUF"):
        depth = 10.6  # metres
    elif (loc == "MPI"):
        depth = 0.3  # metres
    elif (loc == "EDE"):
        depth = 1.0  # metres
    elif (loc == "XENON"):
        depth = 1400  # metres
    elif (loc == "SNOLAB"):
        depth = 2070

    target = loc
    sigmav = v_0 / np.sqrt(2.0)  # 156.0

    # if N_gamma != 11 and args.save_as != None:
    #     args.save_as = 'tmp_' + args.save_as

    print("   ")
    print("    Calculating for...")
    print("        m_x/GeV:", m_x)
    print("        sigma_p/cm^2:", sigma_p)
    print("        detector at :", loc)
    print(f'        considering {N_gamma} angles')
    print(" ")

    # Initialise verne
    verne.core.loadIsotopes(path=path)
    verne.core.loadFFcorrections(m_x, path=path)

    # Calculate the maximum initial speed as a function of incoming angle theta
    Nvals = 1001
    thetavals = np.linspace(0, np.pi, Nvals)


    v_e = np.sqrt(2.0) * sigmav
    VMAX = 1000  # km/s
    VMIN = 1  # km/s

    # Nesc - normalisation constant
    Nesc = (scipy.special.erf(v_esc / (np.sqrt(2.0) * sigmav)) - np.sqrt(2.0 / np.pi) * (
                v_esc / sigmav) * np.exp(-v_esc ** 2 / (2.0 * sigmav ** 2)))
    NNORM = Nesc * sigmav ** 3 * np.sqrt(2.0 * np.pi) * 2.0 * np.pi

    MB.sigmav = sigmav
    MB.vesc = v_esc
    MB.NNORM = NNORM
    MB.Nesc = Nesc

    # Loop over gamma values
    Nv = 30
    gamma_list = np.linspace(0, 1.0, N_gamma)
    vgrid = np.zeros((N_gamma, Nv))
    fgrid = np.zeros((N_gamma, Nv))
    vgrid_average = np.zeros(Nv)
    fgrid_average = np.zeros(Nv)

    def getVelDist(gamma):
        print("        Calculating maximum final speed...")
        a = 1.0
        b = 2 * v_e * (-np.sin(gamma) * np.sin(np.pi - thetavals) + np.cos(gamma) * np.cos(
            np.pi - thetavals))
        c = v_e ** 2 - v_esc ** 2
        v_initial_max = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2.0 * a)

        # Calculate the maximum final speed as a function of incoming angle theta
        v_final_max = 0.0 * v_initial_max
        for i in range(Nvals):
            v_final_max[i] = verne.core.calcVfinal_full(
                v_initial_max[i],
                thetavals[i],
                depth,
                sigma_p,
                m_x,
                target)

        # Calculate interpolation function for max final speed
        vmax = np.max(v_final_max)
        if (vmax < 1.0):
            return np.linspace(0, 1, Nv), np.zeros(Nv)
        vfinal_interp = interp1d(thetavals, v_final_max, kind='linear', bounds_error=False,
                                 fill_value=0)

        print("        Calculating final speed distribution...")

        # Tabulate values of speed distribution
        # v_th = 1.0  # Lowest speed to consider (don't go lower than 1 km/s, other the calculation of derivatives is messed up...)

        # Generate a list of sampling values for v (with some very close to v_th)
        Nv_over_three = Nv // 3
        vlist = np.logspace(np.log10(VMIN), np.log10(0.25 * VMAX), Nv_over_three)  # 10
        vlist = np.append(vlist, np.linspace(0.15 * VMAX, 0.99 * VMAX, Nv - Nv_over_three))  # 20
        vlist = np.sort(vlist)
        f_final = 0.0 * vlist
        for i in range(len(vlist)):
            f_final[i] = verne.core.CalcF(vlist[i], gamma, depth, sigma_p, m_x, target, vfinal_interp)

        return vlist, f_final

    for j in range(N_gamma):
        print("    Calculating for gamma/pi = ", gamma_list[j], "...")
        res = getVelDist(gamma_list[j] * np.pi)
        vgrid[j, :], fgrid[j, :] = res
        vgrid_average += np.array(res[0])
        fgrid_average += np.array(res[1])

    vgrid_average /= N_gamma
    fgrid_average /= N_gamma

    # gamma_rep = np.repeat(gamma_list, Nv)
    # df = pd.DataFrame()
    # df['gamma/pi'] = gamma_rep
    # df['v_[km/s]'] = vgrid.flatten()
    # df['f(v,gamma)_[s/km]'] = fgrid.flatten()

    df_avg = pd.DataFrame()
    df_avg['v_[km/s]'] = vgrid_average
    df_avg['f(v,gamma)_[s/km]'] = fgrid_average

    return df_avg


def save_aververage_df(df, save_name):
    if os.path.exists(save_name):
        # Presumably another instance is also saving the same spectrum or has done so.
        print(f'WARNING:\t{save_name} already exists!')
    else:
        try:
            df.to_csv(save_name, index=False)
        except OSError as e:
            print(
                f'WARNING:\tOSError while writing {save_name} on {time.asctime()}. Taking a nap of 1 minutes and re-trying!')
            print(f'was trying to save df of len(df) = {len(df)}:\n{df}')
            print(f'The error:\n{e}')
            time.sleep(1 * 60)
            if os.path.exists(save_name):
                os.remove(save_name)
                df.to_csv(save_name, index=False)

    def check_on_save():
        df_read = pd.read_csv(save_name)
        return np.shape(df_read) == np.shape(df)

    # assert that the write is correct
    if check_on_save():
        print(f'Dataframe written successfully')
    else:
        print(f'ERROR in writing dataframe trying again in 5 minutes')
        print(f'was trying to save df of len(df) = {len(df)}:\n{df}')
        time.sleep(60)
        if os.path.exists(save_name):
            os.remove(save_name)
        df.to_csv(save_name, index=False)
        if not check_on_save():
            print(f'FATAL ERROR while writing {save_name}, force killing job now')
            sys.exit(-1)


def main():
    print('WELCOME TO THE UPDATED CALCVELDIST 20200417 (v2)')
    # Parse the arguments!
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('-m_x', '--m_x', help='DM mass in GeV', type=float, default=1e5)
    parser.add_argument('-v_0', '--v_0', help='v_0 in km/s', type=float, default=230)
    parser.add_argument('-v_esc', '--v_esc', help='v_esc in km/c', type=float, default=533)
    parser.add_argument('-n_gamma', '--n_gamma', help='number of angles considered', type=int,
                        default=4)
    parser.add_argument('-save_as', '--save_as', default=None,
                        help='name of csv file to save the averaged velocity-distribution',
                        type=str)
    parser.add_argument('-sigma_p', '--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2',
                        type=float, required=True)
    parser.add_argument('-loc', '--location', help='Detector location to consider. `MPI` or `SUF`',
                        type=str, required=True)

    args = parser.parse_args()
    m_x = args.m_x
    sigma_p = args.sigma_p
    loc = args.location
    N_gamma = args.n_gamma
    v_esc = args.v_esc  # 533.0
    v_0 = args.v_0
    save_as = args.save_as
    write_calcveldist(m_x, sigma_p, loc,v_esc, v_0, save_as, N_gamma=4)


if __name__ == '__main__':
    main()
