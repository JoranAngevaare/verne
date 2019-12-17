import numpy as np
from scipy.interpolate import interp1d
import argparse
import pandas as pd
import sys
import scipy



#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float, default = 1e5)
parser.add_argument('-v_0','--v_0', help='v_0 in km/s', type=float, default=230)
parser.add_argument('-v_esc','--v_esc', help='v_esc in km/c', type=float, default=533)
parser.add_argument('-n_gamma','--n_gamma', help='number of angles considered', type=int, default=11)
parser.add_argument('-save_as', '--save_as', default = None,  help='name of csv file to save the averaged velocity-distribution', type=str)
parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
parser.add_argument('-loc','--location', help='Detector location to consider. `MPI` or `SUF`', type=str, required=True)

parser.add_argument('-path','-PATH', help='absolute path of verne', type=str, required=True)

args = parser.parse_args()
m_x = args.m_x
sigma_p = args.sigma_p
loc = args.location
path = args.path

sys.path.insert(1, args.path)
import verne
import MaxwellBoltzmann as MB
results_dir = args.path+"/../results/"

MB.loadPhiInterp(path)

if (loc == "SUF"):
    depth = 10.6 #metres
elif (loc == "MPI"):
    depth = 0.3 #metres
elif (loc == "EDE"):
    depth = 1.0 #metres
elif (loc == "XENON"):
    depth = 1400  # metres

target = loc
N_gamma = args.n_gamma
# if N_gamma != 11 and args.save_as != None:
#     args.save_as = 'tmp_' + args.save_as


print( "   ")
print( "    Calculating for...")
print( "        m_x/GeV:", m_x)
print( "        sigma_p/cm^2:", sigma_p)
# #print "        gamma/pi:", gamma_by_pi)
print( "        detector at :", loc)
print(f'        considering {N_gamma} angles')
print( " ")

#Initialise verne
verne.loadIsotopes(path = path)
verne.loadFFcorrections(m_x, path = path)

#Calculate the maximum initial speed as a function of incoming angle theta
Nvals = 1001
thetavals = np.linspace(0, np.pi, Nvals)

vesc = args.v_esc #533.0
sigmav= args.v_0/np.sqrt(2.0) #156.0
v_e = np.sqrt(2.0)*sigmav
# v_e = np.sqrt(2.0)*MB.sigmav
# vesc = MB.vesc



VMAX = 1000 #km/s
VMIN = 1 #km/s


# Nesc - normalisation constant
Nesc = (scipy.special.erf(vesc/(np.sqrt(2.0)*sigmav)) - np.sqrt(2.0/np.pi)*(vesc/sigmav)*np.exp(-vesc**2/(2.0*sigmav**2)))
NNORM = Nesc*sigmav**3*np.sqrt(2.0*np.pi)*2.0*np.pi


MB.sigmav = sigmav
MB.vesc = vesc
MB.NNORM = NNORM
MB.Nesc = Nesc


def getVelDist(gamma):
    
    print( "        Calculating maximum final speed...")
    a = 1.0
    b = 2*v_e*(-np.sin(gamma)*np.sin(np.pi-thetavals) + np.cos(gamma)*np.cos(np.pi-thetavals))
    c = v_e**2 - vesc**2
    v_initial_max = (-b + np.sqrt(b**2 - 4*a*c))/(2.0*a)

    #Calculate the maximum final speed as a function of incoming angle theta
    v_final_max = 0.0*v_initial_max
    for i in range(Nvals):
        v_final_max[i] = verne.calcVfinal_full(v_initial_max[i], thetavals[i],  depth, sigma_p, m_x, target)

    #Calculate interpolation function for max final speed
    vmax = np.max(v_final_max)
    if (vmax < 1.0):
        return np.linspace(0, 1, 61), np.zeros(61)
    vfinal_interp = interp1d(thetavals, v_final_max, kind='linear', bounds_error=False, fill_value=0)
    
    print( "        Calculating final speed distribution...")

    #Tabulate values of speed distribution
    v_th = 1.0 #Lowest speed to consider (don't go lower than 1 km/s, other the calculation of derivatives is messed up...)

    #Generate a list of sampling values for v (with some very close to v_th)
    # vlist = np.logspace(np.log10(v_th), np.log10(0.25*vmax), 20)    #20
    # vlist = np.append(vlist, np.linspace(0.15*vmax, 0.99*vmax, 40)) #40
    # TODO, I want the same spacing for all angles
    vlist = np.logspace(np.log10(VMIN), np.log10(0.25 * VMAX), 20)  # 20
    vlist = np.append(vlist, np.linspace(0.15*VMAX, 0.99*VMAX, 40)) #40
    vlist = np.sort(vlist)
    f_final = 0.0*vlist
    for i in range(len(vlist)):
        f_final[i] = verne.CalcF(vlist[i], gamma, depth, sigma_p, m_x, target, vfinal_interp)

    #Add on the final point
    # vlist = np.append(vlist, vmax)
    # f_final = np.append(f_final, 0.0)

    return vlist, f_final
    
    
#Loop over gamma values
Nv = 60
gamma_list = np.linspace(0, 1.0, N_gamma)
vgrid = np.zeros((N_gamma, Nv))
fgrid = np.zeros((N_gamma, Nv))
vgrid_average = np.zeros(Nv)
fgrid_average = np.zeros(Nv)
for j in range(N_gamma):
    print( "    Calculating for gamma/pi = ", gamma_list[j],"...")
    res = getVelDist(gamma_list[j]*np.pi)
    vgrid[j,:], fgrid[j,:] = res
    vgrid_average += np.array(res[0])
    fgrid_average += np.array(res[1])

vgrid_average /= N_gamma
fgrid_average /= N_gamma

gamma_rep = np.repeat(gamma_list, Nv)

#Output to file
# fname = results_dir + "veldists/f_" + loc + "_lmx" + '{0:.1f}'.format(np.log10(m_x)) + "_lsig" + '{0:.2f}'.format(np.log10(sigma_p)) + ".txt"
fname = results_dir + "veldists/f_all_%s_%i_%i_%.1f_%.2f.txt"%(loc, args.v_0, args.v_esc, np.log(args.sigma_p), args.m_x)
# headertxt = "mx [GeV]: " + str(m_x) + "\nsigma_p [cm^2]: " + str(sigma_p) + "\ndepth [m]: " + str(depth) + "\nloc: " + target
# headertxt += "\nColumns: gamma/pi, v [km/s], f(v, gamma) [s/km]"
#
# np.savetxt(fname, np.transpose([gamma_rep, vgrid.flatten(), fgrid.flatten()]), header=headertxt)
df = pd.DataFrame()
df['gamma/pi'] = gamma_rep
df['v_[km/s]'] = vgrid.flatten()
df['f(v,gamma)_[s/km]'] = fgrid.flatten()
try:
    df.to_csv(fname[:-3] + 'csv', index=False)
except:
    pass

if args.save_as is not None and type(args.save_as) == str:
    print(f'saving at {args.save_as}')
    fname_avg = args.save_as
else:
    fname_avg = results_dir + "veldists/f_" + loc + "_lmx" + '{0:.1f}'.format(np.log10(m_x)) + "_lsig" + '{0:.2f}'.format(np.log10(sigma_p)) + "_avg" + ".csv"
df_avg = pd.DataFrame()
df_avg['v_[km/s]'] = vgrid_average
df_avg['f(v,gamma)_[s/km]'] = fgrid_average
try:
    df_avg.to_csv(fname_avg, index=False)
except:
    pass