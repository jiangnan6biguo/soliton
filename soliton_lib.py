import matplotlib.pyplot as plt
import scipy.io as scio
import numpy as np
from scipy.fftpack import *
from matplotlib.animation import FuncAnimation

"""
This module is designed to work for soliton data process.
The system is scaler now. May we will develop tools for more complex situation in the future.
Data should match the structure used in XiaoQuan's code.
This module may help user to do some calculation, load data and draw some data.

Original Version
2021.06.01
Richard
"""

x_max = 160.0
x_len = 4096
dt = 0.1
dx = 2*x_max/x_len
suffix = "_norm" #This variable is useful when we want to save multi-file. user can set suffix by him/herself
omega = 0.06




def load_data(data_file,name="dmat"):
    data = scio.loadmat(data_file)[name]
    print("The shape of data is: ",data.shape)

    return data


def load_ground(data_file):
    ground_data = scio.loadmat(data_file)["psip"]

    return ground_data




def convert_to_density(data):
    """
    The data may be the Complex wave function, using convert_to_density will return the density distribution.
    """
    dim = len(data.shape)
    if dim == 1:
        density = []
        for i in range(len(data)):
            density.append((np.abs(data[i])**2))
        density = np.array(density)

    elif dim == 2:
        density = np.zeros(data.shape)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                density[i,j] = np.abs(data[i,j])**2

    else :
        print("The structure of input mat has some probelm, please check! nparray type data is expected, and the dimension should not exceed 2.")


    return density






def find_the_vortex(u,pre_index,bound=0):
    vortex_position = pre_index
    while u[vortex_position+1]<u[vortex_position] or u[vortex_position-1]<u[vortex_position]:
        if u[vortex_position+1]< u[vortex_position]:
            vortex_position = vortex_position+1
            
        elif u[vortex_position-1]<u[vortex_position]:
            vortex_position = vortex_position-1

        if vortex_position >= len(u)-1-bound or vortex_position <= bound:
            vortex_position = 1+bound
            break

    return vortex_position


def mass_center(u):
    weight = 0.0
    position = 0.0
    for i in range(len(u)):
        position = position + u[i]*i
        weight = weight + u[i]

    return position/weight




def soliton_move(data,ground_data,interval=5,pre_position=512):
    """
    This function will find the velocity of soliton moving in scaler BEC.
    The return datas are:
    1. position 2.depth of soliton
    3.velocity of soliton moving
    
    The input data should be density of wave function. Not Complex data.
    The time interval may be too small, so we ignore some media data and the interval is the time index difference.
    """


    position_set = []
    depth_set = []
    for i in range(len(data)):
        pos = find_the_vortex(data[i],pre_position)
        position_set.append(pos)
        depth_set.append(ground_data[pos]-data[i,pos])
        pre_position = pos


    velocity = []
    for i in range(len(position_set)//interval - 1):
        velocity.append((position_set[(i+1)*interval]-position_set[i*interval])/(dt*interval))

    return np.array(position_set),np.array(depth_set),np.array(velocity)*dx






def soliton_energy(data,ground,average=True):
    

    """
    This function can be used to calculate the Energy of soliton. Both the ground data and the exicited data should be given.
    
    The return datas are:
    1. The total energy of the soliton
    2. The kinetic energy of the soliton
    3. The Potential energy of the soliton
    4. The interaction energy of the soliton

    If the average is True, the return will be the average. If false, will be time-dependent.
    

    Noticing! The data shounld be wavefunction, not density.


    """

    ground_kinetic = 0.0

    for i in range(len(ground)-1):
        ground_kinetic = ground_kinetic + np.abs(ground[i+1]-ground[i])**2/dx**2

    ground_interaction = 0.0
    for i in range(len(ground)):
        ground_interaction = ground_interaction + 0.5*np.abs(ground[i])**4





    data_visual = np.zeros((len(data),len(data[0])-1))
    data_kinetic = np.zeros(len(data))
    kine_ground = np.zeros(len(data))+ground_kinetic



    ground_num = 0.0
    for i in range(len(ground)):
        ground_num = ground_num + np.abs(ground[i])**2




    pot_func = np.zeros(len(data[0]))
    for i in range(x_len):
        pot_func[i] = ((float(i)*dx-x_max)**2)*(omega**2)/2

    ground_potential = 0.0

    for i in range(len(ground)):
        ground_potential = ground_potential + pot_func[i]*np.abs(ground[i])**2





    interaction_energy = np.zeros(len(data))
    total_energy = np.zeros(len(data))
    Pot_energy = np.zeros(len(data))
    norm_check = np.zeros(len(data))
    






    for i in range(len(data)):
        potential = 0.0
        inter_energy = 0.0
        num = 0.0
        for  j in range(len(data[0])):
            potential = potential + pot_func[j]*(np.abs(data[i][j])**2)
            inter_energy = inter_energy + abs(data[i,j])**4/2.0

            num = num+abs(data[i,j])**2
        for j in range(len(data[0])-1):
            data_visual[i,j] = (abs(data[i,j+1]-data[i,j])**2)/(dx**2)/2.0

        data_kinetic[i] = sum(data_visual[i])
        interaction_energy[i] = inter_energy
        Pot_energy[i] = potential
        norm_check[i] = num

        total_energy[i] = data_kinetic[i]+interaction_energy[i]+Pot_energy[i]

    factor = norm_check[0]/ground_num
    total_ground = factor*ground_potential + factor*kine_ground+(factor**2)*ground_interaction
    print(total_ground)

    if average == True:
        total_energy = sum(total_energy)/len(total_energy)
        kinetic = sum(data_kinetic)/len(data_kinetic)
        Pot_energy = sum(Pot_energy)/len(Pot_energy)
        inter_energy = sum(interaction_energy)/len(interaction_energy)
        return total_energy,kinetic,Pot_energy,inter_energy

    else :
        return total_energy,data_kinetic,Pot_energy,interaction_energy



def iner_pro(data,fit):
    N = min([len(data),len(fit)])
    res = 0.0
    for i in np.arange(N):
        res = res+data[i]*fit[i]

    return res/N



def cos_fit(dev,ave,omega,t_range):
    return -dev*np.cos(t_range*omega)

def cos_fit_p(dev,ave,omega,t_range):
    return dev*np.cos(t_range*omega)



def circular_fit(data,t_range,omega_mat):
    ave = sum(data)/len(data)
    dev = (max(data)-min(data))/2

    fit_mat = []

    for i in omega_mat:
        fit_mat.append(iner_pro(np.array(data)-ave,cos_fit(dev,ave,i,t_range)))

    return fit_mat
            
 
def circular_fit_p(data,t_range,omega_mat):
    ave = sum(data)/len(data)
    dev = (max(data)-min(data))/2

    fit_mat = []

    for i in omega_mat:
        fit_mat.append(iner_pro(np.array(data)-ave,cos_fit_p(dev,ave,i,t_range)))

    return fit_mat
    




















