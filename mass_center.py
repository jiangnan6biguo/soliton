import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import sys

sys.path.append("/home/richard/文档/copy/code/")
import soliton_lib as soliton



index = "_q_25"
path0 = "/home/richard/文档/copy/code/1d_domain_wall/1DMDW_2/1DMDW/spin1_0"+index+".mat"
pathp = "/home/richard/文档/copy/code/1d_domain_wall/1DMDW_2/1DMDW/spin1_p"+index+".mat"



dt = 0.016

a = dt*np.arange(3000)

p0 = soliton.load_data(path0,"dmat_0")
pp = soliton.load_data(pathp,"dmat_p")

n0 = soliton.convert_to_density(p0)
np = soliton.convert_to_density(pp)

plt.subplot(121)
plt.pcolormesh(n0)
plt.subplot(122)
plt.pcolormesh(np)
plt.show()





"""
mass_center_mat=[]
for i in range(len(n0)):
    mass_center_mat.append(soliton.mass_center(n0[i]+2*np[i]))

plt.plot(a,mass_center_mat)
plt.show()

"""







