import numpy as np
import matplotlib.pyplot as plt 
import netCDF4 as nc

file_nc = 'rigid_rotor_pid.nc'
HUB_ID = 1
BLADE_ID = 10
THETA_0 = 8.0
PID_ID = 1

ds = nc.Dataset(file_nc) 

time = ds.variables['time'][:] 
hub_z = ds.variables[f"node.struct.{HUB_ID}.X"][:,2] 
coll_rotation = ds.variables[f"elem.joint.{BLADE_ID}.Phi"][:,0] 
pid_gain_p = ds.variables[f"elem.loadable.{PID_ID}.P"][:] 


plt.figure()
plt.plot(time, hub_z)
plt.xlabel('Time [s]')
plt.ylabel('Hub Z Position [m]')
plt.title('Hub Z Position Over Time')
plt.grid()
plt.savefig('hub_z_position.png')   

plt.figure()
plt.plot(time, coll_rotation*180/np.pi + THETA_0)       
plt.xlabel('Time [s]')
plt.ylabel('Collective Rotation [deg]')
plt.title('Collective Rotation Over Time')
plt.grid()
plt.savefig('collective_rotation.png')

plt.figure()
plt.plot(time, pid_gain_p)
plt.xlabel('Time [s]')
plt.ylabel('PID Proportional Gain [N/m]')
plt.title('PID Proportional Gain Over Time')
plt.grid()
plt.savefig('pid_proportional_gain.png') 
