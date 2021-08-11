# Rotation velocity
def grad_mw(x,y,z):
	return -accel_mw(x,y,z)

def rot_vel_mw(r):
    return np.sqrt(r*grad_mw(r,0,0)[0])

exec(open("./vel_Sofue13.py").read())


#Plot
r=np.logspace(-3,3,200)

fig = plt.figure(figsize=(21,7))
font = {"size": 15}  
plt.rc('font', **font)
plt.scatter(r,rot_vel_mw(r),s=0.1,marker='o', color='red', label='RAR+barions')

plt.errorbar(v_Sof['r']/1.e3, v_Sof['v'], xerr=v_Sof['err_r']/1.e3 ,yerr=v_Sof['err_v'], fmt='o', color='cyan', label='Sofue+13')
#plt.ylim(-4,2)

plt.xlabel('r')
plt.ylabel('v_r')
plt.grid(True)
plt.xscale('symlog')
