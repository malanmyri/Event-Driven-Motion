from system_class import*
from particle_classes import*
from matplotlib.patches import Rectangle
from scipy import optimize
import matplotlib.pyplot as plt

k_b = 1.38 * 10**(-23)       #Boltzmanns constant
R = 8.3144626181532



#The Gass is a Van der Waals gass
b = 0.03186                    #This is the volume one mole of oxygen atoms take up, https://en.wikipedia.org/wiki/Van_der_Waals_constants_(data_page)

mean_velocity = 461.3          #Mean velocity of oxygen at zero Kelvin, [m/s] https://www.thoughtco.com/root-squmean-velocity-example-problem-607556
mu = 1.66 * 10**(-27)          #Atomic weight [kg]
m_O2 =2*16*mu                  #Weight of O2 [kg]


#Calculating the temperature used in the boltzmann distribution so that it follows the equipartition theorem in two dimensions.

T = mean_velocity**2 * m_O2/2/k_b   #The temperatur set for this simulation [Kelvin]



P = 101325                     #The pressure set for this experiment [pascal]
r_02 = 152*10**(-12)           #Radius of an oxygen atom [meters]
mol =  6.0221415 * 10**(23)    #1 mol



def Van(n): 
    return  n*b + n*R*T/P 

def boltsmann(v): 
    return m_O2*v/k_b /T * np.exp(-m_O2*v**2/2/k_b/T)



# Initial Conditions

num_particles =5000
n = num_particles/mol      #Number of moles
V = Van(n)                        #Volume of the box as calculated from Van der Waals
x_b = y_b = V**(1/3)           #Box boundaries given in meters
elasticity = 1
tc = 0
mass = [m_O2]
radius = [r_02]

system_1 = system(num_particles, elasticity,x_b, y_b,mean_velocity,mass,radius, tc)     # Initializing the system
system_1.uniform_particles()                                                            # Putting the particles in the box with a uniform distribution
system_1.find_collisions()          

vel_start = system_1.return_velocities()
for time in range(100):
    times = 1000
    print(time)
    system_1.update(times)
    system_1.convergence()
vel_end = system_1.return_velocities()

print("Done updating")


fig, ax = plt.subplots()
ax.plot(np.linspace(0,100*1000, 100), system_1.list_of_differences)
fig.suptitle("Convergence")
ax.set_xlabel(f"number of updates")
ax.set_ylabel("var [m²/s²]")
ax.set_title("Variance in the velocity distribution")
plt.savefig("convergence.pdf")



fig, (ax1,ax2) = plt.subplots(1,2, figsize = (15,5), sharex = True)

ax1.hist(vel_start[0],1, label = "vel_start" , density = True, histtype = 'step')
ax1.legend()
ax1.set_xlabel("Velocity [m/s]")
ax1.set_ylabel("Probability")
ax1.set_title("Velocity at the start")

v = np.linspace(0,mean_velocity*3,100)
N = sum(boltsmann(v)*(v[1]-v[0]))                     #Normalisation constant for the Boltzmann distribution

ax2.hist(vel_end,1000, label = "vel_end" , density = True, histtype = 'step')  #1000 is the number of bins we want to plot for 
ax2.plot(v,boltsmann(v)/N, label = "Boltzmann" )
ax2.set_xlabel("Velocity [m/s]")
ax2.set_ylabel("Probability")
ax1.set_title("Velocity at the end")
ax2.legend()


fig.suptitle( "Comparing the velocity distributions")
plt.savefig("Velocity distribution.pdf")





