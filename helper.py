import matplotlib.pyplot as plt
from system_class import*
from particle_classes import*
from matplotlib.patches import Rectangle


mean_velocity = 461.3          #Mean velocity of oxygen at room temperature
mu = 1.66 * 10**(-27)          #Atomic weight
m_02 =2*16*mu                  #Weight of an oxygen atom
mol =  6.0221415 * 10**(23)

#The Gass is a Van der Waals gass
a = 1.382
b = 0.03186
R = 8.3144626181532

r_02 = 152*10**(-12)           #Radius of an oxygen atom
k = 1.38 * 10**(-23)           #Boltzmanns constant
T = 273 + 20                   #Kelvin at roomtemperature
P = 101325                     #atm


def evolution(system_1,times=1):
    x_list = []
    y_list = []
    vx_list = []
    vy_list = []
    for i in range(4):
        x = []
        y = []
        vx = []
        vy = []
        for p in system_1.particles:  
            x.append(p.position.x)
            y.append(p.position.y)
            vx.append(p.velocity.vx)
            vy.append(p.velocity.vy)

        x_list.append(x)
        y_list.append(y)
        vx_list.append(vx)
        vy_list.append(vy)

        system_1.update(times)
    return x_list,y_list,vx_list,vy_list
def boltsmann(v, m):
    return m*v/k / T * np.exp(- m* v**2 / 2 / k / T)
def Boltzmann_comparison(vel_start, vel_end,v, m = m_02):
    '''
    Input:
    vel_start: array with the initial velocity of all the particles.
    vel_end:   array with all the final velocities of the particles
    v: velocity distribution we are calculating the boltsmann distribution over

    Output: 
    Figure comparing the three distributions.
    '''

    fig =plt.figure(figsize = (15,5))

    counts_1, bins_1 = np.histogram(vel_start)
    counts_2, bins_2 = np.histogram(vel_end)

    plt.stairs(counts_1, bins_1, label = "vel_start" )
    plt.stairs(counts_2, bins_2, label = "vel_end" )
    plt.plot(v,boltsmann(v, m), label = "Boltzmann" )
    plt.title( "Comparing the velocity distributions")
    plt.legend()
    return fig


def plotting_evolution(x_list,y_list,vx_list,vy_list):
    fig, ((ax1, ax2),(ax3,ax4)) = plt.subplots(2, 2, figsize = (10,10))
    fig.suptitle('Particle evolution')

    r =1
    s =1
    ax1.scatter(x_list[0], y_list[0], s =r )
    ax1.add_patch(Rectangle((0,0), 1, 1, angle=0.0, fill=True, alpha = 0.3))
    ax1.quiver(x_list[0], y_list[0], vx_list[0], vy_list[0], scale=s)
    ax1.grid()

    ax2.scatter(x_list[1], y_list[1],s =r)
    ax2.add_patch(Rectangle((0,0), 1, 1, angle=0.0, fill=True, alpha = 0.3))
    ax2.quiver(x_list[1], y_list[1], vx_list[1], vy_list[1], scale=s)
    ax2.grid()

    ax3.scatter(x_list[2], y_list[2],s =r)
    ax3.add_patch(Rectangle((0,0), 1, 1, angle=0.0, fill=True, alpha = 0.3))
    ax3.quiver(x_list[2], y_list[2], vx_list[2], vy_list[2], scale=s)
    ax3.grid()

    ax4.scatter(x_list[3], y_list[3],s =r)
    ax4.add_patch(Rectangle((0,0), 1, 1, angle=0.0, fill=True, alpha = 0.3))
    ax4.quiver(x_list[3], y_list[3], vx_list[3], vy_list[3], scale=s)
    ax4.grid()

    return fig

def E(v,m):return m*v**2/2

def Calculating_the_energy(system, number_of_updates, d_collision = 0.001):
    '''
    Input: 
    system: The system we are calculated 
    number_of_updates: 
    d_collision: 

    Output: 


    '''
    vel_m0 = np.zeros((int(system.num_particles/2), number_of_updates))
    vel_4m0 = np.zeros((int(system.num_particles/2), number_of_updates))

    vel = system.return_velocities()
    vel_m0[:,0]=vel[:int(system.num_particles/2)]
    vel_4m0[:,0]=vel[int(system.num_particles/2):]

    average_current = system.average_particle_collision
    average_next = average_current 
    j = 1
    while system.average_particle_collision < 10 and j < number_of_updates:
        while average_next-average_current < d_collision: 
            system.update_step()
            average_next = system.average_particle_collision
        average_current = average_next

        vel = system.return_velocities()
        vel_m0[:,j]=vel[:int(system.num_particles/2)]
        vel_4m0[:,j]=vel[int(system.num_particles/2):]
        j +=1 
    return vel_m0, vel_4m0
def plot_kinetic_energy(vel_m0,vel_4m0):
    fig = plt.figure(figsize = (10,5))
    gs = fig.add_gridspec(1, 3, wspace=0)
    (ax1, ax2, ax3) = gs.subplots(sharey=True)
    ax1.plot(E(vel_m0,m_02).sum(axis = 0))
    ax1.set_title("Ek when m = m0")

    ax2.plot(E(vel_4m0,4*m_02).sum(axis = 0))
    ax2.set_title("Ek when m = 4m0")

    ax3.plot(E(vel_4m0,4*m_02).sum(axis = 0)+E(vel_m0,m_02).sum(axis = 0))
    ax3.set_title("Ek for all m")
    fig.suptitle("The average kinetic energy with elaticity 1")
    return fig