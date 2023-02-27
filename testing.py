from system_class import *
import matplotlib.pyplot as plt

num_particles =3
elasticity= 1
x_boundaries = 1
y_boundaries = 1
mean_velocity = 0.5
mass=[1]
radius=[.001]
tc = 0.001

system_1 = system(num_particles, elasticity,x_boundaries, y_boundaries,mean_velocity,mass,radius, tc)
system_1.uniform_particles()

for p in system_1.particles: 
    system_1.find_collisions_particle(p)

fig, ((ax1, ax2),(ax3,ax4)) = plt.subplots(2, 2)
fig.suptitle('Particle evolution')
from matplotlib.patches import Rectangle

r =10
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
    system_1.update()

ax1.scatter(x_list[0], y_list[0], s =r )
ax1.add_patch(Rectangle((0,0), 1, 1, angle=0.0, fill=True, alpha = 0.3))
ax1.quiver(x_list[0], y_list[0], vx_list[0], vy_list[0], scale=1)
ax1.grid()

ax2.scatter(x_list[1], y_list[1],s =r)
ax2.add_patch(Rectangle((0,0), 1, 1, angle=0.0, fill=True, alpha = 0.3))
ax2.quiver(x_list[1], y_list[1], vx_list[1], vy_list[1], scale=1)
ax2.grid()

ax3.scatter(x_list[2], y_list[2],s =r)
ax3.add_patch(Rectangle((0,0), 1, 1, angle=0.0, fill=True, alpha = 0.3))
ax3.quiver(x_list[2], y_list[2], vx_list[2], vy_list[2], scale=1)
ax3.grid()

ax4.scatter(x_list[3], y_list[3],s =r)
ax4.add_patch(Rectangle((0,0), 1, 1, angle=0.0, fill=True, alpha = 0.3))
ax4.quiver(x_list[3], y_list[3], vx_list[3], vy_list[3], scale=1)
ax4.grid()

plt.show()
print("hello_world")