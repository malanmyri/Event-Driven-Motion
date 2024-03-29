import numpy as np
from heapq import*
import typing
from particle_classes import *
from random import sample

class system: 
    def __init__(self, 
                 num_particles: int, 
                 elasticity: float,
                 x_lim: int, 
                 y_lim: int, 
                 mean_velocity: float,
                 mass: typing.List[float],
                 radius: typing.List[float],
                 tc: float,
                ):
        self.num_particles = num_particles                    #Number of particles in the system
        self.elasticity = elasticity                          #Elasticity constant of the system 
        self.particles = []                                   #List that will hold all the particles in the system
        self.collisions = []                                  #A min heap with all collisions, sorted by the time of incident. 
        self.x_lim = x_lim                                    #x boundaries of the box
        self.y_lim = y_lim                                    #y boundaries of the box
        self.mean_velocity = mean_velocity                    #This is actually not the mean velocity, but the velocity used in the calculations
        self.mass = mass                                      #List containing the masses of the particles. 
        self.radius = radius                                  #List containing the different radiuses of the particles.
        self.tc = tc
        self.current_time = 0
        self.invalid_time =[]
        self.particle_collisions = 0
        self.average_particle_collision = 0


        #Used in the calculation of convergence
        self.vel_curr = []
        self.vel_prior =[self.mean_velocity]
        self.list_of_differences = []
        self.collision_length = []
        self.number_of_vel_mean =[]
        self.tol = 0.1

    
    def uniform_particles(self):
        positions = set([])
        while len(positions) <self.num_particles:
            p_x= np.random.randint(1, int((self.x_lim- self.radius[0])/self.radius[0]), 1)[0]*(self.radius[0])
            p_y= np.random.randint(1, int((self.y_lim- self.radius[0])/self.radius[0]),  1)[0]*(self.radius[0])
            
            positions.add(position(p_x,p_y))

        theta = np.random.uniform(0, 2*np.pi, size=self.num_particles)
        positions = list(positions)
        if len(self.mass) ==1:
            for i in range(self.num_particles):
                
                vel = velocity(np.cos(theta[i])*self.mean_velocity, np.sin(theta[i])*self.mean_velocity)
                self.particles.append(particle(positions[i], vel, self.radius[0], self.mass[0]))
        if len(self.mass) == 2: 
            for i in range(int(self.num_particles/2)):
                vel = velocity(np.cos(theta[i])*self.mean_velocity, np.sin(theta[i])*self.mean_velocity)
                self.particles.append(particle(positions[i], vel, self.radius[0], self.mass[0]))
            for i in range(int(self.num_particles/2), self.num_particles):
                vel = velocity(np.cos(theta[i])*self.mean_velocity, np.sin(theta[i])*self.mean_velocity)
                self.particles.append(particle(positions[i], vel, self.radius[0], self.mass[1]))
    
    def convergence(self):

        if len(self.mass) ==1:
            count = 0
            for p in self.particles:
                if  abs(np.sqrt(p.velocity.vx**2 + p.velocity.vy**2) - self.mean_velocity) < self.tol:

                    count+=1
            
            self.number_of_vel_mean.append(count)

        '''
        if len(self.mass) ==1:
            i = 0
            for p in self.particles:
                if  np.sqrt(p.velocity.vx**2 + p.velocity.vy**2) != self.mean_velocity:
                    self.vel_curr.append(np.sqrt(p.velocity.vx**2 + p.velocity.vy**2))
            
            sqr_prior = sum((np.ones(len(self.vel_prior))*self.mean_velocity- self.vel_prior)**2 /len(self.vel_prior))
            sqr_curr =  sum((np.ones(len(self.vel_curr))*self.mean_velocity- self.vel_curr)**2 /len(self.vel_curr))

            self.vel_prior = self.vel_curr
            self.vel_curr = []
            print(abs(sqr_prior - sqr_curr))
            self.list_of_differences.append(abs(sqr_prior - sqr_curr))
        
        if len(self.mass) == 2:
            i = 0
            for p in self.particles:
                self.vel_curr[i] = np.sqrt(p.velocity.vx**2 + p.velocity.vy**2)
                i+=1
            
            sqr_prior_m_O2 = sum((np.ones(int(self.num_particles/2))*self.mean_velocity- self.vel_prior[:int(self.num_particles/2)])**2 /self.num_particles*2)
            sqr_prior_4m_O2 = sum((np.ones(int(self.num_particles/2))*self.mean_velocity- self.vel_prior[int(self.num_particles/2):])**2 /self.num_particles*2)


            sqr_curr_m_O2 = sum((np.ones(int(self.num_particles/2))*self.mean_velocity- self.vel_curr[:int(self.num_particles/2)])**2 /self.num_particles*2)
            sqr_curr_4m_O2 = sum((np.ones(int(self.num_particles/2))*self.mean_velocity- self.vel_curr[int(self.num_particles/2):])**2 /self.num_particles*2)


            self.vel_prior = self.vel_curr.copy()
            self.list_of_differences.append([abs(sqr_prior_m_O2 - sqr_curr_m_O2),abs(sqr_prior_4m_O2 - sqr_curr_4m_O2)])
        '''



    def calculating_the_energy(self, number_of_updates,m_O2 ,average_max = 0.001, stopping_condition = 10):
        '''
        Input: 
        system: The system we are calculated 
        number_of_updates: The maximum number of updates before we stop the calculation 
        average_max: Defines the step size over which we calculate the averages.
        stopping_condition: the simulation will run untill the average number of particle-particle collisions per article reaches this number.
        m_O2: mass of one oxygen molecule
        Output: 

        average kinetic energy over all particles
        average kinetic energy over all particles weighing m_O2
        average kinetic energy over all particles weighing 4*m_O2
        '''

        if len(self.mass) != 2: 
            print("This function can only be used for systems with particles of different masses.")
            return None

        vel_m0 = np.zeros((int(self.num_particles/2), number_of_updates))
        vel_4m0 = np.zeros((int(self.num_particles/2), number_of_updates))

        vel = self.return_velocities()
        vel_m0[:,0]=vel[:int(self.num_particles/2)]
        vel_4m0[:,0]=vel[int(self.num_particles/2):]

        average_current = self.average_particle_collision
        average_next = average_current 
        j = 1
        while self.average_particle_collision < stopping_condition and j < number_of_updates:
            while average_next-average_current < average_max: 
                self.update_step()
                average_next = self.average_particle_collision
            
            average_current = average_next
            vel = self.return_velocities()
            vel_m0[:,j]=vel[:int(self.num_particles/2)]
            vel_4m0[:,j]=vel[int(self.num_particles/2):]
            j +=1 

        average_m_O2 = (m_O2*vel_m0**2/2).sum(axis = 0)
        average_4m_O2= (4*m_O2*vel_4m0**2/2).sum(axis = 0)
        average_all = average_m_O2 + average_4m_O2
        return average_m_O2 , average_4m_O2 , average_all

    def crater(self):
        for num in range(self.num_particles):
            p = position(np.random.random_integers(self.radius[0], int((self.x_lim- self.radius[0])/self.radius[0]/2), 1)[0]*2*self.radius[0],np.random.random_integers(self.radius[0] ,int((self.y_lim- self.radius[0]-self.y_lim/2)/self.radius[0]/2), 1)[0]*self.radius[0]*2)
            theta = np.random.uniform(0, 2*np.pi, 1)
            vel = velocity(np.cos(theta[0])*self.mean_velocity, np.sin(theta[0])*self.mean_velocity)
            self.particles.append(particle(p, vel, self.radius[0], self.mass[0]))

    def find_collisions_particle(self,p):

        #Wall collisions: 
        delta_tx = np.inf
        if p.velocity.vx>0:
            delta_tx = (self.x_lim-p.radius - p.position.x)/p.velocity.vx
        if p.velocity.vx<0:
            delta_tx = (p.radius - p.position.x)/p.velocity.vx

        delta_ty = np.inf
        if p.velocity.vy>0:
            delta_ty = (self.y_lim-p.radius - p.position.y)/p.velocity.vy
        if p.velocity.vy<0:
            delta_ty = (p.radius - p.position.y)/p.velocity.vy

        if delta_tx<delta_ty: 
            
            heappush(self.collisions, (delta_tx + self.current_time, collision(delta_tx+ self.current_time, [p],p.collision_count, "wall_x")))
        if delta_ty<delta_tx:
            heappush(self.collisions, (delta_ty+ self.current_time, collision(delta_ty+ self.current_time, [p],p.collision_count, "wall_y")))


        #Particle collisions - This for loop can be more efficient 
        for particle_j in self.particles: 
            delta_x = np.array([particle_j.position.x - p.position.x, particle_j.position.y - p.position.y])
            delta_v = np.array([particle_j.velocity.vx - p.velocity.vx, particle_j.velocity.vy - p.velocity.vy])
            d = delta_v.T.dot(delta_x)**2 - delta_v.T.dot(delta_v)*(delta_x.T.dot(delta_x) - (p.radius + particle_j.radius)**2)
            if d > 0 and delta_v.T.dot(delta_x)< 0:
                delta_t = -(delta_v.T.dot(delta_x) + np.sqrt(d))/delta_v.T.dot(delta_v)
                if delta_t > 0: 
                    c = collision(delta_t+ self.current_time, [p, particle_j],p.collision_count + particle_j.collision_count, "particle")
                    heappush(self.collisions, (delta_t+ self.current_time, c))
                else:
                    self.invalid_time.append(self.current_time)
    def find_collisions(self):
        for p in self.particles: 
            self.find_collisions_particle(p)   
    def find_collision_projectile(self):
        self.find_collisions_particle(self.particles[-1]) #Because the last particle is the projectile


    def update_step(self):
        no_valid_collision = True
        particles = []
        while no_valid_collision:
            col_start = list(heappop(self.collisions))[1]
            col = col_start
            while col.time == col_start.time:
                dt = col.time- self.current_time
                if col.type == "wall_x" and col.collision_count == col.entities[0].collision_count:
                    if  col.time - col.entities[0].tc > self.tc:
                        #moving the particle to the wall
                        col.entities[0].position.x +=  col.entities[0].velocity.vx*dt
                        col.entities[0].position.y +=  col.entities[0].velocity.vy*dt

                        #changing the velocities 
                        col.entities[0].velocity.vx = - self.elasticity *col.entities[0].velocity.vx
                        col.entities[0].velocity.vy = self.elasticity *col.entities[0].velocity.vy
                        col.entities[0].tc = col.time 

                    else: #Here we use elasticity equal to 1
                        col.entities[0].position.x +=  col.entities[0].velocity.vx*dt
                        col.entities[0].position.y +=  col.entities[0].velocity.vy*dt
                        col.entities[0].velocity.vx = - col.entities[0].velocity.vx
                        col.entities[0].velocity.vy = col.entities[0].velocity.vy
                        
                    col.entities[0].collision_count += 1
                    no_valid_collision = False
                    particles.append(col.entities[0])
                        
                if col.type == "wall_y"and col.collision_count == col.entities[0].collision_count:

                    if col.time - col.entities[0].tc > self.tc:
                        col.entities[0].position.x +=  col.entities[0].velocity.vx*dt
                        col.entities[0].position.y +=  col.entities[0].velocity.vy*dt
                        col.entities[0].velocity.vx = self.elasticity *col.entities[0].velocity.vx
                        col.entities[0].velocity.vy = -self.elasticity *col.entities[0].velocity.vy
                        col.entities[0].tc = col.time 
                    else:
                        col.entities[0].position.x +=  col.entities[0].velocity.vx*dt
                        col.entities[0].position.y +=  col.entities[0].velocity.vy*dt
                        col.entities[0].velocity.vx = col.entities[0].velocity.vx
                        col.entities[0].velocity.vy = -col.entities[0].velocity.vy 
                        
                    col.entities[0].collision_count += 1
                    no_valid_collision = False
                    particles.append(col.entities[0])

                
                    
                if col.type == "particle" and col.collision_count == (col.entities[0].collision_count + col.entities[1].collision_count):

                    if  col.time - col.entities[0].tc < self.tc or  col.time - col.entities[1].tc < self.tc:
                        self.particle_collision(dt, col.entities[0], col.entities[1], 1)
                    else: 
                        self.particle_collision(dt, col.entities[0], col.entities[1], self.elasticity)
                        col.entities[0].tc = col.time 
                        col.entities[1].tc = col.time 
                    particles.append(col.entities[0])
                    particles.append(col.entities[1])

                    self.particle_collisions += 2
                    no_valid_collision = False

                    #For the odd chance if the heap is empty
                if len(self.collisions) >0:
                    col = list(heappop(self.collisions))[1]
                else:
                    break
            heappush(self.collisions,(col.time, col))

        #Decreasing the delta time of all the future collisions
        col = col_start  
        self.current_time = col.time
             
        
        #Updating positions of all the other particles: 
        for p in self.particles: 
            if p not in particles:
                p.position.x +=  p.velocity.vx*dt
                p.position.y +=  p.velocity.vy*dt
        
        #Finding all the future potential collisions with the new particles. 
        for p in particles: 
            self.find_collisions_particle(p)

        self.average_particle_collision = self.particle_collisions/ self.num_particles


    def particle_collision(self,dt, pj, pi, elasticity):
        delta_v =  np.array([pj.velocity.vx - pi.velocity.vx, pj.velocity.vy - pi.velocity.vy])
        if dt <0: print(f"negativ dt er lik {dt}")
        pi.position.x +=  pi.velocity.vx*dt
        pi.position.y +=  pi.velocity.vy*dt

        pj.position.x +=  pj.velocity.vx*dt
        pj.position.y +=  pj.velocity.vy*dt


        mass_t = pi.mass + pj.mass
        radius_t = pi.radius + pj.radius
        
        delta_x = np.array([pj.position.x - pi.position.x, pj.position.y - pi.position.y])

        vel_i = np.array([pi.velocity.vx,pi.velocity.vy]) + ((1+elasticity)*(pj.mass/mass_t)*(delta_v.dot(delta_x)/(radius_t**2)))*delta_x
        pi.velocity.vx = vel_i[0]
        pi.velocity.vy = vel_i[1]

        vel_j = np.array([pj.velocity.vx,pj.velocity.vy]) - ((1+elasticity)*(pi.mass/mass_t)*(delta_v.dot(delta_x)/(radius_t**2)))*delta_x
        pj.velocity.vx = vel_j[0]
        pj.velocity.vy = vel_j[1]

        #Update Collision count 
        pj.collision_count +=1 
        pi.collision_count +=1


    def return_velocities(self):
        vel = np.zeros(len(self.particles))
        i = 0
        for p in self.particles:
            vel[i] = np.sqrt(p.velocity.vx**2 + p.velocity.vy**2)
            i+=1
        return vel
    
    def energy(self):
        vel = self.return_velocities()
        energy_crater = (vel[:-1]**2).sum()*self.mass[0]/2    #This is the energy of all the small particles
        projectile = self.particles[-1]
        energy_projectile = vel[-1]**2 * projectile.mass /2   #energy of the projectile
        return energy_crater + energy_projectile
    

    def update(self,times):
        for time in range(times): 
            self.update_step()
            self.collision_length.append(len(self.collisions))
            

