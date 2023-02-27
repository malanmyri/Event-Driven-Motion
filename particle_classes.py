from heapq import*
import typing

class position: 
    def __init__(self, 
                 x: float, 
                 y:float,
                 ):
        self.x = x
        self.y = y
class velocity: 
    def __init__(self, 
                 vx: float, 
                 vy:float,
                 ):
        self.vx = vx
        self.vy = vy
class particle: 
    def __init__(self, 
                 position:position, 
                 velocity:velocity, 
                 radius:float,
                 mass:float,
    ):
        self.position = position
        self.velocity = velocity
        self.radius = radius
        self.mass = mass
        self.collision_count = 0
        self.tc = 0
class collision: 
    def __init__(self,
                 time: float,
                 entities:  typing.List[particle], 
                 collision_count:float,
                 type:float,
                 ):
        self.time = time 
        self.entities = entities 
        self.collision_count = collision_count
        self.type = type
    def __lt__(self, other): 
        if self.type == "wall_y" or self.type == "wall_x":
            return False
        else: return True
    def __le__(self, other): 
        if self.type == "wall_y" or self.type == "wall_x":
            return False
        else: return True
