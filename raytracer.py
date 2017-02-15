import pdb
import collections
from math import sqrt

max_depth = 10
vec = collections.namedtuple('vec', 'x y z')
#sphere_center = vec(0.0, 0.0, 1.0)

spheres = []
spheres.append(vec(0.0, 0.45, 1.0))
spheres.append(vec(0.0, -0.45, 1.0))
#spheres.append(vec(0.0, 0.3, 0.8))
#spheres.append(vec(-0.2, -0.3, 1.0))

sphere_colors = {}
sphere_colors[spheres[0]] = (255, 0, 0)
sphere_colors[spheres[1]] = (245, 240, 82)
spec_value = 5.0
spec_power = 100

sphere_radius = 0.5
light_pos = vec(0.5, 0.5, 0.0)
reso = 255

"""
Lambert is defined as dot product between surface normal & vector towards light source
"""
def calculate_lambert(intersect, sphere_center):
    #pdb.set_trace()
    #subtract initial point from terminal point
    #terminal_point - initial_point
    
    #light_pos - intersect = vector from intersection to light start
    light_dir = normalize(sub(light_pos, intersect))
    
    #intersect - sphere_center = vector from sphere center to intersection --> sphere normal
    sphere_normal = normalize(sub(intersect, sphere_center))
   
    lambert = max(0.0, dot(sphere_normal, light_dir))
    if lambert == 0.0:
        return (0, 0, 0)
    else:
        color = sphere_colors[sphere_center]
        color = map(lambda x: x * lambert, color)
        return color


def vector_len(v):
    return sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2)


def normalize(v):
    vlen = vector_len(v)
    return vec(v.x / vlen, v.y / vlen, v.z / vlen)


def pixel_coordinate_to_world_coordinate(coordinate):
    # copy pasted from somewhere, 32 is window height/width and rest are magic
    # plane from -0.5 to 0.5, multiply by 2 makes it -1 .. 1 ?
    return ((coordinate / float(reso)) - 0.5) * 2.0


def dot(v1, v2):
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z


def sub(v1, v2):
    return vec(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z)


def plus(v1, v2):
    return vec(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z)


#in reality this is line intersection instead of vector, thats why we only take positive t!
def ray_intersects_sphere(ray_origin, ray_direction, spheres): 
    """https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    simplified 2nd degree polynomial
    
    D discrimant (term under square root)
    dist -> distance from sphere center to ray origin 
    
    grabs the closest intersection
    """
    point = None
    sphere = None
    t_min = 0.0

    for sphere_center in spheres:
        dist = sub(ray_origin, sphere_center)
        b = dot(ray_direction, dist) 
        c = dot(dist, dist) - sphere_radius ** 2
        D = b * b - c
        if D >= 0:
            t0 = -b - sqrt(D)
            t1 = -b + sqrt(D)
            t = min(t0, t1)
            #negative t goes away from light position -> param d in line equation, distance along line from starting point
            if t > 0.0 and (not point or t_min > t):
                point = sphere_point(ray_origin, ray_direction, t) 
                sphere = sphere_center
                t_min = t
    
    return point, sphere 


def is_shadowed(ray_origin, intersect_sphere):
    other_spheres = list(spheres)
    other_spheres.remove(intersect_sphere)
    ray_direction = normalize(sub(light_pos, ray_origin))
    retval = ray_intersects_sphere(ray_origin, ray_direction, other_spheres)
    return retval[0] != None
 

def sphere_point(ray_origin, ray_direction, t):
    return plus(ray_origin, vec(ray_direction.x * t,
                            ray_direction.y * t,
                            ray_direction.z * t))


def render_image(pixels):
    p = 0
    ray_direction = vec(0.0, 0.0, 1.0)
    for i in range(reso):
        for j in range(reso):
                """if i == 70 and j == 75:
                    pixels[p - 2] = 0xFF
                    pixels[p - 3] = 0xFF
                    pixels[p - 1] = 0xFF
                """
                #if i == 58 and j == 72:
                #    pdb.set_trace()
                #x = j, y = i
                ray_origin = vec(pixel_coordinate_to_world_coordinate(j), pixel_coordinate_to_world_coordinate(i), 0.0)
                initial_color = (pixels[p + 2], pixels[p + 1], pixels[p])
                color = compute_color(initial_color, ray_origin, ray_direction)
                pixels[p] = color[2]
                pixels[p + 1] = color[1]
                pixels[p + 2] = color[0]
                """if i == 58 and j == 72:
                    pixels[p] = 255
                    pixels[p + 1] = 255
                    pixels[p + 2] = 255
                """
                p += 3


""" Returns color as tuple (r,g,b) """
def compute_color(initial_color, ray_origin, ray_direction):
    depth = 1
    reflect_coef = 1.0
    color = (0, 0, 0)
    sphere_self = None

    while depth < max_depth:
        other_spheres = list(spheres)
        if sphere_self:
            other_spheres.remove(sphere_self)
        point, sphere = ray_intersects_sphere(ray_origin, ray_direction, other_spheres)
        if point:
            n = normalize(sub(point, sphere))
            view_dir = normalize(sub(point, ray_origin))
            light_dir = normalize(sub(light_pos, point))
            blinn_dir = normalize(sub(light_dir, view_dir))
            blinn_term = max(0.0, spec_value * pow(dot(blinn_dir, n), spec_power))
            
            reflet = 2.0 * dot(ray_direction, n)
            reflet = vec(*map(lambda x: x * reflet, n))

            ray_origin = point
            #new direction is calculated from old using reflection
            ray_direction = sub(ray_direction, reflet)
            
            if not is_shadowed(point, sphere):
                new_color = calculate_lambert(point, sphere)
                color = map(lambda x, y: min(int(x + y * reflect_coef), 255), color, new_color)
                if blinn_term > 0.0:
                    color = map(lambda x: min(int(x = x + x * blinn_term), 255), color)
            reflect_coef *= 0.6
            sphere_self = sphere
        elif depth == 1:
            color = initial_color
            break
        else:
            break
        depth += 1
   
    return color
        

def targa_file():
    f = open('output.tga', 'wb')
    b = bytearray(reso * reso * 3)
    width = reso
    height = reso
    bpp = 24

    blue = 0
    p = 0

    # 1st pixel is bottom-left
    for i in range(reso):
        green = 0
        for j in range(reso):
            b[p] = blue
            p += 1
            b[p] = green
            p += 1
            b[p] = 0
            p += 1
            green += 255 / reso
        blue += 255 / reso 


    tga_header = bytearray([0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    assert(len(tga_header) == 18)
    tga_header[12] = width & 0xFF 
    tga_header[13] = (width >> 8) & 0xFF
    tga_header[14] = height & 0xFF
    tga_header[15] = (height >> 8) & 0xFF
    tga_header[16] = bpp
    f.write(tga_header)

    render_image(b)
    #pixel x = 64, y = 70
    #i = 128 * 3 * 70 + 3 * 64
    #b[i] = 0xFF
    f.write(b)


targa_file()         

