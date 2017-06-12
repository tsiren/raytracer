import pdb
from collections import namedtuple
from math import sqrt
from math import floor


def normalize(v):
    vlen = v.len()
    if vlen == 0.0:
        return v
    return vec(v.x / vlen, v.y / vlen, v.z / vlen)


def dist(v1, v2):
    return sqrt((v1.x - v2.x)**2 + (v1.y - v2.y)**2 + (v1.z - v2.z)**2)


def mul(v, s):
    return vec(v.x * s, v.y * s, v.z * s)


class vec(namedtuple('vec', "x y z")):

    """Vector / point in 3D space."""

    def __add__(self, other):
        if other.__class__ != self.__class__:
            raise TypeError
        return self.__class__(*[a + b for a, b in zip(self, other)])

    def __sub__(self, other):
        if other.__class__ != self.__class__:
            raise TypeError
        return self.__class__(*[a - b for a, b in zip(self, other)])

    def __mul__(self, other):
        """dot product"""
        if other.__class__ != self.__class__:
            raise TypeError
        return self.x * other.x + self.y * other.y + self.z * other.z

    def len(self):
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)


max_depth = 10

shape = namedtuple('shape', "type point color normal")
shapes = []
shapes.append(shape('sphere', vec(0.0, 0.45, 1.0), (255, 0, 0), 0))
shapes.append(shape('sphere', vec(0.0, -0.45, 1.5), (245, 240, 82), 0))
# y-coord plane, last tuple is plane normal
shapes.append(
    shape('plane', vec(0.0, -1.0, 0.0), (0, 0, 255), vec(0.0, 1.0, 0.0)))

spec_value = 5.0
spec_power = 100

sphere_radius = 0.5
light_pos = vec(0.5, 0.5, 0.0)
reso = 255


def calculate_lambert(intersect, shape):
    """Lambert is defined as
    dot product between surface normal & vector towards light source.
    """

    # light_pos - intersect = vector from intersection to light start
    light_dir = normalize(light_pos - intersect)

    # intersect - sphere_center = vector from sphere center to intersection
    # --> sphere normal
    normal = None
    if shape.type == 'sphere':
        normal = normalize(intersect - shape.point)
    elif shape.type == 'plane':
        normal = shape.normal
    else:
        raise TypeError('invalid type', shape.type)

    lambert = max(0.0, normal * light_dir)
    if lambert == 0.0:
        return (0, 0, 0)
    else:
        color = shape.color
        if shape.type == 'plane':
            color = plane_color(intersect)
        color = map(lambda x: x * lambert, color)
        return color


def plane_color(point):
    tmp = floor(point.x) + floor(point.z)
    return (0, 0, 0) if tmp % 2 == 1 else (255, 255, 255)


def pixel_coordinate_to_world_coordinate(coordinate):
    return ((coordinate / float(reso)) - 0.5) * 2.0


def ray_intersects(ray_origin, ray_direction, shape_self):
    intersect = None
    closest_shape = None
    t_min = 0.0
    t = 0.0
    for shape in shapes:
        tmp = None
        if shape == shape_self:
            continue
        elif shape.type == 'plane':
            tmp, t = ray_intersects_plane(
                ray_origin, ray_direction, shape.point, shape.normal)
        elif shape.type == 'sphere':
            tmp, t = ray_intersects_sphere(
                ray_origin, ray_direction, shape.point)
        else:
            raise TypeError('invalid type', shape.type)
        if tmp and (t_min > t or t_min == 0.0):
            intersect = tmp
            closest_shape = shape
            t_min = t
    return intersect, closest_shape


def ray_intersects_plane(ray_origin, ray_direction, plane_point, plane_normal):
    denom = ray_direction * plane_normal
    if (abs(denom) > 0.00001):
        t = (plane_point - ray_origin) * plane_normal / denom
        if (t >= 0):
            return mul(ray_direction, t), t
    return False, 0.0


def ray_intersects_sphere(ray_origin, ray_direction, sphere_center):
    """https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    simplified 2nd degree polynomial

    D discrimant (term under square root)
    dist -> distance from sphere center to ray origin

    grabs the closest intersection
    this is really line intersection instead of vector,
    thats why we only take positive t!
    """
    point = None

    dist = ray_origin - sphere_center
    b = ray_direction * dist
    c = dist * dist - sphere_radius ** 2
    D = b * b - c
    if D >= 0:
        t0 = -b - sqrt(D)
        t1 = -b + sqrt(D)
        t = min(t0, t1)
        # negative t goes away from light position -> param d in line
        # equation, distance along line from starting point
        if t > 0.0:
            point = sphere_point(ray_origin, ray_direction, t)
            return point, t

    return point, 0.0


def is_shadowed(ray_origin, shape_self):
    ray_direction = normalize(light_pos - ray_origin)
    point = None
    dist_light = dist(ray_origin, light_pos)
    for shape in shapes:
        if shape == shape_self:
            continue
        if shape.type == 'plane':
            point = ray_intersects_plane(
                ray_origin, ray_direction, shape.point, shape.normal)
        elif shape.type == 'sphere':
            point = ray_intersects_sphere(
                ray_origin, ray_direction, shape.point)
        else:
            raise TypeError('invalid type', shape.type)
        if point is not None and point[0] and point[1] < dist_light:
            return True
    return False


def sphere_point(ray_origin, ray_direction, t):
    return ray_origin + vec(ray_direction.x * t,
                            ray_direction.y * t,
                            ray_direction.z * t)


def render_image(pixels):
    p = 0
    ray_origin = vec(0.0, 0.0, 0.0)
    for i in range(reso):
        for j in range(reso):
            ray_direction = normalize(
                vec(pixel_coordinate_to_world_coordinate(j),
                    pixel_coordinate_to_world_coordinate(i),
                    1.0))
            initial_color = (pixels[p + 2], pixels[p + 1], pixels[p])
            color = compute_color(initial_color, ray_origin, ray_direction)
            pixels[p] = color[2]
            pixels[p + 1] = color[1]
            pixels[p + 2] = color[0]
            p += 3


def compute_color(initial_color, ray_origin, ray_direction):
    """Returns color as tuple (r,g,b)."""
    depth = 1
    reflect_coef = 1.0
    color = (0, 0, 0)
    shape_self = None

    while depth < max_depth:
        point, shape = ray_intersects(ray_origin, ray_direction, shape_self)

        if point:
            n = None
            if shape.type == 'plane':
                if depth == 1:
                    return plane_color(point)
                n = shape.normal
            else:
                n = normalize(point - shape.point)
            view_dir = normalize(point - ray_origin)
            light_dir = normalize(light_pos - point)
            blinn_dir = normalize(light_dir - view_dir)
            blinn_term = max(0.0, spec_value * pow(blinn_dir * n, spec_power))

            reflet = 2.0 * (ray_direction * n)
            reflet = vec(*map(lambda x: x * reflet, n))

            ray_origin = point
            # new direction is calculated from old using reflection
            # see http://www.3dkingdoms.com/weekly/weekly.php?a=2
            ray_direction = ray_direction - reflet
            if not is_shadowed(ray_origin, shape):
                new_color = calculate_lambert(point, shape)
                color = map(lambda x, y: min(
                    int(x + y * reflect_coef), 255), color, new_color)
                if blinn_term > 0.0:
                    color = map(lambda x: min(
                        int(x=x + x * blinn_term), 255), color)
            reflect_coef *= 0.6
            shape_self = shape
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

    # 1st pixel is bottom-left, last is top-right
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

    tga_header = bytearray(
        [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert(len(tga_header) == 18)
    tga_header[12] = width & 0xFF
    tga_header[13] = (width >> 8) & 0xFF
    tga_header[14] = height & 0xFF
    tga_header[15] = (height >> 8) & 0xFF
    tga_header[16] = bpp
    f.write(tga_header)

    render_image(b)
    f.write(b)


targa_file()
