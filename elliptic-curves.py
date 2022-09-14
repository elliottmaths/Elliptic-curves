import math
import numpy as np

def sum_of_points(a, b, P, Q):
    # Returns the point P+Q on the elliptic curve y^2 = x^3 + ax + b
    # (given points P and Q on the curve, and where + is the group operation)
    if P == 'id':
        return Q
    if Q == 'id':
        return P
    if P[0] == Q[0]:
        if P[1] == -Q[1]:
            return 'id'
        else:
            s = (3*P[0]**2 + a) / (2*P[1])
            t = P[1] - s * P[0]
            x_3 = s**2 - P[0] - Q[0]
            y_3 = s * x_3 + t
            return (x_3, -y_3)
    s = (Q[1] - P[1]) / (Q[0] - P[0])
    t = P[1] - s * P[0]
    x_3 = s**2 - P[0] - Q[0]
    y_3 = s * x_3 + t
    return (x_3, -y_3)

def discriminant(a, b):
    # Returns discriminant of elliptic curve y^2 = x^3 + ax + b
    return 4 * a**3 + 27 * b ** 2

def divisorGenerator(n):
    if n < 0:
        n = -n
    large_divisors = []
    for i in range(1, int(math.sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i*i != n:
                large_divisors.append(n / i)
    for divisor in reversed(large_divisors):
        yield int(divisor)

def divisors(n):
    # Returns a list containing precisely the divisors of the natural number n
    return list(divisorGenerator(n))

def find_x_given_y(a, b, y):
    # Returns a list of the rational solutions x to the cubic y^2 = x^3 + ax + b
    polynomial = [1, 0, a, b - y**2]
    rational_solns = []
    for root in np.roots(polynomial):
        if root.imag == 0:
            root_int = round(root.real)
            test = root_int**3 + a*root_int + b - y**2
            if test == 0:
                rational_solns.append(root_int)
    return rational_solns

def find_x_given_y_2(a, b, y):
    # Returns a list of the solutions x (rational or otherwise)
    # to the cubic y^2 = x^3 + ax + b
    polynomial = [1, 0, a, b - y**2]
    rational_solns = []
    for root in np.roots(polynomial):
        if root.imag == 0:
            root_int = round(root.real)
            test = root_int**3 + a*root_int + b - y**2
            if test == 0:
                rational_solns.append(root_int)
    return rational_solns

def NagellLutzList(a, b):
    # Returns, for an elliptic curve y^2 = x^3 + ax + b, the list of points (x,y)
    # such that y = 0 or y^2 divides the discriminant
    disc = discriminant(a, b)
    factors = divisors(disc)
    rational_points = []
    x_order_two_roots = np.roots([1, 0, a, b])
    for x in x_order_two_roots:
        if x.imag == 0:
            x_int = round(x.real)
            test = x_int**3 + a*x_int + b
            if test == 0 and (x_int, 0) not in rational_points:
                rational_points.append((x_int, 0))
    possible_y = []
    for factor in factors:
        square_root = math.sqrt(factor)
        if square_root == int(square_root):
            possible_y.append(int(square_root))
            possible_y.append(-int(square_root))
    for y in possible_y:
        possible_x = find_x_given_y(a, b, y)
        for x in possible_x:
            if (x, y) not in rational_points:
                rational_points.append((x, y))
    return rational_points

def multiple_of_P(a, b, n, P):
    # Returns the point nP, for a point P on the elliptic curve
    # y^2 = x^3 + ax + b, and an integer n
    Q = P
    for aux in range(n-1):
        Q = sum_of_points(a, b, P, Q)
    return Q

def isFiniteOrder(a, b, P):
    # Returns True if P is finite order, and false otherwise, where P is
    # a point on the elliptic curve y^2 = x^3 + ax + b
    if P == 'id':
        return True
    Q = P
    for aux in range(11):
        Q = sum_of_points(a, b, P, Q)
        if Q == 'id':
            return True
        if Q[0] != int(Q[0]) or Q[1] != int(Q[1]):
            return False
    return False

def torsion_subgroup(a, b):
    # Returns the torsion subgroup of the group of rational points
    # of the elliptic curve y^2 = x^3 + ax + b
    torsion_points = ['id']
    for candidate in NagellLutzList(a, b):
        if isFiniteOrder(a, b, candidate):
            torsion_points.append(candidate)
    return torsion_points

def order_of_P(a, b, P):
    # Returns the order of a point P on the elliptic curve y^2 = x^3 + ax + b,
    # unless P has infinite order, in which case it returns False
    order = 1
    if P == 'id':
        return order
    Q = P
    for aux in range(11):
        order += 1
        Q = sum_of_points(a, b, P, Q)
        if Q == 'id':
            return order
        if Q[0] != int(Q[0]) or Q[1] != int(Q[1]):
            return False
    return False

def solns_mod_p(a, b, p):
    # Returns the solutions of y^2 = x^3 + ax + b modulo a prime p
    a_red = a % p
    b_red = b % p
    solns = ['id']
    for x_poss in range(p):
        for y_poss in range(p):
            if (y_poss**2) % p == (x_poss**3 + a_red*x_poss + b_red) % p:
                solns.append((x_poss, y_poss))
    return solns

# The below is LMFDB 272.b1
a = -1451
b = 21274
P = (21, 8)
Q = P
