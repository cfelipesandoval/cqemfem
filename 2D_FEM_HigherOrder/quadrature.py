import numpy as np
import sys
import math

"""
GAUSS_QUADRATURE_TRISIMP calculates the integration over an arbitrary
triangle (with straight edges) using the Gaussian quadrature
integration rule
   val = gauss_quadrature_trisimp(p1(1),p1(2),p1(3),p2(1),p2(2),p2(3),p3(1),p3(2),p3(3),f)
   To use, first run the function gauss_quadrature_trisimp_init.  This function
   will then return 'loc' which will be a 3D array of points that
   the function to be integrated needs to be sampled at. After calculating
   these values put them in the same order into a 1D array and use this as
   the 'f' input and run the function.
   Inputs:
   p1(1), the x value of the first node
   p1(2), the y value of the first node
   p1(3), the z value of the first node
   Note: All other node points work similarly, going around the
   triangle in a counterclockwise direction from the first node
   f, the function to be integrated at the required sample points (see
   above) -- for 2D vectorized operation the first dimension should
   correspond to the different quadrature points
   
   Outputs:
   val, the result of the integration
   Accurate to the matlab version by 5 decimal places
"""

def gauss_quadrature_trisimp(p1, p2, p3, f):
    num_pts = len(f)
    # Make the weights array
    match num_pts:
        case 1:
            weights = 1

        case 3:
            weights = np.array([1 / 3, 1 / 3, 1 / 3])

        case 4:
            weights = np.array([-27 / 48, 25 / 48, 25 / 48, 25 / 48])
        case 6:
            weights = np.array(
                [
                    0.22338158967801,
                    0.22338158967801,
                    0.22338158967801,
                    0.10995174365532,
                    0.10995174365532,
                    0.10995174365532,
                ]
            )
        case 7:
            weights = np.array(
                [
                    0.22500000000000,
                    0.13239415278851,
                    0.13239415278851,
                    0.13239415278851,
                    0.12593918054483,
                    0.12593918054483,
                    0.12593918054483,
                ]
            )
        case 12:
            weights = np.array(
                [
                    0.11678627572638,
                    0.11678627572638,
                    0.11678627572638,
                    0.05084490637021,
                    0.05084490637021,
                    0.05084490637021,
                    0.08285107561837,
                    0.08285107561837,
                    0.08285107561837,
                    0.08285107561837,
                    0.08285107561837,
                    0.08285107561837,
                ]
            )
        case 13:
            weights = np.array(
                [
                    -0.14957004446768,
                    0.17561525743321,
                    0.17561525743321,
                    0.17561525743321,
                    0.05334723560884,
                    0.05334723560884,
                    0.05334723560884,
                    0.07711376089026,
                    0.07711376089026,
                    0.07711376089026,
                    0.07711376089026,
                    0.07711376089026,
                    0.07711376089026,
                ]
            )
        case 16:
            weights = np.array(
                [
                    0.14431560767779,
                    0.09509163426728,
                    0.09509163426728,
                    0.09509163426728,
                    0.10321737053472,
                    0.10321737053472,
                    0.10321737053472,
                    0.03245849762320,
                    0.03245849762320,
                    0.03245849762320,
                    0.02723031417443,
                    0.02723031417443,
                    0.02723031417443,
                    0.02723031417443,
                    0.02723031417443,
                    0.02723031417443,
                ]
            )
        case 19:
            weights = np.array(
                [
                    0.097135796282799,
                    0.031334700227139,
                    0.031334700227139,
                    0.031334700227139,
                    0.077827541004774,
                    0.077827541004774,
                    0.077827541004774,
                    0.079647738927210,
                    0.079647738927210,
                    0.079647738927210,
                    0.025577675658698,
                    0.025577675658698,
                    0.025577675658698,
                    0.043283539377289,
                    0.043283539377289,
                    0.043283539377289,
                    0.043283539377289,
                    0.043283539377289,
                    0.043283539377289,
                ]
            )
        case _:
            print(f"Incorrect number of points used")

    val = 0

    Ak = (
        math.sqrt(
            ((p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2]))** 2
            + ((p1[0] - p3[0]) * (p2[2] - p3[2]) - (p2[0] - p3[0]) * (p1[2] - p3[2]))**2
            + ((p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]))**2
        )
        / 2
    )
    # original triangle area
    for i in range(0, len(weights)):


        val += weights[i] * f[i]
        

    # val = (weights*f);
    val *= Ak

    return val

"""Compute Gaussian Quadrature
  @param a Start point
  @param b End point
  @param f Function to be computed 
"""
def gauss_quadrature_1D(start, stop, function):
  num_pts = len(function);

  if num_pts == 2:
      weights = [1,1]
  elif num_pts == 3:
      weights = [0.8888888888888888,0.5555555555555556,0.5555555555555556]
  elif num_pts == 4:
      weights = [0.6521451548625461,0.6521451548625461,0.3478548451374538, 
                0.3478548451374538]
  elif num_pts == 5:
      weights = [0.5688888888888889,0.4786286704993665,0.4786286704993665, 
                0.2369268850561891,0.2369268850561891]
  elif num_pts == 6:
      weights = [0.3607615730481386,0.3607615730481386,0.4679139345726910, 
                0.4679139345726910,0.1713244923791704,0.1713244923791704]
  elif num_pts == 7:
      weights = [0.4179591836734694,0.3818300505051189,0.3818300505051189, 
                0.2797053914892766,0.2797053914892766,0.1294849661688697, 
                0.1294849661688697]
  elif num_pts == 8:
      weights = [0.3626837833783620,0.3626837833783620,0.3137066458778873, 
                0.3137066458778873,0.2223810344533745,0.2223810344533745, 
                0.1012285362903763,0.1012285362903763]
  elif num_pts == 9:
      weights = [0.3302393550012598,0.1806481606948574,0.1806481606948574, 
                0.0812743883615744,0.0812743883615744,0.3123470770400029, 
                0.3123470770400029,0.2606106964029354,0.2606106964029354]
  elif num_pts == 10:
      weights = [0.2955242247147529,0.2955242247147529,0.2692667193099963, 
                0.2692667193099963,0.2190863625159820,0.2190863625159820, 
                0.1494513491505806,0.1494513491505806,0.0666713443086881, 
                0.0666713443086881]
  elif num_pts == 11:
      weights = [0.2729250867779006,0.2628045445102467,0.2628045445102467, 
                0.2331937645919905,0.2331937645919905,0.1862902109277343, 
                0.1862902109277343,0.1255803694649046,0.1255803694649046, 
                0.0556685671161737,0.0556685671161737]

  else:
    print('Error: "N" must be a value between 2 and 11.')
    sys.exit()

  weights = np.array(weights)

  val = ((stop - start) / 2) * weights @ function[:,1]

  return val

"""
GAUSS_QUADRATURE_TRISIMP_INIT calculates the integration points over an 
arbitrary triangle (with straight edges) using the Gaussian quadrature
integration rule
Inputs:
N, the order of the quadrature, gives the exact value of the
integration for a standard interval (triangle of (0,0),(1,0),(0,1)) for polynomials of
order 2N-1 (current version allows N = 1 to 9)
p1(1), the x value of the first node
p1(2), the y value of the first node
p1(3), the z value of the first node
Note: All other node points work similarly, going around the
triangle in a counterclockwise direction from the first node
Outputs:
loc, the points that the function to be integrated at needs to be
sampled at according to the integration rule
Notes:
Using mumpy arrays so everything is faster? (or at least more friendly)
used switch statements instead of if else
return loc array
"""
def gauss_quadrature_trisimp_init(
    N, p1, p2, p3
):  # let p1,p2,p3 be numpy arrays of length 3
    match N:
        case 1:
            alpha = 1 / 3
            beta = 1 / 3
            temp1 = p1[0] * (1 - alpha[i] - beta[i]) + p2[0] * alpha[i] + p3[0] * beta[i]
            temp2 = p1[1] * (1 - alpha[i]- beta[i]) + p2[1] * alpha[i] + p3[1] * beta[i]
            temp3 = p1[2] * (1 - alpha[i] - beta[i]) + p2[2] * alpha[i] + p3[2] * beta[i]

            return np.array([temp1, temp2, temp3])  # loc

        case 2:
            alpha = np.array([1 / 6, 2 / 3, 1 / 6])
            beta = np.array([1 / 6, 1 / 6, 2 / 3])
            loc = np.empty((3, 3))  # empty 3x3 loc array

            for i in range(0, len(alpha)):
                temp1 = p1[0] * (1 - alpha[i] - beta[i]) + p2[0] * alpha[i] + p3[0] * beta[i]
                temp2 = p1[1] * (1 - alpha[i]- beta[i]) + p2[1] * alpha[i] + p3[1] * beta[i]
                temp3 = p1[2] * (1 - alpha[i] - beta[i]) + p2[2] * alpha[i] + p3[2] * beta[i]
                loc[i][0] = temp1
                loc[i][1] = temp2
                loc[i][2] = temp3
            return loc

        case 3:
            alpha = np.array([1 / 3, 1 / 5, 1 / 5, 3 / 5])
            beta = np.array([1 / 3, 3 / 5, 1 / 5, 1 / 5])
            loc = np.empty((4, 3))  # empty 4x3 loc array

            for i in range(0, len(alpha)):
                temp1 = p1[0] * (1 - alpha[i] - beta[i]) + p2[0] * alpha[i] + p3[0] * beta[i]
                temp2 = p1[1] * (1 - alpha[i]- beta[i]) + p2[1] * alpha[i] + p3[1] * beta[i]
                temp3 = p1[2] * (1 - alpha[i] - beta[i]) + p2[2] * alpha[i] + p3[2] * beta[i]
                loc[i][0] = temp1
                loc[i][1] = temp2
                loc[i][2] = temp3
            return loc

        case 4:
            alpha = np.array(
                [
                    0.44594849091597,
                    0.44594849091597,
                    0.10810301816807,
                    0.09157621350977,
                    0.09157621350977,
                    0.81684757298046,
                ]
            )
            beta = np.array(
                [
                    0.44594849091597,
                    0.10810301816807,
                    0.44594849091597,
                    0.09157621350977,
                    0.81684757298046,
                    0.09157621350977,
                ]
            )
            loc = np.empty((6, 3))  # empty 6x3 loc array

            for i in range(0, len(alpha)):
                temp1 = p1[0] * (1 - alpha[i] - beta[i]) + p2[0] * alpha[i] + p3[0] * beta[i]
                temp2 = p1[1] * (1 - alpha[i]- beta[i]) + p2[1] * alpha[i] + p3[1] * beta[i]
                temp3 = p1[2] * (1 - alpha[i] - beta[i]) + p2[2] * alpha[i] + p3[2] * beta[i]
                loc[i][0] = temp1
                loc[i][1] = temp2
                loc[i][2] = temp3
            return loc

        case 5:
            alpha = np.array(
                [
                    0.33333333333333,
                    0.47014206410511,
                    0.47014206410511,
                    0.05971587178977,
                    0.10128650732346,
                    0.10128650732346,
                    0.79742698535309,
                ]
            )
            beta = np.array(
                [
                    0.33333333333333,
                    0.47014206410511,
                    0.05971587178977,
                    0.47014206410511,
                    0.10128650732346,
                    0.79742698535309,
                    0.10128650732346,
                ]
            )
            loc = np.empty((7, 3))  # empty 7x3 loc array

            for i in range(0, len(alpha)):
                temp1 = p1[0] * (1 - alpha[i] - beta[i]) + p2[0] * alpha[i] + p3[0] * beta[i]
                temp2 = p1[1] * (1 - alpha[i]- beta[i]) + p2[1] * alpha[i] + p3[1] * beta[i]
                temp3 = p1[2] * (1 - alpha[i] - beta[i]) + p2[2] * alpha[i] + p3[2] * beta[i]
                loc[i][0] = temp1
                loc[i][1] = temp2
                loc[i][2] = temp3
            return loc

        case 6:
            alpha = np.array(
                [
                    0.24928674517091,
                    0.24928674517091,
                    0.50142650965818,
                    0.06308901449150,
                    0.06308901449150,
                    0.87382197101700,
                    0.31035245103378,
                    0.63650249912140,
                    0.05314504984482,
                    0.63650249912140,
                    0.31035245103378,
                    0.05314504984482,
                ]
            )
            beta = np.array(
                [
                    0.24928674517091,
                    0.50142650965818,
                    0.24928674517091,
                    0.06308901449150,
                    0.87382197101700,
                    0.06308901449150,
                    0.63650249912140,
                    0.05314504984482,
                    0.31035245103378,
                    0.31035245103378,
                    0.05314504984482,
                    0.63650249912140,
                ]
            )
            loc = np.empty((12, 3))  # empty 12x3 loc array

            for i in range(0, len(alpha)):
                temp1 = p1[0] * (1 - alpha[i] - beta[i]) + p2[0] * alpha[i] + p3[0] * beta[i]
                temp2 = p1[1] * (1 - alpha[i]- beta[i]) + p2[1] * alpha[i] + p3[1] * beta[i]
                temp3 = p1[2] * (1 - alpha[i] - beta[i]) + p2[2] * alpha[i] + p3[2] * beta[i]
                loc[i][0] = temp1
                loc[i][1] = temp2
                loc[i][2] = temp3
            return loc

        case 7:
            alpha = np.array(
                [
                    0.33333333333333,
                    0.26034596607904,
                    0.26034596607904,
                    0.47930806784192,
                    0.06513010290222,
                    0.06513010290222,
                    0.86973979419557,
                    0.31286549600487,
                    0.63844418856981,
                    0.04869031542532,
                    0.63844418856981,
                    0.31286549600487,
                    0.04869031542532,
                ]
            )
            beta = np.array(
                [
                    0.33333333333333,
                    0.26034596607904,
                    0.47930806784192,
                    0.26034596607904,
                    0.06513010290222,
                    0.86973979419557,
                    0.06513010290222,
                    0.63844418856981,
                    0.04869031542532,
                    0.31286549600487,
                    0.31286549600487,
                    0.04869031542532,
                    0.63844418856981,
                ]
            )
            loc = np.empty((13, 3))  # empty 13x3 loc array

            for i in range(0, len(alpha)):
                temp1 = p1[0] * (1 - alpha[i] - beta[i]) + p2[0] * alpha[i] + p3[0] * beta[i]
                temp2 = p1[1] * (1 - alpha[i]- beta[i]) + p2[1] * alpha[i] + p3[1] * beta[i]
                temp3 = p1[2] * (1 - alpha[i] - beta[i]) + p2[2] * alpha[i] + p3[2] * beta[i]
                loc[i][0] = temp1
                loc[i][1] = temp2
                loc[i][2] = temp3
            return loc

        case 8:
            alpha = np.array(
                [
                    0.33333333333333,
                    0.45929258829272,
                    0.45929258829272,
                    0.08141482341455,
                    0.17056930775176,
                    0.17056930775176,
                    0.65886138449648,
                    0.05054722831703,
                    0.05054722831703,
                    0.89890554336594,
                    0.26311282963464,
                    0.72849239295540,
                    0.00839477740996,
                    0.72849239295540,
                    0.26311282963464,
                    0.00839477740996,
                ]
            )
            beta = np.array(
                [
                    0.33333333333333,
                    0.45929258829272,
                    0.08141482341455,
                    0.45929258829272,
                    0.17056930775176,
                    0.65886138449648,
                    0.17056930775176,
                    0.05054722831703,
                    0.89890554336594,
                    0.05054722831703,
                    0.72849239295540,
                    0.00839477740996,
                    0.26311282963464,
                    0.26311282963464,
                    0.00839477740996,
                    0.72849239295540,
                ]
            )
            loc = np.empty((16, 3))  # empty 16x3 loc array

            for i in range(0, len(alpha)):
                temp1 = p1[0] * (1 - alpha[i] - beta[i]) + p2[0] * alpha[i] + p3[0] * beta[i]
                temp2 = p1[1] * (1 - alpha[i]- beta[i]) + p2[1] * alpha[i] + p3[1] * beta[i]
                temp3 = p1[2] * (1 - alpha[i] - beta[i]) + p2[2] * alpha[i] + p3[2] * beta[i]
                loc[i][0] = temp1
                loc[i][1] = temp2
                loc[i][2] = temp3
            return loc

        case 9:
            alpha = np.array(
                [
                    0.333333333333333,
                    0.020634961602525,
                    0.489682519198738,
                    0.489682519198738,
                    0.125820817014127,
                    0.437089591492937,
                    0.437089591492937,
                    0.623592928761935,
                    0.188203535619033,
                    0.188203535619033,
                    0.910540973211095,
                    0.044729513394453,
                    0.044729513394453,
                    0.036838412054736,
                    0.221962989160766,
                    0.036838412054736,
                    0.741198598784498,
                    0.221962989160766,
                    0.741198598784498,
                ]
            )
            beta = np.array(
                [
                    0.333333333333333,
                    0.489682519198738,
                    0.020634961602525,
                    0.489682519198738,
                    0.437089591492937,
                    0.125820817014127,
                    0.437089591492937,
                    0.188203535619033,
                    0.623592928761935,
                    0.188203535619033,
                    0.044729513394453,
                    0.910540973211095,
                    0.044729513394453,
                    0.221962989160766,
                    0.036838412054736,
                    0.741198598784498,
                    0.036838412054736,
                    0.741198598784498,
                    0.221962989160766,
                ]
            )
            loc = np.empty((19, 3))  # empty 16x3 loc array

            for i in range(0, len(alpha)):
                temp1 = p1[0] * (1 - alpha[i] - beta[i]) + p2[0] * alpha[i] + p3[0] * beta[i]
                temp2 = p1[1] * (1 - alpha[i]- beta[i]) + p2[1] * alpha[i] + p3[1] * beta[i]
                temp3 = p1[2] * (1 - alpha[i] - beta[i]) + p2[2] * alpha[i] + p3[2] * beta[i]
                loc[i][0] = temp1
                loc[i][1] = temp2
                loc[i][2] = temp3
            return loc

        case _:
            print(f"N must be a value between 1 and 9")
            return