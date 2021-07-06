import sympy
import sympy as sp
from sympy import log
from termcolor import colored, cprint
import math
import sys


def num_after_point(x):
    s = str(x)
    if '.' not in s:
        return 0
    return len(s) - s.index('.') - 1


def checkResult(f, res, epsilon):
    x = round(res, num_after_point(epsilon + 1))
    try:
        abs(f(x))
        return True
    except:
        return False


def mullersMethod(p0, p1, p2, f, epsilon, maxIter):
    """
    :param p0: The first guess
    :param p1: The second guess
    :param p2: The third guess
    :param f: The original function
    :param epsilon: The accepted error
    :param maxIter: The max number of iterations that we are willing to accept
    :return: the amount of the roots found
    """
    countRoots = 0
    start = p0
    end = p2
    s = "The initial guess:\np0 = " + str(p0) + ", p1 = " + str(p1) + ", p2 = " + str(p2) + "\n"
    # Check if one of the initial guess is the root
    flag = False
    maxTries = 10
    move = 0
    while move < maxTries:
        try:
            if f(p0) == 0:
                print(s)
                print(colored("The root is: x = " + str(p0) + ", found after 0 iterations", 'blue'))
                countRoots += 1
                flag = True
            break
        except:
            # print("catch p0")
            move += 1
            p0 += min(0.0001, epsilon)
    if maxTries is move:
        # print("Error - function is undefined in this range :(")
        return countRoots
    move = 0
    while move < maxTries:
        try:
            if f(p1) == 0:
                print(s)
                print(colored("The root is: x = " + str(p1) + ", found after 1 iterations", 'blue'))
                countRoots += 1
                flag = True
            break
        except:
            # print("catch p1")
            move += 1
            p1 += min(0.0001, epsilon)
    if maxTries is move:
        # print("Error - function is undefined in this range :(")
        return countRoots

    move = 0
    while move < maxTries:
        try:
            if f(p2) == 0:
                print(s)
                print(colored("The root is: x = " + str(p2) + ", found after 2 iterations", 'blue'))
                countRoots += 1
                flag = True
            break
        except:
            # print("catch p2")
            move += 1
            p2 -= min(0.0001, epsilon)
    if maxTries is move:
        # print("Error - function is undefined in this range :(")
        return countRoots

    if flag is True:
        return countRoots

    # the iterations
    h1 = p1 - p0
    h2 = p2 - p1
    grad1 = (f(p1) - f(p0)) / h1  # The gradient between p0 and p1
    grad2 = (f(p2) - f(p1)) / h2  # The gradient between p1 and p2
    d = (grad2 - grad1) / (h2 + h1)
    iter = 3  # Start to count the iteration
    while iter <= maxIter:  # while we did not reach the maximum iterations
        b = grad2 + (h2 * d)
        try:
            if (b ** 2 - (
                    4 * f(p2) * d)) < 0:  # if the discriminant is negative - cant do square root on negative number
                # print("*** There is no real root ***")
                return countRoots
        except:
            return countRoots
        D = (b ** 2 - (4 * f(p2) * d)) ** 0.5  # the square root in the denominator
        if abs(b - D) < abs(b + D):  # Save the biggest result
            E = b + D
        else:
            E = b - D
        h = (-2) * f(p2) / E
        p = p2 + h  # The new guess
        try:
            f(p)
        except:
            return countRoots
        if abs(h) < epsilon:  # if the difference is smaller then epsilon
            # if (p > start) and (p < end):
            if (p > start) and (p < end) and (checkResult(f, p, epsilon) is True):  # if the point is in the range
                if abs(f(p)) <= 0 + epsilon:
                    print(s, end="")
                    print(
                        colored("The root is: x = " + str(p) + ", found after " + str(iter) + " iterations\n", 'blue'))
                    countRoots += 1
                    return countRoots
                else:
                    return countRoots
            else:
                # print("*** There is no root in this range ***")
                return countRoots
        s += ("iteration = " + str(iter) + ", x = " + str(p) + ", f(x) = " + str(f(p)) + "\n")
        # Preparing the variable for the next iteration
        p0 = p1
        p1 = p2
        p2 = p
        h1 = p1 - p0
        h2 = p2 - p1
        grad1 = (f(p1) - f(p0)) / h1
        grad2 = (f(p2) - f(p1)) / h2
        d = (grad2 - grad1) / (h2 + h1)
        iter += 1  # increase the number of iteration
    # print("Method failed after " + str(maxIter) + "iterations")
    return countRoots


def rangeDivision(polinom, epsilon, iterNumber, maxRoots=15, start_point=-100, end_point=100):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The accepted error
    :param iterNumber: The max number of iterations that we are willing to accept
    :return: None
    """

    # dividing the range to small parts
    rootsCount = 0
    end = start_point + 0.1
    middle = start_point + 0.05
    while end <= end_point and rootsCount < maxRoots:  # while we don't reach to the end point of the range
        # print("range: [" + str(start_point) + ", " + str(end) + "]")
        rootsCount += mullersMethod(start_point, middle, end, polinom, epsilon,
                                    iterNumber)  # calculate by muller`s method
        start_point = end  # increase the start point of the range
        middle = (start_point + 0.05)  # set the middle point of the range
        end = start_point + 0.1  # set the end point of the range
    print("The number of roots found is: " + str(rootsCount))


def driver():
    x = sp.symbols('x')
    f = x ** 4 - 3 * (x ** 3) + x ** 2 + x + 1
    f = (x ** 4 - 4 * (x ** 2)) / x
    #f = (x ** 3 + 1) ** 15
    # f = ((x - x ** 3) ** 0.5) / ((math.exp(1)) ** x)
    # f = (x) ** 0.5
    # y = sp.lambdify(f, x)
    # print(y)
    # f = (math.exp(1))**x
    f = sympy.cos(x)
    #f = (10 * x * ((math.exp(1) ** x)) - (x ** 3)) / (x**3 - 5)
    #f = log(x ** 4, 10) + x + 1
    epsilon = 10 ** -5
    startPoint = 0.5
    endPoint = 3
    maxIteration = 100
    f = sp.lambdify(x, f)
    rangeDivision(f, epsilon, maxIteration)
    print("program finished :)")


driver()
