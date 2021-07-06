import sympy
import sympy as sp
from sympy import log
from termcolor import colored, cprint
import math
import sys


def num_after_point(x):
    """
    :param x: a number
    :return: the number of digits after the point in the number
    """
    s = str(x)   # cast the number to string
    if '.' not in s:    # search for a dot in the number
        return 0    # no dot was found
    return len(s) - s.index('.') - 1    # return the number of digits after the point


def checkResult(f, res, epsilon):
    """
    :param f: a function
    :param res: an x value
    :param epsilon: an epsilon
    :return: if res (within a range of epsilon) is in the definition range of f
    """
    x = round(res, num_after_point(epsilon + 1))   # calc x by rounding res
    try:
        f(x)   # check if x is in the defenition range of f
        return True   # return true if it does
    except:
        return False   # return false if it doesn't


def mullersMethod(p0, p1, p2, f, epsilon, maxIter):
    """
    :param p0: The first guess
    :param p1: The second guess
    :param p2: The third guess
    :param f: The original function
    :param epsilon: The accepted error
    :param maxIter: The max number of iterations that we are willing to accept
    :return: the amount of the roots found in the range [p0, p2]
    """
    s = "The initial guess:\np0 = " + str(p0) + ", p1 = " + str(p1) + ", p2 = " + str(
        p2) + "\n"  # print the intial guess for the range
    countRoots = 0  # the number of roots that were found in the range
    tempRes = []   # a list that will keep
    start = p0  # will keep the original start point of the range
    end = p2  # will keep the original end point of the range
    # Check if one of the initial guesses is a root of the function
    flag = False  # will save if one of the initial guesses is a root of the function
    maxTries = 10  # the max number of tries to fix an initial guess in case it's not in the definition area
    move = 0  # the number of fixes made for a guess
    while move < maxTries:
        try:
            if f(p0) == 0:
                print(s)
                print(colored("The root is: x = " + str(p0) + ", f(x) = " + str(f(p0)) + ", found after 0 iterations",
                              'blue'))
                countRoots += 1
                tempRes.append(p0)
                flag = True
            break
        except:
            move += 1
            p0 += min(0.0001, epsilon)
    if maxTries is move:
        # print("Error - function is undefined in this range :("
        return countRoots, tempRes
    move = 0
    while move < maxTries:
        try:
            if f(p1) == 0:
                print(s)
                print(colored("The root is: x = " + str(p1) + ", f(x) = " + str(f(p1)) + ", found after 1 iterations",
                              'blue'))
                countRoots += 1
                tempRes.append(p1)
                flag = True
            break
        except:
            move += 1
            p1 += min(0.0001, epsilon)
    if maxTries is move:
        # print("Error - function is undefined in this range :(")
        return countRoots, tempRes

    move = 0
    while move < maxTries:
        try:
            if f(p2) == 0:
                print(s)
                print(colored("The root is: x = " + str(p2) + ", f(x) = " + str(f(p2)) + ", found after 2 iterations",
                              'blue'))
                countRoots += 1
                tempRes.append(p2)
                flag = True
            break
        except:
            move += 1
            p2 -= min(0.0001, epsilon)
    if maxTries is move:
        # print("Error - function is undefined in this range :(")
        return countRoots, tempRes

    if flag is True:
        return countRoots, tempRes

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
                return countRoots, tempRes
        except:
            # print("overflow")
            return countRoots, tempRes
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
            return countRoots, tempRes
        if abs(h) < epsilon:  # if the difference is smaller then epsilon
            # if (p > start) and (p < end):
            if (p > start) and (p < end) and (checkResult(f, p, epsilon) is True):  # if the point is in the range
                if abs(f(p)) <= 0 + epsilon:
                    print(s, end="")
                    print(
                        colored("The root is: x = " + str(p) + ", f(x) = " + str(f(p)) + ", found after " + str(
                            iter) + " iterations\n", 'blue'))
                    countRoots += 1
                    tempRes.append(p)
                    return countRoots, tempRes
                else:
                    return countRoots, tempRes
            else:
                # print("*** There is no root in this range ***")
                return countRoots, tempRes
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
    return countRoots, tempRes


def rangeDivision(polinom, epsilon, iterNumber, maxRoots, start_point, end_point):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The accepted error
    :param iterNumber: The max number of iterations that we are willing to accept
    :param maxRoots: The maximum number of roots the algorithm will find
    :return: a list of the function's roots in the range [start_point, end_point]
    """

    # dividing the range to small parts
    results = []
    rootsCount = 0
    end = start_point + 0.1
    middle = start_point + 0.05
    while end <= end_point and rootsCount < maxRoots:  # while we don't reach to the end point of the range
        # print("range: [" + str(start_point) + ", " + str(end) + "]")
        tempCount, tempRes = mullersMethod(start_point, middle, end, polinom, epsilon,
                                           iterNumber)  # calculate by muller`s method
        rootsCount += tempCount
        results.extend(tempRes)
        start_point = end  # increase the start point of the range
        middle = (start_point + 0.05)  # set the middle point of the range
        end = start_point + 0.1  # set the end point of the range
    print("The number of roots found is: " + str(rootsCount))
    return results



def checkDifferPart3(l, d, name1, name2, epsilon):
    print("check the difference between the roots:")
    # flag = True
    for i in range(len(l)):
        print("x" + str(i + 1) + " - " + name1 + ": " + str(l[i])
              + "\nx" + str(i + 1) + " - " + name2 + ": " + str(d[i]))
        if abs(l[i] - d[i]) > epsilon:
            # flag = False
            print("The difference is bigger than epsilon for some of the components\n")
            return
    print("The difference is smaller than epsilon for all the components\n")


def driver():
    x = sp.symbols('x')
    f = x ** 4 - 3 * (x ** 3) + x ** 2 + x + 1
    f = (x ** 4 - 4 * (x ** 2)) / x
    # f = (x ** 3 + 1) ** 15
    # f = ((x - x ** 3) ** 0.5) / ((math.exp(1)) ** x)
    # f = (x) ** 0.5
    # y = sp.lambdify(f, x)
    # print(y)
    # f = (math.exp(1))**x
    # f = sympy.cos(x)
    # f = (10 * x * ((math.exp(1) ** x)) - (x ** 3)) / (x**3 - 5)
    f = log(x ** 4, 10) + x + 1
    epsilon = 10 ** -5
    maxRoots = 15
    start_point = -100
    end_point = 100
    maxIteration = 100
    f = sp.lambdify(x, f)
    results = rangeDivision(f, epsilon, maxIteration, maxRoots, start_point, end_point)
    checkDifferPart3(results, [0.1, 0.2, 0.3], "Muller's", "Expected", 10 ** -5)
    print("program finished :)")


driver()


