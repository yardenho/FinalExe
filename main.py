import sympy
import sympy as sp
# from colorama import Fore, Back, Style
from sympy import log
from termcolor import colored, cprint
import math


def num_after_point(x):
    """
    :param x: a number
    :return: the number of digits after the point in the number
    """
    s = str(x)  # cast the number to string
    if '.' not in s:  # search for a dot in the number
        return 0  # no dot was found
    return len(s) - s.index('.') - 1  # return the number of digits after the point


def checkResult(f, res, epsilon):
    """
    :param f: a function
    :param res: an x value
    :param epsilon: an epsilon
    :return: if res (within a range of epsilon) is in the definition range of f
    """
    x = round(res, num_after_point(epsilon + 1))  # calc x by rounding res
    try:
        f(x)  # check if x is in the definition range of f
        return True  # return true if it does
    except:
        return False  # return false if it doesn't


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
    tempRes = []  # a list that will keep
    start = p0  # will keep the original start point of the range
    end = p2  # will keep the original end point of the range
    # Check if one of the initial guesses is a root of the function
    flag = False  # will save if one of the initial guesses is a root of the function
    maxTries = 10  # the max number of tries to fix an initial guess in case it's not in the definition area
    move = 0  # the number of fixes made for a guess
    while move < maxTries:  # the max number of tries to fix p0 in case it's not in f definition range
        try:  # try to substitute p0 in f
            if f(p0) == 0:  # if p0 is a root of f
                print(s)  # if it is then print it's a root
                print(colored("The root is: x = " + str(p0) + ", f(x) = " + str(f(p0)) + ", found after 0 iterations",
                              'blue'))
                countRoots += 1  # add 1 to the number of roots that were found in this range
                tempRes.append(p0)  # add p0 to the results list
                flag = True  # change the flag that indicated if one or more of the initial guesses is a root to true
            break  # move on to p1
        except:  # p0 is not in f definition range
            move += 1  # update the number of tries to plus 1
            p0 += min(0.0001, epsilon)  # move p0 to the next try
    if maxTries is move:  # check if p0 is still isn't in the definition range of f
        return countRoots, tempRes  # end checking this range
    move = 0  # update the number of fixes made for a p1 to 0
    while move < maxTries:  # the max number of tries to fix p1 in case it's not in f definition range
        try:  # try to substitute p1 in f
            if f(p1) == 0:  # if p1 is a root of f
                print(s)  # if it is then print it's a root
                print(colored("The root is: x = " + str(p1) + ", f(x) = " + str(f(p1)) + ", found after 1 iterations",
                              'blue'))
                countRoots += 1  # add 1 to the number of roots that were found in this range
                tempRes.append(p1)  # add p1 to the results list
                flag = True  # change the flag that indicated if one or more of the initial guesses is a root to true
            break  # move on to p
        except:  # p1 is not in f definition range
            move += 1  # update the number of tries to plus 1
            p1 += min(0.0001, epsilon)  # move p1 to the next try
    if maxTries is move:  # check if p1 is still isn't in the definition range of f
        return countRoots, tempRes  # finish checking this range

    move = 0  # update the number of fixes made for a p2 to 0
    while move < maxTries:  # the max number of tries to fix p2 in case it's not in f definition range
        try:  # try to substitute p2 in f
            if f(p2) == 0:  # if p2 is a root of f
                print(s)  # if it is then print it's a root
                print(colored("The root is: x = " + str(p2) + ", f(x) = " + str(f(p2)) + ", found after 2 iterations",
                              'blue'))
                countRoots += 1  # add 1 to the number of roots that were found in this range
                tempRes.append(p2)  # add p2 to the results list
                flag = True  # change the flag that indicated if one or more of the initial guesses is a root to true
            break  # finish fixing p2
        except:  # p2 is not in f definition range
            move += 1  # update the number of tries to plus 1
            p2 -= min(0.0001, epsilon)  # move p2 to the next try value
    if maxTries is move:  # check if p2 is still isn't in the definition range of f
        return countRoots, tempRes  # finish checking this range

    if flag is True:  # if at least one of the initial guesses is a root of f
        return countRoots, tempRes  # finish checking this range

    # the iterations
    h1 = p1 - p0  # the distance between the x of p1 and p0
    h2 = p2 - p1  # the distance between the x of p2 and p1
    grad1 = (f(p1) - f(p0)) / h1  # The gradient between p0 and p1
    grad2 = (f(p2) - f(p1)) / h2  # The gradient between p1 and p2
    d = (grad2 - grad1) / (h2 + h1)  # calc d
    iter = 3  # Start to count the iteration from 3
    while iter <= maxIter:  # while we did not reach the maximum number of iterations
        b = grad2 + (h2 * d)  # calc b
        try:
            if (b ** 2 - (
                    4 * f(
                p2) * d)) < 0:  # check if the discriminant is negative(cant do square root on negative number)
                # there is no real root in this range
                return countRoots, tempRes  # move on to the next range
        except:  # overflow
            return countRoots, tempRes  # move on to the next range
        D = (b ** 2 - (4 * f(p2) * d)) ** 0.5  # the square root in the denominator
        if abs(b - D) < abs(b + D):  # Save the biggest result
            E = b + D
        else:
            E = b - D
        h = (-2) * f(p2) / E  # the root of the parabola created from the 3 points
        p = p2 + h  # The new guess
        try:  # check if p is in the definition range of f
            f(p)
        except:
            return countRoots, tempRes  # p is not in the definition range so move to the next range
        if abs(h) < epsilon:  # if the difference is smaller than epsilon
            if (p > start) and (p < end) and (checkResult(f, p, epsilon) is True):  # check if the point is in the range
                if abs(f(p)) <= 0 + epsilon:  # check if p is a root of f
                    print(s, end="")  # if it is then print p as a result
                    print(
                        colored("The root is: x = " + str(p) + ", f(x) = " + str(f(p)) + ", found after " + str(
                            iter) + " iterations\n", 'blue'))
                    countRoots += 1  # add 1 to the number of roots that were found in this range
                    tempRes.append(p)  # add p to the results list
                    return countRoots, tempRes  # move on to the next range
                else:
                    return countRoots, tempRes  # p is not a root of f (handles cases where the functin is constant in this range)
            else:  # there isn't a root in this range
                return countRoots, tempRes  # move to the nect range
        s += ("iteration = " + str(iter) + ", x = " + str(p) + ", f(x) = " + str(
            f(p)) + "\n")  # add the iteration info to the string
        # Preparing the variable for the next iteration
        p0 = p1  # move the points
        p1 = p2  # move the points
        p2 = p  # move the points
        h1 = p1 - p0  # the distance between p1 and p0
        h2 = p2 - p1  # the distance between p2 and p1
        grad1 = (f(p1) - f(p0)) / h1  # The gradient between p1 and p0
        grad2 = (f(p2) - f(p1)) / h2  # The gradient between p2 and p1
        d = (grad2 - grad1) / (h2 + h1)  # calc d
        iter += 1  # increase the number of iteration
    # Method failed after max number of iterations
    return countRoots, tempRes  # move to the next range


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
    results = []  # a list of the all the roots of the function in the range[]
    rootsCount = 0  # the number of roots that were found in the range
    end = start_point + 0.1  # the end point of the first range
    middle = start_point + 0.05  # the middle point of the first range
    while end <= end_point and rootsCount < maxRoots:  # while we are in the range and find the maximum roots number
        tempCount, tempRes = mullersMethod(start_point, middle, end, polinom, epsilon,
                                           iterNumber)  # calculate by muller`s method
        rootsCount += tempCount  # update the number of roots found
        results.extend(tempRes)  # add the new results to the roots list
        start_point = end  # update the start point to the next range
        middle = (start_point + 0.05)  # set the middle point of the new range
        end = start_point + 0.1  # set the end point of the new range
    print("The number of roots found is: " + str(rootsCount))
    return results


def checkDifferPart3(l, d, name1, name2, epsilon):
    """
    :param l: the roots that were received from muller's method
    :param d: the expected roots
    :param name1: the name of the first method
    :param name2: the name of the second method
    :param epsilon: maximum difference between the results
    :return: None
    """
    print("check the difference between the roots:")
    for i in range(len(l)):  # go over all the results
        print("x" + str(i + 1) + " - " + name1 + ": " + str(l[i])
              + "\nx" + str(i + 1) + " - " + name2 + ": " + str(d[i]))
        if abs(l[i] - d[i]) > epsilon:  # check if the difference between the result is bigger than epsilon
            print("The difference is bigger than epsilon for some of the components\n")
            return
    print("The difference is " + colored("smaller", 'yellow') + " than epsilon for all the components\n")


def driver():
    x = sp.symbols('x')
    # f = x ** 4 - 3 * (x ** 3) + x ** 2 + x + 1
    # f = (x ** 4 - 4 * (x ** 2)) / x
    # f = sympy.cos(x)
    epsilon = 10 ** -5
    maxRoots = 15
    start_point = -100
    end_point = 100
    maxIteration = 100
    # ----- 1st run -----
    f = (10 * x * (math.exp(1) ** x) - (x ** 3)) / (x ** 3 - 5)
    f = sp.lambdify(x, f)
    results = rangeDivision(f, epsilon, maxIteration, maxRoots, start_point, end_point)
    checkDifferPart3(results, [-1.49643, 0], "Muller's", "Expected", 10 ** -5)
    # ----- 2nd run -----
    f = log(x ** 4, 10) + x + 1
    f = sp.lambdify(x, f)
    results = rangeDivision(f, epsilon, maxIteration, maxRoots, start_point, end_point)
    checkDifferPart3(results, [-2.76978295, -1, 0.437215962], "Muller's", "Expected", 10 ** -5)
    print("program finished :)")


driver()
