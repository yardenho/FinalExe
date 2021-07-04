import sympy as sp
from termcolor import colored, cprint


# def mullersMethod(p0, p1, p2, f, epsilon, maxIter):
#     """
#     :param p0: The first guess
#     :param p1: The second guess
#     :param p2: The third guess
#     :param f: The original function
#     :param epsilon: The accepted error
#     :param maxIter: The max number of iterations that we are willing to accept
#     :return: None
#     """
#     # Check if one of the initial guess is the root
#     if f(p0) == 0:
#         print("The root is: x = " + str(p0) + ", found after 1 iterations")
#         return
#     if f(p1) == 0:
#         print("The root is: x = " + str(p1) + ", found after 2 iterations")
#         return
#     if f(p2) == 0:
#         print("The root is: x = " + str(p2) + ", found after 3 iterations")
#         return
#
#     print("The initial guess: ")
#     print("p0 = " + str(p0) + ", p1 = " + str(p1) + ", p2 = " + str(p2))
#
#     h1 = p1 - p0  #
#     h2 = p2 - p1  #
#     grad1 = (f(p1) - f(p0)) / h1  # The gradient between p0 and p1
#     grad2 = (f(p2) - f(p1)) / h2  # The gradient between p1 and p2
#     d = (grad2 - grad1) / (h2 + h1)
#     iter = 3  # Start to count the iteration
#     while iter <= maxIter:  # while we did not reach the maximum iterations
#         b = grad2 + (h2 * d)
#         if (b ** 2 - (4 * f(p2) * d)) < 0:  # if the discriminant is negative - cant do square root on negative number
#             print("There is no real root")
#             return
#         D = (b**2 - (4*f(p2)*d)) ** 0.5  # the square root in the denominator
#         if abs(b - D) < abs(b + D):  # Save the biggest result
#             E = b + D
#         else:
#             E = b - D
#         h = (-2) * f(p2) / E
#         p = p2 + h
#         if abs(h) < epsilon:
#             print("The root is: x = " + str(p) + ", found after " + str(iter) + " iterations")
#             return
#         print("iteration = " + str(iter) + ", x = " + str(p) + ", f(x) = " + str(f(p)))
#         p0 = p1
#         p1 = p2
#         p2 = p
#         h1 = p1 - p0
#         h2 = p2 - p1
#         grad1 = (f(p1) - f(p0)) / h1
#         grad2 = (f(p2) - f(p1)) / h2
#         d = (grad2 - grad1) / (h2 + h1)
#         iter += 1
#     print("Method failed after " + str(maxIter) + "iterations")
#     return

def mullersMethod(p0, p1, p2, f, epsilon, maxIter):
    """
    :param p0: The first guess
    :param p1: The second guess
    :param p2: The third guess
    :param f: The original function
    :param epsilon: The accepted error
    :param maxIter: The max number of iterations that we are willing to accept
    :return: None
    """
    start = p0
    end = p2
    s = "The initial guess:\np0 = " + str(p0) + ", p1 = " + str(p1) + ", p2 = " + str(p2) + "\n"
    # Check if one of the initial guess is the root
    flag = False
    if f(p0) == 0:
        print(s)
        print("The root is: x = " + str(p0) + ", found after 1 iterations")
        flag = True
    if f(p1) == 0:
        print(s)
        print("The root is: x = " + str(p1) + ", found after 2 iterations")
        flag = True
    if f(p2) == 0:
        print(s)
        print("The root is: x = " + str(p2) + ", found after 3 iterations")
        flag = True
    if flag is True:
        return

    h1 = p1 - p0
    h2 = p2 - p1
    grad1 = (f(p1) - f(p0)) / h1  # The gradient between p0 and p1
    grad2 = (f(p2) - f(p1)) / h2  # The gradient between p1 and p2
    d = (grad2 - grad1) / (h2 + h1)
    iter = 3  # Start to count the iteration
    while iter <= maxIter:  # while we did not reach the maximum iterations
        b = grad2 + (h2 * d)
        if (b ** 2 - (4 * f(p2) * d)) < 0:  # if the discriminant is negative - cant do square root on negative number
            print("*** There is no real root ***")
            return
        D = (b**2 - (4*f(p2)*d)) ** 0.5  # the square root in the denominator
        if abs(b - D) < abs(b + D):  # Save the biggest result
            E = b + D
        else:
            E = b - D
        h = (-2) * f(p2) / E
        p = p2 + h  # The new guess
        if abs(h) < epsilon:  # if the difference is smaller then epsilon
            if (p > start) and (p < end):  # if the point is in the range
                print(s, end="")
                print(colored("The root is: x = " + str(p) + ", found after " + str(iter) + " iterations\n", 'blue'))
                return
            else:
                print("*** There is no root in this range ***")
                return
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
    print("Method failed after " + str(maxIter) + "iterations")
    return


def rangeDivision(polinom, start_point, end_point, epsilon, iterNumber):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The accepted error
    :param iterNumber: The max number of iterations that we are willing to accept
    :return: None
    """
    # dividing the range to small parts
    end = start_point + 0.1
    middel = start_point + 0.05
    while end <= end_point:  # while we dont reach to the end point of the range
        print("range: [" + str(start_point) + ", " + str(end) + "]")
        mullersMethod(start_point, middel, end, polinom, epsilon, iterNumber)  # calculate by muller`s method
        start_point = end  # increase the start point of the range
        middel = (start_point + 0.05)  # set the middle point of the range
        end = start_point + 0.1  # set the end point of the range

def driver():
    x = sp.symbols('x')
    f = x**4 - 3 * (x ** 3) + x**2 + x + 1
    epsilon = 10**-5
    startPoint = 0.5
    endPoint = 3
    maxIteration = 100
    f = sp.lambdify(x, f)
    rangeDivision(f, startPoint, endPoint, epsilon, maxIteration)

driver()

