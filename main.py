import sympy as sp


def mullersMethod(p0, p1, p2, f, epsilon, maxIter):
    """
    :param p0: The first guess
    :param p1: The second guess
    :param p2: The third guess
    :param f: The original function
    :param epsilon: The accepted error
    :param maxIter: The max number of iterations that we are willing to accept
    :return:
    """
    if f(p0) == 0:
        print("The root is: x = " + str(p0) + ", found after 1 iterations")
        return
    if f(p1) == 0:
        print("The root is: x = " + str(p1) + ", found after 2 iterations")
        return
    if f(p2) == 0:
        print("The root is: x = " + str(p2) + ", found after 3 iterations")
        return
    print("The initial guess: ")
    print("p0 = " + str(p0) + ", p1 = " + str(p1) + ", p2 = " + str(p2))
    h1 = p1 - p0
    h2 = p2 - p1
    grad1 = (f(p1) - f(p0)) / h1
    grad2 = (f(p2) - f(p1)) / h2
    d = (grad2 - grad1) / (h2 + h1)
    iter = 3
    while iter <= maxIter:
        b = grad2 + (h2 * d)
        if (b ** 2 - (4 * f(p2) * d)) < 0:
            print("There is no real root")
            return
        D = (b**2 - (4*f(p2)*d)) ** 0.5
        if abs(b - D) < abs(b + D):
            E = b + D
        else:
            E = b - D
        h = (-2) * f(p2) / E
        p = p2 + h
        if abs(h) < epsilon:
            print("The root is: x = " + str(p) + ", found after " + str(iter) + " iterations")
            return
        print("iteration = " + str(iter) + ", x = " + str(p) + ", f(x) = " + str(f(p)))
        p0 = p1
        p1 = p2
        p2 = p
        h1 = p1 - p0
        h2 = p2 - p1
        grad1 = (f(p1) - f(p0)) / h1
        grad2 = (f(p2) - f(p1)) / h2
        d = (grad2 - grad1) / (h2 + h1)
        iter += 1
    print("Method failed after " + str(maxIter) + "iterations")
    return

def driver():
    x = sp.symbols('x')
    f = x**4 - 3 * x**3 + x**2 + x + 1
    epsilon = 10**-5
    startPoint = 0.5
    middelPoint = 1
    endPoint = 1.5
    maxIteration = 100
    f = sp.lambdify(x, f)
    mullersMethod(startPoint, middelPoint, endPoint, f, epsilon, maxIteration)


driver()