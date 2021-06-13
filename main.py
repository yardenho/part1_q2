"""https://github.com/yardenho/part1_q2.git"""

import sympy as sp
from sympy.utilities.lambdify import lambdify
from numpy import log as ln
import math


def calcDerivative(func):
    """
    :param func: original function
    :return: None
    """
    x = sp.symbols('x')
    f_prime = func.diff(x)
    return f_prime


def simpson(f, startPoint, endPoint, parts):
    """
    :param f: original function
    :param startPoint: start of range
    :param endPoint: end of range
    :param parts: amount of segments
    :return: approximate area of the integral
    """
    if parts % 2 == 1:  # if there is not even numbers of parts
        print("Amount of parts must be even")
        return None
    x = sp.symbols('x')
    func = lambdify(x, f)
    gap = abs(endPoint - startPoint) / parts  # calculate h
    string = "Integral(" + str(startPoint) + ", " + str(endPoint) + ") = 1/3 * " + str(gap) + "[f(" + str(startPoint) + ")"
    appr = func(startPoint)  # placing the start point in the function
    for i in range(1, parts):  # run over the parts
        if i % 2 == 0:  # if is the even place
            string += " + 2 * f(" + str((i * gap) + startPoint) + ")"
            appr += 2 * func((i * gap) + startPoint)
        else:  # if is not the even place
            string += " + 4 * f(" + str((i * gap) + startPoint) + ")"
            appr += 4 * func((i * gap) + startPoint)
        if i % 4 ==0:  # for the printing
            string += "\n"
    string += " * f(" + str(endPoint) + ")]\n"
    print(string)  # print the equation
    appr += func(endPoint)
    appr *= 1 / 3 * gap
    return appr


def rombergMethod(f, a, b, end, epsilon):
    """
    :param f: Original function
    :param a: start of the range
    :param b: end of the range
    :param end: limit of iteration
    :param epsilon: allowed error
    :return: The area in the range
    """
    results = [[0 for i in range(end + 1)] for j in range(end + 1)]  # build matrix
    for k in range(0, end):
        print("R" + str(k+1) + "," + str(1) + " = ", end="")
        res = trapezoidMethod(f, a, b, 2 ** k)  # calculate the values of trapezoid method
        results[k+1][1] = res  # save the value in the matrix
        print(" = " +str(res))  # print the value
    for j in range(2, end + 1):
        for k in range(j, end + 1):
            results[k][j] = results[k][j - 1] + ((1 / ((4 ** (j - 1)) - 1)) * (results[k][j - 1] - results[k - 1][j - 1]))
            print("R" + str(k) + "," + str(j) + " = " + str(results[k][j]))  # print the value
            if abs(results[k][j] - results[k - 1][j]) <= epsilon:  # if the difference is less then epsilon
                return results[k][j]
    return results[k][j]


def trapezoidMethod(f, a, b, n):
    """
    :param f: Original function
    :param a: start of the range
    :param b: end of the range
    :param n: the number of the segments
    :return: The area in the range
    """
    x = sp.symbols('x')
    f = lambdify(x, f)
    h = (b - a) / n
    sum = 0
    save = a
    count = 0
    while a < b:
        sum += 0.5 * ((a + h) - a) * (f(a) + f(a + h))
        count +=1
        if a is not save:
            print(" + ", end="")
        if count is 3:
            print("\n       ", end="")
            count = 0
        print("1/2 * (" + str(b) + " - " + str(a) + ") * (f(" + str(a) + " + f(" + str(a + h) + "))", end="")
        a += h
    return sum


def rangeDivision(polinom, start_point, end_point, epsilon, function):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :param function: The function that needs to be activated
    :return: None
    """
    results = []
    sPoint = start_point
    ePoint = end_point
    flag = False
    temp = start_point + 0.1
    while temp <= end_point:  # while we dont reach to the end point of the range
        print("Range: [" + str(start_point) + ", " + str(temp) + "]")
        res, iter = function(polinom, start_point, temp,
                             epsilon)  # activates the requested function with the original function
        if iter is not None:  # if the return iteration value is not None
            if (res > -epsilon) and (res < epsilon):  # check if the result is very close to 0
                res = 0
            print("The root is " + calcFinalResult(str(res), 10**-4, '13', '18', '41') + "\nNumber of iteration is: " + str(iter))
            results.append(res)
            flag = True
        else:
            print("* There is no Intersection root in range ")
        der = calcDerivative(polinom)  # calculate the derivative
        x = sp.symbols('x')
        res, iter = function(der, start_point, temp,
                             epsilon)  # activates the requested function with the derivative function.
        if iter is not None:  # if the return iteration value is not None
            if (res > -epsilon) and (res < epsilon):  # check if the result is very close to 0
                res = 0
            if (lambdify(x, polinom)(res) > - epsilon) and (
                    lambdify(x, polinom)(res) < epsilon):  # only if the res is root
                print("The root is " + calcFinalResult(str(res), 10**-4, '13', '18', '41') + "\nNumber of iteration is: " + str(iter))
                results.append(res)
                flag = True
            else:
                print("This point is a root of the derivative but is not a root of the function..\n")
        else:
            print("* There is no touching root in range\n")
        start_point = temp  # increase the start point of the range
        temp += 0.1  # increase the end point of the range
    if flag is False:
        print("There is no root in the range " + "[" + str(sPoint) + ", " + str(ePoint) + "]")
    return results


def NewtonRaphson(polinom, start_point, end_point, epsilon):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :return: None
    """
    return rangeDivision(polinom, start_point, end_point, epsilon, calcByNewtonRaphson)


def calcByNewtonRaphson(pol, startPoint, endPoint, epsilon):
    """
    :param pol: Original function
    :param startPoint: int value, the start point of the range
    :param endPoint: int value, the end point of the range
    :param epsilon: The excepted error
    :return: The result and the number of the iteration for getting that result
    """
    Xr = (startPoint + endPoint) / 2  # middle of the range
    iteration = 1
    der = calcDerivative(pol)  # calculate the derivative
    x = sp.symbols('x')
    pol = lambdify(x, pol)
    der = lambdify(x, der)
    if pol(startPoint) * pol(endPoint) > 0:  # check if the is change in the sign in the function
        return None, None
    print("==== Iterations ====")
    while iteration < 100:
        res1 = pol(Xr)
        res2 = der(Xr)
        print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(res1) + ", f`(x) = " + str(res2))
        iteration += 1
        Xnext = Xr - (res1 / res2)
        if abs(Xnext - Xr) < epsilon:
            print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(
                res1) + ", f`(x) = " + str(res2))
            print("==== End of iterations ====\n")
            return Xnext, iteration
        Xr = Xnext
    print("The system does not converge... :(")
    return None, None


def secant_method(polinom, start_point, end_point, epsilon):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :return: None
    """
    return rangeDivision(polinom, start_point, end_point, epsilon, calcBySecant)


def calcBySecant(polinom, start_point, end_point, epsilon):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :return: The result and the number of the iteration for getting that result
    """
    Xr = start_point
    Xnext = end_point
    iteration = 1
    x = sp.symbols('x')
    polinom = lambdify(x, polinom)
    if polinom(start_point) * polinom(end_point) > 0:  # check if the is change in the sign in the function
        return None, None
    print("==== Iterations ====")
    while iteration < 100:
        res1 = polinom(Xr)
        res2 = polinom(Xnext)
        print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(res1) + ", f`(x) = " + str(res2))
        iteration += 1
        temp = Xnext
        Xnext = ((Xr * res2) - (Xnext * res1)) / (res2 - res1)
        Xr = temp
        if abs(Xnext - Xr) < epsilon:
            print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(
                res1) + ", f`(x) = " + str(res2))
            print("==== End of iterations ====\n")
            return Xnext, iteration
    print("The system does not converge... :(")
    return None, None


def checkDiffer(l, d, epsilon):
    print("**** check the difference between the methods: ****")
    flag = True
    for _ in range(len(l)):
        print("Root " + str(_) + ":\nSecant: " + str(l[_]) + ", Newton Raphson: " + str(d[_]))
        if abs(l[_] - d[_]) > epsilon:
            flag = False
            print("The difference is bigger than the epsilon for some of the roots")
            return
    print("The difference is smaller than the epsilon for all the roots")


def calcFinalResult(result, epsilon, day, hour, minutes):
    """
    :param result: the result
    :param epsilon: the epsilon we decided on for the question
    :param day: the day of the calculation
    :param hour: the hour of the calculation
    :param minutes: the minutes of the calculation
    :return: the result by the requested format
    """
    stringRes = str(result)  # cast the result to string
    i = 0
    while stringRes[i] is not ".":  # run over the string while we get to the point
        i += 1  # count how many digits there is before the point
    i += 1
    count = 1
    while epsilon < 1:  # checking how digit needs after the point
        epsilon *= 10
        count += 1
    stringRes = stringRes[:i + count] + "00000" + day + hour + minutes
    return stringRes



def driver():
    print("******* Part 1 - Q2_a *******")
    x = sp.symbols('x')
    f = (sp.cos((x ** 2) + (5 * x) + 6) / (2 * ((math.exp(1)) ** (-x))))
    print("-----------Secant Method -----------")
    l = secant_method(f, -1.5, 2, 0.0001)
    print("-----------Newton Raphson -----------")
    d = NewtonRaphson(f, -1.5, 2, 0.0001)
    checkDiffer(l, d, 10 ** -7)
    print("\n******* Part 1 - Q2_b *******")

    start_point = 0
    end_point = 1
    part = 18
    epsilon = 10 ** -4
    print("-----Simpson Method -----")
    s = simpson(f, start_point, end_point, part)
    print("Final result:\nIntegral(" + str(start_point) + ", " + str(end_point) + ") = " + str(
        calcFinalResult(s, epsilon, '13', '18', '41')))
    print("\n-----Romberg Method -----")
    r = rombergMethod(f, start_point, end_point, 5, epsilon)
    print("\nFinal result:\nIntegral(" + str(start_point) + ", " + str(end_point) + ") = " + str(
        calcFinalResult(r, epsilon, '13', '18', '41')))

    if abs(s - r) <= epsilon:
        print("\n* The difference between the two methods is smaller than the epsilon")
    else:
        print("\n* The difference between the two methods is bigger than the epsilon - needs another method")


driver()



