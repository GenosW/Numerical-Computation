import math as m
import time
import matplotlib.pyplot as plt
import numpy as np

def f_monomial(x, a=2, C=1):
    if a<0:
        return C/x**a
    return C*x**a

def f_exp(x, C=1, b=1):
    return C*m.exp(-b*x)

def f_exp2(x, C=1, b=1):
    return C*m.exp(-b/x)

def knotsOf2(i_max):
    return list(map(lambda x: 2**-x, list(range(0,i_max+1))))

if __name__ == "__main__":
    i_max = 10
    alphas = [0, 0.5, 1, 2]
    knots = knotsOf2(i_max)

    fis = [[f_monomial(x, a=alpha) for x in knots] for alpha in alphas]
    for fi in fis[:]:
        plt.loglog(knots, fi)
    plt.grid()
    plt.title("Monomials f(h) = C*h^alpha")
    legEntries = [str(alpha) for alpha in alphas]
    plt.legend(legEntries)
    plt.savefig("NumComp/Ue4/monomials_in_loglog")
    #plt.show()

    plt.close()
    fis = [[f_exp(x, b=alpha) for x in knots] for alpha in alphas]
    for fi in fis[:]:
        plt.semilogy(knots, fi)
    plt.grid()
    plt.title("Exponential functions f(x) = C*exp(-b*h)")
    legEntries = [str(alpha) for alpha in alphas]
    plt.legend(legEntries)
    plt.savefig("NumComp/Ue4/exp(-x)_in_loglog")
    #plt.show()

    plt.close()
    fis = [[f_exp2(x, b=alpha) for x in knots] for alpha in alphas]
    for fi in fis[:]:
        plt.loglog(knots, fi)
    plt.grid()
    plt.title("Exponential functions f(x) = C*exp(-b/h)")
    legEntries = [str(alpha) for alpha in alphas]
    plt.legend(legEntries)
    plt.savefig("NumComp/Ue4/exp(-1_over_x)_in_loglog")
    #plt.show()

    plt.close()
    fis = [[f_exp2(x, b=alpha) for x in knots] for alpha in alphas]
    for fi in fis[:]:
        plt.semilogx(knots, fi)
    plt.grid()
    plt.title("Exponential functions f(x) = C*exp(-b/h)")
    legEntries = [str(alpha) for alpha in alphas]
    plt.legend(legEntries)
    plt.savefig("NumComp/Ue4/exp(-1_over_x)_in_semilogx")
    #plt.show()

    plt.close()
    fis = [[f_exp2(x, b=alpha) for x in knots] for alpha in alphas]
    for fi in fis[:]:
        plt.semilogy(knots, fi)
    plt.grid()
    plt.title("Exponential functions f(x) = C*exp(-b/h)")
    legEntries = [str(alpha) for alpha in alphas]
    plt.legend(legEntries)
    plt.savefig("NumComp/Ue4/exp(-1_over_x)_in_semilogy")
    #plt.show()

    plt.close()
    fis = [[f_exp2(x, b=alpha) for x in knots] for alpha in alphas]
    for fi in fis[:]:
        plt.plot(knots, fi)
    plt.grid()
    plt.title("Exponential functions f(x) = C*exp(-b/h)")
    legEntries = [str(alpha) for alpha in alphas]
    plt.legend(legEntries)
    plt.savefig("NumComp/Ue4/exp(-1_over_x)_in_plot")
    #plt.show()
    