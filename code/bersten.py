""" Exercise assigned to me by Shri from Bersten (2018) """

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


def load_model():
    dat = np.loadtxt("../data/bersten.txt")
    dt = dat[:,0]
    absmag = dat[:,1]
    order = np.argsort(dt)
    dt = dt[order]
    absmag = absmag[order]
    print(dt[0], dt[-1])
    model = interpolate.interp1d(dt, absmag)
    return model


def plot_model(x):
    model = load_model()
    plt.plot(x, model(x), c='k', lw=2)
    plt.gca().invert_yaxis()
    plt.xscale('log')
    plt.xlabel("Rest-frame days after CC", fontsize=14)
    plt.ylabel("Absolute V mag", fontsize=14)

    plt.show()


def obs_1d():
    """ observe with a 1-day cadence """
    x = np.arange(0,10)


if __name__=="__main__":
    x = np.logspace(np.log10(0.01125), np.log10(536.3), 100000)
    plot_model(x)
