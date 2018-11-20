""" Fit shock breakout model to the rise

Stuff from Piro (2015)
"""

k = 0.1
Msol = 1.9E33 # grams
tp = 23/60/24 # 23 minutes between explosion time and first detection
Lp = 3E44
c = 3E10

# Arnett (1982)
Me = tp**2 * 0.2*c*c / k
# Me = 2.4E-16 solar masses = 4.6E17 grams

# in days
tp = 0.9 * k**(1/2) * E**(-1/4) * (Mc/Msol)**(0.17) * (Me/(0.01*Msol))**(0.57)
