from astropy.constants import G, M_sun, L_sun, R_sun
import numpy as np
import astropy.units as u

def thermal(mass, radius, lum):
	tau=(G*(mass*M_sun)**2)/(radius*R_sun*lum*L_sun)
	return tau.to(u.s)

def dynamical(mass, radius):
	tau=np.sqrt((radius*R_sun)**3/(2*G*mass*M_sun))
	return tau

