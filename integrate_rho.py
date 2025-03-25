#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2023-10-20 at 16:51:05 CEST by David Gaspard <david.gaspard@espci.fr>
## Updated on 2024-02-01 at 17:57:38 CET for the Ebsolve program written in Fortran.
## Python script to integrate the density rho(T) and check if it is normalized to 1.
import sys, os
import numpy as np

def check_meshes(tmesh, rhomesh):
	"""
	Check the given meshes before further computations.
	"""
	## Check the lengths:
	if (tmesh.shape != rhomesh.shape):
		raise ValueError("Invalid meshes. Array lengths do not match.")
	## Check the T mesh:
	if (np.any(tmesh <= 0.)):
		raise ValueError("Invalid T mesh. T must be positive.")
	if (np.any(tmesh >= 1.)):
		raise ValueError("Invalid T mesh. T must be smaller than 1.")
	if (np.any(np.diff(tmesh) < 0.)):
		raise ValueError("Invalid mesh of T. The mesh is not in increasing order.")
	## Check the rho(T) mesh:
	negs = np.sum(np.any(rhomesh < 0.)) ## Number of negative elements.
	if (negs != 0):
		print(f"[INFO] Some values of rho(T) are negative ({negs} values). Absolute values will be assumed.")
	return

def check_rhomesh(rhomesh):
	"""
	Check if the tmesh contains abscissas in increasing order.
	"""
	if (np.any(tmesh <= 0.)):
		raise ValueError("Invalid T mesh. T must be positive.")
	if (np.any(tmesh >= 1.)):
		raise ValueError("Invalid T mesh. T must be smaller than 1.")
	if (np.any(np.diff(tmesh) < 0.)):
		raise ValueError("Invalid mesh of T. The mesh is not in increasing order.")
	return

def rho_integrate_midpoint(tmesh, rhomesh):
	"""
	Computes the integral rho(T) dT over T=[0, 1] represented by the arrays of abscissas 'tmesh' and ordinates 'rhomesh'. Note that the 'tmesh' is generally irregular.
	This function uses midpoint integration.
	"""
	ntm = tmesh.shape[0]
	ti = 0.  ## Initial point over the current interval.
	tf = 0.  ## Final point over the current interval.
	rhoint = 0.
	
	for i in range(ntm):
		ti = tf
		if (i+1 < ntm):
			tf = (tmesh[i] + tmesh[i+1])/2
		else:
			tf = 1.
		rhoint += abs(rhomesh[i]) * (tf-ti)  ## Integrate rho(T) assuming the absolute value.
	
	return rhoint


def rho_mean_midpoint(tmesh, rhomesh):
	"""
	Computes the mean T rho(T) dT over T=[0, 1] represented by the arrays of abscissas 'tmesh' and ordinates 'rhomesh'. Note that the 'tmesh' is generally irregular.
	This function uses midpoint integration.
	"""
	ntm = tmesh.shape[0]
	ti = 0.  ## Initial point over the current interval.
	tf = 0.  ## Final point over the current interval.
	rhoint = 0.
	
	for i in range(ntm):
		ti = tf
		if (i+1 < ntm):
			tf = (tmesh[i] + tmesh[i+1])/2
		else:
			tf = 1.
		rhoint += (ti+tf)/2 * abs(rhomesh[i]) * (tf-ti)  ## Integrate rho(T) assuming the absolute value.
	
	return rhoint

def main(args):
	for fname in args:
		print(f"[INFO] Opening file '{fname}'.")
		if (os.path.isfile(fname)):
			try:
				tmesh, rhomesh, nitermesh = np.loadtxt(fname, delimiter=",", dtype=float, comments=["%", "tm"]).transpose()
			except Exception as exc:
				print(f"[ERROR] Numpy returned exception '{exc}'.") 
				print(f"[ERROR] Input '{fname}' may not be a CSV file.")
				return 1
			#print("[INFO] tmesh =\n", tmesh)
			#print("[INFO] rhomesh =\n", rhomesh)
			check_meshes(tmesh, rhomesh)
			rhoint  = rho_integrate_midpoint(tmesh, rhomesh)
			rhomean = rho_mean_midpoint(tmesh, rhomesh)
			print(f"[INFO] Total niter  = {np.sum(nitermesh)}.")
			print(f"[INFO] rho(T) dT    = {rhoint}.")
			#print(f"[INFO] Relative error     = {abs(rhoint-1)}.")
			print(f"[INFO] T rho(T) dT  = {rhomean}.")
			print(f"[INFO] <T>          = {rhomean/rhoint}.")
		else:
			ValueError(f"Invalid file '{fname}'.")
	return 0

if (__name__ == '__main__'):
	exit(main(sys.argv[1:]))
