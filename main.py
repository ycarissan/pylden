#!/usr/bin/python3

import sys
import numpy as np
import logging
import pyvista as pv
import argparse

import geometry.geometry
import view.view
import pyvista

def main():
    p = pyvista.Plotter()
    geom = geometry.geometry.Geometry("geom.xyz")
    geom.print()
    view.view.addGeometry_to_view(p,geom)
    color=[1,0,0]
    vect_magn=np.array([0.001977, -0.000914, 0.006200])
    view.view.addVector_to_view(p, vect_magn, color)
    color=[0,0,1]
    vect_elec=np.array([0.014972, -0.016617, -0.675772])
    view.view.addVector_to_view(p, vect_elec, color)
    p.show()

main()
