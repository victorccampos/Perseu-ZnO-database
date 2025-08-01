import numpy as np


for nx in range(1,4):
  for ny in range(1,4):
    for nz in range(1,4):
      # indent
      cellsize = np.array([nx, ny, nz])
      kpoints_original = np.array([6,6,6])
      kpoints_cell = kpoints_original // cellsize
      
      print(f'| cellsize = {nx}x{ny}x{nz} | nº átomos = {4*np.prod(cellsize)} | grid-k adaptado = {kpoints_cell} | nº kpontos = {np.prod(kpoints_cell)} ' +
      f'| nelec = {2*(12+6) * np.prod(cellsize)} | nbands >= {int((1/2)*2*(12+6) * np.prod(cellsize))} |')


	    