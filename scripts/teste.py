tamanhos_de_celula = []

for nx in range(1,4):
  for ny in range(1,4):
     for nz in range(1,4):
       tamanhos_de_celula.append((nx,ny,nz)) 
       print(f'{(nx,ny,nz)} - Nº atomos = {4*nx*ny*nz}')



print(tamanhos_de_celula)
