from pertubations import UnitCellCheck, fireworks
from basin_hoping import TiO2

A = TiO2(local=True)
A.add_atoms()
fireworks(A, [12, 13], 5)

c = 0
for jj in range(0, 1000):
    state = UnitCellCheck(A.atoms)
    if state == False:
        print(state)
    else:
        c += 1

print(c)