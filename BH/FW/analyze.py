import numpy as np
import glob

files = glob.glob('npys/*.npy')



cml = 0; cnml = 0
tml = 0; tnml = 0
for file in files:
    data = np.load(file)
    if file[5:7] == 'ML':
        tml += 1
        if data[0] == 1:
            cml += 1
    else:
        tnml += 1
        if data[0] == 1:
            cnml += 1

print('With autobag:')
print('Runs that finished: {}'.format(tml))
print('Runs that converged: {}'.format(cml))

print('Without autobag:')
print('Runs that finished: {}'.format(tnml))
print('Runs that converged: {}'.format(cnml))

