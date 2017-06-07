import numpy as np

def grid3D(min,max,sidelength):
    INT = np.linspace(min,max,sidelength)
    ARRAY = np.empty((sidelength,sidelength,sidelength,3))

    for i in range(0, sidelength):
        for j in range(0, sidelength):
            for k in range(0,sidelength):
                if i < sidelength:
                    if j < sidelength:
                        if k < sidelength:
                            ARRAY[i,j,k,:] = [INT[i], INT[j], INT[k]]
                        else:
                            pass
                    else:
                        pass
                else:
                    pass


    return ARRAY

def binning(length):
    if length < 0.5:
        index = 0
    elif length < 1.5:
        index = 1
    elif length < 2.5:
        index = 2
    elif length < 3.5:
        index = 3
    elif length < 4.5:
        index = 4
    elif length < 5.5:
        index = 5
    elif length < 6.5:
        index = 6
    elif length < 7.5:
        index = 7
    elif length < 8.5:
        index = 8
    elif length < 9.5:
        index = 9
    elif length < 10.5:
        index = 10
    elif length < 11.5:
        index = 11
    elif length < 12.5:
        index = 12
    elif length < 13.5:
        index = 13
    elif length < 14.5:
        index = 14
    elif length < 15.5:
        index = 15
    elif length < 16.5:
        index = 16
    else:
        index = 17
    return index
