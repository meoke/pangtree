def test():
    import numpy as np
    A = np.array([i for i in range(100)])
    B = [34,70,100]
    c = np.intersect1d(A, B)

def test2():
    import numpy as np
    A = np.array([i for i in range(100)])
    B = [34,70,100]
    c = []
    for b in B:
        if b in A:
            c.append(b)

def test3():
    import numpy as np
    def index(a, x):
        i = np.searchsorted(a, x)
        if i != len(a) and a[i] == x:
            return i
        raise ValueError

    A = np.array([i for i in range(100)])
    B = [34, 70, 100]
    c = []
    for b in B:
        try:
            index(A, b)
            c.append(b)
        except:
            pass

if __name__ == '__main__':
    import timeit
    print(timeit.timeit("test3()", setup="from __main__ import test3"))