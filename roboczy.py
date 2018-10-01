import timeit

class Test():
    def __init__(self, a: int):
        self.a = a


d = {i:i+1 for i in range(1000)}

def foo1():
    for i in range(len(d.items())):
        x = d[i]

o = [Test(i) for i in range(len(d.items()))]

def foo2():
    for i in range(len(d.items())):
        x = o[i].a

print(timeit.timeit('foo1()', globals=globals()))
print(timeit.timeit('foo2()', globals=globals()))