from distutils.core import setup

setup(
    name='poapangenome',
    version='0.1dev',
    packages=['poapangenome'],
    install_requires=['biopython', 'numpy', 'jsonpickle', 'ddt',
                      'https://github.com/anialisiecka/Mafgraph/archive/master.zip'],
    license='MIT Licence',
    long_description=open('readme.md').read()
)
