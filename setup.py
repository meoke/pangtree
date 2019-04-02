from distutils.core import setup

setup(
    name='pang',
    version='0.1dev',
    packages=['biopython', 'numpy', 'jsonpickle', 'ddt', 'mafgraph'],
    license='MIT Licence',
    long_description=open('README.md').read(),
)