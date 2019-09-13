import setuptools


setuptools.setup(
    name='pangtreebuild',
    version='0.3dev',
    author="Paulina Dziakdiewicz",
    author_email="pedziadkiewicz@gmail.com",
    description="Multiple sequence alignment analysis with Consensus Tree generation",
    url="https://github.com/meoke/pangtree",
    packages=setuptools.find_packages(exclude=("tests",)),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    install_requires=['biopython', 'numpy', 'jsonpickle', 'ddt', 'networkx', 'six', 'newick'],
    license='MIT Licence',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',

    include_package_data=True
)
