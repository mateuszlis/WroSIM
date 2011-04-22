from setuptools import setup, find_packages
setup(
    name = "pymd2mc",
    version = "0.1",
    packages = find_packages(),
    
    install_requires = ['docutils>=0.3','pymacs>=0.1'],

    package_data = {
        '': ['*.txt', '*.rst'],
        'hello': ['*.msg'],
    },


    author = "Mateusz Lis",
    author_email = "mateusz.lis@pwr.wroc.pl",
    description = "This python project enables you to prepare 2 dimensional Monte Carlo simulation from 3D Molecular Dynamics simulation data. ",
    license = "PSF",
    keywords = "Academic, Python, MonteCarlo, MolecularDynamics, MD, MC, Algorithm, Simulation",
    url = "http://code.google.com/p/pymd2mc/",
    
    entry_points = {
        'console_scripts': [
            'foo = pymd2mc:main_func'
        ],
    }
  
)
