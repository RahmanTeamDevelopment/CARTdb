from setuptools import setup

setup(
    name= 'CARTdb',
    version = '0.1.0',
    description = 'blaCreate CART database for CAVA and gff3 output',
    url = '...',
    author = 'RahmanTeam',
    author_email = 'rahmanlab@icr.ac.uk',
    license = 'MIT',
    packages=['cartdb'],
    scripts=['bin/CARTdb.py','bin/cartdb'],
    zip_safe=False
)
