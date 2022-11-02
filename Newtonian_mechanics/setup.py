import os
from setuptools import find_packages, setup

with open(os.path.join(os.path.dirname(__file__), 'README.md'), encoding='UTF-8') as readme:
    README = readme.read()

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(
    name='Newtonian_Mechanics',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    license='GPL License', 
    description='Newtonian Mechanics for Orbital Motion',
    long_description=README,
    url='',
    author='Tiago',
    author_email='tiagojoaog@gmail.com',
    classifiers=[
        'Intended Audience :: Students',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9', # only tested with these versions
    ],
    install_requires=[
        'keyboard',
        'matplotlib',
        'numpy',
        'vpython',
        'imageio',
        'jupyter-packaging',
        'hyperlink',
        'cryptography', 
        'txaio',
    ]
)
