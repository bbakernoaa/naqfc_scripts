try:
    from setuptools import setup, find_packages
except:
    from distutils.core import setup

setup(
    name='spatial_plot',
    version='1.0',
    author='Barry D. Baker and Patrick Cambell',
    maintainer='Barry Baker',
    packages=find_packages(),
    scripts=['spatial_plot.py'],
    license='GPL',
    long_description=open('README.md').read(),
)
