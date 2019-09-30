# Copyright 2019 Fred Hutchinson Cancer Research Center
# from distutils.core import setup

from setuptools import setup, find_packages


packages = find_packages(exclude=['tests']);


setup(name='locuszoom_plot',
      author='Keith Curtis',
      description='Plotting of genomic regions around a variant',
      license='License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
      version='0.1',
      url='https://github.com/krcurtis/locuszoom-plot',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
          'Programming Language :: Python :: 3.7'
          ],
      keywords='graph dependancy',
      packages=packages,
      install_requires=['matplotlib', 'numpy', 'pandas'],
      python_requires='>=3.7',
      )
