#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name='az',
      version='0.1',
      author="Rory Kirchner",
      description="Plugins for AZ specific sequence analysis",
      url='http://github.com/roryk/bipy',
      author_email='rory.kirchner@gmail.com',
      license='MIT',
      namespace_packages=['az'],
      packages=find_packages(),
      dependency_links=['https://github.com/roryk/bipy/tarball/master#egg=bipy-0.1.0'],
      install_requires=["bipy == 0.1.0"],
      zip_safe=False)
