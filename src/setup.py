#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""PyRAMSES module."""

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages
import os
import re

def read_first_existing(*paths):
    base = os.path.dirname(__file__)
    for path in paths:
        full_path = os.path.join(base, path)
        if os.path.exists(full_path):
            with open(full_path, encoding='utf-8') as f:
                return f.read()
    raise FileNotFoundError("None of the candidate files exist: {}".format(paths))


# Metadata is parsed (not imported) from the package so that setup.py works
# in pip's isolated build environment, where the package's runtime
# dependencies (numpy, scipy, ...) are not installed.
def read_metadata():
    init_path = os.path.join(os.path.dirname(__file__), 'pyramses', '__init__.py')
    with open(init_path, encoding='utf-8') as f:
        content = f.read()
    def grab(name):
        match = re.search(r"^__%s__\s*=\s*['\"]([^'\"]+)['\"]" % name, content, re.M)
        if not match:
            raise RuntimeError("__%s__ not found in %s" % (name, init_path))
        return match.group(1)
    return {name: grab(name) for name in
            ('version', 'author', 'email', 'status', 'url', 'package_name')}

install_requires = ['matplotlib','scipy','numpy>=1.20.0','mkl==2024.2.2','pandapower[all]==3.3.0']

_meta = read_metadata()
__version__ = _meta['version']
__author__ = _meta['author']
__email__ = _meta['email']
__status__ = _meta['status']
__url__ = _meta['url']
__name__ = _meta['package_name']

setup(
    name=__name__,
    version=__version__,
    description='Custom Python library for combining RAMSES dynamic simulator and Pandapower.',
    author=__author__,
    author_email=__email__,
    url=__url__,
    keywords=['RAMSES', 'Power Systems', 'Simulator'],
    license='Apache-2.0',
    long_description=read_first_existing('../README.rst', 'README.rst'),
    long_description_content_type='text/x-rst',
    packages=find_packages(),
    install_requires=install_requires, 
    python_requires=">=3.8",
    package_data={
        'pyramses': ['libs/*.dll', 'libs/*.so', 'libs/*.dylib', 'libs/*.h'],
    },
    classifiers=[
        "Development Status :: " + __status__,
        "Intended Audience :: Developers",
        "Environment :: Console",
        #"License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3"
    ],
    entry_points={
        'console_scripts' : [
            'pystepss-ulg = pyramses.scripts.exec:run',
        ]
    }

)
