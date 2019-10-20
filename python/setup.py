import os
from setuptools import setup, Extension

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

setup(
    name="cmri",
    version='0.0.1',
    author="O Addy",
)
