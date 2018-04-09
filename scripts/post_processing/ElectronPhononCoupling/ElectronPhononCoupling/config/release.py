# coding: utf-8

__all__ = [
    '__version__',
    'name',
    'description',
    'long_description',
    'license',
    '__author__',
    'author',
    'author_email',
    'url',
    ]


__version__ = '4.6.0'

name = "ElectronPhononCoupling"

description = "Python module to analyze electron-phonon related quantities."

long_description = """"
    Compute electron-phonon coupling related quantities, such as:
        - zero-point renormalization
        - temperature dependance of eigenvalues
        - quasiparticle lifetime from the el-ph self-energy
        - frequency-dependent self-energy
        - spectral function
    """

license = 'GPL'

authors = {
    'GA': ('Gabriel Antonius', 'gabriel.antonius at gmail.com'),
    }
        
author = 'The ABINIT group'

author_email = authors['GA'][1]

url = 'https://github.com/GkAntonius/ElectronPhononCoupling'

__author__ = ''
for auth, email in authors.itervalues():
  __author__ += auth + ' <' + email + '>\n'
del auth, email


