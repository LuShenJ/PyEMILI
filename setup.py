version = '1.0.1'

from setuptools import setup, find_packages
import os


with open('requirements.txt') as f:
    requirements = f.read().splitlines()


setup(name='pyemili', 
      version=version,
      description='Spectral lines identifier',
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.rst'),encoding='utf-8').read(),
      author='Zhijun Tu, Xuan Fang, Robert Williams, Jifeng Liu',
      author_email='zjtu@bao.ac.cn',
      url='https://github.com/LuShenJ/PyEMILI',
      install_requires=requirements,
      packages=find_packages(include=['pyemili']),
      package_data={'pyemili':['abundance/*.dat',
                        'Line_dataset/*.dat',
                        'Line_dataset/*.npy',
                        'Line_dataset/*.npz',
                        'recom/*.dat',
                        'eff_reccoe/*.dat',
                        'eff_reccoe/*.npy',]
                    },
)