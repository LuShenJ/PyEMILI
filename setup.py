from setuptools import setup, find_packages
import os
import zipfile


with open('requirements.txt') as f:
    requirements = f.read().splitlines()

linedb_path = os.path.join(os.path.dirname(__file__), 'pyemili','Line_dataset')
linedb = os.path.join(os.path.dirname(__file__), 'pyemili','Line_dataset','Linedb.zip')
with zipfile.ZipFile(linedb, 'r') as zip_ref:
    zip_ref.extractall(linedb_path)


setup(name='pyemili', 
      version='1.0.0',
      description='Spectral lines identifier',
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
      author='Zhijun Tu, Xuan Fang, Robert Williams',
      author_email='zjtu@bao.ac.cn',
      url='https://github.com/LuShenJ/PyEMILI',
      install_requires=requirements,
      packages=find_packages(include=['pyemili']),
      package_data={'pyemili':['abundance/*.dat',
                        'Line_dataset/*.dat',
                        'Line_dataset/*.npy',
                        'recom/*.dat',
                        'eff_reccoe/*.dat',
                        'eff_reccoe/*.npy',]
                    },
)