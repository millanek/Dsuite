from setuptools import setup                                                                                                                                                              

setup(
    name='dtools',
    version='0.1',
    py_modules=['dtools'],
   description='A python module for plotting fbranch',
   author='Hannes Svardal',
   author_email='svardallab@gmail.com',
   scripts=['dtools.py'],
    install_requires=[
        'matplotlib>=3.0.2',
       'pandas>=0.23.4'],
    platforms=[
    'linux-x86_64',
    'macosx-10.10-x86_64'
    ],
    include_package_data=True,
    zip_safe=False
)
