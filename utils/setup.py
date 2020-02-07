from setuptools import setup                                                                                                                                                              

setup(
    name='dtools',
    version='0.1',
    py_modules=['dtools'],
   description='A useful module',
   author='Hannes Svardal',
   author_email='svardallab@gmail.com',
    install_requires=[
        'matplotlib>=3.0.2',
       'pandas>=0.23.4',
       'jupyter>=1.0.0'
    ],
    platforms=[
    'linux-x86_64',
    'macosx-10.10-x86_64'
    ],
    include_package_data=True,
    zip_safe=False
)
