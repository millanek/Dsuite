from setuptools import setup                                                                                                                                                              

setup(
    name='dtools',
    version='0.1',
    py_modules=['dtools'],
    install_requires=[
        'matplotlib>=3.0.2',
        'pandas>=0.23.4'
    ],
    include_package_data=True,
    zip_safe=False
)
