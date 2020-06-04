from setuptools import setup
setup(name='comp_stress_disk',
    version='0.0.1',
    description='Compute stress of rotating annular disk.',
    author='Duy Duc NGUYEN',
    author_email='duyduc.nguyen@protonmail.com',
    packages=['comp_stress_disk'],
    install_requires=['matplotlib',
                      'numpy',
                      'pycodestyle>=2.4.0'])