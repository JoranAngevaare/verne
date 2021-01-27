import setuptools

setuptools.setup(name='verne',
                 packages=[''],
                 package_dir={'': 'src'},
                 version='0.0.0',
                 description='Verne is a python code for calculating the Earth-stopping effect for super-heavy Dark Matter (DM).',
                 author='The Kavanagh, B. J.',
                 url='https://github.com/bradkav/verne',
                 setup_requires=['pytest-runner'],
                 python_requires=">=3.6",
                 zip_safe=False)
