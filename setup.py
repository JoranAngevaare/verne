import setuptools

with open('requirements.txt') as f:
    requires = [
        r.split('/')[-1] if r.startswith('git+') else r
        for r in f.read().splitlines()]

setuptools.setup(name='verne',
                 version='0.0.0',
                 description='Verne is a python code for calculating the Earth-stopping effect for super-heavy Dark Matter (DM).',
                 author='The Kavanagh, B. J.',
                 url='https://github.com/bradkav/verne',
                 install_requires=requires,
                 packages=setuptools.find_packages() + ['data'],
                 package_dir={'verne': 'verne',
                              'data': 'data'},
                 package_data={'verne': ['data/*'],
                               },
                 setup_requires=['pytest-runner'],
                 python_requires=">=3.6",
                 zip_safe=False)
