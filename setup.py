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
                 packages=setuptools.find_packages(),
                 setup_requires=['pytest-runner'],
                 python_requires=">=3.6",
                 zip_safe=False)
#
# setuptools.setup(name='cutax',
#                  version='0.0.1',
#                  description='Peak matching for XENON simulations',
#                  author='The XENON collaboration',
#                  url='https://github.com/XENONnT/peakmatching',
#                  long_description=readme + '\n\n' + history,
#                  long_description_content_type="text/markdown",
#                  setup_requires=['pytest-runner'],
#                  install_requires=requires,
#                  tests_require=requires + [
#                      'pytest',
#                      'hypothesis',
#                      'boltons'],
#                  python_requires=">=3.6",
#                  packages=setuptools.find_packages(),
#                  scripts=['bin/pema_straxer'],
#                  classifiers=[
#                      'Development Status :: 4 - Beta',
#                      'License :: OSI Approved :: BSD License',
#                      'Natural Language :: English',
#                      'Programming Language :: Python :: 3.6',
#                      'Intended Audience :: Science/Research',
#                      'Programming Language :: Python :: Implementation :: CPython',
#                      'Topic :: Scientific/Engineering :: Physics',
#                  ],
#                  zip_safe=False)