from setuptools import setup
import setuptools
import pyplinkqc

# Get the long description from the README file
def load_readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='pyplinkqc',
    version=pyplinkqc.__version__,
    description='genomic data processing and analysis',
    long_description=load_readme(),
    long_description_content_type='text/markdown',
    url='https://github.ibm.com/aur-genomics/pygen.git',
    author='Nathalie Willems',
    author_email='natwille@yahoo.com',
    license='APACHE 2',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Intended Audience :: Science/Research'
    ],
    packages=setuptools.find_packages(),
    zip_safe=False,
    install_requires=[
        'sphinx >= 3.0.3',
        'sphinx_rtd_theme >= 0.4.3',
        'nose >= 1.3.7',
        'coverage >= 5.1',
        'pypi-publisher >= 0.0.4',
        'numpy>=1.18.1',
        'matplotlib>=3.1.3',
        'pandas>=1.0.3',
        'pytest>=5.3.5',
        'mypy>=0.761'
    ],
    python_requires='>=3.7'
)
