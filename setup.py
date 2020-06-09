from setuptools import setup
import setuptools

# Get the long description from the README file
def load_readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='pygen',
    version='0.0.2',
    description='genomic data processing and analysis',
    url='https://github.ibm.com/aur-genomics/pygen.git',
    author='Nathalie Willems',
    author_email='nathalie.willems@ibm.com',
    license='unlicense',
    packages=setuptools.find_packages(include=['pygen', 'pygen.*']),
    zip_safe=False,
    install_requires=[
        'numpy==1.18.1',
        'matplotlib==3.1.3',
        'pandas==1.0.3',
        'pytest==5.3.5',
        'mypy==0.761'
    ],
    python_requires='>=3.7'
)
