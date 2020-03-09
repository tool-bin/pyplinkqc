from setuptools import setup

setup(
    name='pygen',
    version='0.0.1',
    description='genomics processing and analysis',
    url='https://github.ibm.com/aur-genomics/pygen.git',
    author='Isabell Kiral',
    author_email='isa.kiral@gmail.com',
    license='unlicense',
    packages=['pygen'],
    zip_safe=False,
    install_requires=[
        'nutsflow',
        'nutsml',
	    'pandas-plink',
        'numpy',
        'pandas',
        'pytest',
        'mypy'
    ],
)
