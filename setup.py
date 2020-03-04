from setuptools import setup

setup(
    name='nutsgen',
    version='0.0.1',
    description='genomics processing and analysis',
    url='https://github.ibm.com/isabeki/pylinklearn.git',
    author='Isabell Kiral',
    author_email='isa.kiral@gmail.com',
    license='unlicense',
    packages=['nutsgen'],
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
