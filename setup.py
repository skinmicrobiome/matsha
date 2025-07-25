from setuptools import setup, find_packages

setup(
    name='matsha',
    version='0.1.0',
    author='Zeyang Shen',
    author_email='zeyang.shen@wsu.edu',
    description='MATSHA: a microbial GWAS tool for trait association from shotgun metagenomic data',
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    url='https://github.com/skinmicrobiome/matsha',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[],
    entry_points={
        'console_scripts': [
            'matsha=matsha.cli:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.12',
)