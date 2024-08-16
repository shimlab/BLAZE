from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='blaze2',
    version='2.4.0',
    author='Yupei You',
    author_email="youyupei@gmail.com",
    description='Barcode identification from Long reads for AnalyZing single cell gene Expression',
    packages=find_packages(),
    test_suite='test',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={'blaze': ['10X_bc/3M-february-2018.zip', '10X_bc/737K-august-2016.txt', '10X_bc/3M-5pgex-jan-2023.zip', "10X_bc/3M-3pgex-may-2023.zip"]},
    url="https://github.com/shimlab/BLAZE",
    install_requires=["fast-edit-distance==1.2.1", 
                    'matplotlib', 
                    'tqdm', 'numpy', 
                    'pandas'],
    classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
         "Operating System :: OS Independent",
     ],
     entry_points={
        'console_scripts': [
            'blaze = blaze:_pipeline',
        ],
    },
)
