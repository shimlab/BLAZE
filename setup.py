from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='blaze',
    version='2.0',
    author='Yupei You',
    author_email="youyupei@gmail.com",
    description='Barcode identification from Long reads for AnalyZing single cell gene Expression',
    packages=find_packages(),
    long_description=long_description,
    long_description_content_type='text/markdown',
    data_files=[('blaze', ['blaze/10X_bc/3M-february-2018.zip', 'blaze/10X_bc/737K-august-2016.txt'])],
    url="https://github.com/shimlab/BLAZE",
    install_requires=["fast-edit-distance==1.2.1", 
                    'matplotlib', 
                    'tqdm', 'numpy', 
                    'pandas', 
                    'biopython'],
    classifiers=[
         "Programming Language :: Python :: 3",
         "License :: GNU GENERAL PUBLIC LICENSE :: 3",
         "Operating System :: OS Independent",
     ],
     entry_points={
        'console_scripts': [
            'blaze = blaze:_pipeline',
        ],
    },
)
