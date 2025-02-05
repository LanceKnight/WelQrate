from setuptools import setup, find_packages

from setuptools import setup, find_packages

setup(
    name='welqrate', 
    version='0.1.0',  
    description='The official implementation of the WelQrate dataset and benchmark',  
    url='https://github.com/xwang38438/WelQrate.git',  
    author='Xin Wang',  
    author_email='xin.wang.1@vanderbilt.edu',  
    packages=find_packages(include=["welqrate", "welqrate.*"]),  # Only include 'welqrate' package
    install_requires=[],  
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  
    ],
    python_requires='>=3.8',  
)



# setup(
#     name='welqrate', 
#     version='0.1.0',  
#     description='The official implementation of the WelQrate dataset and benchmark',  
#     url='https://github.com/xwang38438/WelQrate.git',  
#     author='Xin Wang',  
#     author_email='xin.wang.1@vanderbilt.edu',  
#     packages=find_packages(),  
#     install_requires=[],  # List dependencies here if any #rdkit
#     classifiers=[
#         'Programming Language :: Python :: 3',
#         'License :: OSI Approved :: MIT License',  
#     ],
#     python_requires='>=3.8',  
# )
