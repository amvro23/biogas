from setuptools import setup


with open("README.md", "r", encoding="utf-8") as fh:
   long_description = fh.read()


setup(
   name='biogas',
   version='0.1.0.rc2',
   description='Biogas Reforming package.',
   long_description=long_description,
   author='Amvrosios Georgiadis',
   author_email='amvro23@gmail.com',
   packages=['biogas'],
   install_requires=[
      'numpy==1.25.2',
      'scipy>=1.13.0',
      'pandas>=2.2.1',
      'matplotlib',
      ],
)
