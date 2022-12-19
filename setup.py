from setuptools import setup


with open("README.md", "r", encoding="utf-8") as fh:
   long_description = fh.read()


setup(
   name='Biogas_reforming',
   version='0.1.0.rc2',
   description='Biogas Reforming package.',
   long_description=long_description,
   author='Amvrosios Georgiadis',
   author_email='amvro23@gmail.com',
   packages=['Biogas_reforming'],
   install_requires=[
      'numpy==1.20.*',
      'scipy>=1.7.*',
      'pandas>=1.1.*',
      'matplotlib',
      ],
)
