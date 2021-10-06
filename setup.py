from setuptools import setup

setup(
   name='kcet',
   version='0.3.1',
   author='Vida Ravanmehr, Peter N Robinson',
   author_email='vida.ravanmehr@jax.org,peter.robinson@jax.org',
   packages=['kcet'],
   scripts=['runRandomForest.py', 'pkpki.py'],
   license='LICENSE',
   description='Prepare protein kinase inhibitor data for ML',
   long_description=open('README.md').read(),
   install_requires=[
       "pandas",
       "nosetest",
       "numpy",
       "pandas",
       "scikit-learn"
   ],
)
