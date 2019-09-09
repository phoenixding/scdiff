#!/usr/bin/env python
from setuptools import setup

setup(  name='scdiff',
		version='1.1.16',
		description='Single Cell Differentiation Model package',
		author='Jun Ding',
		author_email='jund@andrew.cmu.edu',
		url="https://github.com/phoenixding/scdiff",
		license='MIT',
		packages=['scdiff'],
		package_data={'scdiff':['img/logo.gif','tfdata/HumanTFList.txt']},
		entry_points={'console_scripts':['scdiff=scdiff.scdiff:main']},
		install_requires=['scipy>0.13.3','numpy>1.8.2','scikit-learn>=0.20','cloudpickle>0.5,<0.6','pydiffmap>=0.1.1,<0.2.0','matplotlib<3.0','imbalanced_learn<0.5.0'],
		classifiers=[
			'License :: OSI Approved :: MIT License',
			'Programming Language :: Python :: 2',
			'Programming Language :: Python :: 3',
		],
		)
		
