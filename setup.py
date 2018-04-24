from setuptools import setup

setup(  name='scdiff',
		version='1.1.13',
		description='Single Cell Differentiation Model package',
		author='Jun Ding',
		author_email='jund@andrew.cmu.edu',
		url="https://github.com/phoenixding/scdiff",
		license='MIT',
		packages=['scdiff'],
		package_data={'scdiff':['img/logo.gif','img/FlowChart.jpg']},
		entry_points={'console_scripts':['scdiff=scdiff.scdiff:main','scdiff_gui=scdiff.scdiff_gui:main']},
		install_requires=['scipy','numpy','scikit-learn'],
		classifiers=[
			'Development Status :: 3 - Alpha',
			'License :: OSI Approved :: MIT License',
			'Programming Language :: Python :: 2',
			'Programming Language :: Python :: 3',
		],
		)
