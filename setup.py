from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='riboutils',
        version='0.1',
        description="This package contains utilities for other ribosome profiling projects.",
        long_description=readme(),
        keywords="ribosome profiling utilities translation",
        url="",
        author="Brandon Malone",
        author_email="bmmalone@gmail.com",
        license='MIT',
        packages=['riboutils'],
        install_requires = [
            'numpy',
            'pandas',
            'scipy',
            'tqdm',
            'misc[bio]'
        ],
        extras_require = {
        },
        include_package_data=True,
        test_suite='nose.collector',
        tests_require=['nose'],
        entry_points = {
            'console_scripts': [
                                
                               ]
        },
        zip_safe=False
        )
