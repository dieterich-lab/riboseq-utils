from setuptools import setup, find_packages

console_scripts = [
    'extract-metagene-profiles=riboutils.extract_metagene_profiles:main',
    'estimate-metagene-profile-bayes-factors=riboutils.estimate_metagene_profile_bayes_factors:main',
    'select-periodic-offsets=riboutils.select_periodic_offsets:main',
    'bootstrap-ribo-analysis=riboutils.bootstrap_ribo_analysis:main',
    'pickle-stan=riboutils.pickle_stan:main'
]

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='riboutils',
        version='0.2.7',
        description="This package contains utilities for other ribosome profiling projects.",
        long_description=readme(),
        keywords="ribosome profiling utilities translation",
        url="",
        author="Brandon Malone",
        author_email="bmmalone@gmail.com",
        license='MIT',
        packages=find_packages(),
        install_requires = [
            'numpy',
            'pandas',
            'scipy',
            'tqdm',
            'appdirs',
            'statsmodels',
            'pysam',
            'pyyaml',
            'misc==0.2.8',
            'bio-utils==0.2.6',
            'pystan==2.16.0.0'
        ],
        extras_require = {
        },
        include_package_data=True,
        test_suite='nose.collector',
        tests_require=['nose'],
        entry_points = {
            'console_scripts': console_scripts
        },
        zip_safe=False
)
