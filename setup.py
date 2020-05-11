import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="zga",
    version="0.0.3",
    author="Aleksei Korzhenkov",
    author_email="oscypek@ya.ru",
    description="Prokaryotic genome assembly and annotation pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/laxeye/zga",
    packages=setuptools.find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'zga = zga.zga:main',
        ],
    },
    package_data={"zga": ["data/*"]},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires='>=3.6',
    install_requires='biopython'
)
