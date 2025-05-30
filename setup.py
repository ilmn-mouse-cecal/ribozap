from setuptools import setup, find_packages

setup(
    name="ribozap",
    version="0.0.1",
    description="Design probes to deplete rRNA content",
    author="Samuel Bunga",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.10",
    install_requires=[],
    entry_points={
        "console_scripts": [
            "ribozap=ribozap.cli:main",
        ]
    },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
