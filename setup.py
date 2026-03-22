from setuptools import setup, find_packages

setup(
    name="waternes",
    version="0.1",
    packages=find_packages("src"),
    package_dir={"": "src"},
    entry_points={
        "console_scripts": [
            "waternes=waternes.cli:main"
        ]
    },
)
