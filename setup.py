from setuptools import setup, find_packages

setup(
    name="pyrbd_plusplus",
    version="1.0.0",
    packages=find_packages(include=["pyrbd_plusplus", "pyrbd_plusplus.*"]),
    install_requires=[
        "pandas",
        "networkx",
        "tqdm",
        "loguru"
    ],
)