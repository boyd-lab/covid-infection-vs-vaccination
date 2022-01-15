#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

setup(
    author="Maxim Zaslavsky",
    author_email="maxim@maximz.com",
    name="covid_serology",
    description="Covid Serology",
    packages=find_packages(include=["covid_serology", "covid_serology.*"]),
    python_requires=">=3.7",
    version="0.0.1",
)
