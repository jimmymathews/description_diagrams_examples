import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="datacube",
    version="0.0.1",
    author="Jimmy Mathews",
    author_email="mathewj2@mskcc.org",
    description="Feature matrix data analysis and modeling using Steiner graphs.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
