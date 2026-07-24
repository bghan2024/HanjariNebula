from setuptools import setup, find_packages

setup(
    name="hanjari-nebula",
    version="1.0.0",
    description="High-precision macromolecular physical refinement and hydration compilation framework",
    author="Han Byeong-gu",
    author_email="hanbyeonggu@gmail.com",
    url="https://github.com/bghan2024/HanjariNebula",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "gemmi",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
