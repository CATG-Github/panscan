from setuptools import setup, find_packages

setup(
    name="panscan",
    version="1.4.5",
    packages=find_packages(),
    package_data={
        "panscan": [
            "scripts/*",              # Include all files in the scripts directory
            "scripts/perlModules/*"    # Include all files in the perlModules directory
        ]
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "panscan=panscan.__main__:main",  # Assuming you have a `__main__.py` in `panscan/`
        ],
    },
)
