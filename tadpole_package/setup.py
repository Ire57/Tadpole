from setuptools import setup, find_packages

setup(
    name='tadpole',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        # Dependencies on requirements.txt
        'streamlit',
        'matplotlib',
        'scikit-learn',
        'numpy',
        'scipy',
        'biopython',
        'sbol3',
        'graphviz',
        'Pillow',
        'weasyprint',
        'ipython',
        'requests',
        'flask', 
    ],
    entry_points={
        'console_scripts': [
            'tadpole = tadpole.cli:main',
        ],
    },
    author='iGem-UB 2025',
    description='A command-line toolkit for designing RNA switches.',
    long_description='Given a Functional RNA Element (FRE) and a Conformational RNA Element (CRE), our software finds a linker sequence and designs a system that supports two conformations: an OFF state with FRE-CRE hybridisation, and an ON state where the FRE preserves its functional structure.',
    long_description_content_type='text/markdown',
    url='https://gitlab.igem.org/2025/software-tools/barcelona-ub',
)