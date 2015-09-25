from setuptools import setup

setup(name='hpf',
      version='1.0',
      description='NYU Human Proteome Folding Code',
      author='NYU Bonneau Lab',
      license='GPLv2',
      author_email='rbonneau@nyu.edu',
      url='http://bonneaulab.bio.nyu.edu/',
      packages=['hpf', 'hpf.amnh', 'hpf.backup', 'hpf.bionet', 'hpf.blast', 'hpf.enrichment',
                'hpf.examples', 'hpf.fragmentor', 'hpf.function','hpf.hddb','hpf.interpro',
                'hpf.mcm','hpf.ncbi','hpf.pdb','hpf.plot','hpf.rosetta_abinit', 'hpf.scripts',
                'hpf.since_solved','hpf.structure_comparison','hpf.superfunc','hpf.test',
                'hpf.utilities','hpf.yrc'],
      install_requires=[
          'numpy',
          'scipy',
          'networkx',
          'biopython',
          'simplejson',
          'MySQL-python',
          'pymongo>=3.0.1',
          'SQLAlchemy==0.5.8',
      ]
     )