Hier werden alle arbeitsscripts der Membrane Analysis Toolbox verwaltet. Auch plots und specifische scripts für das Carbon Paper werden hier abegelegt


## How to use the scripts?
- (if necessary) log into workstation: `ssh -X st123456@129.69.167.51` (-X is optional to forward plots to the local machine)
- `cd /bigpool/users/st166545/carbon_paper_MAT` to get to the scripts
- `source .venv/bin/activate` to activate the virtual environment in which the library is installed
- (optional) `which python` to ensure that the "python" refers to the virtual environment python under "/bigpool/users/st166545/MembraneAnalysisToolbox/.venv/bin/python"
- run the python files from within the activated virtual environment: `python interactive_PoreAnalysis.py`