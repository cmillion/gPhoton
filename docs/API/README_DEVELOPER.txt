A minimum requirement for ReadTheDocs API to work is that this command must execute without errors, and preferably without warnings as well.  Common problems are addition of imports of non-standard modules (and submodules) that are not Mocked in conf.py, and improper indentation of docstrings:

sphinx-build -b html . _build/html
