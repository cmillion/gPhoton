**Procedure for Regression Testing:**
1. Remove commented block in tests/test_gAperture_regr.py and run nosetests in a local build. Confirm success, or justify failure and change the tests as needed. (Don't commit this change. This block of code takes a while to run and so makes Travis-CI time out / error.)

**Procedure for testing the build:**
1. Make a clean clone / fetch / checkout of the branch.
    git clone https://github.com/cmillion/gPhoton
    cd gPhoton
    git fetch
    git checkout -b BRANCH
    git pull origin BRANCH
2. Run `python setup.py sdist` from the checkout directory to generate a dist file (in ./dist/) with the correct naming convention.
3. Copy that tar.gz file into your Public Dropbox directory and copy the link.
4. Replace download_url in setup.py with that link.
5. Run `python setup.py register -r pypitest`
6. Rerun `python setup.py sdist` and replace the version in the Public Dropbox with the resulting tar.gz.
7. `pip install -i https://testpypi.python.org/pypi gPhoton`
8. Run any tests.
9. `pip uninstall gPhoton`

**Procedure for building:**
Identical to above except that `pypitest` becomes `pypi` everywhere, and the download_url points to MAST, which means that someone needs to put the tar.gz in the appropriate place there.

`python setup.py register -r pypi`

`python setup.py sdist`

`pip install -i https://pypi.python.org/pypi gPhoton`
