**Procedure for Regression Testing:**
1. Rename tests/\_test_gAperture_regr.py and tests/\_test_gMap_regr.py by removing the leading underscore and run nosetests in a local build. Confirm success, or justify failure and change the tests as needed. (Don't commit this change. This block of code takes a while to run and so makes Travis-CI time out / error.)

**Procedure for testing the build:**
1. Make a clean clone / fetch / checkout of the branch.
```
    git clone https://github.com/cmillion/gPhoton
    cd gPhoton
    git fetch
    git checkout -b BRANCH
    git pull origin BRANCH
```
2. Register the test distribution with testPyPI: `python setup.py register -r https://testpypi.python.org/pypi`
3. Build and upload the distribution file to testPyPI: `python setup.py sdist upload -r https://testpypi.python.org/pypi`
4. Install the new version from testPyPI: `pip install -i https://testpypi.python.org/pypi gPhoton`
5. Run any desired tests. A few gAperture commands from the User Guide should be fine.
6. `pip uninstall gPhoton`

**Procedure for building:**
The main PyPI server no longer requires a registration server after the first time. You just upload the file, so...
1. Merge BRANCH into Master.
2. Create a distribution file with: `python setup.py sdist`
3. Upload the file (securely): `twine upload dist/*`
4. Install gPhoton: `pip install gPhoton`

**After building:**
1. Move a copy of the `tar.gz` build file to the appropriate long term home in the MAST archive.
2. Tag this as a "release" on github. Include the build file as an attachment.
