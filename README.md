MAVIS
=======

| branch  | status                                                                                      |
|---------|---------------------------------------------------------------------------------------------|
| master  | ![Build Status](https://www.bcgsc.ca/bamboo/plugins/servlet/wittified/build-status/MAV-MAV27)|
| develop | ![Build Status](https://www.bcgsc.ca/bamboo/plugins/servlet/wittified/build-status/MAV-MAV)|

Install Instructions for Developers
-------------------------------------------

Clone the repository and switch to the development branch

    >>> git clone https://svn.bcgsc.ca/bitbucket/scm/svia/mavis.git
    >>> cd mavis
    >>> git checkout develop

Set up a python virtual environment. If you are developing in python setting up with a virtual environment can be incredibly helpful. 
This can be used to generate the requirements.txt file that pip uses for install. Instructions for setting up the environment
are below

    >>> pip install virtualenv
    >>> virtualenv venv
    >>> source venv/bin/activate
    (venv) >>>

Install the MAVIS python package (currently need to use pip as well due to dependencies stored in svn)

    (venv) >>> python setup.py develop

Run the unit tests and compute code coverage. This will output html into ./coverage

    (venv) >>> python setup.py nosetests 

Make the user manual

    (venv) >>> cd docs
    (venv) >>> make html

The contents of the user manual can then be viewed by opening the build/html/index.html in any available
web browser (i.e. google-chrome, firefox, etc.)

Dependencies
----------------

Other than python3 and the python packages listed in the requirements.txt or setup.py file, the tool also requires samtools (for
sorting and indexing the output bam files) and an aligner. Currently the only aligner supported is blat. Testing
was performed using BLAT v.36x2