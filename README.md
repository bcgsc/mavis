<object type='image/svg+xml' data='docs/source/_static/acronym.svg'>
    <object type='image/svg+xml' data='_static/acronym.svg'>
    	<img src='docs/source/_static/acronym.svg' onerror='this.src="_static/acronym.svg"'>
    </object><br>
</object>

![master branch build Status](https://www.bcgsc.ca/bamboo/plugins/servlet/wittified/build-status/MAV-TEST) 
*(master)* 


![develop branch build status](https://www.bcgsc.ca/bamboo/plugins/servlet/wittified/build-status/MAV-TEST0) 
*(develop)* 

[MAVIS](http://mavis.bcgsc.ca) is python command-line tool for the post-processing of structural variant calls. 
The general [MAVIS](http://mavis.bcgsc.ca) pipeline consists of six main stages
 
- convert
- [cluster](http://mavis.bcgsc.ca/docs/latest/mavis.cluster.html#mavis-cluster)
- [validate](http://mavis.bcgsc.ca/docs/latest/mavis.validate.html#mavis-validate)
- [annotate](http://mavis.bcgsc.ca/docs/latest/mavis.annotate.html#mavis-annotate)
- [pairing](http://mavis.bcgsc.ca/docs/latest/mavis.pairing.html#mavis-pairing)
- [summary](http://mavis.bcgsc.ca/docs/latest/mavis.summary.html#mavis-summary)


## Getting started


There are 3 major steps to setting up and installing [MAVIS](http://mavis.bcgsc.ca)

1. **Install Aligner**

Before [MAVIS](http://mavis.bcgsc.ca) can be installed, the [non-python dependency](http://mavis.bcgsc.ca/docs/latest/about.html#non-python-dependencies) will need to be installed.
After these have been installed MAVIS itself can be installed through pip

The user may configure MAVIS to run using [blat](http://mavis.bcgsc.ca/docs/latest/glossary.html#term-blat) or [bwa mem](http://mavis.bcgsc.ca/docs/latest/glossary.html#term-bwa).

2. **Install MAVIS**

The easiest way to install [MAVIS](http://mavis.bcgsc.ca) is through the python package manager, pip

```
pip install mavis
```

This will install mavis and its python dependencies.

3. **Build reference files**

After [MAVIS](http://mavis.bcgsc.ca) is installed the [reference files](http://mavis.bcgsc.ca/docs/latest/reference.html) must be generated (or downloaded) before it can be run.

Once the above 3 steps are complete [MAVIS](http://mavis.bcgsc.ca) is ready to be run. See [running the pipeline](http://mavis.bcgsc.ca/docs/latest/pipeline.html).


## Help

If you have a question or issue that is not answered in the [FAQs](http://mavis.bcgsc.ca/docs/latest/faqs.html) please submit
an issue to our [github page](https://github.com/bcgsc/mavis/issues) or contact us by email at [mavis@bcgsc.ca](mailto:mavis@bcgsc.ca)
