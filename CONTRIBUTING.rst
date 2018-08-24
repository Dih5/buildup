.. highlight:: shell

============
Contributing
============

Despite this is a rather small project, contributions are very welcome.

You can help this project in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/Dih5/buildupfactor/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Feel free to complete whatever leaves room for improvement: official docs, docstrings, rst files...

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/Dih5/buildupfactor/issues.

If you want to propose a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Do you want to become a volunteer?! Here are some instructions to set up `buildupfactor` for local development.

1. Fork the `buildupfactor` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/buildupfactor.git

3. Install your local copy into a virtualenv. From the repository root you can::

    $ virtualenv venv
    $ source venv/bin/activate
    $ pip install -e .

   You should also install the dev requisites::

    $ pip install -r requirements_dev.txt

   Don't forget to source your virtualenv to work with the development version.

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, you should check at least that they pass the flake8 coding style requirements and
   the tests with your python version::

    $ flake8 buildupfactor tests
    $ python setup.py test or py.test

   If you have multiple python versions installed, you can run the tests in all them with tox, which is more
   recommended::

    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. New tests should be added if they are relevant.
2. If the pull request adds functionality, the docs should be updated. Check the docstrings and the .rst in the docs
   folder (or those linked there).
3. The pull request should work for Python 2.7, 3.4, 3.5 and 3.6. Check
   https://travis-ci.org/Dih5/buildupfactor/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

Run::

    $ make help

in your local repo to find a set of Very Fitting Directives.

To run a subset of tests::

    $ py.test tests/test_foo.py

To run a single test::

    $ py.test tests/test_foo.py::test_bar

Deploying
---------

A reminder for the Volunteers Facing Deployment.
First make sure all relevant changes have been committed (including an entry in HISTORY.rst).
Then run::

$ bumpversion XXX # for XXX in: major / minor / patch
$ git push
$ git push --tags

Travis will then deploy to PyPI if tests pass.
