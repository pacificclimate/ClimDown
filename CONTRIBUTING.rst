Contributing to the ClimDown R package
======================================

Getting Started
---------------

- Create a `Github account`_.
- Fork the repository on Github at https://github.com/pacificclimate/ClimDown.
- Work on the code (see the next section)
- Send us a `pull request`_.

.. _Github account: https://github.com/signup/free
.. _pull request: https://help.github.com/articles/using-pull-requests/

Working on the code
-------------------

We highly recommend using the excellent `devtools`_ package for streamlining development. With `devtools` installed, one should be able to used ``devtools::load_all()``, ``devtools::document()`` and ``devtools::check()`` to load package, build the documentation, or check package validity (including running the tests suite) respectively.

.. _devtools: https://github.com/hadley/devtools

How to build the docs
---------------------

The package documentation is inline in the code. All of the manual pages are built by using ``roxygen2``. Make sure that you have ``roxygen2`` installed and loaded::
   
    $ echo 'roxygen2::roxygenize()' | R --no-save


Or::
   
    $ echo 'devtools::document()' | R --no-save


Don't code? No problem!
-----------------------

Even if you don't program for a living there are plenty of ways to help. Not only is the code open and collaborative, but so is the documentation and issue tracking. Anyone can help with these. If you can't program, consider helping with the following:

- If the documentation doesn't answer your questions, it probably doesn't answer many people's questions. Help us all out and write something that does.
- Take a look through the outstanding `"help wanted" issues`_, and see if you know any of the answers.
- If there are `open bug reports`_, see if you can reproduce the problem and verify that it exists. Having bug reports validated and/or clarified by multiple parties is extremely valuable.
- Tell us your story. If you have used ``ClimDown`` do to some epic climate downscaling, we would love to hear about it. Write a blog post and/or send an e-mail to the `package maintainer`_.

.. _"help wanted" issues: https://github.com/pacificclimate/ClimDown/labels/help%20wanted
.. _open bug reports: https://github.com/pacificclimate/ClimDown/labels/bug
.. _package maintainer: mailto:hiebert@uvic.ca

Releasing
---------

The package maintainer should take the following actions when creating a new release:

1. ``$ git checkout release``
2. ``$ git merge master``
3. Write an entry in the CHANGELOG file.
4. Update the ``Version`` and ``Date`` keys in the DESCRIPTION file.
5. Generate the docs (``$ echo 'roxygen2::roxygenize()' | R --no-save``)
6. git commit all changes.
7. git tag the commit with the version number.
