=================
Contributing code
=================

The following guidelines for developers should help us collaborate on the project and contribute code efficiently.
They are work in progress and suggestions are welcome.


Using pull requests
===================

To add code to the :code:`master` branch of the repository, use the pull request mechanism.
Small bug fixes with no potential to break things and with no need for feedback can be commited directly.
Here is the official documentation on pull requests:

https://help.github.com/articles/using-pull-requests

The idea is to rougly follow what is called the Github Flow, in a "shared repository model" using the terminology above.
This is a Git workflow that uses a single "master" branch and feature branches that are used to create pull requests.
Here is an overview of how it works:

http://guides.github.com/overviews/flow/

For a more in-depth explanation, see for example this article:

http://scottchacon.com/2011/08/31/github-flow.html

While we do not really "deploy", many aspects of this workflow are still very useful.
Note that all the operations can also be performed on the web, as explained here:

https://help.github.com/articles/github-flow-in-the-browser

In principle, one can also use the "fork and pull" model, but we probably have no specific reason to do so in our case.
If it works better for you personally, feel free to use it, of course, but please make sure that you keep the fork private.
From the perspective of the "master" branch, there is no difference between pull requests from the same repository and from a fork.
However, if the feature branch is in the shared main repository, multiple people can contribute to it, either directly or using pull requests.

You can open a pull request early for a branch with an unfinished feature to get feedback.
Mark it as work in progress with the WIP tag so that it is clear it is not ready to merge.

To ask for feedback from a specific person, mention them in the pull request so that they get notified.
If you want a specific person to check and merge the pull request, you can assign it to them.
If you want to indicate that you plan to merge it yourself, you can assign the pull request to yourself.
Please do not merge pull requests that are assigned to someone else.
People with write access to the main repository should probably usually merge their own pull requests, after potential feedback from other in the comments.
Pull requests from people with read-only access need to be merged by someone with write access, after.

If you have seen some else's pull request and think it should be merged, indicate it with a brief comment, as there is no voting mechanism on Github.
You can even use one of the `emoji`_ that Github supports.
You can go to town with these, but at least :+1: and :-1: are a quick way to say "I think this should go in" and "I don't think this should go in".


Meld
====

`Meld`_ is a visual diff and merge tool that can work with files as well as with whole directories.
It can do both two- and three-way comparisons and has support for version control systems.

On Linux, you can get it from you package manager, on Mac, both Homebrew and Macports have it.
They say it should work on Windows as well, but I have never tried.

It is a really cool tool, I have not found anything free/open that would be as good.
There might be commercial tools that are as good or better, but they are usually part of some development environment.

If you set it as your Git difftool, it can be used to show the diff between arbitrary commits or branches.
That can be useful to review pull requests (diff against the remote branch), as the side-by-side diff view is much better for complicated changes.
If anyone is interested, I can share my settings and workflow (Ondrej).


Issues tracker
==============

Use the "Issues" feature of Github to track bugs and enhancement proposals, but also questions or other issues we might want to discuss. Mentioning an open issue in a commit or pull request in one of the ways listed `here`_ can be used to close it.

.. _here: https://help.github.com/articles/closing-issues-via-commit-messages


Code style
==========

Do we want to set (approximate?) code style rules? `PEP 8`_ is probably a pretty good starting point. What about docstrings, do we want to follow `PEP 257`_, something else, or leave it more open?


.. _emoji: http://www.emoji-cheat-sheet.com/
.. _Meld: http://meldmerge.org/
.. _PEP 8: http://www.python.org/dev/peps/pep-0008/
.. _PEP 257: http://www.python.org/dev/peps/pep-0257/
