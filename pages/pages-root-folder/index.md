---
#
# Use the widgets beneath and the content will be
# inserted automagically in the webpage. To make
# this work, you have to use › layout: frontpage
#
layout: page
header:
  image_fullwidth: header-homepage.jpg

#widget1:
#  title: "Blog & Portfolio"
#  url: 'http://phlow.github.io/feeling-responsive/blog/'
#  image: widget-1-302x182.jpg
#  text: 'Every good portfolio website has a blog with fresh news, thoughts and develop&shy;ments of your activities. <em>Feeling Responsive</em> offers you a fully functional blog with an archive page to give readers a quick overview of all your posts.'
#widget2:
#  title: "Why use this theme?"
#  url: 'http://phlow.github.io/feeling-responsive/info/'
#  text: '<em>Feeling Responsive</em> is heavily customizable.<br/>1. Language-Support :)<br/>2. Optimized for speed and it&#39;s responsive.<br/>3. Built on <a href="http://foundation.zurb.com/">Foundation Framework</a>.<br/>4. Seven different Headers.<br/>5. Customizable navigation, footer,...'
#  video: '<a href="#" data-reveal-id="videoModal"><img src="http://phlow.github.io/feeling-responsive/images/start-video-feeling-responsive-302x182.jpg" width="302" height="182" alt=""/></a>'
#widget3:
#  title: "Download Theme"
#  url: 'https://github.com/Phlow/feeling-responsive'
#  image: widget-github-303x182.jpg
#  text: '<em>Feeling Responsive</em> is free and licensed under a MIT License. Make it your own and start building. Grab the <a href="https://github.com/Phlow/feeling-responsive/tree/bare-bones-version">Bare-Bones-Version</a> for a fresh start or learn how to use it with the <a href="https://github.com/Phlow/feeling-responsive/tree/gh-pages">education-version</a> with sample posts and images. Then tell me via Twitter <a href="http://twitter.com/phlow">@phlow</a>.'

#
# Use the call for action to show a button on the frontpage
#
# To make internal links, just use a permalink like this
# url: /getting-started/
#
# To style the button in different colors, use no value
# to use the main color or success, alert or secondary.
# To change colors see sass/_01_settings_colors.scss
#
#
#callforaction:
#  url: https://tinyletter.com/feeling-responsive
#  text: Inform me about new updates and features ›
#  style: alert
#
permalink: /index.html
#
# This is a nasty hack to make the navigation highlight
# this page as active in the topbar navigation
#
homepage: true
---

About i-PI
==========
i-PI is a universal force engine interface
written in Python, designed to be used together with an ab-initio
evaluation of the interactions between the atoms. The main goal is to
decouple the problem of evolving the ionic positions to sample the
appropriate thermodynamic ensemble and the problem of computing the
inter-atomic forces.

The implementation is based on a client-server paradigm, where i-PI
acts as the server and deals with the propagation of the nuclear
motion, whereas the calculation of the potential energy, forces and
the potential energy part of the pressure virial is delegated to one
or more instances of an external code, acting as clients.


i-PI is free software, distributed under a dual MIT/GPLv3 licence. You
are welcome to dowload, use, modify and redistribute it. If you find it
useful for your research, a citation to
[Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 185, 1019-1026 (2014)](http://dx.doi.org/10.1016/j.cpc.2013.10.027)
would be appreciated.

As today, i-PI can be used together
with
[CP2K](https://www.cp2k.org/),
[Lammps](http://lammps.sandia.gov/),
[QuatumEspresso](http://quantum-espresso.org),
[Siesta](http://departments.icmab.es/leem/siesta/) 
