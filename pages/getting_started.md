---
layout: page
title: "Getting Started"
teaser: "Starting i-PI is very easy!"
show_meta: false
header:
   image_fullwidth: "header-homepage.jpg"
permalink: "/resources/getting_started"
---

### Prerequisites

#### Mandatory:
 - A decent version of python >=2.5 (< 3)
 - A decent version of numpy
 
#### Optional:
 - A F90 compiler (gfortran preferred)
 - A program to plot time series and to visualize molecular structures
  (e.g. gnuplot and VMD)
 - Other codes which are compatible with i-PI (Quantum Espresso, CP2K,
   LAMMPS, etc.).
   
### How to install i-PI

---

If you want only to test i-PI you can download our VirtualBox image
which already contains i-PI and several clients properly
configured. 

---

Being written completely in Python, i-PI do not need a real
installation. To keep a low barrier for the new users, we decided to
provide a simple script that sets all the necessary environment
variable in a BASH console. Once you downloaded, you can move the
i-PI root directory where you prefer and then simply source the
```env.sh``` file:

{% highlight bash %}
$ cd <path_to_ipi>
$ source env.sh
{% endhighlight %}

Now i-PI and all his tools are available in the actual console. If you
want to avoid this operation in the future, you can just run the following command:

{% highlight bash %}
$ echo 'source <path_to_ipi>/env.sh' >> ~/.bashrc
{% endhighlight %}


Obviously the ```<path_to_ipi>``` must be replaced with the absolute
path to the i-PI folder.

At this point i-PI should work. To test if i-PI is working just run
the following command:

{% highlight bash %}
$ i-pi
{% endhighlight %}


and check that the output resemble the following:

{% highlight bash %}
Usage: i-pi [options] <input file>

i-pi: error: No input file name provided.
{% endhighlight %}


#### Installing an example driver
To do anything, i-PI needs to communicate with a client code able to
provide forces for a atomic configuration. We provide a simple
suite of potential energy surfaces coded into the "driver" that allows
the users to readily test i-PI.

To compile the driver a F90 compiler is needed. We recommend gfortran.

To prepare the driver executable, just enter the driver directory and write

{% highlight bash %}
$ cd <path_to_ipi>/driver
$ make
{% endhighlight %}


If you have a non standard linux installation or you want to use a
different compiler than gfortran, you can edit the Makefile contained
in this directory.


### First run with i-PI
To run i-PI one must know that i-PI works as a server that only moves
the atoms accordingly to the user input. To provide some valuable
results, we needs a client that connect to i-PI and computes forces,
energy and stresses on the atomic configurations provided by
i-PI. Thus, you will need two terminals: both must have been
configured to work with the i-PI framework (```source <path_to_ipi>/env.sh```).

As a first example move both terminals to the
```examples/driver/zundel-aimd``` directory and run the following:

**Terminal 1**
{% highlight bash %}
$ i-pi input.xml
{% endhighlight %}


After a few lines of output describing the i-PI set-up, you should get
the following text:
{% highlight bash %}
@ForceField: Starting the polling thread main loop
{% endhighlight %}

At this point i-PI is running and is waiting for a client. 
In **Terminal 2** write:

{% highlight bash %}
$ i-pi-driver -m zundel -u -h zundel-md & 
{% endhighlight %}

Now, the **Terminal 1** should report all the i-PI operations.

Congratulation! Your first simulation with i-PI is now running!

Please note that you can stop the client (**Terminal 2**) when you
want (also pressing ```Ctrl+c``` into the terminal directly). Since
you are running a single client, i-PI will stop too but will not
crash. If you decide to restart the client, running the last command
in the **Terminal 2**, the simulation will restart exactly from where
you stopped it!  Wonderful!

---

To go deeper in understanding how i-PI works you can check the
documentation and our tutorials.

---


