
---

# Refer to the libindi library #

**The driver is now included in the [libindi library](http://www.indilib.org) as a 3rd party driver (_indi-eqmod_). I now mainly put new development stuff on the [libindi sourceforge](http://sourceforge.net/p/indi/code/HEAD/tree/trunk/). I use this site to eventually put some specific stuff concerning my home made controller. So do not use it as it is rarely up-to-date**

---


This <a href='http://indilib.org/'>INDI</a> driver interacts with a mount controller using the <a href='http://code.google.com/p/skywatcher/'>Skywatcher Protocol</a> through a serial link. It can directly drive Skywatcher and similar mounts through a serial port <b>after a conversion of signal levels</b> (see <a href='http://eq-mod.sourceforge.net/eqdirect2.htm'>EQDIRECT</a> interface), or any other mount hardware which uses this protocol. This project is similar in its goal to the <a href='http://eq-mod.sourceforge.net/'>EQMOD/EQASCOM</a> project, but for the Linux/INDI world. Thus the name of the INDI device is "EQMod Mount", and the driver executable is <tt>indi_eqmod_telescope</tt>.

<b>The project is in a very early stage</b> since, as of today, <b>it has even not been tested with a real Skywatcher mount</b>, only with a <a href='http://www.flickr.com/photos/87297028@N04/7993769744/'>home-made mount controller</a> driving an <a href='http://www.flickr.com/photos/87297028@N04/7993762889/'>old modded EQ5</a>. However this home-made controller performs well with the EQMOD/EQASCOM driver itself under Windows/Wine. If you plan to test this driver with your own mount, first <b>disassemble your telescope from the mount</b> to prevent potential damages. If you succeed to use the driver, please send me an email.

Current features are:
<ul>
<li>Goto/Slew at maximum microstep resolution</li>
<li>Independent slew speeds for both axes variable between x0.05 to x800 of the sidereal rate (step x0.05)</li>
<li>Sidereal, lunar, solar and custom trackrates</li>
<li>Pulse-guiding</li>
</ul>
Normal syncs and Nearest point alignment are functionnal, N-Star alignment mode is in a testing phase.

**commit(2013/08/26)** _Triangulation used for N-Star alignment (code from [Computational Geometry in C](http://cs.smith.edu/~orourke/books/ftp.html), thanks). Corrected many bugs concerning time/position computations. Polar misalignment estimation (paper from [R. Johnson](http://www.whim.org/nebula/math/pdf/twostar.pdf), thanks. Some mount tools available on [indi-mount-tools](https://code.google.com/p/indi-mount-tools/) (including a GUI Pad enabling the use of a joypad, code from [qjoypad](http://qjoypad.sourceforge.net/), thanks)._

**commit(2013/02/26)** _Preset Speeds. Full Debug and logging support. Log files are in_`/tmp/indi_eqmod_telescope_%TIMESTAMP%.log`. Logging is not included by default (see INSTALL)

**commit(2012/12/25)** _Fixed exceptions handling issues._

**commit(2012/12/14)** _N-Star Alignment code has been written, but **not yet tested** in real conditions. See AlignmentHowto Wiki page for brief explanations._

**commit(2012/12/04)** _Fixed major properties management issues (thanks to R. James), Enabled Alt/Az mount type 114GT for testing purpose, Compliant to Indi version 0.9.6 (no more Req properties)_

**commit(2012/10/20)** _inherits Guider Interface, alignment disabled (no more syncs), park disabled._

A screenshot of the Main Control panel in KStars (<a href='http://www.flickr.com/photos/87297028@N04/sets/72157631552640091/'>other screenshots on flickr</a>)

<a href='http://www.flickr.com/photos/87297028@N04/7993557298/' title='Capture-INDIControlPanel de geehalel, sur Flickr'><img src='http://farm9.staticflickr.com/8041/7993557298_b4b1f559d7.jpg' alt='Capture-INDIControlPanel' width='500' height='491' /></a>
